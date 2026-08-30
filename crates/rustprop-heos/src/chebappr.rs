//! General piecewise-Chebyshev approximation machinery (PLAN.md 4.6 (H,S)) —
//! port of `ChebyshevApproximation1D` and the `detail` helpers of
//! `include/CoolProp/superancillary/superancillary.h`: derivative-expansion
//! coefficients (Mason & Handscomb Eq. 2.52), the transposed colleague
//! matrix, LAPACK-style balancing (arXiv:1401.5766 alg. 3), eigenvalue
//! extrema, monotonic-interval construction, and the general multi-root
//! `get_x_for_y`.
//!
//! Two logged deviations (PLAN.md Decisions):
//! - eigenvalues come from the classic EISPACK/Numerical-Recipes `hqr`
//!   Francis double-shift QR instead of Eigen's `EigenSolver` — the same
//!   algorithm family; extrema agree to conditioning-limited precision, and
//!   every downstream root is re-refined against the expansion itself, so
//!   ULP-level extrema differences cannot propagate;
//! - bracketed roots use bit-tolerance bisection + midpoint in place of
//!   boost TOMS748 (established stand-in, same bracket semantics).

use crate::superancillary::{Cheb1D, OwnedInterval, clenshaw, eval_piecewise};

// ---------------------------------------------------------------------------
// detail:: helpers
// ---------------------------------------------------------------------------

/// Upstream `detail::get_LU_matrices`: L converts node values to
/// coefficients, U converts coefficients to node values.
#[allow(clippy::needless_range_loop)] // symmetric fill writes [j][k] and [k][j]
pub(crate) fn lu_matrices(n: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let nf = n as f64;
    let mut l = vec![vec![0.0; n + 1]; n + 1];
    let mut u = vec![vec![0.0; n + 1]; n + 1];
    for j in 0..=n {
        for k in j..=n {
            let p_j: f64 = if j == 0 || j == n { 2.0 } else { 1.0 };
            let p_k: f64 = if k == 0 || k == n { 2.0 } else { 1.0 };
            let cosjpikn = ((j as f64) * std::f64::consts::PI * (k as f64) / nf).cos();
            l[j][k] = 2.0 / (p_j * p_k * nf) * cosjpikn;
            l[k][j] = l[j][k];
            u[j][k] = cosjpikn;
            u[k][j] = u[j][k];
        }
    }
    (l, u)
}

/// Chebyshev-Lobatto nodes `cos(pi j/N), j = 0..=N` mapped to `[xmin, xmax]`
/// (upstream `ChebyshevExpansion::get_nodes_realworld`; descending in x).
pub(crate) fn nodes_realworld(n: usize, xmin: f64, xmax: f64) -> Vec<f64> {
    (0..=n)
        .map(|j| {
            let node = ((j as f64) * std::f64::consts::PI / (n as f64)).cos();
            ((xmax - xmin) * node + (xmax + xmin)) * 0.5
        })
        .collect()
}

/// Row-dot-product matrix * vector (Eigen gemv stand-in; accumulation order
/// differs from Eigen at ULP level — see the invlnp fit note in PLAN.md).
pub(crate) fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(v).map(|(a, b)| a * b).sum())
        .collect()
}

/// Upstream `ChebyshevExpansion::do_derivs` for one derivative order:
/// coefficients of dy/dx as a Chebyshev expansion on the same interval
/// (Mason & Handscomb p. 34, Eq. 2.52, with the range rescaling).
#[allow(clippy::needless_range_loop)] // index form mirrors the upstream loops
fn deriv_coeffs(coef: &[f64], xmin: f64, xmax: f64) -> Vec<f64> {
    let n = coef.len() - 1; // degree
    let nd = n - 1; // degree of the derivative expansion
    let mut cd = vec![0.0; n];
    for r in 0..=nd {
        for k in (r + 1)..=n {
            if (k - r) % 2 == 1 {
                cd[r] += 2.0 * (k as f64) * coef[k];
            }
        }
        if r == 0 {
            cd[r] /= 2.0;
        }
        cd[r] /= (xmax - xmin) / 2.0;
    }
    cd
}

/// Upstream `detail::companion_matrix_transposed` (colleague matrix of a
/// Chebyshev expansion, transposed; upper Hessenberg).
fn companion_matrix_transposed(coeffs: &[f64]) -> Vec<Vec<f64>> {
    let n = coeffs.len() - 1; // degree
    assert!(
        n >= 2,
        "companion matrix needs a derivative expansion of degree >= 2 (got {} coefficients); \
         upstream would index out of bounds here too",
        coeffs.len()
    );
    let mut a = vec![vec![0.0; n]; n];
    // First col
    a[1][0] = 1.0;
    // Last col
    for i in 0..n {
        a[i][n - 1] = -coeffs[i] / (2.0 * coeffs[n]);
    }
    a[n - 2][n - 1] += 0.5;
    // All the other cols
    for j in 1..(n - 1) {
        a[j - 1][j] = 0.5;
        a[j + 1][j] = 0.5;
    }
    a
}

/// Upstream `detail::balance_matrix` (arXiv:1401.5766 algorithm #3, p = 2,
/// radix 2), including the non-finite/zero row-column guard.
#[allow(clippy::needless_range_loop)] // index form mirrors the upstream loops
fn balance_matrix(a: &mut [Vec<f64>]) {
    let n = a.len();
    let beta = 2.0f64;
    let mut iter = 0;
    loop {
        let mut converged = true;
        for i in 0..n {
            let mut c: f64 = (0..n).map(|k| a[k][i] * a[k][i]).sum::<f64>().sqrt();
            let mut r: f64 = (0..n).map(|k| a[i][k] * a[i][k]).sum::<f64>().sqrt();
            if !c.is_finite() || !r.is_finite() || c == 0.0 || r == 0.0 {
                continue;
            }
            let s = c.powi(2) + r.powi(2);
            let mut f = 1.0;
            while c < r / beta {
                c *= beta;
                r /= beta;
                f *= beta;
            }
            while c >= r * beta {
                c /= beta;
                r *= beta;
                f /= beta;
            }
            // `f` is a product of radix-2 steps, so an extreme row/column
            // norm ratio can drive it to infinity (or to zero the other way)
            // before the two while-loops balance out. Scaling by a non-finite
            // or zero factor would write inf/NaN across a row and a column
            // and hand `hqr` a matrix with no meaningful eigenvalues; leaving
            // the row unscaled is what the existing non-finite c/r guard
            // above already does in the same situation. Real coefficient data
            // is nowhere near this range.
            if c.powi(2) + r.powi(2) < 0.95 * s && f.is_finite() && f > 0.0 {
                converged = false;
                for k in 0..n {
                    a[k][i] *= f;
                }
                for k in 0..n {
                    a[i][k] /= f;
                }
            }
        }
        iter += 1;
        if iter > 50 || converged {
            break;
        }
    }
}

/// Eigenvalues of a real upper Hessenberg matrix by the Francis double-shift
/// QR iteration — the classic EISPACK `hqr` (as given in Numerical Recipes,
/// 2nd ed., §11.6), transcribed 1-based to keep the port mechanical.
/// Upstream evaluates the same balanced matrix through Eigen's RealSchur
/// (the same algorithm family) — see the module docs for the logged
/// deviation. NR's `do {...} while (l < nn-1)` collapses to: deflation
/// branches exit the inner loop (their `l` makes the condition false),
/// the QR-step branch always continues (`l <= nn-2` there).
/// Returns `(re, im)` pairs, or `None` if an eigenvalue fails to converge in
/// 30 iterations.
#[allow(clippy::needless_range_loop)]
fn hqr(h: &[Vec<f64>]) -> Option<Vec<(f64, f64)>> {
    let n = h.len();
    // 1-based working copy: a[1..=n][1..=n]
    let mut a = vec![vec![0.0f64; n + 1]; n + 1];
    for i in 0..n {
        for j in 0..n {
            a[i + 1][j + 1] = h[i][j];
        }
    }
    let mut wr = vec![0.0f64; n + 1];
    let mut wi = vec![0.0f64; n + 1];

    let sign = |a: f64, b: f64| if b >= 0.0 { a.abs() } else { -a.abs() };

    let mut anorm = 0.0;
    for i in 1..=n {
        for j in (i - 1).max(1)..=n {
            anorm += a[i][j].abs();
        }
    }
    let mut nn = n;
    let mut t = 0.0;
    while nn >= 1 {
        let mut its = 0;
        loop {
            // Look for a single small subdiagonal element
            let mut l = 1;
            for cand in (2..=nn).rev() {
                let mut s = a[cand - 1][cand - 1].abs() + a[cand][cand].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if a[cand][cand - 1].abs() + s == s {
                    a[cand][cand - 1] = 0.0;
                    l = cand;
                    break;
                }
            }
            let mut x = a[nn][nn];
            if l == nn {
                // One root found
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
                break;
            }
            let mut y = a[nn - 1][nn - 1];
            let mut w = a[nn][nn - 1] * a[nn - 1][nn];
            if l == nn - 1 {
                // Two roots found
                let p = 0.5 * (y - x);
                let q = p * p + w;
                let mut z = q.abs().sqrt();
                x += t;
                if q >= 0.0 {
                    z = p + sign(z, p);
                    wr[nn - 1] = x + z;
                    wr[nn] = x + z;
                    if z != 0.0 {
                        wr[nn] = x - w / z;
                    }
                    wi[nn - 1] = 0.0;
                    wi[nn] = 0.0;
                } else {
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn] = z;
                    wi[nn - 1] = -z;
                }
                nn -= 2;
                break;
            }
            // No roots found. Continue iteration.
            if its == 30 {
                return None;
            }
            if its == 10 || its == 20 {
                // Exceptional shift
                t += x;
                for i in 1..=nn {
                    a[i][i] -= x;
                }
                let s = a[nn][nn - 1].abs() + a[nn - 1][nn - 2].abs();
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            its += 1;
            // Form shift and look for two consecutive small subdiagonal
            // elements.
            let mut p = 0.0;
            let mut q = 0.0;
            let mut r = 0.0;
            let mut m = nn - 2;
            for cand in (l..=nn - 2).rev() {
                let z = a[cand][cand];
                let rr = x - z;
                let ss = y - z;
                p = (rr * ss - w) / a[cand + 1][cand] + a[cand][cand + 1];
                q = a[cand + 1][cand + 1] - z - rr - ss;
                r = a[cand + 2][cand + 1];
                let s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                m = cand;
                if cand == l {
                    break;
                }
                let u = a[cand][cand - 1].abs() * (q.abs() + r.abs());
                let v =
                    p.abs() * (a[cand - 1][cand - 1].abs() + z.abs() + a[cand + 1][cand + 1].abs());
                if u + v == v {
                    break;
                }
            }
            for i in (m + 2)..=nn {
                a[i][i - 2] = 0.0;
                if i != m + 2 {
                    a[i][i - 3] = 0.0;
                }
            }
            // Double QR step on rows l..=nn and columns m..=nn
            for k in m..=(nn - 1) {
                if k != m {
                    p = a[k][k - 1];
                    q = a[k + 1][k - 1];
                    r = 0.0;
                    if k != nn - 1 {
                        r = a[k + 2][k - 1];
                    }
                    x = p.abs() + q.abs() + r.abs();
                    if x != 0.0 {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }
                let s = sign((p * p + q * q + r * r).sqrt(), p);
                if s != 0.0 {
                    if k == m {
                        if l != m {
                            a[k][k - 1] = -a[k][k - 1];
                        }
                    } else {
                        a[k][k - 1] = -s * x;
                    }
                    p += s;
                    x = p / s;
                    y = q / s;
                    let z = r / s;
                    q /= p;
                    r /= p;
                    // Row modification
                    for j in k..=nn {
                        p = a[k][j] + q * a[k + 1][j];
                        if k != nn - 1 {
                            p += r * a[k + 2][j];
                            a[k + 2][j] -= p * z;
                        }
                        a[k + 1][j] -= p * y;
                        a[k][j] -= p * x;
                    }
                    let mmin = if nn < k + 3 { nn } else { k + 3 };
                    // Column modification
                    for i in l..=mmin {
                        p = x * a[i][k] + y * a[i][k + 1];
                        if k != nn - 1 {
                            p += z * a[i][k + 2];
                            a[i][k + 2] -= p * r;
                        }
                        a[i][k + 1] -= p * q;
                        a[i][k] -= p;
                    }
                }
            }
        }
        if nn == 0 {
            break;
        }
    }
    Some((1..=n).map(|i| (wr[i], wi[i])).collect())
}

// ---------------------------------------------------------------------------
// ChebyshevApproximation1D
// ---------------------------------------------------------------------------

/// Upstream `MonotonicExpansionMatch`.
#[derive(Debug, Clone)]
pub(crate) struct MonotonicExpansionMatch {
    pub idx: usize,
    pub ymin: f64,
    pub ymax: f64,
    pub xmin: f64,
    pub xmax: f64,
}

impl MonotonicExpansionMatch {
    fn contains_y(&self, y: f64) -> bool {
        y >= self.ymin && y <= self.ymax
    }
}

/// Upstream `IntervalMatch`. `xmin`/`xmax` mirror the upstream struct and
/// feed no ported consumer yet.
#[derive(Debug, Clone)]
pub(crate) struct IntervalMatch {
    pub expansioninfo: Vec<MonotonicExpansionMatch>,
    #[allow(dead_code)]
    pub xmin: f64,
    #[allow(dead_code)]
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
}

impl IntervalMatch {
    fn contains_y(&self, y: f64) -> bool {
        y >= self.ymin && y <= self.ymax
    }
}

/// Upstream `ChebyshevApproximation1D` over owned expansions: extrema by
/// colleague-matrix eigenvalues, monotonic intervals, and general
/// multi-root inversion.
pub(crate) struct ChebApprox1d {
    pub expansions: Vec<OwnedInterval>,
    /// Extrema locations (kept for tests and upstream parity).
    #[allow(dead_code)]
    pub x_at_extrema: Vec<f64>,
    monotonic_intervals: Vec<IntervalMatch>,
}

/// Threshold below which an eigenvalue's imaginary part counts as zero
/// (upstream `thresh_imag`).
const THRESH_IMAG: f64 = 1e-15;

impl ChebApprox1d {
    pub fn new(expansions: Vec<OwnedInterval>) -> Self {
        let x_at_extrema = determine_extrema(&expansions, THRESH_IMAG);
        let monotonic_intervals = build_monotonic_intervals(&expansions, &x_at_extrema);
        ChebApprox1d {
            expansions,
            x_at_extrema,
            monotonic_intervals,
        }
    }

    pub fn eval(&self, x: f64) -> f64 {
        eval_piecewise(&self.expansions, x)
    }

    /// Domain lower edge (upstream `xmin()`).
    pub fn xmin(&self) -> f64 {
        self.expansions[0].xmin
    }

    /// Domain upper edge (upstream `xmax()`).
    pub fn xmax(&self) -> f64 {
        self.expansions[self.expansions.len() - 1].xmax
    }

    /// Upstream `get_x_for_y`: all solutions of `y(x) = y` across the
    /// monotonic intervals (each interval piece is monotonic, so endpoint
    /// bracketing is exhaustive).
    pub fn get_x_for_y(&self, y: f64, bits: u32, max_iter: usize, boundsftol: f64) -> Vec<f64> {
        let mut solns = Vec::new();
        for interval in &self.monotonic_intervals {
            if interval.contains_y(y) {
                for ei in &interval.expansioninfo {
                    if ei.contains_y(y) {
                        let e = &self.expansions[ei.idx];
                        solns.push(solve_for_x(
                            e, y, ei.xmin, ei.xmax, bits, max_iter, boundsftol,
                        ));
                    }
                }
            }
        }
        solns
    }
}

/// Upstream `ChebyshevExpansion::solve_for_x_count` (bounds handling +
/// bracketed root at the given bit tolerance).
fn solve_for_x<C: Cheb1D>(
    e: &C,
    y: f64,
    a: f64,
    b: f64,
    bits: u32,
    max_iter: usize,
    boundsytol: f64,
) -> f64 {
    let f = |x: f64| clenshaw(e, x) - y;
    let fa = f(a);
    if fa.abs() < boundsytol {
        return a;
    }
    let fb = f(b);
    if fb.abs() < boundsytol {
        return b;
    }
    bisect_bits(f, a, b, fa, fb, bits, max_iter)
}

/// TOMS748 stand-in: bisect the bracket until its relative width reaches the
/// bit tolerance `2^(1-bits)` (boost `eps_tolerance(bits)` semantics), then
/// return the midpoint, exactly as upstream midpoints the TOMS748 bracket.
pub(crate) fn bisect_bits<F: FnMut(f64) -> f64>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    mut fa: f64,
    fb: f64,
    bits: u32,
    max_iter: usize,
) -> f64 {
    // The bracket check is a development aid for genuine bracketing bugs, so
    // it only applies where "same sign" is meaningful. A non-finite residual
    // reaches here from a non-finite REQUEST — `props_si("Q", "S", NaN,
    // "Hmass", -inf, "R134a")` is one — and `fa * fb` is then NaN, which
    // fails `<= 0.0` and fired the assertion. Release builds compile
    // `debug_assert!` out and return the interval midpoint, so the bare form
    // made a debug build panic where a release build answered: a consumer
    // building in dev profile saw crashes the shipped profile does not have,
    // and `ci.yml`'s debug gate was one seed away from catching it.
    // Non-finite residuals are left to the loop, which handles them by
    // returning a finite midpoint of a finite interval.
    debug_assert!(
        !(fa.is_finite() && fb.is_finite()) || fa * fb <= 0.0,
        "bisect_bits called without a bracket: f(a) = {fa}, f(b) = {fb}"
    );
    let eps = (2.0f64).powi(1 - bits as i32);
    let mut fb = fb;
    for _ in 0..max_iter {
        if (b - a).abs() <= eps * a.abs().max(b.abs()) {
            break;
        }
        let m = 0.5 * (a + b);
        if m <= a.min(b) || m >= a.max(b) {
            break; // machine tightness
        }
        let fm = f(m);
        if fm == 0.0 {
            return m;
        }
        if (fa < 0.0) == (fm < 0.0) {
            a = m;
            fa = fm;
        } else {
            b = m;
            fb = fm;
        }
    }
    let _ = fb;
    0.5 * (a + b)
}

/// Upstream `determine_extrema`: roots of each expansion's derivative via
/// balanced-colleague-matrix eigenvalues, filtered to real values in
/// [-1, 1] and mapped to real-world x; sorted ascending.
fn determine_extrema(expansions: &[OwnedInterval], thresh_im: f64) -> Vec<f64> {
    let mut x_at_extrema = Vec::new();
    for expan in expansions {
        // Derivative coefficients for x in [-1, 1] (the rescale by the
        // interval width only scales all coefficients uniformly and does not
        // move the roots).
        let mut cd = deriv_coeffs(&expan.coef, expan.xmin, expan.xmax);
        // Right-trim coefficients equal to zero, exactly as upstream
        // (head(ilastnonzero) — the index of the last nonzero — dropping
        // that element too; verbatim bug-for-bug).
        let mut ilastnonzero = cd.len() - 1;
        for i in (0..cd.len()).rev() {
            if cd[i].abs() != 0.0 {
                ilastnonzero = i;
                break;
            }
        }
        if ilastnonzero != cd.len() - 1 {
            cd.truncate(ilastnonzero);
        }
        if cd.len() < 3 {
            // Degree < 2 after trimming: no interior extrema representable
            // (upstream would fault on the companion construction here).
            continue;
        }
        let mut a = companion_matrix_transposed(&cd);
        balance_matrix(&mut a);
        let Some(eigs) = hqr(&a) else {
            continue; // non-converged eigensolve: treat as no extrema (loudly rare)
        };
        for (re_n11, im) in eigs {
            if im.abs() < thresh_im && (-1.0..=1.0).contains(&re_n11) {
                x_at_extrema
                    .push(((expan.xmax - expan.xmin) * re_n11 + (expan.xmax + expan.xmin)) / 2.0);
            }
        }
    }
    x_at_extrema.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x_at_extrema
}

/// Upstream `build_monotonic_intervals`.
fn build_monotonic_intervals(
    expansions: &[OwnedInterval],
    x_at_extrema: &[f64],
) -> Vec<IntervalMatch> {
    let sort2 = |x: &mut f64, y: &mut f64| {
        if *x > *y {
            std::mem::swap(x, y);
        }
    };
    let mut intervals = Vec::new();
    let mut newx = x_at_extrema.to_vec();
    newx.insert(0, expansions.first().unwrap().xmin);
    newx.push(expansions.last().unwrap().xmax);

    for j in 0..(newx.len() - 1) {
        let (xmin, xmax) = (newx[j], newx[j + 1]);
        let mut im = IntervalMatch {
            expansioninfo: Vec::new(),
            xmin,
            xmax,
            ymin: 0.0,
            ymax: 0.0,
        };
        for (i, e) in expansions.iter().enumerate() {
            let a = e.xmin.max(xmin);
            let b = e.xmax.min(xmax);
            if a < b {
                let mut ymin = clenshaw(e, a);
                let mut ymax = clenshaw(e, b);
                sort2(&mut ymin, &mut ymax);
                im.expansioninfo.push(MonotonicExpansionMatch {
                    idx: i,
                    ymin,
                    ymax,
                    xmin: a,
                    xmax: b,
                });
            }
        }
        let mut yminoverall = eval_piecewise(expansions, xmin);
        let mut ymaxoverall = eval_piecewise(expansions, xmax);
        sort2(&mut yminoverall, &mut ymaxoverall);
        im.ymin = yminoverall;
        im.ymax = ymaxoverall;
        intervals.push(im);
    }
    intervals
}

#[cfg(test)]
mod tests {
    use super::*;

    /// T_5's interior extrema are exactly cos(k*pi/5), k = 1..=4.
    #[test]
    fn extrema_of_chebyshev_t5() {
        let mut coef = vec![0.0; 6];
        coef[5] = 1.0;
        let exp = OwnedInterval {
            xmin: -1.0,
            xmax: 1.0,
            coef,
        };
        let appr = ChebApprox1d::new(vec![exp]);
        assert_eq!(appr.x_at_extrema.len(), 4, "{:?}", appr.x_at_extrema);
        let mut expected: Vec<f64> = (1..=4)
            .map(|k| ((k as f64) * std::f64::consts::PI / 5.0).cos())
            .collect();
        expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (got, want) in appr.x_at_extrema.iter().zip(&expected) {
            assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
        }
    }

    /// Every reported extremum of a generic expansion must zero the
    /// derivative, and a fine sign-scan of y' must find the same count.
    #[test]
    fn extrema_match_derivative_sign_scan() {
        // A wiggly degree-12 expansion on [2, 7] (fixed pseudo-random coefficients).
        let coef = vec![
            0.7, -1.3, 0.41, 0.9, -0.27, 0.33, -0.12, 0.08, -0.05, 0.03, -0.011, 0.004, -0.0017,
        ];
        let exp = OwnedInterval {
            xmin: 2.0,
            xmax: 7.0,
            coef,
        };
        let dc = deriv_coeffs(&exp.coef, exp.xmin, exp.xmax);
        let dexp = OwnedInterval {
            xmin: 2.0,
            xmax: 7.0,
            coef: dc,
        };
        let appr = ChebApprox1d::new(vec![exp]);
        // Sign-scan the derivative on a fine grid
        let n = 20000;
        let mut scan_count = 0;
        let mut prev = clenshaw(&dexp, 2.0);
        for i in 1..=n {
            let x = 2.0 + 5.0 * (i as f64) / (n as f64);
            let v = clenshaw(&dexp, x);
            if prev != 0.0 && v != 0.0 && (prev < 0.0) != (v < 0.0) {
                scan_count += 1;
            }
            prev = v;
        }
        assert_eq!(appr.x_at_extrema.len(), scan_count);
        for &x in &appr.x_at_extrema {
            let d = clenshaw(&dexp, x);
            assert!(d.abs() < 1e-9, "y'({x}) = {d}");
        }
    }

    /// get_x_for_y on a non-monotonic curve finds every crossing.
    #[test]
    fn get_x_for_y_finds_all_roots() {
        // y = T_3(x) on [-1,1]: y = 0.5 has 3 solutions.
        let mut coef = vec![0.0; 4];
        coef[3] = 1.0;
        let appr = ChebApprox1d::new(vec![OwnedInterval {
            xmin: -1.0,
            xmax: 1.0,
            coef,
        }]);
        let mut solns = appr.get_x_for_y(0.5, 48, 100, 1e-14);
        solns.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(solns.len(), 3, "{solns:?}");
        for &x in &solns {
            let y = 4.0 * x * x * x - 3.0 * x;
            assert!((y - 0.5).abs() < 1e-10, "T3({x}) = {y}");
        }
    }
}
