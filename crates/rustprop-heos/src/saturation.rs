//! Pure-fluid saturation flashes (PLAN.md 4.4) — the v8 superancillary path
//! of `FlashRoutines::QT_flash` / `PQ_flash`: for a pure fluid with a
//! superancillary (ENABLE_SUPERANCILLARIES defaults on), saturation states
//! come from direct superancillary evaluation — the classic Maxwell iteration
//! is only the fallback for fluids without one and is ported when such a
//! fluid arrives.
//!
//! The PQ inverse goes through upstream's fitted `T(ln p)` Chebyshev
//! approximation (`make_invlnp`): Lobatto-node evaluation of `T` from
//! rootfinding on the monotonic p(T) curve, cosine-matrix conversion of node
//! values to coefficients, and dyadic splitting at tol 1e-12 (upstream's
//! parameters, degree = that of the p expansions). Two scoped deviations,
//! logged in PLAN.md:
//! - the saturation p(T) curve is verified strictly monotonic at construction
//!   and the inverse uses per-expansion endpoint brackets — exactly what
//!   upstream's extrema machinery degenerates to when there are no interior
//!   extrema (a non-monotonic curve fails loudly here);
//! - upstream brackets with boost TOMS748 to 64-bit eps tolerance; we bisect
//!   to machine tightness and return the bracket midpoint — the same double
//!   to ~1 ULP.

use crate::superancillary::{Cheb1D, OwnedInterval, eval_piecewise};
use rustprop_core::fluid::SuperAncillaryData;
use rustprop_core::{Error, Result};

/// A saturation state produced by a QT or PQ flash.
#[derive(Debug, Clone, Copy)]
pub struct SatState {
    pub t: f64,
    pub p: f64,
    /// Saturated liquid molar density [mol/m^3]
    pub rho_l: f64,
    /// Saturated vapor molar density [mol/m^3]
    pub rho_v: f64,
    /// Two-phase mixture molar density at the requested quality [mol/m^3]
    pub rhomolar: f64,
    pub q: f64,
}

/// Runtime saturation evaluator for one fluid: the static superancillary data
/// plus the fitted `T(ln p)` inverse (upstream builds it lazily; we build at
/// construction — same fit, no observable difference).
pub struct SaturationSuperAncillary {
    data: &'static SuperAncillaryData,
    invlnp: Vec<OwnedInterval>,
    pmin: f64,
    pmax: f64,
}

impl SaturationSuperAncillary {
    pub fn new(data: &'static SuperAncillaryData) -> Self {
        // Scope guard: the inverse below assumes p(T) strictly increasing.
        // Adjacent expansions evaluate the shared breakpoint a few ULPs
        // apart, so the cross-interval comparison carries ULP-scale slack.
        let mut last = f64::NEG_INFINITY;
        for e in data.p {
            let (ymin, ymax) = (eval_at(e, e.xmin), eval_at(e, e.xmax));
            assert!(
                ymax > ymin && ymin >= last - last.abs() * 1e-12,
                "superancillary p(T) not monotonically increasing; port the upstream extrema machinery"
            );
            last = ymax;
        }
        let first = &data.p[0];
        let last_e = &data.p[data.p.len() - 1];
        let pmin = eval_at(first, first.xmin);
        let pmax = eval_at(last_e, last_e.xmax);
        let invlnp = make_invlnp(data, pmin, pmax);
        SaturationSuperAncillary {
            data,
            invlnp,
            pmin,
            pmax,
        }
    }

    /// Upstream `SuperAncillary::get_pmax` (numerical maximum of the p curve).
    pub fn pmax(&self) -> f64 {
        self.pmax
    }
    pub fn pmin(&self) -> f64 {
        self.pmin
    }

    /// Superancillary branch of upstream `FlashRoutines::QT_flash`.
    pub fn qt_flash(&self, t: f64, q: f64) -> Result<SatState> {
        let tcrit_num = self.data.t_crit_num;
        if t > tcrit_num {
            return Err(Error::Value(format!(
                "Temperature to QT_flash [{t:.8} K] may not be above the numerical critical point of {tcrit_num:.15} K"
            )));
        }
        let rho_l = eval_piecewise(self.data.rho_l, t);
        let rho_v = eval_piecewise(self.data.rho_v, t);
        let p = eval_piecewise(self.data.p, t);
        let rhomolar = 1.0 / (q / rho_v + (1.0 - q) / rho_l);
        Ok(SatState {
            t,
            p,
            rho_l,
            rho_v,
            rhomolar,
            q,
        })
    }

    /// Superancillary branch of upstream `FlashRoutines::PQ_flash`.
    pub fn pq_flash(&self, p: f64, q: f64) -> Result<SatState> {
        if p > self.pmax {
            return Err(Error::Value(format!(
                "Pressure to PQ_flash [{p:.8} Pa] may not be above the numerical critical point of {:.15} Pa",
                self.pmax
            )));
        }
        let t = eval_piecewise(&self.invlnp, p.ln());
        let rho_l = eval_piecewise(self.data.rho_l, t);
        let rho_v = eval_piecewise(self.data.rho_v, t);
        let rhomolar = 1.0 / (q / rho_v + (1.0 - q) / rho_l);
        Ok(SatState {
            t,
            p,
            rho_l,
            rho_v,
            rhomolar,
            q,
        })
    }
}

fn eval_at<C: Cheb1D>(e: &C, x: f64) -> f64 {
    crate::superancillary::clenshaw(e, x)
}

// ---------------------------------------------------------------------------
// Inverse machinery (upstream `get_x_for_y` on the monotonic p curve and
// `make_invlnp`)
// ---------------------------------------------------------------------------

/// Bracketed root of `f` in [a, b] (f(a), f(b) supplied). Upstream uses boost
/// TOMS748 with 64-bit eps tolerance and returns the bracket midpoint;
/// bisection to machine tightness yields the same double to ~1 ULP.
fn bracket_solve<F: Fn(f64) -> f64>(f: F, mut a: f64, mut b: f64, fa_in: f64, fb_in: f64) -> f64 {
    let (mut fa, mut fb) = (fa_in, fb_in);
    debug_assert!(fa * fb <= 0.0);
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if m <= a || m >= b {
            break; // bracket at machine tightness
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
    let _ = (fa, fb);
    0.5 * (a + b)
}

/// Upstream `ChebyshevExpansion::solve_for_x_count` bounds handling +
/// `ChebyshevApproximation1D::get_x_for_y` restricted to a monotonic curve:
/// each expansion is monotonic, so an expansion whose endpoint values bracket
/// `y` contains exactly one solution.
fn get_x_for_y(
    expansions: &[rustprop_core::fluid::ChebyshevInterval],
    y: f64,
    boundsftol: f64,
) -> Vec<f64> {
    let mut solns = Vec::new();
    for e in expansions {
        let (ymin, ymax) = (eval_at(e, e.xmin), eval_at(e, e.xmax));
        let (lo, hi) = if ymin <= ymax {
            (ymin, ymax)
        } else {
            (ymax, ymin)
        };
        if y >= lo && y <= hi {
            let f = |x: f64| eval_at(e, x) - y;
            let fa = f(e.xmin);
            if fa.abs() < boundsftol {
                solns.push(e.xmin);
                continue;
            }
            let fb = f(e.xmax);
            if fb.abs() < boundsftol {
                solns.push(e.xmax);
                continue;
            }
            solns.push(bracket_solve(f, e.xmin, e.xmax, fa, fb));
        }
    }
    solns
}

/// Upstream `detail::get_LU_matrices` — only the L matrix (values -> coeffs).
#[allow(clippy::needless_range_loop)] // symmetric fill writes l[j][k] and l[k][j]
fn l_matrix(n: usize) -> Vec<Vec<f64>> {
    let nf = n as f64;
    let mut l = vec![vec![0.0; n + 1]; n + 1];
    for j in 0..=n {
        for k in j..=n {
            let p_j: f64 = if j == 0 || j == n { 2.0 } else { 1.0 };
            let p_k: f64 = if k == 0 || k == n { 2.0 } else { 1.0 };
            let cosjpikn = ((j as f64) * std::f64::consts::PI * (k as f64) / nf).cos();
            l[j][k] = 2.0 / (p_j * p_k * nf) * cosjpikn;
            l[k][j] = l[j][k];
        }
    }
    l
}

/// Upstream `detail::M_element_norm`: |tail(M)| / |head(M)|.
fn m_element_norm(c: &[f64], m: usize) -> f64 {
    let tail: f64 = c[c.len() - m..].iter().map(|v| v * v).sum::<f64>().sqrt();
    let head: f64 = c[..m].iter().map(|v| v * v).sum::<f64>().sqrt();
    tail / head
}

/// Upstream `SuperAncillary::make_invlnp` + `detail::dyadic_splitting` with
/// upstream's parameters (M = 3, tol = 1e-12, max 26 passes, degree = that of
/// the p expansions).
fn make_invlnp(data: &SuperAncillaryData, pmin: f64, pmax: f64) -> Vec<OwnedInterval> {
    let ndegree = data.p[0].coef.len() - 1;
    let eps = f64::EPSILON;
    let l = l_matrix(ndegree);

    // The node function: T for a given ln(p), by inverting the p curve.
    let func = |i: usize, npts: usize, lnp: f64| -> f64 {
        let p = lnp.exp();
        let solns = get_x_for_y(data.p, p, 1e-8);
        if solns.len() != 1 {
            let edge = i == 0 || i == npts - 1;
            if edge && p > pmin * (1.0 - eps * 1000.0) && p < pmin * (1.0 + eps * 1000.0) {
                return data.p[0].xmin;
            }
            if edge && p > pmax * (1.0 - eps * 1000.0) && p < pmax * (1.0 + eps * 1000.0) {
                return data.p[data.p.len() - 1].xmax;
            }
            panic!(
                "Must have one solution for T(p), got {} for {p} Pa; limits are [{pmin} Pa, {pmax} Pa]; i is {i}",
                solns.len()
            );
        }
        solns[0]
    };

    let builder = |xmin: f64, xmax: f64| -> OwnedInterval {
        // Chebyshev-Lobatto nodes cos(pi*j/N) mapped to [xmin, xmax],
        // j = 0..=N (upstream LinSpaced ordering).
        let n = ndegree;
        let mut fvals = vec![0.0; n + 1];
        for (j, fv) in fvals.iter_mut().enumerate() {
            let node = ((j as f64) * std::f64::consts::PI / (n as f64)).cos();
            let x = ((xmax - xmin) * node + (xmax + xmin)) * 0.5;
            *fv = func(j, n + 1, x);
        }
        let mut coef = vec![0.0; n + 1];
        for (j, c) in coef.iter_mut().enumerate() {
            *c = l[j].iter().zip(&fvals).map(|(lv, fv)| lv * fv).sum();
        }
        assert!(
            coef.iter().all(|c| c.is_finite()),
            "At least one coefficient is non-finite"
        );
        OwnedInterval { xmin, xmax, coef }
    };

    let (m, tol, max_refine_passes) = (3usize, 1e-12, 26);
    let mut expansions = vec![builder(pmin.ln(), pmax.ln())];
    for _pass in 0..max_refine_passes {
        let mut all_converged = true;
        let mut i = expansions.len() as i64 - 1;
        while i >= 0 {
            let idx = i as usize;
            let err = m_element_norm(&expansions[idx].coef, m);
            if err > tol {
                let (xmin, xmax) = (expansions[idx].xmin, expansions[idx].xmax);
                let xmid = (xmin + xmax) / 2.0;
                let newleft = builder(xmin, xmid);
                let newright = builder(xmid, xmax);
                expansions[idx] = newleft;
                expansions.insert(idx + 1, newright);
                all_converged = false;
            }
            i -= 1;
        }
        if all_converged {
            break;
        }
    }
    expansions
}
