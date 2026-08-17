//! Scalar rootfinders ported operation-for-operation from upstream
//! `src/Solvers.cpp` and `src/CPnumerics.cpp` @ v8.0.0 (Halley, Householder4,
//! Brent, solve_cubic), with upstream defaults (omega = 1, xtol_rel = 1e-12).
//! Error conditions map as everywhere else: invalid numbers -> `Error::Value`,
//! iteration caps -> `Error::Solution`.

use rustprop_core::{Error, Result};

/// Residual with derivatives (upstream `FuncWrapper1DWithThreeDerivs`).
/// `call` may mutate cached state (upstream updates the backend state).
pub(crate) trait Resid1D {
    fn call(&mut self, x: f64) -> f64;
    fn deriv(&mut self, x: f64) -> f64;
    fn second_deriv(&mut self, x: f64) -> f64;
    fn third_deriv(&mut self, x: f64) -> f64;
}

/// Upstream `Halley` (omega = 1).
pub(crate) fn halley<R: Resid1D>(f: &mut R, x0: f64, ftol: f64, maxiter: i32) -> Result<f64> {
    halley_omega(f, x0, ftol, maxiter, 1.0)
}

/// Upstream `Halley` with the `omega` damping option (the spinodal finder
/// retries with omega = 0.7, 0.5, 0.3, 0.1).
pub(crate) fn halley_omega<R: Resid1D>(
    f: &mut R,
    x0: f64,
    ftol: f64,
    maxiter: i32,
    omega: f64,
) -> Result<f64> {
    let xtol_rel = 1e-12;
    let mut iter = 0;
    let mut x = x0;
    let mut fval: f64 = 999.0;
    while iter < 2 || fval.abs() > ftol {
        fval = f.call(x);
        let dfdx = f.deriv(x);
        let d2fdx2 = f.second_deriv(x);
        if !fval.is_finite() {
            return Err(Error::Value(
                "Residual function in Halley returned invalid number".into(),
            ));
        }
        if !dfdx.is_finite() {
            return Err(Error::Value(
                "Derivative function in Halley returned invalid number".into(),
            ));
        }
        let dx = -(2.0 * fval * dfdx) / (2.0 * dfdx * dfdx - fval * d2fdx2);
        x += omega * dx;
        if (dx / x).abs() < xtol_rel {
            return Ok(x);
        }
        if iter > maxiter {
            return Err(Error::Solution(
                "Halley reached maximum number of iterations".into(),
            ));
        }
        iter += 1;
    }
    Ok(x)
}

/// Upstream `Householder4` (omega = 1).
pub(crate) fn householder4<R: Resid1D>(f: &mut R, x0: f64, ftol: f64, maxiter: i32) -> Result<f64> {
    let xtol_rel = 1e-12;
    let mut iter = 1;
    let mut x = x0;
    let mut fval: f64 = 999.0;
    while iter < 2 || fval.abs() > ftol {
        fval = f.call(x);
        let dfdx = f.deriv(x);
        let d2fdx2 = f.second_deriv(x);
        let d3fdx3 = f.third_deriv(x);
        if !fval.is_finite() {
            return Err(Error::Value(
                "Residual function in Householder4 returned invalid number".into(),
            ));
        }
        if !dfdx.is_finite() {
            return Err(Error::Value(
                "Derivative function in Householder4 returned invalid number".into(),
            ));
        }
        if !d2fdx2.is_finite() {
            return Err(Error::Value(
                "Second derivative function in Householder4 returned invalid number".into(),
            ));
        }
        if !d3fdx3.is_finite() {
            return Err(Error::Value(
                "Third derivative function in Householder4 returned invalid number".into(),
            ));
        }
        let dx = -fval * (dfdx * dfdx - fval * d2fdx2 / 2.0)
            / (dfdx * dfdx * dfdx - fval * dfdx * d2fdx2 + d3fdx3 * fval * fval / 6.0);
        x += dx;
        if (dx / x).abs() < xtol_rel {
            return Ok(x);
        }
        if iter > maxiter {
            return Err(Error::Solution(
                "Householder4 reached maximum number of iterations".into(),
            ));
        }
        iter += 1;
    }
    Ok(x)
}

/// Upstream `BoundedSecant`: secant with hard bounds — a step past a bound
/// goes half the way to it instead.
pub fn bounded_secant<F: FnMut(f64) -> f64>(
    mut call: F,
    x0: f64,
    xmin: f64,
    xmax: f64,
    dx: f64,
    tol: f64,
    maxiter: i32,
) -> Result<f64> {
    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut y1 = 0.0;
    let mut x;
    let mut fval: f64 = 999.0;
    let mut iter = 1;
    if dx.abs() == 0.0 {
        return Err(Error::Value("dx cannot be zero".into()));
    }
    while iter <= 3 || fval.abs() > tol {
        if iter == 1 {
            x1 = x0;
            x = x1;
        } else if iter == 2 {
            x2 = x0 + dx;
            x = x2;
        } else {
            x = x2;
        }
        fval = call(x);
        if iter == 1 {
            y1 = fval;
        } else {
            let y2 = fval;
            let mut x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
            // Check bounds, go half the way to the limit if limit is exceeded
            if x3 < xmin {
                x3 = (xmin + x2) / 2.0;
            }
            if x3 > xmax {
                x3 = (xmax + x2) / 2.0;
            }
            y1 = y2;
            x1 = x2;
            x2 = x3;
        }
        if iter > maxiter {
            return Err(Error::Solution(format!(
                "BoundedSecant reached maximum number of iterations of {maxiter}"
            )));
        }
        iter += 1;
    }
    Ok(x2)
}

/// Upstream `ExtrapolatingSecant` (`src/Solvers.cpp:429-506`, omega = 1).
///
/// Identical to [`secant`] except in what an INVALID residual means: `Secant`
/// throws ("Residual function in secant returned invalid number",
/// Solvers.cpp:332), while `ExtrapolatingSecant` EXTRAPOLATES through it
/// (Solvers.cpp:472-478) — returning `x0` if the very first call is invalid and
/// otherwise the secant step built from the last two VALID residuals,
/// `x2 - omega*y1/(y1 - y0)*(x2 - x1)`. That is the whole point of the routine
/// and its only upstream caller, `SaturationAncillaryFunction::invert`
/// (`Ancillaries.cpp:111`), relies on it: the pseudo-pure pL/pV curves get
/// inverted at pressures ABOVE their value at `Tmax`, where `evaluate`
/// deliberately returns NaN past the reducing temperature (#1611), and the
/// extrapolated temperature is what the phase determination then uses.
pub fn extrapolating_secant<F: FnMut(f64) -> f64>(
    mut call: F,
    x0: f64,
    dx: f64,
    tol: f64,
    maxiter: i32,
) -> Result<f64> {
    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut x3 = 0.0;
    let mut y0 = 0.0;
    let mut y1 = 0.0;
    let mut x;
    let mut fval: f64 = 999.0;
    let mut iter = 1;
    if dx.abs() == 0.0 {
        return Err(Error::Value("dx cannot be zero".into()));
    }
    while iter <= 2 || fval.abs() > tol {
        if iter == 1 {
            x1 = x0;
            x = x1;
        } else if iter == 2 {
            x2 = x0 + dx;
            x = x2;
        } else {
            x = x2;
        }
        fval = call(x);
        if !fval.is_finite() {
            if iter == 1 {
                return Ok(x);
            }
            return Ok(x2 - y1 / (y1 - y0) * (x2 - x1));
        }
        if iter == 1 {
            y1 = fval;
        }
        if iter > 1 {
            let deltax = x2 - x1;
            if deltax.abs() < 1e-14 {
                return Ok(x);
            }
            let y2 = fval;
            let deltay = y2 - y1;
            if iter > 2 && deltay.abs() < 1e-14 {
                return Ok(x);
            }
            x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
            y0 = y1;
            y1 = y2;
            x1 = x2;
            x2 = x3;
        }
        if iter > maxiter {
            return Err(Error::Solution(
                "Secant reached maximum number of iterations".into(),
            ));
        }
        iter += 1;
    }
    Ok(x3)
}

/// Upstream `Secant` (omega = 1; the `input_not_in_range` hook is unused by
/// the ported callers).
pub fn secant<F: FnMut(f64) -> f64>(
    mut call: F,
    x0: f64,
    dx: f64,
    tol: f64,
    maxiter: i32,
) -> Result<f64> {
    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut x3 = 0.0;
    let mut y1 = 0.0;
    let mut x;
    let mut fval: f64 = 999.0;
    let mut iter = 1;
    if dx.abs() == 0.0 {
        return Err(Error::Value("dx cannot be zero".into()));
    }
    while iter <= 2 || fval.abs() > tol {
        if iter == 1 {
            x1 = x0;
            x = x1;
        } else if iter == 2 {
            x2 = x0 + dx;
            x = x2;
        } else {
            x = x2;
        }
        fval = call(x);
        if !fval.is_finite() {
            return Err(Error::Value(
                "Residual function in secant returned invalid number".into(),
            ));
        }
        if iter == 1 {
            y1 = fval;
        }
        if iter > 1 {
            let deltax = x2 - x1;
            if deltax.abs() < 1e-14 {
                return Ok(x);
            }
            let y2 = fval;
            let deltay = y2 - y1;
            if iter > 2 && deltay.abs() < 1e-14 {
                return Ok(x);
            }
            x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
            y1 = y2;
            x1 = x2;
            x2 = x3;
        }
        if iter > maxiter {
            return Err(Error::Solution(
                "Secant reached maximum number of iterations".into(),
            ));
        }
        iter += 1;
    }
    Ok(x3)
}

/// Upstream `Brent` (the ALGOL-derived variant with all its guards).
pub fn brent<F: FnMut(f64) -> f64>(
    mut call: F,
    mut a: f64,
    mut b: f64,
    macheps: f64,
    t: f64,
    maxiter: i32,
) -> Result<f64> {
    let mut fa = call(a);
    let mut fb = call(b);
    if fb.abs() < t {
        return Ok(b);
    }
    if !fb.is_finite() {
        return Err(Error::Value(format!(
            "Brent's method f(b) is NAN for b = {b}, other input was a = {a}"
        )));
    }
    if fa.abs() < t {
        return Ok(a);
    }
    if !fa.is_finite() {
        return Err(Error::Value(format!(
            "Brent's method f(a) is NAN for a = {a}, other input was b = {b}"
        )));
    }
    if fa * fb > 0.0 {
        return Err(Error::Value(format!(
            "Inputs in Brent [{a:.6e},{b:.6e}] do not bracket the root.  Function values are [{fa:.6e},{fb:.6e}]"
        )));
    }

    let mut c = a;
    let mut fc = fa;
    let mut iter = 1;
    if fc.abs() < fb.abs() {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
    }
    let mut d = b - a;
    let mut e = b - a;
    let mut m = 0.5 * (c - b);
    let mut tol = 2.0 * macheps * b.abs() + t;
    while m.abs() > tol && fb != 0.0 {
        if e.abs() < tol || fa.abs() <= fb.abs() {
            m = 0.5 * (c - b);
            d = m;
            e = m;
        } else {
            let mut p;
            let mut q;
            let mut s = fb / fa;
            if a == c {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                let r = fb / fc;
                m = 0.5 * (c - b);
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            s = e;
            e = d;
            m = 0.5 * (c - b);
            if 2.0 * p < 3.0 * m * q - (tol * q).abs() || p < (0.5 * s * q).abs() {
                d = p / q;
            } else {
                m = 0.5 * (c - b);
                d = m;
                e = m;
            }
        }
        a = b;
        fa = fb;
        if d.abs() > tol {
            b += d;
        } else if m > 0.0 {
            b += tol;
        } else {
            b += -tol;
        }
        fb = call(b);
        if !fb.is_finite() {
            return Err(Error::Value(format!(
                "Brent's method f(t) is NAN for t = {b}"
            )));
        }
        if fb.abs() < macheps {
            return Ok(b);
        }
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        m = 0.5 * (c - b);
        tol = 2.0 * macheps * b.abs() + t;
        iter += 1;
        if iter > maxiter {
            return Err(Error::Solution(format!(
                "Brent's method reached maximum number of steps of {maxiter} "
            )));
        }
        if fb.abs() < 2.0 * macheps * b.abs() {
            return Ok(b);
        }
    }
    Ok(b)
}

/// Upstream `solve_cubic` (CPnumerics.cpp): roots of a*x^3+b*x^2+c*x+d.
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> (i32, f64, f64, f64) {
    if a.abs() < 10.0 * f64::EPSILON {
        if b.abs() < 10.0 * f64::EPSILON {
            return (1, -d / c, f64::NAN, f64::NAN);
        }
        let x0 = (-c + (c * c - 4.0 * b * d).sqrt()) / (2.0 * b);
        let x1 = (-c - (c * c - 4.0 * b * d).sqrt()) / (2.0 * b);
        return (2, x0, x1, f64::NAN);
    }
    let delta = 18.0 * a * b * c * d - 4.0 * b * b * b * d + b * b * c * c
        - 4.0 * a * c * c * c
        - 27.0 * a * a * d * d;
    let p = (3.0 * a * c - b * b) / (3.0 * a * a);
    let q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
    if delta < 0.0 {
        // One real root
        let t0 = if 4.0 * p * p * p + 27.0 * q * q > 0.0 && p < 0.0 {
            -2.0 * q.abs() / q
                * (-p / 3.0).sqrt()
                * (1.0 / 3.0 * (-3.0 * q.abs() / (2.0 * p) * (-3.0 / p).sqrt()).acosh()).cosh()
        } else {
            -2.0 * (p / 3.0).sqrt()
                * (1.0 / 3.0 * (3.0 * q / (2.0 * p) * (3.0 / p).sqrt()).asinh()).sinh()
        };
        let x = t0 - b / (3.0 * a);
        (1, x, x, x)
    } else {
        // Three real roots
        let base = 3.0 * q / (2.0 * p) * (-3.0 / p).sqrt();
        let t0 = 2.0 * (-p / 3.0).sqrt() * (1.0 / 3.0 * base.acos()).cos();
        let t1 = 2.0
            * (-p / 3.0).sqrt()
            * (1.0 / 3.0 * base.acos() - 2.0 * std::f64::consts::PI / 3.0).cos();
        let t2 = 2.0
            * (-p / 3.0).sqrt()
            * (1.0 / 3.0 * base.acos() - 2.0 * 2.0 * std::f64::consts::PI / 3.0).cos();
        (
            3,
            t0 - b / (3.0 * a),
            t1 - b / (3.0 * a),
            t2 - b / (3.0 * a),
        )
    }
}

// ---------------------------------------------------------------------------
// TOMS Algorithm 748 (Boost boost/math/tools/toms748_solve.hpp, verbatim
// port). The mixture sweep flashes need the real interpolating behavior —
// plain bisection can step into pathological pockets of a discontinuous
// residual that the secant/cubic steps hug past (observed: the HSU_P sweep
// walking into the wrong-root pocket near a phase-boundary discontinuity
// that upstream's TOMS748 sails over).
// ---------------------------------------------------------------------------

fn t748_safe_div(num: f64, denom: f64, r: f64) -> f64 {
    if denom.abs() < 1.0 && (denom * f64::MAX).abs() <= num.abs() {
        return r;
    }
    num / denom
}

fn t748_secant_interpolate(a: f64, b: f64, fa: f64, fb: f64) -> f64 {
    let tol = f64::EPSILON * 5.0;
    let c = a - (fa / (fb - fa)) * (b - a);
    if c <= a + a.abs() * tol || c >= b - b.abs() * tol {
        return (a + b) / 2.0;
    }
    c
}

#[allow(clippy::many_single_char_names)]
fn t748_quadratic_interpolate(
    a: f64,
    b: f64,
    d: f64,
    fa: f64,
    fb: f64,
    fd: f64,
    count: u32,
) -> f64 {
    let big_b = t748_safe_div(fb - fa, b - a, f64::MAX);
    let mut big_a = t748_safe_div(fd - fb, d - b, f64::MAX);
    big_a = t748_safe_div(big_a - big_b, d - a, 0.0);

    if big_a == 0.0 {
        return t748_secant_interpolate(a, b, fa, fb);
    }
    let mut c = if big_a.signum() * fa.signum() > 0.0 {
        a
    } else {
        b
    };
    for _ in 1..=count {
        c -= t748_safe_div(
            fa + (big_b + big_a * (c - b)) * (c - a),
            big_b + big_a * (2.0 * c - a - b),
            1.0 + c - a,
        );
    }
    if c <= a || c >= b {
        c = t748_secant_interpolate(a, b, fa, fb);
    }
    c
}

#[allow(clippy::too_many_arguments)]
fn t748_cubic_interpolate(
    a: f64,
    b: f64,
    d: f64,
    e: f64,
    fa: f64,
    fb: f64,
    fd: f64,
    fe: f64,
) -> f64 {
    let q11 = (d - e) * fd / (fe - fd);
    let q21 = (b - d) * fb / (fd - fb);
    let q31 = (a - b) * fa / (fb - fa);
    let d21 = (b - d) * fd / (fd - fb);
    let d31 = (a - b) * fb / (fb - fa);
    let q22 = (d21 - q11) * fb / (fe - fb);
    let q32 = (d31 - q21) * fa / (fd - fa);
    let d32 = (d31 - q21) * fd / (fd - fa);
    let q33 = (d32 - q22) * fa / (fe - fa);
    let c = q31 + q32 + q33 + a;
    if c <= a || c >= b {
        return t748_quadratic_interpolate(a, b, d, fa, fb, fd, 3);
    }
    c
}

/// Boost `detail::bracket`: evaluate at (adjusted) c and shrink [a,b],
/// tracking the ejected third-best point in (d, fd).
#[allow(clippy::too_many_arguments)]
fn t748_bracket<F: FnMut(f64) -> Result<f64>>(
    f: &mut F,
    a: &mut f64,
    b: &mut f64,
    mut c: f64,
    fa: &mut f64,
    fb: &mut f64,
    d: &mut f64,
    fd: &mut f64,
) -> Result<()> {
    let tol = f64::EPSILON * 2.0;
    if (*b - *a) < 2.0 * tol * *a {
        c = *a + (*b - *a) / 2.0;
    } else if c <= *a + a.abs() * tol {
        c = *a + a.abs() * tol;
    } else if c >= *b - b.abs() * tol {
        c = *b - b.abs() * tol;
    }
    let fc = f(c)?;
    if fc == 0.0 {
        *a = c;
        *fa = 0.0;
        *d = 0.0;
        *fd = 0.0;
        return Ok(());
    }
    if fa.signum() * fc.signum() < 0.0 {
        *d = *b;
        *fd = *fb;
        *b = c;
        *fb = fc;
    } else {
        *d = *a;
        *fd = *fa;
        *a = c;
        *fa = fc;
    }
    Ok(())
}

/// `boost::math::tools::toms748_solve` with `eps_tolerance<double>(bits)`;
/// returns the midpoint of the final bracket (as every upstream call site
/// takes `(bracket.first + bracket.second) / 2`).
pub(crate) fn toms748_solve<F: FnMut(f64) -> Result<f64>>(
    f: &mut F,
    ax: f64,
    bx: f64,
    fax: f64,
    fbx: f64,
    bits: u32,
    max_iter: u32,
) -> Result<f64> {
    let eps = (2.0_f64.powi(1 - (bits as i32))).max(4.0 * f64::EPSILON);
    let tol = |a: f64, b: f64| (a - b).abs() <= eps * a.abs().min(b.abs());
    let mu = 0.5;

    let mut count = max_iter;
    let (mut a, mut b, mut fa, mut fb) = (ax, bx, fax, fbx);
    if a >= b {
        return Err(Error::Value("Parameters a and b out of order".into()));
    }
    if tol(a, b) || fa == 0.0 || fb == 0.0 {
        if fa == 0.0 {
            b = a;
        } else if fb == 0.0 {
            a = b;
        }
        return Ok((a + b) / 2.0);
    }
    if fa.signum() * fb.signum() > 0.0 {
        return Err(Error::Value(
            "Parameters a and b do not bracket the root".into(),
        ));
    }
    let (mut d, mut fd, mut e, mut fe) = (0.0_f64, 1e5_f64, 1e5_f64, 1e5_f64);
    let mut c;

    if fa != 0.0 {
        c = t748_secant_interpolate(a, b, fa, fb);
        t748_bracket(f, &mut a, &mut b, c, &mut fa, &mut fb, &mut d, &mut fd)?;
        count -= 1;
        if count > 0 && fa != 0.0 && !tol(a, b) {
            c = t748_quadratic_interpolate(a, b, d, fa, fb, fd, 2);
            e = d;
            fe = fd;
            t748_bracket(f, &mut a, &mut b, c, &mut fa, &mut fb, &mut d, &mut fd)?;
            count -= 1;
        }
    }

    while count > 0 && fa != 0.0 && !tol(a, b) {
        let a0 = a;
        let b0 = b;
        let min_diff = f64::MIN_POSITIVE * 32.0;
        let prof = (fa - fb).abs() < min_diff
            || (fa - fd).abs() < min_diff
            || (fa - fe).abs() < min_diff
            || (fb - fd).abs() < min_diff
            || (fb - fe).abs() < min_diff
            || (fd - fe).abs() < min_diff;
        c = if prof {
            t748_quadratic_interpolate(a, b, d, fa, fb, fd, 2)
        } else {
            t748_cubic_interpolate(a, b, d, e, fa, fb, fd, fe)
        };
        e = d;
        fe = fd;
        t748_bracket(f, &mut a, &mut b, c, &mut fa, &mut fb, &mut d, &mut fd)?;
        count -= 1;
        if count == 0 || fa == 0.0 || tol(a, b) {
            break;
        }

        let prof = (fa - fb).abs() < min_diff
            || (fa - fd).abs() < min_diff
            || (fa - fe).abs() < min_diff
            || (fb - fd).abs() < min_diff
            || (fb - fe).abs() < min_diff
            || (fd - fe).abs() < min_diff;
        c = if prof {
            t748_quadratic_interpolate(a, b, d, fa, fb, fd, 3)
        } else {
            t748_cubic_interpolate(a, b, d, e, fa, fb, fd, fe)
        };
        t748_bracket(f, &mut a, &mut b, c, &mut fa, &mut fb, &mut d, &mut fd)?;
        count -= 1;
        if count == 0 || fa == 0.0 || tol(a, b) {
            break;
        }

        // Double-length secant step
        let (u, fu) = if fa.abs() < fb.abs() {
            (a, fa)
        } else {
            (b, fb)
        };
        c = u - 2.0 * (fu / (fb - fa)) * (b - a);
        if (c - u).abs() > (b - a) / 2.0 {
            c = a + (b - a) / 2.0;
        }
        e = d;
        fe = fd;
        t748_bracket(f, &mut a, &mut b, c, &mut fa, &mut fb, &mut d, &mut fd)?;
        count -= 1;
        if count == 0 || fa == 0.0 || tol(a, b) {
            break;
        }

        if (b - a) < mu * (b0 - a0) {
            continue;
        }
        // Not converging fast enough: bisection
        e = d;
        fe = fd;
        let mid = a + (b - a) / 2.0;
        t748_bracket(f, &mut a, &mut b, mid, &mut fa, &mut fb, &mut d, &mut fd)?;
        count -= 1;
    }

    if fa == 0.0 {
        b = a;
    } else if fb == 0.0 {
        a = b;
    }
    Ok((a + b) / 2.0)
}
