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
        x += dx;
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

/// Upstream `Brent` (the ALGOL-derived variant with all its guards).
pub(crate) fn brent<F: FnMut(f64) -> f64>(
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
pub(crate) fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> (i32, f64, f64, f64) {
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
