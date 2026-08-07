//! Port of upstream `include/CoolProp/svd/Hermite1D.h`.
//!
//! Cubic-Hermite interpolation on one interval. Given values `y0`, `y1` and
//! slopes `m0`, `m1` at the ends of `[x0, x1]` with `h = x1 - x0` and
//! `t = (x - x0)/h`:
//!
//! ```text
//! p(t) = h00(t) y0 + h*h10(t) m0 + h01(t) y1 + h*h11(t) m1
//! ```
//!
//! Slopes are an INPUT — this module has no opinion on where they came from.

/// The four Hermite basis polynomials at `t`.
#[derive(Clone, Copy, Debug, Default)]
pub struct HermiteBasis {
    pub h00: f64,
    pub h10: f64,
    pub h01: f64,
    pub h11: f64,
}

#[inline]
pub fn hermite_basis(t: f64) -> HermiteBasis {
    let t2 = t * t;
    let t3 = t2 * t;
    HermiteBasis {
        h00: 2.0 * t3 - 3.0 * t2 + 1.0,
        h10: t3 - 2.0 * t2 + t,
        h01: -2.0 * t3 + 3.0 * t2,
        h11: t3 - t2,
    }
}

#[inline]
pub fn hermite_eval(y0: f64, y1: f64, m0: f64, m1: f64, h: f64, t: f64) -> f64 {
    let b = hermite_basis(t);
    b.h00 * y0 + b.h10 * h * m0 + b.h01 * y1 + b.h11 * h * m1
}

/// `dp/dx` at parameter `t` (note `dt/dx = 1/h`, hence the division).
#[inline]
pub fn hermite_eval_deriv(y0: f64, y1: f64, m0: f64, m1: f64, h: f64, t: f64) -> f64 {
    let t2 = t * t;
    let dh00 = 6.0 * t2 - 6.0 * t;
    let dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    let dh01 = -6.0 * t2 + 6.0 * t;
    let dh11 = 3.0 * t2 - 2.0 * t;
    let dp_dt = dh00 * y0 + dh10 * h * m0 + dh01 * y1 + dh11 * h * m1;
    dp_dt / h
}
