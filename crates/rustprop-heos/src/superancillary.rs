//! v8 Chebyshev superancillary evaluation (PLAN.md 4.3) — port of the
//! evaluation path of `include/CoolProp/superancillary/superancillary.h`:
//! Clenshaw evaluation of one interval, Knuth-midpoint binary search across
//! the piecewise intervals, and `eval_sat`. The inverse (rootfinding) path
//! arrives with the saturation solver (4.4).

use rustprop_core::fluid::{ChebyshevInterval, SuperAncillaryData};

/// Upstream `ChebyshevExpansion::eval` — Clenshaw's method.
fn clenshaw(e: &ChebyshevInterval, x: f64) -> f64 {
    let xscaled = (2.0 * x - (e.xmax + e.xmin)) / (e.xmax - e.xmin);
    let norder = e.coef.len() - 1;
    let mut u_kp1 = e.coef[norder];
    let mut u_kp2 = 0.0;
    for k in (1..norder).rev() {
        let u_k = 2.0 * xscaled * u_kp1 - u_kp2 + e.coef[k];
        u_kp2 = u_kp1;
        u_kp1 = u_k;
    }
    e.coef[0] + xscaled * u_kp1 - u_kp2
}

/// Upstream `ChebyshevApproximation1D::get_index` (Knuth midpoint).
fn get_index(expansions: &[ChebyshevInterval], x: f64) -> usize {
    let midpoint_knuth = |a: i32, b: i32| (a & b) + ((a ^ b) >> 1);
    let mut il = 0i32;
    let mut ir = expansions.len() as i32 - 1;
    while ir - il > 1 {
        let im = midpoint_knuth(il, ir);
        if x >= expansions[im as usize].xmin {
            il = im;
        } else {
            ir = im;
        }
    }
    if x < expansions[il as usize].xmax {
        il as usize
    } else {
        ir as usize
    }
}

/// Upstream `ChebyshevApproximation1D::eval`.
pub(crate) fn eval_piecewise(expansions: &[ChebyshevInterval], x: f64) -> f64 {
    clenshaw(&expansions[get_index(expansions, x)], x)
}

/// Upstream `SuperAncillary::eval_sat`: `k` selects the property (`'P'` or
/// `'D'`), `q` the phase (0 liquid, 1 vapor; ignored for `'P'`).
pub fn eval_sat(sa: &SuperAncillaryData, t: f64, k: char, q: u8) -> f64 {
    let expansions = match (k, q) {
        ('P', _) => sa.p,
        ('D', 0) => sa.rho_l,
        ('D', 1) => sa.rho_v,
        _ => panic!("bad key {k:?}/Q={q} for superancillary eval_sat"),
    };
    eval_piecewise(expansions, t)
}
