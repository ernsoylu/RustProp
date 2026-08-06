//! Classic saturation ancillaries (PLAN.md 4.3) — port of
//! `SaturationAncillaryFunction::evaluate` from
//! `src/Backends/Helmholtz/Fluids/Ancillaries.cpp` @ v8.0.0 for the forms the
//! fluid documents' `pS`/`rhoL`/`rhoV` blocks use: type `"rhoLnoexp"` is the
//! non-exponential form, every other tag ("pV", "rhoV", ...) the exponential
//! form. The `rational_polynomial` type (hL/hLV/sL/sLV blocks) is Phase 12
//! scope.

use rustprop_core::fluid::SaturationAncillary;

/// Upstream `SaturationAncillaryFunction::evaluate`.
///
/// Above the reducing temperature (theta < 0) this returns NaN, mirroring
/// upstream's deliberate #1611 semantics: callers probing supercritical
/// temperatures rely on well-behaved-but-invalid rather than an exception.
pub fn evaluate(anc: &SaturationAncillary, t: f64) -> f64 {
    let theta = 1.0 - t / anc.t_r;
    if theta < 0.0 {
        return f64::NAN;
    }
    // Upstream fills s[i] then left-folds with std::accumulate — same order.
    let mut summer = 0.0;
    for i in 0..anc.n.len() {
        summer += anc.n[i] * theta.powf(anc.t[i]);
    }
    if anc.anc_type == "rhoLnoexp" {
        anc.reducing_value * (1.0 + summer)
    } else {
        let tau_r_value = if anc.using_tau_r { anc.t_r / t } else { 1.0 };
        anc.reducing_value * (tau_r_value * summer).exp()
    }
}
