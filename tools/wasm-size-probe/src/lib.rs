//! Minimal cdylib whose compiled `.wasm` size CI reports per feature set.
//!
//! The linker strips unreachable code, so the numbers only mean something if
//! `probe()` actually calls into every enabled engine's public entry points —
//! extend it under `#[cfg(feature = ...)]` as facade API lands (first: the
//! IF97 baseline, PLAN.md step 2.6).

// The exported symbol needs `#[unsafe(no_mangle)]` (edition 2024), which the
// workspace-wide `deny(unsafe_code)` would reject; this probe is the one
// deliberate exception and ships nowhere.
#![allow(unsafe_code)]

#[unsafe(no_mangle)]
pub extern "C" fn probe() -> usize {
    rustprop::UPSTREAM_VERSION.len()
}

/// Exercises the full IF97 dispatch so the linker keeps the engine
/// (PLAN.md 2.6 size baseline). Runtime arguments prevent const-folding.
#[cfg(feature = "if97")]
#[unsafe(no_mangle)]
pub extern "C" fn probe_if97(t: f64, p: f64) -> f64 {
    use rustprop_core::Param;
    rustprop::if97_api::props(Param::Hmass, Param::T, t, Param::P, p).unwrap_or(f64::NAN)
}

/// Exercises the HEOS engine over every compiled-in fluid (PT flash +
/// enthalpy) so the linker keeps the engine and all enabled fluid data.
/// With `heos-water` that is one fluid; with `heos-all-fluids`, all 130 —
/// the difference is the per-fluid data cost the modularity story is about.
#[cfg(feature = "heos-data")]
#[unsafe(no_mangle)]
pub extern "C" fn probe_heos(t: f64, p: f64) -> f64 {
    let mut acc = 0.0;
    for (_name, data) in rustprop_data::fluids::all() {
        let flash = rustprop_heos::flash_pt::PtFlash::new(data);
        if let Ok((rho, _phase)) = flash.pt_flash(t, p) {
            acc += flash.eos.hmolar(t, rho);
        }
    }
    acc
}

/// Exercises the cubic engine over the whole 116-fluid table (SRK PT flash
/// + enthalpy) so the linker keeps the engine and the `cubic-fluids` data.
#[cfg(feature = "cubics")]
#[unsafe(no_mangle)]
pub extern "C" fn probe_cubics(t: f64, p: f64) -> f64 {
    let mut acc = 0.0;
    for out in ["D", "Hmolar", "Smolar"] {
        acc += rustprop::props_si(out, "T", t, "P", p, "SRK::n-Propane").unwrap_or(f64::NAN);
        acc += rustprop::props_si(out, "T", t, "P", p, "PR::Water").unwrap_or(f64::NAN);
    }
    acc
}

/// Exercises the incompressible engine across the whole table.
#[cfg(feature = "incompressible")]
#[unsafe(no_mangle)]
pub extern "C" fn probe_incomp(t: f64, p: f64) -> f64 {
    rustprop::props_si("D", "T", t, "P", p, "INCOMP::DowQ").unwrap_or(f64::NAN)
        + rustprop::props_si("H", "T", t, "P", p, "INCOMP::MEG[0.4]").unwrap_or(f64::NAN)
}

/// Exercises the humid-air engine (pulls HEOS Water + Air + IF97).
#[cfg(feature = "humid-air")]
#[unsafe(no_mangle)]
pub extern "C" fn probe_ha(t: f64, p: f64) -> f64 {
    rustprop::ha_props_si("H", "T", t, "P", p, "R", 0.5).unwrap_or(f64::NAN)
}
