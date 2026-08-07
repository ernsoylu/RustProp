//! rustprop — pure-Rust port of CoolProp 8 for WebAssembly-first use.
//!
//! Facade crate: exposes the `PropsSI`-style API and re-exports each
//! calculation engine behind a Cargo feature. Default features are empty by
//! design — applications opt into exactly the engines (and, via
//! `rustprop-data`, the fluids) a calculation needs, keeping compiled WASM
//! minimal.
//!
//! # Example: the `PropsSI` string API
//!
//! ```
//! # #[cfg(feature = "all-backends")] {
//! // Equivalent to PropsSI("Dmolar", "T", 300, "P", 101325, "Water").
//! let d = rustprop::props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water").unwrap();
//! // Golden value from the CoolProp 8.0.0 oracle:
//! assert!(((d - 55317.35277350119) / d).abs() < 1e-8);
//! # }
//! # #[cfg(feature = "if97")] {
//! use rustprop::Param;
//! // Equivalent to PropsSI("H", "T", 300, "P", 101325, "IF97::Water").
//! let h = rustprop::if97_api::props(Param::Hmass, Param::T, 300.0, Param::P, 101325.0).unwrap();
//! assert!(((h - 112665.04341853978) / h).abs() < 1e-11);
//! # }
//! ```

pub use rustprop_core::{Error, InputPair, Param, Result, UPSTREAM_VERSION};

#[cfg(feature = "cubics")]
pub use rustprop_cubics as cubics;
#[cfg(feature = "heos")]
pub use rustprop_heos as heos;
#[cfg(feature = "humid-air")]
pub use rustprop_humid_air as humid_air;
#[cfg(feature = "if97")]
pub mod if97_api;
#[cfg(any(
    feature = "heos",
    feature = "if97",
    feature = "cubics",
    feature = "incompressible"
))]
pub mod props_api;
#[cfg(any(
    feature = "heos",
    feature = "if97",
    feature = "cubics",
    feature = "incompressible"
))]
pub use props_api::props_si;

/// Upstream `HAPropsSI` (humid air). See `rustprop_humid_air` for the
/// deliberate-reproduction notes.
#[cfg(feature = "humid-air")]
#[allow(clippy::too_many_arguments)]
pub fn ha_props_si(
    output: &str,
    n1: &str,
    v1: f64,
    n2: &str,
    v2: f64,
    n3: &str,
    v3: f64,
) -> rustprop_core::Result<f64> {
    use std::sync::OnceLock;
    static HA: OnceLock<rustprop_humid_air::HumidAir> = OnceLock::new();
    let ha = HA.get_or_init(|| {
        rustprop_humid_air::HumidAir::new(
            &rustprop_data::fluids::water::WATER,
            &rustprop_data::fluids::air::AIR,
        )
    });
    ha.ha_props_si(output, n1, v1, n2, v2, n3, v3)
}
#[cfg(feature = "if97")]
pub use rustprop_if97 as if97;
#[cfg(feature = "incompressible")]
pub use rustprop_incompressible as incompressible;
#[cfg(feature = "pcsaft")]
pub use rustprop_pcsaft as pcsaft;
#[cfg(feature = "svdsbtl")]
pub use rustprop_svdsbtl as svdsbtl;
#[cfg(feature = "tabular")]
pub use rustprop_tabular as tabular;

#[cfg(test)]
mod tests {
    #[test]
    fn tracks_coolprop_8() {
        assert_eq!(super::UPSTREAM_VERSION, "8.0.0");
    }
}
