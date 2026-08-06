//! Shared foundations for every rustprop crate.
//!
//! Architecture rule: this crate defines the *types* — parameters, state,
//! errors, engine traits, and fluid-data structures — while `rustprop-data`
//! holds only generated data *contents*. Engines depend on `rustprop-core`
//! alone, so an application links fluid data solely for the fluids it opts
//! into.
//!
//! Ports shared machinery from CoolProp 8 `include/` + `src/`.

/// The upstream CoolProp release this port tracks.
pub const UPSTREAM_VERSION: &str = "8.0.0";
