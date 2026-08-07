//! SVDSBTL engine — SVD-compressed tabular lookup, new in v8 (port of CoolProp 8 src/Backends/SVDSBTL)
//!
//! Not yet ported. Fidelity rule: identical algorithms and data as upstream
//! CoolProp v8.0.0 — see CLAUDE.md.

pub mod artifact;
pub mod hermite;
pub mod region;
pub mod surface;
pub mod svd;

pub use rustprop_core::UPSTREAM_VERSION;
