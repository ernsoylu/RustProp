//! Tabular engine — TTSE and bicubic table interpolation (port of CoolProp 8 src/Backends/Tabular)
//!
//! Not yet ported. Fidelity rule: identical algorithms and data as upstream
//! CoolProp v8.0.0 — see CLAUDE.md.

pub mod bicubic;
pub mod tables;
pub mod ttse;

pub use rustprop_core::UPSTREAM_VERSION;
