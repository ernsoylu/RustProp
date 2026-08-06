//! Generated CoolProp 8 fluid and mixture data.
//!
//! Contents are emitted by `tools/rustprop-datagen` from upstream CoolProp
//! v8.0.0 JSON and must never be edited by hand. Every fluid sits behind its
//! own Cargo feature (`default = []`): data dominates WebAssembly binary
//! size, so nothing is compiled in unless an application opts in.
//!
//! Data *types* live in `rustprop-core`; only *contents* live here.

pub use rustprop_core::UPSTREAM_VERSION;

pub mod fluids;
