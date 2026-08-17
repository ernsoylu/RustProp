//! SVDSBTL engine — SVD-compressed tabular lookup, new in v8 (port of CoolProp 8 src/Backends/SVDSBTL)
//!
//! The EVALUATOR is fully ported ([`surface::SvdSurface`] dispatch over the
//! Hermite-basis rank-r kernel, axis transforms, boundary curves, and region
//! atlas) plus an [`artifact`] reader for flat `.svds` blobs ingested from
//! upstream surfaces by `tools/rustprop-svdgen`; the builder is deliberately
//! not ported (Eigen BDCSVD is not bitwise-reproducible). Like upstream
//! (`available_in_high_level()` is false), this is a LOW-LEVEL API only —
//! see the [`surface`] module docs.

pub mod artifact;
pub mod hermite;
pub mod region;
pub mod surface;
pub mod svd;

pub use rustprop_core::UPSTREAM_VERSION;
