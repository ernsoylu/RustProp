//! HEOS engine — multiparameter Helmholtz EOS for pure fluids and mixtures (port of CoolProp 8 src/Backends/Helmholtz)
//!
//! Phase 4 of PLAN.md. Currently ported: the Helmholtz term machinery
//! ([`alpha`]) for the families Water uses, golden-verified term-by-term
//! against the CoolProp 8.0.0 wheel's `AbstractState` derivative accessors.

pub mod alpha;
pub mod ancillary;
mod chebappr;
pub mod flash_hs;
pub mod flash_pt;
pub mod flash_px;
pub mod melting;
pub mod mixture;
pub mod props;
pub mod saturation;
pub mod solvers;
pub mod superancillary;
pub mod transport;

pub use alpha::{HelmholtzDerivs, HelmholtzEos};
pub use rustprop_core::UPSTREAM_VERSION;
