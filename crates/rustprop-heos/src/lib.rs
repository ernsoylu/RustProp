//! HEOS engine — multiparameter Helmholtz EOS for pure fluids and mixtures (port of CoolProp 8 src/Backends/Helmholtz)
//!
//! Fully ported and golden-verified against the CoolProp 8.0.0 wheel: the
//! Helmholtz term machinery ([`alpha`]), classic and super-ancillaries
//! ([`ancillary`], [`superancillary`]), every `PropsSI`-reachable flash pair
//! for the 130 pure and 6 pseudo-pure fluids ([`flash_pt`], [`flash_px`],
//! [`flash_hs`]), melting lines ([`melting`]), transport and surface tension
//! ([`transport`]), the generic partial-derivative machinery ([`derivs`]),
//! and the mixture model — GERG reducing/departure functions, VLE, Michelsen
//! stability/split, and the sweep pairs ([`mixture`], [`mixture_vle`],
//! [`mixture_flash`], [`mixture_stability`], [`mixture_sweep`]).

pub mod alpha;
pub mod ancillary;
mod chebappr;
pub mod derivs;
pub mod flash_hs;
pub mod flash_pt;
pub mod flash_px;
pub mod melting;
pub mod mixture;
pub mod mixture_flash;
pub mod mixture_stability;
pub mod mixture_sweep;
pub mod mixture_vle;
pub mod props;
pub mod saturation;
pub mod solvers;
pub mod superancillary;
pub mod transport;

pub use alpha::{HelmholtzDerivs, HelmholtzEos};
pub use rustprop_core::UPSTREAM_VERSION;
