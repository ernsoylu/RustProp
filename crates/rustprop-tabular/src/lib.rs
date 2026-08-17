//! Tabular engine — TTSE and bicubic table interpolation (port of CoolProp 8 src/Backends/Tabular)
//!
//! Fully ported: runtime table construction ([`tables`] — saturation table
//! plus the LogPH/LogPT single-phase grids), TTSE and bicubic evaluation
//! ([`ttse`], [`bicubic`]), and the low-level [`state::TabularState`]
//! covering every input pair `TabularBackend::update` serves. Like upstream
//! (`available_in_high_level()` is false), this is a LOW-LEVEL API only —
//! `props_si` rejects `TTSE&HEOS::`/`BICUBIC&HEOS::` verbatim. See the
//! [`state`] module docs for the update semantics and quirks reproduced.

pub mod bicubic;
pub mod state;
pub mod tables;
pub mod ttse;

pub use state::{Scheme, TabularInput, TabularState};

pub use rustprop_core::UPSTREAM_VERSION;
