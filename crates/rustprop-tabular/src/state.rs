//! The low-level tabular state — port of `TabularBackend::update`'s
//! `PT_INPUTS` branch and the output accessors it serves.
//!
//! Upstream's tabular backends are LOW-LEVEL ONLY:
//! `TabularBackend::available_in_high_level()` returns false ("None of the
//! tabular methods are available from the high-level interface",
//! TabularBackends.h:1077), so `PropsSI("...", "TTSE&HEOS::Water")` is
//! rejected before any state update. This type is the equivalent of
//! `AbstractState::factory("TTSE&HEOS", "Water")`.

use crate::bicubic::{self, CellCoeffGrid};
use crate::tables::{GridKind, GriddedTable, SatTable, TransportSource};
use crate::ttse;
use rustprop_core::params::Param;
use rustprop_core::{Error, Result};
use rustprop_heos::flash_pt::PtFlash;

/// Which interpolation scheme a state evaluates with.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scheme {
    /// `TTSE&HEOS` — second-order Taylor about the nearest node.
    Ttse,
    /// `BICUBIC&HEOS` — 16-coefficient cell interpolant.
    Bicubic,
}

/// A tabular state over one fluid's tables.
pub struct TabularState<'a> {
    scheme: Scheme,
    flash: &'a PtFlash,
    logpt: GriddedTable,
    coeffs_pt: CellCoeffGrid,
    sat: SatTable,
    // Live state
    t: f64,
    p: f64,
    cell: (usize, usize),
}

impl<'a> TabularState<'a> {
    /// Build the tables for `flash` and return a state ready to update.
    /// Upstream builds these lazily on first use and caches them to disk;
    /// this port always builds in memory (see the `tables` module docs).
    pub fn new(
        scheme: Scheme,
        flash: &'a PtFlash,
        nx: usize,
        ny: usize,
        transport: Option<&dyn TransportSource>,
    ) -> Result<Self> {
        let logpt = GriddedTable::build(flash, GridKind::LogPT, nx, ny, transport)?;
        let coeffs_pt = CellCoeffGrid::build(&logpt);
        let sat = SatTable::build(flash, SatTable::DEFAULT_N, transport)?;
        Ok(TabularState {
            scheme,
            flash,
            logpt,
            coeffs_pt,
            sat,
            t: f64::NAN,
            p: f64::NAN,
            cell: (0, 0),
        })
    }

    /// Upstream default grids (200x200 single phase, 1000-point saturation).
    pub fn with_defaults(
        scheme: Scheme,
        flash: &'a PtFlash,
        transport: Option<&dyn TransportSource>,
    ) -> Result<Self> {
        Self::new(
            scheme,
            flash,
            GriddedTable::DEFAULT_N,
            GriddedTable::DEFAULT_N,
            transport,
        )
    }

    pub fn table(&self) -> &GriddedTable {
        &self.logpt
    }

    /// `update(PT_INPUTS, p, T)`: range check, two-phase rejection, then
    /// cache the cell/node this scheme will evaluate from.
    pub fn update_pt(&mut self, p: f64, t: f64) -> Result<()> {
        self.p = p;
        self.t = t;
        if !self.logpt.in_bounds(t, p) {
            return Err(Error::Value(format!(
                "inputs are not in range, p={p} Pa, T={t} K"
            )));
        }
        if self.sat.is_inside_pt(p, t)?.is_some() {
            // Upstream's message names TTSE for both tabular backends.
            return Err(Error::Value(
                "P,T with TTSE cannot be two-phase for now".into(),
            ));
        }
        self.cell = match self.scheme {
            Scheme::Ttse => self.logpt.find_native_nearest_good_neighbor(t, p)?,
            Scheme::Bicubic => {
                let (i, j) = self.logpt.find_native_nearest_good_cell(t, p)?;
                let cell = self.coeffs_pt.cell(i, j);
                if cell.valid() {
                    (i, j)
                } else {
                    cell.alternate().ok_or_else(|| {
                        Error::Value(format!(
                            "Cell is invalid and has no good neighbors for x = {t}, y= {p}"
                        ))
                    })?
                }
            }
        };
        Ok(())
    }

    /// A keyed output at the current state.
    pub fn keyed_output(&self, out: Param) -> Result<f64> {
        let mm = self.flash.eos.molar_mass;
        let molar = match out {
            Param::T => return Ok(self.t),
            Param::P => return Ok(self.p),
            Param::Dmolar | Param::Dmass => Param::Dmolar,
            Param::Hmolar | Param::Hmass => Param::Hmolar,
            Param::Smolar | Param::Smass => Param::Smolar,
            Param::Umolar | Param::Umass => Param::Umolar,
            Param::Viscosity | Param::Conductivity => {
                // `calc_viscosity`/`calc_conductivity` interpolate bilinearly
                // over the cell rooted at the CACHED indices — for TTSE that
                // is the nearest good *node*, reused as a cell corner.
                let (i, j) = self.cell;
                return match self.scheme {
                    Scheme::Ttse => ttse::evaluate_single_phase_transport(
                        &self.logpt,
                        out,
                        self.t,
                        self.p,
                        i,
                        j,
                    ),
                    Scheme::Bicubic => bicubic::evaluate_single_phase_transport(
                        &self.logpt,
                        out,
                        self.t,
                        self.p,
                        i,
                        j,
                    ),
                };
            }
            other => {
                return Err(Error::NotImplemented(format!(
                    "output parameter {} is not ported for the tabular backends",
                    other.short_name()
                )));
            }
        };
        let (i, j) = self.cell;
        let val = match self.scheme {
            Scheme::Ttse => ttse::evaluate_single_phase(&self.logpt, molar, self.t, self.p, i, j)?,
            Scheme::Bicubic => bicubic::evaluate_single_phase(
                &self.logpt,
                &self.coeffs_pt,
                molar,
                self.t,
                self.p,
                i,
                j,
            )?,
        };
        Ok(match out {
            Param::Dmass => val * mm,
            Param::Hmass | Param::Smass | Param::Umass => val / mm,
            _ => val,
        })
    }
}
