//! The low-level tabular state — port of `TabularBackend::update` and the
//! output accessors it serves.
//!
//! Upstream's tabular backends are LOW-LEVEL ONLY:
//! `TabularBackend::available_in_high_level()` returns false ("None of the
//! tabular methods are available from the high-level interface",
//! TabularBackends.h:1077), so `PropsSI("...", "TTSE&HEOS::Water")` is
//! rejected before any state update. This type is the equivalent of
//! `AbstractState::factory("TTSE&HEOS", "Water")`.
//!
//! Upstream builds both single-phase tables in `check_tables()` and caches
//! them to disk. This port builds the LogPT table eagerly and the LogPH table
//! on first use — the LogPH grid costs an (h, p) flash per node, two orders
//! of magnitude slower than the LogPT grid's direct evaluations. The choice
//! affects when work happens, never a result.

use crate::bicubic::{self, CellCoeffGrid};
use crate::tables::{GridKind, GriddedTable, HUGE, SatTable, TransportSource};
use crate::ttse;
use rustprop_core::cformat::fmt_g;
use rustprop_core::params::{Param, Phase};
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

/// The input pairs `TabularBackend::update` accepts. Every other pair raises
/// upstream's "Sorry, but this set of inputs is not supported for Tabular
/// backend", which is why this enum has no other variants.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TabularInput {
    /// (hmolar, p)
    HmolarP,
    /// (p, T)
    PT,
    /// (p, umolar)
    PUmolar,
    /// (p, smolar)
    PSmolar,
    /// (rhomolar, p)
    DmolarP,
    /// (smolar, T)
    SmolarT,
    /// (rhomolar, T)
    DmolarT,
    /// (p, Q)
    PQ,
    /// (Q, T)
    QT,
}

/// Upstream `selected_table`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Selected {
    None,
    Ph,
    Pt,
}

/// A tabular state over one fluid's tables.
pub struct TabularState<'a> {
    scheme: Scheme,
    flash: &'a PtFlash,
    transport: Option<&'a dyn TransportSource>,
    nx: usize,
    ny: usize,
    logpt: GriddedTable,
    coeffs_pt: CellCoeffGrid,
    logph: Option<(GriddedTable, CellCoeffGrid)>,
    sat: SatTable,

    // Live state — upstream's `_T`, `_p`, `_hmolar`, `_Q` and the cached
    // indices, with NaN standing in for "not a valid number".
    single_phase: bool,
    selected: Selected,
    t: f64,
    p: f64,
    hmolar: f64,
    q: f64,
    cell: (usize, usize),
    sat_il: usize,
    sat_iv: usize,
    phase: Phase,
}

impl<'a> TabularState<'a> {
    /// Build the LogPT table and the saturation table for `flash`.
    pub fn new(
        scheme: Scheme,
        flash: &'a PtFlash,
        nx: usize,
        ny: usize,
        transport: Option<&'a dyn TransportSource>,
    ) -> Result<Self> {
        let logpt = GriddedTable::build(flash, GridKind::LogPT, nx, ny, transport)?;
        let coeffs_pt = CellCoeffGrid::build(&logpt);
        let sat = SatTable::build(flash, SatTable::DEFAULT_N, transport)?;
        Ok(TabularState {
            scheme,
            flash,
            transport,
            nx,
            ny,
            logpt,
            coeffs_pt,
            logph: None,
            sat,
            single_phase: false,
            selected: Selected::None,
            t: f64::NAN,
            p: f64::NAN,
            hmolar: f64::NAN,
            q: -1000.0,
            cell: (0, 0),
            sat_il: 0,
            sat_iv: 0,
            phase: Phase::Unknown,
        })
    }

    /// Upstream default grids (200x200 single phase, 1000-point saturation).
    pub fn with_defaults(
        scheme: Scheme,
        flash: &'a PtFlash,
        transport: Option<&'a dyn TransportSource>,
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

    /// `calc_phase()`.
    pub fn phase(&self) -> Phase {
        self.phase
    }

    /// Build the LogPH table if this is its first use (see the module note).
    fn ensure_ph(&mut self) -> Result<()> {
        if self.logph.is_none() {
            let table = GriddedTable::build(
                self.flash,
                GridKind::LogPH,
                self.nx,
                self.ny,
                self.transport,
            )?;
            let coeffs = CellCoeffGrid::build(&table);
            self.logph = Some((table, coeffs));
        }
        Ok(())
    }

    fn ph(&self) -> &(GriddedTable, CellCoeffGrid) {
        self.logph.as_ref().expect("LogPH table built")
    }

    /// `find_native_nearest_good_indices`: the nearest good NODE for TTSE,
    /// the nearest good CELL (with the alternate remap) for bicubic.
    fn native_good_indices(
        scheme: Scheme,
        table: &GriddedTable,
        coeffs: &CellCoeffGrid,
        x: f64,
        y: f64,
    ) -> Result<(usize, usize)> {
        match scheme {
            Scheme::Ttse => table.find_native_nearest_good_neighbor(x, y),
            Scheme::Bicubic => {
                let (i, j) = table.find_native_nearest_good_cell(x, y)?;
                bicubic_good_cell(coeffs, i, j, x, y, "y=")
            }
        }
    }

    /// `update(input_pair, val1, val2)`.
    pub fn update(&mut self, pair: TabularInput, val1: f64, val2: f64) -> Result<()> {
        // Upstream clears the state at the top of update().
        self.t = f64::NAN;
        self.p = f64::NAN;
        self.hmolar = f64::NAN;
        self.q = -1000.0;
        self.selected = Selected::None;
        self.phase = Phase::Unknown;

        match pair {
            TabularInput::HmolarP => self.update_hp(val1, val2),
            TabularInput::PT => self.update_pt(val1, val2),
            TabularInput::PUmolar => self.update_ph_pair(val1, Param::Umolar, val2),
            TabularInput::PSmolar => self.update_ph_pair(val1, Param::Smolar, val2),
            TabularInput::DmolarP => self.update_ph_pair(val2, Param::Dmolar, val1),
            TabularInput::SmolarT => self.update_pt_pair(val2, Param::Smolar, val1),
            TabularInput::DmolarT => self.update_pt_pair(val2, Param::Dmolar, val1),
            TabularInput::PQ => self.update_pq(val1, val2),
            TabularInput::QT => self.update_qt(val1, val2),
        }
    }

    /// The `HmolarP_INPUTS` branch.
    fn update_hp(&mut self, hmolar: f64, p: f64) -> Result<()> {
        self.ensure_ph()?;
        self.hmolar = hmolar;
        self.p = p;
        if !self.ph().0.in_bounds(hmolar, p) {
            self.single_phase = false;
            return Err(Error::Value(format!(
                "inputs are not in range, hmolar={}, p={}",
                fmt_g(hmolar),
                fmt_g(p)
            )));
        }
        self.single_phase = true;
        if let Some(inside) = self.sat.is_inside(Param::P, p, Param::Hmolar, hmolar)? {
            self.single_phase = false;
            self.q = (hmolar - inside.y_l) / (inside.y_v - inside.y_l);
            if !(0.0..=1.0).contains(&self.q) {
                return Err(Error::Value(format!(
                    "vapor quality is not in (0,1) for hmolar: {} p: {}, hL: {} hV: {} ",
                    fmt_g(hmolar),
                    fmt_g(p),
                    fmt_g(inside.y_l),
                    fmt_g(inside.y_v)
                )));
            }
            self.sat_il = inside.il;
            self.sat_iv = inside.iv;
            self.phase = Phase::Twophase;
            return Ok(());
        }
        self.selected = Selected::Ph;
        let (table, coeffs) = self.ph();
        self.cell = Self::native_good_indices(self.scheme, table, coeffs, hmolar, p)?;
        self.recalculate_singlephase_phase()
    }

    /// The `PT_INPUTS` branch. The imposed-phase bump ladder is not ported —
    /// `imposed_phase_index` has no entry point in this API.
    pub fn update_pt(&mut self, p: f64, t: f64) -> Result<()> {
        self.p = p;
        self.t = t;
        if !self.logpt.in_bounds(t, p) {
            self.single_phase = false;
            return Err(Error::Value(format!(
                "inputs are not in range, p={} Pa, T={} K",
                fmt_g(p),
                fmt_g(t)
            )));
        }
        self.single_phase = true;
        if self.sat.is_inside_pt(p, t)?.is_some() {
            // Upstream's message names TTSE for both tabular backends.
            self.single_phase = false;
            return Err(Error::Value(
                "P,T with TTSE cannot be two-phase for now".into(),
            ));
        }
        self.selected = Selected::Pt;
        self.cell = Self::native_good_indices(self.scheme, &self.logpt, &self.coeffs_pt, t, p)?;
        self.recalculate_singlephase_phase()
    }

    /// The `PUmolar_INPUTS` / `PSmolar_INPUTS` / `DmolarP_INPUTS` branch:
    /// locate the node from (p, other), then invert the LogPH table's x axis
    /// for hmolar.
    fn update_ph_pair(&mut self, p: f64, otherkey: Param, otherval: f64) -> Result<()> {
        self.p = p;
        self.single_phase = true;
        if let Some(inside) = self.sat.is_inside(Param::P, p, otherkey, otherval)? {
            self.single_phase = false;
            self.q = two_phase_quality(otherkey, otherval, inside.y_l, inside.y_v);
            if !(0.0..=1.0).contains(&self.q) {
                return Err(Error::Value(format!(
                    "vapor quality is not in (0,1) for {}: {} p: {}",
                    otherkey.short_name(),
                    fmt_g(otherval),
                    fmt_g(p)
                )));
            }
            self.sat_il = inside.il;
            self.sat_iv = inside.iv;
            self.phase = Phase::Twophase;
            return Ok(());
        }
        self.ensure_ph()?;
        self.selected = Selected::Ph;
        let scheme = self.scheme;
        let (cell, hmolar) = {
            let (table, coeffs) = self.ph();
            let (i, j) = table.find_nearest_neighbor(Param::P, p, otherkey, otherval)?;
            let (i, j) = match scheme {
                Scheme::Ttse => (i, j),
                Scheme::Bicubic => bicubic_good_cell(coeffs, i, j, p, otherval, "y = ")?,
            };
            let h = match scheme {
                Scheme::Ttse => ttse::invert_single_phase_x(table, otherkey, otherval, p, i, j)?,
                Scheme::Bicubic => {
                    bicubic::invert_single_phase_x(table, coeffs, otherkey, otherval, p, i, j)?
                }
            };
            ((i, j), h)
        };
        self.cell = cell;
        self.hmolar = hmolar;
        self.recalculate_singlephase_phase()
    }

    /// The `SmolarT_INPUTS` / `DmolarT_INPUTS` branch: locate the node from
    /// (T, other), then invert the LogPT table's y axis for p.
    fn update_pt_pair(&mut self, t: f64, otherkey: Param, otherval: f64) -> Result<()> {
        self.t = t;
        self.single_phase = true;
        if let Some(inside) = self.sat.is_inside(Param::T, t, otherkey, otherval)? {
            self.single_phase = false;
            self.q = two_phase_quality(otherkey, otherval, inside.y_l, inside.y_v);
            if !(0.0..=1.0).contains(&self.q) {
                return Err(Error::Value(format!(
                    "vapor quality is not in (0,1) for {}: {} T: {}",
                    otherkey.short_name(),
                    fmt_g(otherval),
                    fmt_g(t)
                )));
            }
            self.sat_il = inside.il;
            self.sat_iv = inside.iv;
            self.p = self
                .sat
                .evaluate(Param::P, t, self.q, inside.il, inside.iv)?;
            self.phase = Phase::Twophase;
            return Ok(());
        }
        self.selected = Selected::Pt;
        let (i, j) = self
            .logpt
            .find_nearest_neighbor(Param::T, t, otherkey, otherval)?;
        let (i, j) = match self.scheme {
            Scheme::Ttse => (i, j),
            Scheme::Bicubic => bicubic_good_cell(&self.coeffs_pt, i, j, t, otherval, "y = ")?,
        };
        self.cell = (i, j);
        self.p = match self.scheme {
            Scheme::Ttse => ttse::invert_single_phase_y(&self.logpt, otherkey, otherval, t, i, j)?,
            Scheme::Bicubic => bicubic::invert_single_phase_y(
                &self.logpt,
                &self.coeffs_pt,
                otherkey,
                otherval,
                t,
                i,
                j,
            )?,
        };
        self.recalculate_singlephase_phase()
    }

    /// The `PQ_INPUTS` branch.
    fn update_pq(&mut self, p: f64, q: f64) -> Result<()> {
        self.p = p;
        self.q = q;
        self.single_phase = false;
        if !(0.0..=1.0).contains(&q) {
            return Err(Error::Value(format!(
                "vapor quality [{}] is not in (0,1)",
                fmt_g(q)
            )));
        }
        let inside = self
            .sat
            .is_inside(Param::P, p, Param::Q, q)?
            .ok_or_else(|| {
                Error::Value("Not possible to determine whether pressure is inside or not".into())
            })?;
        self.t = q * inside.y_v + (1.0 - q) * inside.y_l;
        self.sat_il = inside.il;
        self.sat_iv = inside.iv;
        self.phase = Phase::Twophase;
        Ok(())
    }

    /// The `QT_INPUTS` branch.
    ///
    /// UPSTREAM QUIRK, REPRODUCED: the `is_inside` result is discarded with a
    /// `(void)` cast, so a temperature outside the saturation table leaves
    /// `pL`/`pV` at their `_HUGE` initialisers and `_p` comes out infinite
    /// instead of raising. Only PQ carries the "not possible to determine"
    /// guard.
    fn update_qt(&mut self, q: f64, t: f64) -> Result<()> {
        self.q = q;
        self.t = t;
        self.single_phase = false;
        if !(0.0..=1.0).contains(&q) {
            return Err(Error::Value(format!(
                "vapor quality [{}] is not in (0,1)",
                fmt_g(q)
            )));
        }
        let (p_l, p_v, il, iv) = match self.sat.is_inside(Param::T, t, Param::Q, q)? {
            Some(r) => (r.y_l, r.y_v, r.il, r.iv),
            None => (HUGE, HUGE, 0, 0),
        };
        self.p = q * p_v + (1.0 - q) * p_l;
        self.sat_il = il;
        self.sat_iv = iv;
        self.phase = Phase::Twophase;
        Ok(())
    }

    /// `recalculate_singlephase_phase()`.
    fn recalculate_singlephase_phase(&mut self) -> Result<()> {
        let p = self.keyed_output(Param::P)?;
        let t = self.keyed_output(Param::T)?;
        self.phase = if p > self.flash.p_critical() {
            if t > self.flash.t_critical() {
                Phase::Supercritical
            } else {
                Phase::SupercriticalLiquid
            }
        } else if t > self.flash.t_critical() {
            Phase::SupercriticalGas
        } else if self.keyed_output(Param::Dmolar)? > self.flash.rhomolar_critical() {
            Phase::Liquid
        } else {
            Phase::Gas
        };
        Ok(())
    }

    /// The selected table, its coefficients and the (x, y) it is queried at.
    fn selected_table(&self) -> Result<(&GriddedTable, &CellCoeffGrid, f64, f64)> {
        match self.selected {
            Selected::Pt => Ok((&self.logpt, &self.coeffs_pt, self.t, self.p)),
            Selected::Ph => {
                let (t, c) = self.ph();
                Ok((t, c, self.hmolar, self.p))
            }
            Selected::None => Err(Error::Value("table not selected".into())),
        }
    }

    /// Evaluate one molar output off the selected table and cached indices.
    fn evaluate(&self, out: Param) -> Result<f64> {
        let (i, j) = self.cell;
        let (table, coeffs, x, y) = self.selected_table()?;
        match self.scheme {
            Scheme::Ttse => ttse::evaluate_single_phase(table, out, x, y, i, j),
            Scheme::Bicubic => bicubic::evaluate_single_phase(table, coeffs, out, x, y, i, j),
        }
    }

    /// Transport off the cached indices — for TTSE those are the nearest good
    /// *node*, reused as a cell corner.
    fn evaluate_transport(&self, out: Param) -> Result<f64> {
        let (i, j) = self.cell;
        let (table, _, x, y) = self.selected_table()?;
        match self.scheme {
            Scheme::Ttse => ttse::evaluate_single_phase_transport(table, out, x, y, i, j),
            Scheme::Bicubic => bicubic::evaluate_single_phase_transport(table, out, x, y, i, j),
        }
    }

    /// A keyed output at the current state — the `calc_*` family.
    pub fn keyed_output(&self, out: Param) -> Result<f64> {
        let mm = self.flash.eos.molar_mass;
        let molar = match out {
            Param::P => return Ok(self.p),
            Param::Q => return Ok(self.q),
            Param::T => {
                // The PT table echoes its input; the PH table evaluates T.
                // In two phase `_T` wins whenever the pair set it.
                if self.single_phase {
                    return match self.selected {
                        Selected::Pt => Ok(self.t),
                        _ => self.evaluate(Param::T),
                    };
                }
                if self.t.is_finite() {
                    return Ok(self.t);
                }
                return self
                    .sat
                    .evaluate(Param::T, self.p, self.q, self.sat_il, self.sat_iv);
            }
            Param::Hmolar | Param::Hmass => Param::Hmolar,
            Param::Dmolar | Param::Dmass => Param::Dmolar,
            Param::Smolar | Param::Smass => Param::Smolar,
            Param::Umolar | Param::Umass => Param::Umolar,
            Param::Viscosity | Param::Conductivity => {
                return if self.single_phase {
                    self.evaluate_transport(out)
                } else {
                    self.sat
                        .evaluate(out, self.p, self.q, self.sat_il, self.sat_iv)
                };
            }
            Param::Cpmolar | Param::Cvmolar | Param::SpeedSound if !self.single_phase => {
                return self
                    .sat
                    .evaluate(out, self.p, self.q, self.sat_il, self.sat_iv);
            }
            other => {
                return Err(Error::NotImplemented(format!(
                    "output parameter {} is not ported for the tabular backends",
                    other.short_name()
                )));
            }
        };
        let val = if !self.single_phase {
            self.sat
                .evaluate(molar, self.p, self.q, self.sat_il, self.sat_iv)?
        } else if molar == Param::Hmolar && self.selected == Selected::Ph {
            // The PH table echoes hmolar — the input, or the inverted value.
            self.hmolar
        } else {
            self.evaluate(molar)?
        };
        Ok(match out {
            Param::Dmass => val * mm,
            Param::Hmass | Param::Smass | Param::Umass => val / mm,
            _ => val,
        })
    }
}

/// The invalid-cell remap shared by `find_native_nearest_good_indices` and
/// `find_nearest_neighbor`. Upstream writes the second value with a space
/// (`y = %g`) in one and without (`y= %g`) in the other; `ysep` carries that.
fn bicubic_good_cell(
    coeffs: &CellCoeffGrid,
    i: usize,
    j: usize,
    x: f64,
    y: f64,
    ysep: &str,
) -> Result<(usize, usize)> {
    let cell = coeffs.cell(i, j);
    if cell.valid() {
        return Ok((i, j));
    }
    cell.alternate().ok_or_else(|| {
        Error::Value(format!(
            "Cell is invalid and has no good neighbors for x = {}, {ysep}{}",
            fmt_g(x),
            fmt_g(y)
        ))
    })
}

/// The two-phase quality upstream computes per "other" key: reciprocal
/// (volume) mixing for density, linear for everything else.
fn two_phase_quality(otherkey: Param, otherval: f64, y_l: f64, y_v: f64) -> f64 {
    if otherkey == Param::Dmolar {
        (1.0 / otherval - 1.0 / y_l) / (1.0 / y_v - 1.0 / y_l)
    } else {
        (otherval - y_l) / (y_v - y_l)
    }
}
