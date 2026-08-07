//! Table construction — port of upstream
//! `PureFluidSaturationTableData::build` and
//! `SinglePhaseGriddedTableData::build` (`src/Backends/Tabular/
//! TabularBackends.cpp`).
//!
//! Tables are built AT RUNTIME from a ported source engine; upstream's
//! msgpack/zlib disk cache (`~/.CoolProp/Tabular/...`) is deliberately NOT
//! ported — a WASM-first library has no home directory, and the cache only
//! ever holds what this code can regenerate. Every other observable is
//! reproduced: grid layout, spacing, limits, the two-phase holes, and the
//! nearest-good-neighbor fixups.
//!
//! Holes use `+inf` (upstream's `_HUGE`), and "is this node good?" is
//! upstream's `ValidNumber` — finite and not NaN.

use rustprop_core::params::Param;
use rustprop_core::{Error, Result};
use rustprop_heos::derivs::StateDerivs;
use rustprop_heos::flash_pt::PtFlash;

/// Transport evaluation, injected by the caller: the ECS conformal-state
/// resolver needs the fluid registry, which lives above this crate. `None`
/// from either method leaves the cell a hole — exactly what upstream's
/// `try { visc = AS->viscosity(); } catch {}` does.
pub trait TransportSource {
    fn viscosity(&self, t: f64, rhomolar: f64) -> Option<f64>;
    fn conductivity(&self, t: f64, rhomolar: f64) -> Option<f64>;
}

/// Upstream `_HUGE` — the value an unfilled table cell keeps.
pub const HUGE: f64 = f64::INFINITY;

/// Upstream `ValidNumber`.
#[inline]
pub fn valid_number(v: f64) -> bool {
    v.is_finite()
}

/// One tabulated property with its first and second derivatives on the
/// (x, y) grid — upstream stores these as the flat `X(T) X(dTdx) X(d2Tdx2)`
/// … matrix family; grouping them keeps the same data with less repetition.
#[derive(Clone)]
pub struct PropGrid {
    pub val: Vec<Vec<f64>>,
    pub dx: Vec<Vec<f64>>,
    pub dy: Vec<Vec<f64>>,
    pub dxx: Vec<Vec<f64>>,
    pub dxy: Vec<Vec<f64>>,
    pub dyy: Vec<Vec<f64>>,
}

impl PropGrid {
    fn new(nx: usize, ny: usize) -> Self {
        let m = || vec![vec![HUGE; ny]; nx];
        PropGrid {
            val: m(),
            dx: m(),
            dy: m(),
            dxx: m(),
            dxy: m(),
            dyy: m(),
        }
    }
}

/// Which gridded table layout (upstream `LogPHTable` / `LogPTTable`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GridKind {
    /// x = hmolar (linear), y = p (log)
    LogPH,
    /// x = T (linear), y = p (log)
    LogPT,
}

impl GridKind {
    pub fn xkey(self) -> Param {
        match self {
            GridKind::LogPH => Param::Hmolar,
            GridKind::LogPT => Param::T,
        }
    }
    /// Both layouts use pressure on y, log spaced.
    pub fn ykey(self) -> Param {
        Param::P
    }
}

/// A single-phase gridded table (upstream `SinglePhaseGriddedTableData`).
pub struct GriddedTable {
    pub kind: GridKind,
    pub nx: usize,
    pub ny: usize,
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub xvec: Vec<f64>,
    pub yvec: Vec<f64>,
    pub t: PropGrid,
    pub p: PropGrid,
    pub rhomolar: PropGrid,
    pub hmolar: PropGrid,
    pub smolar: PropGrid,
    pub umolar: PropGrid,
    /// Transport properties carry no derivatives upstream either.
    pub visc: Vec<Vec<f64>>,
    pub cond: Vec<Vec<f64>>,
    /// Nearest good (i, j) for cells whose T is a hole.
    pub nearest_i: Vec<Vec<usize>>,
    pub nearest_j: Vec<Vec<usize>>,
}

impl GriddedTable {
    /// Upstream default grid (`TABULAR_NX` / `TABULAR_NY`, both 200).
    pub const DEFAULT_N: usize = 200;

    /// `set_limits()` for the requested layout, then `build()`.
    pub fn build(
        flash: &PtFlash,
        kind: GridKind,
        nx: usize,
        ny: usize,
        transport: Option<&dyn TransportSource>,
    ) -> Result<Self> {
        let (xmin, xmax, ymin, ymax) = Self::set_limits(flash, kind)?;
        let mut table = GriddedTable {
            kind,
            nx,
            ny,
            xmin,
            xmax,
            ymin,
            ymax,
            xvec: vec![0.0; nx],
            yvec: vec![0.0; ny],
            t: PropGrid::new(nx, ny),
            p: PropGrid::new(nx, ny),
            rhomolar: PropGrid::new(nx, ny),
            hmolar: PropGrid::new(nx, ny),
            smolar: PropGrid::new(nx, ny),
            umolar: PropGrid::new(nx, ny),
            visc: vec![vec![HUGE; ny]; nx],
            cond: vec![vec![HUGE; ny]; nx],
            nearest_i: vec![vec![usize::MAX; ny]; nx],
            nearest_j: vec![vec![usize::MAX; ny]; nx],
        };
        table.fill(flash, transport)?;
        table.make_good_neighbors();
        Ok(table)
    }

    /// `LogPHTable::set_limits` / `LogPTTable::set_limits`.
    fn set_limits(flash: &PtFlash, kind: GridKind) -> Result<(f64, f64, f64, f64)> {
        let data = flash.fluid();
        let tmin = flash.t_triple().max(data.eos.sat_min_liquid.t);
        // Saturated liquid at Tmin sets the low corner for both layouts.
        let sat = flash.qt_state(tmin, 0.0)?;
        let ymin = sat.p();
        let ymax = data.eos.p_max;
        match kind {
            GridKind::LogPH => {
                let xmin = flash.state_hmolar(&sat);
                // Upstream checks BOTH enthalpies on the Tmax isotherm and
                // takes the larger: the ideal-gas-limit one and the pmax one.
                let t_hi = 1.499 * data.eos.t_max;
                let xmax1 = flash.eos.hmolar(t_hi, 1e-10);
                let (rho2, _) = flash.pt_flash(t_hi, ymax)?;
                let xmax2 = flash.eos.hmolar(t_hi, rho2);
                Ok((xmin, xmax1.max(xmax2), ymin, ymax))
            }
            GridKind::LogPT => Ok((tmin, data.eos.t_max * 1.499, ymin, ymax)),
        }
    }

    /// The build loop: upstream computes x and y INSIDE the loop with its own
    /// spacing expressions (note the asymmetry — x uses
    /// `log(xmax) - log(xmin)` while y uses `log(ymax/ymin)`), then writes
    /// them back into xvec/yvec.
    fn fill(&mut self, flash: &PtFlash, transport: Option<&dyn TransportSource>) -> Result<()> {
        for i in 0..self.nx {
            // Both layouts are linear in x (logx = false upstream).
            let x = self.xmin + (self.xmax - self.xmin) / ((self.nx - 1) as f64) * (i as f64);
            self.xvec[i] = x;
            for j in 0..self.ny {
                // logy = true for both layouts.
                let y = (self.ymin.ln()
                    + (self.ymax / self.ymin).ln() / ((self.ny - 1) as f64) * (j as f64))
                    .exp();
                self.yvec[j] = y;

                // Update the state; failures leave the cell as a hole.
                let state = match self.kind {
                    GridKind::LogPH => flash.hmolar_p_state(x, y),
                    GridKind::LogPT => flash.pt_flash(x, y).map(|(rho, phase)| {
                        rustprop_heos::flash_px::HeosState::SinglePhase {
                            t: x,
                            p: y,
                            rhomolar: rho,
                            phase,
                            q: -1.0,
                        }
                    }),
                };
                let Ok(state) = state else { continue };
                let (t, rho) = (state.t(), state.rhomolar());
                if !valid_number(rho) {
                    continue;
                }
                // Two-phase states stay as holes (upstream: Q in [0, 1]).
                let q = state.q();
                if (0.0..=1.0).contains(&q) {
                    continue;
                }

                self.t.val[i][j] = t;
                self.p.val[i][j] = state.p();
                self.rhomolar.val[i][j] = rho;
                self.hmolar.val[i][j] = flash.eos.hmolar(t, rho);
                self.smolar.val[i][j] = flash.eos.smolar(t, rho);
                self.umolar.val[i][j] = flash.eos.umolar(t, rho);

                // Transport failures stay as holes, the rest of the cell stands.
                if let Some(tr) = transport {
                    if let Some(v) = tr.viscosity(t, rho) {
                        self.visc[i][j] = v;
                    }
                    if let Some(c) = tr.conductivity(t, rho) {
                        self.cond[i][j] = c;
                    }
                }

                let (xk, yk) = (self.kind.xkey(), self.kind.ykey());
                // One Helmholtz evaluation serves all 30 derivative queries.
                let sd = StateDerivs::new(&flash.eos, t, rho);
                for (grid, of) in [
                    (&mut self.t, Param::T),
                    (&mut self.p, Param::P),
                    (&mut self.rhomolar, Param::Dmolar),
                    (&mut self.hmolar, Param::Hmolar),
                    (&mut self.smolar, Param::Smolar),
                    (&mut self.umolar, Param::Umolar),
                ] {
                    grid.dx[i][j] = sd.first_partial_deriv(of, xk, yk)?;
                    grid.dy[i][j] = sd.first_partial_deriv(of, yk, xk)?;
                    grid.dxx[i][j] = sd.second_partial_deriv(of, xk, yk, xk, yk)?;
                    grid.dxy[i][j] = sd.second_partial_deriv(of, xk, yk, yk, xk)?;
                    grid.dyy[i][j] = sd.second_partial_deriv(of, yk, xk, yk, xk)?;
                }
            }
        }
        Ok(())
    }

    /// `make_good_neighbors()`: for every cell whose T is a hole, find the
    /// first valid neighbour in upstream's fixed offset order. Note the
    /// bounds test is `0 < iplus && iplus < Nx-1` — strict on BOTH ends, so
    /// the outermost ring is never chosen as a neighbour (upstream's own
    /// asymmetry, reproduced).
    fn make_good_neighbors(&mut self) {
        const XOFF: [isize; 8] = [-1, 1, 0, 0, -1, 1, 1, -1];
        const YOFF: [isize; 8] = [0, 0, 1, -1, -1, -1, 1, 1];
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.nearest_i[i][j] = i;
                self.nearest_j[i][j] = j;
                if !valid_number(self.t.val[i][j]) {
                    for k in 0..8 {
                        let iplus = i as isize + XOFF[k];
                        let jplus = j as isize + YOFF[k];
                        if iplus > 0
                            && (iplus as usize) < self.nx - 1
                            && jplus > 0
                            && (jplus as usize) < self.ny - 1
                            && valid_number(self.t.val[iplus as usize][jplus as usize])
                        {
                            self.nearest_i[i][j] = iplus as usize;
                            self.nearest_j[i][j] = jplus as usize;
                            break;
                        }
                    }
                }
            }
        }
    }
}

/// The pure-fluid saturation table (upstream `PureFluidSaturationTableData`),
/// log-spaced in pressure from the triple point to `0.9999 * p_critical`,
/// with the last point placed exactly at the critical point.
pub struct SatTable {
    pub n: usize,
    pub p_l: Vec<f64>,
    pub t_l: Vec<f64>,
    pub rhomolar_l: Vec<f64>,
    pub hmolar_l: Vec<f64>,
    pub smolar_l: Vec<f64>,
    pub umolar_l: Vec<f64>,
    pub logp_l: Vec<f64>,
    pub logrhomolar_l: Vec<f64>,
    pub cpmolar_l: Vec<f64>,
    pub cvmolar_l: Vec<f64>,
    pub speed_sound_l: Vec<f64>,
    pub visc_l: Vec<f64>,
    pub cond_l: Vec<f64>,
    pub logvisc_l: Vec<f64>,
    pub p_v: Vec<f64>,
    pub t_v: Vec<f64>,
    pub rhomolar_v: Vec<f64>,
    pub hmolar_v: Vec<f64>,
    pub smolar_v: Vec<f64>,
    pub umolar_v: Vec<f64>,
    pub logp_v: Vec<f64>,
    pub logrhomolar_v: Vec<f64>,
    pub cpmolar_v: Vec<f64>,
    pub cvmolar_v: Vec<f64>,
    pub speed_sound_v: Vec<f64>,
    pub visc_v: Vec<f64>,
    pub cond_v: Vec<f64>,
    pub logvisc_v: Vec<f64>,
}

impl SatTable {
    /// Upstream default point count.
    pub const DEFAULT_N: usize = 1000;

    pub fn build(
        flash: &PtFlash,
        n: usize,
        transport: Option<&dyn TransportSource>,
    ) -> Result<Self> {
        let hole = || vec![HUGE; n];
        let mut s = SatTable {
            n,
            p_l: hole(),
            t_l: hole(),
            rhomolar_l: hole(),
            hmolar_l: hole(),
            smolar_l: hole(),
            umolar_l: hole(),
            logp_l: hole(),
            logrhomolar_l: hole(),
            cpmolar_l: hole(),
            cvmolar_l: hole(),
            speed_sound_l: hole(),
            visc_l: hole(),
            cond_l: hole(),
            logvisc_l: hole(),
            p_v: hole(),
            t_v: hole(),
            rhomolar_v: hole(),
            hmolar_v: hole(),
            smolar_v: hole(),
            umolar_v: hole(),
            logp_v: hole(),
            logrhomolar_v: hole(),
            cpmolar_v: hole(),
            cvmolar_v: hole(),
            speed_sound_v: hole(),
            visc_v: hole(),
            cond_v: hole(),
            logvisc_v: hole(),
        };

        let data = flash.fluid();
        let tmin = flash.t_triple().max(data.eos.sat_min_liquid.t);
        let p_triple = flash.qt_state(tmin, 0.0)?.p();
        let pmin = p_triple;
        let pmax = 0.9999 * flash.p_critical();

        for i in 0..(n - 1) {
            // Log spaced in p. (Upstream disables DONT_CHECK_PROPERTY_LIMITS
            // for i == 0 only — this port has no such global limit check.)
            let p = (pmin.ln() + (pmax.ln() - pmin.ln()) / ((n - 1) as f64) * (i as f64)).exp();

            // Saturated liquid; a failure skips the whole point (upstream
            // `continue`s before touching the vapor branch).
            let Ok(state) = flash.pq_state(p, 0.0) else {
                continue;
            };
            let (t, rho) = (state.t(), state.rhomolar());
            s.p_l[i] = p;
            s.t_l[i] = t;
            s.rhomolar_l[i] = rho;
            s.hmolar_l[i] = flash.eos.hmolar(t, rho);
            s.smolar_l[i] = flash.eos.smolar(t, rho);
            s.umolar_l[i] = flash.eos.umolar(t, rho);
            s.logp_l[i] = p.ln();
            s.logrhomolar_l[i] = rho.ln();
            s.cpmolar_l[i] = flash.eos.cpmolar(t, rho);
            s.cvmolar_l[i] = flash.eos.cvmolar(t, rho);
            s.speed_sound_l[i] = flash.eos.speed_sound(t, rho);
            if let Some(tr) = transport {
                if let Some(v) = tr.viscosity(t, rho) {
                    s.visc_l[i] = v;
                    s.logvisc_l[i] = v.ln();
                }
                if let Some(c) = tr.conductivity(t, rho) {
                    s.cond_l[i] = c;
                }
            }

            let Ok(state) = flash.pq_state(p, 1.0) else {
                continue;
            };
            let (t, rho) = (state.t(), state.rhomolar());
            s.p_v[i] = p;
            s.t_v[i] = t;
            s.rhomolar_v[i] = rho;
            s.hmolar_v[i] = flash.eos.hmolar(t, rho);
            s.smolar_v[i] = flash.eos.smolar(t, rho);
            s.umolar_v[i] = flash.eos.umolar(t, rho);
            s.logp_v[i] = p.ln();
            s.logrhomolar_v[i] = rho.ln();
            s.cpmolar_v[i] = flash.eos.cpmolar(t, rho);
            s.cvmolar_v[i] = flash.eos.cvmolar(t, rho);
            s.speed_sound_v[i] = flash.eos.speed_sound(t, rho);
            if let Some(tr) = transport {
                if let Some(v) = tr.viscosity(t, rho) {
                    s.visc_v[i] = v;
                    s.logvisc_v[i] = v.ln();
                }
                if let Some(c) = tr.conductivity(t, rho) {
                    s.cond_v[i] = c;
                }
            }
        }

        // Last point sits at the critical point — upstream takes BOTH
        // branches from the same PQ(p_critical, Q=1) state, so the liquid and
        // vapor rows are identical there (and no cp/cv/w/transport is stored).
        let state = flash.pq_state(flash.p_critical(), 1.0)?;
        let (t, rho) = (state.t(), state.rhomolar());
        let i = n - 1;
        let (h, sm, u) = (
            flash.eos.hmolar(t, rho),
            flash.eos.smolar(t, rho),
            flash.eos.umolar(t, rho),
        );
        s.p_v[i] = state.p();
        s.t_v[i] = t;
        s.rhomolar_v[i] = rho;
        s.hmolar_v[i] = h;
        s.smolar_v[i] = sm;
        s.umolar_v[i] = u;
        s.logp_v[i] = state.p().ln();
        s.logrhomolar_v[i] = rho.ln();
        s.p_l[i] = state.p();
        s.t_l[i] = t;
        s.rhomolar_l[i] = rho;
        s.hmolar_l[i] = h;
        s.smolar_l[i] = sm;
        s.umolar_l[i] = u;
        s.logp_l[i] = state.p().ln();
        s.logrhomolar_l[i] = rho.ln();

        Ok(s)
    }

    /// Bisection for the index bracketing `p` on the (monotonic) log-p grid.
    pub fn bisect_logp(&self, logp: f64) -> Result<usize> {
        if !(self.logp_v[0]..=self.logp_v[self.n - 1]).contains(&logp) {
            return Err(Error::OutOfRange(format!(
                "pressure {} is outside the saturation table",
                logp.exp()
            )));
        }
        let (mut lo, mut hi) = (0usize, self.n - 1);
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if logp < self.logp_v[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Ok(lo)
    }
}

/// Upstream `CubicInterp`: 4-point Lagrange interpolation.
// Upstream's signature, kept argument-for-argument.
#[allow(clippy::too_many_arguments)]
pub fn cubic_interp(
    x0: f64,
    x1: f64,
    x2: f64,
    x3: f64,
    f0: f64,
    f1: f64,
    f2: f64,
    f3: f64,
    x: f64,
) -> f64 {
    let l0 = ((x - x1) * (x - x2) * (x - x3)) / ((x0 - x1) * (x0 - x2) * (x0 - x3));
    let l1 = ((x - x0) * (x - x2) * (x - x3)) / ((x1 - x0) * (x1 - x2) * (x1 - x3));
    let l2 = ((x - x0) * (x - x1) * (x - x3)) / ((x2 - x0) * (x2 - x1) * (x2 - x3));
    let l3 = ((x - x0) * (x - x1) * (x - x2)) / ((x3 - x0) * (x3 - x1) * (x3 - x2));
    l0 * f0 + l1 * f1 + l2 * f2 + l3 * f3
}

fn cubic_interp_vec(
    x: &[f64],
    y: &[f64],
    i0: usize,
    i1: usize,
    i2: usize,
    i3: usize,
    v: f64,
) -> f64 {
    cubic_interp(x[i0], x[i1], x[i2], x[i3], y[i0], y[i1], y[i2], y[i3], v)
}

/// What `PureFluidSaturationTableData::is_inside` reports when the state is
/// inside the dome: the bracketing indices and the saturation values of the
/// "other" variable on each branch.
#[derive(Clone, Copy, Debug)]
pub struct Inside {
    pub il: usize,
    pub iv: usize,
    pub y_l: f64,
    pub y_v: f64,
}

impl SatTable {
    /// The saturated liquid/vapour vectors for the "other" variable
    /// (upstream's lambda at the top of `is_inside`). `iQ` maps to the
    /// temperature vectors, matching upstream.
    fn other_vecs(&self, other: Param) -> Result<(&Vec<f64>, &Vec<f64>)> {
        Ok(match other {
            Param::T | Param::Q => (&self.t_l, &self.t_v),
            Param::Hmolar => (&self.hmolar_l, &self.hmolar_v),
            Param::Smolar => (&self.smolar_l, &self.smolar_v),
            Param::Umolar => (&self.umolar_l, &self.umolar_v),
            Param::Dmolar => (&self.rhomolar_l, &self.rhomolar_v),
            _ => {
                return Err(Error::Value("invalid input for other in is_inside".into()));
            }
        })
    }

    /// `PureFluidSaturationTableData::is_inside`: is the state inside the
    /// two-phase dome, by the saturation table? `main` must be `P` or `T`.
    ///
    /// For `other == Q` upstream returns TRUE unconditionally once the main
    /// variable is in range, after filling `y_l`/`y_v` with the saturation
    /// temperature (for `main == P`) or pressure (for `main == T`).
    pub fn is_inside(
        &self,
        main: Param,
        mainval: f64,
        other: Param,
        val: f64,
    ) -> Result<Option<Inside>> {
        let (yvec_l, yvec_v) = self.other_vecs(other)?;

        // Trivial checks on the main variable's range
        match main {
            Param::P => {
                let (pmax, pmin) = (self.p_v[self.n - 1], self.p_v[0]);
                if mainval > pmax || mainval < pmin {
                    return Ok(None);
                }
            }
            Param::T => {
                let (tmax, tmin) = (self.t_v[self.n - 1], self.t_v[0]);
                if mainval > tmax || mainval < tmin {
                    return Ok(None);
                }
            }
            _ => {
                return Err(Error::Value("invalid input for other in is_inside".into()));
            }
        }

        // Indices bounding the main variable on each branch
        let (iv, il) = match main {
            Param::P => (
                crate::ttse::bisect_vector(&self.p_v, mainval)?,
                crate::ttse::bisect_vector(&self.p_l, mainval)?,
            ),
            _ => (
                crate::ttse::bisect_vector(&self.t_v, mainval)?,
                crate::ttse::bisect_vector(&self.t_l, mainval)?,
            ),
        };
        let mut ivplus = (iv + 1).min(self.n - 1);
        let mut ilplus = (il + 1).min(self.n - 1);

        if other == Param::Q {
            ivplus = ivplus.max(3);
            ilplus = ilplus.max(3);
            let (y_v, y_l) = if main == Param::P {
                let logp = mainval.ln();
                (
                    cubic_interp_vec(
                        &self.logp_v,
                        &self.t_v,
                        ivplus - 3,
                        ivplus - 2,
                        ivplus - 1,
                        ivplus,
                        logp,
                    ),
                    cubic_interp_vec(
                        &self.logp_l,
                        &self.t_l,
                        ilplus - 3,
                        ilplus - 2,
                        ilplus - 1,
                        ilplus,
                        logp,
                    ),
                )
            } else {
                (
                    cubic_interp_vec(
                        &self.t_v,
                        &self.logp_v,
                        ivplus - 3,
                        ivplus - 2,
                        ivplus - 1,
                        ivplus,
                        mainval,
                    )
                    .exp(),
                    cubic_interp_vec(
                        &self.t_l,
                        &self.logp_l,
                        ilplus - 3,
                        ilplus - 2,
                        ilplus - 1,
                        ilplus,
                        mainval,
                    )
                    .exp(),
                )
            };
            return Ok(Some(Inside {
                il: ilplus - 1,
                iv: ivplus - 1,
                y_l,
                y_v,
            }));
        }

        // Bounding values for the other variable across the four nodes
        let ymin = yvec_l[il]
            .min(yvec_l[ilplus])
            .min(yvec_v[iv])
            .min(yvec_v[ivplus]);
        let ymax = yvec_l[il]
            .max(yvec_l[ilplus])
            .max(yvec_v[iv])
            .max(yvec_v[ivplus]);
        if val < ymin || val > ymax {
            return Ok(None);
        }

        // Actually do the "saturation" call using cubic interpolation
        ivplus = ivplus.max(3);
        ilplus = ilplus.max(3);
        let (y_v, y_l) = if main == Param::P {
            let logp = mainval.ln();
            (
                cubic_interp_vec(
                    &self.logp_v,
                    yvec_v,
                    ivplus - 3,
                    ivplus - 2,
                    ivplus - 1,
                    ivplus,
                    logp,
                ),
                cubic_interp_vec(
                    &self.logp_l,
                    yvec_l,
                    ilplus - 3,
                    ilplus - 2,
                    ilplus - 1,
                    ilplus,
                    logp,
                ),
            )
        } else {
            (
                cubic_interp_vec(
                    &self.t_v,
                    yvec_v,
                    ivplus - 3,
                    ivplus - 2,
                    ivplus - 1,
                    ivplus,
                    mainval,
                ),
                cubic_interp_vec(
                    &self.t_l,
                    yvec_l,
                    ilplus - 3,
                    ilplus - 2,
                    ilplus - 1,
                    ilplus,
                    mainval,
                ),
            )
        };

        // `is_in_closed_range(yV, yL, val)` — upstream's helper sorts its
        // bounds, so the branch order does not matter.
        if val < y_v.min(y_l) || val > y_v.max(y_l) {
            Ok(None)
        } else {
            Ok(Some(Inside {
                il: ilplus - 1,
                iv: ivplus - 1,
                y_l,
                y_v,
            }))
        }
    }

    /// `is_inside(iP, p, iT, T, ...)`, the shape `update(PT_INPUTS)` needs.
    pub fn is_inside_pt(&self, p: f64, t: f64) -> Result<Option<(usize, usize)>> {
        Ok(self
            .is_inside(Param::P, p, Param::T, t)?
            .map(|r| (r.il, r.iv)))
    }

    /// `PureFluidSaturationTableData::evaluate`: the two-phase output at
    /// quality `q`, cubic-interpolated along each branch against log(p) (or
    /// against T when the output is pressure). Density and viscosity
    /// interpolate their logs and mix reciprocally; the rest mix linearly.
    pub fn evaluate(
        &self,
        output: Param,
        p_or_t: f64,
        q: f64,
        il: usize,
        iv: usize,
    ) -> Result<f64> {
        let clamp = |mut i: usize| {
            if i <= 2 {
                i = 2;
            } else if i + 1 == self.n {
                i = self.n - 2;
            }
            i
        };
        let (il, iv) = (clamp(il), clamp(iv));
        let logp = p_or_t.ln();
        // The four-point stencils upstream uses, in its order.
        let cl = |x: &Vec<f64>, y: &Vec<f64>, i: usize, v: f64| {
            cubic_interp_vec(x, y, i - 2, i - 1, i, i + 1, v)
        };
        let mix = |v_v: f64, v_l: f64| q * v_v + (1.0 - q) * v_l;
        let checked = |name: &str, v: f64| -> Result<f64> {
            if valid_number(v) {
                Ok(v)
            } else {
                Err(Error::Value(format!("{name} is invalid")))
            }
        };
        Ok(match output {
            Param::P => {
                let logp_v = cl(&self.t_v, &self.logp_v, iv, p_or_t);
                let logp_l = cl(&self.t_l, &self.logp_l, il, p_or_t);
                q * logp_v.exp() + (1.0 - q) * logp_l.exp()
            }
            Param::T => mix(
                cl(&self.logp_v, &self.t_v, iv, logp),
                cl(&self.logp_l, &self.t_l, il, logp),
            ),
            Param::Smolar => mix(
                cl(&self.logp_v, &self.smolar_v, iv, logp),
                cl(&self.logp_l, &self.smolar_l, il, logp),
            ),
            Param::Hmolar => mix(
                cl(&self.logp_v, &self.hmolar_v, iv, logp),
                cl(&self.logp_l, &self.hmolar_l, il, logp),
            ),
            Param::Umolar => mix(
                cl(&self.logp_v, &self.umolar_v, iv, logp),
                cl(&self.logp_l, &self.umolar_l, il, logp),
            ),
            Param::Dmolar => {
                let rho_v = checked(
                    "rhoV",
                    cl(&self.logp_v, &self.logrhomolar_v, iv, logp).exp(),
                )?;
                let rho_l = checked(
                    "rhoL",
                    cl(&self.logp_l, &self.logrhomolar_l, il, logp).exp(),
                )?;
                1.0 / (q / rho_v + (1.0 - q) / rho_l)
            }
            Param::Conductivity => {
                let k_v = checked("kV", cl(&self.logp_v, &self.cond_v, iv, logp))?;
                let k_l = checked("kL", cl(&self.logp_l, &self.cond_l, il, logp))?;
                mix(k_v, k_l)
            }
            Param::Viscosity => {
                let mu_v = checked("muV", cl(&self.logp_v, &self.logvisc_v, iv, logp).exp())?;
                let mu_l = checked("muL", cl(&self.logp_l, &self.logvisc_l, il, logp).exp())?;
                1.0 / (q / mu_v + (1.0 - q) / mu_l)
            }
            Param::Cpmolar => {
                let cp_v = checked("cpV", cl(&self.logp_v, &self.cpmolar_v, iv, logp))?;
                let cp_l = checked("cpL", cl(&self.logp_l, &self.cpmolar_l, il, logp))?;
                mix(cp_v, cp_l)
            }
            Param::Cvmolar => {
                let cv_v = checked("cvV", cl(&self.logp_v, &self.cvmolar_v, iv, logp))?;
                let cv_l = checked("cvL", cl(&self.logp_l, &self.cvmolar_l, il, logp))?;
                mix(cv_v, cv_l)
            }
            Param::SpeedSound => {
                let w_v = checked("wV", cl(&self.logp_v, &self.speed_sound_v, iv, logp))?;
                let w_l = checked("wL", cl(&self.logp_l, &self.speed_sound_l, il, logp))?;
                mix(w_v, w_l)
            }
            other => {
                return Err(Error::Value(format!(
                    "Output key {} is not valid in pure_saturation.evaluate",
                    other.short_name()
                )));
            }
        })
    }
}
