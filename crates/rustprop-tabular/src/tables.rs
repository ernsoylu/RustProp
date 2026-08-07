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

impl SatTable {
    /// `PureFluidSaturationTableData::is_inside(iP, p, iT, T, ...)` — is the
    /// (p, T) state inside the two-phase dome, by the saturation table?
    /// Returns the bracketing indices as upstream does (`iL`/`iV` are set to
    /// `i*plus - 1` when inside).
    pub fn is_inside_pt(&self, p: f64, t: f64) -> Result<Option<(usize, usize)>> {
        // Trivial check on the pressure range
        let pmax = self.p_v[self.n - 1];
        let pmin = self.p_v[0];
        if p > pmax || p < pmin {
            return Ok(None);
        }
        let iv = crate::ttse::bisect_vector(&self.p_v, p)?;
        let il = crate::ttse::bisect_vector(&self.p_l, p)?;
        let mut ivplus = (iv + 1).min(self.n - 1);
        let mut ilplus = (il + 1).min(self.n - 1);

        // Bounding values for T across the four bracketing nodes
        let ymin = self.t_l[il]
            .min(self.t_l[ilplus])
            .min(self.t_v[iv])
            .min(self.t_v[ivplus]);
        let ymax = self.t_l[il]
            .max(self.t_l[ilplus])
            .max(self.t_v[iv])
            .max(self.t_v[ivplus]);
        if t < ymin || t > ymax {
            return Ok(None);
        }

        // Cubic "saturation" call against log(p)
        ivplus = ivplus.max(3);
        ilplus = ilplus.max(3);
        let logp = p.ln();
        let yv = cubic_interp_vec(
            &self.logp_v,
            &self.t_v,
            ivplus - 3,
            ivplus - 2,
            ivplus - 1,
            ivplus,
            logp,
        );
        let yl = cubic_interp_vec(
            &self.logp_l,
            &self.t_l,
            ilplus - 3,
            ilplus - 2,
            ilplus - 1,
            ilplus,
            logp,
        );
        if t < yv.min(yl) || t > yv.max(yl) {
            Ok(None)
        } else {
            Ok(Some((ilplus - 1, ivplus - 1)))
        }
    }
}
