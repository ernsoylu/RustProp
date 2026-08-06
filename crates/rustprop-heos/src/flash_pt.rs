//! PT flash for pure fluids (PLAN.md 4.5) — port of
//! `FlashRoutines::PT_flash`, the superancillary phase determinations
//! (`T_phase_determination_pure_or_pseudopure` /
//! `p_phase_determination_pure_or_pseudopure`, `other == iP`/`iT` paths), and
//! `solver_rho_Tp` with its SRK guess (`solver_rho_Tp_SRK`).
//!

use crate::alpha::HelmholtzEos;
use crate::ancillary;
use crate::saturation::SaturationSuperAncillary;
use crate::solvers::{Resid1D, brent, halley, householder4, solve_cubic};
use rustprop_core::fluid::FluidData;
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};

/// Everything the PT flash needs for one pure fluid.
pub struct PtFlash {
    pub eos: HelmholtzEos,
    /// The saturation superancillary — every pure fluid has one; `None` for
    /// pseudo-pure fluids (upstream skips superancillaries there and the
    /// classic ancillary machinery takes over).
    sat: Option<SaturationSuperAncillary>,
    fluid: &'static FluidData,
    /// Lazily-built caloric superancillaries for the (H,S) flash (upstream
    /// `ensure_caloric_superancillaries` — built on first HS use, cached;
    /// `OnceLock` keeps `PtFlash: Sync` for the facade's per-fluid cache).
    pub(crate) hs_calorics_cell: std::sync::OnceLock<crate::flash_hs::CaloricSa>,
    /// Lazily-built root-enumeration views over the static saturation-density
    /// branches (liquid, vapor) — upstream `get_all_intersections('D', ...)`
    /// walks the same monotonic-interval machinery the caloric SAs use.
    pub(crate) d_approx_cell:
        std::sync::OnceLock<(crate::chebappr::ChebApprox1d, crate::chebappr::ChebApprox1d)>,
    /// Lazily-built melting-curve caloric fits for the HS cascade's leg 4
    /// (upstream `get_melting_caloric_cached`; None = no curve/build failed).
    pub(crate) melting_caloric_cell: std::sync::OnceLock<Option<crate::flash_hs::MeltingCaloric>>,
}

/// Upstream `SolverTPResid`: relative pressure residual with derivatives in
/// density.
struct SolverTpResid<'a> {
    eos: &'a HelmholtzEos,
    t: f64,
    p: f64,
    rhor: f64,
    delta: f64,
}

impl Resid1D for SolverTpResid<'_> {
    fn call(&mut self, rhomolar: f64) -> f64 {
        self.delta = rhomolar / self.rhor; // needed for derivative
        let peos = self.eos.pressure(self.t, rhomolar);
        (peos - self.p) / self.p
    }
    fn deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self
            .eos
            .alphar_all(self.eos.t_reducing / self.t, self.delta);
        self.eos.gas_constant
            * self.t
            * (1.0 + 2.0 * self.delta * d.d10 + self.delta * self.delta * d.d20)
            / self.p
    }
    fn second_deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self
            .eos
            .alphar_all(self.eos.t_reducing / self.t, self.delta);
        self.eos.gas_constant * self.t / self.rhor
            * (2.0 * d.d10 + 4.0 * self.delta * d.d20 + self.delta * self.delta * d.d30)
            / self.p
    }
    fn third_deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self
            .eos
            .alphar_all(self.eos.t_reducing / self.t, self.delta);
        self.eos.gas_constant * self.t / (self.rhor * self.rhor)
            * (6.0 * d.d20 + 6.0 * self.delta * d.d30 + self.delta * self.delta * d.d40)
            / self.p
    }
}

impl PtFlash {
    pub fn new(fluid: &'static FluidData) -> Self {
        let eos = HelmholtzEos::new(fluid);
        let sat = fluid
            .eos
            .superancillary
            .as_ref()
            .map(SaturationSuperAncillary::new);
        PtFlash {
            eos,
            sat,
            fluid,
            hs_calorics_cell: std::sync::OnceLock::new(),
            d_approx_cell: std::sync::OnceLock::new(),
            melting_caloric_cell: std::sync::OnceLock::new(),
        }
    }

    /// The (liquid, vapor) saturation-density branches as root-enumerable
    /// approximations, built once from the static superancillary intervals.
    pub(crate) fn d_approxes(
        &self,
    ) -> &(crate::chebappr::ChebApprox1d, crate::chebappr::ChebApprox1d) {
        self.d_approx_cell.get_or_init(|| {
            let own = |intervals: &[rustprop_core::fluid::ChebyshevInterval]| {
                crate::chebappr::ChebApprox1d::new(
                    intervals
                        .iter()
                        .map(|iv| crate::superancillary::OwnedInterval {
                            xmin: iv.xmin,
                            xmax: iv.xmax,
                            coef: iv.coef.to_vec(),
                        })
                        .collect(),
                )
            };
            let sa = self
                .fluid
                .eos
                .superancillary
                .as_ref()
                .expect("PT flash requires a superancillary fluid");
            (own(sa.rho_l), own(sa.rho_v))
        })
    }

    /// The superancillary — pure fluids only; the pseudo-pure flash paths
    /// never reach the callers of this accessor.
    pub fn sat(&self) -> &SaturationSuperAncillary {
        self.sat
            .as_ref()
            .expect("superancillary paths are pure-fluid only")
    }

    pub fn fluid(&self) -> &'static FluidData {
        self.fluid
    }

    /// Upstream `T_critical()`: for a superancillary fluid, the NUMERICAL
    /// critical temperature (`meta."Tcrittrue / K"`), not the document's
    /// `STATES.critical` — `calc_T_critical` prefers `get_Tcrit_num()`.
    pub fn t_critical(&self) -> f64 {
        match &self.fluid.eos.superancillary {
            Some(sa) => sa.t_crit_num,
            None => self.fluid.states.critical.t,
        }
    }
    /// Upstream `p_critical()`: the superancillary's numerical maximum
    /// pressure (`get_pmax()`), falling back to the document value.
    pub fn p_critical(&self) -> f64 {
        if self.fluid.eos.superancillary.is_some() {
            self.sat().pmax()
        } else {
            self.fluid.states.critical.p
        }
    }
    /// Upstream `rhomolar_critical()`: the superancillary's numerical
    /// critical density (`get_rhocrit_num()`).
    pub fn rhomolar_critical(&self) -> f64 {
        match &self.fluid.eos.superancillary {
            Some(sa) => sa.rho_crit_num,
            None => self.fluid.states.critical.rhomolar,
        }
    }
    /// Upstream `Ttriple()` / `Tmin()`: BOTH resolve to the EOS's
    /// `sat_min_liquid.T` — FluidLibrary overwrites `EOS.Ttriple` and
    /// `limits.Tmin` with it, ignoring the document's chemical `Ttriple`
    /// key (they differ for 27 of the 130 fluids, e.g. CycloPropane
    /// 145.7 K vs 273.0 K).
    pub(crate) fn t_triple(&self) -> f64 {
        self.fluid.eos.sat_min_liquid.t
    }
    fn p_triple(&self) -> f64 {
        self.fluid.eos.sat_min_liquid.p
    }

    /// Upstream `FlashRoutines::PT_flash` for a pure fluid: returns
    /// (rhomolar, phase).
    pub fn pt_flash(&self, t: f64, p: f64) -> Result<(f64, Phase)> {
        let (tc, pc) = (self.t_critical(), self.p_critical());
        // Critical-point short-circuit (upstream #2738)
        if t >= tc * (1.0 - 1e-10)
            && t <= tc * (1.0 + 1e-10)
            && p >= pc * (1.0 - 1e-10)
            && p <= pc * (1.0 + 1e-10)
        {
            return Ok((self.rhomolar_critical(), Phase::CriticalPoint));
        }
        // Upstream phase determinations (both arms): a state below the
        // melting temperature at this pressure is unsupported.
        if let Some(ml) = &self.fluid.ancillaries.melting_line {
            if p > crate::melting::p_min(ml) {
                let tm = crate::melting::t_of_p(ml, p)?;
                if t < tm - 0.001 {
                    return Err(Error::Value(format!(
                        "For now, we don't support T [{t} K] below Tmelt(p) [{tm} K]"
                    )));
                }
            }
        }
        let phase = if t < 0.9 * self.t_triple() + 0.1 * tc {
            self.p_phase_determination_given_t(t, p)?
        } else {
            self.t_phase_determination_given_p(t, p)?
        };
        if phase == Phase::Twophase {
            return Err(Error::Value("twophase not implemented yet".into()));
        }
        let rho = self.solver_rho_tp(t, p, phase)?;
        Ok((rho, phase))
    }

    /// Upstream `T_phase_determination_pure_or_pseudopure`, `other == iP`.
    fn t_phase_determination_given_p(&self, t: f64, p: f64) -> Result<Phase> {
        let (tc, pc) = (self.t_critical(), self.p_critical());
        if t < tc && p > pc {
            Ok(Phase::SupercriticalLiquid)
        } else if (t - tc).abs() < 10.0 * f64::EPSILON {
            if (p - pc).abs() < 10.0 * f64::EPSILON {
                Ok(Phase::CriticalPoint)
            } else if p > pc {
                Ok(Phase::SupercriticalLiquid)
            } else {
                Ok(Phase::SupercriticalGas)
            }
        } else if t < tc {
            // Superancillary phase determination
            let psat = self.sat().qt_flash(t, 1.0)?.p;
            if (psat / p - 1.0).abs() < 1e-6 {
                Err(Error::Value(format!(
                    "Saturation pressure [{psat} Pa] corresponding to T [{t} K] is within 1e-4 % of given p [{p} Pa]"
                )))
            } else if p < psat {
                Ok(Phase::Gas)
            } else {
                Ok(Phase::Liquid)
            }
        } else if t > tc && t > self.t_triple() {
            if p > pc {
                Ok(Phase::Supercritical)
            } else {
                Ok(Phase::SupercriticalGas)
            }
        } else {
            Err(Error::Value("phase determination fell through".into()))
        }
    }

    /// Upstream `p_phase_determination_pure_or_pseudopure`, `other == iT`.
    fn p_phase_determination_given_t(&self, t: f64, p: f64) -> Result<Phase> {
        let psat_max = self.p_critical(); // calc_pmax_sat, pure fluid
        if p > psat_max {
            if t > self.t_critical() {
                Ok(Phase::Supercritical)
            } else {
                Ok(Phase::SupercriticalLiquid)
            }
        } else if p >= self.p_triple() * 0.9999 && p <= psat_max {
            if p > self.sat().pmax() {
                return Err(Error::Value(format!(
                    "Pressure to PQ_flash [{p:.8e} Pa] may not be above the numerical critical point of {:.15} Pa",
                    self.sat().pmax()
                )));
            }
            let tsat = self.sat().pq_flash(p, 0.0)?.t;
            if t < tsat - 100.0 * f64::EPSILON {
                Ok(Phase::Liquid)
            } else if t > tsat + 100.0 * f64::EPSILON {
                Ok(Phase::Gas)
            } else {
                Ok(Phase::Twophase)
            }
        } else if p < self.p_triple() * 0.9999 {
            // Below the triple-point pressure: gas, if the temperature allows
            if t < self.t_triple() {
                Err(Error::NotImplemented(format!(
                    "For now, we don't support p [{p} Pa] below ptriple [{} Pa] when T [{t}] is less than Tmin [{}]",
                    self.p_triple(),
                    self.t_triple()
                )))
            } else {
                Ok(Phase::Gas)
            }
        } else {
            Err(Error::Value(format!(
                "The pressure [{p} Pa] cannot be used in p_phase_determination"
            )))
        }
    }

    /// dp/drho|T (upstream `first_partial_deriv(iP, iDmolar, iT)`).
    fn dpdrho_t(&self, t: f64, rhomolar: f64) -> f64 {
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let d = self.eos.alphar_all(self.eos.t_reducing / t, delta);
        self.eos.gas_constant * t * (1.0 + 2.0 * delta * d.d10 + delta * delta * d.d20)
    }

    /// d2p/drho2|T (upstream `second_partial_deriv(iP, iDmolar, iT, ...)`).
    fn d2pdrho2_t(&self, t: f64, rhomolar: f64) -> f64 {
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let d = self.eos.alphar_all(self.eos.t_reducing / t, delta);
        self.eos.gas_constant * t / self.eos.rhomolar_reducing
            * (2.0 * d.d10 + 4.0 * delta * d.d20 + delta * delta * d.d30)
    }

    /// Upstream `solver_rho_Tp_SRK` for a single component.
    #[allow(clippy::excessive_precision)] // SRK constants verbatim from upstream
    fn solver_rho_tp_srk(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let r_u = self.eos.gas_constant;
        let tci = self.fluid.eos.reducing.t;
        let pci = self.fluid.eos.reducing.p;
        let acentric = self.fluid.eos.acentric;
        let m_i = 0.480 + 1.574 * acentric - 0.176 * acentric.powf(2.0);
        let b = 0.08664034999649577215890158147700 * r_u * tci / pci;
        let a = 0.42748023335403414043900347952220 * (r_u * tci).powf(2.0) / pci
            * (1.0 + m_i * (1.0 - (t / tci).sqrt())).powf(2.0);

        let big_a = a * p / (r_u * t).powf(2.0);
        let big_b = b * p / (r_u * t);

        let (nsolns, z0, z1, z2) =
            solve_cubic(1.0, -1.0, big_a - big_b - big_b * big_b, -big_a * big_b);
        if nsolns == 1 {
            return Ok(p / (z0 * r_u * t));
        }
        let rho0 = p / (z0 * r_u * t);
        let rho1 = p / (z1 * r_u * t);
        let rho2 = p / (z2 * r_u * t);
        if rho0 > 0.0 && rho1 <= 0.0 && rho2 <= 0.0 {
            return Ok(rho0);
        }
        if rho0 <= 0.0 && rho1 > 0.0 && rho2 <= 0.0 {
            return Ok(rho1);
        }
        if rho0 <= 0.0 && rho1 <= 0.0 && rho2 > 0.0 {
            return Ok(rho2);
        }
        match phase {
            Phase::Liquid | Phase::SupercriticalLiquid => Ok(rho0.max(rho1).max(rho2)),
            Phase::Gas | Phase::SupercriticalGas | Phase::Supercritical => {
                Ok(rho0.min(rho1).min(rho2))
            }
            _ => Err(Error::Value("Bad phase to solver_rho_Tp_SRK".into())),
        }
    }

    /// Upstream `solver_rho_Tp` (no caller-provided guess; pure fluid).
    /// Upstream `solver_rho_Tp(T, p, rhomolar_guess)` with a PROVIDED guess
    /// (the `update_TP_guessrho` path): Householder4 from the guess with the
    /// supercritical Brent fallbacks; no phase-specific stability retries
    /// (the caller's phase is not imposed on this path).
    pub fn solver_rho_tp_guessed(&self, t: f64, p: f64, rho_guess: f64) -> Result<f64> {
        let mut resid = SolverTpResid {
            eos: &self.eos,
            t,
            p,
            rhor: self.eos.rhomolar_reducing,
            delta: f64::NAN,
        };
        match householder4(&mut resid, rho_guess, 1e-8, 20) {
            Ok(rho) if rho.is_finite() && rho >= 0.0 => Ok(rho),
            _ => {
                if t > self.t_critical() {
                    let f = |rho: f64| (self.eos.pressure(t, rho) - p) / p;
                    if let Ok(rho) = brent(
                        f,
                        1e-10,
                        5.0 * self.eos.rhomolar_reducing,
                        f64::EPSILON,
                        1e-8,
                        100,
                    ) {
                        return Ok(rho);
                    }
                    return householder4(&mut resid, 3.0 * self.eos.rhomolar_reducing, 1e-8, 100);
                }
                Err(Error::Value(format!(
                    "solver_rho_Tp was unable to find a solution for T={t:10}, p={p:10}, with guess value {rho_guess:10}"
                )))
            }
        }
    }

    pub fn solver_rho_tp(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        let mut resid = SolverTpResid {
            eos: &self.eos,
            t,
            p,
            rhor: self.eos.rhomolar_reducing,
            delta: f64::NAN,
        };

        let mut rhomolar_guess = self.solver_rho_tp_srk(t, p, phase)?;
        if matches!(
            phase,
            Phase::Gas | Phase::SupercriticalGas | Phase::Supercritical
        ) {
            if rhomolar_guess < 0.0 || !rhomolar_guess.is_finite() {
                rhomolar_guess = p / (self.eos.gas_constant * t);
            }
        } else if phase == Phase::Liquid {
            let rho_l_anc = ancillary::evaluate(&self.fluid.ancillaries.rho_l, t);
            // First Halley from the saturated-liquid ancillary, with the
            // upstream stability checks; then a bounded Brent fallback.
            let attempt = halley(&mut resid, rho_l_anc, 1e-8, 100).and_then(|rho| {
                if !rho.is_finite() || self.dpdrho_t(t, rho) < 0.0 || self.d2pdrho2_t(t, rho) < 0.0
                {
                    Err(Error::Value("Liquid density is invalid".into()))
                } else {
                    Ok(rho)
                }
            });
            return match attempt {
                Ok(rho) => Ok(rho),
                Err(_) => {
                    let rho = brent(
                        |x| resid.call(x),
                        rho_l_anc * 0.9,
                        rho_l_anc * 1.3,
                        f64::EPSILON,
                        1e-8,
                        100,
                    )?;
                    if !rho.is_finite() {
                        return Err(Error::Value("Liquid density is invalid".into()));
                    }
                    Ok(rho)
                }
            };
        } else if phase == Phase::SupercriticalLiquid {
            let rho_l_anc = ancillary::evaluate(&self.fluid.ancillaries.rho_l, t);
            let rho_l_triple_anc =
                ancillary::evaluate(&self.fluid.ancillaries.rho_l, self.t_triple());
            let first = brent(
                |x| resid.call(x),
                rho_l_anc * 0.99,
                self.rhomolar_critical() * 4.0,
                f64::EPSILON,
                1e-8,
                100,
            );
            match first {
                Ok(rho) if rho.is_finite() => return Ok(rho),
                _ => {
                    let second = brent(
                        |x| resid.call(x),
                        rho_l_anc * 0.99,
                        rho_l_triple_anc * 1.1,
                        f64::EPSILON,
                        1e-8,
                        100,
                    );
                    match second {
                        Ok(rho) if rho.is_finite() => return Ok(rho),
                        _ => return self.dense_branch_fallback(t, p),
                    }
                }
            }
        }

        // Gas-like phases: Householder4 from the SRK/ideal-gas guess.
        let attempt = householder4(&mut resid, rhomolar_guess, 1e-8, 20).and_then(|rho| {
            if !rho.is_finite() || rho < 0.0 {
                return Err(Error::Value("invalid density".into()));
            }
            if phase == Phase::Gas {
                let dpdrho = self.dpdrho_t(t, rho);
                let d2pdrho2 = self.d2pdrho2_t(t, rho);
                if dpdrho < 0.0 || d2pdrho2 > 0.0 {
                    // Retry with a tiny density to land on the gas branch
                    return householder4(&mut resid, 1e-6, 1e-8, 100);
                }
            }
            Ok(rho)
        });
        match attempt {
            Ok(rho) => Ok(rho),
            Err(e) => {
                if matches!(phase, Phase::Supercritical | Phase::SupercriticalGas) {
                    return brent(
                        |x| resid.call(x),
                        1e-10,
                        3.0 * self.eos.rhomolar_reducing,
                        f64::EPSILON,
                        1e-8,
                        100,
                    );
                }
                if t > self.t_critical() {
                    return match brent(
                        |x| resid.call(x),
                        1e-10,
                        5.0 * self.eos.rhomolar_reducing,
                        f64::EPSILON,
                        1e-8,
                        100,
                    ) {
                        Ok(rho) => Ok(rho),
                        Err(_) => {
                            householder4(&mut resid, 3.0 * self.eos.rhomolar_reducing, 1e-8, 100)
                        }
                    };
                }
                Err(Error::Value(format!(
                    "solver_rho_Tp was unable to find a solution for T={t:.10e}, p={p:.10e}, with guess value {rhomolar_guess:.10e} with error: {e}"
                )))
            }
        }
    }

    /// Upstream's robust dense-branch fallback for supercritical liquid
    /// (bracket walk + TOMS748; bisection stands in for TOMS748 as in the
    /// saturation module).
    fn dense_branch_fallback(&self, t: f64, p_target: f64) -> Result<f64> {
        let mut rho_hi = 5.0 * self.eos.rhomolar_reducing;
        for _ in 0..60 {
            if self.eos.pressure(t, rho_hi) >= p_target {
                break;
            }
            rho_hi *= 1.5;
        }
        let mut rho_lo = rho_hi;
        let mut bracketed = false;
        for _ in 0..400 {
            rho_lo /= 1.05;
            if rho_lo <= 1e-10 {
                break;
            }
            if self.eos.pressure(t, rho_lo) < p_target {
                bracketed = true;
                break;
            }
        }
        if bracketed {
            let f = |rho: f64| self.eos.pressure(t, rho) - p_target;
            let (fa, fb) = (f(rho_lo), f(rho_hi));
            if fa.is_finite() && fb.is_finite() && fa * fb < 0.0 {
                let mut a = rho_lo;
                let mut b = rho_hi;
                let (mut fa, _fb) = (fa, fb);
                for _ in 0..200 {
                    let m = 0.5 * (a + b);
                    if m <= a || m >= b {
                        break;
                    }
                    let fm = f(m);
                    if fm == 0.0 {
                        return Ok(m);
                    }
                    if (fa < 0.0) == (fm < 0.0) {
                        a = m;
                        fa = fm;
                    } else {
                        b = m;
                    }
                }
                let rho = 0.5 * (a + b);
                if rho.is_finite() && rho > 0.0 {
                    return Ok(rho);
                }
            }
        }
        Err(Error::Value(format!(
            "solver_rho_Tp could not find a supercritical_liquid density for T={t}, p={p_target}"
        )))
    }
}
