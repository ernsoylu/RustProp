//! Remaining pure-fluid flash pairs (PLAN.md 4.6): general-quality (T,Q) and
//! (P,Q) states with upstream's two-phase mixing, the (D,T) flash
//! (`DHSU_T_flash` -> superancillary `T_phase_determination`), and the
//! (H,P)/(P,S) flashes (`HSU_P_flash`: superancillary `p_phase_determination`
//! plus the bracketed single-phase solve in T).
//!
//! Numerical notes, logged in PLAN.md:
//! - upstream resolves the single-phase (P,X) temperature with TOMS748 at a
//!   deliberate 30-bit (~1e-9 relative) tolerance, then re-evaluates at the
//!   bracket midpoint; the ported `solvers::toms748_solve` runs at the same
//!   30 bits with the same midpoint re-evaluation, and the inner density
//!   solve carries the previous probe's density warm exactly as upstream
//!   does — golden agreement is bounded by that shared tolerance, hence the
//!   1e-8 policy;
//! - the no-bracket derivative path and the 2-D Newton fallback are unported
//!   (loud error) until a state needs them;
//! - two-phase mixing is upstream's exact `Q*V + (1-Q)*L` with the
//!   DBL_EPSILON endpoint shortcuts.

use crate::alpha::{DerivsMemo, HelmholtzDerivs, HelmholtzEos};
use crate::flash_pt::PtFlash;
use crate::solvers::Resid1D;
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};

/// `LDBL_EPSILON` of the upstream build (x86-64 80-bit long double,
/// 2^-63) — upstream passes it as the Brent `macheps` in
/// `solver_for_rho_given_T_oneof_HSU`.
const LDBL_EPSILON: f64 = 1.084_202_172_485_504_4e-19;

/// Residual s(T, rho) - target with the analytic density derivatives
/// (upstream `solver_resid` in `solver_for_rho_given_T_oneof_HSU`, whose
/// `first_partial_deriv(iSmolar, iDmolar, iT)` values these closed forms
/// reproduce).
struct CaloricTResid<'a> {
    eos: &'a HelmholtzEos,
    t: f64,
    target: f64,
    key: CaloricKey,
    memo: DerivsMemo,
}

impl CaloricTResid<'_> {
    /// The alphar matrix at (tau, delta), computed once per point.
    fn ar(&mut self, tau: f64, delta: f64) -> HelmholtzDerivs {
        let eos = self.eos;
        self.memo
            .get_or_compute(tau, delta, |tau, delta| eos.alphar_all(tau, delta))
    }
}

impl Resid1D for CaloricTResid<'_> {
    fn call(&mut self, rhomolar: f64) -> f64 {
        // `calc_smolar/hmolar/umolar_nocache` off the shared alphar matrix —
        // same tau/delta expressions and arithmetic as the `HelmholtzEos`
        // methods; alpha0 has no second same-point consumer, so it is
        // computed here exactly as before.
        let tau = self.eos.t_reducing / self.t;
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let residual = self.ar(tau, delta);
        let ideal = self.eos.alpha0_all(tau, delta);
        let r = self.eos.gas_constant;
        let v = match self.key {
            CaloricKey::Smolar => r * (tau * (ideal.d01 + residual.d01) - ideal.d00 - residual.d00),
            CaloricKey::Hmolar => {
                r * self.t * (1.0 + tau * (ideal.d01 + residual.d01) + delta * residual.d10)
            }
            CaloricKey::Umolar => r * self.t * tau * (ideal.d01 + residual.d01),
        };
        v - self.target
    }
    /// d(other)/drho|T (upstream `first_partial_deriv(other, iDmolar, iT)`):
    /// S: R*(tau*d11 - 1/delta - d10)/rho_r
    /// H: R*T*(tau*d11 + d10 + delta*d20)/rho_r
    /// U: R*T*tau*d11/rho_r
    fn deriv(&mut self, rhomolar: f64) -> f64 {
        let tau = self.eos.t_reducing / self.t;
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let d = self.ar(tau, delta);
        let r = self.eos.gas_constant;
        match self.key {
            CaloricKey::Smolar => {
                r * (tau * d.d11 - 1.0 / delta - d.d10) / self.eos.rhomolar_reducing
            }
            CaloricKey::Hmolar => {
                r * self.t * (tau * d.d11 + d.d10 + delta * d.d20) / self.eos.rhomolar_reducing
            }
            CaloricKey::Umolar => r * self.t * tau * d.d11 / self.eos.rhomolar_reducing,
        }
    }
    /// d2(other)/drho2|T:
    /// S: R*(tau*d21 + 1/delta^2 - d20)/rho_r^2
    /// H: R*T*(tau*d21 + 2*d20 + delta*d30)/rho_r^2
    /// U: R*T*tau*d21/rho_r^2
    fn second_deriv(&mut self, rhomolar: f64) -> f64 {
        let tau = self.eos.t_reducing / self.t;
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let d = self.ar(tau, delta);
        let r = self.eos.gas_constant;
        let rho_r2 = self.eos.rhomolar_reducing * self.eos.rhomolar_reducing;
        match self.key {
            CaloricKey::Smolar => r * (tau * d.d21 + 1.0 / (delta * delta) - d.d20) / rho_r2,
            CaloricKey::Hmolar => r * self.t * (tau * d.d21 + 2.0 * d.d20 + delta * d.d30) / rho_r2,
            CaloricKey::Umolar => r * self.t * tau * d.d21 / rho_r2,
        }
    }
    fn third_deriv(&mut self, _rhomolar: f64) -> f64 {
        unreachable!("Halley does not use the third derivative")
    }
}

/// Upstream `HSU_P_flash_singlephase_Brent::solver_resid` — the single-phase
/// (p, X) temperature residual `keyed_output(other) - value`.
///
/// The functor carries the last three converged densities so a probe can seed
/// its inner (T, p) density solve from the previous probe's answer
/// (`update_TP_guessrho` -> `solver_rho_Tp(T, p, guess)`) instead of re-running
/// the cold SRK-seeded `update(PT_INPUTS, p, T)`. Upstream's literal gate:
///
/// ```text
/// if (iter < 2 || std::abs(rhomolar1/rhomolar0 - 1) > 0.05 || force_robust_density)
///     HEOS->update(PT_INPUTS, p, T);        // cold
/// else
///     HEOS->update_TP_guessrho(T, p, rhomolar);  // warm
/// ```
///
/// `force_robust_density = (p > p_crit)`: above the critical pressure the
/// carried-rho Newton can hop branches inside the van der Waals loop (below
/// T_crit) or blow up where dp/drho -> 0 (just above it), which makes the outer
/// residual non-monotone; upstream forces every supercritical-pressure probe
/// cold. Upstream's `eos0`/`eos1` companions are not carried: their only
/// consumer is the "above/below the maximum value" diagnostic of the unported
/// derivative-path tail, whose message is a documented deviation.
///
/// NOT ported: upstream's `PXFLASH_DIRECT_EOS` cache-bypass (default ON, but
/// gated on `is_pure()`, so pseudo-pure fluids never see it), which runs a warm
/// probe's density solve as a Householder3 straight off
/// `residual_helmholtz->all` at a 1e-12 relative pressure residual. Upstream
/// calls it "bit-equivalent within ULP" to the cached path and wraps it in a
/// `catch (...)` that falls through to exactly the `update_TP_guessrho` branch
/// ported here — so this is upstream's own fallback, not an invented one. It is
/// a cache-avoidance trick aimed at a cache this port does not have (each
/// residual owns a `DerivsMemo` instead); the price is ULP-scale disagreement
/// on warm probes, measured at a 2.0e-16 median over the (P, caloric) goldens.
struct PxResid<'a> {
    flash: &'a PtFlash,
    p: f64,
    value: f64,
    key: CaloricKey,
    /// The bracket's working phase — upstream's `specify_phase(iphase_liquid |
    /// iphase_gas)` in the functor's constructor.
    phase: Phase,
    /// The phase the state is actually holding (upstream `_phase`), which the
    /// guessed density solve consults where nothing is imposed.
    live_phase: Phase,
    iter: u32,
    rhomolar: f64,
    rhomolar0: f64,
    rhomolar1: f64,
    force_robust_density: bool,
}

impl PxResid<'_> {
    fn call(&mut self, t: f64) -> Result<f64> {
        let rho = if self.iter < 2
            || (self.rhomolar1 / self.rhomolar0 - 1.0).abs() > 0.05
            || self.force_robust_density
        {
            let (rho, live) = self.flash.px_probe_rho(t, self.p, self.phase)?;
            self.live_phase = live;
            rho
        } else {
            self.flash
                .solver_rho_tp_guessed(t, self.p, self.rhomolar, self.live_phase)?
        };
        let eos = self.flash.px_value(self.key, t, rho);
        self.rhomolar = rho;
        let r = eos - self.value;
        if self.iter == 0 {
            self.rhomolar0 = rho;
        } else if self.iter == 1 {
            self.rhomolar1 = rho;
        } else {
            self.rhomolar0 = self.rhomolar1;
            self.rhomolar1 = rho;
        }
        self.iter += 1;
        Ok(r)
    }
}

/// A fully-determined thermodynamic state for one pure fluid.
#[derive(Debug, Clone, Copy)]
pub enum HeosState {
    SinglePhase {
        t: f64,
        p: f64,
        rhomolar: f64,
        phase: Phase,
        /// Upstream `_Q` sentinel: -1 for the flashes that set it so
        /// (PT/PX/DT and the legacy HS path); the superancillary HS cascade
        /// leaves upstream's 10000, observable through `PropsSI("Q")`.
        q: f64,
    },
    TwoPhase {
        t: f64,
        p: f64,
        rhomolar: f64,
        q: f64,
        rho_l: f64,
        rho_v: f64,
        /// Per-branch saturation temperatures — equal to `t` for pure
        /// fluids; a pseudo-pure PQ state carries the temperature GLIDE
        /// (T_bubble != T_dew), and caloric mixes evaluate each branch at
        /// its own temperature exactly as upstream's SatL/SatV sub-states.
        t_l: f64,
        t_v: f64,
    },
}

impl HeosState {
    pub fn t(&self) -> f64 {
        match self {
            HeosState::SinglePhase { t, .. } | HeosState::TwoPhase { t, .. } => *t,
        }
    }
    pub fn p(&self) -> f64 {
        match self {
            HeosState::SinglePhase { p, .. } | HeosState::TwoPhase { p, .. } => *p,
        }
    }
    pub fn rhomolar(&self) -> f64 {
        match self {
            HeosState::SinglePhase { rhomolar, .. } | HeosState::TwoPhase { rhomolar, .. } => {
                *rhomolar
            }
        }
    }
    /// Upstream `_Q`: the vapor quality for two-phase states, the
    /// flash-specific sentinel for single-phase ones.
    pub fn q(&self) -> f64 {
        match self {
            HeosState::SinglePhase { q, .. } | HeosState::TwoPhase { q, .. } => *q,
        }
    }
}

/// Upstream two-phase mixing (`calc_hmolar` etc.): endpoint shortcuts at
/// DBL_EPSILON, otherwise `Q*V + (1-Q)*L`.
fn mix_two_phase(q: f64, liquid: f64, vapor: f64) -> f64 {
    if q.abs() < f64::EPSILON {
        liquid
    } else if (q - 1.0).abs() < f64::EPSILON {
        vapor
    } else {
        q * vapor + (1.0 - q) * liquid
    }
}

impl PtFlash {
    /// Molar enthalpy of a state [J/mol].
    pub fn state_hmolar(&self, s: &HeosState) -> f64 {
        match s {
            HeosState::SinglePhase { t, rhomolar, .. } => self.eos.hmolar(*t, *rhomolar),
            HeosState::TwoPhase {
                q,
                rho_l,
                rho_v,
                t_l,
                t_v,
                ..
            } => mix_two_phase(
                *q,
                self.eos.hmolar(*t_l, *rho_l),
                self.eos.hmolar(*t_v, *rho_v),
            ),
        }
    }
    /// Molar entropy of a state [J/mol/K].
    pub fn state_smolar(&self, s: &HeosState) -> f64 {
        match s {
            HeosState::SinglePhase { t, rhomolar, .. } => self.eos.smolar(*t, *rhomolar),
            HeosState::TwoPhase {
                q,
                rho_l,
                rho_v,
                t_l,
                t_v,
                ..
            } => mix_two_phase(
                *q,
                self.eos.smolar(*t_l, *rho_l),
                self.eos.smolar(*t_v, *rho_v),
            ),
        }
    }
    /// Molar internal energy of a state [J/mol].
    pub fn state_umolar(&self, s: &HeosState) -> f64 {
        match s {
            HeosState::SinglePhase { t, rhomolar, .. } => self.eos.umolar(*t, *rhomolar),
            HeosState::TwoPhase {
                q,
                rho_l,
                rho_v,
                t_l,
                t_v,
                ..
            } => mix_two_phase(
                *q,
                self.eos.umolar(*t_l, *rho_l),
                self.eos.umolar(*t_v, *rho_v),
            ),
        }
    }

    /// General-quality (T,Q) state (superancillary `QT_flash`; the
    /// pseudo-pure branch uses the explicit pL/pV ancillaries).
    pub fn qt_state(&self, t: f64, q: f64) -> Result<HeosState> {
        if self.fluid().eos.pseudo_pure {
            return self.qt_state_pseudo_pure(t, q);
        }
        let sat = self.sat().qt_flash(t, q)?;
        Ok(HeosState::TwoPhase {
            t: sat.t,
            p: sat.p,
            rhomolar: sat.rhomolar,
            q,
            rho_l: sat.rho_l,
            rho_v: sat.rho_v,

            t_l: sat.t,

            t_v: sat.t,
        })
    }

    /// Upstream `QT_flash`'s pseudo-pure branch: the range guard on
    /// [Tmin_sat - 0.1, Tmax_sat], then p straight from the pL (Q=0) or pV
    /// (Q=1) ancillary and the density from a PT solve seeded by the
    /// rhoL/rhoV ancillary. Fractional quality is undefined upstream.
    fn qt_state_pseudo_pure(&self, t: f64, q: f64) -> Result<HeosState> {
        let fluid = self.fluid();
        let anc = &fluid.ancillaries;
        // Upstream `calc_Tmax_sat` (pseudo-pure: max_sat_T.T) and
        // `calc_Tmin_sat` (max of the two sat-min temperatures).
        let tmax_sat = fluid
            .eos
            .max_sat_t
            .as_ref()
            .map_or_else(|| self.t_critical(), |sp| sp.t)
            + 1e-13;
        let tmin_sat = fluid.eos.sat_min_liquid.t.max(fluid.eos.sat_min_vapor.t) - 1e-13;
        if !(tmin_sat - 0.1..=tmax_sat).contains(&t) {
            return Err(Error::Value(format!(
                "Temperature to QT_flash [{t:.8} K] must be in range [{:.8} K, {:.8} K]",
                tmin_sat - 0.1,
                tmax_sat
            )));
        }
        let (p, rho) = if q.abs() < f64::EPSILON {
            let p = crate::ancillary::evaluate(&anc.p_s, t);
            let rho_anc = crate::ancillary::evaluate(&anc.rho_l, t);
            (p, self.solver_rho_tp_guessed(t, p, rho_anc, Phase::Liquid)?)
        } else if (q - 1.0).abs() < f64::EPSILON {
            let pv = anc.p_v_split.as_ref().unwrap_or(&anc.p_s);
            let p = crate::ancillary::evaluate(pv, t);
            let rho_anc = crate::ancillary::evaluate(&anc.rho_v, t);
            (p, self.solver_rho_tp_guessed(t, p, rho_anc, Phase::Gas)?)
        } else {
            return Err(Error::Value(
                "For pseudo-pure fluid, quality must be equal to 0 or 1.  Two-phase quality is not defined"
                    .into(),
            ));
        };
        // Upstream commits the solved branch density as the state density;
        // the other branch is never solved (the exact-0/1 quality zeroes it
        // out of every mix).
        Ok(HeosState::TwoPhase {
            t,
            p,
            rhomolar: rho,
            q,
            rho_l: rho,
            rho_v: rho,
            t_l: t,
            t_v: t,
        })
    }

    /// General-quality (P,Q) state (superancillary `PQ_flash`; the
    /// pseudo-pure branch inverts the pL/pV ancillaries per branch).
    pub fn pq_state(&self, p: f64, q: f64) -> Result<HeosState> {
        if self.fluid().eos.pseudo_pure {
            return self.pq_state_pseudo_pure(p, q);
        }
        let sat = self.sat().pq_flash(p, q)?;
        Ok(HeosState::TwoPhase {
            t: sat.t,
            p: sat.p,
            rhomolar: sat.rhomolar,
            q,
            rho_l: sat.rho_l,
            rho_v: sat.rho_v,

            t_l: sat.t,

            t_v: sat.t,
        })
    }

    /// Upstream `PQ_flash`'s pseudo-pure branch: invert pL and pV for the
    /// bubble/dew temperatures (the GLIDE — T_L != T_V), solve each branch
    /// density from its own ancillary-seeded PT solve, and quality-mix the
    /// outputs. The state carries both branch temperatures so caloric mixes
    /// evaluate each phase at its own T, exactly as upstream's SatL/SatV.
    fn pq_state_pseudo_pure(&self, p: f64, q: f64) -> Result<HeosState> {
        let anc = &self.fluid().ancillaries;
        let pv_anc = anc.p_v_split.as_ref().unwrap_or(&anc.p_s);
        let t_l = crate::ancillary::invert(&anc.p_s, p)?;
        let t_v = crate::ancillary::invert(pv_anc, p)?;
        let rho_l_anc = crate::ancillary::evaluate(&anc.rho_l, t_l);
        let rho_v_anc = crate::ancillary::evaluate(&anc.rho_v, t_v);
        let rho_l = self.solver_rho_tp_guessed(t_l, p, rho_l_anc, Phase::Liquid)?;
        let rho_v = self.solver_rho_tp_guessed(t_v, p, rho_v_anc, Phase::Gas)?;
        Ok(HeosState::TwoPhase {
            t: q * t_v + (1.0 - q) * t_l,
            p,
            rhomolar: 1.0 / (q / rho_v + (1.0 - q) / rho_l),
            q,
            rho_l,
            rho_v,
            t_l,
            t_v,
        })
    }

    /// (Dmolar, T) flash — upstream `DHSU_T_flash(iDmolar)`, which for a pure
    /// fluid is exactly `T_phase_determination_pure_or_pseudopure(iDmolar)`.
    pub fn dmolar_t_state(&self, rhomolar: f64, t: f64) -> Result<HeosState> {
        let tc = self.t_critical();
        let rhoc = self.rhomolar_critical();
        if (t - tc).abs() < 10.0 * f64::EPSILON {
            // Exactly at Tcrit
            let phase = if (rhomolar - rhoc).abs() < 10.0 * f64::EPSILON {
                Phase::CriticalPoint
            } else if rhomolar > rhoc {
                Phase::SupercriticalLiquid
            } else {
                Phase::SupercriticalGas
            };
            return Ok(HeosState::SinglePhase {
                t,
                p: self.eos.pressure(t, rhomolar),
                rhomolar,
                phase,
                q: -1.0,
            });
        }
        if t < tc {
            // Superancillary phase determination
            let sat = self.sat().qt_flash(t, 0.0)?;
            let q = (1.0 / rhomolar - 1.0 / sat.rho_l) / (1.0 / sat.rho_v - 1.0 / sat.rho_l);
            // The single-phase exits relabel by the FINAL pressure — the
            // shipped wheel reclassifies a compressed liquid with p > pc as
            // supercritical_liquid (oracle-probed: Methane 26.77 kmol/m^3 at
            // 119.3 K, p = 20.4 MPa -> phase 3, while p < pc stays liquid).
            if q <= 0.0 || q >= 1.0 {
                let p = self.eos.pressure(t, rhomolar);
                Ok(HeosState::SinglePhase {
                    t,
                    p,
                    rhomolar,
                    phase: self.recalculated_singlephase_phase(t, p, rhomolar),
                    q: -1.0,
                })
            } else {
                Ok(HeosState::TwoPhase {
                    t,
                    p: sat.p,
                    rhomolar,
                    q,
                    rho_l: sat.rho_l,
                    rho_v: sat.rho_v,

                    t_l: t,

                    t_v: t,
                })
            }
        } else if t > tc && t > self.t_triple() {
            // Above Tc the label comes from the pressure quadrant, not the
            // density (wheel: dense DT at 250 K -> supercritical when
            // p > pc, supercritical_gas when p < pc — never
            // supercritical_liquid).
            let p = self.eos.pressure(t, rhomolar);
            Ok(HeosState::SinglePhase {
                t,
                p,
                rhomolar,
                phase: self.recalculated_singlephase_phase(t, p, rhomolar),
                q: -1.0,
            })
        } else {
            Err(Error::Value(
                "temperature is out of range in dmolar_t_state".into(),
            ))
        }
    }

    /// (Smolar, T) flash — upstream `DHSU_T_flash(iSmolar)`.
    pub fn smolar_t_state(&self, smolar: f64, t: f64) -> Result<HeosState> {
        self.caloric_t_state(smolar, t, CaloricKey::Smolar)
    }

    /// (Hmolar, T) flash — upstream `DHSU_T_flash(iHmolar)`.
    pub fn hmolar_t_state(&self, hmolar: f64, t: f64) -> Result<HeosState> {
        self.caloric_t_state(hmolar, t, CaloricKey::Hmolar)
    }

    /// (T, Umolar) flash — upstream `DHSU_T_flash(iUmolar)`.
    pub fn umolar_t_state(&self, umolar: f64, t: f64) -> Result<HeosState> {
        self.caloric_t_state(umolar, t, CaloricKey::Umolar)
    }

    /// Upstream `DHSU_T_flash(other one of iSmolar/iHmolar/iUmolar)`:
    /// superancillary phase determination
    /// (`T_phase_determination_pure_or_pseudopure`), then
    /// `solver_for_rho_given_T_oneof_HSU` for the single-phase branches.
    fn caloric_t_state(&self, value: f64, t: f64, key: CaloricKey) -> Result<HeosState> {
        let tc = self.t_critical();
        if (t - tc).abs() < 10.0 * f64::EPSILON {
            // Upstream supports only iDmolar/iP at exactly Tcrit.
            return Err(Error::Value(
                "T=Tcrit; invalid input for other to T_phase_determination_pure_or_pseudopure"
                    .into(),
            ));
        }
        if t < tc {
            // Superancillary phase determination
            let sat = self.sat().qt_flash(t, 0.0)?;
            let y_l = self.px_value(key, t, sat.rho_l);
            let y_v = self.px_value(key, t, sat.rho_v);
            let q = (value - y_l) / (y_v - y_l);
            if q < 0.0 {
                self.rho_from_caloric_t(t, value, key, Phase::Liquid, sat.rho_l)
            } else if q > 1.0 {
                self.rho_from_caloric_t(t, value, key, Phase::Gas, sat.rho_v)
            } else {
                Ok(HeosState::TwoPhase {
                    t,
                    p: sat.p,
                    rhomolar: 1.0 / (q / sat.rho_v + (1.0 - q) / sat.rho_l),
                    q,
                    rho_l: sat.rho_l,
                    rho_v: sat.rho_v,

                    t_l: t,

                    t_v: t,
                })
            }
        } else if t > tc && t > self.t_triple() {
            self.rho_from_caloric_t_supercritical(t, value, key)
        } else {
            Err(Error::Value(
                "temperature is out of range in caloric_t_state".into(),
            ))
        }
    }

    /// Subcritical single-phase branches of
    /// `solver_for_rho_given_T_oneof_HSU`. `rho_anc` is the superancillary
    /// saturation density of the branch (upstream's `_rhoLanc`/`_rhoVanc`
    /// set by the phase determination).
    fn rho_from_caloric_t(
        &self,
        t: f64,
        value: f64,
        key: CaloricKey,
        phase: Phase,
        rho_anc: f64,
    ) -> Result<HeosState> {
        let mut resid = CaloricTResid {
            eos: &self.eos,
            t,
            target: value,
            key,
            memo: DerivsMemo::default(),
        };
        let rho = match phase {
            Phase::Liquid => {
                let rhomelt = self.fluid().states.triple_liquid.rhomolar;
                let ymelt = self.px_value(key, t, rhomelt);
                let y_l = self.px_value(key, t, rho_anc);
                let guess = (rhomelt - rho_anc) / (ymelt - y_l) * (value - y_l) + rho_anc;
                match crate::solvers::halley(&mut resid, guess, 1e-8, 100) {
                    Ok(rho) => rho,
                    Err(_) => crate::solvers::secant(
                        |rho| self.px_value(key, t, rho) - value,
                        guess,
                        0.0001 * guess,
                        1e-12,
                        100,
                    )?,
                }
            }
            Phase::Gas => {
                let rhomin = 1e-14;
                match crate::solvers::halley(&mut resid, 0.5 * (rhomin + rho_anc), 1e-8, 100) {
                    Ok(rho) => rho,
                    Err(_) => crate::solvers::brent(
                        |rho| self.px_value(key, t, rho) - value,
                        rhomin,
                        rho_anc,
                        LDBL_EPSILON,
                        1e-12,
                        100,
                    )
                    .map_err(|_| Error::Value(String::new()))?,
                }
            }
            _ => {
                return Err(Error::Value(
                    "phase to solver_for_rho_given_T_oneof_HSU is invalid".into(),
                ));
            }
        };
        let p = self.eos.pressure(t, rho);
        // Upstream tail: `_Q = -1` then `recalculate_singlephase_phase`.
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar: rho,
            phase: self.recalculated_singlephase_phase(t, p, rho),
            q: -1.0,
        })
    }

    /// Supercritical branch of `solver_for_rho_given_T_oneof_HSU`.
    fn rho_from_caloric_t_supercritical(
        &self,
        t: f64,
        value: f64,
        key: CaloricKey,
    ) -> Result<HeosState> {
        let mut rhoc = self.rhomolar_critical();
        let rhomin = 1e-10;
        let yc = self.px_value(key, t, rhoc);
        let ymin = self.px_value(key, t, rhomin);
        let y = value;
        let in_closed = |x1: f64, x2: f64, x: f64| x >= x1.min(x2) && x <= x1.max(x2);
        let f = |rho: f64| self.px_value(key, t, rho) - value;
        let rho = if in_closed(yc, ymin, y) {
            crate::solvers::brent(f, rhoc, rhomin, LDBL_EPSILON, 1e-9, 100)?
        } else if y < yc {
            // Increase rhoc until it bounds the solution
            let mut yc2 = yc;
            let mut step_count = 0;
            while !in_closed(ymin, yc2, y) {
                rhoc *= 1.1;
                yc2 = self.px_value(key, t, rhoc);
                if step_count > 30 {
                    return Err(Error::Value(format!(
                        "Even by increasing rhoc, not able to bound input; input {y} is not in range {yc2},{ymin}"
                    )));
                }
                step_count += 1;
            }
            crate::solvers::brent(f, rhomin, rhoc, LDBL_EPSILON, 1e-9, 100)?
        } else {
            return Err(Error::Value(format!(
                "input {y} is not in range {yc},{ymin}"
            )));
        };
        let p = self.eos.pressure(t, rho);
        let phase = if p < self.p_critical() {
            Phase::SupercriticalGas
        } else {
            Phase::Supercritical
        };
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar: rho,
            phase,
            q: -1.0,
        })
    }

    /// Upstream `recalculate_singlephase_phase` (pure-fluid branch).
    pub(crate) fn recalculated_singlephase_phase(&self, t: f64, p: f64, rho: f64) -> Phase {
        if p > self.p_critical() {
            if t > self.t_critical() {
                Phase::Supercritical
            } else {
                Phase::SupercriticalLiquid
            }
        } else if t > self.t_critical() {
            Phase::SupercriticalGas
        } else if rho > self.rhomolar_critical() {
            Phase::Liquid
        } else {
            Phase::Gas
        }
    }

    /// (Hmolar, P) flash — upstream `HSU_P_flash(iHmolar)`.
    pub fn hmolar_p_state(&self, hmolar: f64, p: f64) -> Result<HeosState> {
        self.px_state(p, hmolar, CaloricKey::Hmolar)
    }

    /// (P, Smolar) flash — upstream `HSU_P_flash(iSmolar)`.
    pub fn p_smolar_state(&self, p: f64, smolar: f64) -> Result<HeosState> {
        self.px_state(p, smolar, CaloricKey::Smolar)
    }

    /// (P, Umolar) flash — upstream `HSU_P_flash(iUmolar)`.
    pub fn p_umolar_state(&self, p: f64, umolar: f64) -> Result<HeosState> {
        self.px_state(p, umolar, CaloricKey::Umolar)
    }

    /// (Dmolar, Hmolar) flash — upstream `HSU_D_flash(iHmolar)`.
    pub fn dmolar_hmolar_state(&self, rhomolar: f64, hmolar: f64) -> Result<HeosState> {
        self.hsu_d_state(rhomolar, hmolar, CaloricKey::Hmolar)
    }

    /// (Dmolar, Smolar) flash — upstream `HSU_D_flash(iSmolar)`.
    pub fn dmolar_smolar_state(&self, rhomolar: f64, smolar: f64) -> Result<HeosState> {
        self.hsu_d_state(rhomolar, smolar, CaloricKey::Smolar)
    }

    /// (Dmolar, Umolar) flash — upstream `HSU_D_flash(iUmolar)`.
    pub fn dmolar_umolar_state(&self, rhomolar: f64, umolar: f64) -> Result<HeosState> {
        self.hsu_d_state(rhomolar, umolar, CaloricKey::Umolar)
    }

    /// (Dmolar, Q) flash — upstream `DQ_flash`: for Q at a saturation
    /// boundary (|Q| or |Q-1| < 1e-10) the strict-mode superancillary path
    /// enumerates every T-root of the density on that branch and refuses to
    /// silently pick among several (GitHub #2773/#2834); fractional Q falls
    /// back to a Brent solve of the density-implied quality over T.
    pub fn dmolar_q_state(&self, rhomolar: f64, q: f64) -> Result<HeosState> {
        let q_boundary_tol = 1e-10;
        if q.abs() < q_boundary_tol || (q - 1.0).abs() < q_boundary_tol {
            let q_key = usize::from(q.abs() >= q_boundary_tol);
            let t = self.resolve_t_via_superancillary_d(rhomolar, q_key, "DQ_flash")?;
            let sat = self.sat().qt_flash(t, q)?;
            // Upstream restores the INPUT density and quality after QT_flash.
            return Ok(HeosState::TwoPhase {
                t,
                p: sat.p,
                rhomolar,
                q,
                rho_l: sat.rho_l,
                rho_v: sat.rho_v,

                t_l: t,

                t_v: t,
            });
        }
        // Fallback: Brent over [Tmin + 0.1, Tc - 0.1] on the density-implied
        // quality (the original DQ_flash).
        let eps = 1e-12;
        if rhomolar >= self.rhomolar_critical() + eps && q > eps {
            return Err(Error::OutOfRange(format!(
                "DQ inputs are not defined for density ({rhomolar}) above critical density ({}) and Q>0",
                self.rhomolar_critical()
            )));
        }
        let t_max = self.t_critical() - 0.1;
        let t_min = self.t_triple() + 0.1;
        let f = |t: f64| -> f64 {
            match self.sat().qt_flash(t, 0.0) {
                Ok(sat) => {
                    (1.0 / rhomolar - 1.0 / sat.rho_l) / (1.0 / sat.rho_v - 1.0 / sat.rho_l) - q
                }
                Err(_) => f64::NAN,
            }
        };
        let t = crate::solvers::brent(f, t_min, t_max, f64::EPSILON, 1e-10, 100)?;
        let sat = self.sat().qt_flash(t, q)?;
        Ok(HeosState::TwoPhase {
            t,
            p: sat.p,
            rhomolar,
            q,
            rho_l: sat.rho_l,
            rho_v: sat.rho_v,

            t_l: t,

            t_v: t,
        })
    }

    /// Upstream `resolve_T_via_superancillary` (no-guess mode) for the `'D'`
    /// key: enumerate every T-root on the requested saturation branch, dedup
    /// near-identical roots from adjacent intervals at an extremum, and
    /// return the single root — or throw for zero (out of range) or several
    /// (upstream `MultipleSolutionsError`).
    fn resolve_t_via_superancillary_d(
        &self,
        target: f64,
        q_key: usize,
        fn_name: &str,
    ) -> Result<f64> {
        let (d_l, d_v) = self.d_approxes();
        let approx = if q_key == 0 { d_l } else { d_v };
        let phase_name = if q_key == 0 { "liquid" } else { "vapor" };
        let mut solns = approx.get_x_for_y(target, 64, 100, 1e-10);
        if solns.is_empty() {
            return Err(Error::OutOfRange(format!(
                "{fn_name}: no T-root on saturated {phase_name} for D={target}; superancillary range [{}, {}] K",
                approx.xmin(),
                approx.xmax()
            )));
        }
        solns.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let dedup_tol = 1e-6;
        let mut ts: Vec<f64> = Vec::with_capacity(solns.len());
        for t in solns {
            if ts.last().is_none_or(|last| (t - last).abs() > dedup_tol) {
                ts.push(t);
            }
        }
        if ts.len() > 1 {
            let ts_str = ts
                .iter()
                .map(|t| format!("{t} K"))
                .collect::<Vec<_>>()
                .join(", ");
            return Err(Error::Value(format!(
                "{fn_name}: D={target} on saturated {phase_name} has {} T-roots ({ts_str}); use update_with_guesses with guess.T to pick a branch (see GitHub #2773)",
                ts.len()
            )));
        }
        Ok(ts[0])
    }

    /// Upstream `HSU_D_flash`'s superancillary "happy path" (every bundled
    /// fluid is pure with a superancillary, so the legacy ancillary "sad
    /// path" is only its error fallback — unported, loud error instead).
    ///
    /// Candidate T-intervals are cut at every intersection of the specified
    /// density with either saturation branch (`get_all_intersections`),
    /// classified by dome membership at their midpoint, and solved with the
    /// matching bracketed residual. The two-phase residual here is the
    /// full-EOS `Qo - Qd` form (upstream's `use_ca = false` mode): upstream's
    /// default caloric-superancillary fast path is followed by an EOS polish
    /// (`HSU_D_TWOPHASE_EOS_POLISH`, on by default) that converges to the
    /// same EOS root the direct residual finds — the fast path is a seed
    /// optimization, not a different answer.
    fn hsu_d_state(&self, rhomolar: f64, value: f64, key: CaloricKey) -> Result<HeosState> {
        let sa = self
            .fluid()
            .eos
            .superancillary
            .as_ref()
            .expect("HSU_D requires a superancillary fluid");
        let tcrit = sa.t_crit_num;
        // rhoV and rhoL coalesce at Tcrit and the quality ratio degenerates
        // to 0/0, so keep two-phase brackets just shy of it.
        let tcrit_2phase = tcrit - (1e-6f64).max(1e-9 * tcrit);
        let tmin_sa = sa.rho_l[0].xmin;
        let tmax_1phase = self.fluid().eos.t_max * 1.5;
        let tol = 1e-12;
        // 44 correct bits -> tolerance 2^(1-44) (upstream eps_tolerance<double>(44)).
        let eps44 = (2.0f64).powi(1 - 44);

        let inside_dome = |t: f64| -> bool {
            if t >= tcrit {
                return false;
            }
            let rho_v = crate::superancillary::eval_sat(sa, t, 'D', 1);
            let rho_l = crate::superancillary::eval_sat(sa, t, 'D', 0);
            rhomolar > rho_v && rhomolar < rho_l
        };
        // Bracketed root of `f` on [a, b] with endpoint short-circuits;
        // 44-bit bisection + midpoint stands in for TOMS748 as established.
        let bracket_solve =
            |f: &dyn Fn(f64) -> Result<f64>, a: f64, b: f64| -> Result<Option<f64>> {
                let fa = f(a)?;
                if fa.abs() < tol {
                    return Ok(Some(a));
                }
                let fb = f(b)?;
                if fb.abs() < tol {
                    return Ok(Some(b));
                }
                if fa * fb >= 0.0 {
                    return Ok(None);
                }
                let (mut l, mut r, mut fl) = (a, b, fa);
                for _ in 0..100 {
                    if (r - l) <= eps44 * l.abs().max(r.abs()) {
                        break;
                    }
                    let m = 0.5 * (l + r);
                    let fm = f(m)?;
                    if fm == 0.0 {
                        l = m;
                        r = m;
                        break;
                    }
                    if (fl < 0.0) == (fm < 0.0) {
                        l = m;
                        fl = fm;
                    } else {
                        r = m;
                    }
                }
                Ok(Some(0.5 * (l + r)))
            };

        // Single-phase: keyed value at the fixed density, solved over T.
        let solve_1phase = |a: f64, b: f64| -> Option<HeosState> {
            let f = |t: f64| -> Result<f64> { Ok(self.px_value(key, t, rhomolar) - value) };
            let tconv = bracket_solve(&f, a, b).ok().flatten()?;
            // Reject a converged root inside the dome (metastable); interval
            // classification should never produce one, but the SA and EOS
            // saturation curves differ at the ~1e-8 level near the boundary.
            if inside_dome(tconv) {
                return None;
            }
            let p = self.eos.pressure(tconv, rhomolar);
            Some(HeosState::SinglePhase {
                t: tconv,
                p,
                rhomolar,
                phase: self.recalculated_singlephase_phase(tconv, p, rhomolar),
                // Upstream `finalize_1phase` leaves the _Q = 10000 sentinel.
                q: 10000.0,
            })
        };
        // Two-phase: Qo - Qd residual with SA saturation densities and
        // full-EOS caloric values of each phase.
        let qd_at = |t: f64| -> (f64, f64, f64) {
            let rho_l = crate::superancillary::eval_sat(sa, t, 'D', 0);
            let rho_v = crate::superancillary::eval_sat(sa, t, 'D', 1);
            let qd = (1.0 / rhomolar - 1.0 / rho_l) / (1.0 / rho_v - 1.0 / rho_l);
            (qd, rho_l, rho_v)
        };
        let solve_2phase = |a: f64, b: f64| -> Option<HeosState> {
            let f = |t: f64| -> Result<f64> {
                let (qd, rho_l, rho_v) = qd_at(t);
                let y_l = self.px_value(key, t, rho_l);
                let y_v = self.px_value(key, t, rho_v);
                let qo = (value - y_l) / (y_v - y_l);
                let resid = qo - qd;
                if !resid.is_finite() {
                    return Err(Error::Value(format!(
                        "HSU_D superancillary resid not finite @ T={t} K; Qo={qo}; Qd={qd}"
                    )));
                }
                Ok(resid)
            };
            let tsol = bracket_solve(&f, a, b).ok().flatten()?;
            let (qd_final, _, _) = qd_at(tsol);
            // Reject a spurious crossing: quality outside [0, 1] means the
            // specified density is not inside the two-phase band at Tsol.
            let qeps = 1e-8;
            if !(-qeps..=1.0 + qeps).contains(&qd_final) {
                return None;
            }
            self.qt_state(tsol, qd_final.clamp(0.0, 1.0)).ok()
        };
        // The committed state must reproduce BOTH inputs (upstream
        // `committed_ok`): a solve can converge on a spurious root.
        let committed_ok = |state: &HeosState| -> bool {
            let rho_out = state.rhomolar();
            if !rho_out.is_finite() || (rho_out / rhomolar - 1.0).abs() > 1e-7 {
                return false;
            }
            let x_out = match state {
                HeosState::SinglePhase { t, rhomolar, .. } => self.px_value(key, *t, *rhomolar),
                HeosState::TwoPhase {
                    t, q, rho_l, rho_v, ..
                } => {
                    let y_l = self.px_value(key, *t, *rho_l);
                    let y_v = self.px_value(key, *t, *rho_v);
                    q * y_v + (1.0 - q) * y_l
                }
            };
            x_out.is_finite() && (x_out - value).abs() <= 1e-6 * value.abs() + 1e-3
        };

        // Candidate interval edges: every saturation intersection of the
        // density (both branches), between Tmin and Tcrit, plus the
        // supercritical continuation.
        let (d_l, d_v) = self.d_approxes();
        let mut tsats = d_l.get_x_for_y(rhomolar, 48, 100, 1e-13);
        tsats.extend(d_v.get_x_for_y(rhomolar, 48, 100, 1e-13));
        tsats.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let mut edges = vec![tmin_sa];
        for t in &tsats {
            if *t > tmin_sa && *t < tcrit {
                edges.push(*t);
            }
        }
        edges.push(tcrit);
        edges.push(tmax_1phase);

        for w in edges.windows(2) {
            let (a, b) = (w[0], w[1]);
            if b - a < 1e-10 {
                continue;
            }
            let mid = 0.5 * (a + b);
            let solved = if mid < tcrit && inside_dome(mid) {
                let ub = b.min(tcrit_2phase);
                if ub > a {
                    solve_2phase(a, ub)
                } else {
                    solve_1phase(a, b)
                }
            } else {
                solve_1phase(a, b)
            };
            if let Some(state) = solved {
                if committed_ok(&state) {
                    return Ok(state);
                }
            }
        }
        Err(Error::Value(
            "HSU_D superancillary: no candidate interval reproduced the inputs (the legacy ancillary path is not ported)"
                .into(),
        ))
    }

    fn px_state(&self, p: f64, value: f64, key: CaloricKey) -> Result<HeosState> {
        if self.fluid().eos.pseudo_pure {
            return self.px_state_pseudo_pure(p, value, key);
        }
        let pc = self.p_critical();
        let tc = self.t_critical();
        let rhoc = self.rhomolar_critical();
        let p_triple = self.fluid().eos.sat_min_liquid.p;

        // Upstream `p_phase_determination_pure_or_pseudopure(other=H/S)`
        let phase = if p > pc {
            let crit_value = self.px_value(key, tc, rhoc);
            if value > crit_value {
                Phase::SupercriticalGas
            } else {
                Phase::SupercriticalLiquid
            }
        } else if p >= p_triple * 0.9999 {
            if p > self.sat().pmax() {
                return Err(Error::Value(format!(
                    "Pressure to PQ_flash [{p:.8e} Pa] may not be above the numerical critical point of {:.15} Pa",
                    self.sat().pmax()
                )));
            }
            let sat = self.sat().pq_flash(p, 0.0)?;
            let liq = self.px_value(key, sat.t, sat.rho_l);
            let vap = self.px_value(key, sat.t, sat.rho_v);
            let q = (value - liq) / (vap - liq);
            if q < -1e-9 {
                Phase::Liquid
            } else if q > 1.0 + 1e-9 {
                Phase::Gas
            } else {
                let rhomolar = 1.0 / (q / sat.rho_v + (1.0 - q) / sat.rho_l);
                return Ok(HeosState::TwoPhase {
                    t: sat.t,
                    p,
                    rhomolar,
                    q,
                    rho_l: sat.rho_l,
                    rho_v: sat.rho_v,

                    t_l: sat.t,

                    t_v: sat.t,
                });
            }
        } else {
            // Below the triple-point pressure no saturation exists — the
            // state is gas (upstream's determination; liquid-like inputs
            // fail in the gas-bracket solve exactly as upstream).
            Phase::Gas
        };

        // Single-phase: bracket in T per upstream `HSU_P_flash`, then solve.
        let t_min_fluid = self.t_triple(); // upstream Tmin() == Ttriple() == sat_min_liquid.T
        let (t_min, t_max) = match phase {
            Phase::Gas => {
                let t_max = 1.5 * self.fluid().eos.t_max;
                let t_min = if p < p_triple {
                    t_min_fluid
                } else {
                    self.sat().pq_flash(p, 0.0)?.t
                };
                (t_min, t_max)
            }
            Phase::Liquid => (self.px_t_floor(p)?, self.sat().pq_flash(p, 0.0)?.t),
            Phase::SupercriticalLiquid | Phase::SupercriticalGas | Phase::Supercritical => {
                (self.px_t_floor(p)?, 1.5 * self.fluid().eos.t_max)
            }
            _ => return Err(Error::Value("Not a valid homogeneous state".into())),
        };

        self.px_solve_single_phase(p, value, key, phase, t_min, t_max)
    }

    /// The single-phase (p, X) temperature solve shared by the pure and
    /// pseudo-pure `HSU_P_flash` arms (upstream
    /// `HSU_P_flash_singlephase_Brent`), plus the flash tail.
    ///
    /// Upstream evaluates the residual at both bracket ends, hands the pair to
    /// `boost::math::tools::toms748_solve` with `eps_tolerance<double>(30)` and
    /// `max_iter = 100`, then RE-EVALUATES at the midpoint of the final bracket
    /// with the probe's iteration counter reset — so the state served is that
    /// midpoint evaluation's, taken on the cold density path.
    fn px_solve_single_phase(
        &self,
        p: f64,
        value: f64,
        key: CaloricKey,
        phase: Phase,
        t_min: f64,
        t_max: f64,
    ) -> Result<HeosState> {
        let mut resid = PxResid {
            flash: self,
            p,
            value,
            key,
            phase,
            live_phase: phase,
            iter: 0,
            rhomolar: f64::INFINITY,
            rhomolar0: f64::INFINITY,
            rhomolar1: f64::INFINITY,
            // Upstream `solver_resid::p_crit`, read at construction.
            force_robust_density: p > self.p_critical(),
        };
        let r_min = resid.call(t_min)?;
        let r_max = resid.call(t_max)?;
        if r_min * r_max >= 0.0 {
            return Err(Error::Solution(format!(
                "unable to bracket the (p,X) solution in [{t_min}, {t_max}] (residuals {r_min:e}, {r_max:e}); the derivative path is not ported"
            )));
        }
        resid.iter = 0;
        let t = crate::solvers::toms748_solve(
            &mut |t: f64| resid.call(t),
            t_min,
            t_max,
            r_min,
            r_max,
            30,
            100,
        )?;
        // Upstream re-evaluates at `0.5 * (l + r)` after the solve — the
        // returned midpoint — with `iter` reset, i.e. on the cold density
        // path, and serves THAT state.
        resid.iter = 0;
        resid.call(t)?;
        let rho = resid.rhomolar;
        // Upstream's post-solve range guard. `toms748_solve` returns the
        // midpoint of a bracket that never leaves `[t_min, t_max]`, so this
        // arm is unreachable here; it is carried because it is upstream's.
        if !(t_min..=t_max).contains(&t) {
            return Err(Error::Value(format!(
                "TOMS748 method yielded out of bound T of {t}"
            )));
        }
        // Upstream `HSU_P_flash` tail: `recalculate_singlephase_phase` runs
        // after the solve, so the bracket's working label gives way to the
        // final (T, p, rho) quadrant (wheel: PS steam at 1568 K -> 2).
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar: rho,
            phase: self.recalculated_singlephase_phase(t, p, rho),
            q: -1.0,
        })
    }

    /// Pseudo-pure `HSU_P_flash` (upstream's classic-ancillary arm): the
    /// pseudo-pure phase determination, then the same bracketed single-phase
    /// temperature solve the pure-fluid path uses, with upstream's bracket
    /// selection — `SatV->T()`/`SatL->T()` when the slow VLE ran,
    /// `_TVanc - 0.01`/`_TLanc + 0.01` when the ancillary bands classified.
    fn px_state_pseudo_pure(&self, p: f64, value: f64, key: CaloricKey) -> Result<HeosState> {
        let det = self.p_phase_determination_pseudo_pure(PpOther::Caloric(key), value, p)?;
        let (phase, t_sat, t_anc) = match det {
            PpPhaseDet::TwoPhase(state) => return Ok(state),
            PpPhaseDet::Single {
                phase,
                t_sat,
                t_anc,
            } => (phase, t_sat, t_anc),
        };
        let eos = &self.fluid().eos;
        let p_triple = eos.sat_min_liquid.p;
        let (t_min, t_max) = match phase {
            Phase::Gas => {
                let t_max = 1.5 * eos.t_max;
                let t_min = if p < p_triple {
                    // Upstream `max(Tmin(), Ttriple())` — both resolve to
                    // `sat_min_liquid.T`.
                    self.t_triple()
                } else if let Some((_, t_v_sat)) = t_sat {
                    t_v_sat
                } else {
                    t_anc
                        .expect("subcritical gas classification carries TVanc")
                        .1
                        - 0.01
                };
                (t_min, t_max)
            }
            Phase::Liquid => {
                let t_max = if let Some((t_l_sat, _)) = t_sat {
                    t_l_sat
                } else {
                    t_anc
                        .expect("subcritical liquid classification carries TLanc")
                        .0
                        + 0.01
                };
                (self.px_t_floor(p)?, t_max)
            }
            Phase::SupercriticalLiquid | Phase::SupercriticalGas | Phase::Supercritical => {
                (self.px_t_floor(p)?, 1.5 * eos.t_max)
            }
            _ => return Err(Error::Value("Not a valid homogeneous state".into())),
        };
        self.px_solve_single_phase(p, value, key, phase, t_min, t_max)
    }

    /// Upstream `p_phase_determination_pure_or_pseudopure` for a PSEUDO-PURE
    /// fluid (the classic-ancillary arm; pure fluids take the superancillary
    /// paths inlined in `px_state`/`dmolar_p_state`):
    /// - `p > max_sat_p.p`: supercritical split on the STATES.critical
    ///   calorics (H/S/U) or density (D);
    /// - triple..max_sat_p: invert pL/pV for `_TLanc`/`_TVanc`, then the
    ///   rational-polynomial caloric bands (H/S at their fit error, U at
    ///   1.5x) or the 0.95/1.05 density bands (D), and the slow VLE
    ///   (pseudo-pure PQ at Q=0) with the lever rule when inconclusive;
    /// - below triple: gas outright (`other != iT` asks no questions).
    fn p_phase_determination_pseudo_pure(
        &self,
        other: PpOther,
        value: f64,
        p: f64,
    ) -> Result<PpPhaseDet> {
        let eos = &self.fluid().eos;
        // Upstream `calc_pmax_sat`: `max_sat_p.p` for a pseudo-pure fluid.
        let psat_max = eos
            .max_sat_p
            .as_ref()
            .expect("pseudo-pure fluids carry max_sat_p")
            .p;
        let p_triple = eos.sat_min_liquid.p;
        if p > psat_max {
            let phase = match other {
                PpOther::Dmolar => {
                    if value < self.rhomolar_critical() {
                        Phase::SupercriticalGas
                    } else {
                        Phase::SupercriticalLiquid
                    }
                }
                PpOther::Caloric(key) => {
                    // `calc_{s,h,u}molar_nocache(T_critical, rhomolar_critical)`
                    // — the STATES.critical point for a pseudo-pure fluid.
                    let crit_value =
                        self.px_value(key, self.t_critical(), self.rhomolar_critical());
                    if value > crit_value {
                        Phase::SupercriticalGas
                    } else {
                        Phase::SupercriticalLiquid
                    }
                }
            };
            return Ok(PpPhaseDet::Single {
                phase,
                t_sat: None,
                t_anc: None,
            });
        }
        if p >= p_triple * 0.9999 {
            let anc = &self.fluid().ancillaries;
            let tl_anc = crate::ancillary::invert(&anc.p_s, p)?;
            let tv_anc = crate::ancillary::invert(anc.p_v_split.as_ref().unwrap_or(&anc.p_s), p)?;
            let t_anc = Some((tl_anc, tv_anc));
            let single = |phase: Phase, t_sat: Option<(f64, f64)>| PpPhaseDet::Single {
                phase,
                t_sat,
                t_anc,
            };
            let mut definitely_two_phase = false;
            if let PpOther::Caloric(key) = other {
                // (x_liq, x_liq_error_band, x_vap, x_vap_error_band), None
                // when the gating ancillary is absent (upstream `enabled()`;
                // all six pseudo-pure fluids carry all four curves).
                let bands = match key {
                    CaloricKey::Hmolar => anc.h_l.as_ref().map(|h_l| {
                        let h_lv = anc.h_lv.as_ref().expect("hLV accompanies hL");
                        // Ancillaries are h - h_anchor, so add back h_anchor.
                        let h_liq = crate::ancillary::evaluate_rational_poly(h_l, tl_anc)
                            + eos.hs_anchor.hmolar;
                        let h_liq_band = h_l.max_abs_error;
                        let h_vap = h_liq + crate::ancillary::evaluate_rational_poly(h_lv, tl_anc);
                        (h_liq, h_liq_band, h_vap, h_liq_band + h_lv.max_abs_error)
                    }),
                    CaloricKey::Smolar => anc.s_l.as_ref().map(|s_l| {
                        let s_lv = anc.s_lv.as_ref().expect("sLV accompanies sL");
                        let s_liq = crate::ancillary::evaluate_rational_poly(s_l, tl_anc)
                            + eos.hs_anchor.smolar;
                        let s_liq_band = s_l.max_abs_error;
                        // Upstream evaluates sLV at `_TVanc` (hLV is at
                        // `_TLanc`) — the asymmetry is reproduced verbatim.
                        let s_vap = s_liq + crate::ancillary::evaluate_rational_poly(s_lv, tv_anc);
                        (s_liq, s_liq_band, s_vap, s_liq_band + s_lv.max_abs_error)
                    }),
                    CaloricKey::Umolar => anc.h_l.as_ref().map(|h_l| {
                        // u = h - p/rho off the enthalpy + density
                        // ancillaries; "most of error is in enthalpy", so
                        // the bands are 1.5x the enthalpy bands.
                        let h_lv = anc.h_lv.as_ref().expect("hLV accompanies hL");
                        let h_liq = crate::ancillary::evaluate_rational_poly(h_l, tl_anc)
                            + eos.hs_anchor.hmolar;
                        let h_liq_band = h_l.max_abs_error;
                        let h_vap = h_liq + crate::ancillary::evaluate_rational_poly(h_lv, tl_anc);
                        let h_vap_band = h_liq_band + h_lv.max_abs_error;
                        let rho_vap = crate::ancillary::evaluate(&anc.rho_v, tv_anc);
                        let rho_liq = crate::ancillary::evaluate(&anc.rho_l, tl_anc);
                        (
                            h_liq - p / rho_liq,
                            1.5 * h_liq_band,
                            h_vap - p / rho_vap,
                            1.5 * h_vap_band,
                        )
                    }),
                };
                if let Some((x_liq, x_liq_band, x_vap, x_vap_band)) = bands {
                    if value > x_vap + x_vap_band {
                        return Ok(single(Phase::Gas, None));
                    } else if value < x_liq - x_liq_band {
                        return Ok(single(Phase::Liquid, None));
                    } else if value > x_liq + x_liq_band && value < x_vap - x_vap_band {
                        definitely_two_phase = true;
                    }
                }
            }
            // Upstream's !definitely_two_phase block always evaluates the
            // rhoV/rhoL ancillaries but only the iDmolar case decides —
            // Dmolar never sets the flag, so the guard is kept for shape.
            if !definitely_two_phase {
                if let PpOther::Dmolar = other {
                    let rho_vap = 0.95 * crate::ancillary::evaluate(&anc.rho_v, tv_anc);
                    let rho_liq = 1.05 * crate::ancillary::evaluate(&anc.rho_l, tl_anc);
                    if value < rho_vap {
                        return Ok(single(Phase::Gas, None));
                    } else if value > rho_liq {
                        return Ok(single(Phase::Liquid, None));
                    }
                }
            }
            // The slow full VLE calculation is required (upstream: a fresh
            // backend through `PQ_flash` at Q=0 — the pseudo-pure branch,
            // whose SatL/SatV sit at the bubble/dew temperatures).
            let HeosState::TwoPhase {
                rho_l,
                rho_v,
                t_l,
                t_v,
                ..
            } = self.pq_state_pseudo_pure(p, 0.0)?
            else {
                unreachable!("pseudo-pure PQ always returns a two-phase state")
            };
            let q = match other {
                PpOther::Dmolar => (1.0 / value - 1.0 / rho_l) / (1.0 / rho_v - 1.0 / rho_l),
                PpOther::Caloric(key) => {
                    let y_l = self.px_value(key, t_l, rho_l);
                    let y_v = self.px_value(key, t_v, rho_v);
                    (value - y_l) / (y_v - y_l)
                }
            };
            let t_sat = Some((t_l, t_v));
            if q < -1e-9 {
                return Ok(single(Phase::Liquid, t_sat));
            } else if q > 1.0 + 1e-9 {
                return Ok(single(Phase::Gas, t_sat));
            }
            // Two-phase: upstream loads T/rho straight off the lever rule
            // (raw mixes, no endpoint shortcuts here).
            return Ok(PpPhaseDet::TwoPhase(HeosState::TwoPhase {
                t: q * t_v + (1.0 - q) * t_l,
                p,
                rhomolar: 1.0 / (q / rho_v + (1.0 - q) / rho_l),
                q,
                rho_l,
                rho_v,
                t_l,
                t_v,
            }));
        }
        if p < p_triple * 0.9999 {
            // `other != iT` asks no further questions below the triple point.
            return Ok(PpPhaseDet::Single {
                phase: Phase::Gas,
                t_sat: None,
                t_anc: None,
            });
        }
        // Only reachable for a NaN pressure (all three range tests false).
        Err(Error::Value(format!(
            "The pressure [{p} Pa] cannot be used in p_phase_determination"
        )))
    }

    /// Liquid/supercritical lower T-bound of the (p,X) bracket (upstream
    /// `HSU_P_flash`): the melting temperature at this pressure when a
    /// melting line covers it, else Tmin — both minus 1e-3.
    fn px_t_floor(&self, p: f64) -> Result<f64> {
        if let Some(ml) = &self.fluid().ancillaries.melting_line {
            if p > crate::melting::p_min(ml) {
                return Ok(crate::melting::t_of_p(ml, p)? - 1e-3);
            }
        }
        Ok(self.t_triple() - 1e-3)
    }

    fn px_value(&self, key: CaloricKey, t: f64, rhomolar: f64) -> f64 {
        match key {
            CaloricKey::Hmolar => self.eos.hmolar(t, rhomolar),
            CaloricKey::Smolar => self.eos.smolar(t, rhomolar),
            CaloricKey::Umolar => self.eos.umolar(t, rhomolar),
        }
    }

    /// Inner density for a COLD (p,T) probe — upstream `update(PT_INPUTS, p,
    /// T)`: imposed phase for liquid/gas, full determination for the
    /// supercritical classifications (which can legitimately flip along the
    /// bracket). Returns the density AND the phase the state is left holding
    /// (upstream `_phase`, which `solver_rho_Tp`'s guessed form consults when
    /// no phase is imposed — `clear()` does not reset it).
    fn px_probe_rho(&self, t: f64, p: f64, phase: Phase) -> Result<(f64, Phase)> {
        match phase {
            Phase::Liquid | Phase::Gas => Ok((self.solver_rho_tp(t, p, phase)?, phase)),
            _ => self.pt_flash(t, p),
        }
    }

    /// (Dmolar, P) flash — upstream `DP_flash`: superancillary (p, Dmolar)
    /// phase determination, then a Halley solve for T at fixed density
    /// (Peng-Robinson seed for gas-like phases, saturation/1.1*Tc seeds
    /// otherwise), with the 30-bit bracketed fallback.
    pub fn dmolar_p_state(&self, rhomolar: f64, p: f64) -> Result<HeosState> {
        if self.fluid().eos.pseudo_pure {
            return self.dmolar_p_state_pseudo_pure(rhomolar, p);
        }
        let pc = self.p_critical();
        let tc = self.t_critical();
        let rhoc = self.rhomolar_critical();
        let p_triple = self.fluid().eos.sat_min_liquid.p;

        // The determination's label only seeds the T solve now — the final
        // label is recalculated from the converged state below.
        let (_phase, t0) = if p > pc {
            if rhomolar < rhoc {
                (
                    Phase::SupercriticalGas,
                    self.t_dp_peng_robinson(rhomolar, p),
                )
            } else {
                (Phase::SupercriticalLiquid, 1.1 * tc)
            }
        } else if p >= p_triple * 0.9999 {
            if p > self.sat().pmax() {
                return Err(Error::Value(format!(
                    "Pressure to PQ_flash [{p:.8e} Pa] may not be above the numerical critical point of {:.15} Pa",
                    self.sat().pmax()
                )));
            }
            let sat = self.sat().pq_flash(p, 0.0)?;
            let q = (1.0 / rhomolar - 1.0 / sat.rho_l) / (1.0 / sat.rho_v - 1.0 / sat.rho_l);
            if q < -1e-9 {
                (Phase::Liquid, sat.t)
            } else if q > 1.0 + 1e-9 {
                (Phase::Gas, self.t_dp_peng_robinson(rhomolar, p))
            } else {
                // Upstream recomputes the mixture density from the mixing rule
                // even though density was the input.
                let rhomix = 1.0 / (q / sat.rho_v + (1.0 - q) / sat.rho_l);
                return Ok(HeosState::TwoPhase {
                    t: sat.t,
                    p,
                    rhomolar: rhomix,
                    q,
                    rho_l: sat.rho_l,
                    rho_v: sat.rho_v,

                    t_l: sat.t,

                    t_v: sat.t,
                });
            }
        } else {
            return Err(Error::NotImplemented(
                "(Dmolar,P) flash below the triple-point pressure is not ported yet".into(),
            ));
        };

        if !t0.is_finite() {
            return Err(Error::Value(
                "Starting value of T0 is not valid in DP_flash".into(),
            ));
        }
        let t = self.solve_t_dp(rhomolar, p, t0)?;
        // Upstream `DP_flash` tail: `_Q = -1` then
        // `recalculate_singlephase_phase` — the seed label is only the
        // solver's working hypothesis.
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar,
            phase: self.recalculated_singlephase_phase(t, p, rhomolar),
            q: -1.0,
        })
    }

    /// Pseudo-pure `DP_flash` (upstream's classic-ancillary arm): the
    /// pseudo-pure phase determination, then upstream's T0 seed per phase —
    /// `SatL->T()`/`_TLanc` for liquid, `1.1*T_critical` for supercritical
    /// liquid, Peng-Robinson for the gas-like labels — into the shared
    /// Halley/bracket temperature solve.
    fn dmolar_p_state_pseudo_pure(&self, rhomolar: f64, p: f64) -> Result<HeosState> {
        let det = self.p_phase_determination_pseudo_pure(PpOther::Dmolar, rhomolar, p)?;
        let (phase, t_sat, t_anc) = match det {
            PpPhaseDet::TwoPhase(state) => return Ok(state),
            PpPhaseDet::Single {
                phase,
                t_sat,
                t_anc,
            } => (phase, t_sat, t_anc),
        };
        let t0 = match phase {
            Phase::Liquid => {
                if let Some((t_l_sat, _)) = t_sat {
                    t_l_sat
                } else {
                    t_anc
                        .expect("subcritical liquid classification carries TLanc")
                        .0
                }
            }
            Phase::SupercriticalLiquid => 1.1 * self.t_critical(),
            Phase::Gas | Phase::SupercriticalGas | Phase::Supercritical => {
                self.t_dp_peng_robinson(rhomolar, p)
            }
            _ => return Err(Error::Value("I should never get here".into())),
        };
        if !t0.is_finite() {
            return Err(Error::Value(
                "Starting value of T0 is not valid in DP_flash".into(),
            ));
        }
        let t = self.solve_t_dp(rhomolar, p, t0)?;
        // Upstream `DP_flash` tail: `_Q = -1` then
        // `recalculate_singlephase_phase`.
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar,
            phase: self.recalculated_singlephase_phase(t, p, rhomolar),
            q: -1.0,
        })
    }

    /// Upstream `T_DP_PengRobinson`: PR-based T seed at fixed (rho, p).
    fn t_dp_peng_robinson(&self, rhomolar: f64, p: f64) -> f64 {
        let omega = self.fluid().eos.acentric;
        let tc = self.t_critical();
        let pc = self.p_critical();
        let r = self.eos.gas_constant;
        let v = 1.0 / rhomolar;
        let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        let a = 0.457235 * r * r * tc * tc / pc;
        let b = 0.077796 * r * tc / pc;
        let den = v * v + 2.0 * b * v - b * b;
        let big_a = r * tc / (v - b) - a * kappa * kappa / den;
        let big_b = 2.0 * a * kappa * (1.0 + kappa) / den;
        let big_c = -a * (1.0 + 2.0 * kappa + kappa * kappa) / den - p;
        let sqrt_tr1 = (-big_b + (big_b * big_b - 4.0 * big_a * big_c).sqrt()) / (2.0 * big_a);
        sqrt_tr1 * sqrt_tr1 * tc
    }

    /// dp/dT at constant rho and its T-derivative (for the Halley solve).
    fn dpdt_rho(&self, t: f64, rhomolar: f64) -> (f64, f64) {
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let tau = self.eos.t_reducing / t;
        let d = self.eos.alphar_all(tau, delta);
        let r = self.eos.gas_constant;
        let dpdt = rhomolar * r * (1.0 + delta * d.d10 - delta * tau * d.d11);
        let d2pdt2 = rhomolar * r * delta * tau * tau * d.d12 / t;
        (dpdt, d2pdt2)
    }

    /// Upstream `DP_flash` solver: Halley on (p(T,rho)-p)/p with acceptance
    /// checks, then the 30-bit bracketed fallback on [Tmin, 1.5*Tmax].
    fn solve_t_dp(&self, rhomolar: f64, p: f64, t0: f64) -> Result<f64> {
        struct DpResid<'a> {
            flash: &'a PtFlash,
            rhomolar: f64,
            p: f64,
        }
        impl crate::solvers::Resid1D for DpResid<'_> {
            fn call(&mut self, t: f64) -> f64 {
                (self.flash.eos.pressure(t, self.rhomolar) - self.p) / self.p
            }
            fn deriv(&mut self, t: f64) -> f64 {
                self.flash.dpdt_rho(t, self.rhomolar).0 / self.p
            }
            fn second_deriv(&mut self, t: f64) -> f64 {
                self.flash.dpdt_rho(t, self.rhomolar).1 / self.p
            }
            fn third_deriv(&mut self, _t: f64) -> f64 {
                unreachable!("Halley does not use the third derivative")
            }
        }
        let mut resid = DpResid {
            flash: self,
            rhomolar,
            p,
        };
        let t_hi = 1.5 * self.fluid().eos.t_max;
        use crate::solvers::Resid1D as _;
        if let Ok(t) = crate::solvers::halley(&mut resid, t0, 1e-10, 100) {
            if t.is_finite() && t > 0.0 && t <= t_hi && resid.call(t).abs() < 1e-7 {
                return Ok(t);
            }
        }
        // Bracketed fallback: p(T, rho) is monotone in T at fixed rho.
        let t_lo = self.t_triple();
        let (f_lo, f_hi) = (resid.call(t_lo), resid.call(t_hi));
        if !(f_lo.is_finite() && f_hi.is_finite() && f_lo * f_hi < 0.0) {
            return Err(Error::Solution(format!(
                "DP_flash could not bracket T for rho={rhomolar}, p={p}"
            )));
        }
        let tol = (2.0f64).powi(1 - 30);
        let (mut a, mut b) = (t_lo, t_hi);
        let mut fa = f_lo;
        for _ in 0..200 {
            if (b - a) <= tol * a.abs().max(b.abs()) {
                break;
            }
            let m = 0.5 * (a + b);
            let fm = resid.call(m);
            if fm == 0.0 {
                a = m;
                b = m;
                break;
            }
            if (fa < 0.0) == (fm < 0.0) {
                a = m;
                fa = fm;
            } else {
                b = m;
            }
        }
        Ok(0.5 * (a + b))
    }
}

#[derive(Clone, Copy)]
enum CaloricKey {
    Hmolar,
    Smolar,
    Umolar,
}

/// The `other` parameter of the pseudo-pure
/// `p_phase_determination_pure_or_pseudopure` arm (H/S/U come from
/// `HSU_P_flash`, D from `DP_flash`; T is the already-ported PT path).
#[derive(Clone, Copy)]
enum PpOther {
    Caloric(CaloricKey),
    Dmolar,
}

/// Outcome of the pseudo-pure phase determination.
enum PpPhaseDet {
    /// The determination itself loaded the full two-phase state.
    TwoPhase(HeosState),
    /// Homogeneous phase, plus what the flash's bracket/seed selection
    /// needs (upstream's `saturation_called` + `SatL`/`SatV` or the cached
    /// `_TLanc`/`_TVanc` members).
    Single {
        phase: Phase,
        /// `(SatL->T(), SatV->T())` when the slow VLE ran
        /// (`saturation_called` true).
        t_sat: Option<(f64, f64)>,
        /// `(_TLanc, _TVanc)` when the subcritical band inverted the
        /// pressure ancillaries.
        t_anc: Option<(f64, f64)>,
    },
}
