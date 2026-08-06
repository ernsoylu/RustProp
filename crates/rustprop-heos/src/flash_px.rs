//! Remaining pure-fluid flash pairs (PLAN.md 4.6): general-quality (T,Q) and
//! (P,Q) states with upstream's two-phase mixing, the (D,T) flash
//! (`DHSU_T_flash` -> superancillary `T_phase_determination`), and the
//! (H,P)/(P,S) flashes (`HSU_P_flash`: superancillary `p_phase_determination`
//! plus the bracketed single-phase solve in T).
//!
//! Numerical notes, logged in PLAN.md:
//! - upstream resolves the single-phase (P,X) temperature with TOMS748 at a
//!   deliberate 30-bit (~1e-9 relative) tolerance, then re-evaluates at the
//!   bracket midpoint; we bisect to the same relative tolerance — golden
//!   agreement is bounded by that shared tolerance, hence the 1e-8 policy;
//! - the no-bracket derivative path and the 2-D Newton fallback are unported
//!   (loud error) until a state needs them;
//! - two-phase mixing is upstream's exact `Q*V + (1-Q)*L` with the
//!   DBL_EPSILON endpoint shortcuts.

use crate::alpha::HelmholtzEos;
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
struct SmolarTResid<'a> {
    eos: &'a HelmholtzEos,
    t: f64,
    target: f64,
}

impl Resid1D for SmolarTResid<'_> {
    fn call(&mut self, rhomolar: f64) -> f64 {
        self.eos.smolar(self.t, rhomolar) - self.target
    }
    /// ds/drho|T = R*(tau*d11 - 1/delta - d10)/rho_r
    fn deriv(&mut self, rhomolar: f64) -> f64 {
        let tau = self.eos.t_reducing / self.t;
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let d = self.eos.alphar_all(tau, delta);
        self.eos.gas_constant * (tau * d.d11 - 1.0 / delta - d.d10) / self.eos.rhomolar_reducing
    }
    /// d2s/drho2|T = R*(tau*d21 + 1/delta^2 - d20)/rho_r^2
    fn second_deriv(&mut self, rhomolar: f64) -> f64 {
        let tau = self.eos.t_reducing / self.t;
        let delta = rhomolar / self.eos.rhomolar_reducing;
        let d = self.eos.alphar_all(tau, delta);
        self.eos.gas_constant * (tau * d.d21 + 1.0 / (delta * delta) - d.d20)
            / (self.eos.rhomolar_reducing * self.eos.rhomolar_reducing)
    }
    fn third_deriv(&mut self, _rhomolar: f64) -> f64 {
        unreachable!("Halley does not use the third derivative")
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
                t, q, rho_l, rho_v, ..
            } => mix_two_phase(*q, self.eos.hmolar(*t, *rho_l), self.eos.hmolar(*t, *rho_v)),
        }
    }
    /// Molar entropy of a state [J/mol/K].
    pub fn state_smolar(&self, s: &HeosState) -> f64 {
        match s {
            HeosState::SinglePhase { t, rhomolar, .. } => self.eos.smolar(*t, *rhomolar),
            HeosState::TwoPhase {
                t, q, rho_l, rho_v, ..
            } => mix_two_phase(*q, self.eos.smolar(*t, *rho_l), self.eos.smolar(*t, *rho_v)),
        }
    }
    /// Molar internal energy of a state [J/mol].
    pub fn state_umolar(&self, s: &HeosState) -> f64 {
        match s {
            HeosState::SinglePhase { t, rhomolar, .. } => self.eos.umolar(*t, *rhomolar),
            HeosState::TwoPhase {
                t, q, rho_l, rho_v, ..
            } => mix_two_phase(*q, self.eos.umolar(*t, *rho_l), self.eos.umolar(*t, *rho_v)),
        }
    }

    /// General-quality (T,Q) state (superancillary `QT_flash`).
    pub fn qt_state(&self, t: f64, q: f64) -> Result<HeosState> {
        let sat = self.sat().qt_flash(t, q)?;
        Ok(HeosState::TwoPhase {
            t: sat.t,
            p: sat.p,
            rhomolar: sat.rhomolar,
            q,
            rho_l: sat.rho_l,
            rho_v: sat.rho_v,
        })
    }

    /// General-quality (P,Q) state (superancillary `PQ_flash`).
    pub fn pq_state(&self, p: f64, q: f64) -> Result<HeosState> {
        let sat = self.sat().pq_flash(p, q)?;
        Ok(HeosState::TwoPhase {
            t: sat.t,
            p: sat.p,
            rhomolar: sat.rhomolar,
            q,
            rho_l: sat.rho_l,
            rho_v: sat.rho_v,
        })
    }

    /// (Dmolar, T) flash — upstream `DHSU_T_flash(iDmolar)`, which for a pure
    /// fluid is exactly `T_phase_determination_pure_or_pseudopure(iDmolar)`.
    pub fn dmolar_t_state(&self, rhomolar: f64, t: f64) -> Result<HeosState> {
        let tc = self.fluid().states.critical.t;
        let rhoc = self.fluid().states.critical.rhomolar;
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
            if q <= 0.0 {
                Ok(HeosState::SinglePhase {
                    t,
                    p: self.eos.pressure(t, rhomolar),
                    rhomolar,
                    phase: Phase::Liquid,
                    q: -1.0,
                })
            } else if q >= 1.0 {
                Ok(HeosState::SinglePhase {
                    t,
                    p: self.eos.pressure(t, rhomolar),
                    rhomolar,
                    phase: Phase::Gas,
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
                })
            }
        } else if t > tc && t > self.fluid().eos.t_triple {
            let phase = if rhomolar > rhoc {
                Phase::SupercriticalLiquid
            } else {
                Phase::SupercriticalGas
            };
            Ok(HeosState::SinglePhase {
                t,
                p: self.eos.pressure(t, rhomolar),
                rhomolar,
                phase,
                q: -1.0,
            })
        } else {
            Err(Error::Value(
                "temperature is out of range in dmolar_t_state".into(),
            ))
        }
    }

    /// (Smolar, T) flash — upstream `DHSU_T_flash(iSmolar)`: superancillary
    /// phase determination (`T_phase_determination_pure_or_pseudopure`),
    /// then `solver_for_rho_given_T_oneof_HSU` for the single-phase
    /// branches. Ported for the legacy HS path; it is also the (S,T) input
    /// pair itself.
    pub fn smolar_t_state(&self, smolar: f64, t: f64) -> Result<HeosState> {
        let tc = self.fluid().states.critical.t;
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
            let s_l = self.eos.smolar(t, sat.rho_l);
            let s_v = self.eos.smolar(t, sat.rho_v);
            let q = (smolar - s_l) / (s_v - s_l);
            if q < 0.0 {
                self.rho_from_smolar_t(t, smolar, Phase::Liquid, sat.rho_l)
            } else if q > 1.0 {
                self.rho_from_smolar_t(t, smolar, Phase::Gas, sat.rho_v)
            } else {
                Ok(HeosState::TwoPhase {
                    t,
                    p: sat.p,
                    rhomolar: 1.0 / (q / sat.rho_v + (1.0 - q) / sat.rho_l),
                    q,
                    rho_l: sat.rho_l,
                    rho_v: sat.rho_v,
                })
            }
        } else if t > tc && t > self.fluid().eos.t_triple {
            self.rho_from_smolar_t_supercritical(t, smolar)
        } else {
            Err(Error::Value(
                "temperature is out of range in smolar_t_state".into(),
            ))
        }
    }

    /// Subcritical single-phase branches of
    /// `solver_for_rho_given_T_oneof_HSU(iSmolar)`. `rho_anc` is the
    /// superancillary saturation density of the branch (upstream's
    /// `_rhoLanc`/`_rhoVanc` set by the phase determination).
    fn rho_from_smolar_t(
        &self,
        t: f64,
        smolar: f64,
        phase: Phase,
        rho_anc: f64,
    ) -> Result<HeosState> {
        let mut resid = SmolarTResid {
            eos: &self.eos,
            t,
            target: smolar,
        };
        let rho = match phase {
            Phase::Liquid => {
                let rhomelt = self.fluid().states.triple_liquid.rhomolar;
                let ymelt = self.eos.smolar(t, rhomelt);
                let y_l = self.eos.smolar(t, rho_anc);
                let guess = (rhomelt - rho_anc) / (ymelt - y_l) * (smolar - y_l) + rho_anc;
                match crate::solvers::halley(&mut resid, guess, 1e-8, 100) {
                    Ok(rho) => rho,
                    Err(_) => crate::solvers::secant(
                        |rho| self.eos.smolar(t, rho) - smolar,
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
                        |rho| self.eos.smolar(t, rho) - smolar,
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

    /// Supercritical branch of `solver_for_rho_given_T_oneof_HSU(iSmolar)`.
    fn rho_from_smolar_t_supercritical(&self, t: f64, smolar: f64) -> Result<HeosState> {
        let mut rhoc = self.fluid().states.critical.rhomolar;
        let rhomin = 1e-10;
        let yc = self.eos.smolar(t, rhoc);
        let ymin = self.eos.smolar(t, rhomin);
        let y = smolar;
        let in_closed = |x1: f64, x2: f64, x: f64| x >= x1.min(x2) && x <= x1.max(x2);
        let f = |rho: f64| self.eos.smolar(t, rho) - smolar;
        let rho = if in_closed(yc, ymin, y) {
            crate::solvers::brent(f, rhoc, rhomin, LDBL_EPSILON, 1e-9, 100)?
        } else if y < yc {
            // Increase rhoc until it bounds the solution
            let mut yc2 = yc;
            let mut step_count = 0;
            while !in_closed(ymin, yc2, y) {
                rhoc *= 1.1;
                yc2 = self.eos.smolar(t, rhoc);
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
        let phase = if p < self.fluid().states.critical.p {
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
        let crit = &self.fluid().states.critical;
        if p > crit.p {
            if t > crit.t {
                Phase::Supercritical
            } else {
                Phase::SupercriticalLiquid
            }
        } else if t > crit.t {
            Phase::SupercriticalGas
        } else if rho > crit.rhomolar {
            Phase::Liquid
        } else {
            Phase::Gas
        }
    }

    /// (Hmolar, P) flash — upstream `HSU_P_flash(iHmolar)`.
    pub fn hmolar_p_state(&self, hmolar: f64, p: f64) -> Result<HeosState> {
        self.px_state(p, hmolar, PxKey::Hmolar)
    }

    /// (P, Smolar) flash — upstream `HSU_P_flash(iSmolar)`.
    pub fn p_smolar_state(&self, p: f64, smolar: f64) -> Result<HeosState> {
        self.px_state(p, smolar, PxKey::Smolar)
    }

    fn px_state(&self, p: f64, value: f64, key: PxKey) -> Result<HeosState> {
        let pc = self.fluid().states.critical.p;
        let tc = self.fluid().states.critical.t;
        let rhoc = self.fluid().states.critical.rhomolar;
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
                });
            }
        } else {
            return Err(Error::NotImplemented(
                "(p,X) flash below the triple-point pressure is not ported yet".into(),
            ));
        };

        // Single-phase: bracket in T per upstream `HSU_P_flash`, then solve.
        let t_min_fluid = self.fluid().eos.t_triple; // Tmin == Ttriple for the ported fluids
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
            Phase::Liquid => {
                // Melting line deferred (PLAN.md): lower bound is Tmin - 1e-3.
                (t_min_fluid - 1e-3, self.sat().pq_flash(p, 0.0)?.t)
            }
            Phase::SupercriticalLiquid | Phase::SupercriticalGas | Phase::Supercritical => {
                (t_min_fluid - 1e-3, 1.5 * self.fluid().eos.t_max)
            }
            _ => return Err(Error::Value("Not a valid homogeneous state".into())),
        };

        // Residual: keyed value at the (p,T) state minus the target. For
        // liquid/gas the phase is imposed on the inner density solve
        // (upstream `specify_phase`); supercritical probes re-determine.
        let probe = |t: f64| -> Result<f64> {
            let rho = self.px_probe_rho(t, p, phase)?;
            Ok(self.px_value(key, t, rho))
        };
        let r_min = probe(t_min)? - value;
        let r_max = probe(t_max)? - value;
        if r_min * r_max >= 0.0 {
            return Err(Error::Solution(format!(
                "unable to bracket the (p,X) solution in [{t_min}, {t_max}] (residuals {r_min:e}, {r_max:e}); the derivative path is not ported"
            )));
        }
        // Bisection at upstream's 30-bit relative tolerance, then the bracket
        // midpoint (upstream re-evaluates the state there).
        let tol = (2.0f64).powi(1 - 30);
        let (mut a, mut b) = (t_min, t_max);
        let mut fa = r_min;
        for _ in 0..200 {
            if (b - a) <= tol * a.abs().max(b.abs()) {
                break;
            }
            let m = 0.5 * (a + b);
            let fm = probe(m)? - value;
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
        let t = 0.5 * (a + b);
        let rho = self.px_probe_rho(t, p, phase)?;
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar: rho,
            phase,
            q: -1.0,
        })
    }

    fn px_value(&self, key: PxKey, t: f64, rhomolar: f64) -> f64 {
        match key {
            PxKey::Hmolar => self.eos.hmolar(t, rhomolar),
            PxKey::Smolar => self.eos.smolar(t, rhomolar),
        }
    }

    /// Inner density for a (p,T) probe: imposed phase for liquid/gas,
    /// full determination for the supercritical classifications (which can
    /// legitimately flip along the bracket).
    fn px_probe_rho(&self, t: f64, p: f64, phase: Phase) -> Result<f64> {
        match phase {
            Phase::Liquid | Phase::Gas => self.solver_rho_tp(t, p, phase),
            _ => {
                let (rho, _phase) = self.pt_flash(t, p)?;
                Ok(rho)
            }
        }
    }

    /// (Dmolar, P) flash — upstream `DP_flash`: superancillary (p, Dmolar)
    /// phase determination, then a Halley solve for T at fixed density
    /// (Peng-Robinson seed for gas-like phases, saturation/1.1*Tc seeds
    /// otherwise), with the 30-bit bracketed fallback.
    pub fn dmolar_p_state(&self, rhomolar: f64, p: f64) -> Result<HeosState> {
        let pc = self.fluid().states.critical.p;
        let tc = self.fluid().states.critical.t;
        let rhoc = self.fluid().states.critical.rhomolar;
        let p_triple = self.fluid().eos.sat_min_liquid.p;

        let (phase, t0) = if p > pc {
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
        Ok(HeosState::SinglePhase {
            t,
            p,
            rhomolar,
            phase,
            q: -1.0,
        })
    }

    /// Upstream `T_DP_PengRobinson`: PR-based T seed at fixed (rho, p).
    fn t_dp_peng_robinson(&self, rhomolar: f64, p: f64) -> f64 {
        let omega = self.fluid().eos.acentric;
        let tc = self.fluid().states.critical.t;
        let pc = self.fluid().states.critical.p;
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
        let t_lo = self.fluid().eos.t_triple;
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
enum PxKey {
    Hmolar,
    Smolar,
}
