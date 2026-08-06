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
struct CaloricTResid<'a> {
    eos: &'a HelmholtzEos,
    t: f64,
    target: f64,
    key: CaloricKey,
}

impl Resid1D for CaloricTResid<'_> {
    fn call(&mut self, rhomolar: f64) -> f64 {
        let v = match self.key {
            CaloricKey::Smolar => self.eos.smolar(self.t, rhomolar),
            CaloricKey::Hmolar => self.eos.hmolar(self.t, rhomolar),
            CaloricKey::Umolar => self.eos.umolar(self.t, rhomolar),
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
        let d = self.eos.alphar_all(tau, delta);
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
        let d = self.eos.alphar_all(tau, delta);
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
        } else if t > tc && t > self.t_triple() {
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
                });
            }
        } else {
            return Err(Error::NotImplemented(
                "(p,X) flash below the triple-point pressure is not ported yet".into(),
            ));
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

    fn px_value(&self, key: CaloricKey, t: f64, rhomolar: f64) -> f64 {
        match key {
            CaloricKey::Hmolar => self.eos.hmolar(t, rhomolar),
            CaloricKey::Smolar => self.eos.smolar(t, rhomolar),
            CaloricKey::Umolar => self.eos.umolar(t, rhomolar),
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
        let pc = self.p_critical();
        let tc = self.t_critical();
        let rhoc = self.rhomolar_critical();
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
