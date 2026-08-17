//! Slice 10f: the mixture stability machinery and the full stability-tested
//! PT flash — upstream `StabilityRoutines::StabilityEvaluationClass`
//! (Michelsen TPD, the default algorithm; `check_stability_legacy` is the
//! non-default `MIXTURE_STABILITY_ALGORITHM=0` path and stays unported),
//! `guess_split_from_wilson`, `PTflash_twophase::solve_michelsen`
//! (`solve_legacy` likewise unported), the global lowest-Gibbs density
//! solver, and `PT_flash_mixtures`' full two-phase glue.

#![allow(clippy::needless_range_loop)]

use crate::alpha::{DerivsMemo, HelmholtzDerivs};
use crate::mixture::MixtureModel;
use crate::mixture::XnFlag;
use crate::mixture_flash::MixtureState;
use crate::mixture_vle::{
    SatState, normalize_vector, rachford_rice_beta_bisect, solve_linear, wilson_ln_k_factor,
    x_and_y_from_k,
};
use crate::solvers::{Resid1D, brent, halley_omega};
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};

/// `SRK_covolume` (upstream hardcodes this R, and uses the CRITICAL states —
/// unlike `solver_rho_Tp_SRK`, which uses the reducing states).
fn srk_covolume(model: &MixtureModel, x: &[f64]) -> f64 {
    let mut b = 0.0;
    for i in 0..x.len() {
        let tc = model.crit_t()[i];
        let pc = model.crit_p()[i];
        let r = 8.3144598;
        #[allow(clippy::excessive_precision)]
        let omega_b = 0.08664034999649577215890158147700;
        b += x[i] * omega_b * r * tc / pc;
    }
    b
}

/// `calc_rhomolar_max_bound`.
pub(crate) fn rhomolar_max_bound(model: &MixtureModel, x: &[f64]) -> f64 {
    0.9 / srk_covolume(model, x)
}

/// dp/drho residual with derivatives (upstream `dpdrho_resid`).
struct DpdrhoResid<'a> {
    model: &'a MixtureModel,
    x: &'a [f64],
    t: f64,
    rhor: f64,
    tau: f64,
    delta: f64,
    /// Pressure at the last `call` point (the slow heavy-side probe gates on
    /// upstream's `this->p()` after the state update inside `call`).
    p_last: f64,
    memo: DerivsMemo,
}

impl DpdrhoResid<'_> {
    /// The mixture alphar matrix at (self.tau, delta), computed once per
    /// point.
    fn ar(&mut self, delta: f64) -> HelmholtzDerivs {
        let (model, x, tau) = (self.model, self.x, self.tau);
        self.memo
            .get_or_compute(tau, delta, |tau, delta| model.alphar_all(x, tau, delta))
    }
}

impl Resid1D for DpdrhoResid<'_> {
    fn call(&mut self, rhomolar: f64) -> f64 {
        self.delta = rhomolar / self.rhor;
        let d = self.ar(self.delta);
        self.p_last = rhomolar * self.model.gas_constant() * self.t * (1.0 + self.delta * d.d10);
        self.model.gas_constant()
            * self.t
            * (1.0 + 2.0 * self.delta * d.d10 + self.delta * self.delta * d.d20)
    }
    fn deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self.ar(self.delta);
        self.model.gas_constant() * self.t / self.rhor
            * (2.0 * d.d10 + 4.0 * self.delta * d.d20 + self.delta * self.delta * d.d30)
    }
    fn second_deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self.ar(self.delta);
        self.model.gas_constant() * self.t / (self.rhor * self.rhor)
            * (6.0 * d.d20 + 6.0 * self.delta * d.d30 + self.delta * self.delta * d.d40)
    }
    fn third_deriv(&mut self, _rhomolar: f64) -> f64 {
        unreachable!("Halley never calls third_deriv")
    }
}

enum StationaryPoints {
    Zero,
    One,
    Two { light: f64, heavy: f64 },
}

/// `solver_dpdrho0_Tp`: locate the two densities where dp/drho|T = 0.
fn solver_dpdrho0_tp(
    model: &MixtureModel,
    x: &[f64],
    t: f64,
    rhomax: f64,
) -> Result<StationaryPoints> {
    let rhor = model.reducing.rhormolar(x);
    let tau = model.reducing.tr(x) / t;
    let mut resid = DpdrhoResid {
        model,
        x,
        t,
        rhor,
        tau,
        delta: f64::NAN,
        p_last: f64::NAN,
        memo: DerivsMemo::default(),
    };

    let mut light = -1.0;
    let attempt = halley_omega(&mut resid, 1e-6, 1e-8, 100, 1.0).and_then(|l| {
        if resid.deriv(l) > 0.0 {
            Err(Error::Value("curvature cannot be positive".into()))
        } else {
            Ok(l)
        }
    });
    if let Ok(l) = attempt {
        light = l;
    }
    if light < 0.0 {
        // Slow path: increase density until curvature is positive
        let mut rho = 1e-6;
        for _ in 0..=100 {
            resid.call(rho);
            let curvature = resid.deriv(rho);
            if curvature > 0.0 {
                light = rho;
                break;
            }
            rho *= 2.0;
        }
    }

    let mut heavy = -1.0;
    let mut omega = 0.7;
    while omega > 0.0 {
        let attempt = halley_omega(&mut resid, rhomax, 1e-8, 100, omega).and_then(|h| {
            if resid.deriv(h) < 0.0 {
                Err(Error::Value("curvature cannot be negative".into()))
            } else {
                Ok(h)
            }
        });
        if let Ok(h) = attempt {
            heavy = h;
            break;
        }
        heavy = -1.0;
        omega -= 0.2;
    }
    if heavy < 0.0 {
        // Slow path: decrease density until curvature or pressure goes negative
        let mut rho = rhomax;
        for _ in 0..=100 {
            resid.call(rho);
            let curvature = resid.deriv(rho);
            if curvature < 0.0 || resid.p_last < 0.0 {
                heavy = rho;
                break;
            }
            rho /= 1.1;
        }
    }

    if light > 0.0 && heavy > 0.0 {
        Ok(StationaryPoints::Two { light, heavy })
    } else if light < 0.0 && heavy < 0.0 {
        let dpdrho_min = resid.call(1e-10);
        let dpdrho_max = resid.call(rhomax);
        if dpdrho_max * dpdrho_min > 0.0 {
            Ok(StationaryPoints::Zero)
        } else {
            Err(Error::Value(
                "zero stationary points -- does this make sense?".into(),
            ))
        }
    } else {
        Ok(StationaryPoints::One)
    }
}

/// `solver_rho_Tp_global`: the lowest-Gibbs density root at (T, p).
pub(crate) fn solver_rho_tp_global(
    model: &MixtureModel,
    x: &[f64],
    t: f64,
    p: f64,
    mut rhomolar_max: f64,
) -> Result<f64> {
    let retval = solver_dpdrho0_tp(model, x, t, rhomolar_max)?;
    let presid = |rho: f64| (model.pressure(x, t, rho) - p) / p;

    match retval {
        StationaryPoints::Zero => brent(presid, 1e-10, rhomolar_max, f64::EPSILON, 1e-8, 100),
        StationaryPoints::Two { light, heavy } => {
            let p_at_rhomin_stationary = model.pressure(x, t, light);
            let p_at_rhomax_stationary = model.pressure(x, t, heavy);

            let mut rho_liq = -1.0;
            let mut rho_vap = -1.0;
            if p > p_at_rhomax_stationary {
                for _ in 0..=10 {
                    let p_at_rhomax = model.pressure(x, t, rhomolar_max);
                    if p_at_rhomax < p {
                        rhomolar_max *= 1.05;
                    } else {
                        break;
                    }
                }
                rho_liq =
                    brent(presid, heavy, rhomolar_max, f64::EPSILON, 1e-8, 100).unwrap_or(-1.0);
            }
            if p < p_at_rhomin_stationary {
                rho_vap = brent(presid, light, 1e-10, f64::EPSILON, 1e-8, 100).unwrap_or(-1.0);
            }

            if rho_vap > 0.0 && rho_liq > 0.0 {
                if (rho_vap - rho_liq).abs() < 1e-10 {
                    Ok(rho_vap)
                } else {
                    let gibbs_vap = model.gibbsmolar_nocache(x, t, rho_vap);
                    let gibbs_liq = model.gibbsmolar_nocache(x, t, rho_liq);
                    if gibbs_liq < gibbs_vap {
                        Ok(rho_liq)
                    } else {
                        Ok(rho_vap)
                    }
                }
            } else if rho_vap < 0.0 && rho_liq > 0.0 {
                Ok(rho_liq)
            } else if rho_vap > 0.0 && rho_liq < 0.0 {
                Ok(rho_vap)
            } else {
                Err(Error::Value(format!(
                    "No density solutions for T={t},p={p}"
                )))
            }
        }
        StationaryPoints::One => Err(Error::Value(format!(
            "One stationary point (not good) for T={t},p={p}"
        ))),
    }
}

/// `solve_trial_rho_warm`: warm-started trial-phase density at (T, p) with
/// branch-jump rejection and the global-solver fallback. Updates the state.
/// `phase` is the SatL/SatV constructor imposition (liquid for the satl
/// instance, gas for satv) that arms the warm solve's mechanical-stability
/// retry.
fn solve_trial_rho_warm(
    sat: &mut SatState<'_>,
    t: f64,
    p: f64,
    rho_warm: &mut f64,
    phase: Phase,
) -> Result<f64> {
    if *rho_warm > 0.0 {
        let warm = sat.update_tp_guessrho_result(t, p, *rho_warm, phase);
        if let Ok(r) = warm {
            if r.is_finite() && r > 0.0 && r < 2.0 * *rho_warm && r > 0.5 * *rho_warm {
                *rho_warm = r;
                return Ok(r);
            }
        }
    }
    let model = sat.model();
    let x = sat.x.clone();
    let rg = solver_rho_tp_global(model, &x, t, p, rhomolar_max_bound(model, &x))?;
    sat.set_state_public(t, rg);
    *rho_warm = rg;
    Ok(rg)
}

// ---------------------------------------------------------------------------
// Stability evaluation (Michelsen TPD)
// ---------------------------------------------------------------------------

/// The verdict of `StabilityEvaluationClass`: stable/unstable plus the trial
/// phases when unstable, and the non-conclusive flag.
pub struct StabilityVerdict {
    pub stable: bool,
    pub uncertain: bool,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub rhomolar_liq: f64,
    pub rhomolar_vap: f64,
}

/// `check_stability_michelsen` over one feed z at (T, p).
pub fn check_stability_michelsen(
    model: &MixtureModel,
    z: &[f64],
    the_t: f64,
    the_p: f64,
) -> Result<StabilityVerdict> {
    let n = z.len();
    let mut verdict = StabilityVerdict {
        stable: true,
        uncertain: false,
        x: z.to_vec(),
        y: z.to_vec(),
        rhomolar_liq: -1.0,
        rhomolar_vap: -1.0,
    };
    let mut any_uncertain = false;

    // Feed fugacities d_i = ln(z_i) + ln(phi_i(z))
    let mut satl = SatState::new(model, z.to_vec());
    let rho_b = match solver_rho_tp_global(model, z, the_t, the_p, rhomolar_max_bound(model, z)) {
        Ok(r) => r,
        Err(_) => {
            // Between the spinodal pressures: SRK-seeded gas-phase fallback.
            model.solver_rho_tp(z, the_t, the_p, Phase::Gas)?
        }
    };
    satl.set_state_public(the_t, rho_b);
    let mut ln_f_z = vec![0.0; n];
    for i in 0..n {
        if z[i] > 0.0 {
            ln_f_z[i] = z[i].ln() + satl.fugacity_coefficient(i).ln();
        } else {
            ln_f_z[i] = -1e30;
        }
    }

    // Wilson trial compositions (Michelsen 1982a Eq. 28)
    let mut y_v = vec![0.0; n];
    let mut x_l = vec![0.0; n];
    for i in 0..n {
        let ki = wilson_ln_k_factor(model, the_t, the_p, i).exp();
        y_v[i] = z[i] * ki;
        x_l[i] = z[i] / ki;
    }

    let mut satv = SatState::new(model, z.to_vec());
    let trials = [y_v, x_l];
    for (t_idx, trial) in trials.into_iter().enumerate() {
        let mut y_cap = trial;
        let mut rho_warm = -1.0;

        // --- Phase 1: SS with GDEM acceleration ---
        let max_ss_loops = 4;
        let cntol = 1e-7;
        let mut ss_decided = false;

        'ss: for _loop in 0..max_ss_loops {
            if ss_decided {
                break;
            }
            let mut esq_pair = [0.0_f64; 2];
            let mut err = vec![0.0; n];

            for kk in 0..2 {
                if ss_decided {
                    break;
                }
                let sum_y: f64 = y_cap.iter().sum();
                let y_norm: Vec<f64> = y_cap.iter().map(|v| v / sum_y).collect();

                satv.set_mole_fractions(&y_norm);
                // Upstream runs every trial direction through SatV, whose
                // constructor imposition is gas.
                if solve_trial_rho_warm(&mut satv, the_t, the_p, &mut rho_warm, Phase::Gas).is_err()
                {
                    ss_decided = true;
                    break;
                }

                let mut tm = 1.0;
                let mut gmax = 0.0_f64;
                let mut esq = 0.0;
                for i in 0..n {
                    let ln_phi_y = satv.fugacity_coefficient(i).ln();
                    let ln_y_new = ln_f_z[i] - ln_phi_y;
                    let ln_y_old = y_cap[i].max(1e-300).ln();
                    let diff = ln_y_new - ln_y_old;
                    err[i] = diff;
                    esq += y_cap[i] * diff * diff;
                    gmax = gmax.max(diff.abs());

                    let s_i = ln_y_old + ln_phi_y - ln_f_z[i];
                    tm += y_cap[i] * (s_i - 1.0);

                    y_cap[i] = ln_y_new.exp();
                }

                // Early exit: tm < 0 means unstable
                if tm < -cntol {
                    let s_y: f64 = y_cap.iter().sum();
                    let y_norm: Vec<f64> = y_cap.iter().map(|v| v / s_y).collect();
                    if t_idx == 0 {
                        verdict.y = y_norm;
                        verdict.x = z.to_vec();
                    } else {
                        verdict.x = y_norm;
                        verdict.y = z.to_vec();
                    }
                    verdict.stable = false;
                    return Ok(verdict);
                }
                if gmax < cntol {
                    ss_decided = true;
                    break;
                }
                // Proximity test
                let mut distance_sq = 0.0;
                let mut curvature = 0.0;
                for i in 0..n {
                    let zysq = (y_cap[i] * z[i]).sqrt();
                    distance_sq += y_cap[i] + z[i] - 2.0 * zysq;
                    curvature -= (y_cap[i] - zysq) * err[i];
                }
                if distance_sq < 0.0 {
                    distance_sq = -distance_sq;
                }
                if distance_sq.sqrt() < 0.1 && curvature > 0.0 && tm / curvature > 0.8 {
                    ss_decided = true;
                    break;
                }
                esq_pair[kk] = esq;
            }
            if ss_decided {
                break 'ss;
            }
            // GDEM extrapolation
            if esq_pair[0] > 0.0 {
                let mut ratio = (esq_pair[1] / esq_pair[0]).sqrt();
                if !ratio.is_finite() || !(0.0..0.95).contains(&ratio) {
                    ratio = 0.95;
                }
                let factor = ratio / (1.0 - ratio);
                for i in 0..n {
                    let ln_y = y_cap[i].max(1e-300).ln() + factor * err[i];
                    y_cap[i] = ln_y.exp();
                }
            }
        }

        // --- Phase 2: second-order TPD minimization ---
        let mut trial_unstable = false;
        let trial_ok = minimize_tpd(
            model,
            &mut satv,
            z,
            &mut y_cap,
            &ln_f_z,
            the_t,
            the_p,
            &mut trial_unstable,
        );
        if !trial_ok {
            any_uncertain = true;
        }
        if trial_ok && trial_unstable {
            let s_y: f64 = y_cap.iter().sum();
            let y_norm: Vec<f64> = y_cap.iter().map(|v| v / s_y).collect();
            if t_idx == 0 {
                verdict.y = y_norm;
                verdict.x = z.to_vec();
            } else {
                verdict.x = y_norm;
                verdict.y = z.to_vec();
            }
            verdict.stable = false;
            return Ok(verdict);
        }
    }
    verdict.uncertain = any_uncertain;
    verdict.stable = true;
    Ok(verdict)
}

/// `minimize_tpd`: trust-region quasi-Newton TPD minimization in alpha
/// variables (`alpha_i = 2 sqrt(Y_i)`). Returns false when non-conclusive.
#[allow(clippy::too_many_arguments)]
fn minimize_tpd(
    _model: &MixtureModel,
    satv: &mut SatState<'_>,
    z: &[f64],
    y_cap: &mut [f64],
    ln_f_z: &[f64],
    the_t: f64,
    the_p: f64,
    is_unstable: &mut bool,
) -> bool {
    let n = y_cap.len();
    let cntol = 1e-7;
    let max_iter = 20;
    *is_unstable = false;

    let mut alpha: Vec<f64> = y_cap.iter().map(|y| 2.0 * y.max(1e-300).sqrt()).collect();
    let mut alpha_old = vec![0.0; n];
    let mut trust_radius = 0.25;
    let mut diagonal_shift = 0.0;
    let mut rho_warm = -1.0;

    for _iter in 0..max_iter {
        for i in 0..n {
            y_cap[i] = (0.5 * alpha[i]) * (0.5 * alpha[i]);
            alpha_old[i] = alpha[i];
        }
        let sum_y: f64 = y_cap.iter().sum();
        let y_norm: Vec<f64> = y_cap.iter().map(|v| v / sum_y).collect();

        satv.set_mole_fractions(&y_norm);
        if solve_trial_rho_warm(satv, the_t, the_p, &mut rho_warm, Phase::Gas).is_err() {
            return false;
        }

        let mut scaled_grad = vec![0.0; n];
        let mut grad = vec![0.0; n];
        let mut half_alpha = vec![0.0; n];
        let mut max_gradient = 0.0_f64;
        let mut obj_old = 1.0;
        for i in 0..n {
            half_alpha[i] = alpha[i] * 0.5;
            if z[i] > 0.0 {
                let ln_y = y_cap[i].max(1e-300).ln();
                let ln_phi = satv.fugacity_coefficient(i).ln();
                scaled_grad[i] = ln_y + ln_phi - ln_f_z[i];
            } else {
                scaled_grad[i] = 0.0;
            }
            grad[i] = -scaled_grad[i] * half_alpha[i];
            max_gradient = max_gradient.max(grad[i].abs());
            obj_old += y_cap[i] * (scaled_grad[i] - 1.0);
        }
        if max_gradient < cntol {
            *is_unstable = obj_old < -cntol;
            return true;
        }

        // Hessian (alpha variables); dln_phi/dn via XN_INDEPENDENT
        let mut h = vec![vec![0.0; n]; n];
        for i in 0..n {
            let ahi = half_alpha[i] / sum_y;
            for j in 0..=i {
                let dln_phi_dnj =
                    satv.dln_fugacity_dxj__const_t_p_xi_pub(i, j, XnFlag::Independent);
                let term = ahi * half_alpha[j] * dln_phi_dnj;
                h[i][j] = term;
                h[j][i] = term;
            }
            h[i][i] += 1.0 + 0.5 * scaled_grad[i];
        }

        let max_inner = 20;
        let mut step_accepted = false;
        for _inner in 0..max_inner {
            let mut h_shifted = h.clone();
            for i in 0..n {
                h_shifted[i][i] += diagonal_shift;
            }
            let delta_alpha_v =
                match solve_linear(&mut h_shifted, &grad.iter().map(|g| -g).collect::<Vec<_>>()) {
                    Ok(v) => v,
                    Err(_) => break,
                };
            // solve_linear solves J v = -r with r as passed; upstream solves
            // (H+shift) delta = -grad, so pass r = grad (solve_linear negates).
            let mut delta_alpha = delta_alpha_v;
            let mut step_norm_sq = 0.0;
            for i in 0..n {
                let mut da = delta_alpha[i];
                if da + alpha_old[i] <= 0.0 {
                    da = -0.9 * alpha_old[i];
                }
                delta_alpha[i] = da;
                alpha[i] = alpha_old[i] + da;
                y_cap[i] = (0.5 * alpha[i]) * (0.5 * alpha[i]);
                step_norm_sq += da * da;
            }
            let step_size = step_norm_sq.sqrt();

            if step_size > trust_radius && diagonal_shift == 0.0 {
                diagonal_shift = step_size / trust_radius - 1.0;
                for i in 0..n {
                    alpha[i] = alpha_old[i];
                    y_cap[i] = (0.5 * alpha[i]) * (0.5 * alpha[i]);
                }
                continue;
            }

            let sum_y2: f64 = y_cap.iter().sum();
            let y_norm2: Vec<f64> = y_cap.iter().map(|v| v / sum_y2).collect();
            satv.set_mole_fractions(&y_norm2);
            if solve_trial_rho_warm(satv, the_t, the_p, &mut rho_warm, Phase::Gas).is_err() {
                trust_radius = step_size / 3.0;
                diagonal_shift = 0.0;
                for i in 0..n {
                    alpha[i] = alpha_old[i];
                    y_cap[i] = (0.5 * alpha[i]) * (0.5 * alpha[i]);
                }
                continue;
            }

            let mut obj_new = 1.0;
            for i in 0..n {
                let ln_y = y_cap[i].max(1e-300).ln();
                let ln_phi = satv.fugacity_coefficient(i).ln();
                obj_new += y_cap[i] * (ln_y + ln_phi - ln_f_z[i] - 1.0);
            }
            if obj_new < -cntol {
                *is_unstable = true;
                return true;
            }
            if obj_new > obj_old + 1e-12 {
                trust_radius = step_size / 3.0;
                diagonal_shift = 0.0;
                for i in 0..n {
                    alpha[i] = alpha_old[i];
                    y_cap[i] = (0.5 * alpha[i]) * (0.5 * alpha[i]);
                }
                continue;
            }

            // Trust-region update from actual vs predicted reduction
            let mut hd = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    hd[i] += h[i][j] * delta_alpha[j];
                }
            }
            let mut predicted = 0.0;
            for i in 0..n {
                predicted += 0.5 * delta_alpha[i] * hd[i] - delta_alpha[i] * grad[i];
            }
            let actual = obj_new - obj_old;
            let ratio = if predicted != 0.0 {
                actual / predicted
            } else {
                1.0
            };
            if ratio < 0.25 {
                trust_radius = step_size / 2.0;
            } else if ratio > 0.75 && diagonal_shift > 0.0 {
                trust_radius = step_size * 2.0;
            } else {
                trust_radius = step_size;
            }
            diagonal_shift = 0.0;
            step_accepted = true;
            break;
        }
        if !step_accepted {
            return false;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Wilson-seeded speculative split
// ---------------------------------------------------------------------------

/// `successive_substitution_guessrho`. A density-solve failure PROPAGATES —
/// upstream's `update_TP_guessrho` calls here are unguarded
/// (`VLERoutines.cpp:1845-1847`), so the throw escapes the whole Wilson
/// guess block and the PT flash's catch treats it as not-two-phase.
/// Swallowing it instead kept the PARTIAL seed alive, and `solve_michelsen`
/// could converge it to a degenerate one-component-per-phase "equilibrium"
/// that slips through the verify's y_i < 1e-12 skip (found on
/// Methane[0.5]&Ethane[0.5] at 225 K / 8.75 MPa, where the wheel stays
/// single-phase liquid).
#[allow(clippy::too_many_arguments)]
fn successive_substitution_guessrho(
    satl: &mut SatState<'_>,
    satv: &mut SatState<'_>,
    x: &mut [f64],
    y: &mut [f64],
    rhomolar_liq: &mut f64,
    rhomolar_vap: &mut f64,
    z: &[f64],
    t: f64,
    p: f64,
    num_steps: i32,
) -> Result<()> {
    let n = z.len();
    let tol = 1e-10;
    let mut ln_k = vec![0.0; n];
    let mut k = vec![0.0; n];
    for ss in 0..num_steps {
        let rho_liq_prev = *rhomolar_liq;
        let rho_vap_prev = *rhomolar_vap;
        satl.set_mole_fractions(x);
        satl.update_tp_guessrho(t, p, *rhomolar_liq, Phase::Liquid)?;
        satv.set_mole_fractions(y);
        satv.update_tp_guessrho(t, p, *rhomolar_vap, Phase::Gas)?;
        *rhomolar_liq = satl.rhomolar;
        *rhomolar_vap = satv.rhomolar;

        let mut g0 = 0.0;
        let mut g1 = 0.0;
        let mut max_ln_k_change = 0.0_f64;
        let mut finite = true;
        for i in 0..n {
            let ln_k_new = (satl.fugacity_coefficient(i) / satv.fugacity_coefficient(i)).ln();
            if !ln_k_new.is_finite() {
                finite = false;
                break;
            }
            max_ln_k_change = max_ln_k_change.max((ln_k_new - ln_k[i]).abs());
            ln_k[i] = ln_k_new;
            k[i] = ln_k[i].exp();
            if !k[i].is_finite() {
                finite = false;
                break;
            }
            g0 += z[i] * (k[i] - 1.0);
            g1 += z[i] * (1.0 - 1.0 / k[i]);
        }
        if !finite {
            *rhomolar_liq = rho_liq_prev;
            *rhomolar_vap = rho_vap_prev;
            break;
        }
        let beta = if g0 < 0.0 {
            0.0
        } else if g1 > 0.0 {
            1.0
        } else {
            rachford_rice_beta_bisect(z, &k)
        };
        x_and_y_from_k(beta, &k, z, x, y);
        normalize_vector(x);
        normalize_vector(y);
        if ss > 0 && max_ln_k_change < tol {
            break;
        }
    }
    Ok(())
}

/// `guess_split_from_wilson`: ideal K-factor split seed + SS refinement.
/// Returns None when no usable estimate exists.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn guess_split_from_wilson(
    model: &MixtureModel,
    satl: &mut SatState<'_>,
    satv: &mut SatState<'_>,
    z: &[f64],
    t: f64,
    p: f64,
    num_steps: i32,
    require_bracket: bool,
) -> Option<(Vec<f64>, Vec<f64>, f64, f64)> {
    let n = z.len();
    let mut k = vec![0.0; n];
    let mut g0 = 0.0;
    let mut g1 = 0.0;
    for i in 0..n {
        let ln_ki = wilson_ln_k_factor(model, t, p, i);
        if !ln_ki.is_finite() {
            return None;
        }
        k[i] = ln_ki.exp();
        if !k[i].is_finite() {
            return None;
        }
        g0 += z[i] * (k[i] - 1.0);
        g1 += z[i] * (1.0 - 1.0 / k[i]);
    }
    let bracketed = g0 > 0.0 && g1 < 0.0;
    if require_bracket && !bracketed {
        return None;
    }
    let beta = if bracketed {
        rachford_rice_beta_bisect(z, &k)
    } else {
        0.5
    };
    let mut x = vec![0.0; n];
    let mut y = vec![0.0; n];
    x_and_y_from_k(beta, &k, z, &mut x, &mut y);
    normalize_vector(&mut x);
    normalize_vector(&mut y);

    let mut rhomolar_liq =
        solver_rho_tp_global(model, &x, t, p, rhomolar_max_bound(model, &x)).ok()?;
    let mut rhomolar_vap =
        solver_rho_tp_global(model, &y, t, p, rhomolar_max_bound(model, &y)).ok()?;
    if !rhomolar_liq.is_finite()
        || rhomolar_liq <= 0.0
        || !rhomolar_vap.is_finite()
        || rhomolar_vap <= 0.0
    {
        return None;
    }
    successive_substitution_guessrho(
        satl,
        satv,
        &mut x,
        &mut y,
        &mut rhomolar_liq,
        &mut rhomolar_vap,
        z,
        t,
        p,
        num_steps,
    )
    .ok()?;
    Some((x, y, rhomolar_liq, rhomolar_vap))
}

// ---------------------------------------------------------------------------
// PTflash_twophase (Michelsen)
// ---------------------------------------------------------------------------

/// In/out block of `PTflash_twophase` (upstream `PTflash_twophase_options`).
pub struct PtFlashTwophase {
    pub t: f64,
    pub p: f64,
    pub z: Vec<f64>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub rhomolar_liq: f64,
    pub rhomolar_vap: f64,
    pub beta: f64,
    pub nonconvergence: bool,
}

/// Smallest eigenvalue of a symmetric matrix by cyclic Jacobi rotations
/// (stands in for Eigen's SelfAdjointEigenSolver; N is component count).
fn min_eigenvalue_symmetric(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a.to_vec();
    for _sweep in 0..50 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += m[i][j] * m[i][j];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if m[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
                let t_val = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t_val * t_val + 1.0).sqrt();
                let s = t_val * c;
                for k in 0..n {
                    let mkp = m[k][p];
                    let mkq = m[k][q];
                    m[k][p] = c * mkp - s * mkq;
                    m[k][q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let mpk = m[p][k];
                    let mqk = m[q][k];
                    m[p][k] = c * mpk - s * mqk;
                    m[q][k] = s * mpk + c * mqk;
                }
            }
        }
    }
    (0..n).map(|i| m[i][i]).fold(f64::INFINITY, f64::min)
}

/// `PTflash_twophase::solve_michelsen`.
pub fn ptflash_twophase_solve_michelsen(
    model: &MixtureModel,
    satl: &mut SatState<'_>,
    satv: &mut SatState<'_>,
    io: &mut PtFlashTwophase,
) -> Result<()> {
    let n = io.x.len();
    io.nonconvergence = false;
    let mut rho_warm_l = -1.0;
    let mut rho_warm_v = -1.0;

    let mut ln_k = vec![0.0; n];
    for i in 0..n {
        let ratio = io.y[i] / io.x[i].max(1e-300);
        ln_k[i] = ratio.max(1e-300).ln();
    }
    let mut beta = io.beta;

    for i in 0..n {
        if !ln_k[i].is_finite() || !io.x[i].is_finite() || !io.y[i].is_finite() {
            io.nonconvergence = true;
            return Err(Error::Solution(format!(
                "PTflash_twophase::solve_michelsen got a non-finite seed at T = {} K, p = {} Pa",
                io.t, io.p
            )));
        }
    }

    // Rachford-Rice in log-K space, Newton + bisection safeguards
    macro_rules! solve_rachford_rice {
        () => {{
            let mut beta_min = 0.0;
            let mut beta_max = 1.0;
            if !beta.is_finite() || !(0.0..=1.0).contains(&beta) {
                beta = 0.5;
            }
            for _rr_iter in 0..50 {
                let mut r = 0.0;
                let mut dr = 0.0;
                for i in 0..n {
                    let ki = ln_k[i].min(350.0).exp();
                    let term = ki - 1.0;
                    let denom = 1.0 + beta * term;
                    r += io.z[i] * term / denom;
                    dr -= io.z[i] * term * term / (denom * denom);
                }
                if r > 0.0 {
                    beta_min = beta;
                } else {
                    beta_max = beta;
                }
                if r.abs() < 1e-11 {
                    break;
                }
                let mut beta_new = beta - r / dr;
                if !beta_new.is_finite() || beta_new <= beta_min || beta_new >= beta_max {
                    beta_new = 0.5 * (beta_min + beta_max);
                }
                if (beta_new - beta).abs() < 1e-11 {
                    break;
                }
                beta = beta_new;
            }
            for i in 0..n {
                let ki = ln_k[i].min(350.0).exp();
                io.x[i] = io.z[i] / (1.0 + beta * (ki - 1.0));
                io.y[i] = ki * io.x[i];
            }
            normalize_vector(&mut io.x);
            normalize_vector(&mut io.y);
        }};
    }
    macro_rules! evaluate_phases {
        () => {{
            satl.set_mole_fractions(&io.x);
            let ok_l = solve_trial_rho_warm(satl, io.t, io.p, &mut rho_warm_l, Phase::Liquid);
            let mut ok = false;
            if let Ok(r) = ok_l {
                io.rhomolar_liq = r;
                satv.set_mole_fractions(&io.y);
                if let Ok(rv) = solve_trial_rho_warm(satv, io.t, io.p, &mut rho_warm_v, Phase::Gas)
                {
                    io.rhomolar_vap = rv;
                    ok = true;
                }
            }
            ok
        }};
    }

    // --- Phase 1: SS with GDEM ---
    let max_ss_loops = 4;
    let ss_tol = 1e-7;
    let mut ss_converged = false;
    for _loop in 0..max_ss_loops {
        if ss_converged {
            break;
        }
        let mut esq_pair = [0.0_f64; 2];
        let mut err = vec![0.0; n];
        for kk in 0..2 {
            if ss_converged {
                break;
            }
            solve_rachford_rice!();
            if !evaluate_phases!() {
                return Err(Error::Solution(
                    "PT flash lost a phase density solve during successive substitution".into(),
                ));
            }
            let mut gmax = 0.0_f64;
            let mut esq = 0.0;
            for i in 0..n {
                let ln_k_new =
                    satl.fugacity_coefficient(i).ln() - satv.fugacity_coefficient(i).ln();
                let diff = ln_k_new - ln_k[i];
                err[i] = diff;
                esq += io.z[i] * diff * diff;
                gmax = gmax.max(diff.abs());
                ln_k[i] = ln_k_new;
            }
            esq_pair[kk] = esq;
            if gmax < ss_tol {
                ss_converged = true;
            }
        }
        if ss_converged {
            break;
        }
        if esq_pair[0] > 0.0 {
            let mut ratio = (esq_pair[1] / esq_pair[0]).sqrt();
            if !ratio.is_finite() || !(0.0..0.95).contains(&ratio) {
                ratio = 0.95;
            }
            let factor = ratio / (1.0 - ratio);
            for i in 0..n {
                ln_k[i] += factor * err[i];
            }
        }
    }
    solve_rachford_rice!();
    if !evaluate_phases!() {
        return Err(Error::Solution(
            "PT flash lost a phase density solve after successive substitution".into(),
        ));
    }

    // --- Phase 2: second-order Gibbs minimization ---
    macro_rules! compute_gibbs {
        () => {{
            let mut g = 0.0;
            for i in 0..n {
                if io.x[i] > 0.0 {
                    g +=
                        (1.0 - beta) * io.x[i] * (io.x[i].ln() + satl.fugacity_coefficient(i).ln());
                }
                if io.y[i] > 0.0 {
                    g += beta * io.y[i] * (io.y[i].ln() + satv.fugacity_coefficient(i).ln());
                }
            }
            g
        }};
    }
    let mut g_old = compute_gibbs!();

    let gibbs_tol = 1e-9;
    let max_gibbs_iter = 50;
    let max_restart = 2;
    let max_inner = 40;
    let mut converged = false;

    for restart in 0..max_restart {
        if converged {
            break;
        }
        let mut trust_radius = if restart == 0 { 1.0 } else { 0.2 };

        for _gibbs_iter in 0..max_gibbs_iter {
            if converged {
                break;
            }
            let l_frac = 1.0 - beta;
            let v_frac = beta;

            satl.set_mole_fractions(&io.x);
            satl.set_state_public(io.t, io.rhomolar_liq);
            satv.set_mole_fractions(&io.y);
            satv.set_state_public(io.t, io.rhomolar_vap);

            let mut dl = vec![vec![0.0; n]; n];
            let mut dv = vec![vec![0.0; n]; n];
            for i in 0..n {
                for j in 0..n {
                    dl[i][j] = satl.dln_fugacity_dxj__const_t_p_xi_pub(i, j, XnFlag::Independent);
                    dv[i][j] = satv.dln_fugacity_dxj__const_t_p_xi_pub(i, j, XnFlag::Independent);
                }
            }
            let mut g = vec![0.0; n];
            let mut dia = vec![0.0; n];
            let mut h = vec![vec![0.0; n]; n];
            let mut max_g = 0.0_f64;
            for i in 0..n {
                let mut sum_x_dl = 0.0;
                let mut sum_y_dv = 0.0;
                for k in 0..n {
                    sum_x_dl += io.x[k] * dl[i][k];
                    sum_y_dv += io.y[k] * dv[i][k];
                }
                for j in 0..n {
                    let dln_phi_l_dnj = dl[i][j] - sum_x_dl;
                    let dln_phi_v_dnj = dv[i][j] - sum_y_dv;
                    h[i][j] = v_frac * dln_phi_l_dnj + l_frac * dln_phi_v_dnj - 1.0;
                }
                h[i][i] += (v_frac / io.x[i]) + (l_frac / io.y[i]);
                let l_act = io.x[i].ln() + satl.fugacity_coefficient(i).ln();
                let v_act = io.y[i].ln() + satv.fugacity_coefficient(i).ln();
                g[i] = v_frac * l_frac * (v_act - l_act);
                max_g = max_g.max((v_act - l_act).abs());
                let dia_i = (io.z[i] / (io.x[i] * io.y[i]).max(1e-300)).sqrt();
                dia[i] = if dia_i > 1e-300 { dia_i } else { 1.0 };
            }
            // (upstream also stores max_g into last_max_g here; the final
            // recompute below overwrites it unconditionally — dead store)
            if max_g < gibbs_tol {
                converged = true;
                break;
            }

            // Scaled system
            let mut gs = vec![0.0; n];
            let mut hs = vec![vec![0.0; n]; n];
            for i in 0..n {
                gs[i] = g[i] / dia[i];
                for j in 0..n {
                    hs[i][j] = h[i][j] / (dia[i] * dia[j]);
                }
            }

            let mut diagonal_shift = 0.0;
            let mut step_ok = false;
            for _inner in 0..max_inner {
                let mut hl = hs.clone();
                for i in 0..n {
                    hl[i][i] += diagonal_shift;
                }
                let min_eig = min_eigenvalue_symmetric(&hl);
                if !min_eig.is_finite() {
                    break;
                }
                if min_eig < 1e-8 {
                    diagonal_shift += 1e-8 - min_eig;
                    continue;
                }
                let ds = match solve_linear(&mut hl.clone(), &gs) {
                    Ok(v) => v,
                    Err(_) => break,
                };
                let step_norm = ds.iter().map(|v| v * v).sum::<f64>().sqrt();
                if step_norm > trust_radius {
                    diagonal_shift = if diagonal_shift > 0.0 {
                        diagonal_shift * 3.0
                    } else {
                        (step_norm / trust_radius - 1.0) * min_eig.max(1e-3)
                    };
                    continue;
                }

                let delta_v: Vec<f64> = (0..n).map(|i| ds[i] / dia[i]).collect();
                let mut pos_scale = 1.0_f64;
                for i in 0..n {
                    let v_old = beta * io.y[i];
                    if delta_v[i] > 0.0 && v_old + delta_v[i] > io.z[i] {
                        pos_scale = pos_scale.min(0.99 * (io.z[i] - v_old) / delta_v[i]);
                    }
                    if delta_v[i] < 0.0 && v_old + delta_v[i] < 0.0 {
                        pos_scale = pos_scale.min(0.99 * (-v_old) / delta_v[i]);
                    }
                }
                let mut v_new_total = 0.0;
                let mut l_new_total = 0.0;
                let mut v_new = vec![0.0; n];
                let mut l_new = vec![0.0; n];
                for i in 0..n {
                    v_new[i] = beta * io.y[i] + pos_scale * delta_v[i];
                    l_new[i] = io.z[i] - v_new[i];
                    v_new_total += v_new[i];
                    l_new_total += l_new[i];
                }
                let x_trial: Vec<f64> = (0..n).map(|i| l_new[i] / l_new_total).collect();
                let y_trial: Vec<f64> = (0..n).map(|i| v_new[i] / v_new_total).collect();

                let mut eval_ok = false;
                satl.set_mole_fractions(&x_trial);
                let rl = solver_rho_tp_global(
                    model,
                    &x_trial,
                    io.t,
                    io.p,
                    rhomolar_max_bound(model, &x_trial),
                );
                if let Ok(rl) = rl {
                    satv.set_mole_fractions(&y_trial);
                    let rv = solver_rho_tp_global(
                        model,
                        &y_trial,
                        io.t,
                        io.p,
                        rhomolar_max_bound(model, &y_trial),
                    );
                    if let Ok(rv) = rv {
                        if rl.is_finite() && rv.is_finite() && rl > 0.0 && rv > 0.0 {
                            satl.set_state_public(io.t, rl);
                            satv.set_state_public(io.t, rv);
                            let beta_save = beta;
                            let x_save = io.x.clone();
                            let y_save = io.y.clone();
                            let rl_save = io.rhomolar_liq;
                            let rv_save = io.rhomolar_vap;
                            beta = v_new_total;
                            io.x = x_trial.clone();
                            io.y = y_trial.clone();
                            io.rhomolar_liq = rl;
                            io.rhomolar_vap = rv;
                            let g_new = compute_gibbs!();
                            if g_new < g_old + 1e-12 {
                                g_old = g_new;
                                step_ok = true;
                                eval_ok = true;
                                if step_norm > 0.8 * trust_radius {
                                    trust_radius = (2.0 * trust_radius).min(1e3);
                                }
                            } else {
                                beta = beta_save;
                                io.x = x_save;
                                io.y = y_save;
                                io.rhomolar_liq = rl_save;
                                io.rhomolar_vap = rv_save;
                            }
                        }
                    }
                }
                if eval_ok {
                    break;
                }
                trust_radius = 0.5 * step_norm;
                if trust_radius < 1e-10 {
                    break;
                }
                diagonal_shift = 0.0;
            }
            if !step_ok {
                break;
            }
        }
    }
    io.beta = beta;

    // Recompute the residual on the final published state
    satl.set_mole_fractions(&io.x);
    satl.set_state_public(io.t, io.rhomolar_liq);
    satv.set_mole_fractions(&io.y);
    satv.set_state_public(io.t, io.rhomolar_vap);
    let mut final_max_g = 0.0_f64;
    for i in 0..n {
        let l_act = io.x[i].ln() + satl.fugacity_coefficient(i).ln();
        let v_act = io.y[i].ln() + satv.fugacity_coefficient(i).ln();
        final_max_g = final_max_g.max((v_act - l_act).abs());
    }
    let last_max_g = final_max_g;
    converged = final_max_g < gibbs_tol;

    if !converged {
        let mut spread = 0.0_f64;
        for i in 0..n {
            spread = spread.max((io.x[i] - io.y[i]).abs());
        }
        let genuine = last_max_g.is_finite()
            && last_max_g <= 1e-5
            && spread >= 1e-4
            && beta > 1e-8
            && beta < 1.0 - 1e-8;
        if !genuine {
            io.nonconvergence = true;
            return Err(Error::Solution(format!(
                "PTflash_twophase::solve_michelsen failed to converge: max|ln f_V - ln f_L| = {last_max_g:e} at T = {} K, p = {} Pa",
                io.t, io.p
            )));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// PT_flash_mixtures: full two-phase glue
// ---------------------------------------------------------------------------

/// What the full PT flash publishes.
pub enum PtFlashResult {
    Single(MixtureState),
    TwoPhase {
        t: f64,
        p: f64,
        q: f64,
        rhomolar: f64,
        x: Vec<f64>,
        y: Vec<f64>,
        rhomolar_liq: f64,
        rhomolar_vap: f64,
    },
}

/// Upstream `PT_flash_mixtures` without an imposed phase: stability test,
/// Wilson cross-check on a "stable" verdict, two-phase Michelsen split with
/// verify, trivial-split collapse, single-phase fallback.
pub fn pt_flash_mixtures(model: &MixtureModel, z: &[f64], t: f64, p: f64) -> Result<PtFlashResult> {
    let verdict = check_stability_michelsen(model, z, t, p)?;
    let mut do_twophase = !verdict.stable;
    let mut wilson_seeded = false;

    let mut satl = SatState::new(model, z.to_vec());
    let mut satv = SatState::new(model, z.to_vec());
    let mut io = PtFlashTwophase {
        t,
        p,
        z: z.to_vec(),
        x: verdict.x.clone(),
        y: verdict.y.clone(),
        rhomolar_liq: verdict.rhomolar_liq,
        rhomolar_vap: verdict.rhomolar_vap,
        beta: 0.5,
        nonconvergence: false,
    };

    if !do_twophase {
        // Cross-check the "stable" verdict with a speculative Wilson split
        // (forced when the verdict was non-conclusive).
        let require_bracket = !verdict.uncertain;
        if let Some((x, y, rho_l, rho_v)) =
            guess_split_from_wilson(model, &mut satl, &mut satv, z, t, p, 10, require_bracket)
        {
            io.x = x;
            io.y = y;
            io.rhomolar_liq = rho_l;
            io.rhomolar_vap = rho_v;
            do_twophase = true;
            wilson_seeded = true;
        }
    }

    if do_twophase {
        if wilson_seeded {
            // Speculative: any failure falls back to single phase; a
            // "successful" split must verify as a genuine equilibrium.
            if ptflash_twophase_solve_michelsen(model, &mut satl, &mut satv, &mut io).is_err() {
                do_twophase = false;
            }
            if do_twophase {
                let mut spread = 0.0_f64;
                for i in 0..z.len() {
                    spread = spread.max((io.x[i] - io.y[i]).abs());
                }
                satl.set_mole_fractions(&io.x);
                satl.set_state_public(io.t, io.rhomolar_liq);
                satv.set_mole_fractions(&io.y);
                satv.set_state_public(io.t, io.rhomolar_vap);
                let mut fug_resid = 0.0_f64;
                for i in 0..z.len() {
                    if io.x[i] < 1e-12 || io.y[i] < 1e-12 {
                        continue;
                    }
                    let lnf_l = io.x[i].ln() + satl.fugacity_coefficient(i).ln();
                    let lnf_v = io.y[i].ln() + satv.fugacity_coefficient(i).ln();
                    fug_resid = fug_resid.max((lnf_v - lnf_l).abs());
                }
                let ok = spread >= 1e-6 && fug_resid.is_finite() && fug_resid <= 1e-7;
                if !ok {
                    do_twophase = false;
                }
            }
        } else {
            // Genuine instability: nonconvergence falls back to single phase,
            // any other failure propagates.
            match ptflash_twophase_solve_michelsen(model, &mut satl, &mut satv, &mut io) {
                Ok(()) => {}
                Err(e) => {
                    if !io.nonconvergence {
                        return Err(e);
                    }
                    do_twophase = false;
                }
            }
        }
    }

    if do_twophase {
        if io.beta < 1e-10 {
            let phase = if io.rhomolar_liq < model.reducing.rhormolar(z) {
                Phase::Gas
            } else {
                Phase::Liquid
            };
            Ok(PtFlashResult::Single(MixtureState {
                t,
                p,
                rhomolar: io.rhomolar_liq,
                q: -1.0,
                phase,
            }))
        } else if io.beta > 1.0 - 1e-10 {
            let phase = if io.rhomolar_vap < model.reducing.rhormolar(z) {
                Phase::Gas
            } else {
                Phase::Liquid
            };
            Ok(PtFlashResult::Single(MixtureState {
                t,
                p,
                rhomolar: io.rhomolar_vap,
                q: -1.0,
                phase,
            }))
        } else {
            Ok(PtFlashResult::TwoPhase {
                t,
                p,
                q: io.beta,
                rhomolar: 1.0 / (io.beta / io.rhomolar_vap + (1.0 - io.beta) / io.rhomolar_liq),
                x: io.x,
                y: io.y,
                rhomolar_liq: io.rhomolar_liq,
                rhomolar_vap: io.rhomolar_vap,
            })
        }
    } else {
        model.pt_flash(z, t, p).map(PtFlashResult::Single)
    }
}
