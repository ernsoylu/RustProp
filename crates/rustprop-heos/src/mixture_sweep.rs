//! Slice 10f part 2: the sweep-based mixture flashes — upstream
//! `DHSU_T_flash` / `HSU_P_flash` / `HSU_D_flash` mixture branches. Each
//! sweeps a (T or P) axis with the full stability-tested PT flash inside,
//! with the fast single-phase nocache paths tried first and upstream's
//! post-solve verification gates (never silently return a wrong state).
//! The Boost TOMS748 solver is ported verbatim (`solvers::toms748_solve`) —
//! plain bisection walked into wrong-root pockets of the discontinuous
//! residuals that the real algorithm's interpolating steps sail past.

use crate::mixture::MixtureModel;
use crate::mixture_stability::{PtFlashResult, check_stability_michelsen, pt_flash_mixtures};
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};

/// Which property a sweep matches (upstream `parameters other`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SweepVar {
    Dmolar,
    Hmolar,
    Smolar,
    Umolar,
}

impl SweepVar {
    fn short(self) -> &'static str {
        match self {
            SweepVar::Dmolar => "Dmolar",
            SweepVar::Hmolar => "Hmolar",
            SweepVar::Smolar => "Smolar",
            SweepVar::Umolar => "Umolar",
        }
    }
}

impl PtFlashResult {
    /// `keyed_output` of the published flash state for the sweep residuals.
    pub(crate) fn keyed(&self, model: &MixtureModel, z: &[f64], var: SweepVar) -> f64 {
        match self {
            PtFlashResult::Single(s) => match var {
                SweepVar::Dmolar => s.rhomolar,
                SweepVar::Hmolar => model.hmolar(z, s.t, s.rhomolar),
                SweepVar::Smolar => model.smolar(z, s.t, s.rhomolar),
                SweepVar::Umolar => model.umolar(z, s.t, s.rhomolar),
            },
            PtFlashResult::TwoPhase {
                t,
                q,
                rhomolar,
                x,
                y,
                rhomolar_liq,
                rhomolar_vap,
                ..
            } => match var {
                SweepVar::Dmolar => *rhomolar,
                SweepVar::Hmolar => {
                    q * model.hmolar(y, *t, *rhomolar_vap)
                        + (1.0 - q) * model.hmolar(x, *t, *rhomolar_liq)
                }
                SweepVar::Smolar => {
                    q * model.smolar(y, *t, *rhomolar_vap)
                        + (1.0 - q) * model.smolar(x, *t, *rhomolar_liq)
                }
                SweepVar::Umolar => {
                    q * model.umolar(y, *t, *rhomolar_vap)
                        + (1.0 - q) * model.umolar(x, *t, *rhomolar_liq)
                }
            },
        }
    }
}

/// `toms748_solve(..., eps_tolerance<double>(40), max_iter = 100)` as every
/// upstream sweep call site invokes it.
fn toms748_standin<F: FnMut(f64) -> Result<f64>>(
    f: &mut F,
    a: f64,
    b: f64,
    fa: f64,
    fb: f64,
) -> Result<f64> {
    crate::solvers::toms748_solve(f, a, b, fa, fb, 40, 100)
}

/// Mole-fraction-weighted limits (upstream `calc_Tmin`/`calc_Tmax`/
/// `calc_p_triple`/`calc_pmax`).
fn weighted_tmin(model: &MixtureModel, z: &[f64]) -> f64 {
    model.triple_t().iter().zip(z).map(|(t, x)| t * x).sum()
}
fn weighted_tmax(model: &MixtureModel, z: &[f64]) -> f64 {
    model.t_max().iter().zip(z).map(|(t, x)| t * x).sum()
}
fn weighted_p_triple(model: &MixtureModel, z: &[f64]) -> f64 {
    model.triple_p().iter().zip(z).map(|(p, x)| p * x).sum()
}
fn weighted_pmax(model: &MixtureModel, z: &[f64]) -> f64 {
    model.p_max().iter().zip(z).map(|(p, x)| p * x).sum()
}

fn nocache_value(model: &MixtureModel, z: &[f64], t: f64, rho: f64, var: SweepVar) -> f64 {
    match var {
        SweepVar::Dmolar => rho,
        SweepVar::Hmolar => model.hmolar(z, t, rho),
        SweepVar::Smolar => model.smolar(z, t, rho),
        SweepVar::Umolar => model.umolar(z, t, rho),
    }
}

/// Direct single-phase publish (`update_DmolarT_direct` + `_Q = -1` +
/// `recalculate_singlephase_phase`, mixture branch: reducing-density proxy).
fn direct_state(model: &MixtureModel, z: &[f64], t: f64, rho: f64) -> PtFlashResult {
    let phase = if rho > model.reducing.rhormolar(z) {
        Phase::Liquid
    } else {
        Phase::Gas
    };
    PtFlashResult::Single(crate::mixture_flash::MixtureState {
        t,
        p: model.pressure(z, t, rho),
        rhomolar: rho,
        q: -1.0,
        phase,
    })
}

/// The shared P-sweep at fixed T (upstream DHSU_T slow path): log scan in
/// half-decade steps over [max(p_triple, 100), pmax], TOMS748 refine, final
/// PT flash at the root.
fn p_sweep_at_t(
    model: &MixtureModel,
    z: &[f64],
    t: f64,
    var: SweepVar,
    value: f64,
) -> Result<PtFlashResult> {
    let mut pmin_bound = weighted_p_triple(model, z);
    let pmax_bound = weighted_pmax(model, z);
    if pmin_bound < 100.0 {
        pmin_bound = 100.0;
    }
    let mut p_resid = |p: f64| -> Result<f64> {
        let state = pt_flash_mixtures(model, z, t, p)?;
        Ok(state.keyed(model, z, var) - value)
    };

    let log_pmin = pmin_bound.log10();
    let log_pmax = pmax_bound.log10();
    let dlog_p = 0.5;
    let nsteps = ((log_pmax - log_pmin) / dlog_p) as i32 + 1;

    let mut p_lo = -1.0;
    let mut p_hi = -1.0;
    let mut f_lo = f64::NAN;
    let mut f_hi = f64::NAN;
    let mut f_prev = 0.0;
    let mut p_prev = 0.0;
    let mut have_prev = false;
    for i in 0..=nsteps {
        let mut log_p = log_pmin + f64::from(i) * dlog_p;
        if log_p > log_pmax {
            log_p = log_pmax;
        }
        let p = 10.0_f64.powf(log_p);
        let f = match p_resid(p) {
            Ok(f) => f,
            Err(_) => {
                have_prev = false;
                continue;
            }
        };
        if have_prev && f_prev * f < 0.0 {
            p_lo = p_prev;
            p_hi = p;
            f_lo = f_prev;
            f_hi = f;
            break;
        }
        p_prev = p;
        f_prev = f;
        have_prev = true;
    }

    if p_lo > 0.0 && p_hi > 0.0 {
        let p_sol = toms748_standin(&mut p_resid, p_lo, p_hi, f_lo, f_hi).map_err(|e| {
            Error::Value(format!(
                "DHSU_T_flash P-sweep TOMS748 for mixture failed: T={t}, target={}, value={value}, bracket=[{p_lo}, {p_hi}]: {e}",
                var.short()
            ))
        })?;
        pt_flash_mixtures(model, z, t, p_sol)
    } else {
        Err(Error::Value(format!(
            "DHSU_T_flash P-sweep for mixture: no bracket found scanning P in [{pmin_bound}, {pmax_bound}] at T={t} for target {}={value}",
            var.short()
        )))
    }
}

/// Upstream `DHSU_T_flash` mixture branch: (Dmolar|Hmolar|Smolar|Umolar, T).
pub fn dhsu_t_flash_mixtures(
    model: &MixtureModel,
    z: &[f64],
    t: f64,
    var: SweepVar,
    value: f64,
) -> Result<PtFlashResult> {
    let mut solved: Option<PtFlashResult> = None;

    if var == SweepVar::Dmolar {
        // Fast path: EOS pressure at (T, rho); accept if stability says stable.
        let p_eos = model.pressure(z, t, value);
        if p_eos > 0.0 {
            if let Ok(verdict) = check_stability_michelsen(model, z, t, p_eos) {
                if verdict.stable {
                    solved = Some(direct_state(model, z, t, value));
                }
            }
        }
    } else {
        // Fast path: nocache density sweeps in the gas and liquid ranges.
        let rho_min = 1e-10;
        let rho_max = crate::mixture_stability::rhomolar_max_bound(model, z);
        let rho_reducing = model.reducing.rhormolar(z);

        let is_mechanically_stable = |rho: f64| -> bool {
            let eps_rho = (rho * 1e-6).max(1e-10);
            let p_lo = model.pressure(z, t, rho - eps_rho);
            let p_hi = model.pressure(z, t, rho + eps_rho);
            (p_hi - p_lo) > 0.0
        };
        let mut rho_resid =
            |rho: f64| -> Result<f64> { Ok(nocache_value(model, z, t, rho, var) - value) };

        let mut bracketed_root = |a: f64, b: f64| -> Option<f64> {
            let fa = rho_resid(a).ok()?;
            let fb = rho_resid(b).ok()?;
            if fa * fb > 0.0 {
                return None;
            }
            toms748_standin(&mut rho_resid, a, b, fa, fb).ok()
        };

        let mut rho_gas = -1.0;
        if let Some(rho_cand) = bracketed_root(rho_min, rho_reducing) {
            if model.pressure(z, t, rho_cand) > 0.0 && is_mechanically_stable(rho_cand) {
                rho_gas = rho_cand;
            }
        }
        let mut rho_liq = -1.0;
        if let Some(rho_cand) = bracketed_root(rho_reducing, rho_max) {
            if model.pressure(z, t, rho_cand) > 0.0 && is_mechanically_stable(rho_cand) {
                rho_liq = rho_cand;
            }
        }
        // Fallback: locate the liquid spinodal via the p > 0 crossing
        if rho_liq < 0.0 {
            let mut rho_neg = -1.0;
            let mut rho_pos = -1.0;
            let mut rho_prev = rho_max;
            let mut p_prev = model.pressure(z, t, rho_prev);
            let mut rho_scan = rho_max * 0.95;
            // The 0.95 decay reaches any positive `rho_reducing` in
            // ln(rho_max / rho_reducing) / ln(1/0.95) steps — under 300 even
            // for a millionfold span. A non-positive reducing density is
            // never produced by real mixture data, but if one ever arrived
            // the loop would only stop when `rho_scan` underflowed to zero,
            // ~33,000 full Helmholtz pressure evaluations later. The cap is
            // two orders of magnitude above anything reachable, so it cannot
            // change a scan that currently terminates.
            let mut steps = 0u32;
            while rho_scan > rho_reducing && steps < 10_000 {
                steps += 1;
                let p_scan = model.pressure(z, t, rho_scan);
                if p_prev > 0.0 && p_scan < 0.0 {
                    rho_pos = rho_prev;
                    rho_neg = rho_scan;
                    break;
                }
                rho_prev = rho_scan;
                p_prev = p_scan;
                rho_scan *= 0.95;
            }
            if rho_neg > 0.0 && rho_pos > 0.0 {
                for _ in 0..30 {
                    let rho_mid = 0.5 * (rho_neg + rho_pos);
                    if model.pressure(z, t, rho_mid) > 0.0 {
                        rho_pos = rho_mid;
                    } else {
                        rho_neg = rho_mid;
                    }
                }
                if let Some(rho_cand) = bracketed_root(rho_pos, rho_max) {
                    if is_mechanically_stable(rho_cand) {
                        rho_liq = rho_cand;
                    }
                }
            }
        }

        let rho_cand = if rho_gas > 0.0 && rho_liq > 0.0 {
            let gibbs_gas = model.gibbsmolar_nocache(z, t, rho_gas);
            let gibbs_liq = model.gibbsmolar_nocache(z, t, rho_liq);
            if gibbs_liq <= gibbs_gas {
                rho_liq
            } else {
                rho_gas
            }
        } else if rho_gas > 0.0 {
            rho_gas
        } else {
            rho_liq
        };
        if rho_cand > 0.0 {
            let p_eos = model.pressure(z, t, rho_cand);
            if p_eos > 0.0 {
                if let Ok(verdict) = check_stability_michelsen(model, z, t, p_eos) {
                    if verdict.stable {
                        solved = Some(direct_state(model, z, t, rho_cand));
                    }
                }
            }
        }
    }

    let state = match solved {
        Some(s) => s,
        None => p_sweep_at_t(model, z, t, var, value)?,
    };

    // Verify (H/S/U only; D+T is a direct evaluation).
    if var != SweepVar::Dmolar {
        let resid = state.keyed(model, z, var) - value;
        let scale = value.abs() + 1.0;
        if !resid.is_finite() || resid.abs() > 1e-6 * scale {
            return Err(Error::Value(format!(
                "DHSU_T_flash for mixture did not converge to the specification: residual {resid:e} (target {}={value:e}) at T={t} K -- the (T,p) flash is misclassifying the phase for this mixture",
                var.short()
            )));
        }
    }
    Ok(state)
}

/// Upstream `HSU_P_flash` mixture branch: (Hmolar|Smolar|Umolar, P).
pub fn hsu_p_flash_mixtures(
    model: &MixtureModel,
    z: &[f64],
    p: f64,
    var: SweepVar,
    value: f64,
) -> Result<PtFlashResult> {
    let mut tmin = weighted_tmin(model, z);
    let mut tmax = weighted_tmax(model, z);

    // PQ-based bracket narrowing / exact-saturation shortcut.
    let pq = (|| -> Result<(f64, f64, f64, f64)> {
        let bubble = model.pq_flash(p, 0.0, z)?;
        let t_bubble = bubble.t;
        let val_bubble = two_phase_keyed(model, &bubble, var);
        let dew = model.pq_flash(p, 1.0, z)?;
        let t_dew = dew.t;
        let val_dew = two_phase_keyed(model, &dew, var);
        if t_bubble.is_finite()
            && t_dew.is_finite()
            && val_bubble.is_finite()
            && val_dew.is_finite()
        {
            Ok((t_bubble, val_bubble, t_dew, val_dew))
        } else {
            Err(Error::Value("non-finite PQ".into()))
        }
    })();
    if let Ok((t_bubble, val_bubble, t_dew, val_dew)) = pq {
        let tol_sat = 1e-6 * (value.abs() + 1.0);
        if (value - val_bubble).abs() < tol_sat {
            return model.pq_flash(p, 0.0, z).map(two_phase_result);
        }
        if (value - val_dew).abs() < tol_sat {
            return model.pq_flash(p, 1.0, z).map(two_phase_result);
        }
        let (val_lo, val_hi, t_at_lo, t_at_hi) = if val_bubble < val_dew {
            (val_bubble, val_dew, t_bubble, t_dew)
        } else {
            (val_dew, val_bubble, t_dew, t_bubble)
        };
        if value < val_lo {
            tmax = t_at_lo;
        } else if value > val_hi {
            tmin = t_at_hi;
        }
    }

    let mut resid = |t: f64| -> Result<f64> {
        let state = pt_flash_mixtures(model, z, t, p)?;
        Ok(state.keyed(model, z, var) - value)
    };

    let mut resid_lo = f64::NAN;
    let mut resid_hi = f64::NAN;
    let mut lo_ok = false;
    let mut hi_ok = false;
    if let Ok(r) = resid(tmin) {
        resid_lo = r;
        lo_ok = r.is_finite();
    }
    if let Ok(r) = resid(tmax) {
        resid_hi = r;
        hi_ok = r.is_finite();
    }

    // Binary-search validity recovery from a failing endpoint.
    if !lo_ok && hi_ok {
        let mut a = tmin;
        let mut b = tmax;
        for _ in 0..50 {
            if (b - a) <= 0.01 {
                break;
            }
            let mid = 0.5 * (a + b);
            match resid(mid) {
                Ok(r) if r.is_finite() => {
                    resid_lo = r;
                    lo_ok = true;
                    tmin = mid;
                    b = mid;
                }
                _ => {
                    a = mid;
                }
            }
        }
    } else if !hi_ok && lo_ok {
        let mut a = tmin;
        let mut b = tmax;
        for _ in 0..50 {
            if (b - a) <= 0.01 {
                break;
            }
            let mid = 0.5 * (a + b);
            match resid(mid) {
                Ok(r) if r.is_finite() => {
                    resid_hi = r;
                    hi_ok = true;
                    tmax = mid;
                    a = mid;
                }
                _ => {
                    b = mid;
                }
            }
        }
    }

    let final_state = if lo_ok && hi_ok && resid_lo * resid_hi < 0.0 {
        let t_sol = toms748_standin(&mut resid, tmin, tmax, resid_lo, resid_hi).map_err(|e| {
            Error::Value(format!(
                "HSU_P_flash for mixture failed with Tmin={tmin}, Tmax={tmax}, p={p}: {e}"
            ))
        })?;
        pt_flash_mixtures(model, z, t_sol, p)?
    } else if lo_ok || hi_ok {
        let t_best = if !hi_ok || (lo_ok && resid_lo.abs() <= resid_hi.abs()) {
            tmin
        } else {
            tmax
        };
        pt_flash_mixtures(model, z, t_best, p).map_err(|e| {
            Error::Value(format!(
                "HSU_P_flash for mixture: endpoints do not bracket (resid_lo={resid_lo:e}, resid_hi={resid_hi:e}) and fallback at T={t_best} failed: {e}"
            ))
        })?
    } else {
        return Err(Error::Value(format!(
            "HSU_P_flash for mixture: neither endpoint evaluable (Tmin={tmin}, Tmax={tmax}, p={p})"
        )));
    };

    let resid_final = final_state.keyed(model, z, var) - value;
    let resid_scale = value.abs() + 1.0;
    if !resid_final.is_finite() || resid_final.abs() > 1e-6 * resid_scale {
        return Err(Error::Value(format!(
            "HSU_P_flash for mixture did not converge to the specification: residual {resid_final:e} (target {value:e}) at p={p} Pa -- the (T,p) flash is misclassifying the phase for this mixture"
        )));
    }
    Ok(final_state)
}

fn two_phase_keyed(
    model: &MixtureModel,
    s: &crate::mixture_vle::MixtureTwoPhase,
    var: SweepVar,
) -> f64 {
    match var {
        SweepVar::Dmolar => s.rhomolar,
        SweepVar::Hmolar => {
            s.q * model.hmolar(&s.y_vap, s.t, s.rhomolar_vap)
                + (1.0 - s.q) * model.hmolar(&s.x_liq, s.t, s.rhomolar_liq)
        }
        SweepVar::Smolar => {
            s.q * model.smolar(&s.y_vap, s.t, s.rhomolar_vap)
                + (1.0 - s.q) * model.smolar(&s.x_liq, s.t, s.rhomolar_liq)
        }
        SweepVar::Umolar => {
            s.q * model.umolar(&s.y_vap, s.t, s.rhomolar_vap)
                + (1.0 - s.q) * model.umolar(&s.x_liq, s.t, s.rhomolar_liq)
        }
    }
}

fn two_phase_result(s: crate::mixture_vle::MixtureTwoPhase) -> PtFlashResult {
    PtFlashResult::TwoPhase {
        t: s.t,
        p: s.p,
        q: s.q,
        rhomolar: s.rhomolar,
        x: s.x_liq,
        y: s.y_vap,
        rhomolar_liq: s.rhomolar_liq,
        rhomolar_vap: s.rhomolar_vap,
    }
}

/// Upstream `HSU_D_flash` mixture branch: (Dmolar, Hmolar|Smolar|Umolar).
pub fn hsu_d_flash_mixtures(
    model: &MixtureModel,
    z: &[f64],
    rho_target: f64,
    var: SweepVar,
    value: f64,
) -> Result<PtFlashResult> {
    let mut solved: Option<PtFlashResult> = None;

    // Fast path: T-sweep at fixed rho via nocache evaluations.
    {
        let tmin = weighted_tmin(model, z);
        let tmax = weighted_tmax(model, z);
        let mut nocache_resid =
            |t: f64| -> Result<f64> { Ok(nocache_value(model, z, t, rho_target, var) - value) };
        let is_mechanically_stable = |t: f64| -> bool {
            let eps = (rho_target * 1e-6).max(1e-10);
            let p_lo = model.pressure(z, t, rho_target - eps);
            let p_hi = model.pressure(z, t, rho_target + eps);
            (p_hi - p_lo) > 0.0
        };
        let root = (|| -> Option<f64> {
            let fa = nocache_resid(tmin).ok()?;
            let fb = nocache_resid(tmax).ok()?;
            if fa * fb > 0.0 {
                return None;
            }
            toms748_standin(&mut nocache_resid, tmin, tmax, fa, fb).ok()
        })();
        if let Some(t_cand) = root {
            let p_eos = model.pressure(z, t_cand, rho_target);
            if p_eos > 0.0 && is_mechanically_stable(t_cand) {
                if let Ok(verdict) = check_stability_michelsen(model, z, t_cand, p_eos) {
                    if verdict.stable {
                        solved = Some(direct_state(model, z, t_cand, rho_target));
                    }
                }
            }
        }
    }

    // Slow path: nested 1D (outer T log-scan, inner P-sweep to hit rho).
    if solved.is_none() {
        let tmin = weighted_tmin(model, z);
        let tmax = weighted_tmax(model, z);
        let mut pmin_bound = weighted_p_triple(model, z);
        let pmax_bound = weighted_pmax(model, z);
        if pmin_bound < 100.0 {
            pmin_bound = 100.0;
        }

        // Inner: at T, find P with rho(T,P) = rho_target; return the caloric.
        let solve_for_caloric_at_t = |t: f64| -> Result<(f64, PtFlashResult)> {
            let mut rho_resid = |p: f64| -> Result<f64> {
                let state = pt_flash_mixtures(model, z, t, p)?;
                Ok(state.keyed(model, z, SweepVar::Dmolar) - rho_target)
            };
            let log_pmin = pmin_bound.log10();
            let log_pmax = pmax_bound.log10();
            let dlog_p = 0.5;
            let nsteps = ((log_pmax - log_pmin) / dlog_p) as i32 + 1;
            let mut p_lo = -1.0;
            let mut p_hi = -1.0;
            let mut f_lo = f64::NAN;
            let mut f_hi = f64::NAN;
            let mut f_prev = 0.0;
            let mut p_prev = 0.0;
            let mut have_prev = false;
            for i in 0..=nsteps {
                let mut log_p = log_pmin + f64::from(i) * dlog_p;
                if log_p > log_pmax {
                    log_p = log_pmax;
                }
                let p = 10.0_f64.powf(log_p);
                let f = match rho_resid(p) {
                    Ok(f) => f,
                    Err(_) => {
                        have_prev = false;
                        continue;
                    }
                };
                if have_prev && f_prev * f < 0.0 {
                    p_lo = p_prev;
                    p_hi = p;
                    f_lo = f_prev;
                    f_hi = f;
                    break;
                }
                p_prev = p;
                f_prev = f;
                have_prev = true;
            }
            if p_lo < 0.0 || p_hi < 0.0 {
                return Err(Error::Value("no inner P bracket".into()));
            }
            let p_sol = toms748_standin(&mut rho_resid, p_lo, p_hi, f_lo, f_hi)?;
            let state = pt_flash_mixtures(model, z, t, p_sol)?;
            let x = state.keyed(model, z, var);
            Ok((x, state))
        };

        // Outer: T log-scan (dlogT = 0.1) for a sign change in X - value.
        let log_tmin = tmin.log10();
        let log_tmax = tmax.log10();
        let dlog_t = 0.1;
        let nsteps_t = ((log_tmax - log_tmin) / dlog_t) as i32 + 1;
        let mut t_lo = -1.0;
        let mut t_hi = -1.0;
        let mut f_lo = f64::NAN;
        let mut f_hi = f64::NAN;
        let mut f_prev_t = 0.0;
        let mut t_prev = 0.0;
        let mut have_prev_t = false;
        for i in 0..=nsteps_t {
            let mut log_t = log_tmin + f64::from(i) * dlog_t;
            if log_t > log_tmax {
                log_t = log_tmax;
            }
            let t = 10.0_f64.powf(log_t);
            let x = match solve_for_caloric_at_t(t) {
                Ok((x, _)) => x,
                Err(_) => {
                    have_prev_t = false;
                    continue;
                }
            };
            let f = x - value;
            if have_prev_t && f_prev_t * f < 0.0 {
                t_lo = t_prev;
                t_hi = t;
                f_lo = f_prev_t;
                f_hi = f;
                break;
            }
            t_prev = t;
            f_prev_t = f;
            have_prev_t = true;
        }

        if t_lo > 0.0 && t_hi > 0.0 {
            let mut outer_resid = |t: f64| -> Result<f64> {
                let (x, _) = solve_for_caloric_at_t(t)
                    .map_err(|_| Error::Value("inner P-sweep failed during T-sweep".into()))?;
                Ok(x - value)
            };
            let t_sol = toms748_standin(&mut outer_resid, t_lo, t_hi, f_lo, f_hi)?;
            let (_, state) = solve_for_caloric_at_t(t_sol)?;
            solved = Some(state);
        } else {
            return Err(Error::Value(format!(
                "HSU_D_flash for mixture: no T bracket found in [{tmin}, {tmax}] for target {}={value} at rho={rho_target}",
                var.short()
            )));
        }
    }

    let state = solved.expect("state set on both paths");
    // Verify BOTH the caloric and the density.
    let resid_cal = state.keyed(model, z, var) - value;
    let resid_rho = state.keyed(model, z, SweepVar::Dmolar) - rho_target;
    let scale_cal = value.abs() + 1.0;
    let scale_rho = rho_target.abs() + 1.0;
    if !resid_cal.is_finite()
        || resid_cal.abs() > 1e-6 * scale_cal
        || !resid_rho.is_finite()
        || resid_rho.abs() > 1e-6 * scale_rho
    {
        return Err(Error::Value(format!(
            "HSU_D_flash for mixture did not converge to the specification: caloric residual {resid_cal:e} (target {}={value:e}), density residual {resid_rho:e} (target {rho_target:e}) -- the (T,p) flash is misclassifying the phase",
            var.short()
        )));
    }
    Ok(state)
}
