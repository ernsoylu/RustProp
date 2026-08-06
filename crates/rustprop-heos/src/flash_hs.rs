//! (H,S) flash (PLAN.md 4.6, final pair) — port of upstream
//! `FlashRoutines::HS_flash`'s superancillary "happy path" and its
//! machinery in `src/Backends/Helmholtz/FlashRoutines.cpp` @ v8.0.0:
//!
//! - the caloric superancillaries (`SuperAncillary::add_variable`): Chebyshev
//!   expansions of h and s along both saturation branches, built at runtime
//!   from the EOS at the rho-expansions' own Lobatto nodes (upstream builds
//!   lazily under a mutex; we build lazily via `OnceCell` — same trigger
//!   points, single-threaded);
//! - the fast two-phase screen `hs_two_phase_likely` (no EOS: scan the
//!   Qh==Qs residual on the caloric curves);
//! - the single-phase cascade `hs_cascade` with its three dome-free legs
//!   (saturation anchor, supercritical isentrope, ideal-gas departure via
//!   lambda-continuation) sharing the `(T, ln rho)` homotopy corrector, and
//!   the stability-based acceptance `hs_accept` + `hs_inside_dome` veto;
//! - the EOS-exact two-phase solve `HS_flash_twophase` (Brent on Qh - Qs
//!   over the saturation temperature range);
//! - the legacy TS-scan "sad path" behind the superancillary path (scan on
//!   T with the (Smolar,T) flash inside, first-adjacent-sign-change
//!   bracketing at 30 bits, whole-range Brent fallback).
//!
//! Scoped deviations, logged in PLAN.md:
//! - cascade leg 4 (melting-line caloric anchor for the sub-triple
//!   compressed-liquid corner) is deferred with the melting line itself;
//!   legs 1-3 cover the EOS range above the triple temperature (leg 1
//!   anchors cold liquids at the triple-point saturated liquid);
//! - the reference-state offset "stamp" shift (#2773) is ported as the full
//!   formula, but reference-state mutation does not exist in this port, so
//!   the shift always evaluates to exactly zero;
//! - the corrector's Jacobian entries come from the same closed forms the
//!   departure leg's lambda-model uses at lambda = 1 (mathematically
//!   identical to upstream's `first_partial_deriv` values; both correctors
//!   converge the identical residual to norm < 1e-11, so ULP-level Jacobian
//!   differences cannot move the accepted root beyond that norm).

use crate::alpha::HelmholtzEos;
use crate::chebappr::{ChebApprox1d, bisect_bits, lu_matrices, mat_vec, nodes_realworld};
use crate::flash_pt::PtFlash;
use crate::flash_px::HeosState;
use crate::solvers::brent;
use crate::superancillary::{OwnedInterval, eval_sat};
use rustprop_core::fluid::{Alpha0Term, FluidData, SuperAncillaryData};
use rustprop_core::{Error, Result};

// ---------------------------------------------------------------------------
// Caloric superancillaries (upstream SuperAncillary::add_variable / #2773)
// ---------------------------------------------------------------------------

/// Runtime-built caloric saturation expansions for one fluid: h and s along
/// both branches, on the same intervals as the rho expansions, plus the
/// alpha0-offset stamp recorded at build time.
pub(crate) struct CaloricSa {
    h_l: ChebApprox1d,
    h_v: ChebApprox1d,
    s_l: ChebApprox1d,
    s_v: ChebApprox1d,
    /// (a1, a2) totals of `IdealGasHelmholtzEnthalpyEntropyOffset` at build
    /// time (upstream `get_caloric_alpha0_stamp`).
    stamp: (f64, f64),
}

/// Total (a1, a2) of the document's enthalpy/entropy offset terms — upstream
/// `FlashRoutines::alpha0_offset_total` with the parse-time Core slot only
/// (runtime reference-state offsets don't exist in this port; the alpha0
/// prefactor is 1 for every bundled fluid).
fn alpha0_offset_total(fluid: &FluidData) -> (f64, f64) {
    let (mut a1, mut a2) = (0.0, 0.0);
    for term in fluid.eos.alpha0 {
        if let Alpha0Term::EnthalpyEntropyOffset { a1: x1, a2: x2, .. } = term {
            a1 += x1;
            a2 += x2;
        }
    }
    (a1, a2)
}

/// Upstream `SuperAncillary::add_variable_locked` for H and S (U joins with
/// the HSU_D flashes): rebuild each rho interval's node values through the
/// EOS and refit with the L matrix. The degree-12 L/U matrices are
/// upstream's hard-coded choice; the rhoL/rhoV interval grids are asserted
/// identical, which upstream assumes implicitly.
fn build_calorics(eos: &HelmholtzEos, data: &SuperAncillaryData, fluid: &FluidData) -> CaloricSa {
    let ndeg = 12usize;
    let (lmat, umat) = lu_matrices(ndeg);
    assert_eq!(
        data.rho_l.len(),
        data.rho_v.len(),
        "L&V are not the same size"
    );

    let mut h_l = Vec::new();
    let mut h_v = Vec::new();
    let mut s_l = Vec::new();
    let mut s_v = Vec::new();
    for (el, ev) in data.rho_l.iter().zip(data.rho_v) {
        assert!(
            el.xmin == ev.xmin && el.xmax == ev.xmax && el.coef.len() == ndeg + 1,
            "caloric build requires matching degree-12 rhoL/rhoV intervals"
        );
        let t_nodes = nodes_realworld(ndeg, el.xmin, el.xmax);
        let rho_l_nodes = mat_vec(&umat, el.coef);
        let rho_v_nodes = mat_vec(&umat, ev.coef);
        let build = |vals: Vec<f64>, xmin: f64, xmax: f64| -> OwnedInterval {
            let coef = mat_vec(&lmat, &vals);
            assert!(
                coef.iter().all(|c| c.is_finite()),
                "At least one caloric coefficient is non-finite"
            );
            OwnedInterval { xmin, xmax, coef }
        };
        let hvals = |rhos: &[f64]| -> Vec<f64> {
            t_nodes
                .iter()
                .zip(rhos)
                .map(|(&t, &rho)| eos.hmolar(t, rho))
                .collect()
        };
        let svals = |rhos: &[f64]| -> Vec<f64> {
            t_nodes
                .iter()
                .zip(rhos)
                .map(|(&t, &rho)| eos.smolar(t, rho))
                .collect()
        };
        h_l.push(build(hvals(&rho_l_nodes), el.xmin, el.xmax));
        h_v.push(build(hvals(&rho_v_nodes), ev.xmin, ev.xmax));
        s_l.push(build(svals(&rho_l_nodes), el.xmin, el.xmax));
        s_v.push(build(svals(&rho_v_nodes), ev.xmin, ev.xmax));
    }
    CaloricSa {
        h_l: ChebApprox1d::new(h_l),
        h_v: ChebApprox1d::new(h_v),
        s_l: ChebApprox1d::new(s_l),
        s_v: ChebApprox1d::new(s_v),
        stamp: alpha0_offset_total(fluid),
    }
}

impl CaloricSa {
    /// `SuperAncillary::eval_sat` for the caloric variables.
    fn eval_sat(&self, t: f64, k: char, q: u8) -> f64 {
        match (k, q) {
            ('H', 0) => self.h_l.eval(t),
            ('H', 1) => self.h_v.eval(t),
            ('S', 0) => self.s_l.eval(t),
            ('S', 1) => self.s_v.eval(t),
            _ => panic!("bad caloric key {k:?}/Q={q}"),
        }
    }

    /// `SuperAncillary::get_all_intersections`: every saturation temperature
    /// where either branch of variable `k` equals `val` (liquid solutions
    /// first, then vapor, as upstream concatenates).
    fn get_all_intersections(
        &self,
        k: char,
        val: f64,
        bits: u32,
        max_iter: usize,
        boundsftol: f64,
    ) -> Vec<f64> {
        let (l, v) = match k {
            'H' => (&self.h_l, &self.h_v),
            'S' => (&self.s_l, &self.s_v),
            _ => panic!("bad caloric key {k:?}"),
        };
        let mut solns = l.get_x_for_y(val, bits, max_iter, boundsftol);
        solns.extend(v.get_x_for_y(val, bits, max_iter, boundsftol));
        solns
    }
}

/// Target entropy shifted into the caloric cache's stamped frame (#2773).
/// Reference-state mutation is not ported, so `stamp == current` and the
/// shift is exactly zero — kept as the full formula for fidelity.
fn hs_s_to_cache(eos: &HelmholtzEos, cal: &CaloricSa, current: (f64, f64), s_t: f64) -> f64 {
    s_t - eos.gas_constant * (cal.stamp.0 - current.0)
}

/// Target enthalpy shifted into the caloric cache's stamped frame (#2773).
fn hs_h_to_cache(eos: &HelmholtzEos, cal: &CaloricSa, current: (f64, f64), h_t: f64) -> f64 {
    h_t + eos.gas_constant * eos.t_reducing * (cal.stamp.1 - current.1)
}

// ---------------------------------------------------------------------------
// The lambda-scaled model (upstream hs_leg_departure's lprops); at
// lambda = 1 these are the physical h, s, dp/drho|T and the (T, rho)
// partials the corrector's Newton step needs.
// ---------------------------------------------------------------------------

struct Lp {
    h: f64,
    s: f64,
    prho: f64,
    dh_dt: f64,
    dh_drho: f64,
    ds_dt: f64,
    ds_drho: f64,
}

/// Properties of `alpha = alpha0 + lam*alphar` at (T, rho). A non-finite
/// value aborts the calling leg (upstream throws `ValueError` here).
fn lprops(eos: &HelmholtzEos, t: f64, rho: f64, lam: f64) -> Result<Lp> {
    let rg = eos.gas_constant;
    let tr = eos.t_reducing;
    let rhor = eos.rhomolar_reducing;
    let tau = tr / t;
    let delta = rho / rhor;
    let a0d = eos.alpha0_all(tau, delta);
    let ard_all = eos.alphar_all(tau, delta);
    let (a0, a0t, a0tt) = (a0d.d00, a0d.d01, a0d.d02);
    let (ar, art, ard) = (ard_all.d00, ard_all.d01, ard_all.d10);
    let (artt, ardd, artd) = (ard_all.d02, ard_all.d20, ard_all.d11);
    let at = a0t + lam * art;
    let att = a0tt + lam * artt;
    let hh = 1.0 + tau * at + delta * lam * ard;
    let l = Lp {
        h: rg * t * hh,
        s: rg * (tau * at - (a0 + lam * ar)),
        prho: rg * t * (1.0 + 2.0 * delta * lam * ard + delta * delta * lam * ardd),
        dh_dt: rg * (hh - tau * (at + tau * att + delta * lam * artd)),
        dh_drho: rg * t * (lam * ard + tau * lam * artd + delta * lam * ardd) / rhor,
        ds_dt: -rg * tau * tau * att / t,
        ds_drho: rg * (tau * lam * artd - 1.0 / delta - lam * ard) / rhor,
    };
    if !l.h.is_finite()
        || !l.s.is_finite()
        || !l.prho.is_finite()
        || !l.dh_dt.is_finite()
        || !l.dh_drho.is_finite()
        || !l.ds_dt.is_finite()
        || !l.ds_drho.is_finite()
    {
        return Err(Error::Value(
            "hs lambda model: non-finite property/derivative".into(),
        ));
    }
    Ok(l)
}

// ---------------------------------------------------------------------------
// Homotopy corrector (upstream hs_corrector)
// ---------------------------------------------------------------------------

/// Shared (T, w=ln rho) homotopy corrector: homotope (h,s) linearly from the
/// anchor's values to the target with adaptive subdivision; the
/// dp/drho|_T > 0 guard keeps the corrector on the mechanically stable
/// branch. Returns the converged (T, rho), or None when every subdivision
/// schedule fails (or an evaluation goes non-finite, upstream's throw).
fn hs_corrector(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    t0: f64,
    rho0: f64,
    h_t: f64,
    s_t: f64,
    tlo_override: f64,
) -> Option<(f64, f64)> {
    let rgas = eos.gas_constant;
    let tsc = fluid.states.critical.t;
    let (hscale, sscale) = (rgas * tsc, rgas);
    let tlo = if tlo_override > 0.0 {
        tlo_override
    } else {
        fluid.eos.sat_min_liquid.t * (1.0 - 2e-2)
    };
    let thi = 1.5 * fluid.eos.t_max;
    let l0 = lprops(eos, t0, rho0, 1.0).ok()?;
    let (h0, s0) = (l0.h, l0.s);
    let mut n = 1;
    while n <= 128 {
        let mut t = t0;
        let mut w = rho0.ln();
        let mut failed = false;
        let mut k = 1;
        while k <= n && !failed {
            let lam = f64::from(k) / f64::from(n);
            let ht = h0 + lam * (h_t - h0);
            let st = s0 + lam * (s_t - s0);
            let mut conv = false;
            let mut best_norm = 1e300;
            let mut best_t = t;
            let mut best_w = w;
            for _iter in 0..40 {
                let rho = w.exp();
                let Ok(l) = lprops(eos, t, rho, 1.0) else {
                    return None;
                };
                let rh = l.h - ht;
                let rs = l.s - st;
                let norm = rh.abs() / hscale + rs.abs() / sscale;
                if norm < best_norm {
                    best_norm = norm;
                    best_t = t;
                    best_w = w;
                }
                if norm < 1e-11 {
                    conv = true;
                    break;
                }
                let a11 = l.dh_dt;
                let a12 = rho * l.dh_drho;
                let a21 = l.ds_dt;
                let a22 = rho * l.ds_drho;
                let det = a11 * a22 - a12 * a21;
                if !det.is_finite() || det.abs() < 1e-300 {
                    break;
                }
                let dt = -(a22 * rh - a12 * rs) / det;
                let dw = -(-a21 * rh + a11 * rs) / det;
                let mut f = 1.0;
                while (t + f * dt < tlo || t + f * dt > thi || (f * dw).abs() > 2.0) && f > 1e-6 {
                    f *= 0.5;
                }
                for _g in 0..8 {
                    if f <= 1e-3 {
                        break;
                    }
                    let Ok(trial) = lprops(eos, t + f * dt, (w + f * dw).exp(), 1.0) else {
                        return None;
                    };
                    if trial.prho > 0.0 {
                        break;
                    }
                    f *= 0.5;
                }
                t += f * dt;
                w += f * dw;
            }
            if !conv && best_norm < 1e-8 {
                t = best_t;
                w = best_w;
                conv = true;
            }
            if !conv {
                failed = true;
            }
            k += 1;
        }
        if !failed {
            return Some((t, w.exp()));
        }
        n *= 2;
    }
    None
}

// ---------------------------------------------------------------------------
// Cascade legs
// ---------------------------------------------------------------------------

/// Leg 1: saturation anchor — anchor at the saturated state whose entropy
/// equals the target's (native caloric S inversion), disambiguated by
/// caloric-H closeness, with the dilute-gas / cold-liquid fallbacks.
#[allow(clippy::too_many_arguments)]
fn hs_leg_saturation(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    data: &SuperAncillaryData,
    cal: &CaloricSa,
    h_t: f64,
    s_t: f64,
) -> Option<(f64, f64)> {
    let tc = data.t_crit_num;
    let tmin = data.p[0].xmin;
    let tcap = tc - (1e-4f64).max(1e-6 * tc);
    let (a, b) = (tmin + 1e-3, tcap);
    let current = alpha0_offset_total(fluid);
    let s_cache = hs_s_to_cache(eos, cal, current, s_t);
    let h_cache = hs_h_to_cache(eos, cal, current, h_t);

    let mut t0 = 0.0;
    let mut rho0 = 0.0;
    let mut best_hgap = 1e300;
    let mut found = false;
    let mut near_critical = false;
    for tcand in cal.get_all_intersections('S', s_cache, 48, 100, 1e-12) {
        if tcand >= b {
            near_critical = true;
        }
        if tcand <= a || tcand >= b {
            continue;
        }
        let qroot = if (cal.eval_sat(tcand, 'S', 0) - s_cache).abs()
            <= (cal.eval_sat(tcand, 'S', 1) - s_cache).abs()
        {
            0
        } else {
            1
        };
        let rho = eval_sat(data, tcand, 'D', qroot);
        let hgap = (cal.eval_sat(tcand, 'H', qroot) - h_cache).abs();
        if hgap.is_finite() && hgap < best_hgap {
            best_hgap = hgap;
            t0 = tcand;
            rho0 = rho;
            found = true;
        }
    }
    if !found {
        if near_critical {
            return None; // hand to the supercritical isentrope leg
        }
        let s_crit = eos.smolar(tcap, eval_sat(data, tcap, 'D', 0));
        if s_t > s_crit {
            // dilute gas: T from the monotone 1D enthalpy inversion
            rho0 = data.rho_crit_num * 1e-4;
            let hres = |t: f64| eos.hmolar(t, rho0) - h_t;
            let ta = tmin + 1e-3;
            let tb = 1.5 * fluid.eos.t_max;
            let (fa, fb) = (hres(ta), hres(tb));
            if fa * fb <= 0.0 {
                t0 = bisect_bits(hres, ta, tb, fa, fb, 30, 60);
            } else {
                t0 = if fa.abs() <= fb.abs() { ta } else { tb };
            }
        } else {
            // cold compressed liquid: triple-point saturated liquid
            t0 = a;
            rho0 = eval_sat(data, a, 'D', 0);
        }
    }
    hs_corrector(eos, fluid, t0, rho0, h_t, s_t, -1.0)
}

/// Leg 2: supercritical isentrope — anchor on the dome-free T=Tmax isotherm
/// at the density whose entropy equals s_t.
fn hs_leg_isentrope(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    h_t: f64,
    s_t: f64,
) -> Option<(f64, f64)> {
    let ta = fluid.eos.t_max;
    let rhoc = fluid.states.critical.rhomolar;
    let sres = |rho: f64| eos.smolar(ta, rho) - s_t;
    let ra = 1e-6 * rhoc;
    let ga = sres(ra);
    let mut rb = 6.0 * rhoc;
    let mut gb = sres(rb);
    for _e in 0..12 {
        if ga * gb <= 0.0 {
            break;
        }
        let rb_next = rb * 2.0;
        let gb_next = sres(rb_next);
        if !gb_next.is_finite() {
            break;
        }
        rb = rb_next;
        gb = gb_next;
    }
    let rho0 = if ga * gb > 0.0 {
        if ga.abs() <= gb.abs() { ra } else { rb }
    } else {
        bisect_bits(sres, ra, rb, ga, gb, 40, 80)
    };
    hs_corrector(eos, fluid, ta, rho0, h_t, s_t, -1.0)
}

/// Leg 3: ideal-gas departure — scale residual Helmholtz by lambda; anchor
/// at lambda=0 (no dome at all) and continue lambda 0 -> 1.
fn hs_leg_departure(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    h_t: f64,
    s_t: f64,
) -> Option<(f64, f64)> {
    let rg = eos.gas_constant;
    let tr = eos.t_reducing;
    let rhor = eos.rhomolar_reducing;
    let tmin = fluid.eos.sat_min_liquid.t;
    let tmax = fluid.eos.t_max;

    // lambda=0 anchor: T from h(T) at the reducing density, then rho from
    // s(T0, rho) — two decoupled 1-D solves.
    let hig = |t: f64| lprops(eos, t, rhor, 0.0).map_or(f64::NAN, |l| l.h - h_t);
    let (ta, tb) = (0.3 * tmin, 3.0 * tmax);
    let (fa, fb) = (hig(ta), hig(tb));
    if !fa.is_finite() || !fb.is_finite() || fa * fb > 0.0 {
        return None;
    }
    let t0 = bisect_bits(hig, ta, tb, fa, fb, 40, 80);
    let sig = |rho: f64| lprops(eos, t0, rho, 0.0).map_or(f64::NAN, |l| l.s - s_t);
    let (ra, rb) = (1e-8 * rhor, 50.0 * rhor);
    let (ga, gb) = (sig(ra), sig(rb));
    if !ga.is_finite() || !gb.is_finite() {
        return None;
    }
    let rho0 = if ga * gb > 0.0 {
        if ga.abs() <= gb.abs() { ra } else { rb }
    } else {
        bisect_bits(sig, ra, rb, ga, gb, 40, 80)
    };

    // Continuation lam: 0 -> 1, N-doubling schedule.
    let (hscale, sscale) = (rg * tr, rg);
    let mut n = 1;
    while n <= 256 {
        let mut t = t0;
        let mut w = rho0.ln();
        let mut failed = false;
        let mut k = 1;
        while k <= n && !failed {
            let lam = f64::from(k) / f64::from(n);
            let mut conv = false;
            let mut best_norm = 1e300;
            let mut best_t = t;
            let mut best_w = w;
            for _iter in 0..40 {
                let rho = w.exp();
                let Ok(l) = lprops(eos, t, rho, lam) else {
                    return None;
                };
                let rh = l.h - h_t;
                let rs = l.s - s_t;
                let norm = rh.abs() / hscale + rs.abs() / sscale;
                if norm < best_norm {
                    best_norm = norm;
                    best_t = t;
                    best_w = w;
                }
                if norm < 1e-11 {
                    conv = true;
                    break;
                }
                let a11 = l.dh_dt;
                let a12 = rho * l.dh_drho;
                let a21 = l.ds_dt;
                let a22 = rho * l.ds_drho;
                let det = a11 * a22 - a12 * a21;
                if !det.is_finite() || det.abs() < 1e-300 {
                    break;
                }
                let dt = -(a22 * rh - a12 * rs) / det;
                let dw = -(-a21 * rh + a11 * rs) / det;
                let mut f = 1.0;
                while (t + f * dt <= 0.0 || t + f * dt > 3.0 * tmax || (f * dw).abs() > 2.0)
                    && f > 1e-6
                {
                    f *= 0.5;
                }
                for _g in 0..8 {
                    if f <= 1e-3 {
                        break;
                    }
                    let Ok(trial) = lprops(eos, t + f * dt, (w + f * dw).exp(), lam) else {
                        return None;
                    };
                    if trial.prho > 0.0 {
                        break;
                    }
                    f *= 0.5;
                }
                t += f * dt;
                w += f * dw;
            }
            if !conv && best_norm < 1e-8 {
                t = best_t;
                w = best_w;
                conv = true;
            }
            if !conv {
                failed = true;
            }
            k += 1;
        }
        if !failed {
            return Some((t, w.exp()));
        }
        n *= 2;
    }
    None
}

// ---------------------------------------------------------------------------
// Acceptance and dispatch
// ---------------------------------------------------------------------------

/// dp/drho|_T (upstream `first_partial_deriv(iP, iDmolar, iT)`).
fn dpdrho_t(eos: &HelmholtzEos, t: f64, rho: f64) -> f64 {
    let tau = eos.t_reducing / t;
    let delta = rho / eos.rhomolar_reducing;
    let d = eos.alphar_all(tau, delta);
    eos.gas_constant * t * (1.0 + 2.0 * delta * d.d10 + delta * delta * d.d20)
}

/// Accept only a faithful (h,s) reproduction that is in-range and fully
/// intrinsically stable (dp/drho|_T > 0 AND cv > 0) — over that region the
/// (h,s)->(T,rho) map is injective.
fn hs_accept(eos: &HelmholtzEos, fluid: &FluidData, t: f64, rho: f64, h_t: f64, s_t: f64) -> bool {
    if !t.is_finite() || !rho.is_finite() || rho <= 0.0 {
        return false;
    }
    let rg = eos.gas_constant;
    let tc = fluid.states.critical.t;
    let tmin_eff = fluid.eos.sat_min_liquid.t;
    if t < tmin_eff * (1.0 - 1e-6) || t > fluid.eos.t_max * (1.0 + 1e-6) {
        return false;
    }
    if (eos.hmolar(t, rho) - h_t).abs() > 1e-6 * rg * tc
        || (eos.smolar(t, rho) - s_t).abs() > 1e-6 * rg
    {
        return false;
    }
    let cv = eos.cvmolar(t, rho);
    if !cv.is_finite() || cv <= 0.0 {
        return false;
    }
    dpdrho_t(eos, t, rho) > 0.0
}

/// True if (T, rho) lies strictly inside the two-phase dome (a metastable
/// single-phase extension) — the cascade must reject such roots.
fn hs_inside_dome(data: &SuperAncillaryData, t: f64, rho: f64) -> bool {
    if t >= data.t_crit_num {
        return false;
    }
    rho > eval_sat(data, t, 'D', 1) && rho < eval_sat(data, t, 'D', 0)
}

/// The three-leg cascade (upstream `hs_cascade`; leg 4 deferred with the
/// melting line).
fn hs_cascade(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    data: &SuperAncillaryData,
    cal: &CaloricSa,
    h_t: f64,
    s_t: f64,
) -> Option<(f64, f64)> {
    let good =
        |t: f64, rho: f64| hs_accept(eos, fluid, t, rho, h_t, s_t) && !hs_inside_dome(data, t, rho);
    if let Some((t, rho)) = hs_leg_saturation(eos, fluid, data, cal, h_t, s_t)
        && good(t, rho)
    {
        return Some((t, rho));
    }
    if let Some((t, rho)) = hs_leg_isentrope(eos, fluid, h_t, s_t)
        && good(t, rho)
    {
        return Some((t, rho));
    }
    if let Some((t, rho)) = hs_leg_departure(eos, fluid, h_t, s_t)
        && good(t, rho)
    {
        return Some((t, rho));
    }
    None
}

/// Cheap superancillary two-phase detector (no EOS): scan the Qh==Qs
/// residual on the caloric curves for a sign change with quality strictly
/// inside (0,1). Conservative by design — both error directions are safe.
fn hs_two_phase_likely(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    data: &SuperAncillaryData,
    cal: &CaloricSa,
    h_t: f64,
    s_t: f64,
) -> bool {
    let current = alpha0_offset_total(fluid);
    let hc = hs_h_to_cache(eos, cal, current, h_t);
    let sc = hs_s_to_cache(eos, cal, current, s_t);
    let tmin = data.p[0].xmin;
    let tc = data.t_crit_num;
    let qh_minus_qs = |t: f64| -> f64 {
        let s_l = cal.eval_sat(t, 'S', 0);
        let s_v = cal.eval_sat(t, 'S', 1);
        let h_l = cal.eval_sat(t, 'H', 0);
        let h_v = cal.eval_sat(t, 'H', 1);
        (hc - h_l) / (h_v - h_l) - (sc - s_l) / (s_v - s_l)
    };
    let m = 40;
    let tlo = tmin + 1e-3;
    let thi = tc - (0.5f64).max(1e-3 * tc);
    let mut tprev = tlo;
    let mut fprev = qh_minus_qs(tlo);
    for i in 1..=m {
        let t = tlo + (thi - tlo) * f64::from(i) / f64::from(m);
        let f = qh_minus_qs(t);
        if fprev.is_finite() && f.is_finite() && fprev * f <= 0.0 {
            let tsol = bisect_bits(qh_minus_qs, tprev, t, fprev, f, 40, 60);
            let s_l = cal.eval_sat(tsol, 'S', 0);
            let s_v = cal.eval_sat(tsol, 'S', 1);
            let q = (sc - s_l) / (s_v - s_l);
            if q > 1e-6 && q < 1.0 - 1e-6 {
                return true; // strictly interior => two-phase
            }
        }
        tprev = t;
        fprev = f;
    }
    false
}

impl PtFlash {
    /// Lazily build (once) and return the caloric superancillaries
    /// (upstream `ensure_caloric_superancillaries`).
    pub(crate) fn hs_calorics(&self) -> &CaloricSa {
        self.hs_calorics_cell.get_or_init(|| {
            let data = self
                .fluid()
                .eos
                .superancillary
                .as_ref()
                .expect("HS flash currently requires a superancillary fluid");
            build_calorics(&self.eos, data, self.fluid())
        })
    }

    /// The EOS-exact two-phase solve (upstream `HS_flash_twophase`): Brent on
    /// the quality mismatch Qh - Qs over the saturation temperature range,
    /// then the state at the resolved (Q, T).
    fn hs_flash_twophase(&self, hmolar_spec: f64, smolar_spec: f64) -> Result<HeosState> {
        let fluid = self.fluid();
        let mut last_t = f64::NAN;
        let mut last_qs = f64::NAN;
        let resid = |t: f64| -> f64 {
            let Ok(sat) = self.sat().qt_flash(t, 0.0) else {
                return f64::NAN; // upstream: the inner QT update throws
            };
            let s_l = self.eos.smolar(t, sat.rho_l);
            let s_v = self.eos.smolar(t, sat.rho_v);
            let h_l = self.eos.hmolar(t, sat.rho_l);
            let h_v = self.eos.hmolar(t, sat.rho_v);
            let qs = (smolar_spec - s_l) / (s_v - s_l);
            let qh = (hmolar_spec - h_l) / (h_v - h_l);
            last_t = t;
            last_qs = qs;
            qh - qs
        };
        let tmax_sat = fluid.states.critical.t - 1e-13;
        let tmin_sat = fluid.eos.sat_min_liquid.t.max(fluid.eos.sat_min_vapor.t) - 1e-13;
        brent(resid, tmin_sat, tmax_sat - 0.01, f64::EPSILON, 1e-12, 20)?;
        // Run once more with the final vapor quality (upstream
        // `update(QT_INPUTS, resid.Qs, HEOS.T())` — the T of the last
        // residual evaluation).
        self.qt_state(last_t, last_qs)
    }

    /// (Hmolar, Smolar) flash — upstream `FlashRoutines::HS_flash`: the
    /// superancillary happy path (screen -> cascade -> two-phase fallback)
    /// with the legacy TS-scan sad path behind it.
    pub fn hmolar_smolar_state(&self, h_t: f64, s_t: f64) -> Result<HeosState> {
        let fluid = self.fluid();
        let data = fluid
            .eos
            .superancillary
            .as_ref()
            .expect("HS flash currently requires a superancillary fluid");
        let cal = self.hs_calorics();

        let reproduces = |st: &HeosState| -> bool {
            let hh = self.state_hmolar(st);
            let ss = self.state_smolar(st);
            hh.is_finite()
                && ss.is_finite()
                && (hh - h_t).abs() <= 1e-6 * h_t.abs() + 1e-3
                && (ss - s_t).abs() <= 1e-6 * s_t.abs() + 1e-5
        };

        let mut tried_2ph = false;
        // (0) Fast two-phase screen (no EOS).
        if hs_two_phase_likely(&self.eos, fluid, data, cal, h_t, s_t) {
            tried_2ph = true;
            match self.hs_flash_twophase(h_t, s_t) {
                Ok(st) => {
                    if reproduces(&st) {
                        return Ok(st);
                    }
                }
                // Upstream: the exception aborts the whole superancillary
                // block (skipping the cascade) and falls to the legacy path.
                Err(_) => return self.hs_legacy(h_t, s_t),
            }
        }
        // (1) Single-phase cascade.
        if let Some((t, rho)) = hs_cascade(&self.eos, fluid, data, cal, h_t, s_t) {
            let p = self.eos.pressure(t, rho);
            let st = HeosState::SinglePhase {
                t,
                p,
                rhomolar: rho,
                phase: self.recalculated_singlephase_phase(t, p, rho),
                // Upstream HS_flash leaves `_Q = 10000` on the cascade's
                // single-phase result (oracle-confirmed via PropsSI("Q")).
                q: 10000.0,
            };
            if reproduces(&st) {
                return Ok(st);
            }
        }
        // (2) Two-phase fallback (if the screen did not already try it).
        if !tried_2ph {
            match self.hs_flash_twophase(h_t, s_t) {
                Ok(st) => {
                    if reproduces(&st) {
                        return Ok(st);
                    }
                }
                Err(_) => return self.hs_legacy(h_t, s_t),
            }
        }
        self.hs_legacy(h_t, s_t)
    }

    /// Upstream HS_flash's legacy "sad path": iterate on T with the
    /// (Smolar,T) flash inside — scan upward from Tmin for the FIRST
    /// adjacent sign change of h(SmolarT(s,T)) - h and refine there
    /// (skipping temperatures where the inner solve fails), falling back to
    /// a whole-range Brent.
    fn hs_legacy(&self, hmolar: f64, smolar: f64) -> Result<HeosState> {
        let fluid = self.fluid();
        let eval_resid = |t: f64| -> Result<(f64, HeosState)> {
            let st = self.smolar_t_state(smolar, t)?;
            Ok((self.state_hmolar(&st) - hmolar, st))
        };

        // Find minimum temperature (upstream `Ttriple()` = sat_min_liquid.T)
        let mut tmin = fluid.eos.sat_min_liquid.t;
        let mut rmin = f64::NAN;
        let mut good_tmin = false;
        while !good_tmin {
            if let Ok((r, _)) = eval_resid(tmin) {
                rmin = r;
                good_tmin = true;
            } else {
                tmin += 0.5;
            }
            if tmin > fluid.eos.t_max {
                return Err(Error::Value("Cannot find good Tmin".into()));
            }
        }

        // Find maximum temperature (a little above Tmax, as upstream)
        let mut tmax = fluid.eos.t_max * 1.01;
        let mut rmax = f64::NAN;
        let mut good_tmax = false;
        while !good_tmax {
            if let Ok((r, _)) = eval_resid(tmax) {
                rmax = r;
                good_tmax = true;
            } else {
                tmax -= 0.1;
            }
            if tmax < tmin {
                return Err(Error::Value("Cannot find good Tmax".into()));
            }
        }
        if rmin * rmax > 0.0 && rmax.abs() < rmin.abs() {
            return Err(Error::Value(format!(
                "HS inputs correspond to temperature above maximum temperature of EOS [{} K]",
                fluid.eos.t_max
            )));
        }

        // Scan upward for the FIRST sign change between ADJACENT valid
        // samples; refine that sub-interval at upstream's 30-bit tolerance.
        if rmin.is_finite() {
            let nscan = 50;
            let mut t_prev = tmin;
            let mut f_prev = rmin;
            let mut i_prev = 0;
            for i in 1..=nscan {
                let tt = tmin + (tmax - tmin) * f64::from(i) / f64::from(nscan);
                let Ok((ft, _)) = eval_resid(tt) else {
                    continue;
                };
                if !ft.is_finite() {
                    continue;
                }
                if i == i_prev + 1 && f_prev * ft <= 0.0 {
                    let root = bisect_bits(
                        |t| eval_resid(t).map_or(f64::NAN, |x| x.0),
                        t_prev,
                        tt,
                        f_prev,
                        ft,
                        30,
                        100,
                    );
                    // Leave the state at the converged root (upstream
                    // resid.call on the bracket midpoint).
                    let (_, st) = eval_resid(root)?;
                    return Ok(st);
                }
                t_prev = tt;
                f_prev = ft;
                i_prev = i;
            }
        }
        // No interior sign change found: fall back to bracketing the whole
        // range (upstream Brent(resid, Tmin, Tmax, DBL_EPSILON, 1e-10, 100)).
        let root = brent(
            |t| eval_resid(t).map_or(f64::NAN, |x| x.0),
            tmin,
            tmax,
            f64::EPSILON,
            1e-10,
            100,
        )?;
        let (_, st) = eval_resid(root)?;
        Ok(st)
    }
}
