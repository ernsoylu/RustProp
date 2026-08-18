//! PT-flash goldens (PLAN.md 4.5/4.7): the single-phase density solver and
//! the properties through it must match `PT_INPUTS` states of the wheel at
//! the 1e-9 policy across liquid, vapor, and supercritical states, for
//! every ported fluid.
//!
//! The density is asserted BITWISE — all 218 density records, all 12 fluids —
//! because that is the premise everything else here rests on: upstream assigns
//! the solved root to `_rhomolar` and the port finds the same root. Everything
//! served off the reduced-Helmholtz derivatives instead carries a per-record
//! STALE-CACHE ALLOWANCE — see `stale_cache_allowance` for the upstream
//! mechanism and its bound.

use rustprop_golden_tests::{heos_fluids, load_jsonl};
use rustprop_heos::flash_pt::PtFlash;
use std::path::Path;

/// Upstream's `PT_flash` reports the root but serves properties off the
/// density solver's LAST TRIAL.
///
/// `FlashRoutines.cpp:336` closes with `HEOS._rhomolar = HEOS.solver_rho_Tp(…)`
/// — it assigns the root and nothing else. Every residual evaluation inside
/// that solve went through `SolverTPResid::call`, which mutates the same
/// backend (`HEOS->update_DmolarT_direct(rhomolar, T)`,
/// `HelmholtzEOSMixtureBackend.cpp:2803`) and leaves `_dalphar_dDelta` & co.
/// cached at that trial density. The property calculators then recompute
/// `_delta = _rhomolar / _reducing.rhomolar` FRESH but read the residual
/// derivatives from those caches (`calc_hmolar` at :3200, `calc_cpmolar` at
/// :3326, `calc_speed_sound` at :3368) — so a wheel PT state mixes a
/// root-density delta with last-iterate derivatives and is internally
/// inconsistent. The port has no such cache: it evaluates at the root, which
/// is why its density matches BITWISE while h/s/cp/w do not.
///
/// The displacement is bounded. Upstream's solvers return one Newton step
/// past the iterate whose residual passed, so
/// `|Δln ρ| = |(p_eos − p)/p| / |dln p/dln ρ| < ftol / |dln p/dln ρ|` with
/// `ftol = 1e-8` (`solver_rho_Tp`, :2937/:2944/:2951). Verified exactly on
/// three near-critical records: backing the wheel's cached
/// `dalphar_dtau_constdelta` out to a trial density reproduces its cached
/// `dalphar_ddelta_consttau` to ≤2.7e-13 and gives Δln ρ = +2.215e-8
/// (Fluorine), +1.474e-8 (R125), −1.230e-8 (R22) — each equal to that
/// record's last residual divided by `dln p/dln ρ` to four digits.
///
/// This returns the output's spread over `ρ(1 ± ftol/|dln p/dln ρ|)`, times
/// `MIXING_MARGIN`: the wheel's formulas mix the fresh delta with the stale
/// derivatives, so a pure displacement is a proxy rather than an identity
/// (observed ratio of the true move to this proxy: 0.56 for Fluorine Cp,
/// 1.32 for R22 Cp — the worst over the whole suite is printed on the log).
const MIXING_MARGIN: f64 = 2.0;

fn stale_cache_allowance(flash: &PtFlash, t: f64, rho: f64, eval: impl Fn(f64) -> f64) -> f64 {
    let tau = flash.eos.t_reducing / t;
    let delta = rho / flash.eos.rhomolar_reducing;
    let ar = flash.eos.alphar_all(tau, delta);
    // dln p / dln rho = (1 + 2 δ αr_δ + δ² αr_δδ) / (1 + δ αr_δ)
    let dlnp_dlnrho =
        (1.0 + 2.0 * delta * ar.d10 + delta * delta * ar.d20) / (1.0 + delta * ar.d10);
    let eps = 1e-8 / dlnp_dlnrho.abs();
    let here = eval(rho);
    let spread = (eval(rho * (1.0 + eps)) - here)
        .abs()
        .max((eval(rho * (1.0 - eps)) - here).abs());
    MIXING_MARGIN * spread
}

#[test]
fn pt_states_match_upstream() {
    let mut failures = Vec::new();
    let mut worst_mix: (f64, String) = (0.0, String::new());
    let mut density_records = 0usize;
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_pt.jsonl"));
        let records = load_jsonl(&path);
        let expected_count = if stem == "water" { 100 } else { 90 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let flash = PtFlash::new(fluid);
        let mut worst: (f64, String) = (0.0, String::new());
        for rec in &records {
            assert_eq!((rec.name1.as_str(), rec.name2.as_str()), ("T", "P"));
            let (t, p) = (rec.val1, rec.val2);
            let (rho, _phase) = flash
                .pt_flash(t, p)
                .unwrap_or_else(|e| panic!("{name} {}: {e}", rec.id()));
            let eval = |r: f64| match rec.out.as_str() {
                "Dmolar" => r,
                "Hmolar" => flash.eos.hmolar(t, r),
                "Smolar" => flash.eos.smolar(t, r),
                "Cpmolar" => flash.eos.cpmolar(t, r),
                "A" => flash.eos.speed_sound(t, r),
                other => panic!("unknown output {other}"),
            };
            let actual = eval(rho);
            // The density is the root upstream assigns and the port solves to
            // the same double, so it is asserted BITWISE rather than at a
            // tolerance — if that ever stops holding, the stale-cache
            // allowance below loses the premise it is derived from and this
            // says so directly instead of absorbing it.
            if rec.out == "Dmolar" {
                density_records += 1;
                if actual != rec.expected {
                    failures.push(format!(
                        "{name} {}: density {actual:e} is not BITWISE the wheel's {:e} \
                         (rel {:e}) — the stale-cache allowance assumes it is",
                        rec.id(),
                        rec.expected,
                        ((actual - rec.expected) / rec.expected).abs()
                    ));
                }
                continue;
            }
            let rtol = {
                let allow = stale_cache_allowance(&flash, t, rho, eval) / rec.expected.abs();
                1e-9f64.max(allow)
            };
            let rel = ((actual - rec.expected) / rec.expected).abs();
            if rel > worst.0 {
                worst = (rel, rec.id());
            }
            // How much of the allowance a record actually consumes: the
            // measured justification for MIXING_MARGIN.
            if rtol > 1e-9 {
                let used = MIXING_MARGIN * rel / rtol;
                if used > worst_mix.0 {
                    worst_mix = (used, format!("{name} {}", rec.id()));
                }
            }
            if rel > rtol || rel.is_nan() {
                failures.push(format!(
                    "{name} {}: actual {actual:e}, expected {:e}, rel {rel:e}, rtol {rtol:e}",
                    rec.id(),
                    rec.expected
                ));
            }
        }
        println!("{name}: worst deviation {:e} at {}", worst.0, worst.1);
    }
    println!(
        "{density_records} density records bitwise; worst stale-cache mixing ratio \
         {:.3} (margin {MIXING_MARGIN}) at {}",
        worst_mix.0, worst_mix.1
    );
    // The bitwise premise is only meaningful if the records are actually there.
    assert_eq!(density_records, 218, "density records checked");
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Condition parity: near-saturation PT input errors like upstream, and the
/// critical-point short-circuit engages.
#[test]
fn pt_guards() {
    let water = heos_fluids()[0].2;
    let flash = PtFlash::new(water);
    // psat(400 K) ~ 245769 Pa: within 1e-6 relative must error
    let psat = flash.sat().qt_flash(400.0, 1.0).unwrap().p;
    assert!(flash.pt_flash(400.0, psat * (1.0 + 1e-8)).is_err());
    // Critical short-circuit — upstream returns `rhomolar_critical()`,
    // the superancillary's NUMERICAL critical density.
    let (rho, phase) = flash.pt_flash(647.096, 22064000.0).unwrap();
    assert_eq!(rho, water.eos.superancillary.as_ref().unwrap().rho_crit_num);
    assert_eq!(phase, rustprop_core::params::Phase::CriticalPoint);
}
