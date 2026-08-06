//! All-fluids smoke suite (PLAN.md 4.8): a coarse reduced-coordinate grid —
//! saturation, PT (liquid/gas/supercritical), and one each of the
//! HP/PS/DT/HS pairs — for every pure fluid with a superancillary (130).
//! `#[ignore]`d to keep the default `cargo test` fast; CI runs it in the
//! scheduled/manually-triggered sweep job:
//! `cargo test -p rustprop-golden-tests --test heos_all_smoke -- --ignored`.

use rustprop_golden_tests::load_jsonl;
use rustprop_heos::flash_pt::PtFlash;
use rustprop_heos::flash_px::HeosState;
use std::collections::HashMap;
use std::path::Path;

fn output(flash: &PtFlash, state: &HeosState, out: &str) -> f64 {
    match out {
        "P" => state.p(),
        "T" => state.t(),
        "Dmolar" => state.rhomolar(),
        "Q" => state.q(),
        "Hmolar" => flash.state_hmolar(state),
        "Smolar" => flash.state_smolar(state),
        "Umolar" => flash.state_umolar(state),
        other => panic!("unknown output {other}"),
    }
}

#[test]
#[ignore = "all-fluids sweep: run explicitly or in the scheduled CI job"]
fn all_fluids_smoke_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/heos_all_smoke.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 1922);

    let registry: HashMap<&'static str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    let mut flashes: HashMap<String, PtFlash> = HashMap::new();
    let mut failures = Vec::new();
    let mut fluids_seen = std::collections::HashSet::new();
    for rec in &records {
        fluids_seen.insert(rec.fluid.clone());
        let flash = flashes.entry(rec.fluid.clone()).or_insert_with(|| {
            PtFlash::new(
                registry
                    .get(rec.fluid.as_str())
                    .unwrap_or_else(|| panic!("{} not in registry", rec.fluid)),
            )
        });
        let pair = (rec.name1.as_str(), rec.name2.as_str());
        let state = match pair {
            ("T", "Q") => flash.qt_state(rec.val1, rec.val2),
            ("T", "P") => {
                flash
                    .pt_flash(rec.val1, rec.val2)
                    .map(|(rho, phase)| HeosState::SinglePhase {
                        t: rec.val1,
                        p: rec.val2,
                        rhomolar: rho,
                        phase,
                        q: -1.0,
                    })
            }
            ("Hmolar", "P") => flash.hmolar_p_state(rec.val1, rec.val2),
            ("P", "Smolar") => flash.p_smolar_state(rec.val1, rec.val2),
            ("Dmolar", "T") => flash.dmolar_t_state(rec.val1, rec.val2),
            ("Hmolar", "Smolar") => flash.hmolar_smolar_state(rec.val1, rec.val2),
            other => panic!("unexpected pair {other:?}"),
        };
        let state = match state {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{} {}: flash error {e}", rec.fluid, rec.id()));
                continue;
            }
        };
        let actual = output(flash, &state, &rec.out);
        // Established policy guards (PLAN.md): solver-resolved pairs at 1e-8
        // with Dmolar at 2e-8 through the HS/legacy quality amplification;
        // P measured against the thermal scale rho*R*T of the state.
        let fluid_data = registry[rec.fluid.as_str()];
        let rtol = match (pair, rec.out.as_str()) {
            (("Hmolar", "Smolar"), "Dmolar") => 2e-8,
            _ => 1e-8,
        };
        let denom = if rec.out == "P" {
            rec.expected
                .abs()
                .max(fluid_data.eos.gas_constant * state.t() * state.rhomolar())
        } else {
            rec.expected.abs()
        };
        let rel = (actual - rec.expected).abs() / denom;
        if rel > rtol || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    println!("fluids covered: {}", fluids_seen.len());
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}
