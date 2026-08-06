//! (Hmolar, Smolar) flash goldens (PLAN.md 4.6, final pair): the
//! superancillary HS path — two-phase screen + Qh==Qs solve, single-phase
//! cascade — must match the wheel across subcooled/cold liquid, two-phase,
//! gas, and supercritical states for every ported fluid. Solver-resolved
//! outputs run at the 1e-8 policy; the caloric zero-crossing guard measures
//! near-zero H/S outputs against the thermal scale, as in the flash suite.

use rustprop_golden_tests::{heos_fluids, load_jsonl};
use rustprop_heos::flash_pt::PtFlash;
use std::path::Path;

#[test]
fn hs_flash_matches_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_hs.jsonl"));
        let records = load_jsonl(&path);
        // The wheel's own HS flash rejects a few of the hardest source
        // states for these fluids (try_record filtered them).
        let expected_count = match stem {
            "r134a" | "methanol" => 44,
            "n_heptane" => 40,
            _ => 48,
        };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let flash = PtFlash::new(fluid);
        for rec in &records {
            assert_eq!(
                (rec.name1.as_str(), rec.name2.as_str()),
                ("Hmolar", "Smolar")
            );
            let state = match flash.hmolar_smolar_state(rec.val1, rec.val2) {
                Ok(s) => s,
                Err(e) => {
                    failures.push(format!("{name} {}: flash error {e}", rec.id()));
                    continue;
                }
            };
            let actual = match rec.out.as_str() {
                "T" => state.t(),
                "Dmolar" => state.rhomolar(),
                "P" => state.p(),
                "Q" => state.q(),
                other => panic!("unknown output {other}"),
            };
            // Q = -1 sentinels compare exactly.
            if rec.expected == -1.0 && rec.out == "Q" {
                if actual != -1.0 {
                    failures.push(format!("{name} {}: got {actual}, expected -1", rec.id()));
                }
                continue;
            }
            // Policy adjustments, both logged in PLAN.md:
            // - Dmolar runs at 2e-8: two-phase HS states resolve T at the
            //   shared deliberate 30-bit tolerance and the mixture density
            //   amplifies that by ~1/Q;
            // - P is measured against the thermal scale rho*R*T of the
            //   resolved state: compressed-liquid pressures are the
            //   cancellation residue of terms of that scale (observed: a
            //   0.23 Pa propane state agreeing to 4e-6 Pa absolute).
            let rtol = if rec.out == "Dmolar" { 2e-8 } else { 1e-8 };
            let denom = if rec.out == "P" {
                rec.expected
                    .abs()
                    .max(fluid.eos.gas_constant * state.t() * state.rhomolar())
            } else {
                rec.expected.abs()
            };
            let rel = (actual - rec.expected).abs() / denom;
            if rel > rtol || rel.is_nan() {
                failures.push(format!(
                    "{name} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                    rec.id(),
                    rec.expected
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
