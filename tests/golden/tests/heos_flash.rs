//! Flash-pair goldens (PLAN.md 4.6/4.7): general-quality (T,Q)/(P,Q), (D,T)
//! including two-phase states, and (H,P)/(P,Smolar)/(D,P) including
//! two-phase and backward-T outputs, for every ported fluid. (P,X) pairs run
//! at the 1e-8 policy — upstream itself resolves their T at a 30-bit
//! (~1e-9) tolerance; everything else at 1e-9.

use rustprop_golden_tests::{heos_fluids, load_jsonl};
use rustprop_heos::flash_pt::PtFlash;
use rustprop_heos::flash_px::HeosState;
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
fn flash_pairs_match_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_flash.jsonl"));
        let records = load_jsonl(&path);
        let expected_count = if stem == "water" { 261 } else { 253 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let flash = PtFlash::new(fluid);
        for rec in &records {
            let pair = (rec.name1.as_str(), rec.name2.as_str());
            let (state, rtol) = match pair {
                ("T", "Q") => (flash.qt_state(rec.val1, rec.val2), 1e-9),
                ("P", "Q") => (flash.pq_state(rec.val1, rec.val2), 1e-9),
                ("Dmolar", "T") => (flash.dmolar_t_state(rec.val1, rec.val2), 1e-9),
                ("Hmolar", "P") => (flash.hmolar_p_state(rec.val1, rec.val2), 1e-8),
                ("Dmolar", "P") => (flash.dmolar_p_state(rec.val1, rec.val2), 1e-8),
                ("P", "Smolar") => (flash.p_smolar_state(rec.val1, rec.val2), 1e-8),
                other => panic!("unexpected pair {other:?}"),
            };
            let state = match state {
                Ok(s) => s,
                Err(e) => {
                    failures.push(format!("{name} {}: flash error {e}", rec.id()));
                    continue;
                }
            };
            let actual = output(&flash, &state, &rec.out);
            // Q = -1 sentinels compare exactly; T in the (P,X) pairs may
            // resolve only to the shared 30-bit tolerance.
            if rec.expected == -1.0 && rec.out == "Q" {
                if actual != -1.0 {
                    failures.push(format!("{name} {}: got {actual}, expected -1", rec.id()));
                }
                continue;
            }
            // Caloric outputs of the solver-resolved pairs cross zero at the
            // arbitrary reference state, where pure relative error is
            // ill-posed; measure them against the thermal scale (R*T for
            // H/U, R for S) when |expected| falls below it.
            let solver_pair = matches!(pair, ("Hmolar", "P") | ("P", "Smolar") | ("Dmolar", "P"));
            let denom = match (solver_pair, rec.out.as_str()) {
                (true, "Hmolar" | "Umolar") => {
                    rec.expected.abs().max(fluid.eos.gas_constant * state.t())
                }
                (true, "Smolar") => rec.expected.abs().max(fluid.eos.gas_constant),
                _ => rec.expected.abs(),
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
