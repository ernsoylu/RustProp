//! Flash-pair goldens (PLAN.md 4.6): general-quality (T,Q)/(P,Q), (D,T)
//! including two-phase states, and (H,P)/(P,Smolar) including two-phase and
//! backward-T outputs. (P,X) pairs run at the 1e-8 policy — upstream itself
//! resolves their T at a 30-bit (~1e-9) tolerance; everything else at 1e-9.

use rustprop_data::fluids::water::WATER;
use rustprop_golden_tests::load_jsonl;
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
fn water_flash_pairs_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/heos_water_flash.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 261);

    let flash = PtFlash::new(&WATER);
    let mut failures = Vec::new();
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
                failures.push(format!("{}: flash error {e}", rec.id()));
                continue;
            }
        };
        let actual = output(&flash, &state, &rec.out);
        // Q = -1 sentinels compare exactly; T in the (P,X) pairs may resolve
        // only to the shared 30-bit tolerance.
        if rec.expected == -1.0 && rec.out == "Q" {
            if actual != -1.0 {
                failures.push(format!("{}: got {actual}, expected -1", rec.id()));
            }
            continue;
        }
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > rtol || rel.is_nan() {
            failures.push(format!(
                "{}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.id(),
                rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}
