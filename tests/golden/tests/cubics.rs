//! Cubic-backend goldens (PLAN.md 7.1): `SRK::` / `PR::` PT states
//! (liquid/gas/supercritical), QT at Q = 0/1, PQ across quality, and the
//! trivial outputs, through `props_si` for 12 structurally diverse fluids;
//! plus upstream's error conditions. Direct evaluations reproduce the
//! oracle bit-for-bit; the saturation solves and their downstream outputs
//! run at the 1e-9 policy of the phase's `→ verify` clause.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn cubic_backends_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/cubics.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 1896);

    let mut failures = Vec::new();
    let mut combos = std::collections::HashSet::new();
    for rec in &records {
        combos.insert((rec.backend.clone(), rec.fluid.clone()));
        let actual = match props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("{}::{}", rec.backend, rec.fluid),
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!(
                    "{}::{} {}: error {e}",
                    rec.backend,
                    rec.fluid,
                    rec.id()
                ));
                continue;
            }
        };
        // Caloric outputs cross zero at the reference state; measure against
        // the thermal scale there (established policy: R*Tc for H/U, R for S).
        let denom = match rec.out.as_str() {
            "Hmolar" | "Umolar" => rec.expected.abs().max(8.314 * 300.0),
            "Smolar" => rec.expected.abs().max(8.314),
            _ => rec.expected.abs(),
        };
        let rel = (actual - rec.expected).abs() / denom;
        if rel > 1e-9 || rel.is_nan() {
            failures.push(format!(
                "{}::{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.backend,
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert_eq!(combos.len(), 24, "12 fluids x 2 backends covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Upstream's cubic error conditions.
#[test]
fn cubic_error_conditions() {
    use rustprop::Error;
    // Delegated input pairs die on the fabricated fluid's unset ancillaries.
    match props_si("T", "Hmolar", 1e4, "P", 1e6, "SRK::n-Propane").unwrap_err() {
        Error::Value(msg) => assert!(msg.contains("type not set"), "unexpected message: {msg}"),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Two-phase (D,T) caloric reads reproduce upstream's broken-sub-state
    // throw (P/T/D/Q remain readable — golden-covered).
    match props_si("Hmolar", "Dmolar", 5000.0, "T", 300.0, "SRK::n-Propane").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("calc_alpha0_deriv_nocache returned invalid number"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Unknown cubic fluid.
    assert!(matches!(
        props_si("D", "T", 300.0, "P", 1e5, "SRK::NotAFluid").unwrap_err(),
        Error::Value(_)
    ));
    // Transport is not available for cubics upstream either.
    assert!(props_si("V", "T", 300.0, "P", 1e5, "SRK::n-Propane").is_err());
}
