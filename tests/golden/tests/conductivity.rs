//! Conductivity goldens (PLAN.md 6.1): `L` for all 58 registry fluids with
//! a conductivity model — structured, hardcoded, and ECS — across PT states
//! (including the near-critical region where the Olchowy-Sengers
//! enhancement dominates) and the saturation curve; plus upstream's error
//! conditions.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn structured_conductivity_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/conductivity.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 571);

    let mut failures = Vec::new();
    let mut fluids = std::collections::HashSet::new();
    for rec in &records {
        fluids.insert(rec.fluid.clone());
        let actual = match props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("HEOS::{}", rec.fluid),
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{} {}: error {e}", rec.fluid, rec.id()));
                continue;
            }
        };
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > 1e-8 || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert_eq!(fluids.len(), 58, "all evaluable fluids covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Error-condition parity: every conductivity model class is now ported —
/// the remaining upstream errors are a missing TRANSPORT block and
/// conformal-state-solver failures on deep two-phase mixture states.
#[test]
fn conductivity_error_conditions() {
    use rustprop::Error;
    // No TRANSPORT block at all.
    assert!(matches!(
        props_si("L", "T", 300.0, "P", 1e5, "Acetone").unwrap_err(),
        Error::Value(_)
    ));
    // ECS at a deep two-phase mixture density: upstream's conformal state
    // solver cannot match the reference and throws
    // "Conformal state solver failed; error was Not able to get a solution".
    let err = props_si("L", "T", 243.7975002, "Q", 0.5, "R32").unwrap_err();
    match err {
        Error::Value(msg) => assert!(
            msg.contains("Conformal state solver failed"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
}
