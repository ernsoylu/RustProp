//! Viscosity goldens (PLAN.md 6.1, structured slice): `V` for the 20 fluids
//! whose model is fully structured (dilute + initial-density + higher-order
//! typed families), across PT states and the saturation curve including a
//! two-phase mixture-density evaluation; plus upstream's error conditions.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn structured_viscosity_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/viscosity.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 296);

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
        // The models are direct (T, rho, p) evaluations on top of the
        // 1e-9-verified flash states; 1e-8 policy.
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
    assert_eq!(fluids.len(), 37, "all evaluable fluids covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Error-condition parity for the not-yet-covered classes.
#[test]
fn viscosity_error_conditions() {
    use rustprop::Error;
    // No TRANSPORT block at all (upstream: model not available).
    assert!(matches!(
        props_si("V", "T", 300.0, "P", 1e5, "Acetone").unwrap_err(),
        Error::Value(_)
    ));
    // Model class not ported yet (Chung and rhosr-CS).
    assert!(matches!(
        props_si("V", "T", 350.0, "P", 1e6, "Cyclopentane").unwrap_err(),
        Error::NotImplemented(_)
    ));
    assert!(matches!(
        props_si("V", "T", 300.0, "P", 1e6, "R32").unwrap_err(),
        Error::NotImplemented(_)
    ));
}
