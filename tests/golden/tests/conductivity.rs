//! Conductivity goldens (PLAN.md 6.1, structured slice): `L` for the 15
//! fluids whose dilute/residual/critical trio AND viscosity are fully
//! structured, across PT states (including the near-critical region where
//! the Olchowy-Sengers enhancement dominates) and the saturation curve;
//! plus upstream's error conditions.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn structured_conductivity_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/conductivity.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 306);

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
    assert_eq!(fluids.len(), 34, "all evaluable fluids covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Error-condition parity for the not-yet-covered classes and two-phase.
#[test]
fn conductivity_error_conditions() {
    use rustprop::Error;
    // No TRANSPORT block at all.
    assert!(matches!(
        props_si("L", "T", 300.0, "P", 1e5, "Acetone").unwrap_err(),
        Error::Value(_)
    ));
    // ECS conductivity class not ported yet.
    assert!(matches!(
        props_si("L", "T", 300.0, "P", 1e6, "R11").unwrap_err(),
        Error::NotImplemented(_)
    ));
    // Structured conductivity whose viscosity class is unported blocks the
    // OS enhancement (EthylBenzene's ECS viscosity).
    assert!(matches!(
        props_si("L", "T", 400.0, "P", 1e6, "EthylBenzene").unwrap_err(),
        Error::NotImplemented(_)
    ));
    // Two-phase input: cp (and the enhancement) are undefined.
    assert!(matches!(
        props_si("L", "T", 300.0, "Q", 0.5, "n-Propane").unwrap_err(),
        Error::Value(_)
    ));
}
