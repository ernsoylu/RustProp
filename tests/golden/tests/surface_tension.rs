//! Surface tension goldens (PLAN.md 6.2): `I` along the saturation curve
//! for every pure fluid with a curve, through the string API; plus
//! upstream's error conditions (single-phase states, curveless fluids).

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn surface_tension_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/surface_tension.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 518);

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
        // The curve is a direct polynomial in (1 - T/Tc); QT resolves T
        // exactly, so agreement is essentially bitwise.
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > 1e-12 || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    println!("fluids covered: {}", fluids.len());
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Error-condition parity with upstream `calc_surface_tension`.
#[test]
fn surface_tension_error_conditions() {
    use rustprop::Error;
    // Single-phase input: only two-phase states carry a surface tension.
    assert!(matches!(
        props_si("I", "T", 300.0, "P", 101325.0, "Water").unwrap_err(),
        Error::Value(_)
    ));
    // A fluid without a curve errors even at a two-phase state
    // (upstream SurfaceTensionCorrelation::evaluate on an empty curve).
    assert!(matches!(
        props_si("I", "T", 250.0, "Q", 0.5, "Chlorine").unwrap_err(),
        Error::NotImplemented(_)
    ));
}
