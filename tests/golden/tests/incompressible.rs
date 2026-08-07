//! Incompressible goldens (PLAN.md 8.1): PT states, the DmassP/HmassP/
//! PSmass back-flashes, QT (Q = 0), trivials, and transport, through
//! `props_si` for 12 pure and 10 concentration-bearing fluids; plus
//! upstream's error conditions. Direct evaluations reproduce the oracle
//! bit-for-bit; the Brent back-flashes run at 1e-9.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn incompressible_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/incompressible.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 935);

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
            &format!("INCOMP::{}", rec.fluid),
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{} {}: error {e}", rec.fluid, rec.id()));
                continue;
            }
        };
        // H/S/U cross zero at the hard-coded reference state; measure
        // against a thermal scale there (c*20K ~ 8e4 J/kg for H/U, c/15 for
        // S — use fixed floors matched to the data's magnitudes).
        let denom = match rec.out.as_str() {
            "H" | "U" => rec.expected.abs().max(1e4),
            "S" => rec.expected.abs().max(1e2),
            _ => rec.expected.abs(),
        };
        if actual == rec.expected {
            continue; // exact match (incl. genuine zeros, e.g. LiBr L)
        }
        let rel = (actual - rec.expected).abs() / denom;
        if rel > 1e-9 || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert_eq!(fluids.len(), 22, "all fixture fluids covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Upstream's INCOMP error conditions.
#[test]
fn incompressible_error_conditions() {
    use rustprop::Error;
    // Solution without a bracket gets the PropsSI default x = 1.0, which
    // fails the composition range check.
    match props_si("D", "T", 300.0, "P", 1e5, "INCOMP::MEG").unwrap_err() {
        Error::Value(msg) => assert!(msg.contains("is not between"), "unexpected message: {msg}"),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Out-of-range temperature.
    match props_si("D", "T", 1000.0, "P", 1e5, "INCOMP::DowQ").unwrap_err() {
        Error::Value(msg) => assert!(msg.contains("is not between"), "unexpected message: {msg}"),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Q must be exactly 0.
    match props_si("P", "T", 400.0, "Q", 0.5, "INCOMP::DowQ").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("saturated liquid, Q=0"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Molar outputs are undefined for the mass-based backend.
    assert!(matches!(
        props_si("Dmolar", "T", 300.0, "P", 1e5, "INCOMP::DowQ").unwrap_err(),
        Error::NotImplemented(_)
    ));
    // Unsupported input pair.
    match props_si("T", "P", 1e5, "Umass", 1e4, "INCOMP::DowQ").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("is not yet supported"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Unknown fluid.
    match props_si("D", "T", 300.0, "P", 1e5, "INCOMP::NotAFluid").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("was not found in string_to_index_map"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
}
