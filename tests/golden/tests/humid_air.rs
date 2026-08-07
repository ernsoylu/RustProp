//! Humid-air goldens (PLAN.md 9.1): the full `HAPropsSI` output set over
//! (T, P, humidity) grids including sub-freezing ice paths, plus the
//! T-iterating inverse triples, against the wheel. Direct evaluations run
//! at 1e-9; the finite-difference outputs (cp/cv/speed/isentropic) and the
//! iterative inverses at 1e-8.

use rustprop::ha_props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn humid_air_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/humid_air.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 897);

    let mut failures = Vec::new();
    for rec in &records {
        let actual = match ha_props_si(
            &rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &rec.name3, rec.val3,
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{}: error {e}", rec.id()));
                continue;
            }
        };
        let rtol = match rec.out.as_str() {
            "cp" | "cp_ha" | "CV" | "speed_of_sound" | "isentropic_exponent" => 1e-8,
            _ if rec.name1 != "T" => 1e-8, // inverse triples resolve T iteratively
            _ => 1e-9,
        };
        if actual == rec.expected {
            continue;
        }
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > rtol || rel.is_nan() {
            failures.push(format!(
                "{} [{} {} {}]: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.id(),
                rec.name3,
                rec.val3,
                rec.out,
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

/// Upstream's error conditions (as `Err` — the logged deviation from
/// upstream's swallow-into-inf transport).
#[test]
fn humid_air_error_conditions() {
    use rustprop::Error;
    // No pressure input.
    match ha_props_si("H", "T", 300.0, "W", 0.01, "R", 0.5).unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("Pressure must be one of the inputs"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Duplicate non-pressure inputs.
    match ha_props_si("H", "P", 101325.0, "T", 300.0, "T", 300.0).unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("cannot be the same"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Unknown parameter name.
    match ha_props_si("H", "P", 101325.0, "T", 300.0, "XYZ", 1.0).unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("was not understood to Name2Type"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Two water-content inputs that are not a valid pairing.
    match ha_props_si("T", "P", 101325.0, "W", 0.01, "D", 285.0).unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("cannot provide two inputs that are both water-content"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Out-of-range input.
    assert!(ha_props_si("H", "T", 700.0, "P", 101325.0, "R", 0.5).is_err());
    // Trivial echo happens BEFORE validation (upstream's order).
    assert_eq!(
        ha_props_si("H", "H", 5.0, "T", 300.0, "W", 0.01).unwrap(),
        5.0
    );
}
