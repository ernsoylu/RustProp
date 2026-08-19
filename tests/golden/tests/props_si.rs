//! PropsSI string-API goldens (PLAN.md 5.1/5.2): fluid-key resolution
//! exactly as upstream registers strings, and `props_si` over mass-basis
//! aliases/inputs, swapped pairs, input echo, and trivial outputs — plus
//! the error conditions upstream raises.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[derive(serde::Deserialize)]
struct ResolutionRecord {
    query: String,
    /// Canonical INFO.NAME the wheel resolves to; None = the wheel rejects.
    name: Option<String>,
}

#[test]
fn fluid_resolution_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/fluid_resolution.jsonl");
    let text = std::fs::read_to_string(&path).unwrap();
    let records: Vec<ResolutionRecord> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(records.len(), 639);

    let mut failures = Vec::new();
    for rec in &records {
        match (rustprop::props_api::resolve_fluid(&rec.query), &rec.name) {
            (Ok(data), Some(name)) => {
                if data.name != name {
                    failures.push(format!(
                        "{:?}: resolved to {:?}, wheel says {:?}",
                        rec.query, data.name, name
                    ));
                }
            }
            (Err(_), None) => {}
            (Ok(data), None) => {
                failures.push(format!(
                    "{:?}: resolved to {:?}, wheel rejects",
                    rec.query, data.name
                ));
            }
            (Err(e), Some(name)) => {
                failures.push(format!(
                    "{:?}: rejected ({e}), wheel resolves to {name:?}",
                    rec.query
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} resolution mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn props_si_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/props_si.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 196);

    let mut failures = Vec::new();
    for rec in &records {
        let fluid_string = format!("{}::{}", rec.backend, rec.fluid);
        let actual = match props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &fluid_string,
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{}: error {e}", rec.id()));
                continue;
            }
        };
        // Trivial outputs and input echoes are direct data reads — exact to
        // a multiply; state-derived outputs follow the established solver
        // policies (1e-8 with the (P,X)/HS pairs' documented guards).
        let rtol = 1e-8;
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

/// The bare fluid name defaults to the HEOS backend (upstream backend "?").
#[test]
fn bare_fluid_name_defaults_to_heos() {
    let a = props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water").unwrap();
    let b = props_si("Dmolar", "T", 300.0, "P", 101325.0, "HEOS::Water").unwrap();
    // Bitwise despite the value being an iteratively solved PT density: the
    // two strings differ only in the backend prefix, and after parsing both
    // run the identical flash on the identical inputs. The claim under test
    // is "same route", which is a bit-equality claim.
    assert_eq!(a, b);
}

/// Error conditions mirror upstream's (conditions, not message strings).
#[test]
fn props_si_error_conditions() {
    use rustprop::Error;
    let err = |r: rustprop::Result<f64>| r.unwrap_err();

    // Unknown fluid key
    assert!(matches!(
        err(props_si("T", "P", 1e5, "Q", 0.0, "Watr")),
        Error::Value(_)
    ));
    // Valid params, invalid combination
    assert!(matches!(
        err(props_si("Dmolar", "T", 300.0, "T", 400.0, "Water")),
        Error::Value(_)
    ));
    // Invalid input name with a non-trivial output
    assert!(matches!(
        err(props_si("Dmolar", "nonsense", 300.0, "P", 1e5, "Water")),
        Error::Value(_)
    ));
    // ... but trivial outputs never need the state
    assert!(props_si("Tcrit", "", 0.0, "", 0.0, "Water").is_ok());
    // Vapor quality out of [0, 1]
    assert!(matches!(
        err(props_si("T", "P", 1e5, "Q", 1.5, "Water")),
        Error::OutOfRange(_)
    ));
    // Unknown backend prefix
    assert!(matches!(
        err(props_si("T", "P", 1e5, "Q", 0.0, "REFPROP::Water")),
        Error::Value(_)
    ));
    // Unknown output parameter
    assert!(matches!(
        err(props_si("nonsense", "T", 300.0, "P", 1e5, "Water")),
        Error::Value(_)
    ));
    // IF97 routes through the string API too
    let h = props_si("H", "T", 300.0, "P", 101325.0, "IF97::Water").unwrap();
    assert!(((h - 112665.04341853978) / h).abs() < 1e-11);
}
