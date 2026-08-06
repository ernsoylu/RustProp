//! IF97 golden tests (PLAN.md 2.3): every record in `if97_water.jsonl` —
//! generated from the CoolProp 8.0.0 wheel — must match the port through the
//! facade's `PropsSI`-style dispatch.
//!
//! The algorithm and constants are identical to upstream; residual deviation
//! comes only from libm (exp/ln/pow) differences between the manylinux wheel
//! and Rust std, so the default tolerance is far below the PLAN.md policy.

use rustprop_core::Param;
use rustprop_golden_tests::{GoldenRecord, check, load_jsonl};
use std::path::Path;

// Max observed deviation is 5.6e-12 (surface tension near critical), from
// libm pow differences between the manylinux wheel and Rust std.
const DEFAULT_RTOL: f64 = 1e-11;

fn eval(rec: &GoldenRecord) -> Result<f64, rustprop_core::Error> {
    assert_eq!(rec.backend, "IF97");
    assert_eq!(rec.fluid, "Water");
    let out = Param::parse(&rec.out).unwrap_or_else(|| panic!("bad out {:?}", rec.out));
    let n1 = Param::parse(&rec.name1).unwrap_or_else(|| panic!("bad name1 {:?}", rec.name1));
    let n2 = Param::parse(&rec.name2).unwrap_or_else(|| panic!("bad name2 {:?}", rec.name2));
    rustprop::if97_api::props(out, n1, rec.val1, n2, rec.val2)
}

#[test]
fn if97_water_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/if97_water.jsonl");
    let records = load_jsonl(&path);
    assert!(
        records.len() > 300,
        "expected the full fixture set, got {}",
        records.len()
    );
    let mut failures = Vec::new();
    for rec in &records {
        match eval(rec) {
            Ok(actual) => {
                // Q = -1 sentinel records compare exactly.
                if rec.expected == -1.0 {
                    if actual != -1.0 {
                        failures.push(format!("{}: got {actual}, expected -1", rec.id()));
                    }
                } else if let Err(msg) = check(rec, actual, DEFAULT_RTOL) {
                    failures.push(msg);
                }
            }
            Err(e) => failures.push(format!("{}: error {e}", rec.id())),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} golden records failed:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}
