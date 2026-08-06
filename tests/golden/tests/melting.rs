//! Melting-line goldens (tier-2 deferral): `p(T)`, `T(p)`, and the
//! aggregate limits against the wheel's `AbstractState.melting_line` for
//! all 29 pure fluids with a curve, plus the PT below-melt error condition.
//! Direct evaluations run at 1e-12; the Brent-inverted `T(p)` of the
//! polynomial families at 1e-10.

use rustprop_golden_tests::load_jsonl;
use rustprop_heos::melting;
use std::path::Path;

#[test]
fn melting_line_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/melting.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 290);

    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    let mut failures = Vec::new();
    let mut fluids = std::collections::HashSet::new();
    for rec in &records {
        fluids.insert(rec.fluid.clone());
        let fluid = *registry
            .get(rec.fluid.as_str())
            .unwrap_or_else(|| panic!("{} not in registry", rec.fluid));
        let ml = fluid
            .ancillaries
            .melting_line
            .as_ref()
            .expect("fixture fluid has a melting line");
        let (actual, rtol) = match rec.out.as_str() {
            "melt_pmin" => (Ok(melting::p_min(ml)), 1e-12),
            "melt_pmax" => (Ok(melting::p_max(ml)), 1e-12),
            "melt_p" => (melting::p_of_t(ml, rec.val1), 1e-12),
            "melt_T" => (melting::t_of_p(ml, rec.val1), 1e-10),
            other => panic!("unknown output {other}"),
        };
        let actual = match actual {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{} {}: error {e}", rec.fluid, rec.id()));
                continue;
            }
        };
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > rtol || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert_eq!(fluids.len(), 29, "all melting fluids covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// PT states below the melting temperature raise upstream's error.
#[test]
fn below_melt_error_condition() {
    use rustprop::{Error, props_si};
    for (fluid, t, p) in [
        ("Water", 255.0, 1e5),
        ("Water", 250.0, 1e8),
        ("Nitrogen", 40.0, 1e9),
        ("CarbonDioxide", 216.0, 1e6),
    ] {
        match props_si("D", "T", t, "P", p, &format!("HEOS::{fluid}")).unwrap_err() {
            Error::Value(msg) => assert!(
                msg.contains("below Tmelt(p)"),
                "{fluid}: unexpected message: {msg}"
            ),
            other => panic!("{fluid}: expected Value error, got {other:?}"),
        }
    }
    // Just above the melting line the state is a valid compressed liquid.
    assert!(props_si("D", "T", 270.0, "P", 1e8, "HEOS::Water").is_ok());
}
