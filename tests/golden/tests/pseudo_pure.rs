//! Pseudo-pure flash goldens (tier-2 deferral, final item): PT
//! (liquid/gas/supercritical), QT at the only defined qualities (0/1), and
//! PQ across the temperature glide, through `props_si` for all six
//! pseudo-pure fluids; plus upstream's error conditions. The ancillary
//! curves are exact evaluations and the density solves run at upstream's
//! 1e-8 residuals — 1e-8 policy with caloric outputs measured against the
//! thermal scale where they cross zero.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn pseudo_pure_flashes_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/pseudo_pure.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 330);

    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
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
        // Caloric outputs cross zero at the arbitrary reference state where
        // pure relative error is ill-posed; measure them against the
        // thermal scale (R*Tc for H/U, R for S) when |expected| falls
        // below it (the established policy).
        let fluid = registry[rec.fluid.as_str()];
        let rg = fluid.eos.gas_constant;
        let tc = fluid.states.critical.t;
        let denom = match rec.out.as_str() {
            "Hmolar" | "Umolar" => rec.expected.abs().max(rg * tc),
            "Smolar" => rec.expected.abs().max(rg),
            _ => rec.expected.abs(),
        };
        let rel = (actual - rec.expected).abs() / denom;
        if rel > 1e-8 || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert_eq!(fluids.len(), 6, "all pseudo-pure fluids covered");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Upstream's pseudo-pure error conditions.
#[test]
fn pseudo_pure_error_conditions() {
    use rustprop::Error;
    // QT with fractional quality: "Two-phase quality is not defined".
    match props_si("D", "T", 100.0, "Q", 0.5, "HEOS::Air").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("quality must be equal to 0 or 1"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // PT inside the ancillary band: two-phase inputs unsupported.
    let psat = props_si("P", "T", 100.0, "Q", 0.0, "HEOS::Air").unwrap();
    match props_si("D", "T", 100.0, "P", psat * 0.999, "HEOS::Air").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("Two-phase inputs not supported for pseudo-pure"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Unported pairs stay loud ((D,T) still routes through the guard;
    // (H,P)/(P,S)/(P,U)/(D,P) are served since the classic-ancillary
    // flash port).
    assert!(matches!(
        props_si("P", "Dmolar", 5000.0, "T", 100.0, "HEOS::Air").unwrap_err(),
        Error::NotImplemented(_)
    ));
}
