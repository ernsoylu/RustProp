//! Cubic-superancillary curve goldens (PLAN.md 7.2): saturated p/rhoL/rhoV
//! from the wheel's `update_QT_pure_superanc` across the dome, against the
//! ported Chebyshev tables directly (unit-level — no flash in between).

use rustprop_cubics::{CubicEos, CubicKind};
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn cubic_superancillary_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/cubic_superanc.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 432);

    let mut engines: std::collections::HashMap<(String, String), CubicEos> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &records {
        let key = (rec.backend.clone(), rec.fluid.clone());
        let eos = engines.entry(key).or_insert_with(|| {
            let kind = if rec.backend == "SRK" {
                CubicKind::Srk
            } else {
                CubicKind::PengRobinson
            };
            let upper = rec.fluid.to_uppercase();
            let data = rustprop_data::cubics::CUBIC_FLUIDS
                .iter()
                .find(|f| f.name == upper || f.aliases.contains(&upper.as_str()))
                .expect("fixture fluid in the cubic table");
            CubicEos::new(kind, data)
        });
        let (p, rho_l, rho_v) = match eos.superanc_sat(rec.val1) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!(
                    "{}::{} {}: error {e}",
                    rec.backend,
                    rec.fluid,
                    rec.id()
                ));
                continue;
            }
        };
        let actual = match rec.out.as_str() {
            "sa_p" => p,
            "sa_rhoL" => rho_l,
            "sa_rhoV" => rho_v,
            other => panic!("unknown output {other}"),
        };
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > 1e-12 || rel.is_nan() {
            failures.push(format!(
                "{}::{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.backend,
                rec.fluid,
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
