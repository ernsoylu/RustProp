//! Ancillary goldens (PLAN.md 4.3/4.7): classic pS/rhoL/rhoV fits vs the
//! wheel's `saturation_ancillary`, and each fluid's superancillary vs its
//! embedded extended-precision check points.

use rustprop_golden_tests::heos_fluids;
use std::path::Path;

#[derive(serde::Deserialize)]
struct AncRecord {
    fluid: String,
    t: f64,
    out: String,
    q: u8,
    expected: f64,
}

#[test]
fn classic_ancillaries_match_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("fixtures/heos_{stem}_ancillary.jsonl"));
        let text = std::fs::read_to_string(&path).unwrap();
        let records: Vec<AncRecord> = text
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        let expected_count = if stem == "water" { 33 } else { 36 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        for rec in &records {
            assert_eq!(rec.fluid, name);
            let anc = match (rec.out.as_str(), rec.q) {
                ("P", _) => &fluid.ancillaries.p_s,
                ("Dmolar", 0) => &fluid.ancillaries.rho_l,
                ("Dmolar", 1) => &fluid.ancillaries.rho_v,
                other => panic!("unexpected record {other:?}"),
            };
            let actual = rustprop_heos::ancillary::evaluate(anc, rec.t);
            let rel = ((actual - rec.expected) / rec.expected).abs();
            if rel > 1e-12 || rel.is_nan() {
                failures.push(format!(
                    "{name} {}(Q={}, T={}): actual {actual:e}, expected {:e}, rel {rel:e}",
                    rec.out, rec.q, rec.t, rec.expected
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The superancillary data embeds extended-precision check points together
/// with fastchebpure's own double-precision-eval-over-multiprecision ratios.
/// Our Clenshaw evaluation must reproduce those ratios essentially exactly,
/// for every ported fluid.
#[test]
fn superancillaries_reproduce_embedded_check_points() {
    for (name, _stem, fluid) in heos_fluids() {
        let sa = fluid
            .eos
            .superancillary
            .as_ref()
            .expect("every ported fluid has a superancillary");
        assert!(!sa.check_points.is_empty(), "{name}: no check points");
        for cp in sa.check_points {
            let cases = [
                ('P', 0u8, cp.p, cp.p_ratio, "p"),
                ('D', 0, cp.rho_l, cp.rho_l_ratio, "rhoL"),
                ('D', 1, cp.rho_v, cp.rho_v_ratio, "rhoV"),
            ];
            for (k, q, mp, expected_ratio, what) in cases {
                let ratio = rustprop_heos::superancillary::eval_sat(sa, cp.t, k, q) / mp;
                assert!(
                    (ratio - expected_ratio).abs() <= 1e-14,
                    "{name} {what}(T={}): our ratio {ratio:.16}, fastchebpure {expected_ratio:.16}",
                    cp.t
                );
                assert!(
                    (ratio - 1.0).abs() < 1e-10,
                    "{name} {what}(T={}): ratio {ratio} far from 1",
                    cp.t
                );
            }
        }
    }
}
