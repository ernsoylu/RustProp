//! Ancillary goldens (PLAN.md 4.3): classic pS/rhoL/rhoV fits vs the wheel's
//! `saturation_ancillary`, and the superancillary vs its embedded
//! extended-precision check points.

use rustprop_data::fluids::water::WATER;
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
fn water_classic_ancillaries_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/heos_water_ancillary.jsonl");
    let text = std::fs::read_to_string(&path).unwrap();
    let records: Vec<AncRecord> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(records.len(), 33);

    let mut failures = Vec::new();
    for rec in &records {
        assert_eq!(rec.fluid, "Water");
        let anc = match (rec.out.as_str(), rec.q) {
            ("P", _) => &WATER.ancillaries.p_s,
            ("Dmolar", 0) => &WATER.ancillaries.rho_l,
            ("Dmolar", 1) => &WATER.ancillaries.rho_v,
            other => panic!("unexpected record {other:?}"),
        };
        let actual = rustprop_heos::ancillary::evaluate(anc, rec.t);
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > 1e-12 || rel.is_nan() {
            failures.push(format!(
                "{}(Q={}, T={}): actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.out, rec.q, rec.t, rec.expected
            ));
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
/// Our Clenshaw evaluation must reproduce those ratios essentially exactly.
#[test]
fn water_superancillary_reproduces_embedded_check_points() {
    let sa = WATER
        .eos
        .superancillary
        .as_ref()
        .expect("Water has a superancillary");
    assert_eq!(sa.check_points.len(), 3);
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
                "{what}(T={}): our ratio {ratio:.16}, fastchebpure {expected_ratio:.16}",
                cp.t
            );
            assert!(
                (ratio - 1.0).abs() < 1e-10,
                "{what}(T={}): ratio {ratio} far from 1",
                cp.t
            );
        }
    }
}
