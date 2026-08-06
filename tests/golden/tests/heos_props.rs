//! Property-level HEOS goldens (PLAN.md 4.2/4.7): p, h, s, u, cv, cp, w, g
//! at (T, rhomolar) must match `DmolarT_INPUTS` updates of the wheel at
//! rel <= 1e-9 (policy; observed deviations are far smaller since the terms
//! match at 1e-13/1e-12), for every ported fluid.

use rustprop_golden_tests::heos_fluids;
use rustprop_heos::HelmholtzEos;
use std::path::Path;

#[derive(serde::Deserialize)]
struct PropRecord {
    fluid: String,
    t: f64,
    rhomolar: f64,
    out: String,
    expected: f64,
}

#[test]
fn single_phase_properties_match_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_props.jsonl"));
        let text = std::fs::read_to_string(&path).unwrap();
        let records: Vec<PropRecord> = text
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        // Water: 12 points; parameterized fluids: 13 points x 8 accessors.
        let expected_count = if stem == "water" { 96 } else { 104 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let eos = HelmholtzEos::new(fluid);
        for rec in &records {
            assert_eq!(rec.fluid, name);
            let actual = match rec.out.as_str() {
                "p" => eos.pressure(rec.t, rec.rhomolar),
                "hmolar" => eos.hmolar(rec.t, rec.rhomolar),
                "smolar" => eos.smolar(rec.t, rec.rhomolar),
                "umolar" => eos.umolar(rec.t, rec.rhomolar),
                "cvmolar" => eos.cvmolar(rec.t, rec.rhomolar),
                "cpmolar" => eos.cpmolar(rec.t, rec.rhomolar),
                "speed_sound" => eos.speed_sound(rec.t, rec.rhomolar),
                "gibbsmolar" => eos.gibbsmolar(rec.t, rec.rhomolar),
                other => panic!("unknown property {other}"),
            };
            let rel = ((actual - rec.expected) / rec.expected).abs();
            if rel > 1e-9 || rel.is_nan() {
                failures.push(format!(
                    "{name} {}(T={}, rhomolar={}): actual {actual:e}, expected {:e}, rel {rel:e}",
                    rec.out, rec.t, rec.rhomolar, rec.expected
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} property records failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
