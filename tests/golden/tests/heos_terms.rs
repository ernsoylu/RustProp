//! Term-level HEOS goldens (PLAN.md 4.1): alphar/alpha0 and every tau/delta
//! derivative the wheel exposes must match the port at rel <= 1e-13 on a
//! single-phase (T, rhomolar) grid.

use rustprop_data::fluids::water::WATER;
use rustprop_heos::HelmholtzEos;
use std::path::Path;

#[derive(serde::Deserialize)]
struct TermRecord {
    fluid: String,
    t: f64,
    rhomolar: f64,
    out: String,
    expected: f64,
}

#[test]
fn water_helmholtz_terms_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/heos_water_terms.jsonl");
    let text = std::fs::read_to_string(&path).unwrap();
    let records: Vec<TermRecord> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert!(
        records.len() >= 200,
        "expected the full grid, got {}",
        records.len()
    );

    let eos = HelmholtzEos::new(&WATER);
    let mut failures = Vec::new();
    for rec in &records {
        assert_eq!(rec.fluid, "Water");
        let tau = eos.t_reducing / rec.t;
        let delta = rec.rhomolar / eos.rhomolar_reducing;
        let (residual, ideal) = (eos.alphar_all(tau, delta), eos.alpha0_all(tau, delta));
        let actual = match rec.out.as_str() {
            "alphar" => residual.d00,
            "dalphar_dDelta" => residual.d10,
            "dalphar_dTau" => residual.d01,
            "d2alphar_dDelta2" => residual.d20,
            "d2alphar_dDelta_dTau" => residual.d11,
            "d2alphar_dTau2" => residual.d02,
            "d3alphar_dDelta3" => residual.d30,
            "d3alphar_dDelta2_dTau" => residual.d21,
            "d3alphar_dDelta_dTau2" => residual.d12,
            "d3alphar_dTau3" => residual.d03,
            "alpha0" => ideal.d00,
            "dalpha0_dDelta" => ideal.d10,
            "dalpha0_dTau" => ideal.d01,
            "d2alpha0_dDelta2" => ideal.d20,
            "d2alpha0_dDelta_dTau" => ideal.d11,
            "d2alpha0_dTau2" => ideal.d02,
            "d3alpha0_dDelta3" => ideal.d30,
            "d3alpha0_dDelta2_dTau" => ideal.d21,
            "d3alpha0_dDelta_dTau2" => ideal.d12,
            "d3alpha0_dTau3" => ideal.d03,
            other => panic!("unknown accessor {other}"),
        };
        // Identically-zero derivatives (e.g. alpha0 delta/tau crosses)
        // compare exactly; everything else relatively.
        if rec.expected == 0.0 {
            if actual != 0.0 {
                failures.push(format!(
                    "{}(T={}, rhomolar={}): actual {actual:e}, expected exact 0",
                    rec.out, rec.t, rec.rhomolar
                ));
            }
            continue;
        }
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if !(rel <= 1e-13) {
            failures.push(format!(
                "{}(T={}, rhomolar={}): actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.out, rec.t, rec.rhomolar, rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} term records failed:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}
