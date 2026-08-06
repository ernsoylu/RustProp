//! Term-level HEOS goldens (PLAN.md 4.1/4.7): alphar/alpha0 and every
//! tau/delta derivative the wheel exposes must match the port at
//! rel <= 1e-13 on a single-phase (T, rhomolar) grid, for every ported
//! fluid. The non-water fluids exercise the 4.7 term families: ideal-gas
//! Power (Nitrogen, R134a), PlanckEinsteinFunctionT (Nitrogen),
//! EnthalpyEntropyOffset (CarbonDioxide), and residual GaoB (Ammonia).

use rustprop_golden_tests::heos_fluids;
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
fn helmholtz_terms_match_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_terms.jsonl"));
        let text = std::fs::read_to_string(&path).unwrap();
        let records: Vec<TermRecord> = text
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        // Water: 12 points; parameterized fluids: 13 points x 20 accessors.
        let expected_count = if stem == "water" { 240 } else { 260 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let eos = HelmholtzEos::new(fluid);
        for rec in &records {
            assert_eq!(rec.fluid, name);
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
                        "{name} {}(T={}, rhomolar={}): actual {actual:e}, expected exact 0",
                        rec.out, rec.t, rec.rhomolar
                    ));
                }
                continue;
            }
            // Water's grid holds 1e-13. The reduced-coordinate grids of the
            // other fluids hit two cancellation-limited points (small
            // derivative values from large cancelling term contributions)
            // where BOTH the wheel and the port sit ~1e-14 absolute from the
            // 50-digit true value — see the PLAN.md Decisions log — so those
            // suites run at 1e-12.
            let rtol = if stem == "water" { 1e-13 } else { 1e-12 };
            let rel = ((actual - rec.expected) / rec.expected).abs();
            if rel > rtol || rel.is_nan() {
                failures.push(format!(
                    "{name} {}(T={}, rhomolar={}): actual {actual:e}, expected {:e}, rel {rel:e}",
                    rec.out, rec.t, rec.rhomolar, rec.expected
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} term records failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
