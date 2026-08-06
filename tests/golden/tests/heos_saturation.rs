//! Saturation goldens (PLAN.md 4.4): QT/PQ flashes at Q=0/1 across the dome
//! must match the wheel at the 1e-8 policy (superancillary path — observed
//! deviations are far smaller).

use rustprop_data::fluids::water::WATER;
use rustprop_golden_tests::load_jsonl;
use rustprop_heos::HelmholtzEos;
use rustprop_heos::saturation::{SatState, SaturationSuperAncillary};
use std::path::Path;

fn output(eos: &HelmholtzEos, sat: &SatState, out: &str) -> f64 {
    // Saturated-phase density at the requested quality end
    let rho = if sat.q == 0.0 { sat.rho_l } else { sat.rho_v };
    match out {
        "P" => sat.p,
        "T" => sat.t,
        "Dmolar" => rho,
        "Hmolar" => eos.hmolar(sat.t, rho),
        "Smolar" => eos.smolar(sat.t, rho),
        "Umolar" => eos.umolar(sat.t, rho),
        "Cpmolar" => eos.cpmolar(sat.t, rho),
        "A" => eos.speed_sound(sat.t, rho),
        other => panic!("unknown output {other}"),
    }
}

#[test]
fn water_saturation_states_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/heos_water_sat.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 292);

    let eos = HelmholtzEos::new(&WATER);
    let sa = SaturationSuperAncillary::new(WATER.eos.superancillary.as_ref().unwrap());

    let mut failures = Vec::new();
    let mut worst: (f64, String) = (0.0, String::new());
    for rec in &records {
        let sat = match (rec.name1.as_str(), rec.name2.as_str()) {
            ("T", "Q") => sa.qt_flash(rec.val1, rec.val2).unwrap(),
            ("P", "Q") => sa.pq_flash(rec.val1, rec.val2).unwrap(),
            other => panic!("unexpected pair {other:?}"),
        };
        let actual = output(&eos, &sat, &rec.out);
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > worst.0 {
            worst = (rel, rec.id());
        }
        if rel > 1e-8 || rel.is_nan() {
            failures.push(format!(
                "{}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.id(),
                rec.expected
            ));
        }
    }
    println!("worst deviation: {:e} at {}", worst.0, worst.1);
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Above-critical guards mirror upstream's errors.
#[test]
fn saturation_guards() {
    let sa = SaturationSuperAncillary::new(WATER.eos.superancillary.as_ref().unwrap());
    assert!(sa.qt_flash(650.0, 0.0).is_err());
    assert!(sa.pq_flash(23e6, 0.0).is_err());
    // Mixture density at intermediate quality follows the inverse mixing rule.
    let sat = sa.qt_flash(400.0, 0.5).unwrap();
    let expected = 1.0 / (0.5 / sat.rho_v + 0.5 / sat.rho_l);
    assert_eq!(sat.rhomolar, expected);
}
