//! Saturation goldens (PLAN.md 4.4/4.7): QT/PQ flashes at Q=0/1 across the
//! dome must match the wheel at the 1e-8 policy (superancillary path —
//! observed deviations are far smaller), for every ported fluid.

use rustprop_golden_tests::{heos_fluids, load_jsonl};
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
fn saturation_states_match_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_sat.jsonl"));
        let records = load_jsonl(&path);
        let expected_count = if stem == "water" { 292 } else { 268 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let eos = HelmholtzEos::new(fluid);
        let sa = SaturationSuperAncillary::new(fluid.eos.superancillary.as_ref().unwrap());

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
                    "{name} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                    rec.id(),
                    rec.expected
                ));
            }
        }
        println!("{name}: worst deviation {:e} at {}", worst.0, worst.1);
    }
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
    let water = heos_fluids()[0].2;
    let sa = SaturationSuperAncillary::new(water.eos.superancillary.as_ref().unwrap());
    assert!(sa.qt_flash(650.0, 0.0).is_err());
    assert!(sa.pq_flash(23e6, 0.0).is_err());
    // Mixture density at intermediate quality follows the inverse mixing rule.
    // Bitwise is legitimate: this is an algebraic identity, not an oracle
    // comparison. `qt_flash` evaluates this same expression from the same two
    // doubles it hands back, so anything but bit-equality would mean the
    // mixing rule had changed shape.
    let sat = sa.qt_flash(400.0, 0.5).unwrap();
    let expected = 1.0 / (0.5 / sat.rho_v + 0.5 / sat.rho_l);
    assert_eq!(sat.rhomolar, expected);
}
