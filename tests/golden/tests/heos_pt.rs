//! PT-flash goldens (PLAN.md 4.5/4.7): the single-phase density solver and
//! the properties through it must match `PT_INPUTS` states of the wheel at
//! the 1e-9 policy across liquid, vapor, and supercritical states, for
//! every ported fluid.

use rustprop_golden_tests::{heos_fluids, load_jsonl};
use rustprop_heos::flash_pt::PtFlash;
use std::path::Path;

#[test]
fn pt_states_match_upstream() {
    let mut failures = Vec::new();
    for (name, stem, fluid) in heos_fluids() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("fixtures/heos_{stem}_pt.jsonl"));
        let records = load_jsonl(&path);
        let expected_count = if stem == "water" { 100 } else { 90 };
        assert_eq!(records.len(), expected_count, "{name}: fixture size");

        let flash = PtFlash::new(fluid);
        let mut worst: (f64, String) = (0.0, String::new());
        for rec in &records {
            assert_eq!((rec.name1.as_str(), rec.name2.as_str()), ("T", "P"));
            let (t, p) = (rec.val1, rec.val2);
            let (rho, _phase) = flash
                .pt_flash(t, p)
                .unwrap_or_else(|e| panic!("{name} {}: {e}", rec.id()));
            let actual = match rec.out.as_str() {
                "Dmolar" => rho,
                "Hmolar" => flash.eos.hmolar(t, rho),
                "Smolar" => flash.eos.smolar(t, rho),
                "Cpmolar" => flash.eos.cpmolar(t, rho),
                "A" => flash.eos.speed_sound(t, rho),
                other => panic!("unknown output {other}"),
            };
            // rho/h/s hold 1e-9 everywhere. Cp and speed of sound are
            // second-derivative properties: near-critical they amplify the
            // density solver's upstream stopping residual (ftol 1e-8 —
            // matched deliberately, see PLAN.md), so they run at the
            // solver-dependent 1e-8 tier.
            let rtol = match rec.out.as_str() {
                "Cpmolar" | "A" => 1e-8,
                _ => 1e-9,
            };
            let rel = ((actual - rec.expected) / rec.expected).abs();
            if rel > worst.0 {
                worst = (rel, rec.id());
            }
            if rel > rtol || rel.is_nan() {
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

/// Condition parity: near-saturation PT input errors like upstream, and the
/// critical-point short-circuit engages.
#[test]
fn pt_guards() {
    let water = heos_fluids()[0].2;
    let flash = PtFlash::new(water);
    // psat(400 K) ~ 245769 Pa: within 1e-6 relative must error
    let psat = flash.sat().qt_flash(400.0, 1.0).unwrap().p;
    assert!(flash.pt_flash(400.0, psat * (1.0 + 1e-8)).is_err());
    // Critical short-circuit
    let (rho, phase) = flash.pt_flash(647.096, 22064000.0).unwrap();
    assert_eq!(rho, water.states.critical.rhomolar);
    assert_eq!(phase, rustprop_core::params::Phase::CriticalPoint);
}
