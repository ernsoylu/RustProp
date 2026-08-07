//! PC-SAFT goldens (PLAN.md Phase 11): kernel-level values against the
//! wheel's low-level accessors at imposed-phase (Dmolar, T) states.

use std::path::Path;

fn fluid(name: &str) -> &'static rustprop_core::fluid::PcsaftFluid {
    rustprop_data::pcsaft::PCSAFT_FLUIDS
        .iter()
        .find(|f| f.name == name)
        .unwrap_or_else(|| panic!("no PC-SAFT fluid {name}"))
}

#[test]
fn pcsaft_terms_match_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/pcsaft_terms.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 80);
    let mut failures = Vec::new();
    for rec in &recs {
        let names: Vec<&str> = rec.fluid.split('&').collect();
        let fluids: Vec<_> = names.iter().map(|n| fluid(n)).collect();
        let mut backend = rustprop_pcsaft::PcsaftBackend::new(
            &fluids,
            rustprop_data::pcsaft::PCSAFT_BINARY_PAIRS,
        )
        .expect("backend builds");
        if names.len() == 2 {
            backend.set_mole_fractions(&[rec.val3, 1.0 - rec.val3]);
        } else if names.len() == 3 {
            // Na+ / Cl- / WATER case: x1 = x2, x3 = 1 - 2 x1
            backend.set_mole_fractions(&[rec.val3, rec.val3, 1.0 - 2.0 * rec.val3]);
        }
        backend
            .set_state_dmolar_t(rec.val1, rec.val2)
            .expect("state set");
        let actual = match rec.out.as_str() {
            "P" => backend.calc_pressure(),
            "alphar" => backend.calc_alphar(),
            "Hmolar_residual" => backend.calc_hmolar_residual(),
            "Smolar_residual" => backend.calc_smolar_residual(),
            "Gmolar_residual" => backend.calc_gibbsmolar_residual(),
            other => panic!("unknown out {other}"),
        };
        // Direct evaluations: 1e-12 default; liquid-state P inflates Z's
        // cancellation (Z ~ 1e-3 from O(10) terms), so P carries a
        // per-suite guard at 1e-9 through the fixture rtol field when the
        // generator stamps it (none currently need it).
        if let Err(e) = rustprop_golden_tests::check(rec, actual, 1e-12) {
            failures.push(e);
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} PC-SAFT term records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}
