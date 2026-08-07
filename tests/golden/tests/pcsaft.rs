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

#[test]
fn pcsaft_flash_matches_oracle() {
    // Slice 11c: flashes through the string API. 1e-7 policy: the
    // inside-out loops TARGET maxdif 1e-8 per side (acceptance far looser),
    // so cross-implementation agreement is bounded by ~the loop tolerance —
    // observed up to 1.4e-8 on low-T ACETONE/METHANOL saturation.
    // Near-vacuum states are excluded by the generator's p_sat floor.
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/pcsaft_flash.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 164);
    let mut failures = Vec::new();
    for rec in &recs {
        let fluid = format!("PCSAFT::{}", rec.fluid);
        match rustprop::props_si(&rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &fluid) {
            Ok(actual) => {
                if let Err(e) = rustprop_golden_tests::check(rec, actual, 1e-7) {
                    failures.push(e);
                }
            }
            Err(e) => failures.push(format!("{}: error {e:?}", rec.id())),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} PC-SAFT flash records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn pcsaft_error_conditions() {
    use rustprop::props_si;
    use rustprop_core::Error;
    let err = |r: rustprop_core::Result<f64>| r.unwrap_err();

    // Unsupported input pairs carry upstream's message
    match err(props_si("T", "Hmolar", 1e4, "P", 1e5, "PCSAFT::PROPANE")) {
        Error::Value(m) => assert!(m.contains("is not yet supported"), "got {m}"),
        other => panic!("wrong variant {other:?}"),
    }
    // Absolute calorics are base-class NotImplemented upstream
    match err(props_si("Hmolar", "T", 300.0, "P", 1e5, "PCSAFT::PROPANE")) {
        Error::NotImplemented(_) => {}
        other => panic!("wrong variant {other:?}"),
    }
    // Q range
    match err(props_si("P", "T", 300.0, "Q", 1.5, "PCSAFT::PROPANE")) {
        Error::OutOfRange(m) => {
            assert_eq!(m, "Input vapor quality [Q] must be between 0 and 1")
        }
        other => panic!("wrong variant {other:?}"),
    }
    // Unknown fluid key carries the PCSAFT library message
    match err(props_si("P", "T", 300.0, "Q", 0.0, "PCSAFT::NotAFluid")) {
        Error::Value(m) => assert!(
            m.contains("was not found in string_to_index_map in PCSAFTLibraryClass"),
            "got {m}"
        ),
        other => panic!("wrong variant {other:?}"),
    }
    // Missing kij pair fails at construction
    match err(props_si(
        "P",
        "T",
        300.0,
        "Q",
        0.0,
        "PCSAFT::PROPANE[0.5]&WATER[0.5]",
    )) {
        Error::Value(m) => assert!(m.contains("Could not match the binary pair"), "got {m}"),
        other => panic!("wrong variant {other:?}"),
    }
}
