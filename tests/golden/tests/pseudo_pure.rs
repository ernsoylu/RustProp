//! Pseudo-pure flash goldens (tier-2 deferral, final item): PT
//! (liquid/gas/supercritical), QT at the only defined qualities (0/1), PQ
//! across the temperature glide, and the caloric/density pairs
//! (H,P)/(P,S)/(P,U)/(D,P) over six regions — subcooled liquid, superheated
//! gas, supercritical, mid-glide two-phase, the pcrit..max_sat_p band and
//! sub-triple-pressure gas — through `props_si` for all six pseudo-pure
//! fluids; plus upstream's error conditions and its OWN failure regions — the
//! Air 1-bar low-quality (Hmolar,P) window and the pcrit..max_sat_p band
//! refusals — pinned as error-parity records from both sides. The ancillary
//! curves are exact evaluations and the density solves run at upstream's 1e-8
//! residuals —
//! 1e-8 policy with caloric outputs measured against the thermal scale where
//! they cross zero.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn pseudo_pure_flashes_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/pseudo_pure.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 665);

    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    let mut failures = Vec::new();
    let mut fluids = std::collections::HashSet::new();
    let mut value_records = 0usize;
    for rec in &records {
        fluids.insert(rec.fluid.clone());
        if rec.error.is_some() {
            continue; // pseudo_pure_error_parity_matches_upstream covers these
        }
        value_records += 1;
        let actual = match props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("HEOS::{}", rec.fluid),
        ) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{} {}: error {e}", rec.fluid, rec.id()));
                continue;
            }
        };
        // Caloric outputs cross zero at the arbitrary reference state where
        // pure relative error is ill-posed; measure them against the
        // thermal scale (R*Tc for H/U, R for S) when |expected| falls
        // below it (the established policy).
        let fluid = registry[rec.fluid.as_str()];
        let rg = fluid.eos.gas_constant;
        let tc = fluid.states.critical.t;
        let denom = match rec.out.as_str() {
            "Hmolar" | "Umolar" => rec.expected.abs().max(rg * tc),
            "Smolar" => rec.expected.abs().max(rg),
            _ => rec.expected.abs(),
        };
        let rel = (actual - rec.expected).abs() / denom;
        if rel > 1e-8 || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert_eq!(fluids.len(), 6, "all pseudo-pure fluids covered");
    assert_eq!(value_records, 654, "value records (665 - 11 error-parity)");
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        value_records,
        failures.join("\n")
    );
}

/// ERROR PARITY: the port must refuse exactly the states the oracle refuses and
/// answer exactly the ones it answers. Two upstream failure regions are pinned
/// in the fixture, both of them found by the (H,P)/(P,S)/(P,U)/(D,P) sweep:
///
/// 1. Air at 1 bar, low-quality (Hmolar,P). Below hmolar ~ 206.32 J/mol the
///    rational-polynomial caloric band classifies these two-phase states as
///    single-phase, and the bracketed PY solve then walks out of
///    `[Tmin, SatL->T()]` and throws (`FlashRoutines.cpp:3435`); at 206.5 J/mol
///    and above the same inputs converge on the two-phase branch. Pinned from
///    BOTH sides so the boundary cannot drift in either direction.
/// 2. The pcrit..max_sat_p band, where pL/pV must be inverted above their value
///    at Tmax. Air throws on (Hmolar,P) at the very state whose (P,S)/(P,U)/
///    (D,P) it answers; SES36 throws upstream's `post_update` refusal
///    ("T is not a valid number") on (P,Smolar); R410A cannot even reach the
///    forward PQ state (`solver_rho_Tp` fails, with the same T and density
///    guess the port computes).
///
/// The oracle's verbatim message rides along in each record's `error` field.
/// DIVERGENCE (documented, not fixed): only the `post_update` and
/// `solver_rho_Tp` texts match verbatim; the PY-flash refusals surface as the
/// port's own bracket diagnostic rather than upstream's "unable to solve 1phase
/// PY flash with Tmin=..., Tmax=... due to error: ..." wrapper, because that
/// wrapper quotes `HSU_P_flash_singlephase_Brent`'s guards and the established
/// 30-bit bisection stand-in has none. Refusal vs answer is asserted; the text
/// is not.
#[test]
fn pseudo_pure_error_parity_matches_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/pseudo_pure.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(
        records.iter().filter(|r| r.error.is_some()).count(),
        11,
        "oracle refusals pinned"
    );
    // The Air 1-bar (Hmolar,P) window, both sides of hmolar ~ 206.32 J/mol.
    let window: Vec<_> = records
        .iter()
        .filter(|r| r.fluid == "Air" && r.name1 == "Hmolar" && r.val2 == 1e5)
        .collect();
    assert_eq!(window.len(), 19, "the whole window, both sides");
    assert_eq!(
        window.iter().filter(|r| r.error.is_some()).count(),
        5,
        "refusals inside the window"
    );

    let mut failures = Vec::new();
    for rec in &records {
        let got = props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("HEOS::{}", rec.fluid),
        );
        match (&rec.error, &got) {
            (Some(_), Err(_)) | (None, Ok(_)) => {}
            (Some(msg), Ok(v)) => failures.push(format!(
                "{} {}: port answered {v:e} where the oracle refused ({msg})",
                rec.fluid,
                rec.id()
            )),
            (None, Err(e)) => failures.push(format!(
                "{} {}: port refused ({e}) where the oracle answered",
                rec.fluid,
                rec.id()
            )),
        }
    }
    assert!(
        failures.is_empty(),
        "{} parity failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Upstream's pseudo-pure error conditions.
#[test]
fn pseudo_pure_error_conditions() {
    use rustprop::Error;
    // QT with fractional quality: "Two-phase quality is not defined".
    match props_si("D", "T", 100.0, "Q", 0.5, "HEOS::Air").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("quality must be equal to 0 or 1"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // PT inside the ancillary band: two-phase inputs unsupported.
    let psat = props_si("P", "T", 100.0, "Q", 0.0, "HEOS::Air").unwrap();
    match props_si("D", "T", 100.0, "P", psat * 0.999, "HEOS::Air").unwrap_err() {
        Error::Value(msg) => assert!(
            msg.contains("Two-phase inputs not supported for pseudo-pure"),
            "unexpected message: {msg}"
        ),
        other => panic!("expected Value error, got {other:?}"),
    }
    // Unported pairs stay loud ((D,T) still routes through the guard;
    // (H,P)/(P,S)/(P,U)/(D,P) are served since the classic-ancillary
    // flash port).
    assert!(matches!(
        props_si("P", "Dmolar", 5000.0, "T", 100.0, "HEOS::Air").unwrap_err(),
        Error::NotImplemented(_)
    ));
}
