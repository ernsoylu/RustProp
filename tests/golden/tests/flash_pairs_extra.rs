//! Tier-2 input-pair goldens (PLAN.md 4.6 deferrals): `(Hmolar, T)` and
//! `(T, Umolar)` — upstream `DHSU_T_flash` —, `(P, Umolar)` — upstream
//! `HSU_P_flash` —, and `(Dmolar, Hmolar/Smolar/Umolar)` — upstream
//! `HSU_D_flash`'s superancillary happy path — through `props_si`,
//! including the mass-basis variants, and `(Dmolar, Q)` — upstream
//! `DQ_flash` (superancillary root resolution at the Q boundaries, Brent on
//! the density-implied quality for fractional Q). The T-based pairs solve
//! density directly (1e-9); the (P, X) pair runs at the established 1e-8
//! policy (upstream resolves T at ~30 bits).

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

#[test]
fn extra_flash_pairs_match_upstream() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/flash_pairs_extra.jsonl");
    let records = load_jsonl(&path);
    assert_eq!(records.len(), 693);

    let mut failures = Vec::new();
    for rec in &records {
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
        // (P, X): upstream resolves T at ~30 bits. P outputs of the T-based
        // caloric pairs: the density is Halley-solved to a 1e-8 residual and
        // the liquid's stiff dp/drho amplifies that into P (absolute error
        // stays ~1e-4 Pa) — the established solver-dependent 1e-8 tier.
        let rtol = if rec.name1 == "P" || rec.out == "P" {
            1e-8
        } else {
            1e-9
        };
        // Q = -1 sentinels compare exactly.
        if rec.expected == -1.0 && rec.out == "Q" {
            if actual != -1.0 {
                failures.push(format!(
                    "{} {}: got {actual}, expected -1",
                    rec.fluid,
                    rec.id()
                ));
            }
            continue;
        }
        let rel = ((actual - rec.expected) / rec.expected).abs();
        if rel > rtol || rel.is_nan() {
            failures.push(format!(
                "{} {}: actual {actual:e}, expected {:e}, rel {rel:e}",
                rec.fluid,
                rec.id(),
                rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} failures:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
}

/// Error-condition parity for the pair-level dead ends and the DQ guards.
#[test]
fn input_pair_error_conditions() {
    use rustprop::Error;
    // Upstream's own dead ends: mass_to_molar_inputs never converts these
    // ("This pair of inputs [X] is not yet supported").
    for (n1, v1, n2, v2) in [
        ("Hmass", 2e5, "T", 300.0),
        ("T", 300.0, "Umass", 1.1e5),
        ("Smolar", 80.0, "Umolar", 2e4),
    ] {
        match props_si("D", n1, v1, n2, v2, "HEOS::Water").unwrap_err() {
            Error::Value(msg) => assert!(
                msg.contains("is not yet supported"),
                "unexpected message: {msg}"
            ),
            other => panic!("expected Value error, got {other:?}"),
        }
    }
    // Combinations upstream's generate_update_pair does not know (H+Q, Q+S):
    // INPUT_PAIR_INVALID with non-trivial outputs.
    for (n1, v1, n2, v2) in [("Hmolar", 3e4, "Q", 1.0), ("Q", 0.0, "Smolar", 80.0)] {
        match props_si("D", n1, v1, n2, v2, "HEOS::Water").unwrap_err() {
            Error::Value(msg) => assert!(
                msg.contains("Input pair variable is invalid"),
                "unexpected message: {msg}"
            ),
            other => panic!("expected Value error, got {other:?}"),
        }
    }
    // DQ strict mode: water's liquid-density anomaly gives two T-roots.
    match props_si("T", "Dmolar", 55500.0, "Q", 0.0, "HEOS::Water").unwrap_err() {
        Error::Value(msg) => assert!(msg.contains("2 T-roots"), "unexpected message: {msg}"),
        other => panic!("expected Value error, got {other:?}"),
    }
    // DQ out of range: density between the branches has no root on either.
    match props_si("T", "Dmolar", 10000.0, "Q", 0.0, "HEOS::Water").unwrap_err() {
        Error::OutOfRange(msg) => assert!(msg.contains("no T-root"), "unexpected message: {msg}"),
        other => panic!("expected OutOfRange error, got {other:?}"),
    }
}
