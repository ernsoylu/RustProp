//! Validity parity over ABUSIVE inputs (R11) — states far outside the
//! fluid's range, the dimension the seeded acceptance sweep deliberately does
//! not draw.
//!
//! The question this suite asks is not "is the number right" but "does the
//! port refuse exactly where the wheel refuses". Two upstream gates decide
//! that, and both are pinned here:
//!
//! * `HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache`'s closing
//!   `ValidNumber` (`HelmholtzEOSMixtureBackend.cpp:3621`) — a real C++
//!   throw whose message names the `nTau`/`nDelta` upstream asked for, so the
//!   text identifies WHICH ideal-gas derivative the calculator reached for.
//!   `alpha0_messages_name_the_same_derivative` asserts it verbatim.
//! * `_raise_if_invalid` (`nanobind_interface.cxx:104`), the scalar
//!   PropsSI/HAPropsSI boundary check: a non-finite result is a `ValueError`
//!   carrying the global error string, which is EMPTY when nothing threw
//!   underneath. `Smolar` lands there — `calc_smolar` reads the UNCHECKED
//!   bulk `calc_all_alpha0_derivs_nocache` — and so do `L` and `V`.
//!
//! Before R11 the port answered `Ok(NaN)`/`Ok(inf)` for the second class and
//! a plausible-looking finite speed of sound for part of the first.
//!
//! The `SRK::`/`PR::` records carry a third gate the cubic route was missing
//! outright: `AbstractCubicBackend::update` closes with the inherited
//! `post_update()` (`CubicBackend.cpp:387`), so `PropsSI("Dmolar","T",1000,
//! "P",-1,"SRK::Propane")` is "rhomolar is less than zero", not the negative
//! density the port used to serve.

use rustprop::props_si;
use rustprop_golden_tests::{check, load_jsonl};
use std::path::Path;

fn fixture() -> Vec<rustprop_golden_tests::GoldenRecord> {
    load_jsonl(&Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/validity.jsonl"))
}

/// The wheel's message carries `PropsSI`'s own ` : PropsSI("out",...)` tail
/// (appended by the C++ `PropsSI` catch block); the port has no such
/// wrapper, so comparisons use the message proper.
fn oracle_message(err: &str) -> &str {
    match err.find(" : PropsSI(") {
        Some(i) => &err[..i],
        None => err,
    }
}

#[test]
fn validity_refusal_parity_matches_upstream() {
    let records = fixture();
    assert_eq!(records.len(), 1626);
    let refusals = records.iter().filter(|r| r.error.is_some()).count();
    assert_eq!(refusals, 450, "oracle refusals pinned");

    let mut failures = Vec::new();
    for rec in &records {
        let got = props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("{}::{}", rec.backend, rec.fluid),
        );
        match (&rec.error, &got) {
            (Some(_), Err(_)) => {}
            // 1e-8, not the usual 1e-9: the states here are chosen to be
            // absurd, and four of the answers are dominated by
            // cancellation rather than by the EOS (Water `Z` at 300 K is
            // `1 + delta*dalphar_ddelta` with the bracket ~ -0.99927, so a
            // 1.4e-12 absolute wobble reads as 1.9e-9 relative; `alphar` and
            // `Gmolar_residual` at 1e10 K are the same shape). The
            // classification, not the digit count, is what this suite pins.
            (None, Ok(v)) => {
                if let Err(e) = check(rec, *v, 1e-8) {
                    failures.push(format!("{}: {e}", rec.fluid));
                }
            }
            (Some(msg), Ok(v)) => failures.push(format!(
                "{} {}: port answered {v:e} where the oracle refused ({})",
                rec.fluid,
                rec.id(),
                oracle_message(msg)
            )),
            (None, Err(e)) => failures.push(format!(
                "{} {}: port refused ({e}) where the oracle answered {:e}",
                rec.fluid,
                rec.id(),
                rec.expected
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

/// The alpha0 gate is per-derivative upstream, and the message proves it:
/// `Hmolar` reports `nTau: 1` (`dalpha0_dTau()`), `Cpmolar`/`Cvmolar`/`A`
/// report `nTau: 2` (`d2alpha0_dTau2()`), `Gmolar`/`Helmholtzmolar` report
/// `nTau: 0` (`alpha0()`). A port that gated on "any non-finite alpha0
/// entry" would pass the refusal test above and fail this one.
#[test]
fn alpha0_messages_name_the_same_derivative() {
    let records = fixture();
    let mut checked = 0;
    let mut failures = Vec::new();
    for rec in &records {
        let Some(oracle) = rec.error.as_deref().map(oracle_message) else {
            continue;
        };
        if !oracle.starts_with("calc_alpha0_deriv_nocache") {
            continue;
        }
        checked += 1;
        match props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("{}::{}", rec.backend, rec.fluid),
        ) {
            Err(e) if e.to_string() == oracle => {}
            other => failures.push(format!(
                "{} {}: got {other:?}\n   want {oracle:?}",
                rec.fluid,
                rec.id()
            )),
        }
    }
    assert_eq!(checked, 264, "alpha0-gated refusals pinned");
    assert!(
        failures.is_empty(),
        "{} message mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The boundary-gate half: where nothing threw underneath, upstream's error
/// string is empty and the wheel's `ValueError` carries "". `Smolar` is the
/// canonical case — `calc_smolar` reads the bulk alpha0 path, which has no
/// gate — and the port must refuse with an empty message rather than serve
/// the `inf` the C++ layer computes.
#[test]
fn boundary_gate_refusals_carry_an_empty_message() {
    let records = fixture();
    let empties: Vec<_> = records
        .iter()
        .filter(|r| r.error.as_deref() == Some(""))
        .collect();
    assert_eq!(empties.len(), 25, "empty-message refusals pinned");
    assert!(
        empties.iter().any(|r| r.out == "Smolar"),
        "the Smolar witness is in the fixture"
    );
    for rec in empties {
        let err = props_si(
            &rec.out,
            &rec.name1,
            rec.val1,
            &rec.name2,
            rec.val2,
            &format!("{}::{}", rec.backend, rec.fluid),
        )
        .expect_err("the oracle refuses this state");
        assert_eq!(err.to_string(), "", "{} {}", rec.fluid, rec.id());
    }
}

/// The witness from the R11 handoff, spelled out: one call, three outcomes
/// that used to be `Ok(non-finite)` and one that legitimately answers.
#[test]
fn handoff_witness_is_closed() {
    use rustprop::Error;
    // The handoff's own witness. Water's IAPWS-2011 conductivity reaches
    // `d2alpha0_dTau2()` through cp/cv in the critical term, so the refusal
    // carries upstream's text rather than the boundary gate's empty string.
    match props_si("L", "T", 1e30, "P", 101325.0, "Water").unwrap_err() {
        Error::Value(m) => assert_eq!(
            m,
            "calc_alpha0_deriv_nocache returned invalid number with inputs \
             nTau: 2, nDelta: 0, tau: 6.47096e-28, delta: 6.81824e-31"
        ),
        other => panic!("expected a ValueError, got {other:?}"),
    }
    // Enthalpy names the FIRST tau derivative instead — the gate is
    // per-derivative, exactly as upstream's is.
    match props_si("Hmolar", "T", 1e30, "P", 101325.0, "Water").unwrap_err() {
        Error::Value(m) => assert_eq!(
            m,
            "calc_alpha0_deriv_nocache returned invalid number with inputs \
             nTau: 1, nDelta: 0, tau: 6.47096e-28, delta: 6.81824e-31"
        ),
        other => panic!("expected a ValueError, got {other:?}"),
    }
    // Entropy has no gate to hit: `calc_smolar` reads the unchecked bulk
    // alpha0 path, so the wheel refuses through `_raise_if_invalid` with an
    // empty message and so does the port.
    match props_si("Smolar", "T", 1e30, "P", 101325.0, "Water").unwrap_err() {
        Error::Value(m) => assert_eq!(m, ""),
        other => panic!("expected an empty-message ValueError, got {other:?}"),
    }
    // Speed of sound used to ANSWER here: the non-finite d2alpha0_dTau2 sits
    // in a denominator, so the arithmetic alone produced a plausible number.
    assert!(props_si("A", "T", 1e30, "P", 101325.0, "Water").is_err());
    // Density still answers, exactly as the wheel does.
    let d = props_si("Dmolar", "T", 1e30, "P", 101325.0, "Water").unwrap();
    assert!((d / 1.218673013775591e-26 - 1.0).abs() < 1e-12, "{d:e}");
}

/// The humid-air half. Upstream's scalar `HAPropsSI` binding closes with the
/// same `_raise_if_invalid` (`nanobind_interface.cxx:1104`), and the R11
/// abuse scan (21,424 `HAPropsSI` calls over 4 input triples x abusive
/// ladders) found 800 where the port answered `Ok(inf)`/`Ok(NaN)` and the
/// wheel refused with an empty message. EVERY one of the 800 has a
/// non-finite INPUT, which JSON Lines cannot carry, so this is a literal
/// test rather than a fixture; the wheel's side of each assertion below was
/// read off that scan. Zero rows ran the other way (port refuses, wheel
/// answers), so the gate cannot have over-refused.
#[test]
fn humid_air_boundary_gate_refuses_non_finite() {
    // The echo route is where it bites: the output IS one of the inputs, so
    // the value is handed straight back without a solve and nothing throws.
    for (out, v3) in [
        ("R", f64::INFINITY),
        ("R", f64::NAN),
        ("W", f64::INFINITY),
        ("W", f64::NAN),
        ("H", f64::NAN),
    ] {
        let err = rustprop::ha_props_si(out, "T", 300.0, "P", 101325.0, out, v3)
            .expect_err("the wheel refuses this state");
        assert_eq!(
            err.to_string(),
            "",
            "HAPropsSI({out},T,300,P,101325,{out},{v3})"
        );
    }
    // The control group: the same call with a physical third input answers,
    // and the humid-air suite pins its value.
    assert!(rustprop::ha_props_si("H", "T", 300.0, "P", 101325.0, "R", 0.5).is_ok());
}

/// DIVERGENCE, pinned rather than chased (the 2026-08-17 handoff's candidate
/// #2). `post_update`'s arms are carried in upstream's order and both fire on
/// the right condition, but at this state the port's `(Smass, T)` flash
/// leaves `rhomolar` NaN where upstream's leaves it NEGATIVE, so the two
/// implementations trip DIFFERENT arms of the same gate:
///
/// ```text
/// wheel:    "rhomolar is less than zero"
/// rustprop: "rhomolar is not a valid number"
/// ```
///
/// Refusal-vs-answer agrees, which is the fidelity property; matching the
/// text would mean reproducing a divergent iteration's garbage bit for bit,
/// on a path where neither implementation has an answer. Same family as the
/// R11 scan's residual `a2` rows (the port converging somewhere upstream
/// does not, at sub-triple or negative-temperature inputs).
#[test]
fn post_update_refusal_text_divergence_pinned() {
    use rustprop::Error;
    match props_si("Hmass", "T", 0.0, "Smass", 101325.0, "Water").unwrap_err() {
        Error::Value(m) => assert_eq!(m, "rhomolar is not a valid number"),
        other => panic!("expected a ValueError, got {other:?}"),
    }
}
