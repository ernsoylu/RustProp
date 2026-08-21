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
use rustprop_golden_tests::{GoldenRecord, check, load_jsonl};
use std::path::Path;

fn fixture() -> Vec<GoldenRecord> {
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

/// The input side of one fixture record, for matching against the allowlist.
struct State {
    backend: &'static str,
    fluid: &'static str,
    name1: &'static str,
    val1: f64,
    name2: &'static str,
    val2: f64,
}

/// PLATFORM-SENSITIVE STATES — asserted on Linux/x86-64 ONLY.
///
/// Every fixture in this repo records what a LINUX/x86-64 CoolProp 8.0.0
/// wheel answered. Nothing here can say what CoolProp built against Apple's
/// or MSVC's libm would have answered, because the oracle was never run
/// there. So where the last bits of a libm call decide refusal-versus-answer,
/// asserting that classification off Linux asserts something the fixture
/// cannot support. This list is per-STATE and enumerated by hand; everything
/// else in this 1,626-record fixture — the HEOS alpha0 gates at T=1e30, the
/// melting-line refusals, the negative-pressure `post_update` arm, the cubic
/// states at 300 K and 500 K — keeps asserting on every platform, and the
/// skip counts below are ASSERTED so the set cannot grow unnoticed.
///
/// WHAT MAKES THESE FOUR DIFFERENT. All four are a cubic PT flash at a
/// temperature twenty-plus orders of magnitude outside any physical range,
/// and there the density that comes back is not a converged root at all — it
/// is the rounding residue of a cancellation. `solve_cubic`'s one-real-root
/// branch (`rustprop-heos/src/solvers.rs`) returns `x = t0 - b/(3a)` with
/// `t0 = -2*sqrt(p/3)*sinh(asinh(z)/3)`; for PR::Propane at T = 1e30 K both
/// terms are -713.06352787061..., while the physical root p/(R*T) = 1.2e-26
/// sits THIRTEEN orders of magnitude below ulp(713.06) =
/// 1.1368683772161603e-13. What survives the subtraction is therefore an
/// integer number of ulps, and which integer is decided by the last bits of
/// `sinh`/`asinh` — libm calls, not IEEE-mandated operations. Every other
/// number at the state follows deterministically from that density
/// (`rho_reducing = 1` for the cubics, so `delta` IS the density), which is
/// why the whole state is skipped rather than a hand-picked subset of its 33
/// records.
const PLATFORM_SENSITIVE: &[State] = &[
    // Root = +4 ulp: Dmolar 4.547473508864641e-13 on Linux/x86-64 (wheel and
    // port agree bitwise there), +2 ulp = 2.2737367544323206e-13 on
    // Windows-x64 (run 32468378222). The ideal-gas delta-derivatives inherit
    // the factor exactly, which is the tell:
    //   dalpha0_ddelta   =  1/delta    2.199023255552e12     vs 4.398046511104e12      (2x)
    //   d2alpha0_ddelta2 = -1/delta^2 -4.835703278458517e24  vs -1.9342813113834067e25 (4x)
    //   d3alpha0_ddelta3 =  2/delta^3  2.1267647932558654e37 vs 1.7014118346046923e38  (8x)
    State {
        backend: "PR",
        fluid: "Propane",
        name1: "T",
        val1: 1e30,
        name2: "P",
        val2: 101325.0,
    },
    // Root = 0 ulp on Linux/x86-64: Dmolar is EXACTLY zero, so `delta` is
    // zero, so `calc_alpha0_deriv_nocache` returns +-inf for every
    // nDelta >= 1 and upstream's ValidNumber gate throws — 22 of this state's
    // 33 records are refusals for that reason. Windows-x64 lands one ulp away
    // (1.1368683772161603e-13), the gate never fires, and the port ANSWERS
    // (8.796093022208e12 for dalpha0_ddelta, 1e0 for PIP). One ulp of a
    // -713.0634277389405 cancellation decides the whole state's
    // classification.
    State {
        backend: "PR",
        fluid: "Propane",
        name1: "T",
        val1: 1e20,
        name2: "P",
        val2: 101325.0,
    },
    // Root = +2 ulp of ulp(-3211.738769906704) = 4.547473508864641e-13, i.e.
    // Dmolar 9.094947017729282e-13. Windows-x64 AGREED here (run 32468378222
    // reported PR::Propane and nothing else); listed anyway because the
    // mechanism is the same one, measured: two ulps of a cancellation, not a
    // root.
    State {
        backend: "SRK",
        fluid: "Propane",
        name1: "T",
        val1: 1e30,
        name2: "P",
        val2: 101325.0,
    },
    // Root = -2 ulp: Dmolar -9.094947017729282e-13, a NEGATIVE density, so
    // `post_update` refuses "rhomolar is less than zero" for 31 of this
    // state's 33 records. The SIGN of the residue — and with it the entire
    // state's refuse-versus-answer classification — is two ulps of a
    // -3211.7387504776398 cancellation. Windows-x64 agreed; the margin is two
    // ulps.
    State {
        backend: "SRK",
        fluid: "Propane",
        name1: "T",
        val1: 1e20,
        name2: "P",
        val2: 101325.0,
    },
];

/// True when `rec` sits on an allowlisted state AND this is not the platform
/// the fixture was generated on. Always false on Linux/x86-64, where every
/// record keeps asserting exactly as it did before.
fn platform_sensitive(rec: &GoldenRecord) -> bool {
    !cfg!(all(target_os = "linux", target_arch = "x86_64"))
        && PLATFORM_SENSITIVE.iter().any(|s| {
            s.backend == rec.backend
                && s.fluid == rec.fluid
                && s.name1 == rec.name1
                && s.val1 == rec.val1
                && s.name2 == rec.name2
                && s.val2 == rec.val2
        })
}

/// How many records the allowlist removes on THIS platform: none on
/// Linux/x86-64, the enumerated count anywhere else. Every caller asserts its
/// number, so a skip set that grows — a new knife edge, or a stale entry —
/// turns the suite red instead of quietly shrinking its coverage.
fn expected_skips(off_linux: usize) -> usize {
    if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
        0
    } else {
        off_linux
    }
}

#[test]
fn validity_refusal_parity_matches_upstream() {
    let records = fixture();
    assert_eq!(records.len(), 1626);
    let refusals = records.iter().filter(|r| r.error.is_some()).count();
    assert_eq!(refusals, 450, "oracle refusals pinned");

    let mut failures = Vec::new();
    let mut skipped = 0usize;
    for rec in &records {
        if platform_sensitive(rec) {
            skipped += 1;
            continue;
        }
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
    println!(
        "platform-sensitive records skipped: {skipped} of {}",
        records.len()
    );
    assert_eq!(
        skipped,
        expected_skips(132),
        "the platform-sensitive skip set is exactly the four enumerated states"
    );
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
    let mut skipped = 0usize;
    let mut failures = Vec::new();
    for rec in &records {
        let Some(oracle) = rec.error.as_deref().map(oracle_message) else {
            continue;
        };
        if !oracle.starts_with("calc_alpha0_deriv_nocache") {
            continue;
        }
        if platform_sensitive(rec) {
            skipped += 1;
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
    println!("platform-sensitive alpha0-gated refusals skipped: {skipped}");
    assert_eq!(checked + skipped, 264, "alpha0-gated refusals pinned");
    assert_eq!(
        skipped,
        expected_skips(48),
        "the platform-sensitive skip set is exactly the four enumerated states"
    );
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
    let mut skipped = 0usize;
    for rec in empties {
        if platform_sensitive(rec) {
            skipped += 1;
            continue;
        }
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
    println!("platform-sensitive empty-message refusals skipped: {skipped}");
    assert_eq!(
        skipped,
        expected_skips(6),
        "the platform-sensitive skip set is exactly the four enumerated states"
    );
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
