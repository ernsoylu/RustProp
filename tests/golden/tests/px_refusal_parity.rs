//! (P, caloric) refusal-vs-answer parity (Wave-4 R12).
//!
//! A full-registry scan — all 136 HEOS fluids, every superancillary-derived
//! `(P, Hmass/Smass/Umass)` state over a pressure ladder from just above the
//! triple point to 0.999*pcrit and a 47-point quality ladder, 157,374 states —
//! classified each point answer-vs-refuse on both sides. The wheel refuses at
//! 297 of them, all in three fluids (Air 177, SES36 89, R407C 31), and the
//! port refuses at exactly the same 297 but three: one MethylLinoleate state,
//! drawn once per caloric key. This suite pins both halves of that result.
//!
//! **The Air band** (Wave-3 F8's open question). Upstream's `HSU_P_flash`
//! classifies a low-quality `(p, h)` state as GAS and hands
//! `HSU_P_flash_singlephase_Brent` the bracket [Tsat(p), 1.5*Tmax]. The
//! residual is negative at both ends, so `resid_Tmin * resid_Tmax < 0` fails
//! and TOMS748 never runs; upstream falls to `Halley` from whichever endpoint
//! has the smaller residual, which walks out of range, and then — only when
//! `p > 0.95*pcrit` — to a 2-D Newton that also fails. The port stops at the
//! bracket test, because neither continuation is ported.
//!
//! Bisecting the two classification boundaries on both sides, at 1, 2, 5, 7,
//! 20 bar and 0.8/0.95/0.99/0.999*pcrit: the UPPER edge (where the gas bracket
//! starts to contain the root) is BITWISE identical on both sides at all nine
//! pressures, and the lower edge — the one at the very bottom of the dome,
//! where the classification turns on solver noise — agrees to between 7.4e-13
//! and 7.1e-9 of the dome width in quality (worst at 0.95*pcrit). So the two
//! implementations refuse over the same band for the same reason: the unported
//! ladder never rescues a state, it only changes the message. The MESSAGES are
//! therefore NOT asserted (the PY-flash text divergence in NEXT-STEPS covers
//! that); the classification is.
//!
//! **The MethylLinoleate knife-edge**, pinned by
//! `methyl_linoleate_superanc_extrapolation_divergence_pinned` below.

use rustprop::props_si;
use rustprop_golden_tests::load_jsonl;
use std::path::Path;

fn fixture() -> Vec<rustprop_golden_tests::GoldenRecord> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/px_refusal_parity.jsonl");
    load_jsonl(&path)
}

/// The port answers exactly where the wheel answers and refuses exactly where
/// the wheel refuses — except the three pinned MethylLinoleate rows, which are
/// counted so the pin cannot silently heal or widen.
#[test]
fn px_refusal_parity_matches_upstream() {
    let records = fixture();
    assert_eq!(records.len(), 175);
    let oracle_refusals = records.iter().filter(|r| r.error.is_some()).count();
    assert_eq!(oracle_refusals, 47, "oracle refusals pinned");

    let mut failures = Vec::new();
    let mut pinned = 0usize;
    for rec in &records {
        let spec = format!("HEOS::{}", rec.fluid);
        let got = props_si(&rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &spec);
        match (&rec.error, &got) {
            // Both refuse. The text is not compared: upstream's wrapper quotes
            // whichever of Halley / 2-D Newton failed, neither of which is
            // ported (NEXT-STEPS' PY-flash message row).
            (Some(_), Err(_)) => {}
            (None, Ok(v)) => {
                // (P, X) flashes ride the established 1e-8 solver tier.
                if let Err(e) = rustprop_golden_tests::check(rec, *v, 1e-8) {
                    failures.push(e);
                }
            }
            (Some(msg), Ok(v)) => failures.push(format!(
                "{}: port answers {v:e} where the wheel refuses ({msg})",
                rec.id()
            )),
            (None, Err(e)) => {
                if rec.fluid == "MethylLinoleate" {
                    pinned += 1;
                    assert_eq!(
                        e.to_string(),
                        "rhomolar is less than zero",
                        "{}: the pinned divergence must refuse through upstream's own \
                         post_update arm, not some other path",
                        rec.id()
                    );
                } else {
                    failures.push(format!(
                        "{}: port refuses ({e}) where the wheel answers {:e}",
                        rec.id(),
                        rec.expected
                    ));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} records disagree:\n{}",
        failures.len(),
        records.len(),
        failures.join("\n")
    );
    // Heal detection: if the port ever answers these, delete the pin (and the
    // NEXT-STEPS row) rather than loosening this count.
    assert_eq!(
        pinned, 3,
        "the pinned MethylLinoleate divergence changed shape; see \
         methyl_linoleate_superanc_extrapolation_divergence_pinned"
    );
}

/// DIVERGENCE, pinned rather than chased: the only answer-vs-refuse
/// disagreement in a 157,374-state registry-wide scan.
///
/// The state is MethylLinoleate at p = 1.001 * ptriple = 1.312561985443192e-06
/// Pa, with the caloric input taken from the wheel's own Q = 0 value — i.e.
/// exactly on the saturated-liquid boundary. Three things stack up:
///
/// 1. This fluid's JSON `satminL` block says (T = 260 K, p = 1.3113e-06 Pa),
///    but its superancillary p-curve evaluates p(260 K) = 1.4986e-06 Pa. The
///    two disagree by 14%. Upstream's phase determination gates on the JSON
///    value (`components[0].EOS().ptriple`, itself read from `satminL/p`), so a
///    pressure in between is admitted and then handed to
///    `SuperAncillary::get_T_from_p`, which is `m_invlnp.eval(log(p))` — an
///    EXTRAPOLATION below the fitted domain. Both implementations extrapolate.
/// 2. `make_invlnp` fits to `tol = 1e-12`, and the fitted coefficients differ
///    in their last bits (upstream contracts `L * f` through Eigen; the port
///    sums sequentially). Inside the domain that is worth ~2e-15 relative on
///    T(p) for every fluid measured; extrapolated, here, it is 2.0e-12 — which
///    is 1.4e-12 relative on the saturated-liquid `u`, or 1.3e-6 J/kg.
/// 3. The quality then follows from `(u - u_L)/(u_V - u_L)`, so that 1.3e-6
///    J/kg decides its SIGN: the wheel gets q = 0, the port q = -3.4e-12.
///    Upstream's `q < -1e-9` test keeps both in the two-phase branch, and at
///    1.3 uPa the vapour density is 1.6e-10 mol/m^3 — twelve orders below the
///    liquid — so a quality of -3.4e-12 makes the lever-rule mixture density
///    NEGATIVE. Upstream's own `post_update` gate (ported) then refuses it.
///
/// So the port is running upstream's algorithm, arriving at upstream's guard,
/// and refusing with upstream's message. What differs is which side of q = 0
/// a Chebyshev extrapolation lands on. The control rows in the fixture prove
/// it: 10,000 ulp below the Q = 0 value (q = -3.0e-13) the WHEEL refuses with
/// this exact message, and 10,000 ulp above it both answer.
///
/// Not fixable without reproducing Eigen's reduction order bit for bit — the
/// same ruling as the SVDSBTL `-ffp-contract=fast` row.
#[test]
fn methyl_linoleate_superanc_extrapolation_divergence_pinned() {
    let fl = "HEOS::MethylLinoleate";
    let p = 1.312561985443192e-06;

    // The premise: this pressure is below the superancillary p-curve's own
    // minimum, so T(p) is an extrapolation.
    let p_at_tmin = props_si("P", "T", 260.0, "Q", 0.0, fl).unwrap();
    assert!(
        p < p_at_tmin,
        "the state must sit below p(Tmin) = {p_at_tmin:e} for the pin to describe it"
    );

    // The saturation flash itself is fine and agrees with the wheel to 1e-11.
    let t_sat = props_si("T", "P", p, "Q", 0.0, fl).unwrap();
    assert!(
        ((t_sat - 259.306893811023) / t_sat).abs() < 1e-11,
        "Tsat {t_sat}"
    );

    for key in ["Hmass", "Smass", "Umass"] {
        // The wheel's Q = 0 value, verbatim from the oracle.
        let v0 = match key {
            "Hmass" => -924116.4092398835,
            "Smass" => -2177.730005084435,
            _ => -924116.4092398802,
        };
        let err = props_si("T", "P", p, key, v0, fl)
            .expect_err("the pinned divergence: the wheel answers 259.306893811023 here");
        assert_eq!(err.to_string(), "rhomolar is less than zero", "{key}");
        // Two ulp up is enough headroom for the port's own extrapolated u_L:
        // the disagreement is 3.4e-12 in quality, not in the last bit.
        let up = props_si("T", "P", p, key, v0 * (1.0 - 1e-11), fl)
            .unwrap_or_else(|e| panic!("{key} just above the boundary: {e}"));
        assert!(((up - t_sat) / t_sat).abs() < 1e-6, "{key}: {up}");
    }
}
