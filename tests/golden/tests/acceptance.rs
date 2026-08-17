//! The final acceptance suite (PLAN.md 15.3): a SEEDED randomized sweep of
//! the state space across every engine, replayed against the wheel's answers.
//!
//! Every other suite here samples states a human chose, which means they
//! sample where someone expected trouble. This one draws pseudo-randomly
//! from a committed seed, so its coverage owes nothing to anyone's
//! intuition. `#[ignore]`d and run by the scheduled CI job.

use std::path::Path;

/// Per-engine tolerance, matching what each phase established: direct
/// evaluations are tight, solver-driven pairs get the solver tier.
fn rtol_for(rec: &rustprop_golden_tests::GoldenRecord) -> f64 {
    // A saturation flash whose pressure is a fraction of a pascal is at the
    // numerical floor of any EOS: PC-SAFT propane at 101 K has a vapour
    // density of 2e-6 mol/m^3, and the inner secant's converged point moves
    // by ~1e-5 between implementations. Legitimate states, but not states
    // the phase suites' 1e-7 tier was established on.
    // Detected by either face of the same corner, since a (Q, T) record does
    // not carry its pressure: a sub-pascal pressure anywhere in the record,
    // or a vapour density below a milli-mol per cubic metre.
    let saturation = rec.name1 == "Q" || rec.name2 == "Q";
    let near_vacuum = (rec.name1 == "P" && rec.val1 < 1.0)
        || (rec.name2 == "P" && rec.val2 < 1.0)
        || (rec.out == "P" && rec.expected < 1.0)
        || (rec.out == "D" && rec.expected < 1e-3);
    if saturation && near_vacuum {
        return 1e-4;
    }
    match (rec.backend.as_str(), rec.out.as_str()) {
        // Transport and surface tension carry their own established tier.
        // Prandtl is transport-derived (cp*eta/lambda), so it rides the
        // same tier for every backend, IF97 included.
        (_, "V") | (_, "L") | (_, "I") | (_, "Prandtl") => 1e-8,
        // Derivative-coefficient outputs are Jacobian ratios of first/second
        // partials evaluated at a flash-solved state, so the flash's own
        // 1e-9 tier amplifies through the ratio (observed: beta on an (H,P)
        // n-Propane draw at 2.7e-8). One decade of headroom over the
        // observed worst.
        (_, "isobaric_expansion_coefficient")
        | (_, "isothermal_compressibility")
        | (_, "PIP")
        | (_, "isentropic_expansion_coefficient")
        | (_, "fundamental_derivative_of_gas_dynamics") => 1e-7,
        ("IF97", _) => 1e-11,
        ("INCOMP", _) => 1e-10,
        ("HA", _) => 1e-8,
        ("PCSAFT", _) => 1e-7,
        _ => 1e-8,
    }
}

#[test]
#[ignore]
fn acceptance_sweep_matches_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/acceptance_sweep.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    // 3,720 records from the original 15.3 seed (20260807, bitwise frozen)
    // plus 1,765 from the 2026-08 extension seed (20260816): 865 mixtures +
    // blends, 620 wide outputs, 140 IF97 flash pairs, 40 PC-SAFT Z, 100
    // pseudo-pure transport; plus 535 from the output-tail seed (20260817):
    // 330 HEOS output tail, 55 environmental trivials, 150 cubics output
    // tail.
    assert_eq!(recs.len(), 6020);

    // The one place this port and upstream part company, kept explicit
    // rather than filtered out of the fixture so it cannot quietly widen.
    //
    // Cubic PQ flashes at near-vacuum pressures: upstream's equal-Gibbs
    // secant converges where this port's gives up, at the extreme cold end
    // of the cubic's own saturation range (SRK::CarbonDioxide bottoms out at
    // 91 K / 0.18 Pa — CO2's real triple point is 217 K, so these are pure
    // extrapolation). The seed, step, tolerance and iteration cap are
    // upstream's verbatim; the difference is root conditioning inside the
    // residual, not solver structure. The original draws failed only below
    // one pascal; the 2026-08 extension draws hit the same secant give-up at
    // 1.3 and 1.95 Pa, so the gate is the observed envelope's decade, 10 Pa.
    // The allowance excuses only ERRORS — a converged-but-wrong value in
    // this band still fails — and the heal-detection assert below keeps it
    // honest.
    let known_divergence = |rec: &rustprop_golden_tests::GoldenRecord| {
        (rec.backend == "SRK" || rec.backend == "PR")
            && rec.name1 == "P"
            && rec.name2 == "Q"
            && rec.val1 < 10.0
    };

    // Three mixture records where the WHEEL's recorded answer is demonstrably
    // not its own equilibrium, each triaged to a specific upstream mechanism
    // and pinned to the PORT's answer so the divergence can neither widen nor
    // silently heal. Tuples: (fluid, out, val1, val2, pinned port value).
    //
    // 1. Methane[0.7]&Ethane[0.3] (P,S)->Q and
    // 3. CarbonDioxide[0.3]&Water[0.7] (P,S)->H: upstream's mixture HSU_P
    //    residual sweeps T through PT flashes on the SHARED backend object;
    //    the Tmax-endpoint evaluation leaves stale SatL/SatV state that
    //    disables the two-phase split for every subsequent in-dome PT flash,
    //    so the wheel converges on a corrupted branch and CONTRADICTS ITS OWN
    //    fresh PT flash at the T it returns. Verified both times: a fresh
    //    wheel AbstractState PT flash at the port's converged T reproduces
    //    the port's out/D/S BITWISE. The port's residual is a pure function
    //    of (z, T, p) — the corruption is history-dependence its stateless
    //    architecture deliberately cannot express.
    // 2. CarbonDioxide[0.3]&Water[0.7] (H,P)->D: shallow Michelsen TPD
    //    minimum just inside the dew line at 44.3 MPa; upstream's stability
    //    solve detects the split at only one T in the convergence window
    //    (where port and wheel agree to 8 digits) and otherwise lands the
    //    metastable single-phase gas root; the port detects the split
    //    consistently and lands the lower-Gibbs two-phase branch.
    let mixture_divergences: &[(&str, &str, f64, f64, f64)] = &[
        (
            "Methane[0.7]&Ethane[0.30000000000000004]",
            "Q",
            4264646.527429062,
            2763.1439307662727,
            0.7559507092594944,
        ),
        (
            "CarbonDioxide[0.3]&Water[0.7]",
            "D",
            1486893.968297218,
            44318484.41637476,
            391.5831239366091,
        ),
        (
            "CarbonDioxide[0.3]&Water[0.7]",
            "H",
            69880188.09130834,
            256.87540877887216,
            137867.16749380616,
        ),
    ];

    let mut failures = Vec::new();
    let mut diverged = 0usize;
    let mut diverged_mixture = 0usize;
    let mut by_backend: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    for rec in &recs {
        // Mixtures and predefined blends ride backend "HEOS" (the identity
        // lives in the fluid string); split them out of the coverage tally
        // so the weekly log shows them.
        let tally_key =
            if rec.fluid.contains('&') || rec.fluid.to_ascii_lowercase().ends_with(".mix") {
                format!("{}(mix)", rec.backend)
            } else {
                rec.backend.clone()
            };
        *by_backend.entry(tally_key).or_default() += 1;
        let got = if rec.backend == "HA" {
            rustprop::ha_props_si(
                &rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &rec.name3, rec.val3,
            )
        } else {
            let fluid = format!("{}::{}", rec.backend, rec.fluid);
            rustprop::props_si(&rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &fluid)
        };
        match got {
            Ok(v) => {
                let rtol = rtol_for(rec);
                // The project's standing tolerance policy: caloric outputs
                // are measured against their THERMAL SCALE where they cross
                // zero, because a relative test on a quantity passing
                // through zero measures nothing. R for entropy and heat
                // capacities, R*T for enthalpy and internal energy.
                let scale = match rec.out.as_str() {
                    // The ideal-gas entropies ride the entropy scale: s0 is
                    // reference-offset exactly like S (this fixture's
                    // smallest Smolar_idealgas is 28 J/mol/K, but s0 falls
                    // with ln(rho) and crosses zero at high pressure).
                    "S" | "SMASS" | "SMOLAR" | "C" | "CVMASS" | "CMASS" | "Smolar_residual"
                    | "Smolar_idealgas" | "Smass_idealgas" => Some(8.314_462_618),
                    // Free energies and residual enthalpy cross zero exactly
                    // like H/U do (Gmolar and Helmholtzmolar sit at O(100)
                    // J/mol near ambient), so they take the same thermal
                    // scale. The 20260817 tail joins the family: residual
                    // Gibbs crosses zero in its draws (14 negative / 8
                    // positive, smallest |value| 2.9e-3 J/mol), and the
                    // reference-offset ideal-gas calorics approach it
                    // (Umolar_idealgas bottoms out at 59 J/mol, well under
                    // R*300 — a relative check there measures flash noise,
                    // not agreement).
                    "H" | "U" | "HMASS" | "UMASS" | "G" | "Gmolar" | "Helmholtzmolar"
                    | "Hmolar_residual" | "Helmholtzmass" | "Gmolar_residual"
                    | "Hmolar_idealgas" | "Umolar_idealgas" | "Hmass_idealgas"
                    | "Umass_idealgas" => Some(8.314_462_618 * 300.0),
                    // Dimensionless O(1) quantities that legitimately cross
                    // zero: PIP passes through 0 between its liquid (~-7)
                    // and vapour (~+1) branches; the expansion coefficient
                    // crosses 0 at water's density maximum. Scales are their
                    // natural magnitudes.
                    "PIP" => Some(1.0),
                    "isobaric_expansion_coefficient" => Some(1e-3),
                    // Reduced-Helmholtz values and residual delta/tau
                    // derivatives: dimensionless O(1..10), crossing zero
                    // between attraction- and repulsion-dominated states
                    // (alphar draws split 19 negative / 5 positive with
                    // smallest |value| 3.8e-6; dalphar_ddelta 20/6,
                    // dalphar_dtau 20/1; alpha0 dips to |0.44|). The
                    // ideal-gas delta derivatives are signed powers of
                    // 1/delta and dalpha0_dtau bottoms out at 1.33 in the
                    // fixture, so they need no guard.
                    "alphar" | "alpha0" | "dalphar_ddelta_consttau" | "dalphar_dtau_constdelta" => {
                        Some(1.0)
                    }
                    // Virial coefficients cross zero at their Boyle points
                    // (fixture sign splits: Bvirial 20 negative / 2
                    // positive, Cvirial 3/21, dBvirial_dT 2/20 — quantum
                    // fluids above B's maximum — dCvirial_dT 15/6). Each
                    // scale is its observed median-magnitude decade:
                    // |B| med 4.2e-4, |C| med 5.7e-9, |dB/dT| med 6.6e-7,
                    // |dC/dT| med 2.6e-11.
                    "Bvirial" => Some(1e-4),
                    "Cvirial" => Some(1e-8),
                    "dBvirial_dT" => Some(1e-6),
                    "dCvirial_dT" => Some(1e-10),
                    _ => None,
                };
                let ok = match scale {
                    // A flash converged to its own 1e-9 tier still moves a
                    // near-zero entropy by ~1e-7 of R (dS/dT ~ cp/T), so the
                    // scale-relative tier is the solver's, not the direct
                    // evaluation's.
                    Some(sc) if rec.expected.abs() < sc => (v - rec.expected).abs() <= 1e-7 * sc,
                    _ => rustprop_golden_tests::check(rec, v, rtol).is_ok(),
                };
                if !ok {
                    let pinned = mixture_divergences.iter().find(|(fl, out, v1, v2, _)| {
                        rec.backend == "HEOS"
                            && rec.fluid == *fl
                            && rec.out == *out
                            && rec.val1 == *v1
                            && rec.val2 == *v2
                    });
                    if let Some((_, _, _, _, port_pin)) = pinned {
                        // The excused records must still match the PINNED
                        // port answer — drifting to a THIRD value fails.
                        if ((v - port_pin) / port_pin).abs() <= 1e-6 {
                            diverged_mixture += 1;
                        } else {
                            failures.push(format!(
                                "{}: diverged from BOTH the wheel ({}) and the pinned port answer ({port_pin}): {v}",
                                rec.id(),
                                rec.expected
                            ));
                        }
                    } else {
                        failures.push(
                            rustprop_golden_tests::check(rec, v, rtol)
                                .err()
                                .unwrap_or_else(|| {
                                    format!("{}: {v} vs {}", rec.id(), rec.expected)
                                }),
                        );
                    }
                }
            }
            Err(e) => {
                if known_divergence(rec) {
                    diverged += 1;
                } else {
                    failures.push(format!(
                        "{}: errored where upstream answered: {e}",
                        rec.id()
                    ));
                }
            }
        }
    }
    println!("acceptance sweep coverage: {by_backend:?}");
    // The assert message truncates at 30; a randomized weekly suite wants
    // the complete list on the log for triage.
    for f in &failures {
        println!("FAIL {f}");
    }
    assert!(
        failures.is_empty(),
        "{} of {} records failed:\n{}",
        failures.len(),
        recs.len(),
        failures
            .iter()
            .take(30)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    // If the cubic secant ever gets fixed, this fires and the allowance
    // above should be deleted — a known divergence that silently heals is
    // still a lie in the test.
    assert!(
        diverged > 0,
        "the documented sub-pascal cubic divergence no longer reproduces; remove the allowance"
    );
    // Same heal detection for the three pinned mixture records: if one comes
    // to match the wheel (e.g. after an upstream-side re-generation or a
    // solver change), its excuse never fires and this forces removing it.
    assert_eq!(
        diverged_mixture,
        mixture_divergences.len(),
        "a pinned mixture divergence no longer reproduces; remove its entry"
    );
    println!(
        "acceptance sweep: {} records, {diverged} known cubic divergences, {diverged_mixture} pinned mixture divergences",
        recs.len()
    );
}
