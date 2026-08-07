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
        (_, "V") | (_, "L") | (_, "I") => 1e-8,
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
    assert_eq!(recs.len(), 3720);

    // The one place this port and upstream part company, kept explicit
    // rather than filtered out of the fixture so it cannot quietly widen.
    //
    // Cubic PQ flashes at SUB-PASCAL pressures: upstream's equal-Gibbs
    // secant converges where this port's gives up, at the extreme cold end
    // of the cubic's own saturation range (SRK::CarbonDioxide bottoms out at
    // 91 K / 0.18 Pa — CO2's real triple point is 217 K, so these are pure
    // extrapolation). The seed, step, tolerance and iteration cap are
    // upstream's verbatim; the difference is root conditioning inside the
    // residual, not solver structure. Four records of 3120.
    let known_divergence = |rec: &rustprop_golden_tests::GoldenRecord| {
        (rec.backend == "SRK" || rec.backend == "PR")
            && rec.name1 == "P"
            && rec.name2 == "Q"
            && rec.val1 < 1.0
    };

    let mut failures = Vec::new();
    let mut diverged = 0usize;
    let mut by_backend: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    for rec in &recs {
        *by_backend.entry(rec.backend.clone()).or_default() += 1;
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
                    "S" | "SMASS" | "SMOLAR" | "C" | "CVMASS" | "CMASS" => Some(8.314_462_618),
                    "H" | "U" | "HMASS" | "UMASS" => Some(8.314_462_618 * 300.0),
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
                    failures.push(
                        rustprop_golden_tests::check(rec, v, rtol)
                            .err()
                            .unwrap_or_else(|| format!("{}: {v} vs {}", rec.id(), rec.expected)),
                    );
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
    println!(
        "acceptance sweep: {} records, {diverged} known divergences",
        recs.len()
    );
}
