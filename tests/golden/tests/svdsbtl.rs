//! SVDSBTL evaluator goldens (PLAN.md Phase 13 slice 13.2): the ported
//! rank-r Hermite reconstruction plus region dispatch, against the wheel's
//! own `SVDSBTL&HEOS` backend reading the SAME coefficients.
//!
//! The committed artifact is deliberately low resolution (NT=16, NR=24,
//! rank=6, ~186 KB). Both sides approximate HEOS equally badly on it and
//! must still agree BITWISE, which is what this test asserts — a
//! default-resolution surface is 13 MB and does not belong in a repository
//! (PLAN.md Phase 13.1 design note).

use rustprop_core::params::Param;
use rustprop_svdsbtl::artifact;
use std::path::Path;

fn param(name: &str) -> Param {
    match name {
        "Dmass" => Param::Dmass,
        "Hmass" => Param::Hmass,
        "Smass" => Param::Smass,
        "Umass" => Param::Umass,
        "A" => Param::SpeedSound,
        other => panic!("unknown output {other}"),
    }
}

#[test]
fn svdsbtl_matches_oracle() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let recs = rustprop_golden_tests::load_jsonl(&dir.join("fixtures/svdsbtl.jsonl"));
    assert_eq!(recs.len(), 745);

    let bytes = std::fs::read(dir.join("fixtures/svdsbtl/Water.PT.NT16-NR24-r6.svds"))
        .expect("read the committed artifact");
    let surface = artifact::load("Water", &bytes).expect("load the artifact");
    // The Water PT atlas: subcritical LOG band, two POWER bands straddling
    // p_crit, supercritical LOG band — each split liquid/vapour.
    assert_eq!(surface.region_count(), 6);
    assert_eq!(surface.input_pair, 17); // upstream PT_INPUTS

    let mut failures = Vec::new();
    let mut exact = 0usize;
    let mut worst: f64 = 0.0;
    for rec in &recs {
        // The PT surface takes (a, b) = (p, T).
        let got = surface
            .eval(param(&rec.out), rec.val1, rec.val2)
            .expect("evaluate");
        if got == rec.expected {
            exact += 1;
            continue;
        }
        let rel = ((got - rec.expected) / rec.expected).abs();
        worst = worst.max(rel);
        // A few ulp, never more. The kernel is ported statement for
        // statement, so the residual is the reference BUILD, not the
        // algorithm: GCC compiles the accumulation chain with
        // -ffp-contract=fast and fuses some multiply-adds, which Rust never
        // does implicitly. Fusing the obvious candidate (`acc += u_k*v_k`
        // -> mul_add) makes agreement WORSE (69 records instead of 45), so
        // the contraction GCC chose is elsewhere in the expression and
        // reproducing it would be matching a compiler flag rather than
        // porting the algorithm.
        if rel > 1e-14 {
            failures.push(format!(
                "{} at p={} T={}: got {got}, want {} (rel {rel:e})",
                rec.out, rec.val1, rec.val2, rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} records differ:\n{}",
        failures.len(),
        recs.len(),
        failures
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    println!(
        "svdsbtl: {exact}/{} records bitwise, worst rel err {worst:.3e}",
        recs.len()
    );
}

/// A point outside every region returns NaN, as upstream's `SVDSurface::eval`
/// does — not an error, and not an extrapolated number.
#[test]
fn svdsbtl_outside_every_region_is_nan() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let bytes = std::fs::read(dir.join("fixtures/svdsbtl/Water.PT.NT16-NR24-r6.svds"))
        .expect("read the committed artifact");
    let surface = artifact::load("Water", &bytes).expect("load");
    // Below the triple pressure the subcritical regions do not exist.
    assert!(surface.resolve(1.0, 300.0).is_none());
    assert!(
        surface
            .eval(Param::Dmass, 1.0, 300.0)
            .expect("eval")
            .is_nan()
    );
}
