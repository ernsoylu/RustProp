//! Randomized SVDSBTL acceptance over the three committed low-res PT
//! artifacts (Water plus the n-Propane and CarbonDioxide surfaces added
//! with this suite — genuinely new atlases, boundary splines and
//! coefficient sets; Water stays small because the 745 fixed-grid goldens
//! already blanket its atlas).
//!
//! `#[ignore]`d only for acceptance-family symmetry — it runs in well under
//! a second (artifact load + direct evals; no table builds).

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
#[ignore]
fn acceptance_svdsbtl_matches_oracle() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let recs = rustprop_golden_tests::load_jsonl(&dir.join("fixtures/acceptance_svdsbtl.jsonl"));
    // (40 Water + 80 n-Propane + 80 CarbonDioxide) states x 5 outputs.
    assert_eq!(recs.len(), 1000);

    let mut failures = Vec::new();
    let mut exact = 0usize;
    let mut worst: f64 = 0.0;
    for fluid in ["Water", "n-Propane", "CarbonDioxide"] {
        let bytes =
            std::fs::read(dir.join(format!("fixtures/svdsbtl/{fluid}.PT.NT16-NR24-r6.svds")))
                .expect("committed artifact");
        let surface = artifact::load(fluid, &bytes).expect("load");
        // PT_INPUTS; region_count varies by fluid and stays unasserted.
        assert_eq!(surface.input_pair, 17);
        for rec in recs.iter().filter(|r| r.fluid == fluid) {
            // The PT surface takes (a, b) = (p, T).
            let got = surface
                .eval(param(&rec.out), rec.val1, rec.val2)
                .expect("eval");
            if got == rec.expected {
                exact += 1;
                continue;
            }
            // The wheel only recorded finite answers; a NaN here means the
            // ported region dispatch disagreed on identical doubles — a
            // real defect, never tolerable noise.
            if !got.is_finite() {
                failures.push(format!(
                    "{fluid} {} at p={} T={}: got {got}, want {}",
                    rec.out, rec.val1, rec.val2, rec.expected
                ));
                continue;
            }
            let rel = ((got - rec.expected) / rec.expected).abs();
            worst = worst.max(rel);
            // 13.2 tier: the fixed-grid suite holds 1e-14 (worst observed
            // 1.8e-15, a GCC -ffp-contract=fast residual in the reference
            // build); one decade of headroom for the two new fluids'
            // unexplored coefficient sets.
            if rel > 1e-13 {
                failures.push(format!(
                    "{fluid} {}",
                    rustprop_golden_tests::check(rec, got, 1e-13).unwrap_err()
                ));
            }
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
        "acceptance_svdsbtl: {exact}/{} bitwise, worst rel err {worst:.3e}",
        recs.len()
    );
}
