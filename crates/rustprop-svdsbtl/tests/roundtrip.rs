//! Scratch harness: load a converted artifact from RUSTPROP_SVDS and print
//! evaluations. Ignored by default (needs a locally built upstream table).
use rustprop_core::params::Param;

#[test]
#[ignore]
fn eval_from_env_artifact() {
    let path = std::env::var("RUSTPROP_SVDS").expect("set RUSTPROP_SVDS");
    let bytes = std::fs::read(&path).expect("read artifact");
    let surf = rustprop_svdsbtl::artifact::load("Water", &bytes).expect("load");
    println!(
        "input_pair {} regions {} props {:?}",
        surf.input_pair,
        surf.region_count(),
        surf.properties()
    );
    // PT surface: (a, b) = (p, T)
    for (p, t) in [(101325.0, 400.0), (1.0e6, 500.0), (5.0e7, 300.0)] {
        let r = surf.resolve(p, t);
        println!(
            "p={p} T={t} region={:?} D={:?} H={:?} S={:?}",
            r.map(|r| r.region_idx),
            surf.eval(Param::Dmass, p, t),
            surf.eval(Param::Hmass, p, t),
            surf.eval(Param::Smass, p, t)
        );
    }
}
