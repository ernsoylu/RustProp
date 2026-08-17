//! Regression for the props_api cfg gate: the PCSAFT route must be reachable
//! from a build whose ONLY engine is pcsaft (the gate once listed every other
//! engine but not pcsaft, so a pcsaft-only build silently had no `props_si`).
#![cfg(feature = "pcsaft")]

#[test]
fn pcsaft_feature_alone_serves_props_si() {
    // Oracle: CoolProp 8.0.0 wheel, PropsSI("Z","T",300,"P",101325,"PCSAFT::TOLUENE").
    let z = rustprop::props_si("Z", "T", 300.0, "P", 101325.0, "PCSAFT::TOLUENE").unwrap();
    assert!(((z - 0.004400302347313811) / z).abs() < 1e-9, "Z = {z}");
}
