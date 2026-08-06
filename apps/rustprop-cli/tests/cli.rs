//! End-to-end CLI tests (PLAN.md 2.5): the binary's stdout must match the
//! golden oracle values.

use std::process::Command;

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rustprop-cli"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn props_tsat_matches_golden() {
    let out = run(&["props", "T", "P", "101325", "Q", "0", "IF97::Water"]);
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let v: f64 = String::from_utf8(out.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    // Golden: PropsSI("T", "P", 101325, "Q", 0, "IF97::Water") from the
    // CoolProp 8.0.0 oracle (tests/golden/fixtures/if97_water.jsonl).
    let expected = 373.12430000048056;
    assert!(
        ((v - expected) / expected).abs() < 1e-11,
        "got {v}, expected {expected}"
    );
}

#[test]
fn props_enthalpy_matches_golden() {
    let out = run(&["props", "H", "T", "300", "P", "101325", "IF97::Water"]);
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let v: f64 = String::from_utf8(out.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    let expected = 112665.04341853978;
    assert!(
        ((v - expected) / expected).abs() < 1e-11,
        "got {v}, expected {expected}"
    );
}

#[test]
fn unknown_fluid_fails_cleanly() {
    let out = run(&["props", "T", "P", "101325", "Q", "0", "HEOS::Water"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("unknown fluid"));
}

#[test]
fn unknown_parameter_fails_cleanly() {
    let out = run(&["props", "XYZ", "P", "101325", "Q", "0", "IF97::Water"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("unknown parameter"));
}
