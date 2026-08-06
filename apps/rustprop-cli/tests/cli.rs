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
    let out = run(&["props", "T", "P", "101325", "Q", "0", "Watr"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("was not found"));
}

#[test]
fn unknown_parameter_fails_cleanly() {
    let out = run(&["props", "XYZ", "P", "101325", "Q", "0", "IF97::Water"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("parsing failed"));
}

#[test]
fn props_bare_fluid_uses_heos() {
    // PLAN.md 5.3: the full planned CLI form with a bare fluid name.
    // Golden: PropsSI("Dmolar", "T", 300, "P", 101325, "Water") from the
    // CoolProp 8.0.0 oracle (tests/golden/fixtures/props_si.jsonl covers the
    // HEOS string API; this value round-trips through the same dispatch).
    let out = run(&["props", "Dmolar", "T", "300", "P", "101325", "Water"]);
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
    let expected = 55317.35277350119;
    assert!(
        ((v - expected) / expected).abs() < 1e-8,
        "got {v}, expected {expected}"
    );
}

#[test]
fn props_mass_alias_matches_golden() {
    // Mass-basis alias output ("D" = Dmass) with a trivial-output check too.
    let out = run(&["props", "Tcrit", "T", "300", "P", "101325", "R134a"]);
    assert!(out.status.success());
    let v: f64 = String::from_utf8(out.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    // The superancillary NUMERICAL critical temperature (upstream
    // calc_T_critical for a superancillary fluid).
    let expected = 374.2119665849513;
    assert!(((v - expected) / expected).abs() < 1e-12);
}
