//! Generic partial-derivative goldens (PLAN.md Phase 12 slice 12a): the
//! (T, rho)-basis machinery the Tabular table build consumes, checked
//! against the wheel's own `d(X)/d(Y)|Z` derivative-string outputs.

use rustprop_core::params::Param;
use rustprop_heos::derivs::{first_partial_deriv, second_partial_deriv};
use std::path::Path;

fn parse_first(out: &str) -> Option<(Param, Param, Param)> {
    // "d(T)/d(Hmolar)|P"
    let rest = out.strip_prefix("d(")?;
    let (of, rest) = rest.split_once(")/d(")?;
    let (wrt, constant) = rest.split_once(")|")?;
    Some((
        Param::parse(of)?,
        Param::parse(wrt)?,
        Param::parse(constant)?,
    ))
}

fn parse_second(out: &str) -> Option<(Param, Param, Param, Param, Param)> {
    // "d(d(T)/d(Hmolar)|P)/d(P)|Hmolar"
    let rest = out.strip_prefix("d(d(")?;
    let (of, rest) = rest.split_once(")/d(")?;
    let (wrt1, rest) = rest.split_once(")|")?;
    let (const1, rest) = rest.split_once(")/d(")?;
    let (wrt2, const2) = rest.split_once(")|")?;
    Some((
        Param::parse(of)?,
        Param::parse(wrt1)?,
        Param::parse(const1)?,
        Param::parse(wrt2)?,
        Param::parse(const2)?,
    ))
}

#[test]
fn partial_derivs_match_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/partial_derivs.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 207);
    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();

    let mut failures = Vec::new();
    for rec in &recs {
        let data = registry[rec.fluid.as_str()];
        let flash = rustprop_heos::flash_pt::PtFlash::new(data);
        // (T, P) state -> density through the ported PT flash
        let (rho, _phase) = flash.pt_flash(rec.val1, rec.val2).expect("PT flash");
        let t = rec.val1;
        let is_second = rec.out.starts_with("d(d(");
        let actual = if is_second {
            let (of, w1, c1, w2, c2) = parse_second(&rec.out).expect("second-deriv string");
            second_partial_deriv(&flash.eos, of, w1, c1, w2, c2, t, rho).expect("second deriv")
        } else {
            let (of, wrt, constant) = parse_first(&rec.out).expect("first-deriv string");
            first_partial_deriv(&flash.eos, of, wrt, constant, t, rho).expect("first deriv")
        };
        // The formulas are bitwise; agreement is bounded by the density the
        // PT flash hands both sides (its own ~1e-9 convergence). Second
        // derivatives differentiate that error once more — the compressed-
        // liquid CO2 d2(Smolar)/dT2|P record observes 2.4e-9.
        let rtol = if is_second { 1e-8 } else { 1e-9 };
        if let Err(e) = rustprop_golden_tests::check(rec, actual, rtol) {
            failures.push(e);
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} partial-derivative records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}
