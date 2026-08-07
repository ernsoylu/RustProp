//! Low-level tabular state goldens (PLAN.md Phase 12 slice 12e): upstream's
//! `AbstractState::factory("TTSE&HEOS"/"BICUBIC&HEOS", fluid)` driven through
//! `update(PT_INPUTS, ...)`, covering the state wrapper — mass-basis
//! conversions, echoed inputs, transport accessors — on top of the
//! interpolants that 12c/12d already pin.

use rustprop_core::params::Param;
use rustprop_heos::flash_pt::PtFlash;
use rustprop_tabular::tables::TransportSource;
use rustprop_tabular::{Scheme, TabularState};
use std::path::Path;

fn fluid(name: &str) -> &'static rustprop_core::fluid::FluidData {
    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    registry[name]
}

/// The resolver seam `rustprop-tabular` leaves open: transport needs the
/// fluid registry (ECS conformal references), which lives above that crate.
/// The facade's public string API is exactly that resolution.
struct FacadeTransport(String);

impl TransportSource for FacadeTransport {
    fn viscosity(&self, t: f64, rhomolar: f64) -> Option<f64> {
        rustprop::props_si("V", "T", t, "Dmolar", rhomolar, &self.0).ok()
    }
    fn conductivity(&self, t: f64, rhomolar: f64) -> Option<f64> {
        rustprop::props_si("L", "T", t, "Dmolar", rhomolar, &self.0).ok()
    }
}

fn param(name: &str) -> Param {
    match name {
        "Dmolar" => Param::Dmolar,
        "Dmass" => Param::Dmass,
        "Hmolar" => Param::Hmolar,
        "Hmass" => Param::Hmass,
        "Smolar" => Param::Smolar,
        "Smass" => Param::Smass,
        "Umolar" => Param::Umolar,
        "Umass" => Param::Umass,
        "T" => Param::T,
        "P" => Param::P,
        "Viscosity" => Param::Viscosity,
        "Conductivity" => Param::Conductivity,
        other => panic!("unknown output {other}"),
    }
}

#[test]
fn tabular_state_matches_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/tabular_state.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 552);

    let mut failures = Vec::new();
    let mut worst: f64 = 0.0;
    for fluid_name in ["Water", "n-Propane"] {
        let flash = PtFlash::new(fluid(fluid_name));
        let src = FacadeTransport(format!("HEOS::{fluid_name}"));
        for scheme_name in ["TTSE", "BICUBIC"] {
            let scheme = if scheme_name == "TTSE" {
                Scheme::Ttse
            } else {
                Scheme::Bicubic
            };
            let mut st = TabularState::with_defaults(scheme, &flash, Some(&src)).expect("tables");
            for rec in recs
                .iter()
                .filter(|r| r.fluid == fluid_name && r.backend == scheme_name)
            {
                st.update_pt(rec.val2, rec.val1).expect("PT update");
                let got = st.keyed_output(param(&rec.out)).expect("output");
                if got != rec.expected {
                    worst = worst.max(((got - rec.expected) / rec.expected).abs());
                }
                // 1e-9 matches the tolerance slice 12c established for the
                // TTSE expansion itself; most records here land bitwise.
                if let Err(e) = rustprop_golden_tests::check(rec, got, 1e-9) {
                    failures.push(format!("{scheme_name} {fluid_name} {e}"));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
    println!(
        "tabular_state: {} records, worst rel err {worst:.3e}",
        recs.len()
    );
}
