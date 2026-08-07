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
        "Q" => Param::Q,
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

/// Every input pair `TabularBackend::update` accepts beyond PT (PLAN.md
/// slice 12f), against the wheel's low-level tabular backends.
///
/// `#[ignore]`d for runtime, not for confidence: each state object builds a
/// 200x200 LogPH grid, which costs an (h, p) flash per node — about 100 s
/// per fluid per scheme. The weekly CI job runs it.
#[test]
#[ignore]
fn tabular_pairs_match_oracle() {
    use rustprop_tabular::TabularInput;

    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/tabular_pairs.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 1632);

    let pair_of = |n1: &str, n2: &str| match (n1, n2) {
        ("Hmolar", "P") => TabularInput::HmolarP,
        ("P", "Umolar") => TabularInput::PUmolar,
        ("P", "Smolar") => TabularInput::PSmolar,
        ("Dmolar", "P") => TabularInput::DmolarP,
        ("Smolar", "T") => TabularInput::SmolarT,
        ("Dmolar", "T") => TabularInput::DmolarT,
        ("P", "Q") => TabularInput::PQ,
        ("Q", "T") => TabularInput::QT,
        other => panic!("unknown pair {other:?}"),
    };

    let mut failures = Vec::new();
    let mut worst: f64 = 0.0;
    let mut checked = 0usize;
    for fluid_name in ["Water", "n-Propane"] {
        let flash = PtFlash::new(fluid(fluid_name));
        for scheme_name in ["TTSE", "BICUBIC"] {
            let scheme = if scheme_name == "TTSE" {
                Scheme::Ttse
            } else {
                Scheme::Bicubic
            };
            let mut st = TabularState::with_defaults(scheme, &flash, None).expect("tables");
            let mut last = (f64::NAN, f64::NAN, TabularInput::PT);
            for rec in recs
                .iter()
                .filter(|r| r.fluid == fluid_name && r.backend == scheme_name)
            {
                let pair = pair_of(&rec.name1, &rec.name2);
                if last != (rec.val1, rec.val2, pair) {
                    st.update(pair, rec.val1, rec.val2).unwrap_or_else(|e| {
                        panic!("{scheme_name} {fluid_name} {pair:?} update: {e}")
                    });
                    last = (rec.val1, rec.val2, pair);
                }
                let got = st.keyed_output(param(&rec.out)).expect("output");
                if got != rec.expected && rec.expected != 0.0 {
                    worst = worst.max(((got - rec.expected) / rec.expected).abs());
                }
                checked += 1;
                // The LogPH grid's nodes come from (h, p) FLASHES, which are
                // iterative: my node densities sit 5e-11..3e-10 from the
                // wheel's, both inside solver tolerance. Inverting the
                // bicubic cell for a density target amplifies that — worst
                // observed 2.0e-6, at n-Propane's lowest-pressure corner
                // (rho = 9e-6 mol/m^3) where the cell is most nonlinear. Those
                // records get 1e-5; everything else holds 1e-9.
                let rtol = if scheme == Scheme::Bicubic && rec.name1 == "Dmolar" {
                    1e-5
                } else {
                    1e-9
                };
                if let Err(e) = rustprop_golden_tests::check(rec, got, rtol) {
                    failures.push(format!("{scheme_name} {fluid_name} {e}"));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {checked} records failed:\n{}",
        failures.len(),
        failures
            .iter()
            .take(25)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    println!("tabular_pairs: {checked} records, worst rel err {worst:.3e}");
}
