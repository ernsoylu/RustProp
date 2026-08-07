//! Tabular table-construction goldens (PLAN.md Phase 12 slice 12b): grid
//! limits computed exactly as upstream's `set_limits`, and node values,
//! which are plain source-backend evaluations at the grid coordinates.

use rustprop_heos::flash_pt::PtFlash;
use rustprop_tabular::tables::{GridKind, GriddedTable};
use std::path::Path;

fn fluid(name: &str) -> &'static rustprop_core::fluid::FluidData {
    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    registry[name]
}

#[test]
fn tabular_tables_match_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/tabular_tables.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 140);

    // One coarse LogPT table per fluid covers every node record; the limits
    // are grid-size independent so the same table answers those too.
    let mut tables: std::collections::HashMap<String, (GriddedTable, GriddedTable)> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &recs {
        let entry = tables.entry(rec.fluid.clone()).or_insert_with(|| {
            let flash = PtFlash::new(fluid(&rec.fluid));
            // 2x2 is enough for the limits; the node checks use the 20x20.
            let ph =
                GriddedTable::build(&flash, GridKind::LogPH, 2, 2, None).expect("LogPH limits");
            let pt =
                GriddedTable::build(&flash, GridKind::LogPT, 20, 20, None).expect("LogPT table");
            (ph, pt)
        });
        let (ph, pt) = (&entry.0, &entry.1);

        let actual = match (rec.name1.as_str(), rec.out.as_str()) {
            ("LogPH", "xmin") => ph.xmin,
            ("LogPH", "xmax") => ph.xmax,
            ("LogPH", "ymin") => ph.ymin,
            ("LogPH", "ymax") => ph.ymax,
            ("LogPT", "xmin") => pt.xmin,
            ("LogPT", "xmax") => pt.xmax,
            ("LogPT", "ymin") => pt.ymin,
            ("LogPT", "ymax") => pt.ymax,
            ("LogPT_node", out) => {
                // Locate the node by its (x, y) coordinates.
                let i = pt
                    .xvec
                    .iter()
                    .position(|v| (v - rec.val1).abs() <= 1e-9 * rec.val1.abs().max(1.0))
                    .unwrap_or_else(|| panic!("no x node at {}", rec.val1));
                let j = pt
                    .yvec
                    .iter()
                    .position(|v| (v - rec.val2).abs() <= 1e-9 * rec.val2.abs())
                    .unwrap_or_else(|| panic!("no y node at {}", rec.val2));
                match out {
                    "Dmolar" => pt.rhomolar.val[i][j],
                    "Hmolar" => pt.hmolar.val[i][j],
                    "Smolar" => pt.smolar.val[i][j],
                    "Umolar" => pt.umolar.val[i][j],
                    other => panic!("unknown node output {other}"),
                }
            }
            (kind, out) => panic!("unknown record {kind}/{out}"),
        };
        // Limits and node values are direct evaluations; the only solver in
        // the path is the PT/QT flash the node coordinate implies. Caloric
        // outputs are reference-state anchored and cross zero near the
        // triple point, so they are measured against the thermal scale
        // (R*Tc for H/U, R for S) as everywhere else in this suite.
        let scale = match rec.out.as_str() {
            "Hmolar" | "Umolar" => {
                let flash = PtFlash::new(fluid(&rec.fluid));
                (8.314_462_618_153_24 * flash.t_critical()).max(rec.expected.abs())
            }
            "Smolar" => 8.314_462_618_153_24_f64.max(rec.expected.abs()),
            _ => rec.expected.abs(),
        };
        let err = (actual - rec.expected).abs() / scale;
        // (an inf/NaN error must fail, not silently pass a comparison)
        if err > 1e-9 || err.is_nan() {
            failures.push(format!(
                "{}: actual {actual:e} vs expected {:e} (scaled err {err:e} > 1e-9)",
                rec.id(),
                rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} tabular table records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

/// Upstream's tabular backends are low-level only:
/// `TabularBackend::available_in_high_level()` returns false, so PropsSI
/// rejects them with this exact message (verified against the wheel, which
/// raises the identical string for `PropsSI(..., "TTSE&HEOS::Water")`).
#[test]
fn tabular_is_rejected_by_high_level_api() {
    use rustprop_core::Error;
    for backend in ["TTSE&HEOS::Water", "BICUBIC&HEOS::Water"] {
        match rustprop::props_si("Dmolar", "T", 400.0, "P", 101325.0, backend).unwrap_err() {
            Error::Value(m) => assert_eq!(
                m,
                "This AbstractState derived class cannot be used in the high-level interface; see www.coolprop.org/dev/coolprop/LowLevelAPI.html"
            ),
            other => panic!("wrong variant for {backend}: {other:?}"),
        }
    }
}

/// The low-level path (upstream `AbstractState::factory("TTSE&HEOS", ...)`)
/// works, and rejects out-of-range inputs with upstream's verbatim message.
/// Values are the wheel's, bitwise. (The bulk state goldens live in
/// `tabular_state.rs`; this test pins the errors and the entry point.)
///
/// FIDELITY FINDING, documented rather than asserted: upstream's PT
/// two-phase rejection ("P,T with TTSE cannot be two-phase for now") is
/// unreachable for pure fluids. `SatTable::is_inside` brackets T between the
/// cubic-interpolated liquid and vapour saturation temperatures, which for a
/// pure fluid are the same curve — the "inside" set has ulp width. The wheel
/// confirms it: at exactly Ts(101325 Pa) = 373.12429584766636 K it returns
/// the liquid root (53196.19 mol/m^3) and one nanokelvin above it returns the
/// vapour root (33.17 mol/m^3), with no error in between. The branch is
/// carried here for shape parity with upstream.
#[test]
fn tabular_low_level_state_pt() {
    use rustprop_core::params::Param;
    use rustprop_tabular::{Scheme, TabularState};

    let flash = PtFlash::new(fluid("Water"));
    for (scheme, expected) in [
        (Scheme::Ttse, 30.804_110_745_052_38),
        (Scheme::Bicubic, 30.804_082_705_024_516),
    ] {
        let mut st = TabularState::new(scheme, &flash, 200, 200, None).expect("tables");
        st.update_pt(101_325.0, 400.0).expect("PT update");
        assert_eq!(st.keyed_output(Param::Dmolar).expect("Dmolar"), expected);
        assert_eq!(st.keyed_output(Param::T).expect("T"), 400.0);
        assert_eq!(st.keyed_output(Param::P).expect("P"), 101_325.0);

        match st.update_pt(1.0e-3, 400.0) {
            Err(rustprop_core::Error::Value(m)) => {
                assert_eq!(m, "inputs are not in range, p=0.001 Pa, T=400 K");
            }
            other => panic!("expected the range rejection, got {other:?}"),
        }
    }
}
