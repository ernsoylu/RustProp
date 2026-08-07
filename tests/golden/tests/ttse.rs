//! TTSE evaluation goldens (PLAN.md Phase 12 slice 12c): the second-order
//! Taylor expansion on the LogPT grid, against the wheel's own TTSE backend
//! (which builds the same 200x200 table from the same source engine).

use rustprop_core::params::Param;
use rustprop_heos::flash_pt::PtFlash;
use rustprop_tabular::tables::{GridKind, GriddedTable};
use rustprop_tabular::ttse::evaluate_single_phase;
use std::path::Path;

fn fluid(name: &str) -> &'static rustprop_core::fluid::FluidData {
    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    registry[name]
}

#[test]
fn ttse_matches_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/ttse.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 144);

    let mut tables: std::collections::HashMap<String, GriddedTable> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &recs {
        let table = tables.entry(rec.fluid.clone()).or_insert_with(|| {
            let flash = PtFlash::new(fluid(&rec.fluid));
            GriddedTable::build(&flash, GridKind::LogPT, 200, 200, None).expect("LogPT table")
        });
        let (t, p) = (rec.val1, rec.val2);
        let (i, j) = table
            .find_native_nearest_good_neighbor(t, p)
            .expect("nearest good node");
        let out = match rec.out.as_str() {
            "Dmolar" => Param::Dmolar,
            "Hmolar" => Param::Hmolar,
            "Smolar" => Param::Smolar,
            "Umolar" => Param::Umolar,
            other => panic!("unknown output {other}"),
        };
        let actual = evaluate_single_phase(table, out, t, p, i, j).expect("TTSE eval");
        // The whole path — grid limits, spacing, node selection, expansion —
        // is deterministic on both sides, so this is an exactness check, not
        // an interpolation-accuracy one. Caloric outputs still ride the
        // thermal scale where they cross zero.
        let scale = match rec.out.as_str() {
            "Hmolar" | "Umolar" => {
                let flash = PtFlash::new(fluid(&rec.fluid));
                (8.314_462_618_153_24 * flash.t_critical()).max(rec.expected.abs())
            }
            "Smolar" => 8.314_462_618_153_24_f64.max(rec.expected.abs()),
            _ => rec.expected.abs(),
        };
        let err = (actual - rec.expected).abs() / scale;
        if err > 1e-9 || err.is_nan() {
            failures.push(format!(
                "{}: actual {actual:e} vs expected {:e} (scaled err {err:e})",
                rec.id(),
                rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} TTSE records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn ttse_is_exact_at_nodes() {
    // At a grid node the expansion has zero deltas, so it must return the
    // stored node value bitwise.
    let flash = PtFlash::new(fluid("Water"));
    let table = GriddedTable::build(&flash, GridKind::LogPT, 60, 60, None).expect("table");
    let mut checked = 0;
    for i in [10usize, 25, 44] {
        for j in [5usize, 30, 55] {
            if !table.t.val[i][j].is_finite() {
                continue;
            }
            let (x, y) = (table.xvec[i], table.yvec[j]);
            for (param, grid) in [
                (Param::Dmolar, &table.rhomolar),
                (Param::Hmolar, &table.hmolar),
                (Param::Smolar, &table.smolar),
                (Param::Umolar, &table.umolar),
            ] {
                let v = evaluate_single_phase(&table, param, x, y, i, j).expect("eval");
                assert_eq!(v, grid.val[i][j], "node ({i},{j}) {param:?} not exact");
                checked += 1;
            }
        }
    }
    assert!(checked >= 20, "expected a meaningful number of nodes");
}
