//! Bicubic evaluation goldens (PLAN.md Phase 12 slice 12d): the
//! 16-coefficient cell interpolant on the LogPT grid, against the wheel's
//! own BICUBIC backend.

use rustprop_core::params::Param;
use rustprop_heos::flash_pt::PtFlash;
use rustprop_tabular::bicubic::{CellCoeffGrid, evaluate_single_phase};
use rustprop_tabular::tables::{GridKind, GriddedTable};
use std::path::Path;

fn fluid(name: &str) -> &'static rustprop_core::fluid::FluidData {
    let registry: std::collections::HashMap<&str, &'static rustprop_core::fluid::FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    registry[name]
}

#[test]
fn bicubic_matches_oracle() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/bicubic.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 160);

    let mut built: std::collections::HashMap<String, (GriddedTable, CellCoeffGrid)> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &recs {
        let entry = built.entry(rec.fluid.clone()).or_insert_with(|| {
            let flash = PtFlash::new(fluid(&rec.fluid));
            let table =
                GriddedTable::build(&flash, GridKind::LogPT, 200, 200, None).expect("table");
            let coeffs = CellCoeffGrid::build(&table);
            (table, coeffs)
        });
        let (table, coeffs) = (&entry.0, &entry.1);
        let (t, p) = (rec.val1, rec.val2);
        let (i, j) = table
            .find_native_nearest_good_cell(t, p)
            .expect("cell search");
        // An invalid cell defers to its remapped good neighbour, as upstream's
        // find_native_nearest_good_indices does before evaluating.
        let (i, j) = match coeffs.cell(i, j).alternate() {
            Some(alt) if !coeffs.cell(i, j).valid() => alt,
            _ => (i, j),
        };
        let out = match rec.out.as_str() {
            "Dmolar" => Param::Dmolar,
            "Hmolar" => Param::Hmolar,
            "Smolar" => Param::Smolar,
            "Umolar" => Param::Umolar,
            other => panic!("unknown output {other}"),
        };
        let actual = evaluate_single_phase(table, coeffs, out, t, p, i, j).expect("bicubic eval");
        let scale = match rec.out.as_str() {
            "Hmolar" | "Umolar" => {
                let flash = PtFlash::new(fluid(&rec.fluid));
                (8.314_462_618_153_24 * flash.t_critical()).max(rec.expected.abs())
            }
            "Smolar" => 8.314_462_618_153_24_f64.max(rec.expected.abs()),
            _ => rec.expected.abs(),
        };
        // 1e-8: a query landing exactly on a grid value is mis-located by
        // upstream's `bisect_vector` sign-product bug (documented in
        // ttse.rs and reproduced here), so it evaluates the cubic well
        // OUTSIDE the chosen cell (xhat ~ -0.37). There the polynomial
        // amplifies the last-ulp difference between Eigen's vectorized
        // 16x16 coefficient product and this port's scalar accumulation.
        // Interior states — the other 156 records — agree below 1e-9.
        let err = (actual - rec.expected).abs() / scale;
        if err > 1e-8 || err.is_nan() {
            failures.push(format!(
                "{}: actual {actual:e} vs expected {:e} (scaled err {err:e})",
                rec.id(),
                rec.expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} bicubic records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn bicubic_reproduces_cell_corners() {
    // At a cell corner (xhat = yhat = 0) the interpolant collapses to
    // alpha[0], which the coefficient construction sets to the node value.
    let flash = PtFlash::new(fluid("Water"));
    let table = GriddedTable::build(&flash, GridKind::LogPT, 60, 60, None).expect("table");
    let coeffs = CellCoeffGrid::build(&table);
    let mut checked = 0;
    for i in [12usize, 30, 47] {
        for j in [8usize, 33, 52] {
            if !coeffs.cell(i, j).valid() {
                continue;
            }
            let (x, y) = (table.xvec[i], table.yvec[j]);
            for (param, grid) in [
                (Param::Dmolar, &table.rhomolar),
                (Param::Hmolar, &table.hmolar),
                (Param::Smolar, &table.smolar),
                (Param::Umolar, &table.umolar),
            ] {
                let v = evaluate_single_phase(&table, &coeffs, param, x, y, i, j).expect("eval");
                let want = grid.val[i][j];
                assert!(
                    (v - want).abs() <= 1e-12 * want.abs().max(1.0),
                    "corner ({i},{j}) {param:?}: {v} vs {want}"
                );
                checked += 1;
            }
        }
    }
    assert!(checked >= 20, "expected a meaningful number of corners");
}
