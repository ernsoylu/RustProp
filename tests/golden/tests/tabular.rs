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
