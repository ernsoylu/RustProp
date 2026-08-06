//! Golden tests for the rustprop-core parameter system (PLAN.md step 1.1):
//! every index, name, IO class, unit string, description, and trivial flag
//! must match the dump from the CoolProp 8.0.0 wheel.

use rustprop_core::params::{Param, Phase};
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(name)
}

fn load<T: serde::de::DeserializeOwned>(name: &str) -> Vec<T> {
    let text = std::fs::read_to_string(fixture(name))
        .unwrap_or_else(|e| panic!("cannot read fixture {name}: {e}"));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("bad line in {name}: {e}\n{l}")))
        .collect()
}

#[derive(serde::Deserialize)]
struct ParamRow {
    index: i32,
    short: String,
    io: String,
    units: String,
    long: String,
    trivial: bool,
}

#[derive(serde::Deserialize)]
struct NameRow {
    name: String,
    /// `None` means upstream rejects the name.
    index: Option<i32>,
}

#[test]
fn parameter_table_matches_upstream_dump() {
    let rows: Vec<ParamRow> = load("parameters.jsonl");
    assert_eq!(
        rows.len(),
        Param::ALL.len(),
        "parameter count differs from upstream"
    );
    for row in rows {
        let p = Param::from_index(row.index)
            .unwrap_or_else(|| panic!("no Param for upstream index {} ({})", row.index, row.short));
        assert_eq!(
            p.short_name(),
            row.short,
            "short name of index {}",
            row.index
        );
        assert_eq!(p.io(), row.io, "IO class of {}", row.short);
        assert_eq!(p.units(), row.units, "units of {}", row.short);
        assert_eq!(p.long_desc(), row.long, "long description of {}", row.short);
        assert_eq!(p.is_trivial(), row.trivial, "trivial flag of {}", row.short);
    }
}

#[test]
fn parameter_name_resolution_matches_upstream() {
    let rows: Vec<NameRow> = load("param_aliases.jsonl");
    assert!(
        rows.len() > 150,
        "expected shorts+aliases+case variants, got {}",
        rows.len()
    );
    for row in rows {
        let got = Param::parse(&row.name).map(|p| p.index());
        assert_eq!(got, row.index, "resolution of {:?}", row.name);
    }
}

#[test]
fn phase_indices_match_upstream() {
    let rows: Vec<NameRow> = load("phases.jsonl");
    assert_eq!(
        rows.len(),
        Phase::ALL.len(),
        "phase count differs from upstream"
    );
    for row in rows {
        let got = Phase::parse(&row.name).map(|p| p.index());
        assert_eq!(got, row.index, "phase {}", row.name);
    }
}
