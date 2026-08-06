#!/usr/bin/env python3
"""Generate the committed golden fixtures under tests/golden/fixtures/.

Deterministic by construction: fixed record lists, insertion-ordered keys,
Python repr floats (shortest round-trip). Rerunning must produce
byte-identical files (PLAN.md step 0.4).

Fixture record schema (JSON Lines, one object per line):
    backend, fluid, out, name1, val1, name2, val2, expected[, rtol]
`expected` is the value the oracle (CoolProp 8.0.0 wheel) returns for
PropsSI(out, name1, val1, name2, val2, "backend::fluid").
"""

import json
import platform
from pathlib import Path

import CoolProp
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "golden" / "fixtures"

# Alias list transcribed from ParameterInformation() in
# src/DataStructures.cpp @ v8.0.0; the oracle validates every entry below.
PARAM_ALIASES = [
    "D", "H", "M", "S", "U", "C", "O", "G", "V", "L",
    "pcrit", "Pcrit", "Tcrit", "Ttriple", "ptriple", "rhocrit",
    "Tmin", "Tmax", "pmax", "pmin", "molemass", "molarmass", "A", "I",
]

# Names upstream must reject, so the Rust side must reject them too.
NEGATIVE_NAMES = ["dmolar", "DMoLAR", "t", "nonsense", "QT_INPUTS", ""]

PHASE_NAMES = [
    "phase_liquid", "phase_gas", "phase_twophase", "phase_supercritical",
    "phase_supercritical_gas", "phase_supercritical_liquid",
    "phase_critical_point", "phase_unknown", "phase_not_imposed",
]


def dump_parameters():
    """One row per valid upstream parameter index (contiguous from 1)."""
    rows, idx, misses = [], 1, 0
    while misses < 3:
        try:
            short = CP.get_parameter_information(idx, "short")
        except Exception:
            misses += 1
            idx += 1
            continue
        rows.append({
            "index": idx,
            "short": short,
            "io": CP.get_parameter_information(idx, "IO"),
            "units": CP.get_parameter_information(idx, "units"),
            "long": CP.get_parameter_information(idx, "long"),
            "trivial": bool(CP.is_trivial_parameter(idx)),
        })
        misses = 0
        idx += 1
    return rows


def dump_param_names(param_rows):
    """Resolution of every short name, alias, their uppercase forms, and
    known-invalid names. index=None means upstream rejects the name."""
    candidates = []
    for r in param_rows:
        candidates += [r["short"], r["short"].upper()]
    for a in PARAM_ALIASES:
        candidates += [a, a.upper()]
    candidates += NEGATIVE_NAMES
    rows = []
    for name in dict.fromkeys(candidates):
        try:
            rows.append({"name": name, "index": int(CP.get_parameter_index(name))})
        except Exception:
            rows.append({"name": name, "index": None})
    return rows


def dump_phases():
    return [{"name": n, "index": int(CP.get_phase_index(n))} for n in PHASE_NAMES]


def write_jsonl(name, rows):
    (FIXTURES / name).write_text("".join(json.dumps(r) + "\n" for r in rows))
    print(f"wrote {len(rows):4d} records -> {name}")

# Tiny HEOS::Water smoke set proving generator + harness plumbing (PLAN 0.4).
# Phase-specific suites (IF97, HEOS grids, ...) supersede it for coverage.
WATER_SMOKE = [
    ("HEOS", "Water", "T", "P", 101325.0, "Q", 0.0),
    ("HEOS", "Water", "P", "T", 373.15, "Q", 1.0),
    ("HEOS", "Water", "D", "T", 300.0, "P", 101325.0),
    ("HEOS", "Water", "H", "T", 300.0, "P", 101325.0),
    ("HEOS", "Water", "S", "T", 300.0, "P", 101325.0),
    ("HEOS", "Water", "C", "T", 300.0, "P", 101325.0),
    ("HEOS", "Water", "CVMASS", "T", 300.0, "P", 101325.0),
    ("HEOS", "Water", "A", "T", 300.0, "P", 101325.0),
    ("HEOS", "Water", "D", "T", 400.0, "P", 101325.0),
    ("HEOS", "Water", "D", "T", 700.0, "P", 25e6),
    ("HEOS", "Water", "H", "P", 1e6, "Q", 0.5),
    ("HEOS", "Water", "S", "P", 1e6, "Q", 0.5),
]


def record(backend, fluid, out, name1, val1, name2, val2):
    expected = PropsSI(out, name1, val1, name2, val2, f"{backend}::{fluid}")
    return {
        "backend": backend,
        "fluid": fluid,
        "out": out,
        "name1": name1,
        "val1": val1,
        "name2": name2,
        "val2": val2,
        "expected": expected,
    }


def main():
    FIXTURES.mkdir(parents=True, exist_ok=True)
    write_jsonl("water_propssi_smoke.jsonl", [record(*spec) for spec in WATER_SMOKE])
    param_rows = dump_parameters()
    write_jsonl("parameters.jsonl", param_rows)
    write_jsonl("param_aliases.jsonl", dump_param_names(param_rows))
    write_jsonl("phases.jsonl", dump_phases())

    manifest = {
        "generator": "tools/golden-gen/gen_fixtures.py",
        "coolprop_version": CoolProp.__version__,
        "upstream_tag": "v8.0.0",
        "platform": f"{platform.system()}-{platform.machine()}",
        "files": [
            "param_aliases.jsonl",
            "parameters.jsonl",
            "phases.jsonl",
            "water_propssi_smoke.jsonl",
        ],
    }
    (FIXTURES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
