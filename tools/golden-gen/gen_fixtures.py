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


def try_record(out, n1, v1, n2, v2, backend, fluid):
    try:
        return {
            "backend": backend, "fluid": fluid, "out": out,
            "name1": n1, "val1": v1, "name2": n2, "val2": v2,
            "expected": PropsSI(out, n1, v1, n2, v2, f"{backend}::{fluid}"),
        }
    except Exception:
        return None


def gen_if97_water():
    """IF97::Water goldens over every input pair the Phase-2 facade supports
    (PLAN 2.3). Points span regions 1, 2, 3 (incl. near-critical), the
    saturation dome, and the region boundaries."""
    rows = []
    skipped = 0

    # (T [K], p [Pa]) forward states
    pt_points = [
        (280.0, 1e5), (300.0, 101325.0), (300.0, 3e6), (500.0, 3e6),
        (450.0, 1e5), (700.0, 3500.0), (700.0, 30e6), (650.0, 25.5837018e6),
        (647.3, 22.2e6), (750.0, 78.3095639e6), (863.0, 50e6), (1000.0, 10e6),
        (274.0, 1e5), (623.0, 16e6), (640.0, 20.3e6),
    ]
    for t, p in pt_points:
        for out in ["D", "H", "S", "U", "C", "CVMASS", "A", "V", "L"]:
            r = try_record(out, "T", t, "P", p, "IF97", "Water")
            rows.append(r) if r else (skipped := skipped + 1)

    # (p [Pa], Q) dome states
    for p in [101325.0, 1e6, 1e7, 1.6e7, 2e7, 2.2e7]:
        for q in [0.0, 0.5, 1.0]:
            for out in ["T", "D", "H", "S", "U"]:
                r = try_record(out, "P", p, "Q", q, "IF97", "Water")
                rows.append(r) if r else (skipped := skipped + 1)

    # (Q, T) dome states, including surface tension
    for t in [300.0, 373.124, 500.0, 600.0, 623.5, 645.0]:
        for q in [0.0, 1.0]:
            for out in ["P", "D", "H", "I"]:
                r = try_record(out, "Q", q, "T", t, "IF97", "Water")
                rows.append(r) if r else (skipped := skipped + 1)

    # (h, p) backward states: h from representative forward/dome states
    hp_states = [
        (500e3, 3e6), (1500e3, 80e6), (3000e3, 1e5), (3500e3, 5e6),
        (2700e3, 40e6), (2000e3, 50e6), (2100e3, 1e6), (2500e3, 20e6),
    ]
    for h, p in hp_states:
        for out in ["T", "D", "S", "Q"]:
            r = try_record(out, "H", h, "P", p, "IF97", "Water")
            rows.append(r) if r else (skipped := skipped + 1)

    # (p, s) backward states
    ps_states = [
        (3e6, 0.5e3), (80e6, 3e3), (0.1e6, 7.5e3), (8e6, 6e3),
        (20e6, 5.75e3), (20e6, 3.8e3), (1e6, 6.5e3), (50e6, 4.5e3),
    ]
    for p, s in ps_states:
        for out in ["T", "D", "H", "Q"]:
            r = try_record(out, "P", p, "S", s, "IF97", "Water")
            rows.append(r) if r else (skipped := skipped + 1)

    # (h, s) backward states
    hs_states = [
        (90e3, 0.0), (1500e3, 3.4e3), (2800e3, 6.5e3), (3600e3, 7e3),
        (2800e3, 5.1e3), (2100e3, 4.3e3), (2600e3, 5.1e3), (2400e3, 6e3),
    ]
    for h, s in hs_states:
        for out in ["P", "T"]:
            r = try_record(out, "H", h, "S", s, "IF97", "Water")
            rows.append(r) if r else (skipped := skipped + 1)

    # Trivial outputs (state values are ignored by the oracle)
    for out in ["TCRIT", "PCRIT", "RHOCRIT", "TTRIPLE", "PTRIPLE", "TMIN",
                "TMAX", "PMIN", "PMAX", "M", "ACENTRIC"]:
        r = try_record(out, "T", 300.0, "P", 1e5, "IF97", "Water")
        rows.append(r) if r else (skipped := skipped + 1)

    print(f"if97: {len(rows)} records, {skipped} rejected by the oracle")
    return rows


def gen_heos_water_terms():
    """Term-level goldens (PLAN 4.1): alphar/alpha0 and their tau/delta
    derivatives from AbstractState, on a single-phase (T, rhomolar) grid
    spanning liquid, vapor, supercritical, and near-critical states."""
    state = CP.AbstractState("HEOS", "Water")
    points = [
        (280.0, 55500.0), (300.0, 55000.0), (350.0, 54000.0),
        (400.0, 30.0), (500.0, 100.0), (600.0, 500.0),
        (700.0, 20000.0), (650.0, 17873.72799560906), (630.0, 40000.0),
        (647.2, 17500.0), (648.0, 18500.0), (620.0, 5000.0),
    ]
    accessors = [
        "alphar", "dalphar_dDelta", "dalphar_dTau", "d2alphar_dDelta2",
        "d2alphar_dDelta_dTau", "d2alphar_dTau2", "d3alphar_dDelta3",
        "d3alphar_dDelta2_dTau", "d3alphar_dDelta_dTau2", "d3alphar_dTau3",
        "alpha0", "dalpha0_dDelta", "dalpha0_dTau", "d2alpha0_dDelta2",
        "d2alpha0_dDelta_dTau", "d2alpha0_dTau2", "d3alpha0_dDelta3",
        "d3alpha0_dDelta2_dTau", "d3alpha0_dDelta_dTau2", "d3alpha0_dTau3",
    ]
    rows, skipped = [], []
    for (t, rho) in points:
        state.update(CP.DmolarT_INPUTS, rho, t)
        for name in accessors:
            fn = getattr(state, name, None)
            if fn is None:
                skipped.append(name)
                continue
            rows.append({"fluid": "Water", "t": t, "rhomolar": rho,
                         "out": name, "expected": fn()})
    if skipped:
        print(f"heos terms: accessors missing from wheel: {sorted(set(skipped))}")
    return rows


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
    write_jsonl("if97_water.jsonl", gen_if97_water())
    write_jsonl("heos_water_terms.jsonl", gen_heos_water_terms())
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
            "heos_water_terms.jsonl",
            "if97_water.jsonl",
            "param_aliases.jsonl",
            "parameters.jsonl",
            "phases.jsonl",
            "water_propssi_smoke.jsonl",
        ],
    }
    (FIXTURES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
