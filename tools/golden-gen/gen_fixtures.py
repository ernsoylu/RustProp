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
from CoolProp.CoolProp import PropsSI

FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "golden" / "fixtures"

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
    rows = [record(*spec) for spec in WATER_SMOKE]
    out_path = FIXTURES / "water_propssi_smoke.jsonl"
    out_path.write_text("".join(json.dumps(r) + "\n" for r in rows))

    manifest = {
        "generator": "tools/golden-gen/gen_fixtures.py",
        "coolprop_version": CoolProp.__version__,
        "upstream_tag": "v8.0.0",
        "platform": f"{platform.system()}-{platform.machine()}",
        "files": ["water_propssi_smoke.jsonl"],
    }
    (FIXTURES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {len(rows)} records -> {out_path}")


if __name__ == "__main__":
    main()
