#!/usr/bin/env python3
"""Dump a fluid's final runtime JSON from the CoolProp 8.0.0 oracle wheel into
data/coolprop-json/ (PLAN.md Phase 3 data source decision).

`get_fluid_param_string(fluid, "JSON")` returns the post-injection document
the reference implementation computes with — surface tension, environmental
data, and ancillaries already folded in by upstream's build pipeline, and
numeric values exact (the embedded CBOR carries binary IEEE doubles; the
nlohmann dump prints shortest round-trip decimals).

The string is written verbatim — no re-encoding — so reruns are
byte-identical.

Usage: tools/golden-gen/.venv/bin/python tools/golden-gen/dump_fluid_json.py Water [...]
"""

import sys
from pathlib import Path

from CoolProp.CoolProp import get_fluid_param_string

OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "coolprop-json"


def main():
    fluids = sys.argv[1:]
    if not fluids:
        sys.exit("usage: dump_fluid_json.py <FluidName> [...]")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for fluid in fluids:
        text = get_fluid_param_string(fluid, "JSON")
        path = OUT_DIR / f"{fluid}.json"
        path.write_text(text)
        print(f"wrote {len(text):8d} chars -> {path.relative_to(Path.cwd())}")


if __name__ == "__main__":
    main()
