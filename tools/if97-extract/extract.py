#!/usr/bin/env python3
"""Mechanically extract every coefficient table from the pinned CoolProp/IF97
header (commit 7aaced02, see PLAN.md 2.1) into Rust source.

Output: crates/rustprop-if97/src/tables.rs (fully generated — never edit).
Numeric literals are transferred verbatim, so the Rust f64 values are
bit-identical to the C++ doubles. Entry counts are asserted against the sizes
the upstream classes pass to their constructors.

Usage: python3 tools/if97-extract/extract.py [path-to-IF97-checkout]
"""

import re
import sys
from pathlib import Path

DEFAULT_SRC = Path.home() / "homecloud/dev/IF97"
PINNED_COMMIT = "7aaced024a702f0985474bf293cdaae9c8d06521"
OUT = Path(__file__).resolve().parents[2] / "crates/rustprop-if97/src/tables.rs"

# C++ array name -> (rust name, element kind, expected count)
# kinds: resid (int,int,f64) | ideal (int,f64) | back (f64,f64,f64)
#      | div (int,f64) | satn (int,f64 -> emit n only) | plain (f64 list)
TABLES = {
    "Hresiddata": ("VISC_RESID", "resid", 21),
    "Hidealdata": ("VISC_IDEAL", "ideal", 4),
    "Lresiddata": ("COND_RESID", "resid", 30),
    "Lidealdata": ("COND_IDEAL", "ideal", 5),
    "Region1residdata": ("REGION1_RESID", "resid", 34),
    "Region2residdata": ("REGION2_RESID", "resid", 43),
    "Region2idealdata": ("REGION2_IDEAL", "ideal", 9),
    "Region23data": ("REGION23_N", "plain", 5),
    "Region3residdata": ("REGION3_RESID", "resid", 40),
    "Region5residdata": ("REGION5_RESID", "resid", 6),
    "Region5idealdata": ("REGION5_IDEAL", "ideal", 6),
    "Region2b2cdata": ("REGION2B2C_N", "plain", 5),
    "Region3abdata": ("REGION3AB_N", "plain", 4),
    "Region2abdata": ("REGION2AB_N", "plain", 4),
    "HTmaxdata": ("HTMAX_N", "plain", 4),
    "sat": ("REGION4_N", "satn", 10),
    # Region 3 v(T,p) subregions, SR5-05
    **{f"Region3{c}data": (f"R3_{c}", "resid", n) for c, n in zip(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
        [30, 32, 35, 38, 29, 42, 38, 29, 42, 29, 34, 43, 40, 39, 24, 27, 24,
         27, 29, 33, 38, 39, 35, 36, 20, 23])},
    # Region 3 dividing lines
    **{f"{ab}data": (f"DIV_{ab}", "div", n) for ab, n in [
        ("AB", 5), ("CD", 4), ("GH", 5), ("IJ", 5), ("JK", 5), ("MN", 4),
        ("OP", 5), ("QU", 4), ("RX", 4), ("UV", 4), ("WX", 5)]},
    # Backward T(p,h), T(p,s), p(h,s) and h/s boundary coefficient sets
    **{f"Coeff{s}": (f"COEFF_{s.upper()}", "back", n) for s, n in [
        ("1H", 20), ("1S", 20), ("1HS", 19),
        ("2aH", 34), ("2bH", 38), ("2cH", 23),
        ("2aS", 46), ("2bS", 44), ("2cS", 30),
        ("2aHS", 29), ("2bHS", 33), ("2cHS", 31),
        ("3aH", 31), ("3bH", 33), ("3aS", 33), ("3bS", 28),
        ("3aHS", 33), ("3bHS", 35),
        ("b14HS", 27), ("b3a4HS", 19), ("b2abHS", 30), ("b2c3bHS", 16),
        ("b13HS", 6), ("Tb23HS", 25), ("T4HS", 36)]},
}

STRUCTS = """\
pub struct Resid {
    pub i: i32,
    pub j: i32,
    pub n: f64,
}
pub struct Ideal {
    pub j: i32,
    pub n: f64,
}
pub struct BackResid {
    pub i: f64,
    pub j: f64,
    pub n: f64,
}
pub struct Division {
    pub i: i32,
    pub n: f64,
}
pub struct Table5Row {
    pub region: u8,
    pub t: f64,
    pub p: f64,
    pub v: f64,
}
pub struct Table3Row {
    pub line: &'static str,
    pub p: f64,
    pub t: f64,
}
"""


def f64(tok):
    tok = tok.strip()
    if re.fullmatch(r"[-+]?\d+", tok):
        return tok + ".0"
    return tok


def i32(tok):
    tok = tok.strip()
    assert re.fullmatch(r"[-+]?\d+", tok), f"expected int, got {tok!r}"
    return tok


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def entries(body):
    """Split '{a, b, c}, {d, e, f}, ...' into ['a, b, c', 'd, e, f', ...]."""
    return [m.group(1) for m in re.finditer(r"\{([^{}]*)\}", body)]


def emit(name, kind, body, expected):
    rows = entries(body) if kind != "plain" else None
    out = []
    if kind == "resid":
        out.append(f"pub const {name}: &[Resid] = &[")
        for e in rows:
            i, j, n = e.split(",")
            out.append(f"    Resid {{ i: {i32(i)}, j: {i32(j)}, n: {f64(n)} }},")
    elif kind == "ideal":
        out.append(f"pub const {name}: &[Ideal] = &[")
        for e in rows:
            j, n = e.split(",")
            out.append(f"    Ideal {{ j: {i32(j)}, n: {f64(n)} }},")
    elif kind == "back":
        out.append(f"pub const {name}: &[BackResid] = &[")
        for e in rows:
            i, j, n = e.split(",")
            out.append(f"    BackResid {{ i: {f64(i)}, j: {f64(j)}, n: {f64(n)} }},")
    elif kind == "div":
        out.append(f"pub const {name}: &[Division] = &[")
        for e in rows:
            i, n = e.split(",")
            out.append(f"    Division {{ i: {i32(i)}, n: {f64(n)} }},")
    elif kind == "satn":
        out.append(f"pub const {name}: &[f64] = &[")
        for e in rows:
            i, n = e.split(",")  # index i is implicit (1..=10), keep order only
            out.append(f"    {f64(n)},")
    elif kind == "plain":
        out.append(f"pub const {name}: &[f64] = &[")
        vals = [v for v in body.split(",") if v.strip()]
        rows = vals
        for v in vals:
            out.append(f"    {f64(v)},")
    out.append("];")
    assert len(rows) == expected, f"{name}: {len(rows)} entries, expected {expected}"
    return "\n".join(out)


def main():
    src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SRC
    import subprocess
    head = subprocess.run(["git", "-C", str(src_dir), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    assert head == PINNED_COMMIT, f"IF97 checkout at {head}, expected {PINNED_COMMIT}"
    text = strip_comments((src_dir / "IF97.h").read_text())

    chunks = [
        "//! GENERATED by tools/if97-extract/extract.py — DO NOT EDIT.",
        "//!",
        f"//! Source: CoolProp/IF97 @ {PINNED_COMMIT} (IF97.h),",
        "//! the exact revision CoolProp v8.0.0 pins via CPM. Numeric literals",
        "//! are transferred verbatim, so every f64 here is bit-identical to",
        "//! the corresponding C++ double.",
        "",
        "#![allow(clippy::excessive_precision)] // literals are verbatim upstream",
        "",
        STRUCTS,
    ]

    for cpp_name, (rust_name, kind, expected) in TABLES.items():
        # A named array: `NAME[] = { ... };` (SaturationElement uses `sat[]`)
        m = re.search(re.escape(cpp_name) + r"\s*\[\]\s*=\s*\{(.*?)\};", text, re.S)
        assert m, f"table {cpp_name} not found"
        chunks.append(emit(rust_name, kind, m.group(1), expected))
        chunks.append("")

    # A[6][5] critical-enhancement matrix (row-major as written upstream)
    m = re.search(r"double A\[6\]\[5\] = \{(.*?)\};", text, re.S)
    assert m
    rows = entries(m.group(1))
    assert len(rows) == 6, f"A matrix: {len(rows)} rows"
    a = ["pub const COND_CRIT_A: [[f64; 5]; 6] = ["]
    for r in rows:
        vals = [f64(v) for v in r.split(",") if v.strip()]
        assert len(vals) == 5
        a.append("    [" + ", ".join(vals) + "],")
    a.append("];")
    chunks.append("\n".join(a))
    chunks.append("")

    # Verification tables shipped with the header (ENABLE_CATCH section)
    m = re.search(r"_Table5\[\]\s*=\s*\{(.*?)\};", text, re.S)
    assert m
    rows = entries(m.group(1))
    assert len(rows) == 52, f"Table5: {len(rows)} rows"
    t5 = ["/// SR5-05(2016) Tables 5 & 13: v(T,p) check values per subregion.",
          "pub const TABLE5: &[Table5Row] = &["]
    for r in rows:
        reg, t, p, v = [x.strip() for x in r.split(",")]
        assert re.fullmatch(r"'[A-Z]'", reg)
        t5.append(f"    Table5Row {{ region: b{reg}, t: {f64(t)}, p: {f64(p)}, v: {f64(v)} }},")
    t5.append("];")
    chunks.append("\n".join(t5))
    chunks.append("")

    m = re.search(r"_Table3\[\]\s*=\s*\{(.*?)\};", text, re.S)
    assert m
    rows = entries(m.group(1))
    assert len(rows) == 12, f"Table3: {len(rows)} rows"
    t3 = ["/// SR5-05(2016) Tables 3 & 11: dividing-line T(p) check values.",
          "pub const TABLE3: &[Table3Row] = &["]
    for r in rows:
        line, p, t = [x.strip() for x in r.split(",")]
        line = line.split("LINE_")[-1]
        t3.append(f'    Table3Row {{ line: "{line}", p: {f64(p)}, t: {f64(t)} }},')
    t3.append("];")
    chunks.append("\n".join(t3))
    chunks.append("")

    OUT.write_text("\n".join(chunks))
    n_consts = len(TABLES) + 3
    print(f"wrote {OUT} ({n_consts} tables, all entry counts verified)")


if __name__ == "__main__":
    main()
