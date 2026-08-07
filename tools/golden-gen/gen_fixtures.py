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


def gen_heos_water_props():
    """Property-level goldens (PLAN 4.2): p/h/s/u/cv/cp/w/g at the same
    single-phase (T, rhomolar) grid as the term goldens."""
    state = CP.AbstractState("HEOS", "Water")
    # Liquid densities sit clearly above the saturated-liquid curve; the
    # term-grid values 55500/55000/54000 are a hair inside the dome.
    points = [
        (280.0, 56000.0), (300.0, 55500.0), (350.0, 54500.0),
        (400.0, 30.0), (500.0, 100.0), (600.0, 500.0),
        (700.0, 20000.0), (650.0, 17873.72799560906), (630.0, 40000.0),
        (647.2, 17500.0), (648.0, 18500.0), (620.0, 5000.0),
    ]
    accessors = ["p", "hmolar", "smolar", "umolar", "cvmolar", "cpmolar",
                 "speed_sound", "gibbsmolar"]
    rows = []
    for (t, rho) in points:
        state.update(CP.DmolarT_INPUTS, rho, t)
        for name in accessors:
            rows.append({"fluid": "Water", "t": t, "rhomolar": rho,
                         "out": name, "expected": getattr(state, name)()})
    return rows


def gen_heos_water_ancillary():
    """Classic ancillary goldens (PLAN 4.3) via CP.saturation_ancillary."""
    temps = [275.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 625.0, 640.0, 645.0]
    rows = []
    for t in temps:
        for (out, q) in [("P", 1), ("Dmolar", 0), ("Dmolar", 1)]:
            rows.append({"fluid": "Water", "t": t, "out": out, "q": q,
                         "expected": CP.saturation_ancillary("Water", out, q, "T", t)})
    return rows


def gen_heos_water_sat():
    """Saturation goldens (PLAN 4.4): QT and PQ states at Q=0/1 across the
    dome, near-critical included."""
    rows, skipped = [], 0
    qt_temps = [274.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0,
                620.0, 640.0, 645.0, 646.5, 647.05]
    for t in qt_temps:
        for q in [0.0, 1.0]:
            for out in ["P", "Dmolar", "Hmolar", "Smolar", "Umolar", "Cpmolar", "A"]:
                r = try_record(out, "T", t, "Q", q, "HEOS", "Water")
                rows.append(r) if r else (skipped := skipped + 1)
    pq_pressures = [700.0, 1e3, 1e4, 101325.0, 1e6, 5e6, 1e7, 1.5e7, 2e7,
                    2.15e7, 2.2e7]
    for p_ in pq_pressures:
        for q in [0.0, 1.0]:
            for out in ["T", "Dmolar", "Hmolar", "Smolar", "Umolar"]:
                r = try_record(out, "P", p_, "Q", q, "HEOS", "Water")
                rows.append(r) if r else (skipped := skipped + 1)
    print(f"heos sat: {len(rows)} records, {skipped} rejected by the oracle")
    return rows


def gen_heos_water_pt():
    """PT goldens (PLAN 4.5): liquid, vapor, and the supercritical phases,
    covering both phase-determination paths (T threshold ~310.5 K for Water)
    and the solver strategy branches."""
    points = [
        # liquid (p-path below ~310.5 K, T-path above)
        (280.0, 1e5), (300.0, 101325.0), (300.0, 1e7), (350.0, 5e6),
        (450.0, 5e6), (550.0, 1e7), (600.0, 1.5e7),
        # vapor
        (300.0, 1000.0), (400.0, 1e5), (500.0, 1e6), (600.0, 5e6), (640.0, 1e7),
        # supercritical (T>Tc, p>pc)
        (650.0, 3e7), (660.0, 2.25e7), (700.0, 5e7),
        # supercritical gas (T>Tc, p<pc)
        (700.0, 1e7), (800.0, 2e7),
        # supercritical liquid (T<Tc, p>pc)
        (640.0, 3e7), (300.0, 5e7), (280.0, 8e7),
    ]
    rows, skipped = [], 0
    for (t, p_) in points:
        for out in ["Dmolar", "Hmolar", "Smolar", "Cpmolar", "A"]:
            r = try_record(out, "T", t, "P", p_, "HEOS", "Water")
            rows.append(r) if r else (skipped := skipped + 1)
    print(f"heos pt: {len(rows)} records, {skipped} rejected by the oracle")
    return rows


def gen_heos_water_flash():
    """Flash-pair goldens (PLAN 4.6): general-quality (T,Q)/(P,Q), (D,T)
    incl. two-phase, and (H,P)/(P,S) incl. two-phase and T outputs."""
    rows, skipped = [], 0
    # general-quality QT
    for t in [300.0, 400.0, 500.0, 600.0, 640.0]:
        for q in [0.25, 0.5, 0.75]:
            for out in ["P", "Dmolar", "Hmolar", "Smolar", "Umolar"]:
                r = try_record(out, "T", t, "Q", q, "HEOS", "Water")
                rows.append(r) if r else (skipped := skipped + 1)
    # general-quality PQ
    for p_ in [1e5, 1e6, 5e6, 1e7, 2e7]:
        for q in [0.3, 0.7]:
            for out in ["T", "Dmolar", "Hmolar", "Smolar", "Umolar"]:
                r = try_record(out, "P", p_, "Q", q, "HEOS", "Water")
                rows.append(r) if r else (skipped := skipped + 1)
    # DmolarT: liquid, vapor, supercritical, and two-phase states
    for (rho, t) in [(55000.0, 320.0), (50.0, 400.0), (20000.0, 700.0),
                     (1000.0, 400.0), (30000.0, 550.0), (17873.0, 640.0),
                     (40000.0, 660.0), (5000.0, 680.0)]:
        for out in ["P", "Hmolar", "Smolar", "Q"]:
            r = try_record(out, "Dmolar", rho, "T", t, "HEOS", "Water")
            rows.append(r) if r else (skipped := skipped + 1)
    # HmolarP: liquid, two-phase, gas, supercritical
    for (h, p_) in [(5e3, 1e6), (30e3, 1e6), (55e3, 1e6), (60e3, 1e5),
                    (10e3, 5e7), (40e3, 3e7), (25e3, 2.5e7), (48e3, 5e6)]:
        for out in ["T", "Dmolar", "Smolar", "Q"]:
            r = try_record(out, "Hmolar", h, "P", p_, "HEOS", "Water")
            rows.append(r) if r else (skipped := skipped + 1)
    # PSmolar
    for (p_, s_) in [(1e6, 20.0), (1e6, 70.0), (1e6, 130.0), (5e7, 10.0),
                     (3e7, 90.0), (1e5, 140.0), (2.5e7, 60.0), (5e6, 110.0)]:
        for out in ["T", "Dmolar", "Hmolar", "Q"]:
            r = try_record(out, "P", p_, "Smolar", s_, "HEOS", "Water")
            rows.append(r) if r else (skipped := skipped + 1)
    # DmolarP: liquid, two-phase, gas, and both supercritical classifications
    for (rho, p_) in [(55000.0, 1e6), (54000.0, 1e7), (5000.0, 1e6), (500.0, 1e5),
                      (17873.0, 2.1e7), (30.0, 1e5), (300.0, 1e6),
                      (20000.0, 3e7), (50000.0, 3e7), (5000.0, 2.5e7)]:
        for out in ["T", "Hmolar", "Smolar", "Q"]:
            r = try_record(out, "Dmolar", rho, "P", p_, "HEOS", "Water")
            rows.append(r) if r else (skipped := skipped + 1)
    print(f"heos flash: {len(rows)} records, {skipped} rejected by the oracle")
    return rows


# ---------------------------------------------------------------------------
# Parameterized HEOS suites (PLAN 4.7): the water suites above stay verbatim
# so their fixtures remain byte-identical; every further fluid runs the same
# six suites on grids derived from its own characteristic points (reduced
# coordinates), queried from the pinned oracle wheel — still deterministic.
# ---------------------------------------------------------------------------

HEOS_FLUIDS = ["Nitrogen", "CarbonDioxide", "R134a", "n-Propane", "Ammonia",
               # 4.8 representative new-term-family fluids:
               "Methanol", "R125", "RC318", "R22", "n-Heptane", "Fluorine"]

TERM_ACCESSORS = [
    "alphar", "dalphar_dDelta", "dalphar_dTau", "d2alphar_dDelta2",
    "d2alphar_dDelta_dTau", "d2alphar_dTau2", "d3alphar_dDelta3",
    "d3alphar_dDelta2_dTau", "d3alphar_dDelta_dTau2", "d3alphar_dTau3",
    "alpha0", "dalpha0_dDelta", "dalpha0_dTau", "d2alpha0_dDelta2",
    "d2alpha0_dDelta_dTau", "d2alpha0_dTau2", "d3alpha0_dDelta3",
    "d3alpha0_dDelta2_dTau", "d3alpha0_dDelta_dTau2", "d3alpha0_dTau3",
]

PROP_ACCESSORS = ["p", "hmolar", "smolar", "umolar", "cvmolar", "cpmolar",
                  "speed_sound", "gibbsmolar"]


def module_name(fluid):
    """Mirror of rustprop-datagen's module_name: `n-Propane` -> `n_propane`."""
    s = "".join(c.lower() if c.isalnum() else "_" for c in fluid)
    return ("_" + s) if s[0].isdigit() else s


def gen_heos_fluid_suites(fluid):
    """All six HEOS suites for one fluid; returns {suite: rows}."""
    hf = f"HEOS::{fluid}"
    Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
    Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
    pc = PropsSI("pcrit", "", 0, "", 0, hf)
    pt = PropsSI("ptriple", "", 0, "", 0, hf)
    rhoc = PropsSI("rhomolar_critical", "", 0, "", 0, hf)

    def TL(x):
        return Tt + x * (Tc - Tt)

    def psat(T):
        return PropsSI("P", "T", T, "Q", 1, hf)

    def rhoL(T):
        return PropsSI("Dmolar", "T", T, "Q", 0, hf)

    def rhoV(T):
        return PropsSI("Dmolar", "T", T, "Q", 1, hf)

    def rho_mix(T, q):
        # Backend two-phase density: inverse-volume mixing of the sat curves.
        return 1.0 / (q / rhoV(T) + (1.0 - q) / rhoL(T))

    def prop(out, n1, v1, n2, v2):
        return PropsSI(out, n1, v1, n2, v2, hf)

    suites = {}

    # -- terms + props: shared single-phase (T, rhomolar) grid ---------------
    grid = ([(TL(x), 1.02 * rhoL(TL(x))) for x in (0.1, 0.3, 0.5, 0.98)]
            + [(TL(x), 0.98 * rhoV(TL(x))) for x in (0.3, 0.5, 0.7, 0.85, 0.98)]
            + [(1.05 * Tc, 1.05 * rhoc), (1.02 * Tc, 0.9 * rhoc),
               (1.2 * Tc, 2.0 * rhoc), (1.5 * Tc, 0.3 * rhoc)])
    state = CP.AbstractState("HEOS", fluid)
    term_rows, prop_rows = [], []
    for (t, rho) in grid:
        state.update(CP.DmolarT_INPUTS, rho, t)
        for name in TERM_ACCESSORS:
            term_rows.append({"fluid": fluid, "t": t, "rhomolar": rho,
                              "out": name, "expected": getattr(state, name)()})
        for name in PROP_ACCESSORS:
            prop_rows.append({"fluid": fluid, "t": t, "rhomolar": rho,
                              "out": name, "expected": getattr(state, name)()})
    suites["terms"] = term_rows
    suites["props"] = prop_rows

    # -- classic ancillaries -------------------------------------------------
    anc_rows = []
    for x in (0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98):
        t = TL(x)
        for (out, q) in [("P", 1), ("Dmolar", 0), ("Dmolar", 1)]:
            anc_rows.append({"fluid": fluid, "t": t, "out": out, "q": q,
                             "expected": CP.saturation_ancillary(fluid, out, q, "T", t)})
    suites["ancillary"] = anc_rows

    # -- saturation (QT / PQ at Q = 0/1) -------------------------------------
    rows, skipped = [], 0
    for x in (0.01, 0.1, 0.25, 0.4, 0.55, 0.7, 0.8, 0.9, 0.95, 0.98, 0.995, 0.999):
        for q in [0.0, 1.0]:
            for out in ["P", "Dmolar", "Hmolar", "Smolar", "Umolar", "Cpmolar", "A"]:
                r = try_record(out, "T", TL(x), "Q", q, "HEOS", fluid)
                rows.append(r) if r else (skipped := skipped + 1)
    for y in (0.05, 0.2, 0.4, 0.6, 0.75, 0.85, 0.92, 0.97, 0.99, 0.995):
        p_ = pt * (pc / pt) ** y
        for q in [0.0, 1.0]:
            for out in ["T", "Dmolar", "Hmolar", "Smolar", "Umolar"]:
                r = try_record(out, "P", p_, "Q", q, "HEOS", fluid)
                rows.append(r) if r else (skipped := skipped + 1)
    print(f"{fluid} sat: {len(rows)} records, {skipped} rejected")
    suites["sat"] = rows

    # -- PT flashes ----------------------------------------------------------
    pts = ([(TL(x), 2.5 * psat(TL(x))) for x in (0.1, 0.3, 0.5, 0.7, 0.9)]
           + [(TL(x), 0.5 * psat(TL(x))) for x in (0.15, 0.35, 0.55, 0.75, 0.9)]
           + [(TL(0.5), 1.2 * pt)]
           + [(1.05 * Tc, 1.5 * pc), (1.02 * Tc, 1.02 * pc),
              (1.2 * Tc, 0.5 * pc), (1.5 * Tc, 0.9 * pc),
              (0.995 * Tc, 2.0 * pc), (TL(0.3), 1.5 * pc), (TL(0.6), 3.0 * pc)])
    rows, skipped = [], 0
    for (t, p_) in pts:
        for out in ["Dmolar", "Hmolar", "Smolar", "Cpmolar", "A"]:
            r = try_record(out, "T", t, "P", p_, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    print(f"{fluid} pt: {len(rows)} records, {skipped} rejected")
    suites["pt"] = rows

    # -- flash pairs ---------------------------------------------------------
    rows, skipped = [], 0
    for x in (0.15, 0.4, 0.65, 0.85, 0.95):
        for q in [0.25, 0.5, 0.75]:
            for out in ["P", "Dmolar", "Hmolar", "Smolar", "Umolar"]:
                r = try_record(out, "T", TL(x), "Q", q, "HEOS", fluid)
                rows.append(r) if r else (skipped := skipped + 1)
    for y in (0.2, 0.5, 0.75, 0.9, 0.97):
        p_ = pt * (pc / pt) ** y
        for q in [0.3, 0.7]:
            for out in ["T", "Dmolar", "Hmolar", "Smolar", "Umolar"]:
                r = try_record(out, "P", p_, "Q", q, "HEOS", fluid)
                rows.append(r) if r else (skipped := skipped + 1)
    dt_points = [(1.03 * rhoL(TL(0.2)), TL(0.2)), (0.5 * rhoV(TL(0.5)), TL(0.5)),
                 (1.5 * rhoc, 1.1 * Tc), (rho_mix(TL(0.5), 0.5), TL(0.5)),
                 (rho_mix(TL(0.85), 0.2), TL(0.85)), (0.3 * rhoc, 1.05 * Tc),
                 (rhoc, TL(0.98)), (0.8 * rhoc, 1.3 * Tc)]
    for (rho, t) in dt_points:
        for out in ["P", "Hmolar", "Smolar", "Q"]:
            r = try_record(out, "Dmolar", rho, "T", t, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    T_liq, T_gas, T_mid = TL(0.2), TL(0.9), TL(0.6)
    p_liq, p_gas, p_mid = 2.5 * psat(T_liq), 0.5 * psat(T_gas), 2.0 * psat(T_mid)
    p_2ph_mid, p_2ph_hi = psat(TL(0.5)), psat(TL(0.9))
    xp_pairs = [
        (("T", T_liq, "P", p_liq), p_liq),
        (("P", p_2ph_mid, "Q", 0.4), p_2ph_mid),
        (("P", p_2ph_hi, "Q", 0.7), p_2ph_hi),
        (("T", T_gas, "P", p_gas), p_gas),
        (("T", 1.2 * Tc, "P", 1.5 * pc), 1.5 * pc),
        (("T", 0.995 * Tc, "P", 2.0 * pc), 2.0 * pc),
        (("T", 1.3 * Tc, "P", 0.5 * pc), 0.5 * pc),
        (("T", T_mid, "P", p_mid), p_mid),
    ]
    for (src, p_) in xp_pairs:
        h = prop("Hmolar", *src)
        for out in ["T", "Dmolar", "Smolar", "Q"]:
            r = try_record(out, "Hmolar", h, "P", p_, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    for (src, p_) in xp_pairs:
        s_ = prop("Smolar", *src)
        for out in ["T", "Dmolar", "Hmolar", "Q"]:
            r = try_record(out, "P", p_, "Smolar", s_, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    dp_pairs = [
        (1.03 * rhoL(TL(0.2)), 2.5 * psat(TL(0.2))),
        (rho_mix(TL(0.5), 0.5), psat(TL(0.5))),
        (rho_mix(TL(0.9), 0.3), psat(TL(0.9))),
        (0.5 * rhoV(TL(0.5)), 0.45 * psat(TL(0.5))),
        (1.5 * rhoc, 2.0 * pc),
        (0.3 * rhoc, 0.6 * pc),
        (1.05 * rhoc, 1.05 * pc),
        (2.0 * rhoc, 1.5 * pc),
    ]
    for (rho, p_) in dp_pairs:
        for out in ["T", "Hmolar", "Smolar", "Q"]:
            r = try_record(out, "Dmolar", rho, "P", p_, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    print(f"{fluid} flash: {len(rows)} records, {skipped} rejected")
    suites["flash"] = rows

    return suites


def gen_heos_fluid_hs(fluid):
    """(Hmolar, Smolar) flash suite (PLAN 4.6 final pair): h/s computed by the
    wheel at reduced-coordinate source states spanning subcooled/cold liquid,
    two-phase (low, mid, high T), gas, and the supercritical classes."""
    hf = f"HEOS::{fluid}"
    Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
    Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
    pc = PropsSI("pcrit", "", 0, "", 0, hf)

    def TL(x):
        return Tt + x * (Tc - Tt)

    def psat(T):
        return PropsSI("P", "T", T, "Q", 1, hf)

    sources = [
        ("T", TL(0.2), "P", 2.5 * psat(TL(0.2))),      # subcooled liquid
        ("T", TL(0.05), "P", 10.0 * psat(TL(0.05))),   # cold compressed liquid
        ("P", psat(TL(0.25)), "Q", 0.15),              # two-phase, low T
        ("P", psat(TL(0.5)), "Q", 0.4),                # two-phase, mid T
        ("P", psat(TL(0.9)), "Q", 0.7),                # two-phase, high T
        ("P", psat(TL(0.97)), "Q", 0.5),               # two-phase, near-critical
        ("T", TL(0.9), "P", 0.5 * psat(TL(0.9))),      # gas near saturation
        ("T", TL(0.7), "P", 0.1 * psat(TL(0.7))),      # superheated gas
        ("T", 1.2 * Tc, "P", 1.5 * pc),                # supercritical
        ("T", TL(0.5), "P", 2.0 * pc),                 # supercritical liquid
        ("T", 1.3 * Tc, "P", 0.5 * pc),                # supercritical gas
        ("T", 1.02 * Tc, "P", 1.02 * pc),              # near-critical single phase
    ]
    # Melting-corner sources (HS cascade leg 4): cold compressed liquid just
    # above the melting curve, including sub-triple T where the curve folds
    # back (water). Only fluids with a melting line contribute; states the
    # oracle rejects are dropped by try_record.
    import CoolProp.CoolProp as CPM
    AS = CPM.AbstractState("HEOS", fluid)
    if AS.has_melting_line():
        for pfac in [3.0, 8.0]:
            p_ = pfac * pc
            try:
                Tm = AS.melting_line(CPM.iT, CPM.iP, p_)
                sources.append(("T", Tm * 1.002, "P", p_))
            except ValueError:
                pass
    rows, skipped = [], 0
    for (n1, v1, n2, v2) in sources:
        try:
            h = PropsSI("Hmolar", n1, v1, n2, v2, hf)
            s = PropsSI("Smolar", n1, v1, n2, v2, hf)
        except ValueError:
            skipped += 4
            continue
        for out in ["T", "Dmolar", "P", "Q"]:
            r = try_record(out, "Hmolar", h, "Smolar", s, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    print(f"{fluid} hs: {len(rows)} records, {skipped} rejected")
    return rows


def gen_heos_all_smoke():
    """Coarse smoke grid over every pure fluid with a superancillary
    (PLAN 4.8): saturation, PT (liquid/gas/supercritical), one each of the
    HP/PS/DT/HS pairs, in reduced coordinates. Runs as an `--ignored` suite."""
    import CoolProp.CoolProp as CPP
    all_fluids = CPP.get_global_param_string("fluids_list").split(",")
    pure = []
    for fl in all_fluids:
        d = json.loads(CPP.get_fluid_param_string(fl, "JSON"))[0]
        if not d["EOS"][0].get("pseudo_pure"):
            pure.append(fl)
    rows, skipped = [], 0
    for fluid in pure:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        pc = PropsSI("pcrit", "", 0, "", 0, hf)

        def TL(x):
            return Tt + x * (Tc - Tt)

        def psat(T):
            return PropsSI("P", "T", T, "Q", 1, hf)

        def rec(out, n1, v1, n2, v2):
            nonlocal skipped
            r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)

        try:
            t_mid = TL(0.5)
            p_liq = 2.5 * psat(TL(0.3))
            p_gas = 0.5 * psat(TL(0.8))
        except Exception:
            skipped += 15
            continue
        # saturation
        rec("P", "T", t_mid, "Q", 0.0)
        rec("Dmolar", "T", t_mid, "Q", 0.0)
        rec("Dmolar", "T", t_mid, "Q", 1.0)
        rec("Hmolar", "T", t_mid, "Q", 0.5)
        rec("Smolar", "T", t_mid, "Q", 0.5)
        # PT: liquid, gas, supercritical
        rec("Dmolar", "T", TL(0.3), "P", p_liq)
        rec("Hmolar", "T", TL(0.3), "P", p_liq)
        rec("Dmolar", "T", TL(0.8), "P", p_gas)
        rec("Smolar", "T", TL(0.8), "P", p_gas)
        rec("Dmolar", "T", 1.1 * Tc, "P", 1.5 * pc)
        # HP / PS at wheel-derived single-phase states
        try:
            h_gas = PropsSI("Hmolar", "T", TL(0.8), "P", p_gas, hf)
            s_liq = PropsSI("Smolar", "T", TL(0.3), "P", p_liq, hf)
            rec("T", "Hmolar", h_gas, "P", p_gas)
            rec("T", "P", p_liq, "Smolar", s_liq)
        except Exception:
            skipped += 2
        # DT two-phase
        try:
            rho_mid = PropsSI("Dmolar", "P", psat(t_mid), "Q", 0.5, hf)
            rec("Q", "Dmolar", rho_mid, "T", t_mid)
        except Exception:
            skipped += 1
        # HS: supercritical single-phase and two-phase
        try:
            h_sc = PropsSI("Hmolar", "T", 1.1 * Tc, "P", 1.5 * pc, hf)
            s_sc = PropsSI("Smolar", "T", 1.1 * Tc, "P", 1.5 * pc, hf)
            rec("T", "Hmolar", h_sc, "Smolar", s_sc)
            h_2p = PropsSI("Hmolar", "T", t_mid, "Q", 0.5, hf)
            s_2p = PropsSI("Smolar", "T", t_mid, "Q", 0.5, hf)
            rec("T", "Hmolar", h_2p, "Smolar", s_2p)
        except Exception:
            skipped += 2
    print(f"all-fluids smoke: {len(rows)} records over {len(pure)} fluids, {skipped} skipped")
    return rows


def gen_mixture_helmholtz():
    """Slice 10c: the mixture Helmholtz assembly (corresponding-states sum +
    excess departure term; GERG Table B5 alpha0) probed through the wheel's
    low-level accessors at fixed (T, Dmolar, x). Composition rides name3/val3
    as x1 (x2 = 1 - x1 computed identically on both sides). Pair coverage:
    each departure kind plus the F == 0 placeholder."""
    pairs = [
        ("Methane", "Ethane"),       # GERG-2008 departure (power + gaussian)
        ("Methane", "Nitrogen"),     # GERG-2008 departure
        ("CarbonDioxide", "Water"),  # Exponential departure (Gernert)
        ("R32", "R125"),             # Exponential departure, Lemmon reducing
        ("Helium", "Argon"),         # Gaussian+Exponential departure
        ("Nitrogen", "Oxygen"),      # F = 0 empty departure
    ]
    outs = [
        "alphar", "dalphar_dTau", "dalphar_dDelta",
        "d2alphar_dTau2", "d2alphar_dDelta2", "d2alphar_dDelta_dTau",
        "alpha0", "dalpha0_dTau", "dalpha0_dDelta",
        "d2alpha0_dTau2", "d2alpha0_dDelta2", "d2alpha0_dDelta_dTau",
    ]
    rows = []
    for f1, f2 in pairs:
        state = CP.AbstractState("HEOS", f"{f1}&{f2}")
        # Impose the phase: mixture phase determination re-solves the
        # density (observed delta = 1 + 6e-13 for Helium&Argon), which would
        # bake a perturbed delta into the fixtures. The grid is supercritical
        # by construction, so imposing is truthful and keeps
        # delta = Dmolar/rhor bitwise-reproducible on the Rust side.
        state.specify_phase(CP.iphase_supercritical)
        for x1 in (0.2, 0.5, 0.8, 1.0):
            state.set_mole_fractions([x1, 1.0 - x1])
            tr = state.T_reducing()
            rhor = state.rhomolar_reducing()
            for tfac, dfac in ((1.2, 1.0), (1.2, 0.05), (1.5, 0.6)):
                try:
                    state.update(CP.DmolarT_INPUTS, dfac * rhor, tfac * tr)
                    vals = [(out, getattr(state, out)()) for out in outs]
                except Exception:
                    continue
                for out, val in vals:
                    rows.append({
                        "backend": "HEOS-MIX", "fluid": f"{f1}&{f2}",
                        "out": out,
                        "name1": "T", "val1": tfac * tr,
                        "name2": "Dmolar", "val2": dfac * rhor,
                        "name3": "x1", "val3": x1,
                        "expected": val,
                    })
    return rows


def gen_mixture_pt():
    """Slice 10d: the mixture PT single-phase flash (SRK-seeded lowest-Gibbs
    root selection) + homogeneous mixture properties, via the wheel's real
    PT update (stability machinery included). Points are placed off the dome
    using the wheel's own bubble/dew pressures; any state the wheel labels
    two-phase is skipped (10f territory). x rides name3/val3 as x1."""
    pairs = [
        ("Methane", "Ethane"),
        ("Methane", "Nitrogen"),
        ("CarbonDioxide", "Water"),
        ("R32", "R125"),
        ("Helium", "Argon"),
        ("Nitrogen", "Oxygen"),
    ]
    outs = [
        ("Dmolar", "rhomolar"), ("Hmolar", "hmolar"), ("Smolar", "smolar"),
        ("Umolar", "umolar"), ("Cpmolar", "cpmolar"), ("Cvmolar", "cvmolar"),
        ("speed_of_sound", "speed_sound"), ("Gmolar", "gibbsmolar"),
    ]
    rows = []
    for f1, f2 in pairs:
        state = CP.AbstractState("HEOS", f"{f1}&{f2}")
        for x1 in (0.25, 0.5, 0.75):
            state.set_mole_fractions([x1, 1.0 - x1])
            tr = state.T_reducing()
            # Candidate (T, p) states: compressed liquid above the bubble
            # pressure, superheated gas below the dew pressure, supercritical.
            cands = []
            for tfac in (0.6, 0.75):
                t = tfac * tr
                try:
                    state.update(CP.QT_INPUTS, 0, t)
                    p_bubble = state.p()
                    cands.append((t, 2.0 * p_bubble))
                except Exception:
                    pass
                try:
                    state.update(CP.QT_INPUTS, 1, t)
                    p_dew = state.p()
                    cands.append((t, 0.5 * p_dew))
                except Exception:
                    pass
            cands.append((1.3 * tr, 2.0e7))
            cands.append((1.3 * tr, 1.0e5))
            for t, p in cands:
                try:
                    state.update(CP.PT_INPUTS, p, t)
                    if state.phase() == CP.iphase_twophase:
                        continue
                    vals = [(out, getattr(state, acc)()) for out, acc in outs]
                except Exception:
                    continue
                for out, val in vals:
                    rows.append({
                        "backend": "HEOS-MIX", "fluid": f"{f1}&{f2}",
                        "out": out,
                        "name1": "T", "val1": t,
                        "name2": "P", "val2": p,
                        "name3": "x1", "val3": x1,
                        "expected": val,
                    })
    return rows


def gen_mixture_vle():
    """Slice 10e: blind QT/PQ mixture flashes (Wilson seed -> successive
    substitution -> newton_raphson_saturation). Outputs are the PropsSI-visible
    bulk quantities; x rides name3/val3 as x1."""
    pairs = [
        ("Methane", "Ethane"),
        ("Methane", "Nitrogen"),
        ("R32", "R125"),
        ("Nitrogen", "Oxygen"),
    ]
    rows = []
    for f1, f2 in pairs:
        state = CP.AbstractState("HEOS", f"{f1}&{f2}")
        for x1 in (0.25, 0.5, 0.75):
            state.set_mole_fractions([x1, 1.0 - x1])
            tr = state.T_reducing()
            # QT states
            for tfac in (0.55, 0.65, 0.75):
                for q in (0.0, 0.3, 1.0):
                    t = tfac * tr
                    try:
                        state.update(CP.QT_INPUTS, q, t)
                        vals = [("P", state.p()), ("Dmolar", state.rhomolar()),
                                ("Hmolar", state.hmolar()), ("Smolar", state.smolar())]
                    except Exception:
                        continue
                    for out, val in vals:
                        rows.append({
                            "backend": "HEOS-MIX", "fluid": f"{f1}&{f2}",
                            "out": out, "name1": "T", "val1": t,
                            "name2": "Q", "val2": q,
                            "name3": "x1", "val3": x1,
                            "expected": val,
                        })
            # PQ states seeded from the bubble pressure at 0.65 Tr
            try:
                state.update(CP.QT_INPUTS, 0, 0.65 * tr)
                p_mid = state.p()
            except Exception:
                continue
            for pfac in (0.6, 1.4):
                for q in (0.0, 0.5, 1.0):
                    p = pfac * p_mid
                    try:
                        state.update(CP.PQ_INPUTS, p, q)
                        vals = [("T", state.T()), ("Dmolar", state.rhomolar()),
                                ("Hmolar", state.hmolar()), ("Smolar", state.smolar())]
                    except Exception:
                        continue
                    for out, val in vals:
                        rows.append({
                            "backend": "HEOS-MIX", "fluid": f"{f1}&{f2}",
                            "out": out, "name1": "P", "val1": p,
                            "name2": "Q", "val2": q,
                            "name3": "x1", "val3": x1,
                            "expected": val,
                        })
    return rows


def gen_mixture_propssi():
    """Mixture PropsSI routing goldens: the real string API end to end —
    trivials (weighted limits, reducing values), PT/QT/PQ states with molar
    and mass bases, input echo, swapped order. Error conditions live in the
    Rust test (variant + message asserted, not oracle-recorded)."""
    rows, skipped = [], 0

    def rec(fluid, out, n1, v1, n2, v2):
        nonlocal skipped
        r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
        rows.append(r) if r else (skipped := skipped + 1)

    for f1, f2, x1 in [("Methane", "Ethane", 0.6), ("Methane", "Ethane", 0.25),
                       ("R32", "R125", 0.7), ("Nitrogen", "Oxygen", 0.79)]:
        mix = f"{f1}[{x1}]&{f2}[{1.0 - x1}]"
        hm = f"HEOS::{mix}"
        for out in ["molemass", "gas_constant", "Tmax", "Tmin", "Ttriple",
                    "pmax", "ptriple", "T_reducing", "rhomolar_reducing"]:
            rec(mix, out, "", 0, "", 0)
        tr = PropsSI("T_reducing", "", 0, "", 0, hm)
        # Supercritical + gas PT states (single-phase on both sides)
        for t, p in [(1.3 * tr, 5e6), (1.3 * tr, 1e5)]:
            for out in ["Dmolar", "Hmolar", "Smolar", "Umolar", "Cpmolar",
                        "Cvmolar", "speed_of_sound", "D", "H", "S", "C", "O"]:
                rec(mix, out, "T", t, "P", p)
        # input echo + swapped order
        rec(mix, "T", "T", 1.3 * tr, "P", 1e5)
        rec(mix, "Dmolar", "P", 1e5, "T", 1.3 * tr)
        # mixture transport (log-linear viscosity / linear conductivity of
        # pure components at the bulk state)
        for t, p in [(1.3 * tr, 5e6), (1.3 * tr, 1e5), (0.7 * tr, 2e5)]:
            rec(mix, "viscosity", "T", t, "P", p)
            rec(mix, "conductivity", "T", t, "P", p)
        # QT / PQ two-phase
        t_sat = 0.7 * tr
        for q in (0.0, 1.0, 0.4):
            for out in ["P", "Dmolar", "Hmolar", "Smolar", "D", "H", "S", "Q"]:
                rec(mix, out, "T", t_sat, "Q", q)
        try:
            p_mid = PropsSI("P", "T", t_sat, "Q", 0, hm)
            for q in (0.0, 1.0, 0.5):
                for out in ["T", "Dmolar", "Hmolar", "Smolar"]:
                    rec(mix, out, "P", 0.8 * p_mid, "Q", q)
            # NOTE: mass-basis INPUT records (Hmass&P, P&Smass, Dmass&T)
            # would route into the sweep flashes deferred with slice 10f —
            # they join the fixture when those land.
        except Exception:
            skipped += 1
    print(f"mixture propssi: {len(rows)} records, {skipped} skipped")
    return rows


def gen_mixture_predefined():
    """Predefined mixtures ("<Name>.mix" registry): trivials + PT states for
    binary/ternary/10-component blends, QT/PQ for the refrigerant binaries,
    exercising N > 2 reducing/VLE machinery."""
    rows, skipped = [], 0

    def rec(fluid, out, n1, v1, n2, v2):
        nonlocal skipped
        r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
        rows.append(r) if r else (skipped := skipped + 1)

    blends = ["R410A.mix", "R407C.mix", "R404A.mix", "Air.mix", "Amarillo.mix",
              "R410A.MIX"]
    for mix in blends:
        hm = f"HEOS::{mix}"
        for out in ["molemass", "T_reducing", "rhomolar_reducing", "Tmax", "Ttriple"]:
            rec(mix, out, "", 0, "", 0)
        tr = PropsSI("T_reducing", "", 0, "", 0, hm)
        for t, p in [(1.3 * tr, 5e6), (1.3 * tr, 1e5), (300.0, 1e5)]:
            for out in ["Dmolar", "Hmolar", "Smolar", "Cpmolar", "speed_of_sound"]:
                rec(mix, out, "T", t, "P", p)
        # Two-phase for the refrigerant blends and Air (ternary VLE)
        if mix != "Amarillo.mix":
            t_sat = 0.75 * tr
            for q in (0.0, 1.0):
                for out in ["P", "Dmolar", "Hmolar", "Smolar"]:
                    rec(mix, out, "T", t_sat, "Q", q)
            try:
                p_mid = PropsSI("P", "T", t_sat, "Q", 0, hm)
                for out in ["T", "Dmolar", "Hmolar"]:
                    rec(mix, out, "P", 0.9 * p_mid, "Q", 0.5)
            except Exception:
                skipped += 1
    print(f"mixture predefined: {len(rows)} records, {skipped} skipped")
    return rows


def gen_mixture_pt_twophase():
    """Slice 10f: PT states INSIDE the two-phase dome (stability test +
    Michelsen split) — Q, lever-rule bulk properties. Pressure placed between
    the bubble and dew pressures at each T."""
    pairs = [
        ("Methane", "Ethane"),
        ("Methane", "Nitrogen"),
        ("R32", "R125"),
        ("Nitrogen", "Oxygen"),
    ]
    rows = []
    for f1, f2 in pairs:
        state = CP.AbstractState("HEOS", f"{f1}&{f2}")
        for x1 in (0.3, 0.6):
            state.set_mole_fractions([x1, 1.0 - x1])
            tr = state.T_reducing()
            for tfac in (0.6, 0.72):
                t = tfac * tr
                try:
                    state.update(CP.QT_INPUTS, 0, t)
                    p_bub = state.p()
                    state.update(CP.QT_INPUTS, 1, t)
                    p_dew = state.p()
                except Exception:
                    continue
                for frac in (0.25, 0.6, 0.9):
                    p = p_dew + frac * (p_bub - p_dew)
                    try:
                        state.update(CP.PT_INPUTS, p, t)
                        if state.phase() != CP.iphase_twophase:
                            continue
                        vals = [("Q", state.Q()), ("Dmolar", state.rhomolar()),
                                ("Hmolar", state.hmolar()), ("Smolar", state.smolar())]
                    except Exception:
                        continue
                    for out, val in vals:
                        rows.append({
                            "backend": "HEOS-MIX", "fluid": f"{f1}&{f2}",
                            "out": out, "name1": "T", "val1": t,
                            "name2": "P", "val2": p,
                            "name3": "x1", "val3": x1,
                            "expected": val,
                        })
    return rows


def gen_mixture_sweep():
    """Slice 10f part 2: the sweep-based mixture flashes (DmolarT/HmolarT/
    SmolarT/TUmolar, HmolarP/PSmolar/PUmolar, DmolarHmolar/DmolarSmolar/
    DmolarUmolar). Targets derived from wheel PT states (single-phase gas/
    liquid and in-dome) so both sides solve the same well-posed problems."""
    pairs = [("Methane", "Ethane"), ("R32", "R125")]
    rows, skipped = [], 0

    def rec(mix, out, n1, v1, n2, v2):
        nonlocal skipped
        r = try_record(out, n1, v1, n2, v2, "HEOS", mix)
        rows.append(r) if r else (skipped := skipped + 1)

    for f1, f2 in pairs:
        for x1 in (0.4, 0.7):
            mix = f"{f1}[{x1}]&{f2}[{1.0 - x1}]"
            hm = f"HEOS::{mix}"
            tr = PropsSI("T_reducing", "", 0, "", 0, hm)
            # Anchor states: superheated gas, compressed liquid, in-dome
            anchors = []
            try:
                p_bub = PropsSI("P", "T", 0.7 * tr, "Q", 0, hm)
                p_dew = PropsSI("P", "T", 0.7 * tr, "Q", 1, hm)
                anchors.append((0.7 * tr, 0.5 * p_dew))          # gas
                anchors.append((0.7 * tr, 2.0 * p_bub))          # liquid
                anchors.append((0.7 * tr, 0.5 * (p_bub + p_dew)))  # two-phase
            except Exception:
                skipped += 1
                continue
            for t, p in anchors:
                try:
                    d = PropsSI("Dmolar", "T", t, "P", p, hm)
                    h = PropsSI("Hmolar", "T", t, "P", p, hm)
                    sm = PropsSI("Smolar", "T", t, "P", p, hm)
                    u = PropsSI("Umolar", "T", t, "P", p, hm)
                except Exception:
                    skipped += 1
                    continue
                # DHSU_T
                rec(mix, "P", "Dmolar", d, "T", t)
                rec(mix, "P", "Hmolar", h, "T", t)
                rec(mix, "P", "Smolar", sm, "T", t)
                rec(mix, "P", "T", t, "Umolar", u)
                rec(mix, "Q", "Dmolar", d, "T", t)
                # HSU_P
                rec(mix, "T", "Hmolar", h, "P", p)
                rec(mix, "T", "P", p, "Smolar", sm)
                rec(mix, "T", "P", p, "Umolar", u)
                # HSU_D
                rec(mix, "T", "Dmolar", d, "Hmolar", h)
                rec(mix, "T", "Dmolar", d, "Smolar", sm)
                rec(mix, "T", "Dmolar", d, "Umolar", u)
                rec(mix, "P", "Dmolar", d, "Smolar", sm)
    print(f"mixture sweep: {len(rows)} records, {skipped} skipped")
    return rows


def gen_pcsaft_terms():
    """Phase 11 slice 11b: PC-SAFT kernel values (alphar, Z via P,
    Hmolar_residual, Smolar_residual, Gmolar_residual) at fixed (Dmolar, T)
    with the phase imposed (the wheel's DmolarT phase determination runs
    fragile TQ flashes). Fluid coverage: plain, polar, associating, water
    (runtime sigma), ion system, and a kij mixture."""
    cases = [
        ("TOLUENE", None, [(9033.114359706229, 320.0), (39.44490805826904, 325.0), (8983.377722763931, 325.0)]),
        ("PROPANE", None, [(13000.0, 250.0), (100.0, 300.0)]),
        ("ACETONE", None, [(12000.0, 300.0), (200.0, 350.0)]),
        ("METHANOL", None, [(24000.0, 300.0), (500.0, 400.0)]),
        ("WATER", None, [(50000.0, 300.0), (1000.0, 400.0)]),
        ("METHANOL&CYCLOHEXANE", [0.3, 0.7], [(9000.0, 320.0), (300.0, 350.0)]),
        ("METHANE&N-BUTANE", [0.6, 0.4], [(15000.0, 200.0), (400.0, 300.0)]),
        ("Na+&Cl-&WATER", [0.0907304774758426, 0.0907304774758426, 0.818539045048315],
         [(50000.0, 298.15)]),
    ]
    outs = ["P", "alphar", "Hmolar_residual", "Smolar_residual", "Gmolar_residual"]
    rows, skipped = [], 0
    for names, x, states in cases:
        AS = CP.AbstractState("PCSAFT", names)
        if x is not None:
            AS.set_mole_fractions(x)
        AS.specify_phase(CP.iphase_liquid)
        for rho, t in states:
            try:
                AS.update(CP.DmolarT_INPUTS, rho, t)
                vals = [("P", AS.p()), ("alphar", AS.alphar()),
                        ("Hmolar_residual", AS.hmolar_residual()),
                        ("Smolar_residual", AS.smolar_residual()),
                        ("Gmolar_residual", AS.gibbsmolar_residual())]
            except Exception:
                skipped += 1
                continue
            x1 = x[0] if x is not None else 1.0
            for out, val in vals:
                r = {
                    "backend": "PCSAFT", "fluid": names, "out": out,
                    "name1": "Dmolar", "val1": rho, "name2": "T", "val2": t,
                    "name3": "x1", "val3": x1,
                    "expected": val,
                }
                if out == "P":
                    # P = Z*kb*T*den amplifies Z's cancellation by ~1/Z
                    # (liquid Z ~ 1e-3 from O(10) terms).
                    r["rtol"] = 1e-9
                rows.append(r)
    print(f"pcsaft terms: {len(rows)} records, {skipped} skipped")
    return rows


def gen_pcsaft_flash():
    """Phase 11 slice 11c: PC-SAFT flashes through the string API — QT/PQ
    saturation (incl. the methanol/cyclohexane azeotropic system and the
    NaCl aqueous electrolyte), PT with phase determination, DmolarT."""
    rows, skipped = [], 0

    def rec(fluid, out, n1, v1, n2, v2):
        nonlocal skipped
        r = try_record(out, n1, v1, n2, v2, "PCSAFT", fluid)
        rows.append(r) if r else (skipped := skipped + 1)

    for f in ["PROPANE", "TOLUENE", "N-BUTANE", "ACETONE", "METHANOL", "WATER"]:
        pf = f"PCSAFT::{f}"
        # A subcritical temperature grid via rough Tc-fraction anchors:
        # probe a QT solve and derive states from it.
        for t in (250.0, 300.0, 350.0):
            try:
                p_sat = PropsSI("P", "T", t, "Q", 0, pf)
            except Exception:
                skipped += 1
                continue
            if p_sat < 10.0:
                skipped += 1
                continue
            rec(f, "P", "T", t, "Q", 0)
            rec(f, "Dmolar", "T", t, "Q", 0)
            rec(f, "P", "T", t, "Q", 1)
            rec(f, "Dmolar", "T", t, "Q", 1)
            rec(f, "T", "P", p_sat, "Q", 0)
            if f != "WATER":
                # WATER PT/DT: upstream's phase determination runs on
                # children whose sigma is still the -1 sentinel (quirk 4)
                # and returns physically wrong densities; the port errors
                # loudly instead (documented deviation).
                rec(f, "Dmolar", "T", t, "P", 2.0 * p_sat)   # liquid
                rec(f, "Dmolar", "T", t, "P", 0.5 * p_sat)   # gas
                rec(f, "Q", "Dmolar", 100.0, "T", t)
                rec(f, "Hmolar_residual", "T", t, "P", 2.0 * p_sat)
                rec(f, "Smolar_residual", "T", t, "P", 2.0 * p_sat)
    # Mixtures
    rec("METHANOL[0.3]&CYCLOHEXANE[0.7]", "P", "T", 327.48, "Q", 0)
    rec("METHANOL[0.3]&CYCLOHEXANE[0.7]", "Dmolar", "T", 327.48, "Q", 0)
    rec("METHANE[0.5]&N-BUTANE[0.5]", "P", "T", 250.0, "Q", 0)
    rec("Na+[0.0907304774758426]&Cl-[0.0907304774758426]&WATER[0.818539045048315]",
        "P", "T", 298.15, "Q", 0)
    print(f"pcsaft flash: {len(rows)} records, {skipped} skipped")
    return rows


def gen_partial_derivs():
    """Phase 12 slice 12a: the generic (T,rho)-basis partial derivatives —
    upstream first_partial_deriv / second_partial_deriv, which the Tabular
    table build consumes. Probed through PropsSI derivative strings at
    single-phase states. name3/val3 carry nothing here."""
    rows, skipped = [], 0
    # (Of, Wrt, Constant) triples the LogPH / LogPT table builds use, plus
    # the common thermodynamic ones.
    firsts = [
        ("T", "Hmolar", "P"), ("T", "P", "Hmolar"),
        ("P", "Hmolar", "P"), ("Dmolar", "Hmolar", "P"), ("Dmolar", "P", "Hmolar"),
        ("Hmolar", "T", "P"), ("Hmolar", "P", "T"),
        ("Smolar", "T", "P"), ("Smolar", "P", "T"),
        ("Umolar", "T", "P"), ("Umolar", "P", "T"),
        ("Dmolar", "T", "P"), ("Dmolar", "P", "T"),
        ("P", "T", "Dmolar"), ("P", "Dmolar", "T"),
    ]
    seconds = [
        ("T", "Hmolar", "P", "Hmolar", "P"),
        ("T", "Hmolar", "P", "P", "Hmolar"),
        ("T", "P", "Hmolar", "P", "Hmolar"),
        ("Dmolar", "Hmolar", "P", "Hmolar", "P"),
        ("Dmolar", "P", "Hmolar", "P", "Hmolar"),
        ("Smolar", "T", "P", "T", "P"),
        ("Hmolar", "T", "P", "T", "P"),
        ("P", "Dmolar", "T", "Dmolar", "T"),
    ]
    for fluid in ["Water", "n-Propane", "CarbonDioxide"]:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        t_liq = Tt + 0.4 * (Tc - Tt)
        p_liq = 3.0 * PropsSI("P", "T", t_liq, "Q", 0, hf)
        t_gas = Tt + 0.85 * (Tc - Tt)
        p_gas = 0.4 * PropsSI("P", "T", t_gas, "Q", 1, hf)
        for (t, p) in [(t_liq, p_liq), (t_gas, p_gas), (1.3 * Tc, 5e6)]:
            for (of, wrt, const) in firsts:
                out = f"d({of})/d({wrt})|{const}"
                r = try_record(out, "T", t, "P", p, "HEOS", fluid)
                rows.append(r) if r else (skipped := skipped + 1)
            for (of, w1, c1, w2, c2) in seconds:
                out = f"d(d({of})/d({w1})|{c1})/d({w2})|{c2}"
                r = try_record(out, "T", t, "P", p, "HEOS", fluid)
                rows.append(r) if r else (skipped := skipped + 1)
    print(f"partial derivs: {len(rows)} records, {skipped} skipped")
    return rows


def gen_tabular_tables():
    """Phase 12 slice 12b: the table GRIDS — limits computed exactly as
    upstream's LogPHTable/LogPTTable::set_limits, plus node values at
    sampled grid coordinates (which are plain source-backend evaluations).
    out names: "xmin"/"xmax"/"ymin"/"ymax" carry the limits with the grid
    kind in name1; node records carry the node's (x, y) in val1/val2."""
    rows = []
    for fluid in ["Water", "n-Propane", "CarbonDioxide"]:
        hf = f"HEOS::{fluid}"
        AS = CP.AbstractState("HEOS", fluid)
        Tmin = max(PropsSI("Ttriple", "", 0, "", 0, hf), PropsSI("Tmin", "", 0, "", 0, hf))
        Tmax = PropsSI("Tmax", "", 0, "", 0, hf)
        pmax = PropsSI("pmax", "", 0, "", 0, hf)
        AS.update(CP.QT_INPUTS, 0, Tmin)
        h_satL, p_triple = AS.hmolar(), AS.p()

        # LogPH limits
        AS.update(CP.DmolarT_INPUTS, 1e-10, 1.499 * Tmax)
        xmax1 = AS.hmolar()
        AS.update(CP.PT_INPUTS, pmax, 1.499 * Tmax)
        xmax2 = AS.hmolar()
        ph = {"xmin": h_satL, "xmax": max(xmax1, xmax2), "ymin": p_triple, "ymax": pmax}
        pt = {"xmin": Tmin, "xmax": Tmax * 1.499, "ymin": p_triple, "ymax": pmax}
        for kind, lim in (("LogPH", ph), ("LogPT", pt)):
            for k, v in lim.items():
                rows.append({
                    "backend": "TABULAR", "fluid": fluid, "out": k,
                    "name1": kind, "val1": 0.0, "name2": "", "val2": 0.0,
                    "expected": v,
                })
        # Node values on a coarse 20x20 LogPT grid (x linear in T, y log in p)
        nx = ny = 20
        for i in (0, 7, 13, 19):
            x = pt["xmin"] + (pt["xmax"] - pt["xmin"]) / (nx - 1) * i
            for j in (0, 9, 19):
                import math
                y = math.exp(math.log(pt["ymin"]) + math.log(pt["ymax"] / pt["ymin"]) / (ny - 1) * j)
                try:
                    # Nodes sitting ON the saturation line are decided by the
                    # +-100*DBL_EPSILON band around Tsat(p) — i.e. by ulp-level
                    # inversion noise, not by the algorithm. Both codes carry
                    # the identical band; skip the knife edge.
                    AS.update(CP.PQ_INPUTS, y, 0)
                    if abs(AS.T() - x) / x < 1e-12:
                        continue
                except Exception:
                    pass
                try:
                    AS.update(CP.PT_INPUTS, y, x)
                    if 0.0 <= AS.Q() <= 1.0:
                        continue
                    vals = [("Dmolar", AS.rhomolar()), ("Hmolar", AS.hmolar()),
                            ("Smolar", AS.smolar()), ("Umolar", AS.umolar())]
                except Exception:
                    continue
                for out, v in vals:
                    rows.append({
                        "backend": "TABULAR", "fluid": fluid, "out": out,
                        "name1": "LogPT_node", "val1": x, "name2": "P", "val2": y,
                        "expected": v,
                    })
    print(f"tabular tables: {len(rows)} records")
    return rows


def gen_ttse():
    """Phase 12 slice 12c: TTSE evaluation on the LogPT table through
    PT inputs, against the wheel's own TTSE backend (which builds the same
    200x200 grid). Node records use the exact grid coordinates so the
    expansion must reproduce the stored node value."""
    import math
    rows = []
    for fluid in ["Water", "n-Propane"]:
        AS = CP.AbstractState("TTSE&HEOS", fluid)
        hf = f"HEOS::{fluid}"
        Tmin = max(PropsSI("Ttriple", "", 0, "", 0, hf), PropsSI("Tmin", "", 0, "", 0, hf))
        Tmax = PropsSI("Tmax", "", 0, "", 0, hf)
        pmax = PropsSI("pmax", "", 0, "", 0, hf)
        ref = CP.AbstractState("HEOS", fluid)
        ref.update(CP.QT_INPUTS, 0, Tmin)
        pmin = ref.p()
        xmin, xmax = Tmin, Tmax * 1.499
        # Interior states, plus states landing exactly on grid nodes.
        nx = ny = 200
        states = []
        for fi in (0.13, 0.37, 0.61, 0.88):
            for fj in (0.2, 0.5, 0.8):
                states.append((xmin + (xmax - xmin) * fi,
                               math.exp(math.log(pmin) + math.log(pmax / pmin) * fj)))
        for i in (40, 111, 175):
            for j in (30, 140):
                states.append((xmin + (xmax - xmin) / (nx - 1) * i,
                               math.exp(math.log(pmin) + math.log(pmax / pmin) / (ny - 1) * j)))
        for (t, p) in states:
            try:
                AS.update(CP.PT_INPUTS, p, t)
                vals = [("Dmolar", AS.rhomolar()), ("Hmolar", AS.hmolar()),
                        ("Smolar", AS.smolar()), ("Umolar", AS.umolar())]
            except Exception:
                continue
            for out, v in vals:
                rows.append({
                    "backend": "TTSE", "fluid": fluid, "out": out,
                    "name1": "T", "val1": t, "name2": "P", "val2": p,
                    "expected": v,
                })
    print(f"ttse: {len(rows)} records")
    return rows


def gen_fluid_resolution():
    """5.1 registry goldens: for every pure fluid, how the wheel resolves its
    canonical name, CAS, aliases, and upper(aliases); plus negatives. The
    canonical answer is INFO.NAME of the resolved document."""
    import CoolProp.CoolProp as CPP
    all_fluids = CPP.get_global_param_string("fluids_list").split(",")
    pure, pseudo = [], set()
    for fl in all_fluids:
        d = json.loads(CPP.get_fluid_param_string(fl, "JSON"))[0]
        (pseudo.add(fl) if d["EOS"][0].get("pseudo_pure") else pure.append(fl))
    rows = []
    seen = set()

    def resolve(q):
        try:
            return json.loads(CPP.get_fluid_param_string(q, "JSON"))[0]["INFO"]["NAME"]
        except Exception:
            return None

    for fl in pure:
        d = json.loads(CPP.get_fluid_param_string(fl, "JSON"))[0]
        queries = [fl, d["INFO"]["CAS"]]
        for a in d["INFO"]["ALIASES"]:
            queries += [a, a.upper()]
        for q in queries:
            if q in seen:
                continue
            seen.add(q)
            name = resolve(q)
            if name in pseudo:
                # A pure fluid's string claimed by a pseudo-pure fluid (not
                # in the ported registry) — record for visibility.
                print(f"  NOTE: query {q!r} resolves to pseudo-pure {name!r}")
            rows.append({"query": q, "name": name})
    for q in ["Watr", "", "R134A", "H2O!"]:
        if q not in seen:
            rows.append({"query": q, "name": resolve(q)})
    print(f"fluid resolution: {len(rows)} queries")
    return rows


def gen_props_si():
    """5.2 string-API goldens: mass-basis aliases and inputs, molar forms,
    swapped pair order, input echo, and trivial outputs (empty-name and
    state-input forms), over four fluids."""
    rows, skipped = [], 0

    def rec(fluid, out, n1, v1, n2, v2):
        nonlocal skipped
        r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
        rows.append(r) if r else (skipped := skipped + 1)

    for fluid in ["Water", "CarbonDioxide", "R134a", "Ammonia"]:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)

        def TL(x):
            return Tt + x * (Tc - Tt)

        t_liq = TL(0.3)
        p_liq = 2.5 * PropsSI("P", "T", t_liq, "Q", 1, hf)
        t_gas = TL(0.8)
        p_gas = 0.5 * PropsSI("P", "T", t_gas, "Q", 1, hf)
        # mass-basis output aliases at PT states
        for (t, p) in [(t_liq, p_liq), (t_gas, p_gas)]:
            for out in ["D", "H", "S", "U", "C", "O", "A", "G",
                        "Dmolar", "Hmolar", "Cvmolar"]:
                rec(fluid, out, "T", t, "P", p)
        # swapped order + input echo
        rec(fluid, "Dmolar", "P", p_gas, "T", t_gas)
        rec(fluid, "T", "T", t_gas, "P", p_gas)
        rec(fluid, "P", "T", t_gas, "P", p_gas)
        # mass-basis inputs (values from the wheel at the same states)
        d_mass = PropsSI("D", "T", t_gas, "P", p_gas, hf)
        h_mass = PropsSI("H", "T", t_gas, "P", p_gas, hf)
        s_mass = PropsSI("S", "T", t_gas, "P", p_gas, hf)
        s_liq_mass = PropsSI("S", "T", t_liq, "P", p_liq, hf)
        rec(fluid, "P", "Dmass", d_mass, "T", t_gas)
        rec(fluid, "T", "Hmass", h_mass, "P", p_gas)
        rec(fluid, "T", "P", p_liq, "Smass", s_liq_mass)
        rec(fluid, "T", "Dmass", d_mass, "P", p_gas)
        rec(fluid, "T", "Hmass", h_mass, "Smass", s_mass)
        rec(fluid, "Dmolar", "Smass", s_mass, "T", t_gas)
        # Q pairs incl. mass-basis outputs
        rec(fluid, "H", "T", TL(0.5), "Q", 0.5)
        rec(fluid, "D", "P", PropsSI("P", "T", TL(0.5), "Q", 1, hf), "Q", 0.3)
        rec(fluid, "Q", "T", TL(0.5), "Q", 0.25)
        # trivial outputs: empty-name form and state-input form
        for out in ["Tcrit", "pcrit", "rhocrit", "rhomolar_critical",
                    "Ttriple", "ptriple", "Tmin", "Tmax", "pmax", "M",
                    "acentric", "gas_constant", "T_reducing"]:
            rec(fluid, out, "", 0.0, "", 0.0)
        rec(fluid, "Tcrit", "T", t_gas, "P", p_gas)
        rec(fluid, "M", "T", t_liq, "P", p_liq)
    print(f"props_si: {len(rows)} records, {skipped} rejected")
    return rows


def gen_surface_tension():
    """6.2 goldens: surface tension `I` along the saturation curve for every
    pure fluid with a curve (the wheel errors for fluids without one —
    try_record filters them, and the Rust side asserts the error condition
    separately)."""
    import CoolProp.CoolProp as CPP
    all_fluids = CPP.get_global_param_string("fluids_list").split(",")
    pure = []
    for fl in all_fluids:
        d = json.loads(CPP.get_fluid_param_string(fl, "JSON"))[0]
        if not d["EOS"][0].get("pseudo_pure"):
            pure.append(fl)
    rows, skipped = [], 0
    for fluid in pure:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        for x in (0.05, 0.3, 0.6, 0.9, 0.99):
            t = Tt + x * (Tc - Tt)
            r = try_record("I", "T", t, "Q", 0.5, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)
    print(f"surface tension: {len(rows)} records, {skipped} rejected (incl. curveless fluids)")
    return rows


def gen_flash_pairs_extra():
    """Tier-2 input pairs (PLAN 4.6 deferrals): (Hmolar,T), (T,Umolar) —
    upstream DHSU_T_flash — and (P,Umolar) — upstream HSU_P_flash — across
    liquid/gas/supercritical/two-phase states, plus mass-basis variants."""
    rows, skipped = [], 0
    for fluid in HEOS_FLUIDS:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        pc = PropsSI("pcrit", "", 0, "", 0, hf)

        def TL(x):
            return Tt + x * (Tc - Tt)

        def rec(out, n1, v1, n2, v2):
            nonlocal skipped
            r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)

        # State inventory: (T, p) tuples for the single-phase probes.
        T_liq = TL(0.4)
        p_liq = 2.0 * PropsSI("P", "T", T_liq, "Q", 0, hf)
        T_gas = TL(0.7)
        p_gas = 0.5 * PropsSI("P", "T", T_gas, "Q", 1, hf)
        T_sc = 1.2 * Tc
        p_sc = 1.5 * pc
        for (T, p_) in [(T_liq, p_liq), (T_gas, p_gas), (T_sc, p_sc)]:
            h = PropsSI("Hmolar", "T", T, "P", p_, hf)
            u = PropsSI("Umolar", "T", T, "P", p_, hf)
            for out in ["Dmolar", "P"]:
                rec(out, "Hmolar", h, "T", T)
                rec(out, "T", T, "Umolar", u)
            for out in ["T", "Dmolar"]:
                rec(out, "P", p_, "Umolar", u)
        # Two-phase state at Q=0.3 (DHSU_T two-phase branch; upstream's
        # HSU_P two-phase branch covered via (P,U) at the same state).
        T2 = TL(0.5)
        h2 = PropsSI("Hmolar", "T", T2, "Q", 0.3, hf)
        u2 = PropsSI("Umolar", "T", T2, "Q", 0.3, hf)
        p2 = PropsSI("P", "T", T2, "Q", 0.3, hf)
        for out in ["Dmolar", "Q", "P"]:
            rec(out, "Hmolar", h2, "T", T2)
            rec(out, "T", T2, "Umolar", u2)
        for out in ["T", "Dmolar", "Q"]:
            rec(out, "P", p2, "Umolar", u2)
        # (D, H/S/U) — upstream HSU_D_flash's superancillary happy path —
        # at the same liquid/gas/supercritical/two-phase states.
        for (T, p_) in [(T_liq, p_liq), (T_gas, p_gas), (T_sc, p_sc)]:
            d = PropsSI("Dmolar", "T", T, "P", p_, hf)
            for k in ["Hmolar", "Smolar", "Umolar"]:
                x = PropsSI(k, "T", T, "P", p_, hf)
                for out in ["T", "P"]:
                    rec(out, "Dmolar", d, k, x)
        d2 = PropsSI("Dmolar", "T", T2, "Q", 0.3, hf)
        for k in ["Hmolar", "Smolar", "Umolar"]:
            x = PropsSI(k, "T", T2, "Q", 0.3, hf)
            for out in ["T", "Q", "P"]:
                rec(out, "Dmolar", d2, k, x)
        # (D, Q) — upstream DQ_flash: superancillary root resolution at the
        # Q boundaries and the Brent fallback for fractional quality.
        dl = PropsSI("Dmolar", "T", T2, "Q", 0.0, hf)
        dv = PropsSI("Dmolar", "T", T2, "Q", 1.0, hf)
        dmix = PropsSI("Dmolar", "T", T2, "Q", 0.4, hf)
        for out in ["T", "P"]:
            rec(out, "Dmolar", dl, "Q", 0.0)
            rec(out, "Dmolar", dv, "Q", 1.0)
            rec(out, "Dmolar", dmix, "Q", 0.4)
        rec("T", "Dmass", dl * PropsSI("molemass", "", 0, "", 0, hf), "Q", 0.0)
        # Mass-basis variants exercise mass_to_molar_inputs. (Hmass,T) and
        # (T,Umass) are upstream dead ends — "not yet supported" — asserted
        # Rust-side instead of recorded.
        um = PropsSI("Umass", "T", T_gas, "P", p_gas, hf)
        dm = PropsSI("Dmass", "T", T_liq, "P", p_liq, hf)
        sm = PropsSI("Smass", "T", T_liq, "P", p_liq, hf)
        rec("T", "P", p_gas, "Umass", um)
        rec("T", "Dmass", dm, "Smass", sm)
        # (P, X) below the triple-point pressure: gas states (no saturation).
        p_sub = 0.65 * PropsSI("ptriple", "", 0, "", 0, hf)
        T_sub = 1.05 * Tt
        try:
            h_sub = PropsSI("Hmolar", "T", T_sub, "P", p_sub, hf)
            s_sub = PropsSI("Smolar", "T", T_sub, "P", p_sub, hf)
            u_sub = PropsSI("Umolar", "T", T_sub, "P", p_sub, hf)
            rec("T", "Hmolar", h_sub, "P", p_sub)
            rec("T", "P", p_sub, "Smolar", s_sub)
            rec("T", "P", p_sub, "Umolar", u_sub)
        except ValueError:
            skipped += 3
    print(f"flash pairs extra: {len(rows)} records, {skipped} rejected")
    return rows


def gen_melting():
    """Melting-line goldens (tier-2 deferral): T(p) and p(T) via the wheel's
    `AbstractState.melting_line` for every pure fluid with a curve, plus the
    aggregate limits."""
    import CoolProp.CoolProp as CPCP
    rows, skipped = [], 0
    all_fluids = CPCP.get_global_param_string("fluids_list").split(",")
    pure = []
    for fl in all_fluids:
        d = json.loads(CPCP.get_fluid_param_string(fl, "JSON"))[0]
        if d["EOS"][0].get("pseudo_pure", False):
            continue
        pure.append(fl)
    for fluid in sorted(pure):
        AS = CPCP.AbstractState("HEOS", fluid)
        if not AS.has_melting_line():
            continue
        pmin = AS.melting_line(CPCP.iP_min, -1, -1)
        pmax = AS.melting_line(CPCP.iP_max, -1, -1)
        tmin = AS.melting_line(CPCP.iT_min, -1, -1)
        tmax = AS.melting_line(CPCP.iT_max, -1, -1)
        def rec(out, n1, v1, expected):
            rows.append({"backend": "HEOS", "fluid": fluid, "out": out,
                         "name1": n1, "val1": v1, "name2": "", "val2": 0.0,
                         "expected": expected})
        rec("melt_pmin", "", 0.0, pmin)
        rec("melt_pmax", "", 0.0, pmax)
        # p(T) and T(p) on a grid interior to the fit range
        for x in [0.05, 0.3, 0.6, 0.9]:
            T = tmin + x * (tmax - tmin)
            try:
                rec("melt_p", "T", T, AS.melting_line(CPCP.iP, CPCP.iT, T))
            except ValueError:
                skipped += 1
            p_ = pmin + x * (pmax - pmin)
            try:
                rec("melt_T", "P", p_, AS.melting_line(CPCP.iT, CPCP.iP, p_))
            except ValueError:
                skipped += 1
    print(f"melting: {len(rows)} records, {skipped} rejected")
    return rows


CUBIC_FLUIDS = ["n-Propane", "Water", "CarbonDioxide", "Nitrogen", "R134a",
                "Methane", "Ammonia", "n-Decane", "Benzene", "R32",
                "Isopentane", "MD2M"]


def gen_cubics():
    """Cubic-backend goldens (PLAN 7.1): SRK:: and PR:: PT states
    (liquid/gas/supercritical), QT at Q=0/1, PQ across quality, and the
    trivial outputs, for a structurally diverse fluid set."""
    rows, skipped = [], 0
    for be in ["SRK", "PR"]:
        for fluid in CUBIC_FLUIDS:
            bf = f"{be}::{fluid}"
            Tc = PropsSI("Tcrit", "", 0, "", 0, bf)
            pc = PropsSI("pcrit", "", 0, "", 0, bf)

            def rec(out, n1, v1, n2, v2):
                nonlocal skipped
                r = try_record(out, n1, v1, n2, v2, be, fluid)
                rows.append(r) if r else (skipped := skipped + 1)

            # Trivials
            for out in ["Tcrit", "pcrit", "acentric", "rhomolar_critical",
                        "molemass", "gas_constant"]:
                rec(out, "", 0.0, "", 0.0)
            # QT at the dome
            for x in [0.55, 0.75, 0.9]:
                T = x * Tc
                for q in [0.0, 1.0]:
                    for out in ["P", "Dmolar", "Hmolar", "Smolar"]:
                        rec(out, "T", T, "Q", q)
            # PQ across quality
            psat_mid = PropsSI("P", "T", 0.7 * Tc, "Q", 0, bf)
            for q in [0.0, 0.35, 1.0]:
                for out in ["T", "Dmolar", "Hmolar", "Umolar"]:
                    rec(out, "P", psat_mid, "Q", q)
            # DmolarT (7.2, superancillary route): liquid, gas,
            # supercritical, and in-dome states — upstream's two-phase (D,T)
            # caloric reads throw (broken sub-states), so try_record keeps
            # P/Q there and drops H/S automatically.
            rhoc_kaz = PropsSI("rhomolar_critical", "", 0, "", 0, bf)
            for (rho, T) in [(2.2 * rhoc_kaz, 0.62 * Tc), (0.02 * rhoc_kaz, 0.8 * Tc),
                             (1.5 * rhoc_kaz, 1.3 * Tc), (0.8 * rhoc_kaz, 0.7 * Tc)]:
                for out in ["P", "Q", "Hmolar", "Smolar"]:
                    rec(out, "Dmolar", rho, "T", T)
            # PT: liquid, gas, supercritical(s)
            psat_l = PropsSI("P", "T", 0.6 * Tc, "Q", 0, bf)
            for (T, p_) in [(0.6 * Tc, 3.0 * psat_l), (0.85 * Tc, 0.3 * psat_l),
                            (1.2 * Tc, 1.5 * pc), (1.2 * Tc, 0.5 * pc),
                            (0.7 * Tc, 2.0 * pc)]:
                for out in ["Dmolar", "Hmolar", "Smolar", "Cpmolar", "A"]:
                    rec(out, "T", T, "P", p_)
    print(f"cubics: {len(rows)} records, {skipped} rejected")
    return rows


def gen_cubic_superanc():
    """Cubic-superancillary curve goldens (PLAN 7.2): p/rhoL/rhoV from the
    wheel's `update_QT_pure_superanc` across the dome for both backends."""
    import CoolProp.CoolProp as CPC
    rows = []
    for be in ["SRK", "PR"]:
        for fluid in CUBIC_FLUIDS:
            AS = CPC.AbstractState(be, fluid)
            Tc = AS.T_critical()
            for x in [0.35, 0.5, 0.7, 0.85, 0.95, 0.99]:
                T = x * Tc
                try:
                    AS.update_QT_pure_superanc(0.0, T)
                except Exception:
                    continue
                rows.append({"backend": be, "fluid": fluid, "out": "sa_p",
                             "name1": "T", "val1": T, "name2": "", "val2": 0.0,
                             "expected": AS.p()})
                rows.append({"backend": be, "fluid": fluid, "out": "sa_rhoL",
                             "name1": "T", "val1": T, "name2": "", "val2": 0.0,
                             "expected": AS.saturated_liquid_keyed_output(CPC.iDmolar)})
                rows.append({"backend": be, "fluid": fluid, "out": "sa_rhoV",
                             "name1": "T", "val1": T, "name2": "", "val2": 0.0,
                             "expected": AS.saturated_vapor_keyed_output(CPC.iDmolar)})
    print(f"cubic superanc: {len(rows)} records")
    return rows


INCOMP_PURE = ["DowQ", "DowJ", "Water", "TVP1", "T72", "HC10", "AS10",
               "PMS1", "NaK", "DEB", "FoodWater", "TCO"]
INCOMP_SOLUTIONS = ["MEG[0.3]", "MEG[0.5]", "MPG[0.4]", "LiBr[0.3]", "MAM[0.2]",
                    "MKC[0.3]", "AEG[0.25]", "VNA[0.4]", "MNA[0.2]", "ZM[0.5]"]


def gen_incompressible():
    """Incompressible goldens (PLAN 8.1): PT states across the range plus
    the DmassP/HmassP/PSmass back-flashes and QT (Q=0, psat) for pure and
    concentration-bearing fluids; trivials; V/L where defined."""
    rows, skipped = [], 0
    for fluid in INCOMP_PURE + INCOMP_SOLUTIONS:
        bf = f"INCOMP::{fluid}"
        tmin = PropsSI("Tmin", "", 0, "", 0, bf)
        tmax = PropsSI("Tmax", "", 0, "", 0, bf)

        def TL(f_):
            return tmin + f_ * (tmax - tmin)

        def rec(out, n1, v1, n2, v2):
            nonlocal skipped
            r = try_record(out, n1, v1, n2, v2, "INCOMP", fluid)
            rows.append(r) if r else (skipped := skipped + 1)

        rec("Tmin", "", 0.0, "", 0.0)
        rec("Tmax", "", 0.0, "", 0.0)
        rec("T_freeze", "", 0.0, "", 0.0)
        for x in [0.15, 0.5, 0.85]:
            T = TL(x)
            for p_ in [1e5, 5e5]:
                for out in ["D", "H", "S", "U", "C", "V", "L", "Prandtl"]:
                    rec(out, "T", T, "P", p_)
            # Back-flashes from the state's own values
            try:
                h = PropsSI("H", "T", T, "P", 1e5, bf)
                s_ = PropsSI("S", "T", T, "P", 1e5, bf)
                d = PropsSI("D", "T", T, "P", 1e5, bf)
                rec("T", "H", h, "P", 1e5)
                rec("T", "P", 1e5, "S", s_)
                rec("T", "D", d, "P", 1e5)
            except ValueError:
                skipped += 3
        # Saturated liquid (works only where a psat curve exists and
        # T > TminPsat).
        rec("P", "T", TL(0.9), "Q", 0.0)
    print(f"incompressible: {len(rows)} records, {skipped} rejected")
    return rows


def gen_humid_air():
    """Humid-air goldens (PLAN 9.1): the full output set over (T, P, x)
    grids incl. sub-freezing ice paths, plus the T-iterating input triples
    (P,W,H), (P,R,W), (P,Tdp,R), (P,W,Twb)."""
    from CoolProp.CoolProp import HAPropsSI
    rows, skipped = [], 0

    def rec(out, n1, v1, n2, v2, n3, v3):
        nonlocal skipped
        try:
            rows.append({"backend": "HA", "fluid": "", "out": out,
                         "name1": n1, "val1": v1, "name2": n2, "val2": v2,
                         "name3": n3, "val3": v3,
                         "expected": HAPropsSI(out, n1, v1, n2, v2, n3, v3)})
        except Exception:
            skipped += 1

    outputs = ["W", "psi_w", "Tdp", "Twb", "H", "Hha", "U", "S", "Sha",
               "V", "Vha", "mu", "k", "cp", "cp_ha", "CV", "P_w", "Z",
               "speed_of_sound", "isentropic_exponent"]
    for T in [253.15, 273.15, 298.15, 333.15, 393.15]:
        for p_ in [101325.0, 5e5, 2e6]:
            for (n3, v3) in [("R", 0.2), ("R", 0.85), ("W", 0.005)]:
                for out in outputs:
                    rec(out, "T", T, "P", p_, n3, v3)
    # Inverse triples (no T given)
    for (n2, v2, n3, v3) in [("W", 0.01, "H", 60000.0), ("R", 0.5, "W", 0.008),
                             ("Tdp", 285.0, "R", 0.6), ("W", 0.01, "B", 295.0),
                             ("R", 0.4, "H", 45000.0), ("W", 0.005, "S", 150.0)]:
        for out in ["T", "W", "R", "H", "V"]:
            rec(out, "P", 101325.0, n2, v2, n3, v3)
    print(f"humid air: {len(rows)} records, {skipped} rejected")
    return rows


PSEUDO_PURE = ["Air", "R404A", "R407C", "R410A", "R507A", "SES36"]


def gen_pseudo_pure():
    """Pseudo-pure flash goldens: PT (liquid/gas/supercritical), QT at the
    only defined qualities (0/1), and PQ across the glide, for all six
    pseudo-pure fluids."""
    rows, skipped = [], 0
    for fluid in PSEUDO_PURE:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        pc = PropsSI("pcrit", "", 0, "", 0, hf)

        def TL(x):
            return Tt + x * (Tc - Tt)

        def rec(out, n1, v1, n2, v2):
            nonlocal skipped
            r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)

        T_mid = TL(0.5)
        p_mid = PropsSI("P", "T", T_mid, "Q", 0, hf)
        # QT at the defined qualities
        for q in [0.0, 1.0]:
            for out in ["P", "Dmolar", "Hmolar", "Smolar"]:
                rec(out, "T", T_mid, "Q", q)
                rec(out, "T", TL(0.3), "Q", q)
        # PQ across the glide (fractional quality is defined for PQ)
        for q in [0.0, 0.4, 1.0]:
            for out in ["T", "Dmolar", "Hmolar", "Umolar"]:
                rec(out, "P", p_mid, "Q", q)
                rec(out, "P", 3.0 * p_mid, "Q", q)
        # PT: liquid (above 1.02*pL), gas (below 0.98*pV), supercritical
        pL = PropsSI("P", "T", T_mid, "Q", 0, hf)
        pV = PropsSI("P", "T", T_mid, "Q", 1, hf)
        for (T, p_) in [(T_mid, 2.0 * pL), (T_mid, 0.5 * pV),
                        (1.2 * Tc, 1.5 * pc), (1.2 * Tc, 0.5 * pc),
                        (TL(0.4), 1.5 * pc)]:
            for out in ["Dmolar", "Hmolar", "Smolar"]:
                rec(out, "T", T, "P", p_)
    print(f"pseudo-pure: {len(rows)} records, {skipped} rejected")
    return rows


# The 20 fluids whose viscosity is fully structured (dilute/initial_density/
# higher_order with typed families only) — the 6.1 structured slice.
VISCOSITY_STRUCTURED = [
    "Ammonia", "Argon", "DimethylEther", "Ethanol", "HydrogenSulfide",
    "IsoButane", "Methane", "Nitrogen", "Oxygen", "R123", "R125", "R134a",
    "SulfurHexafluoride", "n-Butane", "n-Decane", "n-Dodecane", "n-Nonane",
    "n-Octane", "n-Pentane", "n-Propane",
    # fully-hardcoded models (slice 3):
    "Water", "HeavyWater", "Helium", "R23", "Methanol", "m-Xylene",
    "o-Xylene", "p-Xylene",
    # section-hardcoded parts (slice 3):
    "CarbonDioxide", "Ethane", "CycloHexane", "Benzene", "Hydrogen",
    "ParaHydrogen", "Toluene", "n-Hexane", "n-Heptane",
    # Chung + rhosr-CS (slice 4):
    "Cyclopentane", "Isopentane",
    "R1234yf", "R1234ze(E)", "R124", "R152A", "R22", "R245fa", "R32",
    # ECS (slice 5) — references Propane/R134a/Nitrogen:
    "EthylBenzene", "Propylene", "R11", "R116", "R12", "R13", "R14",
    "R141b", "R142b", "R143a", "R218", "R227EA", "R236EA", "R236FA",
    "RC318",
]


def gen_viscosity():
    """6.1 structured-viscosity goldens: V at PT liquid/gas/supercritical
    states and along the saturation curve (incl. a two-phase mixture-density
    state, which upstream evaluates verbatim)."""
    rows, skipped = [], 0
    for fluid in VISCOSITY_STRUCTURED:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        pc = PropsSI("pcrit", "", 0, "", 0, hf)

        def TL(x):
            return Tt + x * (Tc - Tt)

        def psat(T):
            return PropsSI("P", "T", T, "Q", 1, hf)

        def rec(out, n1, v1, n2, v2):
            nonlocal skipped
            r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)

        rec("V", "T", TL(0.3), "P", 2.5 * psat(TL(0.3)))
        rec("V", "T", TL(0.6), "P", 2.0 * psat(TL(0.6)))
        rec("V", "T", TL(0.8), "P", 0.5 * psat(TL(0.8)))
        rec("V", "T", 1.1 * Tc, "P", 1.5 * pc)
        rec("V", "T", TL(0.5), "Q", 0.0)
        rec("V", "T", TL(0.5), "Q", 1.0)
        rec("V", "T", TL(0.5), "Q", 0.5)
        rec("viscosity", "T", TL(0.7), "P", 3.0 * psat(TL(0.7)))
    print(f"viscosity: {len(rows)} records, {skipped} rejected")
    return rows


# The 15 fluids whose conductivity trio AND viscosity are fully structured
# (the Olchowy-Sengers enhancement consumes the fluid's viscosity).
CONDUCTIVITY_STRUCTURED = [
    "Argon", "Ethanol", "IsoButane", "Nitrogen", "Oxygen", "R125", "R134a",
    "SulfurHexafluoride", "n-Butane", "n-Decane", "n-Dodecane", "n-Nonane",
    "n-Octane", "n-Pentane", "n-Propane",
    # fully-hardcoded models (slice 3):
    "Water", "HeavyWater", "Helium", "R23", "Methane",
    # structured fluids unlocked by slice-3 viscosity / hardcoded sections:
    "Ammonia", "R123", "CarbonDioxide", "Ethane", "Benzene", "Methanol",
    "Hydrogen", "ParaHydrogen", "Toluene", "m-Xylene", "o-Xylene",
    "p-Xylene", "n-Hexane", "n-Heptane",
    # unlocked by slice-4 viscosity (their conductivity trio is structured):
    "Cyclopentane", "Isopentane", "R1234yf", "R1234ze(E)", "R152A",
    # ECS conductivity (slice 5):
    "Propylene", "R11", "R116", "R12", "R124", "R13", "R14", "R141b",
    "R142b", "R143a", "R218", "R22", "R227EA", "R236EA", "R236FA",
    "R245fa", "R32", "RC318",
    # structured conductivity unlocked by slice-5 ECS viscosity (OS term):
    "EthylBenzene",
]


def gen_conductivity():
    """6.1 structured-conductivity goldens: L at PT states (incl. the
    near-critical region where the Olchowy-Sengers enhancement dominates)
    and along the saturation curve."""
    rows, skipped = [], 0
    for fluid in CONDUCTIVITY_STRUCTURED:
        hf = f"HEOS::{fluid}"
        Tc = PropsSI("Tcrit", "", 0, "", 0, hf)
        Tt = PropsSI("Ttriple", "", 0, "", 0, hf)
        pc = PropsSI("pcrit", "", 0, "", 0, hf)

        def TL(x):
            return Tt + x * (Tc - Tt)

        def psat(T):
            return PropsSI("P", "T", T, "Q", 1, hf)

        def rec(out, n1, v1, n2, v2):
            nonlocal skipped
            r = try_record(out, n1, v1, n2, v2, "HEOS", fluid)
            rows.append(r) if r else (skipped := skipped + 1)

        rec("L", "T", TL(0.3), "P", 2.5 * psat(TL(0.3)))
        rec("L", "T", TL(0.6), "P", 2.0 * psat(TL(0.6)))
        rec("L", "T", TL(0.8), "P", 0.5 * psat(TL(0.8)))
        rec("L", "T", 1.1 * Tc, "P", 1.5 * pc)
        rec("L", "T", 1.02 * Tc, "P", 1.02 * pc)
        rec("L", "T", TL(0.5), "Q", 0.0)
        rec("L", "T", TL(0.5), "Q", 1.0)
        # Strictly two-phase: upstream returns dilute+residual when the OS
        # numerator gate short-circuits (else throws via cpmolar) —
        # try_record keeps whichever states evaluate.
        rec("L", "T", TL(0.5), "Q", 0.5)
        rec("L", "T", TL(0.98), "Q", 0.0)
        rec("conductivity", "T", TL(0.7), "P", 3.0 * psat(TL(0.7)))
    print(f"conductivity: {len(rows)} records, {skipped} rejected")
    return rows


WRITTEN = []


def write_jsonl(name, rows):
    (FIXTURES / name).write_text("".join(json.dumps(r) + "\n" for r in rows))
    WRITTEN.append(name)
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
    write_jsonl("heos_water_props.jsonl", gen_heos_water_props())
    write_jsonl("heos_water_ancillary.jsonl", gen_heos_water_ancillary())
    write_jsonl("heos_water_sat.jsonl", gen_heos_water_sat())
    write_jsonl("heos_water_pt.jsonl", gen_heos_water_pt())
    write_jsonl("heos_water_flash.jsonl", gen_heos_water_flash())
    for fluid in HEOS_FLUIDS:
        module = module_name(fluid)
        for suite, rows in gen_heos_fluid_suites(fluid).items():
            write_jsonl(f"heos_{module}_{suite}.jsonl", rows)
    for fluid in ["Water"] + HEOS_FLUIDS:
        write_jsonl(f"heos_{module_name(fluid)}_hs.jsonl", gen_heos_fluid_hs(fluid))
    write_jsonl("heos_all_smoke.jsonl", gen_heos_all_smoke())
    write_jsonl("fluid_resolution.jsonl", gen_fluid_resolution())
    write_jsonl("props_si.jsonl", gen_props_si())
    write_jsonl("surface_tension.jsonl", gen_surface_tension())
    write_jsonl("viscosity.jsonl", gen_viscosity())
    write_jsonl("flash_pairs_extra.jsonl", gen_flash_pairs_extra())
    write_jsonl("melting.jsonl", gen_melting())
    write_jsonl("pseudo_pure.jsonl", gen_pseudo_pure())
    write_jsonl("cubics.jsonl", gen_cubics())
    write_jsonl("cubic_superanc.jsonl", gen_cubic_superanc())
    write_jsonl("incompressible.jsonl", gen_incompressible())
    write_jsonl("humid_air.jsonl", gen_humid_air())
    write_jsonl("mixture_helmholtz.jsonl", gen_mixture_helmholtz())
    write_jsonl("mixture_pt.jsonl", gen_mixture_pt())
    write_jsonl("mixture_vle.jsonl", gen_mixture_vle())
    write_jsonl("mixture_propssi.jsonl", gen_mixture_propssi())
    write_jsonl("mixture_predefined.jsonl", gen_mixture_predefined())
    write_jsonl("mixture_pt_twophase.jsonl", gen_mixture_pt_twophase())
    write_jsonl("mixture_sweep.jsonl", gen_mixture_sweep())
    write_jsonl("pcsaft_terms.jsonl", gen_pcsaft_terms())
    write_jsonl("pcsaft_flash.jsonl", gen_pcsaft_flash())
    write_jsonl("partial_derivs.jsonl", gen_partial_derivs())
    write_jsonl("tabular_tables.jsonl", gen_tabular_tables())
    write_jsonl("ttse.jsonl", gen_ttse())
    write_jsonl("conductivity.jsonl", gen_conductivity())
    param_rows = dump_parameters()
    write_jsonl("parameters.jsonl", param_rows)
    write_jsonl("param_aliases.jsonl", dump_param_names(param_rows))
    write_jsonl("phases.jsonl", dump_phases())

    manifest = {
        "generator": "tools/golden-gen/gen_fixtures.py",
        "coolprop_version": CoolProp.__version__,
        "upstream_tag": "v8.0.0",
        "platform": f"{platform.system()}-{platform.machine()}",
        "files": sorted(WRITTEN),
    }
    (FIXTURES / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
