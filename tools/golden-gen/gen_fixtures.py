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
    rows, skipped = [], 0
    for (n1, v1, n2, v2) in sources:
        h = PropsSI("Hmolar", n1, v1, n2, v2, hf)
        s = PropsSI("Smolar", n1, v1, n2, v2, hf)
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
