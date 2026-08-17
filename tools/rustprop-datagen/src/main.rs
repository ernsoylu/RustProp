//! Data pipeline (PLAN.md 3.2): reads fluid runtime-JSON dumps from
//! `data/coolprop-json/` (produced by `tools/golden-gen/dump_fluid_json.py`
//! from the CoolProp 8.0.0 oracle wheel) and emits feature-gated Rust data
//! modules into `crates/rustprop-data/src/fluids/`. JSON parsing lives only
//! here — never in shipped crates.
//!
//! Unknown `alpha0`/`alphar` term types fail loudly; new families are added
//! to `rustprop_core::fluid` as the fluids that need them are ported.
//!
//! Usage: `cargo run -p rustprop-datagen [-- Fluid ...]` (default: every
//! dump found). Output is deterministic; CI regenerates and diffs.

use serde::Deserialize;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

// ---------------------------------------------------------------------------
// Serde model of the ported subset of the fluid document
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct Doc {
    #[serde(rename = "INFO")]
    info: InfoJson,
    /// Parsed lazily: only `EOS[0]` is evaluated by upstream's backend (the
    /// `EOS()` accessor), and alternate blocks may carry term families the
    /// active EOS never uses (e.g. Methanol's `ResidualHelmholtzAssociating`
    /// alternate). Strict term parsing applies to `EOS[0]` alone.
    #[serde(rename = "EOS")]
    eos: Vec<serde_json::Value>,
    #[serde(rename = "ANCILLARIES")]
    ancillaries: AncJson,
    #[serde(rename = "STATES")]
    states: StatesJson,
    /// Classified in code: the structured viscosity form parses strictly;
    /// ECS/Chung/rhosr-CS/fully-hardcoded classes are not yet ported.
    #[serde(rename = "TRANSPORT")]
    transport: Option<serde_json::Value>,
}

#[derive(Deserialize)]
struct InfoJson {
    #[serde(rename = "NAME")]
    name: String,
    #[serde(rename = "CAS")]
    cas: String,
    #[serde(rename = "ALIASES")]
    aliases: Vec<String>,
    /// Absent for 11 documents; upstream then keeps the default-constructed
    /// `EnvironmentalFactorsStruct` (every field `_HUGE`).
    #[serde(rename = "ENVIRONMENTAL")]
    environmental: Option<EnvJson>,
}

/// `INFO.ENVIRONMENTAL` — exactly the seven doubles upstream's
/// `parse_environmental` (FluidLibrary.h) reads with the throwing
/// `get_double`, so all are required here too. The block's `ASHRAE34` and
/// `Name` strings are not carried (ASHRAE34 is PropsSI-unreachable; Name is
/// unread by upstream's parse).
#[derive(Deserialize)]
struct EnvJson {
    #[serde(rename = "GWP20")]
    gwp20: f64,
    #[serde(rename = "GWP100")]
    gwp100: f64,
    #[serde(rename = "GWP500")]
    gwp500: f64,
    #[serde(rename = "ODP")]
    odp: f64,
    #[serde(rename = "HH")]
    hh: f64,
    #[serde(rename = "PH")]
    ph: f64,
    #[serde(rename = "FH")]
    fh: f64,
}

#[derive(Deserialize)]
struct SuperAncJson {
    jexpansions_p: Vec<ChebJson>,
    #[serde(rename = "jexpansions_rhoL")]
    jexpansions_rho_l: Vec<ChebJson>,
    #[serde(rename = "jexpansions_rhoV")]
    jexpansions_rho_v: Vec<ChebJson>,
    meta: MetaJson,
    check_points: Vec<CheckPointJson>,
}

#[derive(Deserialize)]
struct ChebJson {
    xmin: f64,
    xmax: f64,
    coef: Vec<f64>,
}

#[derive(Deserialize)]
struct MetaJson {
    #[serde(rename = "Tcrittrue / K")]
    t_crit_num: f64,
    #[serde(rename = "rhocrittrue / mol/m^3")]
    rho_crit_num: f64,
}

#[derive(Deserialize)]
struct CheckPointJson {
    #[serde(rename = "T / K")]
    t: f64,
    #[serde(rename = "p(mp) / Pa")]
    p: f64,
    #[serde(rename = "rho'(mp) / mol/m^3")]
    rho_l: f64,
    #[serde(rename = "rho''(mp) / mol/m^3")]
    rho_v: f64,
    #[serde(rename = "p(SA)/p(mp)")]
    p_ratio: f64,
    #[serde(rename = "rho'(SA)/rho'(mp)")]
    rho_l_ratio: f64,
    #[serde(rename = "rho''(SA)/rho''(mp)")]
    rho_v_ratio: f64,
}

#[derive(Deserialize)]
struct EosJson {
    gas_constant: f64,
    molar_mass: f64,
    p_max: f64,
    #[serde(rename = "T_max")]
    t_max: f64,
    #[serde(rename = "Ttriple")]
    t_triple: f64,
    acentric: f64,
    pseudo_pure: bool,
    #[serde(rename = "STATES")]
    states: EosStatesJson,
    alpha0: Vec<Alpha0Json>,
    alphar: Vec<AlpharJson>,
    #[serde(rename = "SUPERANCILLARY")]
    superancillary: Option<SuperAncJson>,
}

#[derive(Deserialize)]
struct EosStatesJson {
    reducing: Sp,
    sat_min_liquid: Sp,
    sat_min_vapor: Sp,
    hs_anchor: Sp,
    temperature_max_sat: Option<Sp>,
    pressure_max_sat: Option<Sp>,
}

/// `hmolar`/`smolar` default to NaN: the `sat_min_liquid`/`sat_min_vapor`
/// states of a few documents (R1130(E), R1243zf) omit them, and upstream's
/// FluidLibrary never reads them from those states either.
#[derive(Deserialize)]
struct Sp {
    #[serde(rename = "T")]
    t: f64,
    p: f64,
    rhomolar: f64,
    #[serde(default = "f64_nan")]
    hmolar: f64,
    #[serde(default = "f64_nan")]
    smolar: f64,
}

fn f64_nan() -> f64 {
    f64::NAN
}

#[derive(Deserialize)]
#[serde(tag = "type")]
enum Alpha0Json {
    #[serde(rename = "IdealGasHelmholtzLead")]
    Lead { a1: f64, a2: f64 },
    #[serde(rename = "IdealGasHelmholtzLogTau")]
    LogTau { a: f64 },
    #[serde(rename = "IdealGasHelmholtzPlanckEinstein")]
    PlanckEinstein { n: Vec<f64>, t: Vec<f64> },
    #[serde(rename = "IdealGasHelmholtzPlanckEinsteinFunctionT")]
    PlanckEinsteinFunctionT {
        n: Vec<f64>,
        v: Vec<f64>,
        #[serde(rename = "Tcrit")]
        tcrit: f64,
    },
    #[serde(rename = "IdealGasHelmholtzEnthalpyEntropyOffset")]
    EnthalpyEntropyOffset { a1: f64, a2: f64, reference: String },
    #[serde(rename = "IdealGasHelmholtzPower")]
    Power { n: Vec<f64>, t: Vec<f64> },
    #[serde(rename = "IdealGasHelmholtzPlanckEinsteinGeneralized")]
    PlanckEinsteinGeneralized {
        n: Vec<f64>,
        t: Vec<f64>,
        c: Vec<f64>,
        d: Vec<f64>,
    },
    #[serde(rename = "IdealGasHelmholtzCP0Constant")]
    Cp0Constant {
        #[serde(rename = "cp_over_R")]
        cp_over_r: f64,
        #[serde(rename = "Tc")]
        tc: f64,
        #[serde(rename = "T0")]
        t0: f64,
    },
    #[serde(rename = "IdealGasHelmholtzCP0PolyT")]
    Cp0PolyT {
        c: Vec<f64>,
        t: Vec<f64>,
        #[serde(rename = "Tc")]
        tc: f64,
        #[serde(rename = "T0")]
        t0: f64,
    },
    #[serde(rename = "IdealGasHelmholtzCP0AlyLee")]
    Cp0AlyLee {
        c: Vec<f64>,
        #[serde(rename = "Tc")]
        tc: f64,
        #[serde(rename = "T0")]
        t0: f64,
    },
}

#[derive(Deserialize)]
#[serde(tag = "type")]
enum AlpharJson {
    #[serde(rename = "ResidualHelmholtzPower")]
    Power {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        l: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzGaussian")]
    Gaussian {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        eta: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
        epsilon: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzNonAnalytic")]
    NonAnalytic {
        n: Vec<f64>,
        a: Vec<f64>,
        b: Vec<f64>,
        beta: Vec<f64>,
        #[serde(rename = "A")]
        big_a: Vec<f64>,
        #[serde(rename = "B")]
        big_b: Vec<f64>,
        #[serde(rename = "C")]
        big_c: Vec<f64>,
        #[serde(rename = "D")]
        big_d: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzExponential")]
    Exponential {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        g: Vec<f64>,
        l: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzDoubleExponential")]
    DoubleExponential {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        gd: Vec<f64>,
        ld: Vec<f64>,
        gt: Vec<f64>,
        lt: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzLemmon2005")]
    Lemmon2005 {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        l: Vec<f64>,
        m: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzGaoB")]
    GaoB {
        n: Vec<f64>,
        t: Vec<f64>,
        d: Vec<f64>,
        eta: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
        epsilon: Vec<f64>,
        b: Vec<f64>,
    },
}

#[derive(Deserialize)]
struct AncJson {
    /// Pure fluids carry a single `pS` curve; pseudo-pure fluids carry
    /// separate `pL`/`pV` curves (upstream `parse_ancillaries`).
    #[serde(rename = "pS")]
    p_s: Option<SatAncJson>,
    #[serde(rename = "pL")]
    p_l: Option<SatAncJson>,
    #[serde(rename = "pV")]
    p_v: Option<SatAncJson>,
    #[serde(rename = "rhoL")]
    rho_l: SatAncJson,
    #[serde(rename = "rhoV")]
    rho_v: SatAncJson,
    surface_tension: Option<SurfTensJson>,
    melting_line: Option<MeltJson>,
}

#[derive(Deserialize)]
struct SurfTensJson {
    a: Vec<f64>,
    n: Vec<f64>,
    #[serde(rename = "Tc")]
    tc: f64,
}

#[derive(Deserialize)]
struct MeltJson {
    #[serde(rename = "T_m")]
    t_m: f64,
    #[serde(rename = "type")]
    melt_type: String,
    parts: Vec<MeltPartJson>,
}

#[derive(Deserialize)]
struct MeltPartJson {
    #[serde(rename = "T_0")]
    t_0: f64,
    #[serde(rename = "p_0")]
    p_0: f64,
    #[serde(rename = "T_min")]
    t_min: f64,
    #[serde(rename = "T_max")]
    t_max: f64,
    // Simon segments
    a: Option<serde_json::Value>,
    c: Option<f64>,
    // Polynomial segments
    t: Option<Vec<f64>>,
}

#[derive(Deserialize)]
struct SatAncJson {
    #[serde(rename = "type")]
    anc_type: String,
    n: Vec<f64>,
    t: Vec<f64>,
    #[serde(rename = "T_r")]
    t_r: f64,
    reducing_value: f64,
    using_tau_r: bool,
    #[serde(rename = "Tmin")]
    t_min: f64,
    #[serde(rename = "Tmax")]
    t_max: f64,
}

#[derive(Deserialize)]
struct StatesJson {
    critical: Sp,
    triple_liquid: Sp,
    triple_vapor: Sp,
}

// ---------------------------------------------------------------------------
// Emission
// ---------------------------------------------------------------------------

/// Shortest round-trip f64 literal (Rust `{:?}` guarantees exact re-parse).
/// NaN (a document field upstream never reads, e.g. a sat_min caloric)
/// emits as the `f64::NAN` path expression; +inf (upstream's `_HUGE`
/// default for a missing ENVIRONMENTAL block) as `f64::INFINITY`.
fn f(v: f64) -> String {
    if v.is_nan() {
        return "f64::NAN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 {
            "f64::INFINITY"
        } else {
            "f64::NEG_INFINITY"
        }
        .into();
    }
    format!("{v:?}")
}

fn slice(vals: &[f64]) -> String {
    let items: Vec<String> = vals.iter().map(|v| f(*v)).collect();
    format!("&[{}]", items.join(", "))
}

/// `Water` -> `water`, `1-Butene` -> `_1_butene` (a valid Rust identifier).
fn module_name(fluid: &str) -> String {
    let mut s: String = fluid
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect();
    if s.starts_with(|c: char| c.is_ascii_digit()) {
        s.insert(0, '_');
    }
    s
}

fn state_point(sp: &Sp, indent: &str) -> String {
    format!(
        "StatePoint {{\n{indent}    t: {},\n{indent}    p: {},\n{indent}    rhomolar: {},\n{indent}    hmolar: {},\n{indent}    smolar: {},\n{indent}}}",
        f(sp.t),
        f(sp.p),
        f(sp.rhomolar),
        f(sp.hmolar),
        f(sp.smolar),
    )
}

fn sat_anc(a: &SatAncJson, indent: &str) -> String {
    format!(
        "SaturationAncillary {{\n{indent}    anc_type: {:?},\n{indent}    n: {},\n{indent}    t: {},\n{indent}    t_r: {},\n{indent}    reducing_value: {},\n{indent}    using_tau_r: {},\n{indent}    t_min: {},\n{indent}    t_max: {},\n{indent}}}",
        a.anc_type,
        slice(&a.n),
        slice(&a.t),
        f(a.t_r),
        f(a.reducing_value),
        a.using_tau_r,
        f(a.t_min),
        f(a.t_max),
    )
}

fn cheb_intervals(name: &str, blocks: &[ChebJson], out: &mut String) {
    writeln!(out, "        {name}: &[").unwrap();
    for b in blocks {
        writeln!(
            out,
            "            ChebyshevInterval {{ xmin: {}, xmax: {}, coef: {} }},",
            f(b.xmin),
            f(b.xmax),
            slice(&b.coef)
        )
        .unwrap();
    }
    writeln!(out, "        ],").unwrap();
}

fn emit_superancillary(sa: &SuperAncJson, w: &mut String) {
    writeln!(w, "        superancillary: Some(SuperAncillaryData {{").unwrap();
    let mut inner = String::new();
    cheb_intervals("p", &sa.jexpansions_p, &mut inner);
    cheb_intervals("rho_l", &sa.jexpansions_rho_l, &mut inner);
    cheb_intervals("rho_v", &sa.jexpansions_rho_v, &mut inner);
    // shift the inner indentation by four spaces to sit inside Some(..)
    for line in inner.lines() {
        writeln!(w, "    {line}").unwrap();
    }
    writeln!(w, "            t_crit_num: {},", f(sa.meta.t_crit_num)).unwrap();
    writeln!(w, "            rho_crit_num: {},", f(sa.meta.rho_crit_num)).unwrap();
    writeln!(w, "            check_points: &[").unwrap();
    for c in &sa.check_points {
        writeln!(
            w,
            "                SuperAncCheckPoint {{ t: {}, p: {}, rho_l: {}, rho_v: {}, p_ratio: {}, rho_l_ratio: {}, rho_v_ratio: {} }},",
            f(c.t), f(c.p), f(c.rho_l), f(c.rho_v), f(c.p_ratio), f(c.rho_l_ratio), f(c.rho_v_ratio)
        )
        .unwrap();
    }
    writeln!(w, "            ],").unwrap();
    writeln!(w, "        }}),").unwrap();
}

fn render_alpha0(w: &mut String, terms: &[Alpha0Json], indent: &str) {
    writeln!(w, "{indent}alpha0: &[").unwrap();
    for term in terms {
        match term {
            Alpha0Json::Lead { a1, a2 } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::Lead {{ a1: {}, a2: {} }},",
                    f(*a1),
                    f(*a2)
                )
                .unwrap();
            }
            Alpha0Json::LogTau { a } => {
                writeln!(w, "{indent}    Alpha0Term::LogTau {{ a: {} }},", f(*a)).unwrap();
            }
            Alpha0Json::PlanckEinstein { n, t } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::PlanckEinstein {{ n: {}, t: {} }},",
                    slice(n),
                    slice(t)
                )
                .unwrap();
            }
            Alpha0Json::PlanckEinsteinFunctionT { n, v, tcrit } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::PlanckEinsteinFunctionT {{ n: {}, v: {}, tcrit: {} }},",
                    slice(n),
                    slice(v),
                    f(*tcrit)
                )
                .unwrap();
            }
            Alpha0Json::EnthalpyEntropyOffset { a1, a2, reference } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::EnthalpyEntropyOffset {{ a1: {}, a2: {}, reference: {:?} }},",
                    f(*a1),
                    f(*a2),
                    reference
                )
                .unwrap();
            }
            Alpha0Json::Power { n, t } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::Power {{ n: {}, t: {} }},",
                    slice(n),
                    slice(t)
                )
                .unwrap();
            }
            Alpha0Json::PlanckEinsteinGeneralized { n, t, c, d } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::PlanckEinsteinGeneralized {{ n: {}, t: {}, c: {}, d: {} }},",
                    slice(n),
                    slice(t),
                    slice(c),
                    slice(d)
                )
                .unwrap();
            }
            Alpha0Json::Cp0Constant { cp_over_r, tc, t0 } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::Cp0Constant {{ cp_over_r: {}, tc: {}, t0: {} }},",
                    f(*cp_over_r),
                    f(*tc),
                    f(*t0)
                )
                .unwrap();
            }
            Alpha0Json::Cp0PolyT { c, t, tc, t0 } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::Cp0PolyT {{ c: {}, t: {}, tc: {}, t0: {} }},",
                    slice(c),
                    slice(t),
                    f(*tc),
                    f(*t0)
                )
                .unwrap();
            }
            Alpha0Json::Cp0AlyLee { c, tc, t0 } => {
                writeln!(
                    w,
                    "{indent}    Alpha0Term::Cp0AlyLee {{ c: {}, tc: {}, t0: {} }},",
                    slice(c),
                    f(*tc),
                    f(*t0)
                )
                .unwrap();
            }
        }
    }
    writeln!(w, "{indent}],").unwrap();
}

fn emit(doc: &Doc, eos: &EosJson, source_file: &str) -> String {
    let _ = &doc.eos;
    let mut out = String::new();
    let w = &mut out;
    let static_name = module_name(&doc.info.name).to_uppercase();
    writeln!(
        w,
        "//! GENERATED by rustprop-datagen from {source_file} — DO NOT EDIT."
    )
    .unwrap();
    writeln!(w, "//! Regenerate: cargo run -p rustprop-datagen").unwrap();
    writeln!(w).unwrap();
    writeln!(w, "#![cfg_attr(rustfmt, rustfmt::skip)]").unwrap();
    // Verbatim upstream coefficients may approximate mathematical constants
    // (e.g. a fitted exponent of 3.14); the data is bit-faithful by mandate.
    writeln!(w, "#![allow(clippy::approx_constant)]").unwrap();
    writeln!(w).unwrap();
    let surf = if doc.ancillaries.surface_tension.is_some() {
        ", SurfaceTension"
    } else {
        ""
    };
    let sa = if eos.superancillary.is_some() {
        ", ChebyshevInterval, SuperAncCheckPoint, SuperAncillaryData"
    } else {
        ""
    };
    let melt = match &doc.ancillaries.melting_line {
        Some(ml) if ml.melt_type == "Simon" => ", MeltingLine, MeltingLineKind, SimonMeltPart",
        Some(_) => ", MeltingLine, MeltingLineKind, PolyMeltPart",
        None => "",
    };
    let transport_rendered = emit_transport(doc.transport.as_ref());
    let mut timp = String::new();
    for ty in [
        "Transport",
        "TransportModel",
        "ViscosityModel",
        "ConductivityModel",
        "Viscosity",
        "ViscosityDilute",
        "ViscosityInitialDensity",
        "ViscosityHigherOrder",
        "Conductivity",
        "ConductivityDilute",
        "ConductivityResidual",
        "ConductivityCritical",
    ] {
        if transport_rendered.contains(&format!("{ty} {{"))
            || transport_rendered.contains(&format!("{ty}::"))
        {
            timp.push_str(", ");
            timp.push_str(ty);
        }
    }
    writeln!(
        w,
        "use rustprop_core::fluid::{{Alpha0Term, AlpharTerm, Ancillaries, Environmental, Eos, FluidData, SaturationAncillary, StatePoint, States{sa}{surf}{melt}{timp}}};"
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(w, "pub static {static_name}: FluidData = FluidData {{").unwrap();
    writeln!(w, "    name: {:?},", doc.info.name).unwrap();
    writeln!(w, "    cas: {:?},", doc.info.cas).unwrap();
    let aliases: Vec<String> = doc.info.aliases.iter().map(|a| format!("{a:?}")).collect();
    writeln!(w, "    aliases: &[{}],", aliases.join(", ")).unwrap();
    // Missing block -> upstream's default-constructed struct: every field
    // _HUGE (FluidLibrary.cpp, "Environmental data are missing" branch).
    const HUGE: f64 = f64::INFINITY;
    let env = doc.info.environmental.as_ref();
    writeln!(w, "    environmental: Environmental {{").unwrap();
    for (field, v) in [
        ("gwp20", env.map_or(HUGE, |e| e.gwp20)),
        ("gwp100", env.map_or(HUGE, |e| e.gwp100)),
        ("gwp500", env.map_or(HUGE, |e| e.gwp500)),
        ("odp", env.map_or(HUGE, |e| e.odp)),
        ("hh", env.map_or(HUGE, |e| e.hh)),
        ("ph", env.map_or(HUGE, |e| e.ph)),
        ("fh", env.map_or(HUGE, |e| e.fh)),
    ] {
        writeln!(w, "        {field}: {},", f(v)).unwrap();
    }
    writeln!(w, "    }},").unwrap();
    writeln!(w, "    eos: Eos {{").unwrap();
    writeln!(w, "        gas_constant: {},", f(eos.gas_constant)).unwrap();
    writeln!(w, "        molar_mass: {},", f(eos.molar_mass)).unwrap();
    writeln!(w, "        p_max: {},", f(eos.p_max)).unwrap();
    writeln!(w, "        t_max: {},", f(eos.t_max)).unwrap();
    writeln!(w, "        t_triple: {},", f(eos.t_triple)).unwrap();
    writeln!(w, "        acentric: {},", f(eos.acentric)).unwrap();
    writeln!(w, "        pseudo_pure: {},", eos.pseudo_pure).unwrap();
    writeln!(
        w,
        "        reducing: {},",
        state_point(&eos.states.reducing, "        ")
    )
    .unwrap();
    writeln!(
        w,
        "        sat_min_liquid: {},",
        state_point(&eos.states.sat_min_liquid, "        ")
    )
    .unwrap();
    writeln!(
        w,
        "        sat_min_vapor: {},",
        state_point(&eos.states.sat_min_vapor, "        ")
    )
    .unwrap();
    writeln!(
        w,
        "        hs_anchor: {},",
        state_point(&eos.states.hs_anchor, "        ")
    )
    .unwrap();
    for (field, sp) in [
        ("max_sat_t", &eos.states.temperature_max_sat),
        ("max_sat_p", &eos.states.pressure_max_sat),
    ] {
        match sp {
            Some(sp) => {
                writeln!(w, "        {field}: Some({}),", state_point(sp, "        ")).unwrap()
            }
            None => writeln!(w, "        {field}: None,").unwrap(),
        }
    }
    render_alpha0(w, &eos.alpha0, "        ");
    writeln!(w, "        alphar: &[").unwrap();
    for term in &eos.alphar {
        match term {
            AlpharJson::Power { n, d, t, l } => {
                writeln!(
                    w,
                    "            AlpharTerm::Power {{ n: {}, d: {}, t: {}, l: {} }},",
                    slice(n),
                    slice(d),
                    slice(t),
                    slice(l)
                )
                .unwrap();
            }
            AlpharJson::Gaussian {
                n,
                d,
                t,
                eta,
                beta,
                gamma,
                epsilon,
            } => {
                writeln!(
                    w,
                    "            AlpharTerm::Gaussian {{ n: {}, d: {}, t: {}, eta: {}, beta: {}, gamma: {}, epsilon: {} }},",
                    slice(n), slice(d), slice(t), slice(eta), slice(beta), slice(gamma), slice(epsilon)
                )
                .unwrap();
            }
            AlpharJson::NonAnalytic {
                n,
                a,
                b,
                beta,
                big_a,
                big_b,
                big_c,
                big_d,
            } => {
                writeln!(
                    w,
                    "            AlpharTerm::NonAnalytic {{ n: {}, a: {}, b: {}, beta: {}, big_a: {}, big_b: {}, big_c: {}, big_d: {} }},",
                    slice(n), slice(a), slice(b), slice(beta), slice(big_a), slice(big_b), slice(big_c), slice(big_d)
                )
                .unwrap();
            }
            AlpharJson::Exponential { n, d, t, g, l } => {
                writeln!(
                    w,
                    "            AlpharTerm::Exponential {{ n: {}, d: {}, t: {}, g: {}, l: {} }},",
                    slice(n),
                    slice(d),
                    slice(t),
                    slice(g),
                    slice(l)
                )
                .unwrap();
            }
            AlpharJson::DoubleExponential {
                n,
                d,
                t,
                gd,
                ld,
                gt,
                lt,
            } => {
                writeln!(
                    w,
                    "            AlpharTerm::DoubleExponential {{ n: {}, d: {}, t: {}, gd: {}, ld: {}, gt: {}, lt: {} }},",
                    slice(n), slice(d), slice(t), slice(gd), slice(ld), slice(gt), slice(lt)
                )
                .unwrap();
            }
            AlpharJson::Lemmon2005 { n, d, t, l, m } => {
                writeln!(
                    w,
                    "            AlpharTerm::Lemmon2005 {{ n: {}, d: {}, t: {}, l: {}, m: {} }},",
                    slice(n),
                    slice(d),
                    slice(t),
                    slice(l),
                    slice(m)
                )
                .unwrap();
            }
            AlpharJson::GaoB {
                n,
                t,
                d,
                eta,
                beta,
                gamma,
                epsilon,
                b,
            } => {
                writeln!(
                    w,
                    "            AlpharTerm::GaoB {{ n: {}, t: {}, d: {}, eta: {}, beta: {}, gamma: {}, epsilon: {}, b: {} }},",
                    slice(n), slice(t), slice(d), slice(eta), slice(beta), slice(gamma), slice(epsilon), slice(b)
                )
                .unwrap();
            }
        }
    }
    writeln!(w, "        ],").unwrap();
    match &eos.superancillary {
        Some(sa) => emit_superancillary(sa, w),
        None => writeln!(w, "        superancillary: None,").unwrap(),
    }
    writeln!(w, "    }},").unwrap();
    writeln!(w, "    ancillaries: Ancillaries {{").unwrap();
    // Upstream parse_ancillaries: pL+pV present (pseudo-pure) fill the two
    // slots; else a single pS fills both. We alias the pure case through
    // `p_v_split: None` instead of duplicating the curve.
    let (p_l_slot, p_v_slot) = match (
        &doc.ancillaries.p_s,
        &doc.ancillaries.p_l,
        &doc.ancillaries.p_v,
    ) {
        (None, Some(pl), Some(pv)) => (pl, Some(pv)),
        (Some(ps), None, None) => (ps, None),
        _ => panic!("Pressure ancillary curves are missing or invalid"),
    };
    writeln!(w, "        p_s: {},", sat_anc(p_l_slot, "        ")).unwrap();
    match p_v_slot {
        Some(pv) => {
            writeln!(w, "        p_v_split: Some({}),", sat_anc(pv, "        ")).unwrap();
        }
        None => writeln!(w, "        p_v_split: None,").unwrap(),
    }
    writeln!(
        w,
        "        rho_l: {},",
        sat_anc(&doc.ancillaries.rho_l, "        ")
    )
    .unwrap();
    writeln!(
        w,
        "        rho_v: {},",
        sat_anc(&doc.ancillaries.rho_v, "        ")
    )
    .unwrap();
    match &doc.ancillaries.surface_tension {
        Some(st) => writeln!(
            w,
            "        surface_tension: Some(SurfaceTension {{ a: {}, n: {}, tc: {} }}),",
            slice(&st.a),
            slice(&st.n),
            f(st.tc)
        )
        .unwrap(),
        None => writeln!(w, "        surface_tension: None,").unwrap(),
    }
    match &doc.ancillaries.melting_line {
        Some(ml) => {
            writeln!(w, "        melting_line: Some(MeltingLine {{").unwrap();
            writeln!(w, "            t_m: {},", f(ml.t_m)).unwrap();
            let (variant, is_simon) = match ml.melt_type.as_str() {
                "Simon" => ("Simon", true),
                "polynomial_in_Tr" => ("PolynomialInTr", false),
                "polynomial_in_Theta" => ("PolynomialInTheta", false),
                other => panic!("unknown melting_line type {other:?}"),
            };
            writeln!(
                w,
                "            kind: MeltingLineKind::{variant} {{ parts: &["
            )
            .unwrap();
            for part in &ml.parts {
                if is_simon {
                    let a = part
                        .a
                        .as_ref()
                        .and_then(serde_json::Value::as_f64)
                        .expect("Simon melting part requires scalar a");
                    writeln!(
                        w,
                        "                SimonMeltPart {{ t_0: {}, a: {}, c: {}, p_0: {}, t_min: {}, t_max: {} }},",
                        f(part.t_0),
                        f(a),
                        f(part.c.expect("Simon melting part requires c")),
                        f(part.p_0),
                        f(part.t_min),
                        f(part.t_max)
                    )
                    .unwrap();
                } else {
                    let a: Vec<f64> = part
                        .a
                        .as_ref()
                        .and_then(|v| v.as_array().cloned())
                        .expect("polynomial melting part requires array a")
                        .iter()
                        .map(|x| x.as_f64().expect("numeric a"))
                        .collect();
                    writeln!(
                        w,
                        "                PolyMeltPart {{ t_0: {}, p_0: {}, t_min: {}, t_max: {}, a: {}, t: {} }},",
                        f(part.t_0),
                        f(part.p_0),
                        f(part.t_min),
                        f(part.t_max),
                        slice(&a),
                        slice(part.t.as_ref().expect("polynomial melting part requires t"))
                    )
                    .unwrap();
                }
            }
            writeln!(w, "            ] }},").unwrap();
            writeln!(w, "        }}),").unwrap();
        }
        None => writeln!(w, "        melting_line: None,").unwrap(),
    }
    writeln!(w, "    }},").unwrap();
    writeln!(w, "    states: States {{").unwrap();
    writeln!(
        w,
        "        critical: {},",
        state_point(&doc.states.critical, "        ")
    )
    .unwrap();
    writeln!(
        w,
        "        triple_liquid: {},",
        state_point(&doc.states.triple_liquid, "        ")
    )
    .unwrap();
    writeln!(
        w,
        "        triple_vapor: {},",
        state_point(&doc.states.triple_vapor, "        ")
    )
    .unwrap();
    writeln!(w, "    }},").unwrap();
    w.push_str(&emit_transport(doc.transport.as_ref()));
    writeln!(w, "}};").unwrap();
    out
}

/// Classify + emit `TRANSPORT`. Returns the rendered `transport:` field.
/// Per-property tri-state: Absent (no key -> upstream "model not
/// available"), Unported (ECS/Chung/rhosr lists/fully-hardcoded), or the
/// strictly-parsed structured Model.
fn emit_transport(tr: Option<&serde_json::Value>) -> String {
    let Some(tr) = tr else {
        return "    transport: None,\n".into();
    };
    let slot = |key: &str,
                model: &str,
                render: &dyn Fn(&serde_json::Map<String, serde_json::Value>) -> String|
     -> String {
        // Upstream parse_viscosity/parse_conductivity: an array uses its
        // FIRST entry.
        let entry = match tr.get(key) {
            None => return "TransportModel::Absent".into(),
            Some(v) if v.is_array() => &v.as_array().unwrap()[0],
            Some(v) => v,
        };
        let v = entry.as_object().expect("transport entry is an object");
        if let Some(h) = v.get("hardcoded") {
            // Fully-hardcoded per-fluid formulation.
            return format!(
                "TransportModel::Model({model}::Hardcoded {{ name: {:?} }})",
                h.as_str().unwrap()
            );
        }
        match v.get("type").and_then(serde_json::Value::as_str) {
            Some("Chung") if model == "ViscosityModel" => format!(
                "TransportModel::Model(ViscosityModel::Chung {{ rhomolar_critical: {}, acentric: {}, molar_mass: {}, t_critical: {}, dipole_moment_d: {}, kappa: {} }})",
                f(jnum(entry, "rhomolar_critical")),
                f(jnum(entry, "acentric")),
                f(jnum(entry, "molar_mass")),
                f(jnum(entry, "T_critical")),
                f(jnum(entry, "dipole_moment_D")),
                f(jnum(entry, "kappa"))
            ),
            Some("rhosr-CS") if model == "ViscosityModel" => format!(
                "TransportModel::Model(ViscosityModel::RhosrCs {{ c: {}, c_liq: {}, c_vap: {}, rhosr_critical: {}, x_crossover: {} }})",
                f(jnum(entry, "C")),
                slice(&jarr(entry, "c_liq")),
                slice(&jarr(entry, "c_vap")),
                f(jnum(entry, "rhosr_critical")),
                f(jnum(entry, "x_crossover"))
            ),
            Some("ECS") if model == "ViscosityModel" => {
                let psi = entry.get("psi").expect("ECS viscosity has psi");
                format!(
                    "TransportModel::Model(ViscosityModel::Ecs {{ reference_fluid: {:?}, psi_a: {}, psi_t: {}, psi_rhomolar_reducing: {}, sigma_eta: {}, epsilon_over_k: {} }})",
                    entry
                        .get("reference_fluid")
                        .and_then(serde_json::Value::as_str)
                        .unwrap(),
                    slice(&jarr(psi, "a")),
                    slice(&jarr(psi, "t")),
                    f(jnum(psi, "rhomolar_reducing")),
                    f(jnum_opt(entry, "sigma_eta", f64::NAN)),
                    f(jnum_opt(entry, "epsilon_over_k", f64::NAN))
                )
            }
            Some("ECS") => {
                let psi = entry.get("psi").expect("ECS conductivity has psi");
                let fint = entry.get("f_int").expect("ECS conductivity has f_int");
                format!(
                    "TransportModel::Model(ConductivityModel::Ecs {{ reference_fluid: {:?}, psi_a: {}, psi_t: {}, psi_rhomolar_reducing: {}, f_int_a: {}, f_int_t: {}, f_int_t_reducing: {} }})",
                    entry
                        .get("reference_fluid")
                        .and_then(serde_json::Value::as_str)
                        .unwrap(),
                    slice(&jarr(psi, "a")),
                    slice(&jarr(psi, "t")),
                    f(jnum(psi, "rhomolar_reducing")),
                    slice(&jarr(fint, "a")),
                    slice(&jarr(fint, "t")),
                    f(jnum(fint, "T_reducing"))
                )
            }
            Some(_) => "TransportModel::Unported".into(),
            None => render(v),
        }
    };
    let visc = slot("viscosity", "ViscosityModel", &|v| {
        format!(
            "TransportModel::Model(ViscosityModel::Structured(Viscosity {{\n{}        }}))",
            render_viscosity(v)
        )
    });
    let cond = slot("conductivity", "ConductivityModel", &|c| {
        format!(
            "TransportModel::Model(ConductivityModel::Structured(Conductivity {{\n{}        }}))",
            render_conductivity(c)
        )
    });
    format!(
        "    transport: Some(Transport {{\n        viscosity: {visc},\n        conductivity: {cond},\n    }}),\n"
    )
}

fn render_conductivity(c: &serde_json::Map<String, serde_json::Value>) -> String {
    let mut out = String::new();
    let w = &mut out;
    let dil = c.get("dilute").expect("structured conductivity has dilute");
    let dilute = if let Some(h) = dil.get("hardcoded") {
        format!(
            "ConductivityDilute::Hardcoded {{ name: {:?} }}",
            h.as_str().unwrap()
        )
    } else {
        match dil.get("type").and_then(serde_json::Value::as_str).unwrap() {
            "ratio_of_polynomials" => format!(
                "ConductivityDilute::RatioOfPolynomials {{ a: {}, n: {}, b: {}, m: {}, t_reducing: {} }}",
                slice(&jarr(dil, "A")),
                slice(&jarr(dil, "n")),
                slice(&jarr(dil, "B")),
                slice(&jarr(dil, "m")),
                f(jnum(dil, "T_reducing"))
            ),
            "eta0_and_poly" => format!(
                "ConductivityDilute::Eta0AndPoly {{ a: {}, t: {} }}",
                slice(&jarr(dil, "A")),
                slice(&jarr(dil, "t"))
            ),
            other => panic!("unknown dilute conductivity type {other}"),
        }
    };
    writeln!(w, "            dilute: {dilute},").unwrap();

    let res = c
        .get("residual")
        .expect("structured conductivity has residual");
    let residual = match res.get("type").and_then(serde_json::Value::as_str).unwrap() {
        "polynomial" => format!(
            "ConductivityResidual::Polynomial {{ b: {}, t: {}, d: {}, t_reducing: {}, rhomass_reducing: {} }}",
            slice(&jarr(res, "B")),
            slice(&jarr(res, "t")),
            slice(&jarr(res, "d")),
            f(jnum(res, "T_reducing")),
            f(jnum(res, "rhomass_reducing"))
        ),
        "polynomial_and_exponential" => format!(
            "ConductivityResidual::PolynomialAndExponential {{ a: {}, t: {}, d: {}, gamma: {}, l: {} }}",
            slice(&jarr(res, "A")),
            slice(&jarr(res, "t")),
            slice(&jarr(res, "d")),
            slice(&jarr(res, "gamma")),
            slice(&jarr(res, "l"))
        ),
        other => panic!("unknown residual conductivity type {other}"),
    };
    writeln!(w, "            residual: {residual},").unwrap();

    match c.get("critical") {
        None => writeln!(w, "            critical: None,").unwrap(),
        Some(cr) => {
            let body = if let Some(h) = cr.get("hardcoded") {
                format!(
                    "ConductivityCritical::Hardcoded {{ name: {:?} }}",
                    h.as_str().unwrap()
                )
            } else {
                match cr.get("type").and_then(serde_json::Value::as_str).unwrap() {
                    "simplified_Olchowy_Sengers" => format!(
                        "ConductivityCritical::SimplifiedOlchowySengers {{ k: {}, r0: {}, gamma: {}, nu: {}, big_gamma: {}, zeta0: {}, qd: {}, t_ref: {} }}",
                        f(jnum_opt(cr, "k", 1.3806488e-23)),
                        f(jnum_opt(cr, "R0", 1.03)),
                        f(jnum_opt(cr, "gamma", 1.239)),
                        f(jnum_opt(cr, "nu", 0.63)),
                        f(jnum_opt(cr, "GAMMA", 0.0496)),
                        f(jnum_opt(cr, "zeta0", 1.94e-10)),
                        f(jnum_opt(cr, "qD", 2e9)),
                        f(jnum_opt(cr, "T_ref", f64::NAN))
                    ),
                    other => panic!("unknown critical conductivity type {other}"),
                }
            };
            writeln!(w, "            critical: Some({body}),").unwrap();
        }
    }
    out
}

fn jnum(v: &serde_json::Value, key: &str) -> f64 {
    v.get(key)
        .and_then(serde_json::Value::as_f64)
        .unwrap_or_else(|| panic!("missing/invalid numeric key {key}"))
}

fn jarr(v: &serde_json::Value, key: &str) -> Vec<f64> {
    v.get(key)
        .and_then(serde_json::Value::as_array)
        .unwrap_or_else(|| panic!("missing array key {key}"))
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

/// Optional array; empty when absent (upstream's empty-vector channels).
fn jarr_opt(v: &serde_json::Value, key: &str) -> Vec<f64> {
    match v.get(key) {
        Some(a) => a
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect(),
        None => Vec::new(),
    }
}

fn jnum_opt(v: &serde_json::Value, key: &str, default: f64) -> f64 {
    v.get(key)
        .and_then(serde_json::Value::as_f64)
        .unwrap_or(default)
}

fn render_viscosity(v: &serde_json::Map<String, serde_json::Value>) -> String {
    let mut out = String::new();
    let w = &mut out;
    let eps = v.get("epsilon_over_k").and_then(serde_json::Value::as_f64);
    let sig = v.get("sigma_eta").and_then(serde_json::Value::as_f64);
    writeln!(
        w,
        "            epsilon_over_k: {},",
        eps.map_or("f64::NAN".into(), f)
    )
    .unwrap();
    writeln!(
        w,
        "            sigma_eta: {},",
        sig.map_or("f64::NAN".into(), f)
    )
    .unwrap();

    let dil = v.get("dilute").expect("structured viscosity has dilute");
    let dilute = if let Some(h) = dil.get("hardcoded") {
        format!(
            "ViscosityDilute::Hardcoded {{ name: {:?} }}",
            h.as_str().unwrap()
        )
    } else {
        match dil.get("type").and_then(serde_json::Value::as_str).unwrap() {
            "kinetic_theory" => "ViscosityDilute::KineticTheory".into(),
            "collision_integral" => format!(
                "ViscosityDilute::CollisionIntegral {{ a: {}, t: {}, c: {}, molar_mass: {} }}",
                slice(&jarr(dil, "a")),
                slice(&jarr(dil, "t")),
                f(jnum(dil, "C")),
                f(jnum(dil, "molar_mass"))
            ),
            "powers_of_T" => format!(
                "ViscosityDilute::PowersOfT {{ a: {}, t: {} }}",
                slice(&jarr(dil, "a")),
                slice(&jarr(dil, "t"))
            ),
            "powers_of_Tr" => format!(
                "ViscosityDilute::PowersOfTr {{ a: {}, t: {}, t_reducing: {} }}",
                slice(&jarr(dil, "a")),
                slice(&jarr(dil, "t")),
                f(jnum(dil, "T_reducing"))
            ),
            "collision_integral_powers_of_Tstar" => format!(
                "ViscosityDilute::CollisionIntegralPowersOfTstar {{ a: {}, t: {}, c: {}, t_reducing: {} }}",
                slice(&jarr(dil, "a")),
                slice(&jarr(dil, "t")),
                f(jnum(dil, "C")),
                f(jnum(dil, "T_reducing"))
            ),
            other => panic!("unknown dilute viscosity type {other}"),
        }
    };
    writeln!(w, "            dilute: {dilute},").unwrap();

    match v.get("initial_density") {
        None => writeln!(w, "            initial_density: None,").unwrap(),
        Some(id) => {
            let body = match id.get("type").and_then(serde_json::Value::as_str).unwrap() {
                "Rainwater-Friend" => format!(
                    "ViscosityInitialDensity::RainwaterFriend {{ b: {}, t: {} }}",
                    slice(&jarr(id, "b")),
                    slice(&jarr(id, "t"))
                ),
                "empirical" => format!(
                    "ViscosityInitialDensity::Empirical {{ n: {}, d: {}, t: {}, t_reducing: {}, rhomolar_reducing: {} }}",
                    slice(&jarr(id, "n")),
                    slice(&jarr(id, "d")),
                    slice(&jarr(id, "t")),
                    f(jnum(id, "T_reducing")),
                    f(jnum(id, "rhomolar_reducing"))
                ),
                other => panic!("unknown initial_density viscosity type {other}"),
            };
            writeln!(w, "            initial_density: Some({body}),").unwrap();
        }
    }

    let ho = v
        .get("higher_order")
        .expect("structured viscosity has higher_order");
    let higher = if let Some(h) = ho.get("hardcoded") {
        format!(
            "ViscosityHigherOrder::Hardcoded {{ name: {:?} }}",
            h.as_str().unwrap()
        )
    } else {
        match ho.get("type").and_then(serde_json::Value::as_str).unwrap() {
            "modified_Batschinski_Hildebrand" => format!(
                "ViscosityHigherOrder::ModifiedBatschinskiHildebrand {{ a: {}, d1: {}, t1: {}, gamma: {}, l: {}, f: {}, d2: {}, t2: {}, g: {}, h: {}, p: {}, q: {}, t_reduce: {}, rhomolar_reduce: {} }}",
                slice(&jarr(ho, "a")),
                slice(&jarr(ho, "d1")),
                slice(&jarr(ho, "t1")),
                slice(&jarr(ho, "gamma")),
                slice(&jarr(ho, "l")),
                slice(&jarr(ho, "f")),
                slice(&jarr(ho, "d2")),
                slice(&jarr(ho, "t2")),
                slice(&jarr(ho, "g")),
                slice(&jarr(ho, "h")),
                slice(&jarr(ho, "p")),
                slice(&jarr(ho, "q")),
                f(jnum(ho, "T_reduce")),
                f(jnum(ho, "rhomolar_reduce"))
            ),
            "friction_theory" => {
                // Upstream: exactly one of Arr/Adrdr; Aii and Arrr/Aaaa optional.
                format!(
                    "ViscosityHigherOrder::FrictionTheory {{ ai: {}, aa: {}, ar: {}, aaa: {}, arr: {}, adrdr: {}, aii: {}, arrr: {}, aaaa: {}, na: {}, naa: {}, nr: {}, nrr: {}, nii: {}, nrrr: {}, naaa: {}, c1: {}, c2: {}, t_reduce: {} }}",
                    slice(&jarr(ho, "Ai")),
                    slice(&jarr(ho, "Aa")),
                    slice(&jarr(ho, "Ar")),
                    slice(&jarr(ho, "Aaa")),
                    slice(&jarr_opt(ho, "Arr")),
                    slice(&jarr_opt(ho, "Adrdr")),
                    slice(&jarr_opt(ho, "Aii")),
                    slice(&jarr_opt(ho, "Arrr")),
                    slice(&jarr_opt(ho, "Aaaa")),
                    f(jnum(ho, "Na")),
                    f(jnum(ho, "Naa")),
                    f(jnum(ho, "Nr")),
                    f(jnum(ho, "Nrr")),
                    f(jnum_opt(ho, "Nii", 0.0)),
                    f(jnum_opt(ho, "Nrrr", 0.0)),
                    f(jnum_opt(ho, "Naaa", 0.0)),
                    f(jnum(ho, "c1")),
                    f(jnum(ho, "c2")),
                    f(jnum(ho, "T_reduce"))
                )
            }
            other => panic!("unknown higher_order viscosity type {other}"),
        }
    };
    writeln!(w, "            higher_order: {higher},").unwrap();
    out
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct CubicFluidJson {
    name: String,
    #[serde(rename = "CAS")]
    cas: String,
    aliases: Vec<String>,
    #[serde(rename = "Tc")]
    tc: f64,
    pc: f64,
    rhomolarc: f64,
    acentric: f64,
    molemass: f64,
    alpha0: Vec<Alpha0Json>,
}

/// Emit `crates/rustprop-data/src/cubics.rs` from
/// `data/cubics/all_cubic_fluids.json` (upstream `CubicsLibrary` contents),
/// behind the single `cubic-fluids` feature: the whole 116-fluid table costs
/// about one HEOS fluid's worth of rodata, so per-fluid gating buys nothing.
fn emit_cubics(root: &Path) -> Result<(), String> {
    let path = root.join("data/cubics/all_cubic_fluids.json");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let fluids: Vec<CubicFluidJson> =
        serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut w = String::from(
        "//! GENERATED by rustprop-datagen from data/cubics/all_cubic_fluids.json — DO NOT EDIT.\n//! Regenerate: cargo run -p rustprop-datagen\n\n#![cfg_attr(rustfmt, rustfmt::skip)]\n#![allow(clippy::approx_constant)]\n\nuse rustprop_core::fluid::{Alpha0Term, CubicFluid};\n\n",
    );
    writeln!(
        w,
        "/// Upstream `CubicsLibrary` contents, in file order.\npub static CUBIC_FLUIDS: &[CubicFluid] = &["
    )
    .unwrap();
    for fl in &fluids {
        writeln!(w, "    CubicFluid {{").unwrap();
        writeln!(w, "        name: {:?},", fl.name).unwrap();
        writeln!(w, "        cas: {:?},", fl.cas).unwrap();
        let aliases: Vec<String> = fl.aliases.iter().map(|a| format!("{a:?}")).collect();
        writeln!(w, "        aliases: &[{}],", aliases.join(", ")).unwrap();
        writeln!(w, "        tc: {},", f(fl.tc)).unwrap();
        writeln!(w, "        pc: {},", f(fl.pc)).unwrap();
        writeln!(w, "        rhomolarc: {},", f(fl.rhomolarc)).unwrap();
        writeln!(w, "        acentric: {},", f(fl.acentric)).unwrap();
        writeln!(w, "        molemass: {},", f(fl.molemass)).unwrap();
        render_alpha0(&mut w, &fl.alpha0, "        ");
        writeln!(w, "    }},").unwrap();
    }
    writeln!(w, "];").unwrap();
    std::fs::write(root.join("crates/rustprop-data/src/cubics.rs"), w)
        .map_err(|e| e.to_string())?;
    println!(
        "generated crates/rustprop-data/src/cubics.rs ({} fluids)",
        fluids.len()
    );
    Ok(())
}

#[derive(Deserialize)]
struct IncompFluidJson {
    name: String,
    #[serde(rename = "Tmin")]
    tmin: f64,
    #[serde(rename = "Tmax")]
    tmax: f64,
    xmin: Option<f64>,
    xmax: Option<f64>,
    xid: String,
    #[serde(rename = "TminPsat")]
    tmin_psat: Option<f64>,
    #[serde(rename = "Tbase")]
    tbase: Option<f64>,
    xbase: Option<f64>,
    density: IncompBlockJson,
    specific_heat: IncompBlockJson,
    conductivity: IncompBlockJson,
    viscosity: IncompBlockJson,
    saturation_pressure: IncompBlockJson,
    #[serde(rename = "T_freeze")]
    t_freeze: IncompBlockJson,
    mass2input: IncompBlockJson,
    mole2input: IncompBlockJson,
    volume2input: IncompBlockJson,
}

#[derive(Deserialize)]
struct IncompBlockJson {
    #[serde(rename = "type")]
    block_type: String,
    coeffs: serde_json::Value,
}

/// Render one incompressible property block (upstream `parse_coefficients`:
/// the five known tags; anything else on a non-vital block is `NotSet`).
fn incomp_block(b: &IncompBlockJson, vital: bool, ctx: &str) -> String {
    let arr2d = |v: &serde_json::Value| -> String {
        let rows: Vec<String> = v
            .as_array()
            .unwrap_or_else(|| panic!("{ctx}: coeffs not 2D"))
            .iter()
            .map(|row| {
                let vals: Vec<f64> = row
                    .as_array()
                    .unwrap_or_else(|| panic!("{ctx}: ragged coeffs"))
                    .iter()
                    .map(|x| x.as_f64().unwrap())
                    .collect();
                slice(&vals).to_string()
            })
            .collect();
        format!("&[{}]", rows.join(", "))
    };
    let arr1d = |v: &serde_json::Value| -> String {
        let vals: Vec<f64> = v
            .as_array()
            .unwrap_or_else(|| panic!("{ctx}: coeffs not 1D"))
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        slice(&vals).to_string()
    };
    match b.block_type.as_str() {
        "polynomial" => format!("IncompData::Polynomial({})", arr2d(&b.coeffs)),
        "exppolynomial" => format!("IncompData::ExpPolynomial({})", arr2d(&b.coeffs)),
        "exponential" => format!("IncompData::Exponential({})", arr1d(&b.coeffs)),
        "logexponential" => format!("IncompData::LogExponential({})", arr1d(&b.coeffs)),
        "polyoffset" => format!("IncompData::PolyOffset({})", arr1d(&b.coeffs)),
        other if !vital => {
            assert_eq!(
                other, "notdefined",
                "{ctx}: unexpected non-vital tag {other:?}"
            );
            "IncompData::NotSet".to_string()
        }
        other => panic!("{ctx}: vital block has unknown type {other:?}"),
    }
}

/// Emit `crates/rustprop-data/src/incompressible.rs` from the 126 documents
/// in data/incompressible/, behind the single `incompressible-fluids`
/// feature.
fn emit_incompressibles(root: &Path) -> Result<(), String> {
    let dir = root.join("data/incompressible");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .map_err(|e| e.to_string())?
        .filter_map(|e| {
            let p = e.ok()?.path();
            (p.extension()? == "json").then_some(p)
        })
        .collect();
    files.sort();
    let mut w = String::from(
        "//! GENERATED by rustprop-datagen from data/incompressible/*.json — DO NOT EDIT.\n//! Regenerate: cargo run -p rustprop-datagen\n\n#![cfg_attr(rustfmt, rustfmt::skip)]\n#![allow(clippy::approx_constant)]\n\nuse rustprop_core::fluid::{IncompData, IncompFluid, IncompFrac};\n\n/// Upstream `JSONIncompressibleLibrary` contents, in file order.\npub static INCOMP_FLUIDS: &[IncompFluid] = &[\n",
    );
    let mut n = 0;
    for path in &files {
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let fl: IncompFluidJson =
            serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
        let xid = match fl.xid.as_str() {
            "pure" => "Pure",
            "mass" => "Mass",
            "volume" => "Volume",
            other => panic!("{}: unknown xid {other:?}", path.display()),
        };
        writeln!(w, "    IncompFluid {{").unwrap();
        writeln!(w, "        name: {:?},", fl.name).unwrap();
        writeln!(w, "        tmin: {},", f(fl.tmin)).unwrap();
        writeln!(w, "        tmax: {},", f(fl.tmax)).unwrap();
        writeln!(w, "        xmin: {},", f(fl.xmin.unwrap_or(0.0))).unwrap();
        writeln!(w, "        xmax: {},", f(fl.xmax.unwrap_or(1.0))).unwrap();
        writeln!(w, "        xid: IncompFrac::{xid},").unwrap();
        writeln!(w, "        tmin_psat: {},", f(fl.tmin_psat.unwrap_or(0.0))).unwrap();
        writeln!(w, "        tbase: {},", f(fl.tbase.unwrap_or(0.0))).unwrap();
        writeln!(w, "        xbase: {},", f(fl.xbase.unwrap_or(0.0))).unwrap();
        let name = &fl.name;
        writeln!(
            w,
            "        density: {},",
            incomp_block(&fl.density, true, name)
        )
        .unwrap();
        writeln!(
            w,
            "        specific_heat: {},",
            incomp_block(&fl.specific_heat, true, name)
        )
        .unwrap();
        writeln!(
            w,
            "        conductivity: {},",
            incomp_block(&fl.conductivity, false, name)
        )
        .unwrap();
        writeln!(
            w,
            "        viscosity: {},",
            incomp_block(&fl.viscosity, false, name)
        )
        .unwrap();
        writeln!(
            w,
            "        saturation_pressure: {},",
            incomp_block(&fl.saturation_pressure, false, name)
        )
        .unwrap();
        writeln!(
            w,
            "        t_freeze: {},",
            incomp_block(&fl.t_freeze, false, name)
        )
        .unwrap();
        writeln!(
            w,
            "        mass2input: {},",
            incomp_block(&fl.mass2input, false, name)
        )
        .unwrap();
        writeln!(
            w,
            "        mole2input: {},",
            incomp_block(&fl.mole2input, false, name)
        )
        .unwrap();
        writeln!(
            w,
            "        volume2input: {},",
            incomp_block(&fl.volume2input, false, name)
        )
        .unwrap();
        writeln!(w, "    }},").unwrap();
        n += 1;
    }
    w.push_str("];\n");
    std::fs::write(root.join("crates/rustprop-data/src/incompressible.rs"), w)
        .map_err(|e| e.to_string())?;
    println!("generated crates/rustprop-data/src/incompressible.rs ({n} fluids)");
    Ok(())
}

#[derive(Deserialize)]
struct MixPairJson {
    #[serde(rename = "CAS1")]
    cas1: String,
    #[serde(rename = "CAS2")]
    cas2: String,
    #[serde(rename = "F")]
    f: f64,
    #[serde(rename = "betaT")]
    beta_t: Option<f64>,
    #[serde(rename = "gammaT")]
    gamma_t: Option<f64>,
    #[serde(rename = "betaV")]
    beta_v: Option<f64>,
    #[serde(rename = "gammaV")]
    gamma_v: Option<f64>,
    xi: Option<f64>,
    zeta: Option<f64>,
    function: Option<String>,
}

#[derive(Deserialize)]
struct MixDepJson {
    #[serde(rename = "Name")]
    name: String,
    #[serde(rename = "type")]
    dep_type: String,
    #[serde(rename = "Npower")]
    npower: Option<usize>,
    n: Vec<f64>,
    d: Vec<f64>,
    t: Vec<f64>,
    l: Option<Vec<f64>>,
    eta: Option<Vec<f64>>,
    epsilon: Option<Vec<f64>>,
    beta: Option<Vec<f64>>,
    gamma: Option<Vec<f64>>,
}

/// Emit `crates/rustprop-data/src/mixtures.rs` from data/mixtures/, behind
/// the `mixture-data` feature. The six Lemmon xi/zeta records are converted
/// to GERG form here exactly as upstream's load-time
/// `LemmonAirHFCReducingFunction::convert_to_GERG` (using the two fluids'
/// `EOS.reducing` states, read from the runtime fluid dumps).
fn emit_mixtures(root: &Path) -> Result<(), String> {
    // Reducing states by CAS, from the runtime fluid dumps.
    let mut reduce: std::collections::HashMap<String, (f64, f64)> =
        std::collections::HashMap::new();
    for entry in std::fs::read_dir(root.join("data/coolprop-json")).map_err(|e| e.to_string())? {
        let p = entry.map_err(|e| e.to_string())?.path();
        if p.extension().map(|e| e != "json").unwrap_or(true) {
            continue;
        }
        let text = std::fs::read_to_string(&p).map_err(|e| e.to_string())?;
        let doc: serde_json::Value = serde_json::from_str(&text).map_err(|e| e.to_string())?;
        let d = &doc[0];
        let cas = d["INFO"]["CAS"].as_str().unwrap().to_string();
        let red = &d["EOS"][0]["STATES"]["reducing"];
        reduce.insert(
            cas,
            (
                red["T"].as_f64().unwrap(),
                red["rhomolar"].as_f64().unwrap(),
            ),
        );
    }

    let pairs_text = std::fs::read_to_string(root.join("data/mixtures/mixture_binary_pairs.json"))
        .map_err(|e| e.to_string())?;
    let pairs: Vec<MixPairJson> = serde_json::from_str(&pairs_text).map_err(|e| e.to_string())?;
    let deps_text =
        std::fs::read_to_string(root.join("data/mixtures/mixture_departure_functions.json"))
            .map_err(|e| e.to_string())?;
    let deps: Vec<MixDepJson> = serde_json::from_str(&deps_text).map_err(|e| e.to_string())?;

    let mut w = String::from(
        "//! GENERATED by rustprop-datagen from data/mixtures/*.json — DO NOT EDIT.\n//! Regenerate: cargo run -p rustprop-datagen\n\n#![cfg_attr(rustfmt, rustfmt::skip)]\n#![allow(clippy::approx_constant)]\n\nuse rustprop_core::fluid::{MixBinaryPair, MixDepartureFn, MixDepartureKind, PredefinedMixture};\n\n/// Upstream `mixture_binary_pairs.json`, CAS keys stored AS SHIPPED (already\n/// sorted in the source); Lemmon records converted to GERG form.\npub static MIX_BINARY_PAIRS: &[MixBinaryPair] = &[\n",
    );
    let mut n_lemmon = 0;
    for pr in &pairs {
        let (beta_t, gamma_t, beta_v, gamma_v);
        if let (Some(xi), Some(zeta)) = (pr.xi, pr.zeta) {
            // Lemmon xi/zeta -> GERG (upstream convert_to_GERG).
            let (t1, rho1) = reduce
                .get(&pr.cas1)
                .unwrap_or_else(|| panic!("no reducing state for CAS {}", pr.cas1));
            let (t2, rho2) = reduce
                .get(&pr.cas2)
                .unwrap_or_else(|| panic!("no reducing state for CAS {}", pr.cas2));
            beta_t = 1.0;
            beta_v = 1.0;
            gamma_t = (t1 + t2 + xi) / (2.0 * (t1 * t2).sqrt());
            let v1 = 1.0 / rho1;
            let v2 = 1.0 / rho2;
            let one_third = 1.0 / 3.0;
            gamma_v = (v1 + v2 + zeta) / (0.25 * (v1.powf(one_third) + v2.powf(one_third)).powi(3));
            n_lemmon += 1;
        } else {
            beta_t = pr.beta_t.expect("betaT");
            gamma_t = pr.gamma_t.expect("gammaT");
            beta_v = pr.beta_v.expect("betaV");
            gamma_v = pr.gamma_v.expect("gammaV");
        }
        let function = match (&pr.function, pr.f.abs() > f64::EPSILON) {
            (Some(name), true) => format!("Some({name:?})"),
            _ => "None".to_string(),
        };
        writeln!(
            w,
            "    MixBinaryPair {{ cas1: {:?}, cas2: {:?}, beta_t: {}, gamma_t: {}, beta_v: {}, gamma_v: {}, f: {}, function: {} }},",
            pr.cas1, pr.cas2, f(beta_t), f(gamma_t), f(beta_v), f(gamma_v), f(pr.f), function
        )
        .unwrap();
    }
    writeln!(w, "];").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "/// Upstream `mixture_departure_functions.json`, in file order."
    )
    .unwrap();
    writeln!(w, "pub static MIX_DEPARTURE_FNS: &[MixDepartureFn] = &[").unwrap();
    for dp in &deps {
        let (kind, npower) = match dp.dep_type.as_str() {
            "GERG-2008" => ("Gerg2008", dp.npower.expect("Npower")),
            "Exponential" => ("Exponential", dp.n.len()),
            "Gaussian+Exponential" => ("GaussianExponential", dp.npower.expect("Npower")),
            other => panic!("unknown departure type {other:?}"),
        };
        let arr = |v: &Option<Vec<f64>>| -> String {
            match v {
                Some(v) => slice(v).to_string(),
                None => "&[]".to_string(),
            }
        };
        writeln!(w, "    MixDepartureFn {{").unwrap();
        writeln!(w, "        name: {:?},", dp.name).unwrap();
        writeln!(w, "        kind: MixDepartureKind::{kind},").unwrap();
        writeln!(w, "        npower: {npower},").unwrap();
        writeln!(w, "        n: {},", slice(&dp.n)).unwrap();
        writeln!(w, "        d: {},", slice(&dp.d)).unwrap();
        writeln!(w, "        t: {},", slice(&dp.t)).unwrap();
        writeln!(w, "        l: {},", arr(&dp.l)).unwrap();
        writeln!(w, "        eta: {},", arr(&dp.eta)).unwrap();
        writeln!(w, "        epsilon: {},", arr(&dp.epsilon)).unwrap();
        writeln!(w, "        beta: {},", arr(&dp.beta)).unwrap();
        writeln!(w, "        gamma: {},", arr(&dp.gamma)).unwrap();
        writeln!(w, "    }},").unwrap();
    }
    writeln!(w, "];").unwrap();

    let predef_text = std::fs::read_to_string(root.join("data/mixtures/predefined_mixtures.json"))
        .map_err(|e| e.to_string())?;
    let predefs: Vec<PredefinedMixJson> =
        serde_json::from_str(&predef_text).map_err(|e| e.to_string())?;
    writeln!(w).unwrap();
    writeln!(
        w,
        "/// Upstream `predefined_mixtures.json`, in file order (154 blends)."
    )
    .unwrap();
    writeln!(w, "pub static MIX_PREDEFINED: &[PredefinedMixture] = &[").unwrap();
    for pd in &predefs {
        let fluid_list = pd
            .fluids
            .iter()
            .map(|n| format!("{n:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(
            w,
            "    PredefinedMixture {{ name: {:?}, fluids: &[{fluid_list}], mole_fractions: {} }},",
            pd.name,
            slice(&pd.mole_fractions)
        )
        .unwrap();
    }
    writeln!(w, "];").unwrap();
    std::fs::write(root.join("crates/rustprop-data/src/mixtures.rs"), w)
        .map_err(|e| e.to_string())?;
    println!(
        "generated crates/rustprop-data/src/mixtures.rs ({} pairs incl. {n_lemmon} Lemmon-converted, {} departure fns, {} predefined)",
        pairs.len(),
        deps.len(),
        predefs.len()
    );
    Ok(())
}

#[derive(serde::Deserialize)]
struct PredefinedMixJson {
    name: String,
    fluids: Vec<String>,
    mole_fractions: Vec<f64>,
}

fn emit_pcsaft(root: &Path) -> Result<(), String> {
    #[derive(serde::Deserialize)]
    struct PcsaftFluidJson {
        name: String,
        #[serde(rename = "CAS")]
        cas: String,
        aliases: Vec<String>,
        m: f64,
        sigma: f64,
        u: f64,
        molemass: f64,
        #[serde(rename = "uAB")]
        u_ab: Option<f64>,
        #[serde(rename = "volA")]
        vol_a: Option<f64>,
        #[serde(rename = "assocScheme")]
        assoc_scheme: Option<Vec<String>>,
        dipm: Option<f64>,
        dipnum: Option<f64>,
        charge: Option<f64>,
    }
    #[derive(serde::Deserialize)]
    struct PcsaftPairJson {
        #[serde(rename = "CAS1")]
        cas1: String,
        #[serde(rename = "CAS2")]
        cas2: String,
        kij: f64,
        #[serde(rename = "kijT")]
        kij_t: Option<f64>,
    }

    let fluids_text = std::fs::read_to_string(root.join("data/pcsaft/all_pcsaft_fluids.json"))
        .map_err(|e| e.to_string())?;
    let fluids: Vec<PcsaftFluidJson> =
        serde_json::from_str(&fluids_text).map_err(|e| e.to_string())?;
    let pairs_text =
        std::fs::read_to_string(root.join("data/pcsaft/mixture_binary_pairs_pcsaft.json"))
            .map_err(|e| e.to_string())?;
    let pairs: Vec<PcsaftPairJson> =
        serde_json::from_str(&pairs_text).map_err(|e| e.to_string())?;

    let mut w = String::from(
        "//! GENERATED by rustprop-datagen from data/pcsaft/*.json — DO NOT EDIT.\n//! Regenerate: cargo run -p rustprop-datagen\n\n#![cfg_attr(rustfmt, rustfmt::skip)]\n#![allow(clippy::approx_constant)]\n\nuse rustprop_core::fluid::{PcsaftBinaryPair, PcsaftFluid};\n\n/// Upstream `all_pcsaft_fluids.json` (180 fluids), in file order; optional\n/// keys default to 0/empty exactly as `cpjson::make_pcsaft_fluid`.\npub static PCSAFT_FLUIDS: &[PcsaftFluid] = &[\n",
    );
    for fl in &fluids {
        let aliases = fl
            .aliases
            .iter()
            .map(|a| format!("{a:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        let schemes = fl
            .assoc_scheme
            .as_ref()
            .map(|v| {
                v.iter()
                    .map(|a| format!("{a:?}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            })
            .unwrap_or_default();
        writeln!(
            w,
            "    PcsaftFluid {{ name: {:?}, cas: {:?}, aliases: &[{aliases}], m: {}, sigma: {}, u: {}, molemass: {}, u_ab: {}, vol_a: {}, assoc_scheme: &[{schemes}], dipm: {}, dipnum: {}, z: {} }},",
            fl.name,
            fl.cas,
            f(fl.m),
            f(fl.sigma),
            f(fl.u),
            f(fl.molemass),
            f(fl.u_ab.unwrap_or(0.0)),
            f(fl.vol_a.unwrap_or(0.0)),
            f(fl.dipm.unwrap_or(0.0)),
            f(fl.dipnum.unwrap_or(0.0)),
            f(fl.charge.unwrap_or(0.0)),
        )
        .unwrap();
    }
    writeln!(w, "];").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "/// Upstream `mixture_binary_pairs_pcsaft.json` (140 records), in file\n/// order with the CAS pair SORTED as `PCSAFTLibrary::load_from_JSON` does;\n/// absent `kijT` is 0 (the upstream getter's fallback)."
    )
    .unwrap();
    writeln!(
        w,
        "pub static PCSAFT_BINARY_PAIRS: &[PcsaftBinaryPair] = &["
    )
    .unwrap();
    for pr in &pairs {
        let (c1, c2) = if pr.cas1 <= pr.cas2 {
            (&pr.cas1, &pr.cas2)
        } else {
            (&pr.cas2, &pr.cas1)
        };
        writeln!(
            w,
            "    PcsaftBinaryPair {{ cas1: {c1:?}, cas2: {c2:?}, kij: {}, kij_t: {} }},",
            f(pr.kij),
            f(pr.kij_t.unwrap_or(0.0)),
        )
        .unwrap();
    }
    writeln!(w, "];").unwrap();
    std::fs::write(root.join("crates/rustprop-data/src/pcsaft.rs"), w)
        .map_err(|e| e.to_string())?;
    println!(
        "generated crates/rustprop-data/src/pcsaft.rs ({} fluids, {} kij pairs)",
        fluids.len(),
        pairs.len()
    );
    Ok(())
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .unwrap()
}

fn run() -> Result<(), String> {
    let root = workspace_root();
    let json_dir = root.join("data/coolprop-json");
    let out_dir = root.join("crates/rustprop-data/src/fluids");

    let mut fluids: Vec<String> = std::env::args().skip(1).collect();
    if fluids.is_empty() {
        let mut found: Vec<String> = std::fs::read_dir(&json_dir)
            .map_err(|e| format!("cannot read {}: {e}", json_dir.display()))?
            .filter_map(|entry| {
                let name = entry.ok()?.file_name().into_string().ok()?;
                name.strip_suffix(".json").map(str::to_owned)
            })
            .collect();
        found.sort();
        fluids = found;
    }
    if fluids.is_empty() {
        return Err("no fluid dumps found in data/coolprop-json/".into());
    }

    std::fs::create_dir_all(&out_dir).map_err(|e| e.to_string())?;
    let data_manifest = std::fs::read_to_string(root.join("crates/rustprop-data/Cargo.toml"))
        .map_err(|e| e.to_string())?;

    let mut mod_rs = String::from(
        "//! GENERATED by rustprop-datagen — DO NOT EDIT.\n//! One feature-gated module per fluid dump in data/coolprop-json/.\n\n#![cfg_attr(rustfmt, rustfmt::skip)]\n// The registry's cfg-gated pushes require the init-then-push shape.\n#![allow(clippy::vec_init_then_push)]\n\n",
    );
    let mut registry = String::from(
        "\n/// Every compiled-in fluid: (upstream name, data), in dump-file order.\npub fn all() -> Vec<(&'static str, &'static rustprop_core::fluid::FluidData)> {\n    #[allow(unused_mut)]\n    let mut v: Vec<(&'static str, &'static rustprop_core::fluid::FluidData)> = Vec::new();\n",
    );
    for fluid in &fluids {
        let path = json_dir.join(format!("{fluid}.json"));
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let docs: Vec<Doc> =
            serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
        let doc = docs
            .first()
            .ok_or_else(|| format!("{fluid}: empty document"))?;
        if doc.info.name != *fluid {
            return Err(format!("{fluid}: document NAME is {:?}", doc.info.name));
        }
        let module = module_name(fluid);
        if !data_manifest.contains(&format!("{module} = []")) {
            return Err(format!(
                "feature `{module} = []` missing from crates/rustprop-data/Cargo.toml"
            ));
        }
        let eos: EosJson = serde_json::from_value(
            doc.eos
                .first()
                .ok_or_else(|| format!("{fluid}: no EOS entries"))?
                .clone(),
        )
        .map_err(|e| format!("{fluid}: EOS[0]: {e}"))?;
        let rendered = emit(doc, &eos, &format!("data/coolprop-json/{fluid}.json"));
        let out_path = out_dir.join(format!("{module}.rs"));
        std::fs::write(&out_path, rendered).map_err(|e| e.to_string())?;
        writeln!(mod_rs, "#[cfg(feature = \"{module}\")]\npub mod {module};").unwrap();
        writeln!(
            registry,
            "    #[cfg(feature = \"{module}\")]\n    v.push(({:?}, &{module}::{}));",
            doc.info.name,
            module.to_uppercase()
        )
        .unwrap();
        println!("generated crates/rustprop-data/src/fluids/{module}.rs");
    }
    registry.push_str("    v\n}\n");
    mod_rs.push_str(&registry);
    std::fs::write(out_dir.join("mod.rs"), mod_rs).map_err(|e| e.to_string())?;
    emit_cubics(&root)?;
    emit_incompressibles(&root)?;
    emit_mixtures(&root)?;
    emit_pcsaft(&root)?;
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("rustprop-datagen: {e}");
            ExitCode::FAILURE
        }
    }
}
