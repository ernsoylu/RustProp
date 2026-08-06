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
}

#[derive(Deserialize)]
struct InfoJson {
    #[serde(rename = "NAME")]
    name: String,
    #[serde(rename = "CAS")]
    cas: String,
    #[serde(rename = "ALIASES")]
    aliases: Vec<String>,
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
    #[serde(rename = "pS")]
    p_s: SatAncJson,
    #[serde(rename = "rhoL")]
    rho_l: SatAncJson,
    #[serde(rename = "rhoV")]
    rho_v: SatAncJson,
    surface_tension: Option<SurfTensJson>,
}

#[derive(Deserialize)]
struct SurfTensJson {
    a: Vec<f64>,
    n: Vec<f64>,
    #[serde(rename = "Tc")]
    tc: f64,
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
/// emits as the `f64::NAN` path expression.
fn f(v: f64) -> String {
    if v.is_nan() {
        return "f64::NAN".into();
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
    writeln!(
        w,
        "use rustprop_core::fluid::{{Alpha0Term, AlpharTerm, Ancillaries, ChebyshevInterval, Eos, FluidData, SaturationAncillary, StatePoint, States, SuperAncCheckPoint, SuperAncillaryData{surf}}};"
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(w, "pub static {static_name}: FluidData = FluidData {{").unwrap();
    writeln!(w, "    name: {:?},", doc.info.name).unwrap();
    writeln!(w, "    cas: {:?},", doc.info.cas).unwrap();
    let aliases: Vec<String> = doc.info.aliases.iter().map(|a| format!("{a:?}")).collect();
    writeln!(w, "    aliases: &[{}],", aliases.join(", ")).unwrap();
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
    writeln!(w, "        alpha0: &[").unwrap();
    for term in &eos.alpha0 {
        match term {
            Alpha0Json::Lead { a1, a2 } => {
                writeln!(
                    w,
                    "            Alpha0Term::Lead {{ a1: {}, a2: {} }},",
                    f(*a1),
                    f(*a2)
                )
                .unwrap();
            }
            Alpha0Json::LogTau { a } => {
                writeln!(w, "            Alpha0Term::LogTau {{ a: {} }},", f(*a)).unwrap();
            }
            Alpha0Json::PlanckEinstein { n, t } => {
                writeln!(
                    w,
                    "            Alpha0Term::PlanckEinstein {{ n: {}, t: {} }},",
                    slice(n),
                    slice(t)
                )
                .unwrap();
            }
            Alpha0Json::PlanckEinsteinFunctionT { n, v, tcrit } => {
                writeln!(
                    w,
                    "            Alpha0Term::PlanckEinsteinFunctionT {{ n: {}, v: {}, tcrit: {} }},",
                    slice(n),
                    slice(v),
                    f(*tcrit)
                )
                .unwrap();
            }
            Alpha0Json::EnthalpyEntropyOffset { a1, a2, reference } => {
                writeln!(
                    w,
                    "            Alpha0Term::EnthalpyEntropyOffset {{ a1: {}, a2: {}, reference: {:?} }},",
                    f(*a1),
                    f(*a2),
                    reference
                )
                .unwrap();
            }
            Alpha0Json::Power { n, t } => {
                writeln!(
                    w,
                    "            Alpha0Term::Power {{ n: {}, t: {} }},",
                    slice(n),
                    slice(t)
                )
                .unwrap();
            }
            Alpha0Json::PlanckEinsteinGeneralized { n, t, c, d } => {
                writeln!(
                    w,
                    "            Alpha0Term::PlanckEinsteinGeneralized {{ n: {}, t: {}, c: {}, d: {} }},",
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
                    "            Alpha0Term::Cp0Constant {{ cp_over_r: {}, tc: {}, t0: {} }},",
                    f(*cp_over_r),
                    f(*tc),
                    f(*t0)
                )
                .unwrap();
            }
            Alpha0Json::Cp0PolyT { c, t, tc, t0 } => {
                writeln!(
                    w,
                    "            Alpha0Term::Cp0PolyT {{ c: {}, t: {}, tc: {}, t0: {} }},",
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
                    "            Alpha0Term::Cp0AlyLee {{ c: {}, tc: {}, t0: {} }},",
                    slice(c),
                    f(*tc),
                    f(*t0)
                )
                .unwrap();
            }
        }
    }
    writeln!(w, "        ],").unwrap();
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
    writeln!(
        w,
        "        p_s: {},",
        sat_anc(&doc.ancillaries.p_s, "        ")
    )
    .unwrap();
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
    writeln!(w, "}};").unwrap();
    out
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

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
