//! Data fidelity test (PLAN.md 3.3/4.7): every ported field of each
//! generated fluid's data must equal the committed upstream JSON **exactly**
//! (bitwise — the emitter writes shortest round-trip literals, so `==` must
//! hold).
//!
//! This walker parses the JSON independently of the datagen tool (raw
//! `serde_json::Value` with its own path mapping), so a key misread in the
//! generator cannot cancel out. It is also a completeness gate: any JSON key
//! that is neither ported nor on the explicit skip list fails the test, so
//! nothing gets dropped silently.

use rustprop_core::fluid::ChebyshevInterval;
use rustprop_core::fluid::{
    Alpha0Term, AlpharTerm, Conductivity, ConductivityCritical, ConductivityDilute,
    ConductivityModel, ConductivityResidual, FluidData, SaturationAncillary, StatePoint,
    TransportModel, Viscosity, ViscosityDilute, ViscosityHigherOrder, ViscosityInitialDensity,
    ViscosityModel,
};
use serde_json::Value;
use std::path::Path;

struct Walker {
    mismatches: Vec<String>,
}

impl Walker {
    fn num(&mut self, rust: f64, json: &Value, path: &str) {
        match json.as_f64() {
            Some(j) if rust == j => {}
            Some(j) => self
                .mismatches
                .push(format!("{path}: rust {rust:?} != json {j:?}")),
            None => self
                .mismatches
                .push(format!("{path}: not a number in JSON")),
        }
    }
    fn nums(&mut self, rust: &[f64], json: &Value, path: &str) {
        let Some(arr) = json.as_array() else {
            self.mismatches
                .push(format!("{path}: not an array in JSON"));
            return;
        };
        if arr.len() != rust.len() {
            self.mismatches.push(format!(
                "{path}: rust len {} != json len {}",
                rust.len(),
                arr.len()
            ));
            return;
        }
        for (i, (r, j)) in rust.iter().zip(arr).enumerate() {
            self.num(*r, j, &format!("{path}[{i}]"));
        }
    }
    fn string(&mut self, rust: &str, json: &Value, path: &str) {
        match json.as_str() {
            Some(j) if rust == j => {}
            other => self
                .mismatches
                .push(format!("{path}: rust {rust:?} != json {other:?}")),
        }
    }
    fn boolean(&mut self, rust: bool, json: &Value, path: &str) {
        if json.as_bool() != Some(rust) {
            self.mismatches
                .push(format!("{path}: rust {rust} != json {json}"));
        }
    }
    /// Completeness gate: every key must be ported or explicitly skipped.
    fn keys(&mut self, json: &Value, path: &str, ported: &[&str], skipped: &[&str]) {
        let Some(obj) = json.as_object() else {
            self.mismatches.push(format!("{path}: not an object"));
            return;
        };
        for key in obj.keys() {
            let known = ported.contains(&key.as_str())
                || skipped.contains(&key.as_str())
                || key.ends_with("_units")
                || key.starts_with('_');
            if !known {
                self.mismatches
                    .push(format!("{path}.{key}: key neither ported nor skip-listed"));
            }
        }
        for key in ported {
            if !obj.contains_key(*key) {
                self.mismatches
                    .push(format!("{path}.{key}: ported key missing from JSON"));
            }
        }
    }
    fn state_point(&mut self, rust: &StatePoint, json: &Value, path: &str) {
        self.keys(json, path, &["T", "p", "rhomolar", "hmolar", "smolar"], &[]);
        self.num(rust.t, &json["T"], &format!("{path}.T"));
        self.num(rust.p, &json["p"], &format!("{path}.p"));
        self.num(
            rust.rhomolar,
            &json["rhomolar"],
            &format!("{path}.rhomolar"),
        );
        self.num(rust.hmolar, &json["hmolar"], &format!("{path}.hmolar"));
        self.num(rust.smolar, &json["smolar"], &format!("{path}.smolar"));
    }
    /// `sat_min_liquid`/`sat_min_vapor`: a few documents omit the caloric
    /// fields (upstream never reads them from these states); the generated
    /// data must then hold NaN.
    fn sat_min_state_point(&mut self, rust: &StatePoint, json: &Value, path: &str) {
        self.keys(json, path, &["T", "p", "rhomolar"], &["hmolar", "smolar"]);
        self.num(rust.t, &json["T"], &format!("{path}.T"));
        self.num(rust.p, &json["p"], &format!("{path}.p"));
        self.num(
            rust.rhomolar,
            &json["rhomolar"],
            &format!("{path}.rhomolar"),
        );
        for (val, key) in [(rust.hmolar, "hmolar"), (rust.smolar, "smolar")] {
            if json.get(key).is_some() {
                self.num(val, &json[key], &format!("{path}.{key}"));
            } else if !val.is_nan() {
                self.mismatches
                    .push(format!("{path}.{key}: absent in JSON but rust {val:?}"));
            }
        }
    }
    fn sat_ancillary(&mut self, rust: &SaturationAncillary, json: &Value, path: &str) {
        self.keys(
            json,
            path,
            &[
                "type",
                "n",
                "t",
                "T_r",
                "reducing_value",
                "using_tau_r",
                "Tmin",
                "Tmax",
            ],
            &["description", "max_abserror_percentage"],
        );
        self.string(rust.anc_type, &json["type"], &format!("{path}.type"));
        self.nums(rust.n, &json["n"], &format!("{path}.n"));
        self.nums(rust.t, &json["t"], &format!("{path}.t"));
        self.num(rust.t_r, &json["T_r"], &format!("{path}.T_r"));
        self.num(
            rust.reducing_value,
            &json["reducing_value"],
            &format!("{path}.reducing_value"),
        );
        self.boolean(
            rust.using_tau_r,
            &json["using_tau_r"],
            &format!("{path}.using_tau_r"),
        );
        self.num(rust.t_min, &json["Tmin"], &format!("{path}.Tmin"));
        self.num(rust.t_max, &json["Tmax"], &format!("{path}.Tmax"));
    }
}

fn check_fluid(fluid: &FluidData, json_file: &str) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/coolprop-json")
        .join(json_file);
    let text = std::fs::read_to_string(&path).unwrap();
    let docs: Value = serde_json::from_str(&text).unwrap();
    let doc = &docs[0];
    let mut w = Walker {
        mismatches: Vec::new(),
    };

    w.keys(
        doc,
        "$",
        &["INFO", "EOS", "ANCILLARIES", "STATES"],
        &["TRANSPORT"],
    );

    // INFO
    let info = &doc["INFO"];
    w.keys(
        info,
        "INFO",
        &["NAME", "CAS", "ALIASES"],
        &[
            "2DPNG_URL",
            "CHEMSPIDER_ID",
            "ENVIRONMENTAL",
            "FORMULA",
            "INCHI_KEY",
            "INCHI_STRING",
            "REFPROP_NAME",
            "SMILES",
        ],
    );
    w.string(fluid.name, &info["NAME"], "INFO.NAME");
    w.string(fluid.cas, &info["CAS"], "INFO.CAS");
    let aliases = info["ALIASES"].as_array().unwrap();
    assert_eq!(fluid.aliases.len(), aliases.len(), "INFO.ALIASES length");
    for (i, (r, j)) in fluid.aliases.iter().zip(aliases).enumerate() {
        w.string(r, j, &format!("INFO.ALIASES[{i}]"));
    }

    // EOS[0] — some documents (Ammonia) carry alternate historical EOS
    // entries after index 0; upstream's backend evaluates EOSVector[0]
    // exclusively (the `EOS()` accessor), so only EOS[0] is ported.
    let eos_arr = doc["EOS"].as_array().unwrap();
    assert!(!eos_arr.is_empty(), "at least one EOS expected");
    let eos = &eos_arr[0];
    w.keys(
        eos,
        "EOS[0]",
        &[
            "gas_constant",
            "molar_mass",
            "p_max",
            "T_max",
            "Ttriple",
            "acentric",
            "pseudo_pure",
            "STATES",
            "alpha0",
            "alphar",
            "SUPERANCILLARY",
        ],
        &["BibTeX_CP0", "BibTeX_EOS", "critical_region_splines"],
    );
    w.num(
        fluid.eos.gas_constant,
        &eos["gas_constant"],
        "EOS[0].gas_constant",
    );
    w.num(
        fluid.eos.molar_mass,
        &eos["molar_mass"],
        "EOS[0].molar_mass",
    );
    w.num(fluid.eos.p_max, &eos["p_max"], "EOS[0].p_max");
    w.num(fluid.eos.t_max, &eos["T_max"], "EOS[0].T_max");
    w.num(fluid.eos.t_triple, &eos["Ttriple"], "EOS[0].Ttriple");
    w.num(fluid.eos.acentric, &eos["acentric"], "EOS[0].acentric");
    w.boolean(
        fluid.eos.pseudo_pure,
        &eos["pseudo_pure"],
        "EOS[0].pseudo_pure",
    );

    let eos_states = &eos["STATES"];
    w.keys(
        eos_states,
        "EOS[0].STATES",
        &["reducing", "sat_min_liquid", "sat_min_vapor", "hs_anchor"],
        &[],
    );
    w.state_point(
        &fluid.eos.reducing,
        &eos_states["reducing"],
        "EOS[0].STATES.reducing",
    );
    w.sat_min_state_point(
        &fluid.eos.sat_min_liquid,
        &eos_states["sat_min_liquid"],
        "EOS[0].STATES.sat_min_liquid",
    );
    w.sat_min_state_point(
        &fluid.eos.sat_min_vapor,
        &eos_states["sat_min_vapor"],
        "EOS[0].STATES.sat_min_vapor",
    );
    w.state_point(
        &fluid.eos.hs_anchor,
        &eos_states["hs_anchor"],
        "EOS[0].STATES.hs_anchor",
    );

    // alpha0 terms
    let alpha0 = eos["alpha0"].as_array().unwrap();
    assert_eq!(fluid.eos.alpha0.len(), alpha0.len(), "alpha0 term count");
    for (i, (term, json)) in fluid.eos.alpha0.iter().zip(alpha0).enumerate() {
        let path = format!("EOS[0].alpha0[{i}]");
        match term {
            Alpha0Term::Lead { a1, a2 } => {
                w.keys(json, &path, &["type", "a1", "a2"], &[]);
                w.string(
                    "IdealGasHelmholtzLead",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.num(*a1, &json["a1"], &format!("{path}.a1"));
                w.num(*a2, &json["a2"], &format!("{path}.a2"));
            }
            Alpha0Term::LogTau { a } => {
                w.keys(json, &path, &["type", "a"], &[]);
                w.string(
                    "IdealGasHelmholtzLogTau",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.num(*a, &json["a"], &format!("{path}.a"));
            }
            Alpha0Term::PlanckEinstein { n, t } => {
                w.keys(json, &path, &["type", "n", "t"], &[]);
                w.string(
                    "IdealGasHelmholtzPlanckEinstein",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(t, &json["t"], &format!("{path}.t"));
            }
            Alpha0Term::PlanckEinsteinFunctionT { n, v, tcrit } => {
                // Some documents carry informational `R`/`T0` here; upstream's
                // parse reads neither.
                w.keys(json, &path, &["type", "n", "v", "Tcrit"], &["R", "T0"]);
                w.string(
                    "IdealGasHelmholtzPlanckEinsteinFunctionT",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(v, &json["v"], &format!("{path}.v"));
                w.num(*tcrit, &json["Tcrit"], &format!("{path}.Tcrit"));
            }
            Alpha0Term::EnthalpyEntropyOffset { a1, a2, reference } => {
                w.keys(json, &path, &["type", "a1", "a2", "reference"], &[]);
                w.string(
                    "IdealGasHelmholtzEnthalpyEntropyOffset",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.num(*a1, &json["a1"], &format!("{path}.a1"));
                w.num(*a2, &json["a2"], &format!("{path}.a2"));
                w.string(reference, &json["reference"], &format!("{path}.reference"));
            }
            Alpha0Term::Power { n, t } => {
                w.keys(json, &path, &["type", "n", "t"], &[]);
                w.string(
                    "IdealGasHelmholtzPower",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(t, &json["t"], &format!("{path}.t"));
            }
            Alpha0Term::PlanckEinsteinGeneralized { n, t, c, d } => {
                w.keys(json, &path, &["type", "n", "t", "c", "d"], &[]);
                w.string(
                    "IdealGasHelmholtzPlanckEinsteinGeneralized",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(c, &json["c"], &format!("{path}.c"));
                w.nums(d, &json["d"], &format!("{path}.d"));
            }
            Alpha0Term::Cp0Constant { cp_over_r, tc, t0 } => {
                w.keys(json, &path, &["type", "cp_over_R", "Tc", "T0"], &[]);
                w.string(
                    "IdealGasHelmholtzCP0Constant",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.num(*cp_over_r, &json["cp_over_R"], &format!("{path}.cp_over_R"));
                w.num(*tc, &json["Tc"], &format!("{path}.Tc"));
                w.num(*t0, &json["T0"], &format!("{path}.T0"));
            }
            Alpha0Term::Cp0PolyT { c, t, tc, t0 } => {
                // `R` is informational; upstream's parse does not read it.
                w.keys(json, &path, &["type", "c", "t", "Tc", "T0"], &["R"]);
                w.string(
                    "IdealGasHelmholtzCP0PolyT",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(c, &json["c"], &format!("{path}.c"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.num(*tc, &json["Tc"], &format!("{path}.Tc"));
                w.num(*t0, &json["T0"], &format!("{path}.T0"));
            }
            Alpha0Term::Cp0AlyLee { c, tc, t0 } => {
                w.keys(json, &path, &["type", "c", "Tc", "T0"], &[]);
                w.string(
                    "IdealGasHelmholtzCP0AlyLee",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(c, &json["c"], &format!("{path}.c"));
                w.num(*tc, &json["Tc"], &format!("{path}.Tc"));
                w.num(*t0, &json["T0"], &format!("{path}.T0"));
            }
        }
    }

    // alphar terms
    let alphar = eos["alphar"].as_array().unwrap();
    assert_eq!(fluid.eos.alphar.len(), alphar.len(), "alphar term count");
    for (i, (term, json)) in fluid.eos.alphar.iter().zip(alphar).enumerate() {
        let path = format!("EOS[0].alphar[{i}]");
        match term {
            AlpharTerm::Power { n, d, t, l } => {
                w.keys(json, &path, &["type", "n", "d", "t", "l"], &[]);
                w.string(
                    "ResidualHelmholtzPower",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(d, &json["d"], &format!("{path}.d"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(l, &json["l"], &format!("{path}.l"));
            }
            AlpharTerm::Gaussian {
                n,
                d,
                t,
                eta,
                beta,
                gamma,
                epsilon,
            } => {
                w.keys(
                    json,
                    &path,
                    &["type", "n", "d", "t", "eta", "beta", "gamma", "epsilon"],
                    &[],
                );
                w.string(
                    "ResidualHelmholtzGaussian",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(d, &json["d"], &format!("{path}.d"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(eta, &json["eta"], &format!("{path}.eta"));
                w.nums(beta, &json["beta"], &format!("{path}.beta"));
                w.nums(gamma, &json["gamma"], &format!("{path}.gamma"));
                w.nums(epsilon, &json["epsilon"], &format!("{path}.epsilon"));
            }
            AlpharTerm::NonAnalytic {
                n,
                a,
                b,
                beta,
                big_a,
                big_b,
                big_c,
                big_d,
            } => {
                w.keys(
                    json,
                    &path,
                    &["type", "n", "a", "b", "beta", "A", "B", "C", "D"],
                    &[],
                );
                w.string(
                    "ResidualHelmholtzNonAnalytic",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(a, &json["a"], &format!("{path}.a"));
                w.nums(b, &json["b"], &format!("{path}.b"));
                w.nums(beta, &json["beta"], &format!("{path}.beta"));
                w.nums(big_a, &json["A"], &format!("{path}.A"));
                w.nums(big_b, &json["B"], &format!("{path}.B"));
                w.nums(big_c, &json["C"], &format!("{path}.C"));
                w.nums(big_d, &json["D"], &format!("{path}.D"));
            }
            AlpharTerm::Exponential { n, d, t, g, l } => {
                w.keys(json, &path, &["type", "n", "d", "t", "g", "l"], &[]);
                w.string(
                    "ResidualHelmholtzExponential",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(d, &json["d"], &format!("{path}.d"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(g, &json["g"], &format!("{path}.g"));
                w.nums(l, &json["l"], &format!("{path}.l"));
            }
            AlpharTerm::DoubleExponential {
                n,
                d,
                t,
                gd,
                ld,
                gt,
                lt,
            } => {
                w.keys(
                    json,
                    &path,
                    &["type", "n", "d", "t", "gd", "ld", "gt", "lt"],
                    &[],
                );
                w.string(
                    "ResidualHelmholtzDoubleExponential",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(d, &json["d"], &format!("{path}.d"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(gd, &json["gd"], &format!("{path}.gd"));
                w.nums(ld, &json["ld"], &format!("{path}.ld"));
                w.nums(gt, &json["gt"], &format!("{path}.gt"));
                w.nums(lt, &json["lt"], &format!("{path}.lt"));
            }
            AlpharTerm::Lemmon2005 { n, d, t, l, m } => {
                w.keys(json, &path, &["type", "n", "d", "t", "l", "m"], &[]);
                w.string(
                    "ResidualHelmholtzLemmon2005",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(d, &json["d"], &format!("{path}.d"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(l, &json["l"], &format!("{path}.l"));
                w.nums(m, &json["m"], &format!("{path}.m"));
            }
            AlpharTerm::GaoB {
                n,
                t,
                d,
                eta,
                beta,
                gamma,
                epsilon,
                b,
            } => {
                w.keys(
                    json,
                    &path,
                    &[
                        "type", "n", "t", "d", "eta", "beta", "gamma", "epsilon", "b",
                    ],
                    &[],
                );
                w.string(
                    "ResidualHelmholtzGaoB",
                    &json["type"],
                    &format!("{path}.type"),
                );
                w.nums(n, &json["n"], &format!("{path}.n"));
                w.nums(t, &json["t"], &format!("{path}.t"));
                w.nums(d, &json["d"], &format!("{path}.d"));
                w.nums(eta, &json["eta"], &format!("{path}.eta"));
                w.nums(beta, &json["beta"], &format!("{path}.beta"));
                w.nums(gamma, &json["gamma"], &format!("{path}.gamma"));
                w.nums(epsilon, &json["epsilon"], &format!("{path}.epsilon"));
                w.nums(b, &json["b"], &format!("{path}.b"));
            }
        }
    }

    // SUPERANCILLARY
    let sa_json = &eos["SUPERANCILLARY"];
    let sa = fluid
        .eos
        .superancillary
        .as_ref()
        .expect("every ported fluid has a superancillary");
    w.keys(
        sa_json,
        "EOS[0].SUPERANCILLARY",
        &[
            "jexpansions_p",
            "jexpansions_rhoL",
            "jexpansions_rhoV",
            "meta",
            "check_points",
        ],
        &["crit_anc", "source_eos_hash"],
    );
    let cheb_sets: [(&str, &[ChebyshevInterval]); 3] = [
        ("jexpansions_p", sa.p),
        ("jexpansions_rhoL", sa.rho_l),
        ("jexpansions_rhoV", sa.rho_v),
    ];
    for (name, rust_side) in cheb_sets {
        let arr = sa_json[name].as_array().unwrap();
        assert_eq!(rust_side.len(), arr.len(), "{name} interval count");
        for (i, (r, j)) in rust_side.iter().zip(arr).enumerate() {
            let path = format!("EOS[0].SUPERANCILLARY.{name}[{i}]");
            w.keys(j, &path, &["xmin", "xmax", "coef"], &[]);
            w.num(r.xmin, &j["xmin"], &format!("{path}.xmin"));
            w.num(r.xmax, &j["xmax"], &format!("{path}.xmax"));
            w.nums(r.coef, &j["coef"], &format!("{path}.coef"));
        }
    }
    let meta = &sa_json["meta"];
    w.keys(
        meta,
        "EOS[0].SUPERANCILLARY.meta",
        &["Tcrittrue / K", "rhocrittrue / mol/m^3"],
        &[
            "BrhoL / mol/m^3",
            "BrhoV / mol/m^3",
            "Tcrit / K",
            "Treducing / K",
            "Ttriple / K",
            "gas_constant / J/mol/K",
        ],
    );
    w.num(sa.t_crit_num, &meta["Tcrittrue / K"], "meta.Tcrittrue");
    w.num(
        sa.rho_crit_num,
        &meta["rhocrittrue / mol/m^3"],
        "meta.rhocrittrue",
    );
    let cps = sa_json["check_points"].as_array().unwrap();
    assert_eq!(sa.check_points.len(), cps.len(), "check_points count");
    for (i, (r, j)) in sa.check_points.iter().zip(cps).enumerate() {
        let path = format!("EOS[0].SUPERANCILLARY.check_points[{i}]");
        w.keys(
            j,
            &path,
            &[
                "T / K",
                "p(mp) / Pa",
                "rho'(mp) / mol/m^3",
                "rho''(mp) / mol/m^3",
                "p(SA)/p(mp)",
                "rho'(SA)/rho'(mp)",
                "rho''(SA)/rho''(mp)",
            ],
            &[],
        );
        w.num(r.t, &j["T / K"], &format!("{path}.T"));
        w.num(r.p, &j["p(mp) / Pa"], &format!("{path}.p"));
        w.num(r.rho_l, &j["rho'(mp) / mol/m^3"], &format!("{path}.rhoL"));
        w.num(r.rho_v, &j["rho''(mp) / mol/m^3"], &format!("{path}.rhoV"));
        w.num(r.p_ratio, &j["p(SA)/p(mp)"], &format!("{path}.p_ratio"));
        w.num(
            r.rho_l_ratio,
            &j["rho'(SA)/rho'(mp)"],
            &format!("{path}.rhoL_ratio"),
        );
        w.num(
            r.rho_v_ratio,
            &j["rho''(SA)/rho''(mp)"],
            &format!("{path}.rhoV_ratio"),
        );
    }

    // ANCILLARIES
    let anc = &doc["ANCILLARIES"];
    w.keys(
        anc,
        "ANCILLARIES",
        &["pS", "rhoL", "rhoV"],
        &["hL", "hLV", "sL", "sLV", "melting_line", "surface_tension"],
    );
    w.sat_ancillary(&fluid.ancillaries.p_s, &anc["pS"], "ANCILLARIES.pS");
    w.sat_ancillary(&fluid.ancillaries.rho_l, &anc["rhoL"], "ANCILLARIES.rhoL");
    w.sat_ancillary(&fluid.ancillaries.rho_v, &anc["rhoV"], "ANCILLARIES.rhoV");
    // surface_tension (Phase 6.2): ported when present.
    match (
        &fluid.ancillaries.surface_tension,
        anc.get("surface_tension"),
    ) {
        (Some(st), Some(json)) => {
            let path = "ANCILLARIES.surface_tension";
            w.keys(json, path, &["a", "n", "Tc"], &["BibTeX", "description"]);
            w.nums(st.a, &json["a"], &format!("{path}.a"));
            w.nums(st.n, &json["n"], &format!("{path}.n"));
            w.num(st.tc, &json["Tc"], &format!("{path}.Tc"));
        }
        (None, None) => {}
        (rust, json) => w.mismatches.push(format!(
            "ANCILLARIES.surface_tension: rust present={} json present={}",
            rust.is_some(),
            json.is_some()
        )),
    }

    // STATES
    let states = &doc["STATES"];
    w.keys(
        states,
        "STATES",
        &["critical", "triple_liquid", "triple_vapor"],
        &[],
    );
    w.state_point(
        &fluid.states.critical,
        &states["critical"],
        "STATES.critical",
    );
    w.state_point(
        &fluid.states.triple_liquid,
        &states["triple_liquid"],
        "STATES.triple_liquid",
    );
    w.state_point(
        &fluid.states.triple_vapor,
        &states["triple_vapor"],
        "STATES.triple_vapor",
    );

    // TRANSPORT (6.1, viscosity slice): mirror datagen's classification.
    check_transport(&mut w, fluid, doc.get("TRANSPORT"));

    assert!(
        w.mismatches.is_empty(),
        "{} fidelity mismatches:\n{}",
        w.mismatches.len(),
        w.mismatches.join("\n")
    );
}

/// Mirror of datagen's TRANSPORT classification. Per-property tri-state:
/// key absent -> `Absent`; unported class (top-level type/hardcoded, or a
/// rhosr list) -> `Unported`; structured -> full field comparison.
fn check_transport(w: &mut Walker, fluid: &FluidData, json: Option<&Value>) {
    let Some(tr) = json else {
        if fluid.transport.is_some() {
            w.mismatches
                .push("TRANSPORT: absent in JSON but rust Some".into());
        }
        return;
    };
    let Some(rust_tr) = &fluid.transport else {
        w.mismatches
            .push("TRANSPORT: present in JSON but rust None".into());
        return;
    };
    check_slot(
        w,
        &rust_tr.viscosity,
        tr.get("viscosity"),
        "TRANSPORT.viscosity",
        &|w, m, v, path| match m {
            ViscosityModel::Structured(rv) => check_viscosity(w, rv, v, path),
            ViscosityModel::Hardcoded { name } => {
                w.string(name, &v["hardcoded"], &format!("{path}.hardcoded"));
            }
            ViscosityModel::Chung {
                rhomolar_critical,
                acentric,
                molar_mass,
                t_critical,
                dipole_moment_d,
                kappa,
            } => {
                w.string("Chung", &v["type"], &format!("{path}.type"));
                w.num(
                    *rhomolar_critical,
                    &v["rhomolar_critical"],
                    &format!("{path}.rhomolar_critical"),
                );
                w.num(*acentric, &v["acentric"], &format!("{path}.acentric"));
                w.num(*molar_mass, &v["molar_mass"], &format!("{path}.molar_mass"));
                w.num(*t_critical, &v["T_critical"], &format!("{path}.T_critical"));
                w.num(
                    *dipole_moment_d,
                    &v["dipole_moment_D"],
                    &format!("{path}.dipole_moment_D"),
                );
                w.num(*kappa, &v["kappa"], &format!("{path}.kappa"));
            }
            ViscosityModel::RhosrCs {
                c,
                c_liq,
                c_vap,
                rhosr_critical,
                x_crossover,
            } => {
                w.string("rhosr-CS", &v["type"], &format!("{path}.type"));
                w.num(*c, &v["C"], &format!("{path}.C"));
                w.nums(c_liq, &v["c_liq"], &format!("{path}.c_liq"));
                w.nums(c_vap, &v["c_vap"], &format!("{path}.c_vap"));
                w.num(
                    *rhosr_critical,
                    &v["rhosr_critical"],
                    &format!("{path}.rhosr_critical"),
                );
                w.num(
                    *x_crossover,
                    &v["x_crossover"],
                    &format!("{path}.x_crossover"),
                );
            }
        },
    );
    check_slot(
        w,
        &rust_tr.conductivity,
        tr.get("conductivity"),
        "TRANSPORT.conductivity",
        &|w, m, c, path| match m {
            ConductivityModel::Structured(rc) => check_conductivity(w, rc, c, path),
            ConductivityModel::Hardcoded { name } => {
                w.string(name, &c["hardcoded"], &format!("{path}.hardcoded"));
            }
        },
    );
}

/// Classify one property slot exactly as datagen does and dispatch.
fn check_slot<T>(
    w: &mut Walker,
    rust: &TransportModel<T>,
    json: Option<&Value>,
    path: &str,
    check: &dyn Fn(&mut Walker, &T, &Value, &str),
) {
    // Upstream parse: an array uses its FIRST entry.
    let json = json.map(|v| {
        if v.is_array() {
            &v.as_array().unwrap()[0]
        } else {
            v
        }
    });
    let class = match json {
        None => "absent",
        Some(v) => {
            let o = v.as_object().unwrap();
            if o.contains_key("hardcoded") {
                "model" // fully-hardcoded per-fluid formulation
            } else {
                match o.get("type").and_then(Value::as_str) {
                    Some("Chung" | "rhosr-CS") => "model",
                    Some(_) => "unported", // ECS
                    None => "model",       // structured
                }
            }
        }
    };
    match (rust, class) {
        (TransportModel::Absent, "absent") | (TransportModel::Unported, "unported") => {}
        (TransportModel::Model(m), "model") => check(w, m, json.unwrap(), path),
        (rust, class) => {
            let rust_class = match rust {
                TransportModel::Absent => "absent",
                TransportModel::Unported => "unported",
                TransportModel::Model(_) => "structured",
            };
            w.mismatches
                .push(format!("{path}: rust {rust_class} vs json {class}"));
        }
    }
}

fn check_viscosity(w: &mut Walker, rv: &Viscosity, v: &Value, path: &str) {
    // epsilon_over_k / sigma_eta: NaN encodes "absent".
    for (val, key) in [
        (rv.epsilon_over_k, "epsilon_over_k"),
        (rv.sigma_eta, "sigma_eta"),
    ] {
        match v.get(key) {
            Some(j) => w.num(val, j, &format!("{path}.{key}")),
            None if val.is_nan() => {}
            None => w
                .mismatches
                .push(format!("{path}.{key}: absent in JSON but rust {val:?}")),
        }
    }
    check_viscosity_dilute(w, &rv.dilute, &v["dilute"], &format!("{path}.dilute"));
    match (&rv.initial_density, v.get("initial_density")) {
        (None, None) => {}
        (Some(id), Some(j)) => {
            check_viscosity_initial(w, id, j, &format!("{path}.initial_density"));
        }
        (rust, json) => w.mismatches.push(format!(
            "{path}.initial_density: rust present={} json present={}",
            rust.is_some(),
            json.is_some()
        )),
    }
    check_viscosity_higher(
        w,
        &rv.higher_order,
        &v["higher_order"],
        &format!("{path}.higher_order"),
    );
}

fn check_conductivity(w: &mut Walker, rc: &Conductivity, c: &Value, path: &str) {
    let dp = format!("{path}.dilute");
    match &rc.dilute {
        ConductivityDilute::RatioOfPolynomials {
            a,
            n,
            b,
            m,
            t_reducing,
        } => {
            let j = &c["dilute"];
            w.string("ratio_of_polynomials", &j["type"], &format!("{dp}.type"));
            w.nums(a, &j["A"], &format!("{dp}.A"));
            w.nums(n, &j["n"], &format!("{dp}.n"));
            w.nums(b, &j["B"], &format!("{dp}.B"));
            w.nums(m, &j["m"], &format!("{dp}.m"));
            w.num(*t_reducing, &j["T_reducing"], &format!("{dp}.T_reducing"));
        }
        ConductivityDilute::Eta0AndPoly { a, t } => {
            let j = &c["dilute"];
            w.string("eta0_and_poly", &j["type"], &format!("{dp}.type"));
            w.nums(a, &j["A"], &format!("{dp}.A"));
            w.nums(t, &j["t"], &format!("{dp}.t"));
        }
        ConductivityDilute::Hardcoded { name } => {
            w.string(name, &c["dilute"]["hardcoded"], &format!("{dp}.hardcoded"));
        }
    }
    let rp = format!("{path}.residual");
    match &rc.residual {
        ConductivityResidual::Polynomial {
            b,
            t,
            d,
            t_reducing,
            rhomass_reducing,
        } => {
            let j = &c["residual"];
            w.string("polynomial", &j["type"], &format!("{rp}.type"));
            w.nums(b, &j["B"], &format!("{rp}.B"));
            w.nums(t, &j["t"], &format!("{rp}.t"));
            w.nums(d, &j["d"], &format!("{rp}.d"));
            w.num(*t_reducing, &j["T_reducing"], &format!("{rp}.T_reducing"));
            w.num(
                *rhomass_reducing,
                &j["rhomass_reducing"],
                &format!("{rp}.rhomass_reducing"),
            );
        }
        ConductivityResidual::PolynomialAndExponential { a, t, d, gamma, l } => {
            let j = &c["residual"];
            w.string(
                "polynomial_and_exponential",
                &j["type"],
                &format!("{rp}.type"),
            );
            w.nums(a, &j["A"], &format!("{rp}.A"));
            w.nums(t, &j["t"], &format!("{rp}.t"));
            w.nums(d, &j["d"], &format!("{rp}.d"));
            w.nums(gamma, &j["gamma"], &format!("{rp}.gamma"));
            w.nums(l, &j["l"], &format!("{rp}.l"));
        }
    }
    let cp = format!("{path}.critical");
    match (&rc.critical, c.get("critical")) {
        (None, None) => {}
        (
            Some(ConductivityCritical::SimplifiedOlchowySengers {
                k,
                r0,
                gamma,
                nu,
                big_gamma,
                zeta0,
                qd,
                t_ref,
            }),
            Some(j),
        ) => {
            w.string(
                "simplified_Olchowy_Sengers",
                &j["type"],
                &format!("{cp}.type"),
            );
            // Absent keys carry upstream's defaults in the generated data.
            for (val, key, default) in [
                (*k, "k", 1.380_648_8e-23),
                (*r0, "R0", 1.03),
                (*gamma, "gamma", 1.239),
                (*nu, "nu", 0.63),
                (*big_gamma, "GAMMA", 0.0496),
                (*zeta0, "zeta0", 1.94e-10),
                (*qd, "qD", 2e9),
            ] {
                match j.get(key) {
                    Some(jv) => w.num(val, jv, &format!("{cp}.{key}")),
                    None if val == default => {}
                    None => w.mismatches.push(format!(
                        "{cp}.{key}: absent in JSON but rust {val:?} != default {default:?}"
                    )),
                }
            }
            match j.get("T_ref") {
                Some(jv) => w.num(*t_ref, jv, &format!("{cp}.T_ref")),
                None if t_ref.is_nan() => {}
                None => w
                    .mismatches
                    .push(format!("{cp}.T_ref: absent in JSON but rust {t_ref:?}")),
            }
        }
        (Some(ConductivityCritical::Hardcoded { name }), Some(j)) => {
            w.string(name, &j["hardcoded"], &format!("{cp}.hardcoded"));
        }
        (rust, json) => w.mismatches.push(format!(
            "{cp}: rust present={} json present={}",
            rust.is_some(),
            json.is_some()
        )),
    }
}

fn check_viscosity_dilute(w: &mut Walker, rust: &ViscosityDilute, json: &Value, path: &str) {
    match rust {
        ViscosityDilute::KineticTheory => {
            w.string("kinetic_theory", &json["type"], &format!("{path}.type"));
        }
        ViscosityDilute::CollisionIntegral {
            a,
            t,
            c,
            molar_mass,
        } => {
            w.string("collision_integral", &json["type"], &format!("{path}.type"));
            w.nums(a, &json["a"], &format!("{path}.a"));
            w.nums(t, &json["t"], &format!("{path}.t"));
            w.num(*c, &json["C"], &format!("{path}.C"));
            w.num(
                *molar_mass,
                &json["molar_mass"],
                &format!("{path}.molar_mass"),
            );
        }
        ViscosityDilute::PowersOfT { a, t } => {
            w.string("powers_of_T", &json["type"], &format!("{path}.type"));
            w.nums(a, &json["a"], &format!("{path}.a"));
            w.nums(t, &json["t"], &format!("{path}.t"));
        }
        ViscosityDilute::PowersOfTr { a, t, t_reducing } => {
            w.string("powers_of_Tr", &json["type"], &format!("{path}.type"));
            w.nums(a, &json["a"], &format!("{path}.a"));
            w.nums(t, &json["t"], &format!("{path}.t"));
            w.num(
                *t_reducing,
                &json["T_reducing"],
                &format!("{path}.T_reducing"),
            );
        }
        ViscosityDilute::CollisionIntegralPowersOfTstar {
            a,
            t,
            c,
            t_reducing,
        } => {
            w.string(
                "collision_integral_powers_of_Tstar",
                &json["type"],
                &format!("{path}.type"),
            );
            w.nums(a, &json["a"], &format!("{path}.a"));
            w.nums(t, &json["t"], &format!("{path}.t"));
            w.num(*c, &json["C"], &format!("{path}.C"));
            w.num(
                *t_reducing,
                &json["T_reducing"],
                &format!("{path}.T_reducing"),
            );
        }
        ViscosityDilute::Hardcoded { name } => {
            w.string(name, &json["hardcoded"], &format!("{path}.hardcoded"));
        }
    }
}

fn check_viscosity_initial(
    w: &mut Walker,
    rust: &ViscosityInitialDensity,
    json: &Value,
    path: &str,
) {
    match rust {
        ViscosityInitialDensity::RainwaterFriend { b, t } => {
            w.string("Rainwater-Friend", &json["type"], &format!("{path}.type"));
            w.nums(b, &json["b"], &format!("{path}.b"));
            w.nums(t, &json["t"], &format!("{path}.t"));
        }
        ViscosityInitialDensity::Empirical {
            n,
            d,
            t,
            t_reducing,
            rhomolar_reducing,
        } => {
            w.string("empirical", &json["type"], &format!("{path}.type"));
            w.nums(n, &json["n"], &format!("{path}.n"));
            w.nums(d, &json["d"], &format!("{path}.d"));
            w.nums(t, &json["t"], &format!("{path}.t"));
            w.num(
                *t_reducing,
                &json["T_reducing"],
                &format!("{path}.T_reducing"),
            );
            w.num(
                *rhomolar_reducing,
                &json["rhomolar_reducing"],
                &format!("{path}.rhomolar_reducing"),
            );
        }
    }
}

#[allow(clippy::too_many_lines)]
fn check_viscosity_higher(w: &mut Walker, rust: &ViscosityHigherOrder, json: &Value, path: &str) {
    match rust {
        ViscosityHigherOrder::ModifiedBatschinskiHildebrand {
            a,
            d1,
            t1,
            gamma,
            l,
            f,
            d2,
            t2,
            g,
            h,
            p,
            q,
            t_reduce,
            rhomolar_reduce,
        } => {
            w.string(
                "modified_Batschinski_Hildebrand",
                &json["type"],
                &format!("{path}.type"),
            );
            w.nums(a, &json["a"], &format!("{path}.a"));
            w.nums(d1, &json["d1"], &format!("{path}.d1"));
            w.nums(t1, &json["t1"], &format!("{path}.t1"));
            w.nums(gamma, &json["gamma"], &format!("{path}.gamma"));
            w.nums(l, &json["l"], &format!("{path}.l"));
            w.nums(f, &json["f"], &format!("{path}.f"));
            w.nums(d2, &json["d2"], &format!("{path}.d2"));
            w.nums(t2, &json["t2"], &format!("{path}.t2"));
            w.nums(g, &json["g"], &format!("{path}.g"));
            w.nums(h, &json["h"], &format!("{path}.h"));
            w.nums(p, &json["p"], &format!("{path}.p"));
            w.nums(q, &json["q"], &format!("{path}.q"));
            w.num(*t_reduce, &json["T_reduce"], &format!("{path}.T_reduce"));
            w.num(
                *rhomolar_reduce,
                &json["rhomolar_reduce"],
                &format!("{path}.rhomolar_reduce"),
            );
        }
        ViscosityHigherOrder::FrictionTheory {
            ai,
            aa,
            ar,
            aaa,
            arr,
            adrdr,
            aii,
            arrr,
            aaaa,
            na,
            naa,
            nr,
            nrr,
            nii,
            nrrr,
            naaa,
            c1,
            c2,
            t_reduce,
        } => {
            w.string("friction_theory", &json["type"], &format!("{path}.type"));
            w.nums(ai, &json["Ai"], &format!("{path}.Ai"));
            w.nums(aa, &json["Aa"], &format!("{path}.Aa"));
            w.nums(ar, &json["Ar"], &format!("{path}.Ar"));
            w.nums(aaa, &json["Aaa"], &format!("{path}.Aaa"));
            for (vals, key) in [
                (arr, "Arr"),
                (adrdr, "Adrdr"),
                (aii, "Aii"),
                (arrr, "Arrr"),
                (aaaa, "Aaaa"),
            ] {
                match json.get(key) {
                    Some(j) => w.nums(vals, j, &format!("{path}.{key}")),
                    None if vals.is_empty() => {}
                    None => w
                        .mismatches
                        .push(format!("{path}.{key}: absent in JSON but rust non-empty")),
                }
            }
            w.num(*na, &json["Na"], &format!("{path}.Na"));
            w.num(*naa, &json["Naa"], &format!("{path}.Naa"));
            w.num(*nr, &json["Nr"], &format!("{path}.Nr"));
            w.num(*nrr, &json["Nrr"], &format!("{path}.Nrr"));
            for (val, key) in [(*nii, "Nii"), (*nrrr, "Nrrr"), (*naaa, "Naaa")] {
                match json.get(key) {
                    Some(j) => w.num(val, j, &format!("{path}.{key}")),
                    None if val == 0.0 => {}
                    None => w
                        .mismatches
                        .push(format!("{path}.{key}: absent in JSON but rust {val:?}")),
                }
            }
            w.num(*c1, &json["c1"], &format!("{path}.c1"));
            w.num(*c2, &json["c2"], &format!("{path}.c2"));
            w.num(*t_reduce, &json["T_reduce"], &format!("{path}.T_reduce"));
        }
        ViscosityHigherOrder::Hardcoded { name } => {
            w.string(name, &json["hardcoded"], &format!("{path}.hardcoded"));
        }
    }
}

#[test]
fn every_fluid_data_matches_upstream_json_exactly() {
    let fluids = rustprop_data::fluids::all();
    assert_eq!(fluids.len(), 130, "all fluid features enabled");
    for (name, data) in fluids {
        check_fluid(data, &format!("{name}.json"));
    }
}
