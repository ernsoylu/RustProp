//! Data fidelity test (PLAN.md 3.3): every ported field of the generated
//! Water data must equal the committed upstream JSON **exactly** (bitwise —
//! the emitter writes shortest round-trip literals, so `==` must hold).
//!
//! This walker parses the JSON independently of the datagen tool (raw
//! `serde_json::Value` with its own path mapping), so a key misread in the
//! generator cannot cancel out. It is also a completeness gate: any JSON key
//! that is neither ported nor on the explicit skip list fails the test, so
//! nothing gets dropped silently.

use rustprop_core::fluid::ChebyshevInterval;
use rustprop_core::fluid::{Alpha0Term, AlpharTerm, SaturationAncillary, StatePoint};
use rustprop_data::fluids::water::WATER;
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

#[test]
fn water_data_matches_upstream_json_exactly() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/coolprop-json/Water.json");
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
    w.string(WATER.name, &info["NAME"], "INFO.NAME");
    w.string(WATER.cas, &info["CAS"], "INFO.CAS");
    let aliases = info["ALIASES"].as_array().unwrap();
    assert_eq!(WATER.aliases.len(), aliases.len(), "INFO.ALIASES length");
    for (i, (r, j)) in WATER.aliases.iter().zip(aliases).enumerate() {
        w.string(r, j, &format!("INFO.ALIASES[{i}]"));
    }

    // EOS[0]
    let eos_arr = doc["EOS"].as_array().unwrap();
    assert_eq!(eos_arr.len(), 1, "exactly one EOS expected");
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
        WATER.eos.gas_constant,
        &eos["gas_constant"],
        "EOS[0].gas_constant",
    );
    w.num(
        WATER.eos.molar_mass,
        &eos["molar_mass"],
        "EOS[0].molar_mass",
    );
    w.num(WATER.eos.p_max, &eos["p_max"], "EOS[0].p_max");
    w.num(WATER.eos.t_max, &eos["T_max"], "EOS[0].T_max");
    w.num(WATER.eos.t_triple, &eos["Ttriple"], "EOS[0].Ttriple");
    w.num(WATER.eos.acentric, &eos["acentric"], "EOS[0].acentric");
    w.boolean(
        WATER.eos.pseudo_pure,
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
        &WATER.eos.reducing,
        &eos_states["reducing"],
        "EOS[0].STATES.reducing",
    );
    w.state_point(
        &WATER.eos.sat_min_liquid,
        &eos_states["sat_min_liquid"],
        "EOS[0].STATES.sat_min_liquid",
    );
    w.state_point(
        &WATER.eos.sat_min_vapor,
        &eos_states["sat_min_vapor"],
        "EOS[0].STATES.sat_min_vapor",
    );
    w.state_point(
        &WATER.eos.hs_anchor,
        &eos_states["hs_anchor"],
        "EOS[0].STATES.hs_anchor",
    );

    // alpha0 terms
    let alpha0 = eos["alpha0"].as_array().unwrap();
    assert_eq!(WATER.eos.alpha0.len(), alpha0.len(), "alpha0 term count");
    for (i, (term, json)) in WATER.eos.alpha0.iter().zip(alpha0).enumerate() {
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
        }
    }

    // alphar terms
    let alphar = eos["alphar"].as_array().unwrap();
    assert_eq!(WATER.eos.alphar.len(), alphar.len(), "alphar term count");
    for (i, (term, json)) in WATER.eos.alphar.iter().zip(alphar).enumerate() {
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
        }
    }

    // SUPERANCILLARY
    let sa_json = &eos["SUPERANCILLARY"];
    let sa = WATER
        .eos
        .superancillary
        .as_ref()
        .expect("Water has a superancillary");
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
    w.sat_ancillary(&WATER.ancillaries.p_s, &anc["pS"], "ANCILLARIES.pS");
    w.sat_ancillary(&WATER.ancillaries.rho_l, &anc["rhoL"], "ANCILLARIES.rhoL");
    w.sat_ancillary(&WATER.ancillaries.rho_v, &anc["rhoV"], "ANCILLARIES.rhoV");

    // STATES
    let states = &doc["STATES"];
    w.keys(
        states,
        "STATES",
        &["critical", "triple_liquid", "triple_vapor"],
        &[],
    );
    w.state_point(
        &WATER.states.critical,
        &states["critical"],
        "STATES.critical",
    );
    w.state_point(
        &WATER.states.triple_liquid,
        &states["triple_liquid"],
        "STATES.triple_liquid",
    );
    w.state_point(
        &WATER.states.triple_vapor,
        &states["triple_vapor"],
        "STATES.triple_vapor",
    );

    assert!(
        w.mismatches.is_empty(),
        "{} fidelity mismatches:\n{}",
        w.mismatches.len(),
        w.mismatches.join("\n")
    );
}
