//! Mixture goldens (PLAN.md Phase 10, slice 10b): the GERG-2008 reducing
//! function against the wheel's `T_reducing()`/`rhomolar_reducing()` for
//! GERG- and Lemmon-converted binary pairs across compositions.

use rustprop_core::fluid::FluidData;
use rustprop_heos::mixture::Gerg2008Reducing;

fn fluid(name: &str) -> &'static FluidData {
    let registry: std::collections::HashMap<&str, &'static FluidData> =
        rustprop_data::fluids::all().into_iter().collect();
    registry[name]
}

#[test]
fn reducing_function_matches_upstream() {
    // (fluid1, fluid2, x1, Tr_oracle, rhor_oracle) — from the wheel.
    let cases = [
        (
            "Methane",
            "Ethane",
            0.5,
            250.5718599087635,
            8205.78588373207,
        ),
        (
            "Methane",
            "Ethane",
            0.25,
            278.43022073906417,
            7479.200563404327,
        ),
        (
            "Nitrogen",
            "CarbonDioxide",
            0.3,
            251.7309025964105,
            10619.912488821874,
        ),
        ("R32", "R125", 0.6973, 353.7083335545, 6773.742371291401),
        (
            "Methane",
            "n-Propane",
            0.8,
            230.67616922000263,
            8429.138696443493,
        ),
    ];
    for (f1, f2, x1, tr_exp, rhor_exp) in cases {
        let comps = [fluid(f1), fluid(f2)];
        let red = Gerg2008Reducing::new(&comps, rustprop_data::mixtures::MIX_BINARY_PAIRS)
            .expect("pair present");
        let x = [x1, 1.0 - x1];
        let tr = red.tr(&x);
        let rhor = red.rhormolar(&x);
        assert!(
            ((tr - tr_exp) / tr_exp).abs() < 1e-12,
            "{f1}&{f2} Tr: {tr} vs {tr_exp}"
        );
        assert!(
            ((rhor - rhor_exp) / rhor_exp).abs() < 1e-12,
            "{f1}&{f2} rhor: {rhor} vs {rhor_exp}"
        );
    }
}

/// The composition derivatives satisfy their finite-difference identities
/// (the wheel exposes no direct accessors for these; the flashes that
/// consume them get golden-verified end to end in slices 10d/10e).
#[test]
fn reducing_derivatives_consistent() {
    use rustprop_heos::mixture::XnFlag;
    let comps = [fluid("Methane"), fluid("Ethane")];
    let red = Gerg2008Reducing::new(&comps, rustprop_data::mixtures::MIX_BINARY_PAIRS).unwrap();
    let x = [0.4, 0.6];
    let h = 1e-7;
    for i in 0..2 {
        let mut xp = x;
        xp[i] += h;
        let mut xm = x;
        xm[i] -= h;
        let fd_t = (red.tr(&xp) - red.tr(&xm)) / (2.0 * h);
        let an_t = red.dtrdxi__constxj(&x, i, XnFlag::Independent);
        assert!(
            ((fd_t - an_t) / an_t).abs() < 1e-6,
            "dTr/dx{i}: fd {fd_t} vs analytic {an_t}"
        );
        let fd_v = (1.0 / red.rhormolar(&xp) - 1.0 / red.rhormolar(&xm)) / (2.0 * h);
        let an_v = red.dvrmolardxi__constxj(&x, i, XnFlag::Independent);
        assert!(
            ((fd_v - an_v) / an_v).abs() < 1e-6,
            "dvr/dx{i}: fd {fd_v} vs analytic {an_v}"
        );
        for j in 0..2 {
            let fd2 = (red.dtrdxi__constxj(&xp, j, XnFlag::Independent)
                - red.dtrdxi__constxj(&xm, j, XnFlag::Independent))
                / (2.0 * h);
            let an2 = red.d2trdxidxj(&x, j, i, XnFlag::Independent);
            assert!(
                (fd2 - an2).abs() / an2.abs().max(1.0) < 1e-5,
                "d2Tr/dx{j}dx{i}: fd {fd2} vs analytic {an2}"
            );
        }
    }
}

#[test]
fn mixture_helmholtz_matches_oracle() {
    // Slice 10c: corresponding-states + excess alphar and Table B5 alpha0
    // against the wheel's low-level accessors. Direct evaluations: 1e-12.
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/mixture_helmholtz.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 864);

    let mut models: std::collections::HashMap<String, rustprop_heos::mixture::MixtureModel> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &recs {
        let (f1, f2) = rec.fluid.split_once('&').expect("pair fluid string");
        let model = models.entry(rec.fluid.clone()).or_insert_with(|| {
            rustprop_heos::mixture::MixtureModel::new(
                &[fluid(f1), fluid(f2)],
                rustprop_data::mixtures::MIX_BINARY_PAIRS,
                rustprop_data::mixtures::MIX_DEPARTURE_FNS,
            )
            .expect("model builds")
        });
        let x = [rec.val3, 1.0 - rec.val3];
        let tr = model.reducing.tr(&x);
        let rhor = model.reducing.rhormolar(&x);
        let tau = tr / rec.val1;
        let delta = rec.val2 / rhor;
        let actual = if rec.out.starts_with("alphar") || rec.out.contains("alphar_") {
            let ar = model.alphar_all(&x, tau, delta);
            match rec.out.as_str() {
                "alphar" => ar.d00,
                "dalphar_dTau" => ar.d01,
                "dalphar_dDelta" => ar.d10,
                "d2alphar_dTau2" => ar.d02,
                "d2alphar_dDelta2" => ar.d20,
                "d2alphar_dDelta_dTau" => ar.d11,
                other => panic!("unknown out {other}"),
            }
        } else {
            let a0 = model.alpha0_all(&x, tau, delta, tr, rhor);
            match rec.out.as_str() {
                "alpha0" => a0.d00,
                "dalpha0_dTau" => a0.d01,
                "dalpha0_dDelta" => a0.d10,
                "d2alpha0_dTau2" => a0.d02,
                "d2alpha0_dDelta2" => a0.d20,
                "d2alpha0_dDelta_dTau" => a0.d11,
                other => panic!("unknown out {other}"),
            }
        };
        if let Err(e) = rustprop_golden_tests::check(rec, actual, 1e-12) {
            failures.push(format!("{} x1={}: {e}", rec.fluid, rec.val3));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} mixture Helmholtz records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn mixture_pt_flash_matches_oracle() {
    // Slice 10d: PT single-phase flash (SRK seed, lowest-Gibbs root pick) and
    // homogeneous mixture properties. Solver tier: 1e-9.
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/mixture_pt.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 696);

    let mut models: std::collections::HashMap<String, rustprop_heos::mixture::MixtureModel> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &recs {
        let (f1, f2) = rec.fluid.split_once('&').expect("pair fluid string");
        let model = models.entry(rec.fluid.clone()).or_insert_with(|| {
            rustprop_heos::mixture::MixtureModel::new(
                &[fluid(f1), fluid(f2)],
                rustprop_data::mixtures::MIX_BINARY_PAIRS,
                rustprop_data::mixtures::MIX_DEPARTURE_FNS,
            )
            .expect("model builds")
        });
        let x = [rec.val3, 1.0 - rec.val3];
        let (t, p) = (rec.val1, rec.val2);
        let state = match model.pt_flash(&x, t, p) {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!(
                    "{} x1={}: {} flash failed: {e:?}",
                    rec.fluid,
                    rec.val3,
                    rec.id()
                ));
                continue;
            }
        };
        let rho = state.rhomolar;
        let actual = match rec.out.as_str() {
            "Dmolar" => rho,
            "Hmolar" => model.hmolar(&x, t, rho),
            "Smolar" => model.smolar(&x, t, rho),
            "Umolar" => model.umolar(&x, t, rho),
            "Cpmolar" => model.cpmolar(&x, t, rho),
            "Cvmolar" => model.cvmolar(&x, t, rho),
            "speed_of_sound" => model.speed_sound(&x, t, rho),
            "Gmolar" => model.gibbsmolar_nocache(&x, t, rho),
            other => panic!("unknown out {other}"),
        };
        if let Err(e) = rustprop_golden_tests::check(rec, actual, 1e-9) {
            failures.push(format!("{} x1={}: {e}", rec.fluid, rec.val3));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} mixture PT records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn mixture_vle_matches_oracle() {
    // Slice 10e: blind QT/PQ flashes. Solver tier: 1e-8 (NR converges the
    // ln-fugacity residuals to 1e-7 rms on both sides from the same seed).
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/mixture_vle.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 720);

    let mut models: std::collections::HashMap<String, rustprop_heos::mixture::MixtureModel> =
        std::collections::HashMap::new();
    let mut failures = Vec::new();
    for rec in &recs {
        let (f1, f2) = rec.fluid.split_once('&').expect("pair fluid string");
        let model = models.entry(rec.fluid.clone()).or_insert_with(|| {
            rustprop_heos::mixture::MixtureModel::new(
                &[fluid(f1), fluid(f2)],
                rustprop_data::mixtures::MIX_BINARY_PAIRS,
                rustprop_data::mixtures::MIX_DEPARTURE_FNS,
            )
            .expect("model builds")
        });
        let z = [rec.val3, 1.0 - rec.val3];
        let state = match rec.name1.as_str() {
            "T" => model.qt_flash(rec.val2, rec.val1, &z),
            "P" => model.pq_flash(rec.val1, rec.val2, &z),
            other => panic!("unknown imposed variable {other}"),
        };
        let state = match state {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!(
                    "{} x1={}: {} flash failed: {e:?}",
                    rec.fluid,
                    rec.val3,
                    rec.id()
                ));
                continue;
            }
        };
        let actual = match rec.out.as_str() {
            "T" => state.t,
            "P" => state.p,
            "Dmolar" => state.rhomolar,
            "Hmolar" => state.hmolar(),
            "Smolar" => state.smolar(),
            other => panic!("unknown out {other}"),
        };
        if let Err(e) = rustprop_golden_tests::check(rec, actual, 1e-8) {
            failures.push(format!("{} x1={}: {e}", rec.fluid, rec.val3));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} mixture VLE records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

/// Same identities under the XN_DEPENDENT convention (x1 = 1 - x0): the
/// newton_raphson_saturation Jacobian consumes exactly these — a missing
/// Dependent branch in d2Yrdxidxj once sent the 10e NR to a wrong region.
#[test]
fn reducing_derivatives_consistent_dependent() {
    use rustprop_heos::mixture::XnFlag;
    let comps = [fluid("Methane"), fluid("Ethane")];
    let red = Gerg2008Reducing::new(&comps, rustprop_data::mixtures::MIX_BINARY_PAIRS).unwrap();
    let flag = XnFlag::Dependent;
    let x = [0.4_f64, 0.6];
    let h = 1e-7;
    let xp = [x[0] + h, x[1] - h];
    let xm = [x[0] - h, x[1] + h];

    let fd_t = (red.tr(&xp) - red.tr(&xm)) / (2.0 * h);
    let an_t = red.dtrdxi__constxj(&x, 0, flag);
    assert!(
        ((fd_t - an_t) / an_t).abs() < 1e-6,
        "dTr/dx0 dep: fd {fd_t} vs analytic {an_t}"
    );
    let fd_r = (red.rhormolar(&xp) - red.rhormolar(&xm)) / (2.0 * h);
    let an_r = red.drhormolardxi__constxj(&x, 0, flag);
    assert!(
        ((fd_r - an_r) / an_r).abs() < 1e-6,
        "drhor/dx0 dep: fd {fd_r} vs analytic {an_r}"
    );
    let fd2_t = (red.dtrdxi__constxj(&xp, 0, flag) - red.dtrdxi__constxj(&xm, 0, flag)) / (2.0 * h);
    let an2_t = red.d2trdxidxj(&x, 0, 0, flag);
    assert!(
        (fd2_t - an2_t).abs() / an2_t.abs().max(1.0) < 1e-5,
        "d2Tr/dx0dx0 dep: fd {fd2_t} vs analytic {an2_t}"
    );
    let fd2_r = (red.drhormolardxi__constxj(&xp, 0, flag)
        - red.drhormolardxi__constxj(&xm, 0, flag))
        / (2.0 * h);
    let an2_r = red.d2rhormolardxidxj(&x, 0, 0, flag);
    assert!(
        (fd2_r - an2_r).abs() / an2_r.abs().max(1.0) < 1e-5,
        "d2rhor/dx0dx0 dep: fd {fd2_r} vs analytic {an2_r}"
    );
    let fd_nd =
        (red.ndtrdni__constnj(&xp, 0, flag) - red.ndtrdni__constnj(&xm, 0, flag)) / (2.0 * h);
    let an_nd = red.d_ndtrdni_dxj__constxi(&x, 0, 0, flag);
    assert!(
        (fd_nd - an_nd).abs() / an_nd.abs().max(1.0) < 1e-5,
        "d_ndTrdni_dxj dep: fd {fd_nd} vs analytic {an_nd}"
    );
    let fd_ndr = (red.ndrhorbardni__constnj(&xp, 0, flag)
        - red.ndrhorbardni__constnj(&xm, 0, flag))
        / (2.0 * h);
    let an_ndr = red.d_ndrhorbardni_dxj__constxi(&x, 0, 0, flag);
    assert!(
        (fd_ndr - an_ndr).abs() / an_ndr.abs().max(1.0) < 1e-5,
        "d_ndrhorbardni_dxj dep: fd {fd_ndr} vs analytic {an_ndr}"
    );
}

#[test]
fn mixture_propssi_matches_oracle() {
    // PropsSI mixture-string routing end to end: trivials at 1e-12 (direct
    // arithmetic), states at 1e-8 (solver tier).
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/mixture_propssi.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 284);
    let mut failures = Vec::new();
    for rec in &recs {
        let fluid = format!("HEOS::{}", rec.fluid);
        match rustprop::props_si(&rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &fluid) {
            Ok(actual) => {
                let rtol = if rec.name1.is_empty() { 1e-12 } else { 1e-8 };
                if let Err(e) = rustprop_golden_tests::check(rec, actual, rtol) {
                    failures.push(e);
                }
            }
            Err(e) => failures.push(format!("{}: error {e:?}", rec.id())),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} mixture PropsSI records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn mixture_propssi_error_conditions() {
    use rustprop::props_si;
    use rustprop_core::Error;
    let err = |r: rustprop_core::Result<f64>| r.unwrap_err();
    let mix = "HEOS::Methane[0.6]&Ethane[0.4]";

    // Fraction parsing errors (verbatim upstream messages)
    match err(props_si(
        "T",
        "P",
        1e5,
        "Q",
        0.0,
        "HEOS::Methane[0.6]&Ethane",
    )) {
        Error::Value(m) => assert_eq!(m, "Fluid entry [Ethane] must end with ']' character"),
        other => panic!("wrong variant {other:?}"),
    }
    match err(props_si(
        "T",
        "P",
        1e5,
        "Q",
        0.0,
        "HEOS::Methane[1.2]&Ethane[-0.2]",
    )) {
        Error::Value(m) => assert_eq!(
            m,
            "fraction [1.2] was not converted to a value between 0 and 1 inclusive"
        ),
        other => panic!("wrong variant {other:?}"),
    }
    match err(props_si(
        "T",
        "P",
        1e5,
        "Q",
        0.0,
        "HEOS::Methane[abc]&Ethane[0.4]",
    )) {
        Error::Value(m) => assert_eq!(m, "fraction [abc] was not converted fully"),
        other => panic!("wrong variant {other:?}"),
    }
    // No fractions with '&': the factory-size mismatch
    match err(props_si("T", "P", 1e5, "Q", 0.0, "HEOS::Methane&Ethane")) {
        Error::Value(m) => assert!(m.contains(
            "size of mole fraction vector [1] does not equal that of component vector [2]"
        )),
        other => panic!("wrong variant {other:?}"),
    }
    // Upstream's own mixture dead ends (verbatim)
    match err(props_si("T", "Dmolar", 100.0, "P", 1e5, mix)) {
        Error::NotImplemented(m) => assert_eq!(m, "DP_flash not ready for mixtures"),
        other => panic!("wrong variant {other:?}"),
    }
    match err(props_si("T", "Dmolar", 100.0, "Q", 0.5, mix)) {
        Error::NotImplemented(m) => assert_eq!(m, "DQ_flash not ready for mixtures"),
        other => panic!("wrong variant {other:?}"),
    }
    // HQ/QS: no update-pair row exists (string-API dead end, same as pure)
    match err(props_si("T", "Hmolar", 1e4, "Q", 0.5, mix)) {
        Error::Value(m) => assert!(m.contains("Input pair variable is invalid")),
        other => panic!("wrong variant {other:?}"),
    }
    // Q out of range
    match err(props_si("P", "T", 200.0, "Q", 1.5, mix)) {
        Error::OutOfRange(m) => {
            assert_eq!(m, "Input vapor quality [Q] must be between 0 and 1")
        }
        other => panic!("wrong variant {other:?}"),
    }
    // Critical-family trivials fail with PropsSI's wrapper message
    match err(props_si("Tcrit", "", 0.0, "", 0.0, mix)) {
        Error::Value(m) => assert_eq!(m, "No outputs were able to be calculated"),
        other => panic!("wrong variant {other:?}"),
    }
    // Unmatched binary pair carries upstream's message
    match err(props_si(
        "T",
        "P",
        1e5,
        "Q",
        0.0,
        "HEOS::Methanol[0.5]&Novec649[0.5]",
    )) {
        Error::Value(m) => assert!(m.contains("Could not match the binary pair"), "got: {m}"),
        other => panic!("wrong variant {other:?}"),
    }
    // Pure-fluid collapse: single component after parsing behaves as pure
    let pure = props_si("Dmolar", "T", 300.0, "P", 1e5, "HEOS::Water[1.0]").unwrap();
    let direct = props_si("Dmolar", "T", 300.0, "P", 1e5, "Water").unwrap();
    assert_eq!(pure, direct);
    // Deferred sweep pairs error loudly (documented 10f deviation)
    match err(props_si("P", "Dmolar", 100.0, "T", 300.0, mix)) {
        Error::NotImplemented(m) => assert!(m.contains("deferred with slice 10f"), "got: {m}"),
        other => panic!("wrong variant {other:?}"),
    }
}

#[test]
fn mixture_predefined_matches_oracle() {
    // Predefined "<Name>.mix" blends through the string API — binary,
    // ternary (Air, R404A/R407C), and the 10-component Amarillo natural gas.
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/mixture_predefined.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 175);
    let mut failures = Vec::new();
    for rec in &recs {
        let fluid = format!("HEOS::{}", rec.fluid);
        match rustprop::props_si(&rec.out, &rec.name1, rec.val1, &rec.name2, rec.val2, &fluid) {
            Ok(actual) => {
                let rtol = if rec.name1.is_empty() { 1e-12 } else { 1e-8 };
                if let Err(e) = rustprop_golden_tests::check(rec, actual, rtol) {
                    failures.push(e);
                }
            }
            Err(e) => failures.push(format!("{}: error {e:?}", rec.id())),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} predefined-mixture records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}

#[test]
fn mixture_predefined_case_sensitivity() {
    use rustprop::props_si;
    // Exact and uppercase keys resolve; other casings fall through to the
    // pure registry and fail its lookup (same as upstream's factory error
    // path, whose outer message wraps the same missing-key condition).
    assert!(props_si("molemass", "", 0.0, "", 0.0, "HEOS::R410A.mix").is_ok());
    assert!(props_si("molemass", "", 0.0, "", 0.0, "HEOS::R410A.MIX").is_ok());
    assert!(props_si("molemass", "", 0.0, "", 0.0, "HEOS::r410a.mix").is_err());
    assert!(props_si("molemass", "", 0.0, "", 0.0, "HEOS::R410a.mix").is_err());
    // The predefined R410A blend and the pseudo-pure R410A fluid are
    // different models: same name stem, different route.
    let blend = props_si("Dmolar", "T", 300.0, "P", 1e5, "HEOS::R410A.mix").unwrap();
    let pseudo = props_si("Dmolar", "T", 300.0, "P", 1e5, "HEOS::R410A").unwrap();
    assert!(blend != pseudo, "blend {blend} vs pseudo-pure {pseudo}");
}

#[test]
fn mixture_pt_twophase_matches_oracle() {
    // Slice 10f: in-dome PT — stability test + Michelsen split. Both sides
    // converge the same equilibrium to ~1e-9 fugacity residuals; 1e-6 policy
    // for the split-derived outputs: Q's condition number against that
    // residual scales with 1/(composition spread), and the near-azeotropic
    // R32&R125 (spread ~0.03, K ~ 1) observes up to ~5e-7.
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/mixture_pt_twophase.jsonl");
    let recs = rustprop_golden_tests::load_jsonl(&path);
    assert_eq!(recs.len(), 192);
    let mut failures = Vec::new();
    for rec in &recs {
        let (f1, f2) = rec.fluid.split_once('&').expect("pair fluid string");
        let fluid = format!("HEOS::{f1}[{}]&{f2}[{}]", rec.val3, 1.0 - rec.val3);
        match rustprop::props_si(&rec.out, "T", rec.val1, "P", rec.val2, &fluid) {
            Ok(actual) => {
                if let Err(e) = rustprop_golden_tests::check(rec, actual, 1e-6) {
                    failures.push(e);
                }
            }
            Err(e) => failures.push(format!("{}: error {e:?}", rec.id())),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} in-dome PT records failed:\n{}",
        failures.len(),
        recs.len(),
        failures.join("\n")
    );
}
