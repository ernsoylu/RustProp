//! IAPWS published check tables (PLAN.md 2.2/2.3 verify) — reference values
//! transcribed from the upstream verification driver IF97.cpp @ 7aaced02,
//! which reproduces them from IAPWS R7-97(2012), SR5-05(2016), SR3-03(2014),
//! SR2-01(2014), SR4-04(2014), R12-08, R15-11, and R1-76(2014).
//!
//! Driver values are in IAPWS units (MPa, kJ); tests convert to SI at the
//! boundary. Region-3 rows of the forward table are checked through the
//! Newton-refined density (`rhomass_iterated`) exactly as the upstream driver
//! does with `REGION3_ITERATE` — the production (CoolProp-config) direct path
//! is covered by the SR5-05 Table 5/13 check at 1e-12 and by the golden
//! fixtures.

use rustprop_if97 as if97;

const MPA: f64 = 1e6;
const KJ: f64 = 1e3;

fn assert_rel(actual: f64, expected: f64, rtol: f64, what: &str) {
    let rel = ((actual - expected) / expected).abs();
    assert!(
        rel <= rtol,
        "{what}: actual {actual:e}, expected {expected:e}, rel {rel:e} > {rtol:e}"
    );
}

/// R7-97 Tables 5, 15, 33, 42 — regions 1, 2, 3, 5, three points each.
#[allow(clippy::type_complexity)]
#[rustfmt::skip]
const FORWARD: [(f64, f64, f64, f64, f64, f64, f64, f64); 12] = [
    // T [K], p [MPa], v [m^3/kg], h [kJ/kg], u [kJ/kg], s [kJ/kg-K], cp [kJ/kg-K], w [m/s]
    (300.0, 3.0, 0.00100215168, 115.331273, 112.324818, 0.392294792, 4.17301218, 1507.73921),
    (300.0, 80.0, 0.000971180894, 184.142828, 106.448356, 0.368563852, 4.01008987, 1634.69054),
    (500.0, 3.0, 0.001202418, 975.542239, 971.934985, 2.58041912, 4.65580682, 1240.71337),
    (300.0, 0.0035, 39.4913866, 2549.91145, 2411.6916, 8.52238967, 1.91300162, 427.920172),
    (700.0, 0.0035, 92.3015898, 3335.68375, 3012.62819, 10.1749996, 2.08141274, 644.289068),
    (700.0, 30.0, 0.00542946619, 2631.49474, 2468.61076, 5.17540298, 10.3505092, 480.386523),
    (650.0, 25.5837018, 0.002, 1863.43019, 1812.26279, 4.05427273, 13.8935717, 502.005554),
    (650.0, 22.2930643, 0.005, 2375.12401, 2263.65868, 4.85438792, 44.6579342, 383.444594),
    (750.0, 78.3095639, 0.002, 2258.68845, 2102.06932, 4.46971906, 6.34165359, 760.696041),
    (1500.0, 0.5, 1.3845509, 5219.76855, 4527.4931, 9.65408875, 2.61609445, 917.06869),
    (1500.0, 30.0, 0.0230761299, 5167.23514, 4474.95124, 7.72970133, 2.72724317, 928.548002),
    (2000.0, 30.0, 0.0311385219, 6571.22604, 5637.07038, 8.53640523, 2.88569882, 1067.36948),
];

#[test]
fn forward_tables_5_15_42_gibbs_regions() {
    for (i, &(t, p_mpa, v, h, u, s, cp, w)) in FORWARD.iter().enumerate() {
        if (6..=8).contains(&i) {
            continue; // region 3 rows handled in the iterated test below
        }
        let p = p_mpa * MPA;
        assert_rel(
            1.0 / if97::rhomass_tp(t, p).unwrap(),
            v,
            1e-8,
            &format!("v row {i}"),
        );
        assert_rel(
            if97::hmass_tp(t, p).unwrap(),
            h * KJ,
            1e-8,
            &format!("h row {i}"),
        );
        assert_rel(
            if97::umass_tp(t, p).unwrap(),
            u * KJ,
            1e-8,
            &format!("u row {i}"),
        );
        assert_rel(
            if97::smass_tp(t, p).unwrap(),
            s * KJ,
            1e-8,
            &format!("s row {i}"),
        );
        assert_rel(
            if97::cpmass_tp(t, p).unwrap(),
            cp * KJ,
            1e-8,
            &format!("cp row {i}"),
        );
        assert_rel(
            if97::speed_sound_tp(t, p).unwrap(),
            w,
            1e-8,
            &format!("w row {i}"),
        );
    }
}

#[test]
fn forward_table_33_region3_via_newton_refined_density() {
    for (i, &(t, p_mpa, v, h, u, s, cp, w)) in FORWARD.iter().enumerate() {
        if !(6..=8).contains(&i) {
            continue;
        }
        let p = p_mpa * MPA;
        let region = if97::region3::region_determination(t, p).unwrap();
        let rho0 = 1.0 / if97::region3::v_tp(region, t, p).unwrap();
        let rho = if97::region3::rhomass_iterated(t, p, rho0).unwrap();
        assert_rel(1.0 / rho, v, 5e-8, &format!("v row {i}"));
        assert_rel(
            if97::region3::hmass(t, rho),
            h * KJ,
            5e-8,
            &format!("h row {i}"),
        );
        assert_rel(
            if97::region3::umass(t, rho),
            u * KJ,
            5e-8,
            &format!("u row {i}"),
        );
        assert_rel(
            if97::region3::smass(t, rho),
            s * KJ,
            5e-8,
            &format!("s row {i}"),
        );
        assert_rel(
            if97::region3::cpmass(t, rho),
            cp * KJ,
            2e-7,
            &format!("cp row {i}"),
        );
        assert_rel(
            if97::region3::speed_sound(t, rho),
            w,
            5e-8,
            &format!("w row {i}"),
        );
    }
}

/// SR5-05(2016) Tables 5 & 13: direct v(T,p) per subregion, and subregion
/// determination — upstream's own check compares at 1e-12 absolute.
#[test]
fn sr5_table5_direct_v_and_subregion_determination() {
    for row in if97::tables::TABLE5 {
        let v = if97::region3::v_tp(row.region, row.t, row.p).unwrap();
        assert!(
            (v - row.v).abs() <= 1e-12,
            "v({}, {}, {}): {v:e} vs {:e}",
            row.region as char,
            row.t,
            row.p,
            row.v
        );
        let det = if97::region3::region_determination(row.t, row.p).unwrap();
        assert_eq!(det, row.region, "subregion at T={} p={}", row.t, row.p);
    }
}

/// SR5-05(2016) Tables 3 & 11: dividing-line temperatures at 1e-7 absolute.
#[test]
fn sr5_table3_dividing_lines() {
    for row in if97::tables::TABLE3 {
        let line = if97::region3::Line::parse(row.line).unwrap();
        let t = if97::region3::dividing_line(line, row.p);
        assert!(
            (t - row.t).abs() <= 1e-7,
            "line {} at p={}: {t} vs {}",
            row.line,
            row.p,
            row.t
        );
    }
}

/// R7-97 Tables 35 & 36: saturation pressure and temperature.
#[test]
fn region4_tables_35_36() {
    let ts = [300.0, 500.0, 600.0];
    let pact = [0.353658941E-2, 0.263889776E+1, 0.123443146E+2];
    for (t, p_mpa) in ts.iter().zip(pact) {
        assert_rel(
            if97::psat97(*t).unwrap(),
            p_mpa * MPA,
            5e-9,
            &format!("psat({t})"),
        );
    }
    let ps = [0.1, 1.0, 10.0];
    let tact = [0.372755919E+3, 0.453035632E+3, 0.584149488E+3];
    for (p_mpa, t) in ps.iter().zip(tact) {
        assert_rel(
            if97::tsat97(p_mpa * MPA).unwrap(),
            t,
            5e-9,
            &format!("Tsat({p_mpa} MPa)"),
        );
    }
}

/// R7-97 Tables 7 & 24 and SR3-03(2014) Tables 5 & 12: T(p,h) and T(p,s),
/// plus the backward region integer.
#[test]
fn backward_t_ph_and_t_ps() {
    let expected_region = [1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3];
    #[rustfmt::skip]
    let ph: [(f64, f64, f64); 18] = [
        (3.0, 500.0, 391.798509), (80.0, 500.0, 378.108626), (80.0, 1500.0, 611.041229),
        (0.001, 3000.0, 534.433241), (3.0, 3000.0, 575.373370), (3.0, 4000.0, 1010.77577),
        (5.0, 3500.0, 801.299102), (5.0, 4000.0, 1015.31583), (25.0, 3500.0, 875.279054),
        (40.0, 2700.0, 743.056411), (60.0, 2700.0, 791.137067), (60.0, 3200.0, 882.756860),
        (20.0, 1700.0, 629.3083892), (50.0, 2000.0, 690.5718338), (100.0, 2100.0, 733.6163014),
        (20.0, 2500.0, 641.8418053), (50.0, 2400.0, 735.1848618), (100.0, 2700.0, 842.0460876),
    ];
    for (i, &(p_mpa, h_kj, t_expected)) in ph.iter().enumerate() {
        let (p, h) = (p_mpa * MPA, h_kj * KJ);
        assert_rel(
            if97::t_phmass(p, h).unwrap(),
            t_expected,
            5e-9,
            &format!("T(p,h) row {i}"),
        );
        assert_eq!(
            if97::region_ph(p, h).unwrap(),
            expected_region[i],
            "region_ph row {i}"
        );
    }
    #[rustfmt::skip]
    let ps: [(f64, f64, f64); 18] = [
        (3.0, 0.5, 307.842258), (80.0, 0.5, 309.979785), (80.0, 3.0, 565.899909),
        (0.1, 7.5, 399.517097), (0.1, 8.0, 514.127081), (2.5, 8.0, 1039.84917),
        (8.0, 6.0, 600.48404), (8.0, 7.5, 1064.95556), (90.0, 6.0, 1038.01126),
        (20.0, 5.75, 697.992849), (80.0, 5.25, 854.011484), (80.0, 5.75, 949.017998),
        (20.0, 3.8, 628.2959869), (50.0, 3.6, 629.7158726), (100.0, 4.0, 705.6880237),
        (20.0, 5.0, 640.1176443), (50.0, 4.5, 716.3687517), (100.0, 5.0, 847.4332825),
    ];
    for (i, &(p_mpa, s_kj, t_expected)) in ps.iter().enumerate() {
        let (p, s) = (p_mpa * MPA, s_kj * KJ);
        assert_rel(
            if97::t_psmass(p, s).unwrap(),
            t_expected,
            5e-9,
            &format!("T(p,s) row {i}"),
        );
    }
}

/// SR2-01(2014) Tables 3 & 9 and SR4-04(2014) Tables 5 & 29: p(h,s) and
/// Tsat(h,s).
#[test]
fn backward_p_hs_and_tsat_hs() {
    #[rustfmt::skip]
    let hs_p: [(f64, f64, f64); 18] = [
        (0.001, 0.0, 9.800980614e-4), (90.0, 0.0, 91.92954727), (1500.0, 3.4, 58.68294423),
        (2800.0, 6.5, 1.371012767), (2800.0, 9.5, 1.879743844e-3), (4100.0, 9.5, 0.1024788997),
        (2800.0, 6.0, 4.793911442), (3600.0, 6.0, 83.95519209), (3600.0, 7.0, 7.527161441),
        (2800.0, 5.1, 94.3920206), (2800.0, 5.8, 8.414574124), (3400.0, 5.8, 83.76903879),
        (1700.0, 3.8, 25.55703246), (2000.0, 4.2, 45.40873468), (2100.0, 4.3, 60.7812334),
        (2600.0, 5.1, 34.34999263), (2400.0, 4.7, 63.63924887), (2700.0, 5.0, 88.39043281),
    ];
    for (i, &(h_kj, s_kj, p_mpa)) in hs_p.iter().enumerate() {
        let p = if97::p_hsmass(h_kj * KJ, s_kj * KJ).unwrap();
        assert_rel(p, p_mpa * MPA, 1e-9, &format!("p(h,s) row {i}"));
    }
    let hs_t: [(f64, f64, f64); 3] = [
        (1800.0, 5.3, 346.8475498),
        (2400.0, 6.0, 425.1373305),
        (2500.0, 5.5, 522.5579013),
    ];
    for (i, &(h_kj, s_kj, t_expected)) in hs_t.iter().enumerate() {
        let t = if97::t_hsmass(h_kj * KJ, s_kj * KJ).unwrap();
        assert_rel(t, t_expected, 1e-9, &format!("Tsat(h,s) row {i}"));
    }
}

/// R12-08 Table 4: viscosity from (T, rho). The last six points sit near the
/// critical point where the industrial formulation (mu2 = 1) deviates by
/// design — upstream documents errors up to 8.417e-2 there.
#[test]
fn viscosity_r12_table4() {
    #[rustfmt::skip]
    let rows: [(f64, f64, f64); 17] = [
        (298.15, 998.0, 889.7351), (298.15, 1200.0, 1437.649467), (373.15, 1000.0, 307.883622),
        (433.15, 1.0, 14.538324), (433.15, 1000.0, 217.685358), (873.15, 1.0, 32.619287),
        (873.15, 100.0, 35.802262), (873.15, 600.0, 77.430195), (1173.15, 1.0, 44.217245),
        (1173.15, 100.0, 47.640433), (1173.15, 400.0, 64.154608),
        (647.35, 122.0, 25.520677), (647.35, 222.0, 31.337589), (647.35, 272.0, 36.228143),
        (647.35, 322.0, 42.961579), (647.35, 372.0, 45.688204), (647.35, 422.0, 49.436256),
    ];
    for (i, &(t, rho, mu_upas)) in rows.iter().enumerate() {
        let mu = if97::visc_trho(t, rho) * 1.0e6;
        let rtol = if i < 11 { 1e-7 } else { 0.09 };
        assert_rel(mu, mu_upas, rtol, &format!("visc row {i}"));
    }
}

/// R15-11 Tables 7, 8, 9: thermal conductivity via tcond_tp. Near-critical
/// rows use the simplified industrial enhancement and, in our CoolProp
/// configuration, the direct region-3 density — both looser by design.
#[test]
fn conductivity_r15_tables789() {
    #[rustfmt::skip]
    let rows: [(f64, f64, f64); 13] = [
        // T [K], p [MPa], lambda [mW/m/K]
        (620.0, 20.0, 0.481485195 * 1000.0), (620.0, 50.0, 0.54503894 * 1000.0),
        (650.0, 0.3, 0.0522311024 * 1000.0), (800.0, 50.0, 0.177709914 * 1000.0),
        (647.35, 21.98406271345, 0.36687941 * 1000.0), (647.35, 22.1321600249828, 1.24182415 * 1000.0),
        (647.35, 0.297422657, 51.9298924), (647.35, 19.45771946, 130.922885),
        (647.35, 21.98406271, 367.787459), (647.35, 22.11526557, 757.959776),
        (647.35, 22.13216002, 1443.75556), (647.35, 22.15298122, 650.319402),
        (647.35, 22.33268694, 448.883487),
    ];
    for (i, &(t, p_mpa, k_mw)) in rows.iter().enumerate() {
        let k = if97::tcond_tp(t, p_mpa * MPA).unwrap() * 1000.0;
        let rtol = if i < 4 { 1e-7 } else { 0.15 };
        assert_rel(k, k_mw, rtol, &format!("tcond row {i}"));
    }
}

/// R1-76(2014) Table 1: surface tension (table values are 2-decimal mN/m).
#[test]
fn surface_tension_r1_table1() {
    #[rustfmt::skip]
    let t_c: [f64; 75] = [
        0.01, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0,
        75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0,
        140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0,
        205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0,
        270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0,
        335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0,
    ];
    #[rustfmt::skip]
    let sigma_mn: [f64; 75] = [
        75.64, 74.94, 74.23, 73.49, 72.74, 71.98, 71.19, 70.41, 69.59, 68.78, 67.93, 67.09, 66.24,
        65.36, 64.47, 63.57, 62.68, 61.76, 60.82, 59.88, 58.92, 57.95, 56.97, 55.98, 54.97, 53.96,
        52.94, 51.9, 50.86, 49.81, 48.75, 47.67, 46.58, 45.49, 44.4, 43.3, 42.19, 41.07, 39.95,
        38.82, 37.68, 36.54, 35.4, 34.24, 33.09, 31.92, 30.76, 29.58, 28.4, 27.22, 26.05, 24.86,
        23.66, 22.46, 21.29, 20.14, 18.93, 17.76, 16.6, 15.45, 14.3, 13.18, 12.04, 10.92, 9.81,
        8.73, 7.66, 6.61, 5.59, 4.6, 3.64, 2.74, 1.89, 1.12, 0.45,
    ];
    for (i, (t_c, sig)) in t_c.iter().zip(sigma_mn).enumerate() {
        let sigma = if97::sigma97(t_c + 273.15).unwrap() * 1000.0;
        if *t_c > 360.0 {
            // sigma -> 0 at the critical point: the table's last rows differ
            // from the normative correlation by up to 0.07 mN/m absolute.
            assert!(
                (sigma - sig).abs() <= 0.07,
                "sigma row {i} (T={t_c} C): {sigma} vs {sig}"
            );
        } else {
            let rtol = if *t_c < 260.0 { 2e-3 } else { 2e-2 };
            assert_rel(sigma, sig, rtol, &format!("sigma row {i} (T={t_c} C)"));
        }
    }
}

/// Trivial accessors and a few sanity anchors.
#[test]
fn trivials() {
    assert_eq!(if97::get_tcrit(), 647.096);
    assert_eq!(if97::get_pcrit(), 22.064e6);
    assert_eq!(if97::get_rhocrit(), 322.0);
    assert_eq!(if97::get_mw(), 0.018015268);
    let acentric = if97::get_acentric().unwrap();
    assert!((0.34..0.35).contains(&acentric), "acentric {acentric}");
}
