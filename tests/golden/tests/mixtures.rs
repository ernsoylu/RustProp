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
