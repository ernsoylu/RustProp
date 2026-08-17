//! Classic saturation ancillaries (PLAN.md 4.3) — port of
//! `SaturationAncillaryFunction::evaluate` from
//! `src/Backends/Helmholtz/Fluids/Ancillaries.cpp` @ v8.0.0 for the forms the
//! fluid documents' `pS`/`rhoL`/`rhoV` blocks use: type `"rhoLnoexp"` is the
//! non-exponential form, every other tag ("pV", "rhoV", ...) the exponential
//! form. The `rational_polynomial` type (the pseudo-pure hL/hLV/sL/sLV
//! caloric blocks) is `evaluate_rational_poly`.

use rustprop_core::Result;
use rustprop_core::fluid::{RationalPolyAncillary, SaturationAncillary};

/// Upstream `SaturationAncillaryFunction::evaluate`.
///
/// Above the reducing temperature (theta < 0) this returns NaN, mirroring
/// upstream's deliberate #1611 semantics: callers probing supercritical
/// temperatures rely on well-behaved-but-invalid rather than an exception.
pub fn evaluate(anc: &SaturationAncillary, t: f64) -> f64 {
    let theta = 1.0 - t / anc.t_r;
    if theta < 0.0 {
        return f64::NAN;
    }
    // Upstream fills s[i] then left-folds with std::accumulate — same order.
    let mut summer = 0.0;
    for i in 0..anc.n.len() {
        summer += anc.n[i] * theta.powf(anc.t[i]);
    }
    if anc.anc_type == "rhoLnoexp" {
        anc.reducing_value * (1.0 + summer)
    } else {
        let tau_r_value = if anc.using_tau_r { anc.t_r / t } else { 1.0 };
        anc.reducing_value * (tau_r_value * summer).exp()
    }
}

/// Upstream `SurfaceTensionCorrelation::evaluate`:
/// `sigma = sum a_i * (1 - T/Tc)^n_i` [N/m]; `T > Tc` errors.
/// Upstream `SaturationAncillaryFunction::invert` (Ancillaries.cpp:82-113): T
/// for a given output value — Brent on [Tmin - 0.01, Tmax], then
/// `ExtrapolatingSecant` from Tmax with a -0.01 step. The fallback must be the
/// EXTRAPOLATING secant, not plain `Secant`: inverting the pseudo-pure pL/pV
/// curves at a pressure above their value at Tmax (the pcrit..max_sat_p band,
/// where Brent cannot bracket) walks the residual past the reducing
/// temperature, where `evaluate` returns NaN by design, and upstream
/// extrapolates through that instead of failing.
pub fn invert(anc: &SaturationAncillary, value: f64) -> Result<f64> {
    let f = |t: f64| evaluate(anc, t) - value;
    match crate::solvers::brent(f, anc.t_min - 0.01, anc.t_max, f64::EPSILON, 1e-10, 100) {
        Ok(t) => Ok(t),
        Err(_) => crate::solvers::extrapolating_secant(f, anc.t_max, -0.01, 1e-12, 100),
    }
}

/// Upstream `SaturationAncillaryFunction::evaluate` for
/// `TYPE_RATIONAL_POLYNOMIAL` (Ancillaries.cpp:46-48): the pseudo-pure
/// hL/hLV/sL/sLV caloric curves, `poly(A, T) / poly(B, T)`.
pub fn evaluate_rational_poly(anc: &RationalPolyAncillary, t: f64) -> f64 {
    poly_eval(anc.a, t) / poly_eval(anc.b, t)
}

/// Upstream `Polynomial2D::evaluate` (PolyMath.cpp:159-163), which is
/// `Eigen::poly_eval` (Eigen 5.0.1
/// `unsupported/Eigen/src/Polynomials/PolynomialUtils.h:47-61`) on the
/// degree-ascending coefficients: for `x^2 <= 1` plain Horner from the top
/// coefficient; otherwise Eigen's "stabilized" scheme — fold `1/x` upward
/// from `c[0]` and rescale by `x^(n-1)`. Kelvin inputs always take the
/// second branch; both are replicated for exactness.
fn poly_eval(c: &[f64], x: f64) -> f64 {
    if x * x <= 1.0 {
        let mut val = c[c.len() - 1];
        for &ci in c[..c.len() - 1].iter().rev() {
            val = val * x + ci;
        }
        val
    } else {
        let mut val = c[0];
        let inv_x = 1.0 / x;
        for &ci in &c[1..] {
            val = val * inv_x + ci;
        }
        x.powf((c.len() - 1) as f64) * val
    }
}

pub fn surface_tension(
    st: &rustprop_core::fluid::SurfaceTension,
    t: f64,
) -> rustprop_core::Result<f64> {
    if t > st.tc {
        return Err(rustprop_core::Error::Value(
            "Must be saturated state : T <= Tc".into(),
        ));
    }
    let theta = 1.0 - t / st.tc;
    Ok(st.a.iter().zip(st.n).map(|(a, n)| a * theta.powf(*n)).sum())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Expected values computed independently with the oracle venv's Python
    /// (3.12) by replicating Eigen 5.0.1 `poly_eval` on the A/B arrays taken
    /// verbatim from `data/coolprop-json/{Air,R410A}.json` — CPython and
    /// Rust perform the same IEEE-754 double ops (no FMA contraction in
    /// either), so agreement must be bitwise. The wheel exposes no direct
    /// hL/sL evaluation route (`calc_saturation_ancillary` serves only
    /// pL/pV/rhoL/rhoV/sigma), hence the replica oracle.
    #[test]
    fn rational_poly_matches_python_oracle_bitwise() {
        let air_h_l = RationalPolyAncillary {
            a: &[
                138123.1143472345,
                -10990.781133976443,
                353.76801616323996,
                -6.196252980236228,
                0.06433161035294266,
                -0.00039662318978531325,
                1.3454787957202308e-6,
                -1.93892599649897e-9,
            ],
            b: &[1.0, -0.00752748779665181],
            max_abs_error: 29.91964643871188,
            t_min: 59.75,
            t_max: 132.4306,
        };
        let air_s_lv = RationalPolyAncillary {
            a: &[
                -13823.264935008381,
                1102.4746662606901,
                -36.608903508494514,
                0.6620654722822006,
                -0.007061778753904959,
                4.4478254997636e-5,
                -1.533052652255337e-7,
                2.232451590190225e-10,
            ],
            b: &[1.0, -0.005030285100436906],
            max_abs_error: 4.036921827207072,
            t_min: 59.75,
            t_max: 132.4306,
        };
        let r410a_h_lv = RationalPolyAncillary {
            a: &[
                13181208.205563573,
                -353816.44393662986,
                4053.373798329682,
                -25.65820688339148,
                0.0969058419286629,
                -0.000218383941900713,
                2.7193316296232616e-7,
                -1.4435015378161098e-10,
            ],
            b: &[1.0, -0.0026896548966614153],
            max_abs_error: 789.9373183003213,
            t_min: 200.0,
            t_max: 344.394,
        };
        let r410a_s_l = RationalPolyAncillary {
            a: &[
                70603.1839083014,
                -1877.4786213374375,
                21.187576149997447,
                -0.13170006347714705,
                0.0004871198326707739,
                -1.0721367131653266e-6,
                1.3001891159789996e-9,
                -6.701928888007127e-13,
            ],
            b: &[1.0, -0.0019352055158899544],
            max_abs_error: 2.3257113031897028,
            t_min: 200.0,
            t_max: 344.394,
        };
        let cases: [(&RationalPolyAncillary, f64, f64); 12] = [
            (&air_h_l, 60.0, -6344.7615144200545),
            (&air_h_l, 100.0, -4094.961154892782),
            (&air_h_l, 132.0, -1223.449233328278),
            (&air_s_lv, 60.0, 110.68156097433149),
            (&air_s_lv, 100.0, 49.427685581843164),
            (&air_s_lv, 132.0, 12.577271697090302),
            (&r410a_h_lv, 180.0, 25110.785544304796),
            (&r410a_h_lv, 280.0, 15451.502502177267),
            (&r410a_h_lv, 340.0, 5046.893978405953),
            (&r410a_s_l, 180.0, -76.12749488960816),
            (&r410a_s_l, 280.0, -44.16020812742055),
            (&r410a_s_l, 340.0, -17.06481847454429),
        ];
        for (anc, t, expected) in cases {
            let got = evaluate_rational_poly(anc, t);
            assert_eq!(
                got.to_bits(),
                expected.to_bits(),
                "T={t}: got {got:?}, expected {expected:?}"
            );
        }
    }

    /// The `x^2 <= 1` Horner branch (unreachable for kelvin inputs, but part
    /// of the replicated `Eigen::poly_eval`) — synthetic coefficients,
    /// expected value from the same Python replica.
    #[test]
    fn rational_poly_horner_branch_bitwise() {
        let anc = RationalPolyAncillary {
            a: &[1.1, -2.2, 3.3],
            b: &[0.7, 0.3],
            max_abs_error: 0.0,
            t_min: 0.0,
            t_max: 1.0,
        };
        let got = evaluate_rational_poly(&anc, 0.5);
        assert_eq!(got.to_bits(), 0.9705882352941176_f64.to_bits());
    }
}
