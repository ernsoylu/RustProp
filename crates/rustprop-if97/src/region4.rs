//! Region 4 (saturation curve) and surface tension — port of `Region4` from
//! IF97.h, including the loop-restructured evaluation of the quadratic
//! saturation equations (kept verbatim for bit-order fidelity).

use crate::tables::REGION4_N;
use crate::{P_FACT, PCRIT, PMIN, TCRIT, TMIN, TTRIP, powi};
use rustprop_core::Error;

/// Upstream indexes coefficients 1..=10 with `n[0] = 0`.
fn n(i: usize) -> f64 {
    if i == 0 { 0.0 } else { REGION4_N[i - 1] }
}

const P_STAR: f64 = 1.0 * P_FACT;
const T_STAR: f64 = 1.0;

/// Saturation pressure [Pa] from temperature [K] (upstream `Region4::p_T`).
pub(crate) fn p_t(t: f64) -> Result<f64, Error> {
    if t < TMIN || t > TCRIT {
        return Err(Error::OutOfRange("Temperature out of range".into()));
    }
    let theta = t / T_STAR + n(9) / (t / T_STAR - n(10));
    let mut abc = [1.0, n(3), n(6)];
    for j in 1..3 {
        for x in abc.iter_mut() {
            *x *= theta;
        }
        for (i, x) in abc.iter_mut().enumerate() {
            *x += n(i * 3 + j);
        }
    }
    let (aa, bb, cc) = (abc[0], abc[1], abc[2]);
    Ok(P_STAR * powi(2.0 * cc / (-bb + (bb * bb - 4.0 * aa * cc).sqrt()), 4))
}

/// Saturation temperature [K] from pressure [Pa] (upstream `Region4::T_p`).
pub(crate) fn t_p(p: f64) -> Result<f64, Error> {
    if p < PMIN || p > PCRIT {
        return Err(Error::OutOfRange("Pressure out of range".into()));
    }
    let beta2 = (p / P_STAR).sqrt();
    let beta = beta2.sqrt();
    let mut efg = [1.0, n(1), n(2)];
    for x in efg.iter_mut() {
        *x *= beta;
    }
    for (i, x) in efg.iter_mut().enumerate() {
        *x += n(i + 3);
    }
    for x in efg.iter_mut() {
        *x *= beta;
    }
    for (i, x) in efg.iter_mut().enumerate() {
        *x += n(i + 6);
    }
    let (e, f, g) = (efg[0], efg[1], efg[2]);
    let d = 2.0 * g / (-f - (f * f - 4.0 * e * g).sqrt());
    let n10pd = n(10) + d;
    Ok(T_STAR * 0.5 * (n10pd - (n10pd * n10pd - 4.0 * (n(9) + n(10) * d)).sqrt()))
}

/// Surface tension [N/m] from temperature [K], IAPWS R1-76(2014)
/// (upstream `Region4::sigma_t`). Extrapolates down to 25 K below triple.
pub(crate) fn sigma_t(t: f64) -> Result<f64, Error> {
    if t < (TTRIP - 25.0) || t > TCRIT {
        return Err(Error::OutOfRange("Temperature out of range".into()));
    }
    let tau = 1.0 - t / TCRIT;
    let b_cap = 235.8 / 1000.0; // published in mN/m; SI here in all cases
    let b = -0.625;
    let mu = 1.256;
    Ok(b_cap * tau.powf(mu) * (1.0 + b * tau))
}
