//! Region 3 — Helmholtz formulation plus the SR5-05(2016) backward v(T,p)
//! machinery (26 subregions, dividing lines, subregion determination), ported
//! from `Region3`, `Region3Backwards`, and friends in IF97.h.
//!
//! Configuration note (PLAN.md): CoolProp includes IF97.h with neither
//! `IAPWS_UNITS` nor `REGION3_ITERATE` defined, so [`output`] uses the direct
//! (non-iterated) SR5-05 density. The Newton refinement upstream keeps behind
//! `REGION3_ITERATE` is ported as [`rhomass_iterated`] and exercised by the
//! IAPWS check-table tests, which need its tighter accuracy.

use crate::tables::{
    DIV_AB, DIV_CD, DIV_GH, DIV_IJ, DIV_JK, DIV_MN, DIV_OP, DIV_QU, DIV_RX, DIV_UV, DIV_WX,
    Division, R3_A, R3_B, R3_C, R3_D, R3_E, R3_F, R3_G, R3_H, R3_I, R3_J, R3_K, R3_L, R3_M, R3_N,
    R3_O, R3_P, R3_Q, R3_R, R3_S, R3_T, R3_U, R3_V, R3_W, R3_X, R3_Y, R3_Z, REGION3_RESID, Resid,
};
use crate::{P_FACT, Prop, R_FACT, RGAS, RHOCRIT, SatState, TCRIT, powi, region4, transport};
use rustprop_core::Error;

// ---------------------------------------------------------------------------
// Helmholtz phi(delta, tau) machinery (upstream `Region3` methods).
// The first term is n[0]*ln(delta); sums run over terms 1..40.
// ---------------------------------------------------------------------------

fn nr(i: usize) -> &'static Resid {
    &REGION3_RESID[i]
}

pub(crate) fn phi(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = nr(0).n * delta.ln();
    for e in &REGION3_RESID[1..40] {
        summer += e.n * powi(delta, e.i) * powi(tau, e.j);
    }
    summer
}

fn dphi_ddelta(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = nr(0).n / delta;
    for e in &REGION3_RESID[1..40] {
        summer += e.n * f64::from(e.i) * powi(delta, e.i - 1) * powi(tau, e.j);
    }
    summer
}

fn d2phi_ddelta2(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = -nr(0).n / (delta * delta);
    for e in &REGION3_RESID[1..40] {
        summer +=
            e.n * f64::from(e.i) * (f64::from(e.i) - 1.0) * powi(delta, e.i - 2) * powi(tau, e.j);
    }
    summer
}

fn delta_dphi_ddelta(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = nr(0).n;
    for e in &REGION3_RESID[1..40] {
        summer += e.n * f64::from(e.i) * powi(delta, e.i) * powi(tau, e.j);
    }
    summer
}

fn tau_dphi_dtau(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = 0.0;
    for e in &REGION3_RESID[1..40] {
        summer += e.n * f64::from(e.j) * powi(delta, e.i) * powi(tau, e.j);
    }
    summer
}

fn delta2_d2phi_ddelta2(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = -nr(0).n;
    for e in &REGION3_RESID[1..40] {
        summer += e.n * f64::from(e.i) * f64::from(e.i - 1) * powi(delta, e.i) * powi(tau, e.j);
    }
    summer
}

fn tau2_d2phi_dtau2(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = 0.0;
    for e in &REGION3_RESID[1..40] {
        summer += e.n * f64::from(e.j) * f64::from(e.j - 1) * powi(delta, e.i) * powi(tau, e.j);
    }
    summer
}

fn deltatau_d2phi_ddelta_dtau(t: f64, rho: f64) -> f64 {
    let (delta, tau) = (rho / RHOCRIT, TCRIT / t);
    let mut summer = 0.0;
    for e in &REGION3_RESID[1..40] {
        summer += e.n * f64::from(e.j) * f64::from(e.i) * powi(delta, e.i) * powi(tau, e.j);
    }
    summer
}

// ---------------------------------------------------------------------------
// Properties from (T, rho)
// ---------------------------------------------------------------------------

pub fn p(t: f64, rho: f64) -> f64 {
    rho * RGAS * t * delta_dphi_ddelta(t, rho) * (P_FACT / 1000.0 / R_FACT)
}
pub fn umass(t: f64, rho: f64) -> f64 {
    RGAS * t * tau_dphi_dtau(t, rho)
}
pub fn smass(t: f64, rho: f64) -> f64 {
    RGAS * (tau_dphi_dtau(t, rho) - phi(t, rho))
}
pub fn hmass(t: f64, rho: f64) -> f64 {
    RGAS * t * (tau_dphi_dtau(t, rho) + delta_dphi_ddelta(t, rho))
}
pub fn cpmass(t: f64, rho: f64) -> f64 {
    RGAS * (-tau2_d2phi_dtau2(t, rho)
        + powi(
            delta_dphi_ddelta(t, rho) - deltatau_d2phi_ddelta_dtau(t, rho),
            2,
        ) / (2.0 * delta_dphi_ddelta(t, rho) + delta2_d2phi_ddelta2(t, rho)))
}
pub fn cvmass(t: f64, rho: f64) -> f64 {
    RGAS * (-tau2_d2phi_dtau2(t, rho))
}
pub fn speed_sound(t: f64, rho: f64) -> f64 {
    let rhs = 2.0 * delta_dphi_ddelta(t, rho) + delta2_d2phi_ddelta2(t, rho)
        - powi(
            delta_dphi_ddelta(t, rho) - deltatau_d2phi_ddelta_dtau(t, rho),
            2,
        ) / tau2_d2phi_dtau2(t, rho);
    (RGAS * (1000.0 / R_FACT) * t * rhs).sqrt()
}
/// From IAPWS Revised Advisory Note No. 3 — takes rho, not p.
pub fn drhodp(t: f64, rho: f64) -> f64 {
    (rho / p(t, rho)) / (2.0 + delta2_d2phi_ddelta2(t, rho) / delta_dphi_ddelta(t, rho))
}

fn tcond(t: f64, rho: f64) -> f64 {
    let lambda2 = transport::lambda2_region3(
        t,
        rho,
        cpmass(t, rho),
        cvmass(t, rho),
        transport::visc(t, rho),
        drhodp(t, rho),
    );
    0.001 * (transport::lambda0(t) * transport::lambda1(t, rho) + lambda2)
}

// ---------------------------------------------------------------------------
// Newton-Raphson refinement (upstream `REGION3_ITERATE` block)
// ---------------------------------------------------------------------------

fn f(t: f64, p: f64, rho0: f64) -> f64 {
    1.0 / (rho0 * rho0)
        - RGAS * t * dphi_ddelta(t, rho0) / (p * RHOCRIT) * (P_FACT / 1000.0 / R_FACT)
}
fn df(t: f64, p: f64, rho0: f64) -> f64 {
    let rho_c2 = 322.0 * 322.0;
    -2.0 / (rho0 * rho0 * rho0)
        - RGAS * t * d2phi_ddelta2(t, rho0) / (p * rho_c2) * (P_FACT / 1000.0 / R_FACT)
}

/// Newton-Raphson solution of p(T,rho) for rho, seeded with `rho0` (the
/// direct SR5-05 value). Not used by [`output`] — see the module doc.
pub fn rhomass_iterated(t: f64, p: f64, mut rho0: f64) -> Result<f64, Error> {
    let mut iter = 100;
    let mut f_val = f(t, p, rho0);
    while f_val.abs() > 1.0e-14 {
        rho0 -= f_val / df(t, p, rho0);
        iter -= 1;
        if iter == 0 {
            return Err(Error::Solution("Failed to converge!".into()));
        }
        f_val = f(t, p, rho0);
    }
    Ok(rho0)
}

// ---------------------------------------------------------------------------
// SR5-05 backward v(T,p): subregion evaluation
// ---------------------------------------------------------------------------

struct SubRegion {
    v_star: f64,
    p_star: f64,
    t_star: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    data: &'static [Resid],
    /// Subregion N uses `exp(sum)` instead of `sum^e`.
    exp_form: bool,
}

impl SubRegion {
    fn v(&self, t: f64, p: f64) -> f64 {
        let (pi, theta) = (p / self.p_star, t / self.t_star);
        if self.exp_form {
            let mut summer = 0.0;
            for e in self.data {
                summer += e.n * powi(pi - self.a, e.i) * powi(theta - self.b, e.j);
            }
            return summer.exp() * self.v_star;
        }
        let mut summer = 0.0;
        for el in self.data {
            summer += el.n
                * powi((pi - self.a).powf(self.c), el.i)
                * powi((theta - self.b).powf(self.d), el.j);
        }
        summer.powf(self.e) * self.v_star
    }
}

#[rustfmt::skip]
const SUBREGIONS: [SubRegion; 26] = [
    SubRegion { v_star: 0.0024, p_star: 100.0 * P_FACT, t_star: 760.0, a: 0.085, b: 0.817, c: 1.0, d: 1.0, e: 1.0, data: R3_A, exp_form: false },
    SubRegion { v_star: 0.0041, p_star: 100.0 * P_FACT, t_star: 860.0, a: 0.280, b: 0.779, c: 1.0, d: 1.0, e: 1.0, data: R3_B, exp_form: false },
    SubRegion { v_star: 0.0022, p_star: 40.0 * P_FACT, t_star: 690.0, a: 0.259, b: 0.903, c: 1.0, d: 1.0, e: 1.0, data: R3_C, exp_form: false },
    SubRegion { v_star: 0.0029, p_star: 40.0 * P_FACT, t_star: 690.0, a: 0.559, b: 0.939, c: 1.0, d: 1.0, e: 4.0, data: R3_D, exp_form: false },
    SubRegion { v_star: 0.0032, p_star: 40.0 * P_FACT, t_star: 710.0, a: 0.587, b: 0.918, c: 1.0, d: 1.0, e: 1.0, data: R3_E, exp_form: false },
    SubRegion { v_star: 0.0064, p_star: 40.0 * P_FACT, t_star: 730.0, a: 0.587, b: 0.891, c: 0.5, d: 1.0, e: 4.0, data: R3_F, exp_form: false },
    SubRegion { v_star: 0.0027, p_star: 25.0 * P_FACT, t_star: 660.0, a: 0.872, b: 0.971, c: 1.0, d: 1.0, e: 4.0, data: R3_G, exp_form: false },
    SubRegion { v_star: 0.0032, p_star: 25.0 * P_FACT, t_star: 660.0, a: 0.898, b: 0.983, c: 1.0, d: 1.0, e: 4.0, data: R3_H, exp_form: false },
    SubRegion { v_star: 0.0041, p_star: 25.0 * P_FACT, t_star: 660.0, a: 0.910, b: 0.984, c: 0.5, d: 1.0, e: 4.0, data: R3_I, exp_form: false },
    SubRegion { v_star: 0.0054, p_star: 25.0 * P_FACT, t_star: 670.0, a: 0.875, b: 0.964, c: 0.5, d: 1.0, e: 4.0, data: R3_J, exp_form: false },
    SubRegion { v_star: 0.0077, p_star: 25.0 * P_FACT, t_star: 680.0, a: 0.802, b: 0.935, c: 1.0, d: 1.0, e: 1.0, data: R3_K, exp_form: false },
    SubRegion { v_star: 0.0026, p_star: 24.0 * P_FACT, t_star: 650.0, a: 0.908, b: 0.989, c: 1.0, d: 1.0, e: 4.0, data: R3_L, exp_form: false },
    SubRegion { v_star: 0.0028, p_star: 23.0 * P_FACT, t_star: 650.0, a: 1.0, b: 0.997, c: 1.0, d: 0.25, e: 1.0, data: R3_M, exp_form: false },
    SubRegion { v_star: 0.0031, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.976, b: 0.997, c: 0.0, d: 0.0, e: 0.0, data: R3_N, exp_form: true },
    SubRegion { v_star: 0.0034, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.974, b: 0.996, c: 0.5, d: 1.0, e: 1.0, data: R3_O, exp_form: false },
    SubRegion { v_star: 0.0041, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.972, b: 0.997, c: 0.5, d: 1.0, e: 1.0, data: R3_P, exp_form: false },
    SubRegion { v_star: 0.0022, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.848, b: 0.983, c: 1.0, d: 1.0, e: 4.0, data: R3_Q, exp_form: false },
    SubRegion { v_star: 0.0054, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.874, b: 0.982, c: 1.0, d: 1.0, e: 1.0, data: R3_R, exp_form: false },
    SubRegion { v_star: 0.0022, p_star: 21.0 * P_FACT, t_star: 640.0, a: 0.886, b: 0.990, c: 1.0, d: 1.0, e: 4.0, data: R3_S, exp_form: false },
    SubRegion { v_star: 0.0088, p_star: 20.0 * P_FACT, t_star: 650.0, a: 0.803, b: 1.02, c: 1.0, d: 1.0, e: 1.0, data: R3_T, exp_form: false },
    SubRegion { v_star: 0.0026, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.902, b: 0.988, c: 1.0, d: 1.0, e: 1.0, data: R3_U, exp_form: false },
    SubRegion { v_star: 0.0031, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.960, b: 0.995, c: 1.0, d: 1.0, e: 1.0, data: R3_V, exp_form: false },
    SubRegion { v_star: 0.0039, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.959, b: 0.995, c: 1.0, d: 1.0, e: 4.0, data: R3_W, exp_form: false },
    SubRegion { v_star: 0.0049, p_star: 23.0 * P_FACT, t_star: 650.0, a: 0.910, b: 0.988, c: 1.0, d: 1.0, e: 1.0, data: R3_X, exp_form: false },
    SubRegion { v_star: 0.0031, p_star: 22.0 * P_FACT, t_star: 650.0, a: 0.996, b: 0.994, c: 1.0, d: 1.0, e: 4.0, data: R3_Y, exp_form: false },
    SubRegion { v_star: 0.0038, p_star: 22.0 * P_FACT, t_star: 650.0, a: 0.993, b: 0.994, c: 1.0, d: 1.0, e: 4.0, data: R3_Z, exp_form: false },
];

/// Upstream `Region3Backwards::Region3_v_TP`.
pub fn v_tp(region: u8, t: f64, p: f64) -> Result<f64, Error> {
    if !region.is_ascii_uppercase() {
        return Err(Error::OutOfRange("Unable to match region".into()));
    }
    Ok(SUBREGIONS[(region - b'A') as usize].v(t, p))
}

// ---------------------------------------------------------------------------
// Dividing lines between subregions
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Line {
    AB,
    CD,
    EF,
    GH,
    IJ,
    JK,
    MN,
    OP,
    QU,
    RX,
    UV,
    WX,
}

impl Line {
    pub fn parse(name: &str) -> Option<Line> {
        Some(match name {
            "AB" => Line::AB,
            "CD" => Line::CD,
            "EF" => Line::EF,
            "GH" => Line::GH,
            "IJ" => Line::IJ,
            "JK" => Line::JK,
            "MN" => Line::MN,
            "OP" => Line::OP,
            "QU" => Line::QU,
            "RX" => Line::RX,
            "UV" => Line::UV,
            "WX" => Line::WX,
            _ => return None,
        })
    }
}

fn poly_t_p(data: &[Division], p: f64, log_pi: bool) -> f64 {
    let pi = p / (1.0 * P_FACT);
    let x = if log_pi { pi.ln() } else { pi };
    let mut summer = 0.0;
    for e in data {
        summer += e.n * powi(x, e.i);
    }
    summer * 1.0 // sum is multiplied by T* = 1.0 [K]
}

/// Upstream `Region3Backwards::DividingLine`: T [K] on the line at p [Pa].
pub fn dividing_line(line: Line, p: f64) -> f64 {
    match line {
        Line::AB => poly_t_p(DIV_AB, p, true),
        Line::CD => poly_t_p(DIV_CD, p, false),
        Line::EF => {
            let pi = p / (1.0 * P_FACT);
            3.727888004 * (pi - 22.064) + 647.096
        }
        Line::GH => poly_t_p(DIV_GH, p, false),
        Line::IJ => poly_t_p(DIV_IJ, p, false),
        Line::JK => poly_t_p(DIV_JK, p, false),
        Line::MN => poly_t_p(DIV_MN, p, false),
        Line::OP => poly_t_p(DIV_OP, p, true),
        Line::QU => poly_t_p(DIV_QU, p, false),
        Line::RX => poly_t_p(DIV_RX, p, false),
        Line::UV => poly_t_p(DIV_UV, p, false),
        Line::WX => poly_t_p(DIV_WX, p, true),
    }
}

// ---------------------------------------------------------------------------
// Subregion determination
// ---------------------------------------------------------------------------

/// Upstream `BackwardsRegion3SubRegionDetermination` — the near-critical mess.
fn sub_region_determination(t: f64, p: f64) -> Result<u8, Error> {
    let line = |l| dividing_line(l, p);
    if p > 22.5 * P_FACT {
        Err(Error::OutOfRange("Out of range".into()))
    } else if 22.11 * P_FACT < p && p <= 22.5 * P_FACT {
        // Supercritical
        if line(Line::QU) < t && t <= line(Line::UV) {
            Ok(b'U')
        } else if line(Line::UV) < t && t <= line(Line::EF) {
            Ok(b'V')
        } else if line(Line::EF) < t && t <= line(Line::WX) {
            Ok(b'W')
        } else if line(Line::WX) < t && t <= line(Line::RX) {
            Ok(b'X')
        } else {
            Ok(b'?')
        }
    } else if 22.064 * P_FACT < p && p <= 22.11 * P_FACT {
        // Supercritical
        if line(Line::QU) < t && t <= line(Line::UV) {
            Ok(b'U')
        } else if line(Line::UV) < t && t <= line(Line::EF) {
            Ok(b'Y')
        } else if line(Line::EF) < t && t <= line(Line::WX) {
            Ok(b'Z')
        } else if line(Line::WX) < t && t <= line(Line::RX) {
            Ok(b'X')
        } else {
            Ok(b'?')
        }
    } else if t <= region4::t_p(p)? {
        if 21.93161551 * P_FACT < p && p <= 22.064 * P_FACT {
            // Sub-critical
            if line(Line::QU) < t && t <= line(Line::UV) {
                Ok(b'U')
            } else if line(Line::UV) < t {
                Ok(b'Y')
            } else {
                Ok(b'?')
            }
        } else {
            Ok(b'U')
        }
    } else if 21.90096265 * P_FACT < p && p <= 22.064 * P_FACT {
        // Sub-critical
        if t <= line(Line::WX) {
            Ok(b'Z')
        } else if line(Line::WX) < t && t <= line(Line::RX) {
            Ok(b'X')
        } else {
            Ok(b'?')
        }
    } else {
        Ok(b'X')
    }
}

/// Upstream `BackwardsRegion3RegionDetermination`.
#[allow(clippy::excessive_precision)] // boundary literals are verbatim upstream
pub fn region_determination(t: f64, p: f64) -> Result<u8, Error> {
    let line = |l| dividing_line(l, p);
    if p > 100.0 * P_FACT {
        Err(Error::OutOfRange("pressure out of range".into()))
    } else if p > 40.0 * P_FACT && p <= 100.0 * P_FACT {
        if t <= line(Line::AB) {
            Ok(b'A')
        } else {
            Ok(b'B')
        }
    } else if p > 25.0 * P_FACT && p <= 40.0 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= line(Line::AB) {
            Ok(b'D')
        } else if line(Line::AB) < t && t <= line(Line::EF) {
            Ok(b'E')
        } else {
            Ok(b'F')
        }
    } else if p > 23.5 * P_FACT && p <= 25.0 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= line(Line::GH) {
            Ok(b'G')
        } else if line(Line::GH) < t && t <= line(Line::EF) {
            Ok(b'H')
        } else if line(Line::EF) < t && t <= line(Line::IJ) {
            Ok(b'I')
        } else if line(Line::IJ) < t && t <= line(Line::JK) {
            Ok(b'J')
        } else {
            Ok(b'K')
        }
    } else if p > 23.0 * P_FACT && p <= 23.5 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= line(Line::GH) {
            Ok(b'L')
        } else if line(Line::GH) < t && t <= line(Line::EF) {
            Ok(b'H')
        } else if line(Line::EF) < t && t <= line(Line::IJ) {
            Ok(b'I')
        } else if line(Line::IJ) < t && t <= line(Line::JK) {
            Ok(b'J')
        } else {
            Ok(b'K')
        }
    } else if p > 22.5 * P_FACT && p <= 23.0 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= line(Line::GH) {
            Ok(b'L')
        } else if line(Line::GH) < t && t <= line(Line::MN) {
            Ok(b'M')
        } else if line(Line::MN) < t && t <= line(Line::EF) {
            Ok(b'N')
        } else if line(Line::EF) < t && t <= line(Line::OP) {
            Ok(b'O')
        } else if line(Line::OP) < t && t <= line(Line::IJ) {
            Ok(b'P')
        } else if line(Line::IJ) < t && t <= line(Line::JK) {
            Ok(b'J')
        } else {
            Ok(b'K')
        }
    } else if p > 21.04336732 * P_FACT && p <= 22.5 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= line(Line::QU) {
            Ok(b'Q')
        } else if line(Line::RX) < t && t <= line(Line::JK) {
            Ok(b'R')
        } else if t > line(Line::JK) {
            Ok(b'K')
        } else {
            sub_region_determination(t, p)
        }
    } else if p > 20.5 * P_FACT && p <= 21.04336732 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= region4::t_p(p)? {
            Ok(b'S')
        } else if region4::t_p(p)? < t && t <= line(Line::JK) {
            Ok(b'R')
        } else if t > line(Line::JK) {
            Ok(b'K')
        } else {
            Ok(b'?')
        }
    } else if p > 19.00881189173929 * P_FACT && p <= 20.5 * P_FACT {
        if t <= line(Line::CD) {
            Ok(b'C')
        } else if line(Line::CD) < t && t <= region4::t_p(p)? {
            Ok(b'S')
        } else if region4::t_p(p)? < t {
            Ok(b'T')
        } else {
            Ok(b'?')
        }
    } else if p > 16.529164252604481 * P_FACT && p <= 19.00881189173929 * P_FACT {
        if t <= region4::t_p(p)? {
            Ok(b'C')
        } else {
            Ok(b'T')
        }
    } else {
        Ok(b'?')
    }
}

/// Upstream `Region3::SatSubRegionAdjust`: force the subregion to the correct
/// side of the saturation curve when a saturated state is requested.
fn sat_sub_region_adjust(state: SatState, p: f64, subregion: u8) -> u8 {
    match state {
        SatState::Vapor => match subregion {
            b'C' => b'T',
            b'S' => {
                if p < 20.5 * P_FACT {
                    b'T'
                } else {
                    b'R'
                }
            }
            b'U' => {
                if p < 21.90096265 * P_FACT {
                    b'X'
                } else {
                    b'Z'
                }
            }
            b'Y' => b'Z',
            _ => subregion,
        },
        SatState::Liquid => match subregion {
            b'Z' => {
                if p > 21.93161551 * P_FACT {
                    b'Y'
                } else {
                    b'U'
                }
            }
            b'X' => b'U',
            b'R' | b'K' => b'S',
            b'T' => {
                if p > 19.00881189173929 * P_FACT {
                    b'S'
                } else {
                    b'C'
                }
            }
            _ => subregion,
        },
        SatState::None => subregion,
    }
}

/// Upstream `Region3::output` in CoolProp's configuration: density from the
/// direct SR5-05 backward v(T,p), no Newton refinement.
pub(crate) fn output(key: Prop, t: f64, p_in: f64, state: SatState) -> Result<f64, Error> {
    let region = region_determination(t, p_in)?;
    let region = sat_sub_region_adjust(state, p_in, region);
    let rho = 1.0 / v_tp(region, t, p_in)?;
    match key {
        Prop::Dmass => Ok(rho),
        Prop::Hmass => Ok(hmass(t, rho)),
        Prop::Smass => Ok(smass(t, rho)),
        Prop::Umass => Ok(umass(t, rho)),
        Prop::Cpmass => Ok(cpmass(t, rho)),
        Prop::Cvmass => Ok(cvmass(t, rho)),
        Prop::W => Ok(speed_sound(t, rho)),
        Prop::Mu => Ok(transport::visc(t, rho)),
        Prop::K => Ok(tcond(t, rho)),
        Prop::DrhoDp => Ok(drhodp(t, rho)),
        _ => Err(Error::Input("Bad key to output".into())),
    }
}
