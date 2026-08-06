//! Transport properties — viscosity (IAPWS R12-08, industrial: mu2 = 1) and
//! thermal conductivity (IAPWS R15-11, simplified critical enhancement),
//! ported from the shared machinery in IF97.h's `BaseRegion` and `Region3`.
//!
//! The two `lambda2` variants deliberately differ exactly as upstream does:
//! the Gibbs-region version uses the literal `PI = 3.141592654` and gates on
//! `delChi > 0`; the region-3 version uses exact pi (`2*acos(0)`) and clamps
//! `zeta` instead. `del_tr` likewise exists in both power-evaluation orders.

use crate::tables::{COND_CRIT_A, COND_IDEAL, COND_RESID, VISC_IDEAL, VISC_RESID};
use crate::{PCRIT, R_FACT, RHOCRIT, TCRIT, powi};

fn tr_term(t: f64) -> f64 {
    TCRIT / t - 1.0
}
fn rhor_term(rho: f64) -> f64 {
    rho / RHOCRIT - 1.0
}

fn mu0(t: f64) -> f64 {
    let t_bar = t / TCRIT;
    let summer: f64 = VISC_IDEAL.iter().map(|e| e.n / powi(t_bar, e.j)).sum();
    100.0 * t_bar.sqrt() / summer
}

fn mu1(t: f64, rho: f64) -> f64 {
    let rho_bar = rho / RHOCRIT;
    let summer: f64 = VISC_RESID
        .iter()
        .map(|e| rho_bar * powi(tr_term(t), e.i) * e.n * powi(rhor_term(rho), e.j))
        .sum();
    summer.exp()
}

/// Viscosity [Pa-s] from T [K] and rho [kg/m^3]; valid in every region.
pub(crate) fn visc(t: f64, rho: f64) -> f64 {
    let mu_star = 1.0e-6; // Reference viscosity [Pa-s]
    let mu2 = 1.0; // For Industrial Formulation (IF97), mu2 = 1.0
    mu_star * mu0(t) * mu1(t, rho) * mu2
}

pub(crate) fn lambda0(t: f64) -> f64 {
    let t_bar = t / TCRIT;
    let summer: f64 = COND_IDEAL.iter().map(|e| e.n / powi(t_bar, e.j)).sum();
    t_bar.sqrt() / summer
}

pub(crate) fn lambda1(t: f64, rho: f64) -> f64 {
    let rho_bar = rho / RHOCRIT;
    let summer: f64 = COND_RESID
        .iter()
        .map(|e| rho_bar * powi(tr_term(t), e.i) * e.n * powi(rhor_term(rho), e.j))
        .sum();
    summer.exp()
}

/// drhodp correlation at the reducing temperature — Gibbs-region variant
/// (`BaseRegion::delTr`, powers via `powi`).
fn del_tr_gibbs(rho: f64) -> f64 {
    let rhobar = rho / RHOCRIT;
    let j = del_tr_column(rhobar);
    let mut summer = 0.0;
    for (i, row) in COND_CRIT_A.iter().enumerate() {
        summer += row[j] * powi(rhobar, i as i32);
    }
    1.0 / summer
}

/// drhodp correlation at the reducing temperature — region-3 variant
/// (`Region3::delTr`, powers via a running product).
pub(crate) fn del_tr_region3(rho: f64) -> f64 {
    let rhobar = rho / RHOCRIT;
    let j = del_tr_column(rhobar);
    let mut summer = 0.0;
    let mut pow_rhobar = 1.0;
    for row in COND_CRIT_A.iter() {
        summer += row[j] * pow_rhobar;
        pow_rhobar *= rhobar;
    }
    1.0 / summer
}

fn del_tr_column(rhobar: f64) -> usize {
    if rhobar <= 0.310559006 {
        0
    } else if rhobar <= 0.776397516 {
        1
    } else if rhobar <= 1.242236025 {
        2
    } else if rhobar <= 1.863354037 {
        3
    } else {
        4
    }
}

/// Critical enhancement for regions 1 and 2 (`BaseRegion::lambda2`).
#[allow(clippy::approx_constant)] // upstream uses the truncated 3.141592654
pub(crate) fn lambda2_gibbs(t: f64, rho: f64, cp: f64, cv: f64, mu: f64, drhodp: f64) -> f64 {
    let rhobar = rho / RHOCRIT;
    let lambda_cap = 177.8514;
    let q_d = 1.0 / 0.40;
    let tr = 1.5 * TCRIT;
    let xi0 = 0.13;
    let nu = 0.630;
    let gam = 1.239;
    let gamma0 = 0.06;
    let pi = 3.141592654; // upstream literal
    let cp_star = 0.46151805 * R_FACT; // slightly lower than IF97 Rgas
    let mut cpbar = cp / cp_star;
    if cpbar < 0.0 || cpbar > 1.0e13 {
        cpbar = 1.0e13;
    }
    let k = cp / cv;
    let mubar = mu / 1.0e-6;
    let del_chi = rhobar * (PCRIT / RHOCRIT * drhodp - del_tr_gibbs(rho) * tr / t);
    // At low (T,p), delChi can go negative, making y imaginary from the
    // nth-root; limit to delChi > 0.
    let y = if del_chi > 0.0 {
        q_d * xi0 * (del_chi / gamma0).powf(nu / gam)
    } else {
        0.0
    };
    let z = if y < 1.2e-7 {
        0.0
    } else {
        2.0 / pi / y
            * (((1.0 - 1.0 / k) * y.atan() + y / k)
                - (1.0 - (-1.0 / (1.0 / y + y * y / (3.0 * rhobar * rhobar))).exp()))
    };
    lambda_cap * rhobar * cpbar * t / (TCRIT * mubar) * z
}

/// Critical enhancement for region 3 (`Region3::lambda2`).
pub(crate) fn lambda2_region3(t: f64, rho: f64, cp: f64, cv: f64, mu: f64, drhodp: f64) -> f64 {
    let rhobar = rho / RHOCRIT;
    let lambda_cap = 177.8514;
    let q_d = 1.0 / 0.40;
    let tr = 1.5 * TCRIT;
    let xi0 = 0.13;
    let nu = 0.630;
    let gam = 1.239;
    let gamma0 = 0.06;
    let pi = std::f64::consts::PI; // upstream: 2*acos(0)
    let cp_star = 0.46151805 * R_FACT;
    let mut cpbar = cp / cp_star;
    if cpbar < 0.0 || cpbar > 1.0e13 {
        cpbar = 1.0e13;
    }
    let k = cp / cv;
    let mubar = mu / 1.0e-6;
    let mut zeta = PCRIT / RHOCRIT * drhodp;
    if zeta < 0.0 || zeta > 1.0e13 {
        zeta = 1.0e13;
    }
    let del_chi = rhobar * (zeta - del_tr_region3(rho) * tr / t);
    let y = q_d * xi0 * (del_chi / gamma0).powf(nu / gam);
    let z = if y < 1.2e-7 {
        0.0
    } else {
        2.0 / (pi * y)
            * (((1.0 - 1.0 / k) * y.atan() + y / k)
                - (1.0 - (-1.0 / (1.0 / y + y * y / (3.0 * rhobar * rhobar))).exp()))
    };
    lambda_cap * rhobar * cpbar * t / (TCRIT * mubar) * z
}
