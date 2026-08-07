//! HEOS mixture machinery (PLAN.md Phase 10) — the GERG-2008 reducing
//! function (upstream `GERG2008ReducingFunction`, used for every HEOS
//! mixture), through the second-order composition derivatives the PT/PQ/QT
//! flashes consume. Third-order machinery (`d3Yr*`, `PSI_*`) is envelope/
//! critical-locus territory and deliberately unported.
//!
//! Upstream conventions reproduced:
//! - The constructor cross terms use each fluid's `EOS().reduce` state (NOT
//!   `crit`): `T_c[i][j] = sqrt(Tr_i Tr_j)`,
//!   `v_c[i][j] = 1/8 (rhor_i^{-1/3} + rhor_j^{-1/3})^3`.
//! - `beta` is antisymmetric under component swap (`beta[j][i] =
//!   1/beta[i][j]`), `gamma` symmetric.
//! - The quadratic form is `Y_r = sum_i x_i^2 Y_ci + cross terms` (the
//!   upstream Doxygen writes `x_i Y_ci^2` — a doc typo; the code is
//!   authoritative).

// Function names mirror upstream's `dYrdxi__constxj` convention (the double
// underscore separates the derivative from its held-constant list).
#![allow(non_snake_case)]

use rustprop_core::fluid::FluidData;
use rustprop_core::{Error, Result};

/// Which mole fraction convention the composition derivatives use
/// (upstream `x_N_dependency_flag`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum XnFlag {
    /// Kunz-Wagner 2012 Table B9: all N fractions independent.
    Independent,
    /// Gernert 2014 Table S1: `x_N = 1 - sum x_i`.
    Dependent,
}

/// The GERG-2008 reducing function for one component set.
pub struct Gerg2008Reducing {
    n: usize,
    /// `beta_T[i][j]`, antisymmetric.
    beta_t: Vec<f64>,
    gamma_t: Vec<f64>,
    beta_v: Vec<f64>,
    gamma_v: Vec<f64>,
    /// `T_c[i][j] = sqrt(Tr_i Tr_j)`.
    t_c: Vec<f64>,
    /// `v_c[i][j]`.
    v_c: Vec<f64>,
    /// Diagonal `Yc` vectors.
    yc_t: Vec<f64>,
    yc_v: Vec<f64>,
}

impl Gerg2008Reducing {
    /// Build from the component list and the CAS-keyed binary-pair table
    /// (upstream `MixtureParameters::set_mixture_parameters`): every i != j
    /// pair must exist or construction fails with upstream's message.
    pub fn new(
        components: &[&'static FluidData],
        pairs: &[rustprop_core::fluid::MixBinaryPair],
    ) -> Result<Self> {
        let n = components.len();
        let mut beta_t = vec![1.0; n * n];
        let mut gamma_t = vec![1.0; n * n];
        let mut beta_v = vec![1.0; n * n];
        let mut gamma_v = vec![1.0; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let (a, b) = (components[i].cas, components[j].cas);
                // The table ships in document order; match either
                // orientation (upstream sorts the lookup key, we search).
                let pair = pairs
                    .iter()
                    .find(|p| (p.cas1 == a && p.cas2 == b) || (p.cas1 == b && p.cas2 == a))
                    .ok_or_else(|| {
                        let mut sorted = [a, b];
                        sorted.sort_unstable();
                        Error::Value(format!(
                            "Could not match the binary pair [{},{}] - for now this is an error.",
                            sorted[0], sorted[1]
                        ))
                    })?;
                // beta inverts when the live component order is the swap of
                // the stored orientation.
                let stored_forward = a == pair.cas1;
                let (bt, bv) = if stored_forward {
                    (pair.beta_t, pair.beta_v)
                } else {
                    (1.0 / pair.beta_t, 1.0 / pair.beta_v)
                };
                beta_t[i * n + j] = bt;
                beta_t[j * n + i] = 1.0 / bt;
                gamma_t[i * n + j] = pair.gamma_t;
                gamma_t[j * n + i] = pair.gamma_t;
                beta_v[i * n + j] = bv;
                beta_v[j * n + i] = 1.0 / bv;
                gamma_v[i * n + j] = pair.gamma_v;
                gamma_v[j * n + i] = pair.gamma_v;
            }
        }
        let mut t_c = vec![0.0; n * n];
        let mut v_c = vec![0.0; n * n];
        let mut yc_t = vec![0.0; n];
        let mut yc_v = vec![0.0; n];
        for i in 0..n {
            yc_t[i] = components[i].eos.reducing.t;
            yc_v[i] = 1.0 / components[i].eos.reducing.rhomolar;
            for j in 0..n {
                t_c[i * n + j] =
                    (components[i].eos.reducing.t * components[j].eos.reducing.t).sqrt();
                v_c[i * n + j] = 1.0 / 8.0
                    * (components[i].eos.reducing.rhomolar.powf(-1.0 / 3.0)
                        + components[j].eos.reducing.rhomolar.powf(-1.0 / 3.0))
                    .powi(3);
            }
        }
        Ok(Gerg2008Reducing {
            n,
            beta_t,
            gamma_t,
            beta_v,
            gamma_v,
            t_c,
            v_c,
            yc_t,
            yc_v,
        })
    }

    fn c_y_ij(&self, i: usize, j: usize, beta: &[f64], gamma: &[f64], y_c: &[f64]) -> f64 {
        2.0 * beta[i * self.n + j] * gamma[i * self.n + j] * y_c[i * self.n + j]
    }

    fn f_y_ij(x: &[f64], i: usize, j: usize, beta_y: f64) -> f64 {
        let (xi, xj) = (x[i], x[j]);
        xi * xj * (xi + xj) / (beta_y * beta_y * xi + xj)
    }

    /// The generalized quadratic form (upstream `Yr`).
    fn yr(&self, x: &[f64], beta: &[f64], gamma: &[f64], y_c_ij: &[f64], yc: &[f64]) -> f64 {
        let n = self.n;
        let mut yr = 0.0;
        for i in 0..n {
            yr += x[i] * x[i] * yc[i];
            for j in (i + 1)..n {
                yr +=
                    self.c_y_ij(i, j, beta, gamma, y_c_ij) * Self::f_y_ij(x, i, j, beta[i * n + j]);
            }
        }
        yr
    }

    pub fn tr(&self, x: &[f64]) -> f64 {
        self.yr(x, &self.beta_t, &self.gamma_t, &self.t_c, &self.yc_t)
    }

    pub fn rhormolar(&self, x: &[f64]) -> f64 {
        1.0 / self.yr(x, &self.beta_v, &self.gamma_v, &self.v_c, &self.yc_v)
    }

    // --- f_Y_ij building blocks (upstream ReducingFunctions.cpp:550-609) ---

    fn dfykidxi__constxk(x: &[f64], k: usize, i: usize, beta_y: f64) -> f64 {
        let (xk, xi) = (x[k], x[i]);
        let b2 = beta_y * beta_y;
        xk * (xk + xi) / (b2 * xk + xi)
            + xk * xi / (b2 * xk + xi) * (1.0 - (xk + xi) / (b2 * xk + xi))
    }
    fn dfyikdxi__constxk(x: &[f64], i: usize, k: usize, beta_y: f64) -> f64 {
        let (xi, xk) = (x[i], x[k]);
        let b2 = beta_y * beta_y;
        xk * (xi + xk) / (b2 * xi + xk)
            + xi * xk / (b2 * xi + xk) * (1.0 - b2 * (xi + xk) / (b2 * xi + xk))
    }
    fn d2fyikdxi2__constxk(x: &[f64], i: usize, k: usize, beta_y: f64) -> f64 {
        let (xi, xk) = (x[i], x[k]);
        let b2 = beta_y * beta_y;
        1.0 / (b2 * xi + xk)
            * (1.0 - b2 * (xi + xk) / (b2 * xi + xk))
            * (2.0 * xk - 2.0 * b2 * xi * xk / (b2 * xi + xk))
    }
    fn d2fykidxi2__constxk(x: &[f64], k: usize, i: usize, beta_y: f64) -> f64 {
        let (xk, xi) = (x[k], x[i]);
        let b2 = beta_y * beta_y;
        1.0 / (b2 * xk + xi)
            * (1.0 - (xk + xi) / (b2 * xk + xi))
            * (2.0 * xk - 2.0 * xk * xi / (b2 * xk + xi))
    }
    fn d2fyijdxidxj(x: &[f64], i: usize, j: usize, beta_y: f64) -> f64 {
        let (xi, xj) = (x[i], x[j]);
        let b2 = beta_y * beta_y;
        (xi + xj) / (b2 * xi + xj)
            + xj / (b2 * xi + xj) * (1.0 - (xi + xj) / (b2 * xi + xj))
            + xi / (b2 * xi + xj) * (1.0 - b2 * (xi + xj) / (b2 * xi + xj))
            - xi * xj / (b2 * xi + xj).powi(2) * (1.0 + b2 - 2.0 * b2 * (xi + xj) / (b2 * xi + xj))
    }

    // --- first derivatives (upstream Yr composition derivatives) ---

    #[allow(clippy::too_many_arguments)]
    fn dyrdxi__constxj(
        &self,
        x: &[f64],
        i: usize,
        beta: &[f64],
        gamma: &[f64],
        y_c_ij: &[f64],
        yc: &[f64],
        flag: XnFlag,
    ) -> f64 {
        let n = self.n;
        match flag {
            XnFlag::Independent => {
                let mut d = 2.0 * x[i] * yc[i];
                for k in 0..i {
                    d += self.c_y_ij(k, i, beta, gamma, y_c_ij)
                        * Self::dfykidxi__constxk(x, k, i, beta[k * n + i]);
                }
                for k in (i + 1)..n {
                    d += self.c_y_ij(i, k, beta, gamma, y_c_ij)
                        * Self::dfyikdxi__constxk(x, i, k, beta[i * n + k]);
                }
                d
            }
            XnFlag::Dependent => {
                if i == n - 1 {
                    return 0.0;
                }
                let xn = x[n - 1];
                let mut d = 2.0 * x[i] * yc[i] - 2.0 * xn * yc[n - 1];
                for k in 0..i {
                    d += self.c_y_ij(k, i, beta, gamma, y_c_ij)
                        * Self::dfykidxi__constxk(x, k, i, beta[k * n + i]);
                }
                for k in (i + 1)..(n - 1) {
                    d += self.c_y_ij(i, k, beta, gamma, y_c_ij)
                        * Self::dfyikdxi__constxk(x, i, k, beta[i * n + k]);
                }
                let b_in = beta[i * n + (n - 1)];
                let b2 = b_in * b_in;
                d += self.c_y_ij(i, n - 1, beta, gamma, y_c_ij)
                    * (xn * (x[i] + xn) / (b2 * x[i] + xn)
                        + (1.0 - b2) * x[i] * xn * xn / (b2 * x[i] + xn).powi(2));
                for k in 0..(n - 1) {
                    let b_kn = beta[k * n + (n - 1)];
                    let bk2 = b_kn * b_kn;
                    d += self.c_y_ij(k, n - 1, beta, gamma, y_c_ij)
                        * (-x[k] * (x[k] + xn) / (bk2 * x[k] + xn)
                            + (1.0 - bk2) * xn * x[k] * x[k] / (bk2 * x[k] + xn).powi(2));
                }
                d
            }
        }
    }

    pub fn dtrdxi__constxj(&self, x: &[f64], i: usize, flag: XnFlag) -> f64 {
        self.dyrdxi__constxj(
            x,
            i,
            &self.beta_t,
            &self.gamma_t,
            &self.t_c,
            &self.yc_t,
            flag,
        )
    }
    pub fn dvrmolardxi__constxj(&self, x: &[f64], i: usize, flag: XnFlag) -> f64 {
        self.dyrdxi__constxj(
            x,
            i,
            &self.beta_v,
            &self.gamma_v,
            &self.v_c,
            &self.yc_v,
            flag,
        )
    }
    /// GERG 2004 Eq. 7.57: `drhor/dxi = -rhor^2 dvr/dxi`.
    pub fn drhormolardxi__constxj(&self, x: &[f64], i: usize, flag: XnFlag) -> f64 {
        let rhor = self.rhormolar(x);
        -rhor * rhor * self.dvrmolardxi__constxj(x, i, flag)
    }

    // --- second derivatives ---

    #[allow(clippy::too_many_arguments)]
    fn d2yrdxidxj(
        &self,
        x: &[f64],
        i: usize,
        j: usize,
        beta: &[f64],
        gamma: &[f64],
        y_c_ij: &[f64],
        yc: &[f64],
        flag: XnFlag,
    ) -> f64 {
        let n = self.n;
        if i == j {
            // d2Yr/dxi2 (upstream d2Yrdxi2__constxj), XN_INDEPENDENT form.
            let mut d = 2.0 * yc[i];
            for k in 0..i {
                d += self.c_y_ij(k, i, beta, gamma, y_c_ij)
                    * Self::d2fykidxi2__constxk(x, k, i, beta[k * n + i]);
            }
            for k in (i + 1)..n {
                d += self.c_y_ij(i, k, beta, gamma, y_c_ij)
                    * Self::d2fyikdxi2__constxk(x, i, k, beta[i * n + k]);
            }
            d
        } else {
            let _ = flag;
            self.c_y_ij(i.min(j), i.max(j), beta, gamma, y_c_ij)
                * Self::d2fyijdxidxj(x, i.min(j), i.max(j), beta[i.min(j) * n + i.max(j)])
        }
    }

    pub fn d2trdxidxj(&self, x: &[f64], i: usize, j: usize, flag: XnFlag) -> f64 {
        self.d2yrdxidxj(
            x,
            i,
            j,
            &self.beta_t,
            &self.gamma_t,
            &self.t_c,
            &self.yc_t,
            flag,
        )
    }
    pub fn d2vrmolardxidxj(&self, x: &[f64], i: usize, j: usize, flag: XnFlag) -> f64 {
        self.d2yrdxidxj(
            x,
            i,
            j,
            &self.beta_v,
            &self.gamma_v,
            &self.v_c,
            &self.yc_v,
            flag,
        )
    }
    /// GERG 2004 Eqns 7.58/7.59.
    pub fn d2rhormolardxidxj(&self, x: &[f64], i: usize, j: usize, flag: XnFlag) -> f64 {
        let rhor = self.rhormolar(x);
        let dvi = self.dvrmolardxi__constxj(x, i, flag);
        let dvj = self.dvrmolardxi__constxj(x, j, flag);
        2.0 * rhor.powi(3) * dvi * dvj - rhor * rhor * self.d2vrmolardxidxj(x, i, j, flag)
    }

    // --- composite n d/dn forms (upstream base class) ---

    /// GERG Eq. 7.54.
    pub fn ndtrdni__constnj(&self, x: &[f64], i: usize, flag: XnFlag) -> f64 {
        let lim = match flag {
            XnFlag::Independent => self.n,
            XnFlag::Dependent => self.n - 1,
        };
        let mut summer = 0.0;
        for j in 0..lim {
            summer += x[j] * self.dtrdxi__constxj(x, j, flag);
        }
        self.dtrdxi__constxj(x, i, flag) - summer
    }
    pub fn ndrhorbardni__constnj(&self, x: &[f64], i: usize, flag: XnFlag) -> f64 {
        let lim = match flag {
            XnFlag::Independent => self.n,
            XnFlag::Dependent => self.n - 1,
        };
        let mut summer = 0.0;
        for j in 0..lim {
            summer += x[j] * self.drhormolardxi__constxj(x, j, flag);
        }
        self.drhormolardxi__constxj(x, i, flag) - summer
    }

    /// GERG Eq. 7.56 / Gernert A28.
    pub fn d_ndtrdni_dxj__constxi(&self, x: &[f64], i: usize, j: usize, flag: XnFlag) -> f64 {
        match flag {
            XnFlag::Independent => {
                let mut summer = 0.0;
                for k in 0..self.n {
                    summer += x[k] * self.d2trdxidxj(x, j, k, flag);
                }
                self.d2trdxidxj(x, i, j, flag) - self.dtrdxi__constxj(x, j, flag) - summer
            }
            XnFlag::Dependent => {
                if j == self.n - 1 {
                    return 0.0;
                }
                let mut summer = 0.0;
                for k in 0..(self.n - 1) {
                    summer += x[k] * self.d2trdxidxj(x, k, j, flag);
                }
                self.d2trdxidxj(x, j, i, flag) - self.dtrdxi__constxj(x, j, flag) - summer
            }
        }
    }
    pub fn d_ndrhorbardni_dxj__constxi(&self, x: &[f64], i: usize, j: usize, flag: XnFlag) -> f64 {
        // Upstream has no xN branch here.
        let mut summer = 0.0;
        for k in 0..self.n {
            summer += x[k] * self.d2rhormolardxidxj(x, j, k, flag);
        }
        self.d2rhormolardxidxj(x, j, i, flag) - self.drhormolardxi__constxj(x, j, flag) - summer
    }
}
