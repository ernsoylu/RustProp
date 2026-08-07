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
                    .powf(3.0);
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
            * (2.0 * xk - xi * xk * 2.0 * b2 / (b2 * xi + xk))
    }
    fn d2fykidxi2__constxk(x: &[f64], k: usize, i: usize, beta_y: f64) -> f64 {
        let (xk, xi) = (x[k], x[i]);
        let b2 = beta_y * beta_y;
        1.0 / (b2 * xk + xi)
            * (1.0 - (xk + xi) / (b2 * xk + xi))
            * (2.0 * xk - xk * xi * 2.0 / (b2 * xk + xi))
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
    /// Upstream `d2Yrdxi2__constxj` — Kunz-Wagner Table B9 (Independent) /
    /// Gernert Table S1 (Dependent).
    #[allow(clippy::too_many_arguments)]
    fn d2yrdxi2__constxj(
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
            }
            XnFlag::Dependent => {
                if i == n - 1 {
                    return 0.0;
                }
                let xn = x[n - 1];
                let mut d = 2.0 * yc[i] + 2.0 * yc[n - 1];
                for k in 0..i {
                    d += self.c_y_ij(k, i, beta, gamma, y_c_ij)
                        * Self::d2fykidxi2__constxk(x, k, i, beta[k * n + i]);
                }
                for k in (i + 1)..(n - 1) {
                    d += self.c_y_ij(i, k, beta, gamma, y_c_ij)
                        * Self::d2fyikdxi2__constxk(x, i, k, beta[i * n + k]);
                }
                let b_in = beta[i * n + (n - 1)];
                let b2 = b_in * b_in;
                d += 2.0
                    * self.c_y_ij(i, n - 1, beta, gamma, y_c_ij)
                    * (-(x[i] + xn) / (b2 * x[i] + xn)
                        + (1.0 - b2)
                            * (xn * xn / (b2 * x[i] + xn).powi(2)
                                + ((1.0 - b2) * x[i] * xn * xn - b2 * x[i] * x[i] * xn)
                                    / (b2 * x[i] + xn).powi(3)));
                for k in 0..(n - 1) {
                    let b_kn = beta[k * n + (n - 1)];
                    let bk2 = b_kn * b_kn;
                    d += 2.0
                        * self.c_y_ij(k, n - 1, beta, gamma, y_c_ij)
                        * x[k]
                        * x[k]
                        * (1.0 - bk2)
                        / (bk2 * x[k] + xn).powi(2)
                        * (xn / (bk2 * x[k] + xn) - 1.0);
                }
                d
            }
        }
    }

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
        match flag {
            XnFlag::Independent => {
                if i == j {
                    self.d2yrdxi2__constxj(x, i, beta, gamma, y_c_ij, yc, flag)
                } else {
                    // See Table B9 from Kunz Wagner 2012 (GERG 2008)
                    self.c_y_ij(i, j, beta, gamma, y_c_ij)
                        * Self::d2fyijdxidxj(x, i, j, beta[i * n + j])
                }
            }
            XnFlag::Dependent => {
                // Table S1 from Gernert, 2014, supplemental information
                if j == n - 1 || i == n - 1 {
                    return 0.0;
                }
                if i == j {
                    return self.d2yrdxi2__constxj(x, i, beta, gamma, y_c_ij, yc, flag);
                }
                let mut d = 2.0 * yc[n - 1];
                d += self.c_y_ij(i, j, beta, gamma, y_c_ij)
                    * Self::d2fyijdxidxj(x, i, j, beta[i * n + j]);
                for k in 0..(n - 1) {
                    d += self.c_y_ij(k, n - 1, beta, gamma, y_c_ij)
                        * Self::d2fykidxi2__constxk(x, k, n - 1, beta[k * n + (n - 1)]);
                }
                d -= self.c_y_ij(i, n - 1, beta, gamma, y_c_ij)
                    * Self::d2fyijdxidxj(x, i, n - 1, beta[i * n + (n - 1)]);
                d -= self.c_y_ij(j, n - 1, beta, gamma, y_c_ij)
                    * Self::d2fyijdxidxj(x, j, n - 1, beta[j * n + (n - 1)]);
                d
            }
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

// ---------------------------------------------------------------------------
// Mixture Helmholtz model (slice 10c)
// ---------------------------------------------------------------------------

/// Upstream `R_U_CODATA` (Configuration.h) — `calc_gas_constant` returns it
/// for every mixture because `NORMALIZE_GAS_CONSTANTS` defaults to true.
pub const R_U_CODATA: f64 = 8.314_462_618_153_24;

/// One `i < j` pair of the excess term: the scaling factor `F[i][j]` and the
/// departure function (upstream `ExcessTerm`'s F matrix +
/// `DepartureFunctionMatrix`).
struct ExcessPair {
    i: usize,
    j: usize,
    f: f64,
    dep: crate::alpha::DepartureEval,
}

/// The Helmholtz-energy model of one mixture: per-component pure containers
/// (corresponding-states part), the GERG-2008 reducing function, and the
/// excess term (upstream `HelmholtzEOSMixtureBackend`'s `residual_helmholtz`
/// + `Reducing` + the Table B5 ideal part).
pub struct MixtureModel {
    /// Per-component full pure-fluid containers evaluated at the MIXTURE
    /// (tau, delta) in the corresponding-states sum.
    pure: Vec<crate::alpha::HelmholtzEos>,
    /// `STATES.critical` T [K] per component (upstream `iT_critical`) — the
    /// Table B5 alpha0 scales, NOT the reducing state.
    crit_t: Vec<f64>,
    /// `STATES.critical` rhomolar [mol/m^3] per component.
    crit_rhomolar: Vec<f64>,
    /// `STATES.critical` p [Pa] per component (upstream `iP_critical`).
    crit_p: Vec<f64>,
    /// `EOS().sat_min_liquid.T` per component (upstream `iT_triple`).
    triple_t: Vec<f64>,
    /// `EOS().sat_min_liquid.p` per component (upstream `iP_triple`).
    triple_p: Vec<f64>,
    /// `EOS().R_u` per component.
    r_component: Vec<f64>,
    /// `(EOS().reduce.T, EOS().reduce.p, EOS().acentric)` per component —
    /// the SRK seed constants.
    srk_consts: Vec<(f64, f64, f64)>,
    /// `EOS().molar_mass` per component.
    molar_mass: Vec<f64>,
    pub reducing: Gerg2008Reducing,
    excess: Vec<ExcessPair>,
}

impl MixtureModel {
    /// Build from the component list and the datagen tables (the caller
    /// passes `rustprop_data::mixtures::{MIX_BINARY_PAIRS, MIX_DEPARTURE_FNS}`
    /// — engines only depend on core). Mirrors upstream
    /// `MixtureParameters::set_mixture_parameters`.
    pub fn new(
        components: &[&'static FluidData],
        pairs: &[rustprop_core::fluid::MixBinaryPair],
        departures: &[rustprop_core::fluid::MixDepartureFn],
    ) -> Result<Self> {
        let reducing = Gerg2008Reducing::new(components, pairs)?;
        let n = components.len();
        let mut excess = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let (a, b) = (components[i].cas, components[j].cas);
                // Gerg2008Reducing::new already proved the pair exists.
                let pair = pairs
                    .iter()
                    .find(|p| (p.cas1 == a && p.cas2 == b) || (p.cas1 == b && p.cas2 == a))
                    .expect("pair vetted by Gerg2008Reducing::new");
                // Upstream: |F| < DBL_EPSILON gets the empty departure
                // function that just returns 0.
                let dep = if pair.f.abs() < f64::EPSILON {
                    crate::alpha::DepartureEval::zero()
                } else {
                    let name = pair.function.unwrap_or("");
                    let dep_fn = departures.iter().find(|d| d.name == name).ok_or_else(|| {
                        Error::Value(format!(
                            "Departure function name [{name}] seems to be invalid"
                        ))
                    })?;
                    crate::alpha::DepartureEval::new(dep_fn)
                };
                excess.push(ExcessPair {
                    i,
                    j,
                    f: pair.f,
                    dep,
                });
            }
        }
        Ok(MixtureModel {
            pure: components
                .iter()
                .map(|c| crate::alpha::HelmholtzEos::new(c))
                .collect(),
            crit_t: components.iter().map(|c| c.states.critical.t).collect(),
            crit_rhomolar: components
                .iter()
                .map(|c| c.states.critical.rhomolar)
                .collect(),
            crit_p: components.iter().map(|c| c.states.critical.p).collect(),
            triple_t: components.iter().map(|c| c.eos.sat_min_liquid.t).collect(),
            triple_p: components.iter().map(|c| c.eos.sat_min_liquid.p).collect(),
            r_component: components.iter().map(|c| c.eos.gas_constant).collect(),
            srk_consts: components
                .iter()
                .map(|c| (c.eos.reducing.t, c.eos.reducing.p, c.eos.acentric))
                .collect(),
            molar_mass: components.iter().map(|c| c.eos.molar_mass).collect(),
            reducing,
            excess,
        })
    }

    pub(crate) fn srk_component_consts(&self) -> &[(f64, f64, f64)] {
        &self.srk_consts
    }

    pub(crate) fn n_components(&self) -> usize {
        self.pure.len()
    }

    /// `iT_critical` / `iP_critical` / `irhomolar_critical` /
    /// `iacentric_factor` / `iT_triple` / `iP_triple` per component.
    pub(crate) fn crit_t(&self) -> &[f64] {
        &self.crit_t
    }
    pub(crate) fn crit_p(&self) -> &[f64] {
        &self.crit_p
    }
    pub(crate) fn crit_rhomolar(&self) -> &[f64] {
        &self.crit_rhomolar
    }
    pub(crate) fn acentric(&self, i: usize) -> f64 {
        self.srk_consts[i].2
    }
    pub(crate) fn triple_t(&self) -> &[f64] {
        &self.triple_t
    }
    pub(crate) fn triple_p(&self) -> &[f64] {
        &self.triple_p
    }

    /// Component i's pure alphar derivatives at the mixture (tau, delta).
    pub(crate) fn component_alphar_all(
        &self,
        i: usize,
        tau: f64,
        delta: f64,
    ) -> crate::alpha::HelmholtzDerivs {
        self.pure[i].alphar_all(tau, delta)
    }

    /// The excess pairs as `(i, j, F, departure)` views.
    pub(crate) fn excess_pairs(
        &self,
    ) -> impl Iterator<Item = (usize, usize, f64, &crate::alpha::DepartureEval)> {
        self.excess.iter().map(|p| (p.i, p.j, p.f, &p.dep))
    }

    pub(crate) fn molar_masses(&self) -> &[f64] {
        &self.molar_mass
    }

    /// `gas_constant()` — R for every mixture property relation.
    pub fn gas_constant(&self) -> f64 {
        R_U_CODATA
    }

    /// All residual derivatives at the mixture (tau, delta):
    /// corresponding-states sum plus excess term (upstream
    /// `ResidualHelmholtz::all` without the convenience products).
    pub fn alphar_all(&self, x: &[f64], tau: f64, delta: f64) -> crate::alpha::HelmholtzDerivs {
        let mut summer = crate::alpha::HelmholtzDerivs::default();
        for (i, eos) in self.pure.iter().enumerate() {
            let derivs = eos.alphar_all(tau, delta);
            summer.add_scaled(&derivs, x[i]);
        }
        for pair in &self.excess {
            let term = pair.dep.all(tau, delta);
            summer.add_scaled(&term, x[pair.i] * x[pair.j] * pair.f);
        }
        summer
    }

    /// The six ideal-gas derivatives through second order at the mixture
    /// (tau, delta) — Table B5, GERG 2008 (upstream
    /// `calc_all_alpha0_derivs_nocache`, mixture branch). Components are
    /// evaluated at their shifted `tau_i = T_ci tau / Tr`,
    /// `delta_i = delta rhor / rho_ci` and rescaled by `R_i / R_mix`; the
    /// `x_i ln x_i` entropy-of-mixing piece rides only the (0,0) derivative.
    /// (Upstream also calls `set_Tred(Tr)` here for GERG-2004 sinh/cosh
    /// alpha0 terms; no fluid document in the ported set carries them.)
    pub fn alpha0_all(
        &self,
        x: &[f64],
        tau: f64,
        delta: f64,
        tr: f64,
        rhor: f64,
    ) -> crate::alpha::HelmholtzDerivs {
        let r_mix = self.gas_constant();
        let mut ders = crate::alpha::HelmholtzDerivs::default();
        for (i, eos) in self.pure.iter().enumerate() {
            let rho_ci = self.crit_rhomolar[i];
            let t_ci = self.crit_t[i];
            let tau_i = t_ci * tau / tr;
            let delta_i = delta * rhor / rho_ci;
            let rratio = self.r_component[i] / r_mix;

            let pure = eos.alpha0_all(tau_i, delta_i);
            let logxi = if x[i].abs() > f64::EPSILON {
                x[i].ln()
            } else {
                0.0
            };
            ders.d00 += x[i] * rratio * (pure.d00 + logxi);
            ders.d10 += x[i] * rratio * rhor / rho_ci * pure.d10;
            ders.d01 += x[i] * rratio * t_ci / tr * pure.d01;
            ders.d20 += x[i] * rratio * (rhor / rho_ci).powi(2) * pure.d20;
            ders.d11 += x[i] * rratio * rhor / rho_ci * t_ci / tr * pure.d11;
            ders.d02 += x[i] * rratio * (t_ci / tr).powi(2) * pure.d02;
        }
        ders
    }
}
