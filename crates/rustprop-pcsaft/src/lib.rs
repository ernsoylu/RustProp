//! PC-SAFT EOS engine — operation-for-operation port of CoolProp 8
//! `src/Backends/PCSAFT/PCSAFTBackend.cpp` (hard-chain + dispersion +
//! Gross-Vrabec polar + Huang-Radosz association + ePC-SAFT ion terms).
//!
//! Slice 11b scope: the backend state, `calc_alphar`, `calc_dadt`,
//! `calc_compressibility_factor`, the residual caloric assemblies, and the
//! association helpers (`XA_find`, `dXAdt_find`, the damped-SS iterations
//! with upstream's per-kernel tolerances: 1e-15 in alphar/dadt, 1e-14 in Z).
//! Upstream's structural duplication (the ~150-line prep block appearing in
//! each kernel) collapses into one shared `Prep`, arithmetic order preserved.
//!
//! Deliberately reproduced upstream quirks (see the port survey):
//! - WATER ships `sigma = -1`; `calc_water_sigma(T)` replaces it (guards at
//!   t > 473.16 / t < 273 with messages quoting 473.15/273.15).
//! - `dielc_water` range-checks its ARGUMENT but evaluates the polynomials
//!   at the CURRENT `_T` member — callers passing a perturbed T get the
//!   dielectric constant at the stale live temperature.
//! - `dielc` is only assigned when water is present; the ion term of a
//!   water-free ionic system reads 0 here (Rust initializes; upstream reads
//!   an uninitialized double — UB we cannot and do not reproduce).
//! - `reduced_to_molar` ignores the ion T-independent-diameter rule.
//! - The association iteration has no convergence check (100 damped SS
//!   steps, then whatever the last iterate is).

#![allow(clippy::needless_range_loop)]
#![allow(clippy::too_many_arguments)]
// Upstream-verbatim expression shapes: `-1 * x` products, the PI literal
// (upstream defines its own), and the multi-value assoc-block return.
#![allow(clippy::neg_multiply)]
#![allow(clippy::approx_constant)]
#![allow(clippy::type_complexity)]
#![allow(clippy::assign_op_pattern)]

use rustprop_core::fluid::{PcsaftBinaryPair, PcsaftFluid};
use rustprop_core::{Error, Result};

pub use rustprop_core::UPSTREAM_VERSION;

// Constants (upstream PCSAFTBackend.h:16-20)
const KB: f64 = 1.380649e-23;
const PI: f64 = 3.141592653589793;
const N_AV: f64 = 6.02214076e23;
const E_CHRG: f64 = 1.6021766208e-19;
const PERM_VAC: f64 = 8.854187817e-22;

// Universal dispersion constants (duplicated verbatim in four upstream fns)
const A0: [f64; 7] = [
    0.9105631445,
    0.6361281449,
    2.6861347891,
    -26.547362491,
    97.759208784,
    -159.59154087,
    91.297774084,
];
const A1: [f64; 7] = [
    -0.3084016918,
    0.1860531159,
    -2.5030047259,
    21.419793629,
    -65.255885330,
    83.318680481,
    -33.746922930,
];
const A2: [f64; 7] = [
    -0.0906148351,
    0.4527842806,
    0.5962700728,
    -1.7241829131,
    -4.1302112531,
    13.776631870,
    -8.6728470368,
];
const B0: [f64; 7] = [
    0.7240946941,
    2.2382791861,
    -4.0025849485,
    -21.003576815,
    26.855641363,
    206.55133841,
    -355.60235612,
];
const B1: [f64; 7] = [
    -0.5755498075,
    0.6995095521,
    3.8925673390,
    -17.215471648,
    192.67226447,
    -161.82646165,
    -165.20769346,
];
const B2: [f64; 7] = [
    0.0976883116,
    -0.2557574982,
    -9.1558561530,
    20.642075974,
    -38.804430052,
    93.626774077,
    -29.666905585,
];

// Gross-Vrabec dipole tables
const A0DIP: [f64; 5] = [0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308];
const A1DIP: [f64; 5] = [0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135];
const A2DIP: [f64; 5] = [-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575];
const B0DIP: [f64; 5] = [0.2187939, -1.1896431, 1.1626889, 0.0, 0.0];
const B1DIP: [f64; 5] = [-0.5873164, 1.2489132, -0.5085280, 0.0, 0.0];
const B2DIP: [f64; 5] = [3.4869576, -14.915974, 15.372022, 0.0, 0.0];
const C0DIP: [f64; 5] = [-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0.0];
const C1DIP: [f64; 5] = [-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0.0];
const C2DIP: [f64; 5] = [-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0.0];
/// Conversion factor, see the note below Table 2 in Gross and Vrabec 2006.
const CONV: f64 = 7242.702976750923;

/// A live (mutable) component: WATER's sigma is assigned at runtime.
#[derive(Clone)]
struct Component {
    cas: &'static str,
    m: f64,
    sigma: f64,
    u: f64,
    u_ab: f64,
    vol_a: f64,
    assoc_scheme: &'static [&'static str],
    dipm: f64,
    dipnum: f64,
    z: f64,
    molemass: f64,
}

/// The PC-SAFT backend state (upstream `PCSAFTBackend` minus the flash
/// scaffolding, which lands with slice 11c).
pub struct PcsaftBackend {
    n: usize,
    components: Vec<Component>,
    pub mole_fractions: Vec<f64>,
    /// Flat N*N; EMPTY for pure fluids (the `k_ij.empty()` branches test it).
    k_ij: Vec<f64>,
    k_ij_t: Vec<f64>,
    ion_term: bool,
    polar_term: bool,
    assoc_term: bool,
    water_present: bool,
    water_idx: usize,
    /// Sites per component (from the association schemes).
    assoc_num: Vec<i32>,
    /// Flat num_sites x num_sites site-compatibility matrix.
    assoc_matrix: Vec<i32>,
    /// Only ever assigned when water is present (upstream reads an
    /// uninitialized double otherwise; Rust zero-initializes).
    dielc: f64,
    pub t: f64,
    pub rhomolar: f64,
}

/// `get_scheme_index` site charges (DataStructures.cpp:438-493, case-sensitive).
fn scheme_charges(scheme: &str) -> Result<&'static [i32]> {
    Ok(match scheme {
        "1" => &[0],
        "2A" => &[0, 0],
        "2B" => &[-1, 1],
        "3A" => &[0, 0, 0],
        "3B" => &[-1, -1, 1],
        "4A" => &[0, 0, 0, 0],
        "4B" => &[1, 1, 1, -1],
        "4C" => &[-1, -1, 1, 1],
        other => {
            return Err(Error::Value(format!(
                "{other} is not a valid association type."
            )));
        }
    })
}

impl PcsaftBackend {
    /// Build from resolved fluids and the kij table (upstream ctor; every
    /// i != j pair must exist in the table or construction fails).
    pub fn new(fluids: &[&'static PcsaftFluid], pairs: &[PcsaftBinaryPair]) -> Result<Self> {
        let n = fluids.len();
        let mut ion_term = false;
        let mut polar_term = false;
        let mut assoc_term = false;
        let mut water_present = false;
        let mut water_idx = 0;
        let components: Vec<Component> = fluids
            .iter()
            .map(|f| Component {
                cas: f.cas,
                m: f.m,
                sigma: f.sigma,
                u: f.u,
                u_ab: f.u_ab,
                vol_a: f.vol_a,
                assoc_scheme: f.assoc_scheme,
                dipm: f.dipm,
                dipnum: f.dipnum,
                z: f.z,
                molemass: f.molemass,
            })
            .collect();
        for (i, c) in components.iter().enumerate() {
            if c.z != 0.0 {
                ion_term = true;
            }
            if c.dipm != 0.0 {
                polar_term = true;
            }
            if c.vol_a != 0.0 {
                assoc_term = true;
            }
            if c.cas == "7732-18-5" {
                water_present = true;
                water_idx = i;
            }
        }

        // Association scheme -> per-component site counts + site matrix
        let mut assoc_num = Vec::new();
        let mut assoc_matrix = Vec::new();
        if assoc_term {
            let mut charge: Vec<i32> = Vec::new();
            for c in &components {
                let mut num_sites = 0;
                for scheme in c.assoc_scheme {
                    let ch = scheme_charges(scheme)?;
                    charge.extend_from_slice(ch);
                    num_sites += ch.len() as i32;
                }
                assoc_num.push(num_sites);
            }
            for &c1 in &charge {
                for &c2 in &charge {
                    let ok = c1 == 0 || c2 == 0 || (c1 == 1 && c2 == -1) || (c1 == -1 && c2 == 1);
                    assoc_matrix.push(i32::from(ok));
                }
            }
        }

        // kij load (pure fluids keep the vectors EMPTY)
        let mut k_ij = Vec::new();
        let mut k_ij_t = Vec::new();
        if n > 1 {
            k_ij = vec![0.0; n * n];
            k_ij_t = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    if i != j {
                        let (a, b) = (components[i].cas, components[j].cas);
                        let pair = pairs
                            .iter()
                            .find(|p| {
                                (p.cas1 == a && p.cas2 == b) || (p.cas1 == b && p.cas2 == a)
                            })
                            .ok_or_else(|| {
                                let mut sorted = [a, b];
                                sorted.sort_unstable();
                                Error::Value(format!(
                                    "Could not match the binary pair [{},{}] - for now this is an error.",
                                    sorted[0], sorted[1]
                                ))
                            })?;
                        k_ij[i * n + j] = pair.kij;
                        k_ij_t[i * n + j] = pair.kij_t;
                    }
                }
            }
        }

        let mole_fractions = if n == 1 { vec![1.0] } else { Vec::new() };
        Ok(PcsaftBackend {
            n,
            components,
            mole_fractions,
            k_ij,
            k_ij_t,
            ion_term,
            polar_term,
            assoc_term,
            water_present,
            water_idx,
            assoc_num,
            assoc_matrix,
            dielc: 0.0,
            t: f64::NAN,
            rhomolar: f64::NAN,
        })
    }

    pub fn n_components(&self) -> usize {
        self.n
    }

    pub fn set_mole_fractions(&mut self, x: &[f64]) {
        self.mole_fractions = x.to_vec();
    }

    /// `calc_molar_mass`: x-weighted component masses.
    pub fn molar_mass(&self) -> f64 {
        self.components
            .iter()
            .zip(&self.mole_fractions)
            .map(|(c, x)| x * c.molemass)
            .sum()
    }

    /// `PCSAFTFluid::calc_water_sigma` — replaces WATER's -1 sigma sentinel.
    /// Guard bounds (473.16, 273) deliberately mismatch the message text.
    pub fn calc_water_sigma(&mut self, t: f64) -> Result<()> {
        if t > 473.16 {
            return Err(Error::Value(
                "The current function for sigma for water is only valid for temperatures below 473.15 K.".into(),
            ));
        } else if t < 273.0 {
            return Err(Error::Value(
                "The current function for sigma for water is only valid for temperatures above 273.15 K.".into(),
            ));
        }
        self.components[self.water_idx].sigma =
            3.8395 + 1.2828 * (-0.0074944 * t).exp() - 1.3939 * (-0.00056029 * t).exp();
        Ok(())
    }

    /// `dielc_water(t)` — range checks use `t`, the POLYNOMIALS use the live
    /// `_T` member (upstream bug, reproduced).
    pub fn dielc_water(&self, t: f64) -> Result<f64> {
        if t < 263.15 {
            Err(Error::Value(
                "The current function for the dielectric constant for water is only valid for temperatures above 263.15 K.".into(),
            ))
        } else if t <= 368.15 {
            Ok(7.6555618295E-04 * self.t * self.t - 8.1783881423E-01 * self.t + 2.5419616803E+02)
        } else if t <= 443.15 {
            Ok(0.0005003272124 * self.t * self.t - 0.6285556029 * self.t + 220.4467027)
        } else {
            Err(Error::Value(
                "The current function for the dielectric constant for water is only valid for temperatures less than 443.15 K.".into(),
            ))
        }
    }

    /// Set the live (T, rhomolar) state, refreshing water sigma/dielc as the
    /// `DmolarT_INPUTS` update path does.
    pub fn set_state_dmolar_t(&mut self, rhomolar: f64, t: f64) -> Result<()> {
        self.t = t;
        self.rhomolar = rhomolar;
        if self.water_present {
            self.calc_water_sigma(t)?;
            let d = self.dielc_water(t)?;
            self.dielc = d;
        }
        Ok(())
    }

    /// `reduced_to_molar(nu, T)` — always uses the T-dependent diameter,
    /// even for ions (upstream quirk 5).
    pub fn reduced_to_molar(&self, nu: f64, t: f64) -> f64 {
        let mut summ = 0.0;
        for i in 0..self.n {
            let d =
                self.components[i].sigma * (1.0 - 0.12 * (-3.0 * self.components[i].u / t).exp());
            summ += self.mole_fractions[i] * self.components[i].m * d.powi(3);
        }
        6.0 / PI * nu / summ * 1.0e30 / N_AV
    }

    // -- shared prep (the block upstream duplicates in each kernel) --------

    fn prep(&self) -> Prep {
        let ncomp = self.n;
        let t = self.t;
        let x = &self.mole_fractions;
        let mut d = vec![0.0; ncomp];
        let mut dd_dt = vec![0.0; ncomp];
        for i in 0..ncomp {
            let c = &self.components[i];
            d[i] = c.sigma * (1.0 - 0.12 * (-3.0 * c.u / t).exp());
            dd_dt[i] = c.sigma * -3.0 * c.u / t / t * 0.12 * (-3.0 * c.u / t).exp();
        }
        if self.ion_term {
            for i in 0..ncomp {
                if self.components[i].z != 0.0 {
                    d[i] = self.components[i].sigma * (1.0 - 0.12);
                    dd_dt[i] = 0.0;
                }
            }
        }

        let den = self.rhomolar * N_AV / 1.0e30;

        let mut zeta = [0.0; 4];
        for i in 0..4 {
            let mut summ = 0.0;
            for j in 0..ncomp {
                summ += x[j] * self.components[j].m * d[j].powi(i as i32);
            }
            zeta[i] = PI / 6.0 * den * summ;
        }
        let mut dzeta_dt = [0.0; 4];
        for i in 1..4 {
            let mut summ = 0.0;
            for j in 0..ncomp {
                summ +=
                    x[j] * self.components[j].m * (i as f64) * dd_dt[j] * d[j].powi(i as i32 - 1);
            }
            dzeta_dt[i] = PI / 6.0 * den * summ;
        }

        let eta = zeta[3];
        let mut m_avg = 0.0;
        for i in 0..ncomp {
            m_avg += x[i] * self.components[i].m;
        }

        let mut ghs = vec![0.0; ncomp * ncomp];
        let mut dghs_dt = vec![0.0; ncomp * ncomp];
        let mut denghs = vec![0.0; ncomp * ncomp];
        let mut e_ij = vec![0.0; ncomp * ncomp];
        let mut s_ij = vec![0.0; ncomp * ncomp];
        let mut m2es3 = 0.0;
        let mut m2e2s3 = 0.0;
        let mut idx = 0usize;
        for i in 0..ncomp {
            for j in 0..ncomp {
                let ci = &self.components[i];
                let cj = &self.components[j];
                s_ij[idx] = (ci.sigma + cj.sigma) / 2.0;
                let dispersion_allowed = if self.ion_term {
                    // like-charge pairs keep e_ij at zero
                    ci.z * cj.z <= 0.0
                } else {
                    true
                };
                if dispersion_allowed {
                    if self.k_ij.is_empty() {
                        e_ij[idx] = (ci.u * cj.u).sqrt();
                    } else {
                        e_ij[idx] =
                            (ci.u * cj.u).sqrt() * (1.0 - (self.k_ij[idx] + self.k_ij_t[idx] * t));
                    }
                }
                m2es3 += x[i] * x[j] * ci.m * cj.m * e_ij[idx] / t * s_ij[idx].powi(3);
                m2e2s3 += x[i] * x[j] * ci.m * cj.m * (e_ij[idx] / t).powi(2) * s_ij[idx].powi(3);
                ghs[idx] = 1.0 / (1.0 - zeta[3])
                    + (d[i] * d[j] / (d[i] + d[j])) * 3.0 * zeta[2]
                        / (1.0 - zeta[3])
                        / (1.0 - zeta[3])
                    + (d[i] * d[j] / (d[i] + d[j])).powi(2) * 2.0 * zeta[2] * zeta[2]
                        / (1.0 - zeta[3]).powi(3);
                let ddij_dt = (d[i] * d[j] / (d[i] + d[j]))
                    * (dd_dt[i] / d[i] + dd_dt[j] / d[j] - (dd_dt[i] + dd_dt[j]) / (d[i] + d[j]));
                dghs_dt[idx] = dzeta_dt[3] / (1.0 - zeta[3]).powi(2)
                    + 3.0 * (ddij_dt * zeta[2] + (d[i] * d[j] / (d[i] + d[j])) * dzeta_dt[2])
                        / (1.0 - zeta[3]).powi(2)
                    + 4.0
                        * (d[i] * d[j] / (d[i] + d[j]))
                        * zeta[2]
                        * (1.5 * dzeta_dt[3]
                            + ddij_dt * zeta[2]
                            + (d[i] * d[j] / (d[i] + d[j])) * dzeta_dt[2])
                        / (1.0 - zeta[3]).powi(3)
                    + 6.0 * ((d[i] * d[j] / (d[i] + d[j])) * zeta[2]).powi(2) * dzeta_dt[3]
                        / (1.0 - zeta[3]).powi(4);
                denghs[idx] = zeta[3] / (1.0 - zeta[3]) / (1.0 - zeta[3])
                    + (d[i] * d[j] / (d[i] + d[j]))
                        * (3.0 * zeta[2] / (1.0 - zeta[3]) / (1.0 - zeta[3])
                            + 6.0 * zeta[2] * zeta[3] / (1.0 - zeta[3]).powi(3))
                    + (d[i] * d[j] / (d[i] + d[j])).powi(2)
                        * (4.0 * zeta[2] * zeta[2] / (1.0 - zeta[3]).powi(3)
                            + 6.0 * zeta[2] * zeta[2] * zeta[3] / (1.0 - zeta[3]).powi(4));
                idx += 1;
            }
        }

        // Dispersion a/b coefficients + I1/I2 + C1/C2
        let mut a = [0.0; 7];
        let mut b = [0.0; 7];
        for i in 0..7 {
            a[i] = A0[i]
                + (m_avg - 1.0) / m_avg * A1[i]
                + (m_avg - 1.0) / m_avg * (m_avg - 2.0) / m_avg * A2[i];
            b[i] = B0[i]
                + (m_avg - 1.0) / m_avg * B1[i]
                + (m_avg - 1.0) / m_avg * (m_avg - 2.0) / m_avg * B2[i];
        }
        let c1 = 1.0
            / (1.0
                + m_avg * (8.0 * eta - 2.0 * eta * eta) / (1.0 - eta).powi(4)
                + (1.0 - m_avg)
                    * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta.powi(3) - 2.0 * eta.powi(4))
                    / ((1.0 - eta) * (2.0 - eta)).powf(2.0));
        let c2 = -1.0
            * c1
            * c1
            * (m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / (1.0 - eta).powi(5)
                + (1.0 - m_avg) * (2.0 * eta.powi(3) + 12.0 * eta * eta - 48.0 * eta + 40.0)
                    / ((1.0 - eta) * (2.0 - eta)).powf(3.0));

        Prep {
            d,
            dd_dt,
            den,
            zeta,
            dzeta_dt,
            eta,
            m_avg,
            ghs,
            dghs_dt,
            denghs,
            e_ij,
            s_ij,
            m2es3,
            m2e2s3,
            a,
            b,
            c1,
            c2,
        }
    }

    fn dipm_sq(&self) -> Vec<f64> {
        (0..self.n)
            .map(|i| {
                let c = &self.components[i];
                c.dipm.powf(2.0) / (c.m * c.u * c.sigma.powf(3.0)) * CONV
            })
            .collect()
    }

    /// Association-site bookkeeping: (num_sites, iA, x_assoc).
    fn assoc_sites(&self) -> (usize, Vec<usize>, Vec<f64>) {
        let mut num_sites = 0usize;
        let mut i_a = Vec::new();
        for (comp, &cnt) in self.assoc_num.iter().enumerate() {
            num_sites += cnt as usize;
            for _ in 0..cnt {
                i_a.push(comp);
            }
        }
        let x_assoc: Vec<f64> = i_a.iter().map(|&i| self.mole_fractions[i]).collect();
        (num_sites, i_a, x_assoc)
    }

    /// `XA_find`: one successive-substitution step.
    fn xa_find(&self, xa_guess: &[f64], delta_ij: &[f64], den: f64, x: &[f64]) -> Vec<f64> {
        let num_sites = xa_guess.len();
        let mut xa = xa_guess.to_vec();
        let mut idxij = 0usize;
        for i in 0..num_sites {
            let mut summ = 0.0;
            for j in 0..num_sites {
                summ += den * x[j] * xa_guess[j] * delta_ij[idxij];
                idxij += 1;
            }
            xa[i] = 1.0 / (1.0 + summ);
        }
        xa
    }

    /// The damped-SS XA iteration (no convergence check after the cap).
    fn iterate_xa(&self, delta_ij: &[f64], den: f64, x_assoc: &[f64], tol: f64, xa: &mut [f64]) {
        let num_sites = xa.len();
        let mut ctr = 0;
        let mut dif = 1000.0;
        let mut xa_old = xa.to_vec();
        while ctr < 100 && dif > tol {
            ctr += 1;
            let xa_new = self.xa_find(&xa_old, delta_ij, den, x_assoc);
            dif = 0.0;
            for i in 0..num_sites {
                dif += (xa_new[i] - xa_old[i]).abs();
            }
            for i in 0..num_sites {
                xa[i] = xa_new[i];
                xa_old[i] = (xa_new[i] + xa_old[i]) / 2.0;
            }
        }
    }

    /// `dXAdt_find`: dense LU solve for dXA/dT.
    fn dxadt_find(
        &self,
        delta_ij: &[f64],
        den: f64,
        xa: &[f64],
        ddelta_dt: &[f64],
        x: &[f64],
    ) -> Vec<f64> {
        let num_sites = xa.len();
        let mut b = vec![0.0; num_sites];
        let mut a = vec![vec![0.0; num_sites]; num_sites];
        let mut ij = 0usize;
        for i in 0..num_sites {
            let mut summ = 0.0;
            for j in 0..num_sites {
                b[i] -= x[j] * xa[j] * ddelta_dt[ij];
                a[i][j] = x[j] * delta_ij[ij];
                summ += x[j] * xa[j] * delta_ij[ij];
                ij += 1;
            }
            a[i][i] = (1.0 + den * summ).powf(2.0) / den;
        }
        solve_dense(&mut a, &b)
    }

    /// XA at the live state for one kernel (shared Δ construction).
    fn assoc_block(
        &self,
        prep: &Prep,
        with_dt: bool,
        tol: f64,
    ) -> (usize, Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let ncomp = self.n;
        let t = self.t;
        let (num_sites, i_a, x_assoc) = self.assoc_sites();
        let mut xa = vec![0.0; num_sites];
        let mut delta_ij = vec![0.0; num_sites * num_sites];
        let mut ddelta_dt = vec![0.0; num_sites * num_sites];
        let mut idxa = 0usize;
        for i in 0..num_sites {
            let idxi = i_a[i] * ncomp + i_a[i];
            for j in 0..num_sites {
                let idxj = i_a[j] * ncomp + i_a[j];
                if self.assoc_matrix[idxa] != 0 {
                    let e_abij =
                        (self.components[i_a[i]].u_ab + self.components[i_a[j]].u_ab) / 2.0;
                    let vol_abij = (self.components[i_a[i]].vol_a * self.components[i_a[j]].vol_a)
                        .sqrt()
                        * ((prep.s_ij[idxi] * prep.s_ij[idxj]).sqrt()
                            / (0.5 * (prep.s_ij[idxi] + prep.s_ij[idxj])))
                            .powi(3);
                    delta_ij[idxa] = prep.ghs[i_a[i] * ncomp + i_a[j]]
                        * ((e_abij / t).exp() - 1.0)
                        * prep.s_ij[i_a[i] * ncomp + i_a[j]].powi(3)
                        * vol_abij;
                    if with_dt {
                        ddelta_dt[idxa] = prep.s_ij[idxj].powi(3)
                            * vol_abij
                            * (-e_abij / t.powi(2)
                                * (e_abij / t).exp()
                                * prep.ghs[i_a[i] * ncomp + i_a[j]]
                                + prep.dghs_dt[i_a[i] * ncomp + i_a[j]]
                                    * ((e_abij / t).exp() - 1.0));
                    }
                }
                idxa += 1;
            }
            xa[i] = (-1.0 + (1.0 + 8.0 * prep.den * delta_ij[i * num_sites + i]).sqrt())
                / (4.0 * prep.den * delta_ij[i * num_sites + i]);
            if !xa[i].is_finite() {
                xa[i] = 0.02;
            }
        }
        self.iterate_xa(&delta_ij, prep.den, &x_assoc, tol, &mut xa);
        (num_sites, i_a, x_assoc, xa, delta_ij, ddelta_dt)
    }

    // -- calc_alphar -------------------------------------------------------

    /// `calc_alphar()` at the live (T, rhomolar, x).
    pub fn calc_alphar(&self) -> f64 {
        let ncomp = self.n;
        let t = self.t;
        let x = &self.mole_fractions;
        let p = self.prep();

        let ares_hs = 1.0 / p.zeta[0]
            * (3.0 * p.zeta[1] * p.zeta[2] / (1.0 - p.zeta[3])
                + p.zeta[2].powi(3) / (p.zeta[3] * (1.0 - p.zeta[3]).powi(2))
                + (p.zeta[2].powi(3) / p.zeta[3].powi(2) - p.zeta[0]) * (1.0 - p.zeta[3]).ln());
        let mut summ = 0.0;
        for i in 0..ncomp {
            summ += x[i] * (self.components[i].m - 1.0) * p.ghs[i * ncomp + i].ln();
        }
        let ares_hc = p.m_avg * ares_hs - summ;

        let mut i1 = 0.0;
        let mut i2 = 0.0;
        for i in 0..7 {
            i1 += p.a[i] * p.eta.powi(i as i32);
            i2 += p.b[i] * p.eta.powi(i as i32);
        }
        let ares_disp =
            -2.0 * PI * p.den * i1 * p.m2es3 - PI * p.den * p.m_avg * p.c1 * i2 * p.m2e2s3;

        // Polar
        let mut ares_polar = 0.0;
        if self.polar_term {
            let dipm_sq = self.dipm_sq();
            let mut a2t = 0.0;
            let mut a3t = 0.0;
            for i in 0..ncomp {
                for j in 0..ncomp {
                    let mut m_ij = (self.components[i].m * self.components[j].m).sqrt();
                    if m_ij > 2.0 {
                        m_ij = 2.0;
                    }
                    let mut j2 = 0.0;
                    for l in 0..5 {
                        let adip = A0DIP[l]
                            + (m_ij - 1.0) / m_ij * A1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * A2DIP[l];
                        let bdip = B0DIP[l]
                            + (m_ij - 1.0) / m_ij * B1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * B2DIP[l];
                        j2 += (adip + bdip * p.e_ij[i * ncomp + j] / t) * p.eta.powi(l as i32);
                    }
                    a2t += x[i] * x[j] * p.e_ij[i * ncomp + i] / t * p.e_ij[j * ncomp + j] / t
                        * p.s_ij[i * ncomp + i].powi(3)
                        * p.s_ij[j * ncomp + j].powi(3)
                        / p.s_ij[i * ncomp + j].powi(3)
                        * self.components[i].dipnum
                        * self.components[j].dipnum
                        * dipm_sq[i]
                        * dipm_sq[j]
                        * j2;

                    for k in 0..ncomp {
                        let mut m_ijk =
                            (self.components[i].m * self.components[j].m * self.components[k].m)
                                .powf(1.0 / 3.0);
                        if m_ijk > 2.0 {
                            m_ijk = 2.0;
                        }
                        let mut j3 = 0.0;
                        for l in 0..5 {
                            let cdip = C0DIP[l]
                                + (m_ijk - 1.0) / m_ijk * C1DIP[l]
                                + (m_ijk - 1.0) / m_ijk * (m_ijk - 2.0) / m_ijk * C2DIP[l];
                            j3 += cdip * p.eta.powi(l as i32);
                        }
                        a3t += x[i] * x[j] * x[k] * p.e_ij[i * ncomp + i] / t
                            * p.e_ij[j * ncomp + j]
                            / t
                            * p.e_ij[k * ncomp + k]
                            / t
                            * p.s_ij[i * ncomp + i].powi(3)
                            * p.s_ij[j * ncomp + j].powi(3)
                            * p.s_ij[k * ncomp + k].powi(3)
                            / p.s_ij[i * ncomp + j]
                            / p.s_ij[i * ncomp + k]
                            / p.s_ij[j * ncomp + k]
                            * self.components[i].dipnum
                            * self.components[j].dipnum
                            * self.components[k].dipnum
                            * dipm_sq[i]
                            * dipm_sq[j]
                            * dipm_sq[k]
                            * j3;
                    }
                }
            }
            let a2 = -PI * p.den * a2t;
            let a3 = -4.0 / 3.0 * PI * PI * p.den * p.den * a3t;
            if a2 != 0.0 {
                ares_polar = a2 / (1.0 - a3 / a2);
            }
        }

        // Association
        let mut ares_assoc = 0.0;
        if self.assoc_term {
            let (num_sites, i_a, _x_assoc, xa, _delta, _ddt) = self.assoc_block(&p, false, 1e-15);
            for i in 0..num_sites {
                ares_assoc += x[i_a[i]] * (xa[i].ln() - 0.5 * xa[i] + 0.5);
            }
        }

        // Ion
        let mut ares_ion = 0.0;
        if self.ion_term {
            let q: Vec<f64> = self.components.iter().map(|c| c.z * E_CHRG).collect();
            let mut summ = 0.0;
            for i in 0..ncomp {
                summ += self.components[i].z * self.components[i].z * x[i];
            }
            let kappa = (p.den * E_CHRG * E_CHRG / KB / t / (self.dielc * PERM_VAC) * summ).sqrt();
            if kappa != 0.0 {
                let mut summ2 = 0.0;
                for i in 0..ncomp {
                    let chi = 3.0 / (kappa * p.d[i]).powi(3)
                        * (1.5 + (1.0 + kappa * p.d[i]).ln() - 2.0 * (1.0 + kappa * p.d[i])
                            + 0.5 * (1.0 + kappa * p.d[i]).powi(2));
                    summ2 += x[i] * q[i] * q[i] * chi * kappa;
                }
                ares_ion = -1.0 / 12.0 / PI / KB / t / (self.dielc * PERM_VAC) * summ2;
            }
        }

        ares_hc + ares_disp + ares_polar + ares_assoc + ares_ion
    }

    // -- calc_dadt ---------------------------------------------------------

    /// `calc_dadt()`: d(alphar)/dT at constant rho, x.
    pub fn calc_dadt(&self) -> f64 {
        let ncomp = self.n;
        let t = self.t;
        let x = &self.mole_fractions;
        let p = self.prep();

        let dadt_hs = 1.0 / p.zeta[0]
            * (3.0 * (p.dzeta_dt[1] * p.zeta[2] + p.zeta[1] * p.dzeta_dt[2]) / (1.0 - p.zeta[3])
                + 3.0 * p.zeta[1] * p.zeta[2] * p.dzeta_dt[3] / (1.0 - p.zeta[3]).powi(2)
                + 3.0 * p.zeta[2].powi(2) * p.dzeta_dt[2] / p.zeta[3] / (1.0 - p.zeta[3]).powi(2)
                + p.zeta[2].powi(3) * p.dzeta_dt[3] * (3.0 * p.zeta[3] - 1.0)
                    / p.zeta[3].powi(2)
                    / (1.0 - p.zeta[3]).powi(3)
                + (3.0 * p.zeta[2].powi(2) * p.dzeta_dt[2] * p.zeta[3]
                    - 2.0 * p.zeta[2].powi(3) * p.dzeta_dt[3])
                    / p.zeta[3].powi(3)
                    * (1.0 - p.zeta[3]).ln()
                + (p.zeta[0] - p.zeta[2].powi(3) / p.zeta[3].powi(2)) * p.dzeta_dt[3]
                    / (1.0 - p.zeta[3]));

        let mut i1 = 0.0;
        let mut i2 = 0.0;
        let mut di1_dt = 0.0;
        let mut di2_dt = 0.0;
        for i in 0..7 {
            i1 += p.a[i] * p.eta.powi(i as i32);
            i2 += p.b[i] * p.eta.powi(i as i32);
            di1_dt += p.a[i] * p.dzeta_dt[3] * (i as f64) * p.eta.powi(i as i32 - 1);
            di2_dt += p.b[i] * p.dzeta_dt[3] * (i as f64) * p.eta.powi(i as i32 - 1);
        }
        let dc1_dt = p.c2 * p.dzeta_dt[3];

        let mut summ = 0.0;
        for i in 0..ncomp {
            summ += x[i] * (self.components[i].m - 1.0) * p.dghs_dt[i * ncomp + i]
                / p.ghs[i * ncomp + i];
        }
        let dadt_hc = p.m_avg * dadt_hs - summ;
        let dadt_disp = -2.0 * PI * p.den * (di1_dt - i1 / t) * p.m2es3
            - PI * p.den * p.m_avg * (dc1_dt * i2 + p.c1 * di2_dt - 2.0 * p.c1 * i2 / t) * p.m2e2s3;

        // Polar
        let mut dadt_polar = 0.0;
        if self.polar_term {
            let dipm_sq = self.dipm_sq();
            let mut a2t = 0.0;
            let mut a3t = 0.0;
            let mut da2_dt_t = 0.0;
            let mut da3_dt_t = 0.0;
            for i in 0..ncomp {
                for j in 0..ncomp {
                    let mut m_ij = (self.components[i].m * self.components[j].m).sqrt();
                    if m_ij > 2.0 {
                        m_ij = 2.0;
                    }
                    let mut j2 = 0.0;
                    let mut dj2_dt = 0.0;
                    for l in 0..5 {
                        let adip = A0DIP[l]
                            + (m_ij - 1.0) / m_ij * A1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * A2DIP[l];
                        let bdip = B0DIP[l]
                            + (m_ij - 1.0) / m_ij * B1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * B2DIP[l];
                        j2 += (adip + bdip * p.e_ij[i * ncomp + j] / t) * p.eta.powi(l as i32);
                        // NOTE upstream uses e_ij[j*ncomp+j] here (the jj
                        // diagonal), unlike J2's i*ncomp+j — reproduced.
                        dj2_dt += adip * (l as f64) * p.eta.powi(l as i32 - 1) * p.dzeta_dt[3]
                            + bdip
                                * p.e_ij[j * ncomp + j]
                                * (1.0 / t * (l as f64) * p.eta.powi(l as i32 - 1) * p.dzeta_dt[3]
                                    - 1.0 / t.powf(2.0) * p.eta.powi(l as i32));
                    }
                    a2t += x[i] * x[j] * p.e_ij[i * ncomp + i] / t * p.e_ij[j * ncomp + j] / t
                        * p.s_ij[i * ncomp + i].powi(3)
                        * p.s_ij[j * ncomp + j].powi(3)
                        / p.s_ij[i * ncomp + j].powi(3)
                        * self.components[i].dipnum
                        * self.components[j].dipnum
                        * dipm_sq[i]
                        * dipm_sq[j]
                        * j2;
                    da2_dt_t += x[i]
                        * x[j]
                        * p.e_ij[i * ncomp + i]
                        * p.e_ij[j * ncomp + j]
                        * p.s_ij[i * ncomp + i].powi(3)
                        * p.s_ij[j * ncomp + j].powi(3)
                        / p.s_ij[i * ncomp + j].powi(3)
                        * self.components[i].dipnum
                        * self.components[j].dipnum
                        * dipm_sq[i]
                        * dipm_sq[j]
                        * (dj2_dt / t.powi(2) - 2.0 * j2 / t.powi(3));

                    for k in 0..ncomp {
                        let mut m_ijk =
                            (self.components[i].m * self.components[j].m * self.components[k].m)
                                .powf(1.0 / 3.0);
                        if m_ijk > 2.0 {
                            m_ijk = 2.0;
                        }
                        let mut j3 = 0.0;
                        let mut dj3_dt = 0.0;
                        for l in 0..5 {
                            let cdip = C0DIP[l]
                                + (m_ijk - 1.0) / m_ijk * C1DIP[l]
                                + (m_ijk - 1.0) / m_ijk * (m_ijk - 2.0) / m_ijk * C2DIP[l];
                            j3 += cdip * p.eta.powi(l as i32);
                            dj3_dt += cdip * (l as f64) * p.eta.powi(l as i32 - 1) * p.dzeta_dt[3];
                        }
                        a3t += x[i] * x[j] * x[k] * p.e_ij[i * ncomp + i] / t
                            * p.e_ij[j * ncomp + j]
                            / t
                            * p.e_ij[k * ncomp + k]
                            / t
                            * p.s_ij[i * ncomp + i].powi(3)
                            * p.s_ij[j * ncomp + j].powi(3)
                            * p.s_ij[k * ncomp + k].powi(3)
                            / p.s_ij[i * ncomp + j]
                            / p.s_ij[i * ncomp + k]
                            / p.s_ij[j * ncomp + k]
                            * self.components[i].dipnum
                            * self.components[j].dipnum
                            * self.components[k].dipnum
                            * dipm_sq[i]
                            * dipm_sq[j]
                            * dipm_sq[k]
                            * j3;
                        da3_dt_t += x[i]
                            * x[j]
                            * x[k]
                            * p.e_ij[i * ncomp + i]
                            * p.e_ij[j * ncomp + j]
                            * p.e_ij[k * ncomp + k]
                            * p.s_ij[i * ncomp + i].powi(3)
                            * p.s_ij[j * ncomp + j].powi(3)
                            * p.s_ij[k * ncomp + k].powi(3)
                            / p.s_ij[i * ncomp + j]
                            / p.s_ij[i * ncomp + k]
                            / p.s_ij[j * ncomp + k]
                            * self.components[i].dipnum
                            * self.components[j].dipnum
                            * self.components[k].dipnum
                            * dipm_sq[i]
                            * dipm_sq[j]
                            * dipm_sq[k]
                            * (-3.0 * j3 / t.powi(4) + dj3_dt / t.powi(3));
                    }
                }
            }
            let a2 = -PI * p.den * a2t;
            let a3 = -4.0 / 3.0 * PI * PI * p.den * p.den * a3t;
            let da2_dt = -PI * p.den * da2_dt_t;
            let da3_dt = -4.0 / 3.0 * PI * PI * p.den * p.den * da3_dt_t;
            if a2 != 0.0 {
                dadt_polar = (da2_dt - 2.0 * a3 / a2 * da2_dt + da3_dt) / (1.0 - a3 / a2).powi(2);
            }
        }

        // Association
        let mut dadt_assoc = 0.0;
        if self.assoc_term {
            let (num_sites, i_a, x_assoc, xa, delta_ij, ddelta_dt) =
                self.assoc_block(&p, true, 1e-15);
            let dxa_dt = self.dxadt_find(&delta_ij, p.den, &xa, &ddelta_dt, &x_assoc);
            for i in 0..num_sites {
                dadt_assoc += x[i_a[i]] * (1.0 / xa[i] - 0.5) * dxa_dt[i];
            }
        }

        // Ion
        let mut dadt_ion = 0.0;
        if self.ion_term {
            let q: Vec<f64> = self.components.iter().map(|c| c.z * E_CHRG).collect();
            let mut summ = 0.0;
            for i in 0..ncomp {
                summ += self.components[i].z * self.components[i].z * x[i];
            }
            let kappa = (p.den * E_CHRG * E_CHRG / KB / t / (self.dielc * PERM_VAC) * summ).sqrt();
            if kappa != 0.0 {
                let mut chi = vec![0.0; ncomp];
                let mut dchikap_dk = vec![0.0; ncomp];
                let mut summ2 = 0.0;
                for i in 0..ncomp {
                    chi[i] = 3.0 / (kappa * p.d[i]).powi(3)
                        * (1.5 + (1.0 + kappa * p.d[i]).ln() - 2.0 * (1.0 + kappa * p.d[i])
                            + 0.5 * (1.0 + kappa * p.d[i]).powi(2));
                    dchikap_dk[i] = -2.0 * chi[i] + 3.0 / (1.0 + kappa * p.d[i]);
                    summ2 += x[i] * self.components[i].z * self.components[i].z;
                }
                let dkappa_dt =
                    -0.5 * p.den * E_CHRG * E_CHRG / KB / t / t / (self.dielc * PERM_VAC) * summ2
                        / kappa;
                let mut summ3 = 0.0;
                for i in 0..ncomp {
                    summ3 += x[i]
                        * q[i]
                        * q[i]
                        * (dchikap_dk[i] * dkappa_dt / t - kappa * chi[i] / t / t);
                }
                dadt_ion = -1.0 / 12.0 / PI / KB / (self.dielc * PERM_VAC) * summ3;
            }
        }

        dadt_hc + dadt_disp + dadt_assoc + dadt_polar + dadt_ion
    }

    // -- calc_compressibility_factor --------------------------------------

    /// `calc_compressibility_factor()`: Z at the live state (XA tol 1e-14
    /// here — upstream's inconsistency with alphar's 1e-15, reproduced).
    pub fn calc_compressibility_factor(&self) -> f64 {
        let ncomp = self.n;
        let t = self.t;
        let x = &self.mole_fractions;
        let p = self.prep();

        let zhs = p.zeta[3] / (1.0 - p.zeta[3])
            + 3.0 * p.zeta[1] * p.zeta[2] / p.zeta[0] / (1.0 - p.zeta[3]) / (1.0 - p.zeta[3])
            + (3.0 * p.zeta[2].powi(3) - p.zeta[3] * p.zeta[2].powi(3))
                / p.zeta[0]
                / (1.0 - p.zeta[3]).powi(3);

        let mut det_i1_det = 0.0;
        let mut det_i2_det = 0.0;
        let mut i2 = 0.0;
        for i in 0..7 {
            det_i1_det += p.a[i] * ((i + 1) as f64) * p.eta.powi(i as i32);
            det_i2_det += p.b[i] * ((i + 1) as f64) * p.eta.powi(i as i32);
            i2 += p.b[i] * p.eta.powi(i as i32);
        }

        let mut summ = 0.0;
        for i in 0..ncomp {
            summ += x[i] * (self.components[i].m - 1.0) / p.ghs[i * ncomp + i]
                * p.denghs[i * ncomp + i];
        }
        let zid = 1.0;
        let zhc = p.m_avg * zhs - summ;
        let zdisp = -2.0 * PI * p.den * det_i1_det * p.m2es3
            - PI * p.den * p.m_avg * (p.c1 * det_i2_det + p.c2 * p.eta * i2) * p.m2e2s3;

        // Polar (separate triple loop, as upstream structures this kernel)
        let mut zpolar = 0.0;
        if self.polar_term {
            let dipm_sq = self.dipm_sq();
            let mut a2t = 0.0;
            let mut da2_det_t = 0.0;
            for i in 0..ncomp {
                for j in 0..ncomp {
                    let mut m_ij = (self.components[i].m * self.components[j].m).sqrt();
                    if m_ij > 2.0 {
                        m_ij = 2.0;
                    }
                    let mut j2 = 0.0;
                    let mut det_j2_det = 0.0;
                    for l in 0..5 {
                        let adip = A0DIP[l]
                            + (m_ij - 1.0) / m_ij * A1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * A2DIP[l];
                        let bdip = B0DIP[l]
                            + (m_ij - 1.0) / m_ij * B1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * B2DIP[l];
                        j2 += (adip + bdip * p.e_ij[i * ncomp + j] / t) * p.eta.powi(l as i32);
                        det_j2_det += (adip + bdip * p.e_ij[i * ncomp + j] / t)
                            * ((l + 1) as f64)
                            * p.eta.powi(l as i32);
                    }
                    let common = x[i] * x[j] * p.e_ij[i * ncomp + i] / t * p.e_ij[j * ncomp + j]
                        / t
                        * p.s_ij[i * ncomp + i].powi(3)
                        * p.s_ij[j * ncomp + j].powi(3)
                        / p.s_ij[i * ncomp + j].powi(3)
                        * self.components[i].dipnum
                        * self.components[j].dipnum
                        * dipm_sq[i]
                        * dipm_sq[j];
                    a2t += common * j2;
                    da2_det_t += common * det_j2_det;
                }
            }
            let mut a3t = 0.0;
            let mut da3_det_t = 0.0;
            for i in 0..ncomp {
                for j in 0..ncomp {
                    for k in 0..ncomp {
                        let mut m_ijk =
                            (self.components[i].m * self.components[j].m * self.components[k].m)
                                .powf(1.0 / 3.0);
                        if m_ijk > 2.0 {
                            m_ijk = 2.0;
                        }
                        let mut j3 = 0.0;
                        let mut det_j3_det = 0.0;
                        for l in 0..5 {
                            let cdip = C0DIP[l]
                                + (m_ijk - 1.0) / m_ijk * C1DIP[l]
                                + (m_ijk - 1.0) / m_ijk * (m_ijk - 2.0) / m_ijk * C2DIP[l];
                            j3 += cdip * p.eta.powi(l as i32);
                            det_j3_det += cdip * ((l + 2) as f64) * p.eta.powi(l as i32 + 1);
                        }
                        let common = x[i] * x[j] * x[k] * p.e_ij[i * ncomp + i] / t
                            * p.e_ij[j * ncomp + j]
                            / t
                            * p.e_ij[k * ncomp + k]
                            / t
                            * p.s_ij[i * ncomp + i].powi(3)
                            * p.s_ij[j * ncomp + j].powi(3)
                            * p.s_ij[k * ncomp + k].powi(3)
                            / p.s_ij[i * ncomp + j]
                            / p.s_ij[i * ncomp + k]
                            / p.s_ij[j * ncomp + k]
                            * self.components[i].dipnum
                            * self.components[j].dipnum
                            * self.components[k].dipnum
                            * dipm_sq[i]
                            * dipm_sq[j]
                            * dipm_sq[k];
                        a3t += common * j3;
                        da3_det_t += common * det_j3_det;
                    }
                }
            }
            let a2 = -PI * p.den * a2t;
            let a3 = -4.0 / 3.0 * PI * PI * p.den * p.den * a3t;
            let da2_det = -PI * p.den / p.eta * da2_det_t;
            let da3_det = -4.0 / 3.0 * PI * PI * p.den / p.eta * p.den / p.eta * da3_det_t;
            if a2 != 0.0 {
                zpolar = p.eta
                    * ((da2_det * (1.0 - a3 / a2) + (da3_det * a2 - a3 * da2_det) / a2)
                        / (1.0 - a3 / a2)
                        / (1.0 - a3 / a2));
            }
        }

        // Association
        let mut zassoc = 0.0;
        if self.assoc_term {
            let (num_sites, i_a, x_assoc, xa, delta_ij, _ddt) = self.assoc_block(&p, false, 1e-14);
            // ddelta_dx (per component k)
            let mut ddelta_dx = vec![0.0; num_sites * num_sites * ncomp];
            let mut idx_ddelta = 0usize;
            for k in 0..ncomp {
                let mut idxa = 0usize;
                for i in 0..num_sites {
                    let idxi = i_a[i] * ncomp + i_a[i];
                    for j in 0..num_sites {
                        let idxj = i_a[j] * ncomp + i_a[j];
                        if self.assoc_matrix[idxa] != 0 {
                            let e_abij =
                                (self.components[i_a[i]].u_ab + self.components[i_a[j]].u_ab) / 2.0;
                            let vol_abij = (self.components[i_a[i]].vol_a
                                * self.components[i_a[j]].vol_a)
                                .sqrt()
                                * ((p.s_ij[idxi] * p.s_ij[idxj]).sqrt()
                                    / (0.5 * (p.s_ij[idxi] + p.s_ij[idxj])))
                                    .powi(3);
                            let dghsd_dx = PI / 6.0
                                * self.components[k].m
                                * (p.d[k].powi(3) / (1.0 - p.zeta[3]) / (1.0 - p.zeta[3])
                                    + 3.0 * p.d[i_a[i]] * p.d[i_a[j]]
                                        / (p.d[i_a[i]] + p.d[i_a[j]])
                                        * (p.d[k] * p.d[k]
                                            / (1.0 - p.zeta[3])
                                            / (1.0 - p.zeta[3])
                                            + 2.0 * p.d[k].powi(3) * p.zeta[2]
                                                / (1.0 - p.zeta[3]).powi(3))
                                    + 2.0
                                        * (p.d[i_a[i]] * p.d[i_a[j]]
                                            / (p.d[i_a[i]] + p.d[i_a[j]]))
                                            .powi(2)
                                        * (2.0 * p.d[k] * p.d[k] * p.zeta[2]
                                            / (1.0 - p.zeta[3]).powi(3)
                                            + 3.0
                                                * (p.d[k].powi(3) * p.zeta[2] * p.zeta[2]
                                                    / (1.0 - p.zeta[3]).powi(4))));
                            ddelta_dx[idx_ddelta] = dghsd_dx
                                * ((e_abij / t).exp() - 1.0)
                                * p.s_ij[i_a[i] * ncomp + i_a[j]].powi(3)
                                * vol_abij;
                        }
                        idx_ddelta += 1;
                        idxa += 1;
                    }
                }
            }

            let dxa_dx = self.dxadx_find(&delta_ij, p.den, &xa, &ddelta_dx, &x_assoc);
            let mut summ = 0.0;
            let mut ij = 0usize;
            for i in 0..ncomp {
                for j in 0..num_sites {
                    summ += x[i] * p.den * x[i_a[j]] * (1.0 / xa[j] - 0.5) * dxa_dx[ij];
                    ij += 1;
                }
            }
            zassoc = summ;
        }

        // Ion
        let mut zion = 0.0;
        if self.ion_term {
            let q: Vec<f64> = self.components.iter().map(|c| c.z * E_CHRG).collect();
            let mut summ = 0.0;
            for i in 0..ncomp {
                summ += self.components[i].z.powf(2.0) * x[i];
            }
            let kappa = (p.den * E_CHRG * E_CHRG / KB / t / (self.dielc * PERM_VAC) * summ).sqrt();
            if kappa != 0.0 {
                let mut summ2 = 0.0;
                for i in 0..ncomp {
                    let chi = 3.0 / (kappa * p.d[i]).powi(3)
                        * (1.5 + (1.0 + kappa * p.d[i]).ln() - 2.0 * (1.0 + kappa * p.d[i])
                            + 0.5 * (1.0 + kappa * p.d[i]).powi(2));
                    let sigma_k = -2.0 * chi + 3.0 / (1.0 + kappa * p.d[i]);
                    summ2 += q[i] * q[i] * x[i] * sigma_k;
                }
                zion = -1.0 * kappa / 24.0 / PI / KB / t / (self.dielc * PERM_VAC) * summ2;
            }
        }

        zid + zhc + zdisp + zpolar + zassoc + zion
    }

    /// `dXAdx_find`: dense LU solve for dXA/d(rho_i).
    fn dxadx_find(
        &self,
        delta_ij: &[f64],
        den: f64,
        xa: &[f64],
        ddelta_dx: &[f64],
        x: &[f64],
    ) -> Vec<f64> {
        let num_sites = xa.len();
        let ncomp = self.assoc_num.len();
        let dim = num_sites * ncomp;
        let mut a = vec![vec![0.0; dim]; dim];
        let mut b = vec![0.0; dim];

        let mut idx1 = 0usize;
        let mut ij = 0usize;
        for i in 0..ncomp {
            for j in 0..num_sites {
                let mut sum1 = 0.0;
                for k in 0..num_sites {
                    sum1 += den
                        * x[k]
                        * (xa[k] * ddelta_dx[i * num_sites * num_sites + j * num_sites + k]);
                    a[ij][i * num_sites + k] =
                        xa[j] * xa[j] * den * x[k] * delta_ij[j * num_sites + k];
                }
                let mut sum2 = 0.0;
                for l in 0..(self.assoc_num[i] as usize) {
                    sum2 += xa[idx1 + l] * delta_ij[idx1 * num_sites + l * num_sites + j];
                }
                a[ij][ij] += 1.0;
                b[ij] = -1.0 * xa[j] * xa[j] * (sum1 + sum2);
                ij += 1;
            }
            idx1 += self.assoc_num[i] as usize;
        }
        solve_dense(&mut a, &b)
    }

    // -- derived quantities ------------------------------------------------

    /// `calc_pressure` [Pa].
    pub fn calc_pressure(&self) -> f64 {
        let den = self.rhomolar * N_AV / 1.0e30;
        self.calc_compressibility_factor() * KB * self.t * den * 1.0e30
    }

    /// `calc_hmolar_residual` (Gross & Sadowski 2001 Eq. A.46).
    pub fn calc_hmolar_residual(&self) -> f64 {
        let z = self.calc_compressibility_factor();
        let dares_dt = self.calc_dadt();
        (-self.t * dares_dt + (z - 1.0)) * KB * N_AV * self.t
    }

    /// `calc_smolar_residual`.
    pub fn calc_smolar_residual(&self) -> f64 {
        let dares_dt = self.calc_dadt();
        let ares = self.calc_alphar();
        KB * N_AV * (-self.t * dares_dt - ares)
    }

    /// `calc_fugacity_coefficients()`: ln phi_i via composition
    /// derivatives; XA at 1e-15 here (as alphar/dadt).
    pub fn calc_fugacity_coefficients(&self) -> Vec<f64> {
        let ncomp = self.n;
        let t = self.t;
        let x = &self.mole_fractions;
        let p = self.prep();

        let ares_hs = 1.0 / p.zeta[0]
            * (3.0 * p.zeta[1] * p.zeta[2] / (1.0 - p.zeta[3])
                + p.zeta[2].powi(3) / (p.zeta[3] * (1.0 - p.zeta[3]).powi(2))
                + (p.zeta[2].powi(3) / p.zeta[3].powi(2) - p.zeta[0]) * (1.0 - p.zeta[3]).ln());
        let zhs = p.zeta[3] / (1.0 - p.zeta[3])
            + 3.0 * p.zeta[1] * p.zeta[2] / p.zeta[0] / (1.0 - p.zeta[3]) / (1.0 - p.zeta[3])
            + (3.0 * p.zeta[2].powi(3) - p.zeta[3] * p.zeta[2].powi(3))
                / p.zeta[0]
                / (1.0 - p.zeta[3]).powi(3);

        let mut det_i1_det = 0.0;
        let mut det_i2_det = 0.0;
        let mut i1 = 0.0;
        let mut i2 = 0.0;
        for i in 0..7 {
            det_i1_det += p.a[i] * ((i + 1) as f64) * p.eta.powi(i as i32);
            det_i2_det += p.b[i] * ((i + 1) as f64) * p.eta.powi(i as i32);
            i2 += p.b[i] * p.eta.powi(i as i32);
            i1 += p.a[i] * p.eta.powi(i as i32);
        }

        let mut summ = 0.0;
        for i in 0..ncomp {
            summ += x[i] * (self.components[i].m - 1.0) * p.ghs[i * ncomp + i].ln();
        }
        let ares_hc = p.m_avg * ares_hs - summ;
        let ares_disp =
            -2.0 * PI * p.den * i1 * p.m2es3 - PI * p.den * p.m_avg * p.c1 * i2 * p.m2e2s3;

        let mut summ = 0.0;
        for i in 0..ncomp {
            summ += x[i] * (self.components[i].m - 1.0) / p.ghs[i * ncomp + i]
                * p.denghs[i * ncomp + i];
        }
        let zhc = p.m_avg * zhs - summ;
        let zdisp = -2.0 * PI * p.den * det_i1_det * p.m2es3
            - PI * p.den * p.m_avg * (p.c1 * det_i2_det + p.c2 * p.eta * i2) * p.m2e2s3;

        // Composition derivatives of the hs/hc/dispersion parts
        let mut dghsii_dx = vec![0.0; ncomp * ncomp];
        let mut dahs_dx = vec![0.0; ncomp];
        let mut idx = 0usize;
        for i in 0..ncomp {
            let mut dzeta_dx = [0.0; 4];
            for l in 0..4 {
                dzeta_dx[l] = PI / 6.0 * p.den * self.components[i].m * p.d[i].powi(l as i32);
            }
            for j in 0..ncomp {
                dghsii_dx[idx] = dzeta_dx[3] / (1.0 - p.zeta[3]) / (1.0 - p.zeta[3])
                    + (p.d[j] * p.d[j] / (p.d[j] + p.d[j]))
                        * (3.0 * dzeta_dx[2] / (1.0 - p.zeta[3]) / (1.0 - p.zeta[3])
                            + 6.0 * p.zeta[2] * dzeta_dx[3] / (1.0 - p.zeta[3]).powi(3))
                    + (p.d[j] * p.d[j] / (p.d[j] + p.d[j])).powi(2)
                        * (4.0 * p.zeta[2] * dzeta_dx[2] / (1.0 - p.zeta[3]).powi(3)
                            + 6.0 * p.zeta[2] * p.zeta[2] * dzeta_dx[3]
                                / (1.0 - p.zeta[3]).powi(4));
                idx += 1;
            }
            dahs_dx[i] = -dzeta_dx[0] / p.zeta[0] * ares_hs
                + 1.0 / p.zeta[0]
                    * (3.0 * (dzeta_dx[1] * p.zeta[2] + p.zeta[1] * dzeta_dx[2])
                        / (1.0 - p.zeta[3])
                        + 3.0 * p.zeta[1] * p.zeta[2] * dzeta_dx[3]
                            / (1.0 - p.zeta[3])
                            / (1.0 - p.zeta[3])
                        + 3.0 * p.zeta[2] * p.zeta[2] * dzeta_dx[2]
                            / p.zeta[3]
                            / (1.0 - p.zeta[3])
                            / (1.0 - p.zeta[3])
                        + p.zeta[2].powi(3) * dzeta_dx[3] * (3.0 * p.zeta[3] - 1.0)
                            / p.zeta[3]
                            / p.zeta[3]
                            / (1.0 - p.zeta[3]).powi(3)
                        + (1.0 - p.zeta[3]).ln()
                            * ((3.0 * p.zeta[2] * p.zeta[2] * dzeta_dx[2] * p.zeta[3]
                                - 2.0 * p.zeta[2].powi(3) * dzeta_dx[3])
                                / p.zeta[3].powi(3)
                                - dzeta_dx[0])
                        + (p.zeta[0] - p.zeta[2].powi(3) / p.zeta[3] / p.zeta[3]) * dzeta_dx[3]
                            / (1.0 - p.zeta[3]));
        }

        let mut dadisp_dx = vec![0.0; ncomp];
        let mut dahc_dx = vec![0.0; ncomp];
        for i in 0..ncomp {
            let dzeta3_dx = PI / 6.0 * p.den * self.components[i].m * p.d[i].powi(3);
            let mut di1_dx = 0.0;
            let mut di2_dx = 0.0;
            let mut dm2es3_dx = 0.0;
            let mut dm2e2s3_dx = 0.0;
            for l in 0..7 {
                let daa_dx = self.components[i].m / p.m_avg / p.m_avg * A1[l]
                    + self.components[i].m / p.m_avg / p.m_avg * (3.0 - 4.0 / p.m_avg) * A2[l];
                let db_dx = self.components[i].m / p.m_avg / p.m_avg * B1[l]
                    + self.components[i].m / p.m_avg / p.m_avg * (3.0 - 4.0 / p.m_avg) * B2[l];
                di1_dx += p.a[l] * (l as f64) * dzeta3_dx * p.eta.powi(l as i32 - 1)
                    + daa_dx * p.eta.powi(l as i32);
                di2_dx += p.b[l] * (l as f64) * dzeta3_dx * p.eta.powi(l as i32 - 1)
                    + db_dx * p.eta.powi(l as i32);
            }
            for j in 0..ncomp {
                dm2es3_dx += x[j]
                    * self.components[j].m
                    * (p.e_ij[i * ncomp + j] / t)
                    * p.s_ij[i * ncomp + j].powi(3);
                dm2e2s3_dx += x[j]
                    * self.components[j].m
                    * (p.e_ij[i * ncomp + j] / t).powi(2)
                    * p.s_ij[i * ncomp + j].powi(3);
                dahc_dx[i] += x[j] * (self.components[j].m - 1.0) / p.ghs[j * ncomp + j]
                    * dghsii_dx[i * ncomp + j];
            }
            let dm2es3_dx = dm2es3_dx * 2.0 * self.components[i].m;
            let dm2e2s3_dx = dm2e2s3_dx * 2.0 * self.components[i].m;
            dahc_dx[i] = self.components[i].m * ares_hs + p.m_avg * dahs_dx[i]
                - dahc_dx[i]
                - (self.components[i].m - 1.0) * p.ghs[i * ncomp + i].ln();
            let dc1_dx = p.c2 * dzeta3_dx
                - p.c1
                    * p.c1
                    * (self.components[i].m * (8.0 * p.eta - 2.0 * p.eta * p.eta)
                        / (1.0 - p.eta).powi(4)
                        - self.components[i].m
                            * (20.0 * p.eta - 27.0 * p.eta * p.eta + 12.0 * p.eta.powi(3)
                                - 2.0 * p.eta.powi(4))
                            / ((1.0 - p.eta) * (2.0 - p.eta)).powi(2));

            dadisp_dx[i] = -2.0 * PI * p.den * (di1_dx * p.m2es3 + i1 * dm2es3_dx)
                - PI * p.den
                    * ((self.components[i].m * p.c1 * i2
                        + p.m_avg * dc1_dx * i2
                        + p.m_avg * p.c1 * di2_dx)
                        * p.m2e2s3
                        + p.m_avg * p.c1 * i2 * dm2e2s3_dx);
        }

        let mut mu_hc = vec![0.0; ncomp];
        let mut mu_disp = vec![0.0; ncomp];
        for i in 0..ncomp {
            for j in 0..ncomp {
                mu_hc[i] += x[j] * dahc_dx[j];
                mu_disp[i] += x[j] * dadisp_dx[j];
            }
            mu_hc[i] = ares_hc + zhc + dahc_dx[i] - mu_hc[i];
            mu_disp[i] = ares_disp + zdisp + dadisp_dx[i] - mu_disp[i];
        }

        // Polar
        let mut mu_polar = vec![0.0; ncomp];
        if self.polar_term {
            let dipm_sq = self.dipm_sq();
            let mut a2t = 0.0;
            let mut a3t = 0.0;
            let mut da2_det_t = 0.0;
            let mut da3_det_t = 0.0;
            let mut da2_dx = vec![0.0; ncomp];
            let mut da3_dx = vec![0.0; ncomp];
            for i in 0..ncomp {
                for j in 0..ncomp {
                    let mut m_ij = (self.components[i].m * self.components[j].m).sqrt();
                    if m_ij > 2.0 {
                        m_ij = 2.0;
                    }
                    let mut j2 = 0.0;
                    let mut dj2_det = 0.0;
                    let mut det_j2_det = 0.0;
                    for l in 0..5 {
                        let adip = A0DIP[l]
                            + (m_ij - 1.0) / m_ij * A1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * A2DIP[l];
                        let bdip = B0DIP[l]
                            + (m_ij - 1.0) / m_ij * B1DIP[l]
                            + (m_ij - 1.0) / m_ij * (m_ij - 2.0) / m_ij * B2DIP[l];
                        let term = adip + bdip * p.e_ij[i * ncomp + j] / t;
                        j2 += term * p.eta.powi(l as i32);
                        dj2_det += term * (l as f64) * p.eta.powi(l as i32 - 1);
                        det_j2_det += term * ((l + 1) as f64) * p.eta.powi(l as i32);
                    }
                    let base2 = p.e_ij[i * ncomp + i] / t * p.e_ij[j * ncomp + j] / t
                        * p.s_ij[i * ncomp + i].powi(3)
                        * p.s_ij[j * ncomp + j].powi(3)
                        / p.s_ij[i * ncomp + j].powi(3)
                        * self.components[i].dipnum
                        * self.components[j].dipnum
                        * dipm_sq[i]
                        * dipm_sq[j];
                    a2t += x[i] * x[j] * base2 * j2;
                    da2_det_t += x[i] * x[j] * base2 * det_j2_det;
                    if i == j {
                        da2_dx[i] += base2
                            * (x[i] * x[j] * dj2_det * PI / 6.0
                                * p.den
                                * self.components[i].m
                                * p.d[i].powi(3)
                                + 2.0 * x[j] * j2);
                    } else {
                        da2_dx[i] += base2
                            * (x[i] * x[j] * dj2_det * PI / 6.0
                                * p.den
                                * self.components[i].m
                                * p.d[i].powi(3)
                                + x[j] * j2);
                    }

                    for k in 0..ncomp {
                        let mut m_ijk =
                            (self.components[i].m * self.components[j].m * self.components[k].m)
                                .powf(1.0 / 3.0);
                        if m_ijk > 2.0 {
                            m_ijk = 2.0;
                        }
                        let mut j3 = 0.0;
                        let mut dj3_det = 0.0;
                        let mut det_j3_det = 0.0;
                        for l in 0..5 {
                            let cdip = C0DIP[l]
                                + (m_ijk - 1.0) / m_ijk * C1DIP[l]
                                + (m_ijk - 1.0) / m_ijk * (m_ijk - 2.0) / m_ijk * C2DIP[l];
                            j3 += cdip * p.eta.powi(l as i32);
                            dj3_det += cdip * (l as f64) * p.eta.powi(l as i32 - 1);
                            det_j3_det += cdip * ((l + 2) as f64) * p.eta.powi(l as i32 + 1);
                        }
                        let base3 = p.e_ij[i * ncomp + i] / t * p.e_ij[j * ncomp + j] / t
                            * p.e_ij[k * ncomp + k]
                            / t
                            * p.s_ij[i * ncomp + i].powi(3)
                            * p.s_ij[j * ncomp + j].powi(3)
                            * p.s_ij[k * ncomp + k].powi(3)
                            / p.s_ij[i * ncomp + j]
                            / p.s_ij[i * ncomp + k]
                            / p.s_ij[j * ncomp + k]
                            * self.components[i].dipnum
                            * self.components[j].dipnum
                            * self.components[k].dipnum
                            * dipm_sq[i]
                            * dipm_sq[j]
                            * dipm_sq[k];
                        a3t += x[i] * x[j] * x[k] * base3 * j3;
                        da3_det_t += x[i] * x[j] * x[k] * base3 * det_j3_det;
                        if i == j && i == k {
                            da3_dx[i] += base3
                                * (x[i] * x[j] * x[k] * dj3_det * PI / 6.0
                                    * p.den
                                    * self.components[i].m
                                    * p.d[i].powi(3)
                                    + 3.0 * x[j] * x[k] * j3);
                        } else if i == j || i == k {
                            da3_dx[i] += base3
                                * (x[i] * x[j] * x[k] * dj3_det * PI / 6.0
                                    * p.den
                                    * self.components[i].m
                                    * p.d[i].powi(3)
                                    + 2.0 * x[j] * x[k] * j3);
                        } else {
                            da3_dx[i] += base3
                                * (x[i] * x[j] * x[k] * dj3_det * PI / 6.0
                                    * p.den
                                    * self.components[i].m
                                    * p.d[i].powi(3)
                                    + x[j] * x[k] * j3);
                        }
                    }
                }
            }
            let a2 = -PI * p.den * a2t;
            let a3 = -4.0 / 3.0 * PI * PI * p.den * p.den * a3t;
            let da2_det = -PI * p.den / p.eta * da2_det_t;
            let da3_det = -4.0 / 3.0 * PI * PI * p.den / p.eta * p.den / p.eta * da3_det_t;
            for i in 0..ncomp {
                da2_dx[i] = -PI * p.den * da2_dx[i];
                da3_dx[i] = -4.0 / 3.0 * PI * PI * p.den * p.den * da3_dx[i];
            }
            let mut dapolar_dx = vec![0.0; ncomp];
            for i in 0..ncomp {
                dapolar_dx[i] = (da2_dx[i] * (1.0 - a3 / a2)
                    + (da3_dx[i] * a2 - a3 * da2_dx[i]) / a2)
                    / (1.0 - a3 / a2).powi(2);
            }
            if a2 != 0.0 {
                let ares_polar = a2 / (1.0 - a3 / a2);
                let zpolar = p.eta
                    * ((da2_det * (1.0 - a3 / a2) + (da3_det * a2 - a3 * da2_det) / a2)
                        / (1.0 - a3 / a2)
                        / (1.0 - a3 / a2));
                for i in 0..ncomp {
                    for j in 0..ncomp {
                        mu_polar[i] += x[j] * dapolar_dx[j];
                    }
                    mu_polar[i] = ares_polar + zpolar + dapolar_dx[i] - mu_polar[i];
                }
            }
        }

        // Association
        let mut mu_assoc = vec![0.0; ncomp];
        if self.assoc_term {
            let (num_sites, i_a, x_assoc, xa, delta_ij, _ddt) = self.assoc_block(&p, false, 1e-15);
            let mut ddelta_dx = vec![0.0; num_sites * num_sites * ncomp];
            let mut idx_ddelta = 0usize;
            for k in 0..ncomp {
                let mut idxa = 0usize;
                for i in 0..num_sites {
                    let idxi = i_a[i] * ncomp + i_a[i];
                    for j in 0..num_sites {
                        let idxj = i_a[j] * ncomp + i_a[j];
                        if self.assoc_matrix[idxa] != 0 {
                            let e_abij =
                                (self.components[i_a[i]].u_ab + self.components[i_a[j]].u_ab) / 2.0;
                            let vol_abij = (self.components[i_a[i]].vol_a
                                * self.components[i_a[j]].vol_a)
                                .sqrt()
                                * ((p.s_ij[idxi] * p.s_ij[idxj]).sqrt()
                                    / (0.5 * (p.s_ij[idxi] + p.s_ij[idxj])))
                                    .powi(3);
                            let dghsd_dx = PI / 6.0
                                * self.components[k].m
                                * (p.d[k].powi(3) / (1.0 - p.zeta[3]) / (1.0 - p.zeta[3])
                                    + 3.0 * p.d[i_a[i]] * p.d[i_a[j]]
                                        / (p.d[i_a[i]] + p.d[i_a[j]])
                                        * (p.d[k] * p.d[k]
                                            / (1.0 - p.zeta[3])
                                            / (1.0 - p.zeta[3])
                                            + 2.0 * p.d[k].powi(3) * p.zeta[2]
                                                / (1.0 - p.zeta[3]).powi(3))
                                    + 2.0
                                        * (p.d[i_a[i]] * p.d[i_a[j]]
                                            / (p.d[i_a[i]] + p.d[i_a[j]]))
                                            .powi(2)
                                        * (2.0 * p.d[k] * p.d[k] * p.zeta[2]
                                            / (1.0 - p.zeta[3]).powi(3)
                                            + 3.0
                                                * (p.d[k].powi(3) * p.zeta[2] * p.zeta[2]
                                                    / (1.0 - p.zeta[3]).powi(4))));
                            ddelta_dx[idx_ddelta] = dghsd_dx
                                * ((e_abij / t).exp() - 1.0)
                                * p.s_ij[i_a[i] * ncomp + i_a[j]].powi(3)
                                * vol_abij;
                        }
                        idx_ddelta += 1;
                        idxa += 1;
                    }
                }
            }
            let dxa_dx = self.dxadx_find(&delta_ij, p.den, &xa, &ddelta_dx, &x_assoc);
            let mut ij = 0usize;
            for i in 0..ncomp {
                for j in 0..num_sites {
                    mu_assoc[i] += x[i_a[j]] * p.den * dxa_dx[ij] * (1.0 / xa[j] - 0.5);
                    ij += 1;
                }
            }
            for i in 0..num_sites {
                mu_assoc[i_a[i]] += xa[i].ln() - 0.5 * xa[i] + 0.5;
            }
        }

        // Ion
        let mut mu_ion = vec![0.0; ncomp];
        if self.ion_term {
            let q: Vec<f64> = self.components.iter().map(|c| c.z * E_CHRG).collect();
            let mut summ = 0.0;
            for i in 0..ncomp {
                summ += self.components[i].z * self.components[i].z * x[i];
            }
            let kappa = (p.den * E_CHRG * E_CHRG / KB / t / (self.dielc * PERM_VAC) * summ).sqrt();
            if kappa != 0.0 {
                let mut chi = vec![0.0; ncomp];
                let mut sigma_k = vec![0.0; ncomp];
                let mut summ1 = 0.0;
                let mut summ2 = 0.0;
                for i in 0..ncomp {
                    chi[i] = 3.0 / (kappa * p.d[i]).powi(3)
                        * (1.5 + (1.0 + kappa * p.d[i]).ln() - 2.0 * (1.0 + kappa * p.d[i])
                            + 0.5 * (1.0 + kappa * p.d[i]).powi(2));
                    sigma_k[i] = -2.0 * chi[i] + 3.0 / (1.0 + kappa * p.d[i]);
                    summ1 += q[i] * q[i] * x[i] * sigma_k[i];
                    summ2 += x[i] * q[i] * q[i];
                }
                for i in 0..ncomp {
                    mu_ion[i] = -q[i] * q[i] * kappa / 24.0 / PI / KB / t / (self.dielc * PERM_VAC)
                        * (2.0 * chi[i] + summ1 / summ2);
                }
            }
        }

        let z = self.calc_compressibility_factor();
        let mut fugcoef = vec![0.0; ncomp];
        for i in 0..ncomp {
            let mu = mu_hc[i] + mu_disp[i] + mu_polar[i] + mu_assoc[i] + mu_ion[i];
            fugcoef[i] = (mu - z.ln()).exp();
        }
        fugcoef
    }

    /// `calc_gibbsmolar_residual` — (T, V) basis (upstream #1943).
    pub fn calc_gibbsmolar_residual(&self) -> f64 {
        let z = self.calc_compressibility_factor();
        let ares = self.calc_alphar();
        (ares + (z - 1.0)) * KB * N_AV * self.t
    }
}

impl PcsaftBackend {
    /// `update_DmolarT(rho)`: sets the live density, returns the pressure.
    pub fn update_dmolar_t(&mut self, rho: f64) -> f64 {
        self.rhomolar = rho;
        self.calc_pressure()
    }

    /// `solver_rho_Tp(T, p, phase)`: reduced-density bracket scan + Brent.
    /// Mutates the live `_rhomolar` while probing, exactly as upstream.
    /// Returns `f64::MAX`-like garbage only via upstream's own no-else hole
    /// (a phase that is neither liquid- nor gas-like with 2-3 brackets).
    pub fn solver_rho_tp(&mut self, t: f64, p: f64, phase: PcsaftPhase) -> Result<f64> {
        self.t = t;
        // bracket scan over reduced density nu
        let mut x_lo: Vec<f64> = Vec::new();
        let mut x_hi: Vec<f64> = Vec::new();
        let num_pts = 20;
        let mut limit_lower = -8.0_f64;
        let mut limit_upper = -1.0_f64;
        let mut rho_guess = 1e-13_f64;
        let mut rho_guess_prev = rho_guess;
        let rm = self.reduced_to_molar(rho_guess, t);
        let mut err_prev = (self.update_dmolar_t(rm) - p) / p;
        for i in 0..num_pts {
            rho_guess = 10.0_f64
                .powf((limit_upper - limit_lower) / (num_pts as f64) * (i as f64) + limit_lower);
            let rm = self.reduced_to_molar(rho_guess, t);
            let err = (self.update_dmolar_t(rm) - p) / p;
            if err * err_prev < 0.0 {
                x_lo.push(rho_guess_prev);
                x_hi.push(rho_guess);
            }
            err_prev = err;
            rho_guess_prev = rho_guess;
        }
        limit_lower = 0.1;
        limit_upper = 0.7405;
        for i in 0..num_pts {
            rho_guess = (limit_upper - limit_lower) / (num_pts as f64) * (i as f64) + limit_lower;
            let rm = self.reduced_to_molar(rho_guess, t);
            let err = (self.update_dmolar_t(rm) - p) / p;
            if err * err_prev < 0.0 {
                x_lo.push(rho_guess_prev);
                x_hi.push(rho_guess);
            }
            err_prev = err;
            rho_guess_prev = rho_guess;
        }

        let mut brent_resid = |rho: f64, this: &mut Self| -> f64 {
            let peos = this.update_dmolar_t(rho);
            let cost = (peos - p) / p;
            if cost.is_finite() { cost } else { 1.0e20 }
        };

        let mut rho = f64::MAX; // upstream _HUGE
        if x_lo.len() == 1 {
            let a = self.reduced_to_molar(x_lo[0], t);
            let b = self.reduced_to_molar(x_hi[0], t);
            rho = brent_mut(self, &mut brent_resid, a, b, f64::EPSILON, 1e-8, 200)?;
        } else if x_lo.len() <= 3 && !x_lo.is_empty() {
            if matches!(
                phase,
                PcsaftPhase::Liquid | PcsaftPhase::SupercriticalLiquid
            ) {
                let a = self.reduced_to_molar(*x_lo.last().unwrap(), t);
                let b = self.reduced_to_molar(*x_hi.last().unwrap(), t);
                rho = brent_mut(self, &mut brent_resid, a, b, f64::EPSILON, 1e-8, 200)?;
            } else if matches!(
                phase,
                PcsaftPhase::Gas | PcsaftPhase::SupercriticalGas | PcsaftPhase::Supercritical
            ) {
                let a = self.reduced_to_molar(x_lo[0], t);
                let b = self.reduced_to_molar(x_hi[0], t);
                rho = brent_mut(self, &mut brent_resid, a, b, f64::EPSILON, 1e-8, 200)?;
            }
        } else if x_lo.len() > 3 {
            // Minimum-Gibbs root selection (Privat/Gani/Jaubert 2010), on
            // the (T, P) basis via the local -ln(Z)*RT conversion.
            let mut g_min = 1e60;
            for i in 0..x_lo.len() {
                let a = self.reduced_to_molar(x_lo[i], t);
                let b = self.reduced_to_molar(x_hi[i], t);
                let rho_i = brent_mut(self, &mut brent_resid, a, b, f64::EPSILON, 1e-8, 200)?;
                let rho_original = self.rhomolar;
                self.rhomolar = rho_i;
                let z_i = self.calc_compressibility_factor();
                let g_i = self.calc_gibbsmolar_residual() - z_i.ln() * KB * N_AV * self.t;
                self.rhomolar = rho_original;
                if g_i < g_min {
                    g_min = g_i;
                    rho = rho_i;
                }
            }
        } else {
            // 0 brackets: minimum-|error| scan (may be far from a root)
            let num_pts = 25;
            let mut err_min = 1e40;
            let mut rho_min = f64::NAN;
            for i in 0..num_pts {
                let rho_guess = (0.7405 - 1e-8) / (num_pts as f64) * (i as f64) + 1e-8;
                let rm = self.reduced_to_molar(rho_guess, t);
                let err = (self.update_dmolar_t(rm) - p) / p;
                if err.abs() < err_min {
                    err_min = err.abs();
                    rho_min = self.reduced_to_molar(rho_guess, t);
                }
            }
            rho = rho_min;
        }
        Ok(rho)
    }
}

/// Upstream `phases` values the PC-SAFT paths distinguish.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PcsaftPhase {
    Liquid,
    Gas,
    Supercritical,
    SupercriticalGas,
    SupercriticalLiquid,
    TwoPhase,
    Unknown,
}

/// Upstream `Brent` (Solvers.cpp) over a state-mutating residual —
/// operation-for-operation identical to the verified rustprop-heos port.
fn brent_mut<F: FnMut(f64, &mut PcsaftBackend) -> f64>(
    state: &mut PcsaftBackend,
    call: &mut F,
    mut a: f64,
    mut b: f64,
    macheps: f64,
    t: f64,
    maxiter: i32,
) -> Result<f64> {
    let mut fa = call(a, state);
    let mut fb = call(b, state);
    if fb.abs() < t {
        return Ok(b);
    }
    if !fb.is_finite() {
        return Err(Error::Value(format!(
            "Brent's method f(b) is NAN for b = {b}, other input was a = {a}"
        )));
    }
    if fa.abs() < t {
        return Ok(a);
    }
    if !fa.is_finite() {
        return Err(Error::Value(format!(
            "Brent's method f(a) is NAN for a = {a}, other input was b = {b}"
        )));
    }
    if fa * fb > 0.0 {
        return Err(Error::Value(format!(
            "Inputs in Brent [{a:.6e},{b:.6e}] do not bracket the root.  Function values are [{fa:.6e},{fb:.6e}]"
        )));
    }

    let mut c = a;
    let mut fc = fa;
    let mut iter = 1;
    if fc.abs() < fb.abs() {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
    }
    let mut d = b - a;
    let mut e = b - a;
    let mut m = 0.5 * (c - b);
    let mut tol = 2.0 * macheps * b.abs() + t;
    while m.abs() > tol && fb != 0.0 {
        if e.abs() < tol || fa.abs() <= fb.abs() {
            m = 0.5 * (c - b);
            d = m;
            e = m;
        } else {
            let mut p;
            let mut q;
            let mut s = fb / fa;
            if a == c {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                let r = fb / fc;
                m = 0.5 * (c - b);
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            s = e;
            e = d;
            m = 0.5 * (c - b);
            if 2.0 * p < 3.0 * m * q - (tol * q).abs() || p < (0.5 * s * q).abs() {
                d = p / q;
            } else {
                m = 0.5 * (c - b);
                d = m;
                e = m;
            }
        }
        a = b;
        fa = fb;
        if d.abs() > tol {
            b += d;
        } else if m > 0.0 {
            b += tol;
        } else {
            b += -tol;
        }
        fb = call(b, state);
        if !fb.is_finite() {
            return Err(Error::Value(format!(
                "Brent's method f(t) is NAN for t = {b}"
            )));
        }
        if fb.abs() < macheps {
            return Ok(b);
        }
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        m = 0.5 * (c - b);
        tol = 2.0 * macheps * b.abs() + t;
        iter += 1;
        if iter > maxiter {
            return Err(Error::Solution(format!(
                "Brent's method reached maximum number of steps of {maxiter} "
            )));
        }
        if fb.abs() < 2.0 * macheps * b.abs() {
            return Ok(b);
        }
    }
    Ok(b)
}

/// One shared prep of the quantities upstream recomputes in each kernel.
struct Prep {
    d: Vec<f64>,
    #[allow(dead_code)]
    dd_dt: Vec<f64>,
    den: f64,
    zeta: [f64; 4],
    dzeta_dt: [f64; 4],
    eta: f64,
    m_avg: f64,
    ghs: Vec<f64>,
    dghs_dt: Vec<f64>,
    denghs: Vec<f64>,
    e_ij: Vec<f64>,
    s_ij: Vec<f64>,
    m2es3: f64,
    m2e2s3: f64,
    a: [f64; 7],
    b: [f64; 7],
    c1: f64,
    c2: f64,
}

/// Dense partial-pivot LU solve (stands in for Eigen's `A.lu().solve(B)`).
fn solve_dense(a: &mut [Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut aug: Vec<Vec<f64>> = (0..n)
        .map(|row| {
            let mut v = a[row].clone();
            v.push(b[row]);
            v
        })
        .collect();
    for col in 0..n {
        let mut piv = col;
        for row in (col + 1)..n {
            if aug[row][col].abs() > aug[piv][col].abs() {
                piv = row;
            }
        }
        aug.swap(col, piv);
        for row in (col + 1)..n {
            let factor = aug[row][col] / aug[col][col];
            for k in col..=n {
                aug[row][k] -= factor * aug[col][k];
            }
        }
    }
    let mut v = vec![0.0; n];
    for row in (0..n).rev() {
        let mut sum = aug[row][n];
        for col in (row + 1)..n {
            sum -= aug[row][col] * v[col];
        }
        v[row] = sum / aug[row][row];
    }
    v
}
