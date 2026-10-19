//! Mixture critical points — upstream
//! `HelmholtzEOSMixtureBackend::calc_all_critical_points` (the Deiters & Bell
//! criticality-contour tracer) behind `calc_T_critical` / `calc_p_critical` /
//! `calc_rhomolar_critical` for mixtures, @ v8.0.0.
//!
//! The method: walk the `L1* = det(L*) = 0` contour in (tau, delta) from
//! `delta = 0.5`, and wherever `M1* = det(M*)` changes sign along it, solve
//! `L1* = M1* = 0` by a 2-D Newton step. `L*` is the matrix of
//! `n (d ln f_i / d n_j)` and `M*` is `L*` with its last row replaced by the
//! traces of `adj(L*) * n^2 d L*/d n_i`. Everything below is the
//! `MixtureDerivatives` layer that builds them (Kunz & Wagner 2012 Table B,
//! Bell & Jäger 2016), restricted to the `XN_INDEPENDENT` branch upstream
//! calls with.
//!
//! Two arithmetic divergences, both below the 1e-10 Newton tolerance the
//! answer converges to: Eigen's dynamic `determinant()` (partial-pivot LU)
//! and `trace(A*B)` reduce in vectorized order, and the 2-D Newton step
//! solves the 2x2 system directly where upstream calls
//! `colPivHouseholderQr`.

// Function names mirror upstream's `MixtureDerivatives`; loops keep its
// index shapes for reviewability against the C++.
#![allow(non_snake_case)]
#![allow(clippy::needless_range_loop)]

use crate::alpha::HelmholtzDerivs;
use crate::mixture::{MixtureModel, XnFlag};
use crate::solvers::{Resid1D, brent, halley};
use rustprop_core::{Error, Result};
use std::f64::consts::PI;

/// One critical point (upstream `CriticalState`, without the unused
/// `stable` flag — `calc_T_critical` and friends never read it).
#[derive(Debug, Clone, Copy)]
pub struct CriticalPoint {
    pub t: f64,
    pub p: f64,
    pub rhomolar: f64,
}

// ---------------------------------------------------------------------------
// Dense helpers (Eigen stand-ins)
// ---------------------------------------------------------------------------

/// Row-major square matrix.
#[derive(Clone)]
struct Mat {
    n: usize,
    a: Vec<f64>,
}

impl Mat {
    fn zeros(n: usize) -> Self {
        Mat {
            n,
            a: vec![0.0; n * n],
        }
    }
    fn at(&self, i: usize, j: usize) -> f64 {
        self.a[i * self.n + j]
    }
    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] = v;
    }
    /// Build a symmetric matrix from its upper triangle, upstream's fill.
    fn symmetric(n: usize, mut f: impl FnMut(usize, usize) -> f64) -> Self {
        let mut m = Mat::zeros(n);
        for i in 0..n {
            for j in i..n {
                m.set(i, j, f(i, j));
            }
        }
        for i in 0..n {
            for j in 0..i {
                let v = m.at(j, i);
                m.set(i, j, v);
            }
        }
        m
    }
    /// `(self * other).trace()`.
    fn trace_of_product(&self, other: &Mat) -> f64 {
        let n = self.n;
        let mut t = 0.0;
        for i in 0..n {
            let mut s = 0.0;
            for k in 0..n {
                s += self.at(i, k) * other.at(k, i);
            }
            t += s;
        }
        t
    }
    /// Eigen `determinant()` on a dynamic matrix: partial-pivot LU, the
    /// product of the pivots and the transposition sign.
    fn determinant(&self) -> f64 {
        let n = self.n;
        if n == 0 {
            return 1.0;
        }
        let mut lu = self.a.clone();
        let mut transpositions = 0;
        for k in 0..n {
            let mut row = k;
            let mut biggest = lu[k * n + k].abs();
            for r in (k + 1)..n {
                let v = lu[r * n + k].abs();
                if v > biggest {
                    biggest = v;
                    row = r;
                }
            }
            if biggest != 0.0 {
                if row != k {
                    for c in 0..n {
                        lu.swap(k * n + c, row * n + c);
                    }
                    transpositions += 1;
                }
                let pivot = lu[k * n + k];
                for r in (k + 1)..n {
                    lu[r * n + k] /= pivot;
                }
            }
            for r in (k + 1)..n {
                let l = lu[r * n + k];
                for c in (k + 1)..n {
                    lu[r * n + c] -= l * lu[k * n + c];
                }
            }
        }
        let mut det = 1.0;
        for k in 0..n {
            det *= lu[k * n + k];
        }
        if transpositions % 2 == 1 { -det } else { det }
    }
    /// Upstream `minor_matrix`: drop row `i` and column `j`.
    fn minor(&self, i: usize, j: usize) -> Mat {
        let n = self.n;
        let mut m = Mat::zeros(n - 1);
        let mut rr = 0;
        for r in 0..n {
            if r == i {
                continue;
            }
            let mut cc = 0;
            for c in 0..n {
                if c == j {
                    continue;
                }
                m.set(rr, cc, self.at(r, c));
                cc += 1;
            }
            rr += 1;
        }
        m
    }
    /// Upstream `adjugate`.
    fn adjugate(&self) -> Mat {
        let n = self.n;
        if n == 1 {
            let mut m = Mat::zeros(1);
            m.set(0, 0, 1.0);
            return m;
        }
        let mut adj = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
                adj.set(i, j, sign * self.minor(j, i).determinant());
            }
        }
        adj
    }
    /// Upstream `adjugate_derivative`.
    fn adjugate_derivative(&self, d: &Mat) -> Mat {
        let n = self.n;
        let mut adj = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
                let v = self.minor(j, i).adjugate().trace_of_product(&d.minor(j, i));
                adj.set(i, j, sign * v);
            }
        }
        adj
    }
}

/// `HelmholtzDerivs` by (delta order, tau order).
fn g(d: &HelmholtzDerivs, nd: u8, nt: u8) -> f64 {
    match (nd, nt) {
        (0, 0) => d.d00,
        (1, 0) => d.d10,
        (0, 1) => d.d01,
        (2, 0) => d.d20,
        (1, 1) => d.d11,
        (0, 2) => d.d02,
        (3, 0) => d.d30,
        (2, 1) => d.d21,
        (1, 2) => d.d12,
        (0, 3) => d.d03,
        (4, 0) => d.d40,
        (3, 1) => d.d31,
        (2, 2) => d.d22,
        (1, 3) => d.d13,
        (0, 4) => d.d04,
        _ => unreachable!("derivative order above four"),
    }
}

/// Composition-only reducing quantities, computed once: the composition is
/// fixed for the whole critical search.
struct Reducing {
    tr: f64,
    rhor: f64,
    dtr: Vec<f64>,
    drho: Vec<f64>,
    ndtr: Vec<f64>,
    ndrho: Vec<f64>,
    d_ndtr: Vec<f64>,
    d_ndrho: Vec<f64>,
    psi_t: Vec<f64>,
    psi_rho: Vec<f64>,
    d_psi_t: Vec<f64>,
    d_psi_rho: Vec<f64>,
    d2_psi_t: Vec<f64>,
    d2_psi_rho: Vec<f64>,
}

// ---------------------------------------------------------------------------
// The state the derivatives are taken at (upstream's `critical_state`)
// ---------------------------------------------------------------------------

struct CritState<'m> {
    model: &'m MixtureModel,
    x: Vec<f64>,
    n: usize,
    red: Reducing,
    t: f64,
    rhomolar: f64,
    tau: f64,
    delta: f64,
    ar: HelmholtzDerivs,
    comp: Vec<HelmholtzDerivs>,
    /// `F_ij * departure_ij` derivatives, symmetric, `None` with no pair.
    pair: Vec<Option<HelmholtzDerivs>>,
}

impl<'m> CritState<'m> {
    fn new(model: &'m MixtureModel, x: &[f64]) -> Self {
        let n = x.len();
        let r = &model.reducing;
        let f = XnFlag::Independent;
        let idx2 = |i: usize, j: usize| i * n + j;
        let mut red = Reducing {
            tr: r.tr(x),
            rhor: r.rhormolar(x),
            dtr: (0..n).map(|i| r.dtrdxi__constxj(x, i, f)).collect(),
            drho: (0..n).map(|i| r.drhormolardxi__constxj(x, i, f)).collect(),
            ndtr: (0..n).map(|i| r.ndtrdni__constnj(x, i, f)).collect(),
            ndrho: (0..n).map(|i| r.ndrhorbardni__constnj(x, i, f)).collect(),
            d_ndtr: vec![0.0; n * n],
            d_ndrho: vec![0.0; n * n],
            psi_t: (0..n).map(|i| r.psi_t(x, i)).collect(),
            psi_rho: (0..n).map(|i| r.psi_rho(x, i)).collect(),
            d_psi_t: vec![0.0; n * n],
            d_psi_rho: vec![0.0; n * n],
            d2_psi_t: vec![0.0; n * n * n],
            d2_psi_rho: vec![0.0; n * n * n],
        };
        for i in 0..n {
            for j in 0..n {
                red.d_ndtr[idx2(i, j)] = r.d_ndtrdni_dxj__constxi(x, i, j, f);
                red.d_ndrho[idx2(i, j)] = r.d_ndrhorbardni_dxj__constxi(x, i, j, f);
                red.d_psi_t[idx2(i, j)] = r.d_psi_t_dxj(x, i, j);
                red.d_psi_rho[idx2(i, j)] = r.d_psi_rho_dxj(x, i, j);
                for k in 0..n {
                    red.d2_psi_t[(i * n + j) * n + k] = r.d2_psi_t_dxj_dxk(x, i, j, k);
                    red.d2_psi_rho[(i * n + j) * n + k] = r.d2_psi_rho_dxj_dxk(x, i, j, k);
                }
            }
        }
        CritState {
            model,
            x: x.to_vec(),
            n,
            red,
            t: f64::NAN,
            rhomolar: f64::NAN,
            tau: f64::NAN,
            delta: f64::NAN,
            ar: HelmholtzDerivs::default(),
            comp: Vec::new(),
            pair: vec![None; n * n],
        }
    }

    /// Upstream `update_DmolarT_direct` (and the imposed-gas `DmolarT`
    /// update, which reaches the same state): tau and delta are recomputed
    /// from (T, rho), exactly as upstream does.
    fn update(&mut self, rhomolar: f64, t: f64) -> Result<()> {
        if rhomolar < 0.0 {
            return Err(Error::Value(format!(
                "The molar density of {rhomolar:.6} mol/m3 is below the minimum of 0.000000 mol/m3"
            )));
        }
        if t < 0.0 {
            return Err(Error::Value(format!(
                "The temperature of {t:.6} K is below the minimum of 0.000000 K"
            )));
        }
        self.t = t;
        self.rhomolar = rhomolar;
        self.tau = self.red.tr / t;
        self.delta = rhomolar / self.red.rhor;
        let (tau, delta) = (self.tau, self.delta);
        self.ar = self.model.alphar_all(&self.x, tau, delta);
        self.comp = (0..self.n)
            .map(|i| self.model.component_alphar_all(i, tau, delta))
            .collect();
        self.pair = vec![None; self.n * self.n];
        for (i, j, f, dep) in self.model.excess_pairs() {
            let mut d = HelmholtzDerivs::default();
            d.add_scaled(&dep.all(tau, delta), f);
            self.pair[i * self.n + j] = Some(d);
            self.pair[j * self.n + i] = Some(d);
        }
        Ok(())
    }

    fn p(&self) -> f64 {
        self.rhomolar * self.model.gas_constant() * self.t * (1.0 + self.delta * self.ar.d10)
    }

    // --- mixture alphar and its residual composition derivatives ---

    fn a(&self, nd: u8, nt: u8) -> f64 {
        g(&self.ar, nd, nt)
    }
    /// `residual_helmholtz->d^(nd+nt+1) alphar / d x_i ...` (CS + Excess).
    fn rx(&self, i: usize, nd: u8, nt: u8) -> f64 {
        let mut ex = 0.0;
        for k in 0..self.n {
            if i != k {
                if let Some(p) = &self.pair[i * self.n + k] {
                    ex += self.x[k] * g(p, nd, nt);
                }
            }
        }
        g(&self.comp[i], nd, nt) + ex
    }
    /// The second composition derivative family (excess only).
    fn rxx(&self, i: usize, j: usize, nd: u8, nt: u8) -> f64 {
        if i == j {
            return 0.0;
        }
        self.pair[i * self.n + j]
            .as_ref()
            .map_or(0.0, |p| g(p, nd, nt))
    }

    // --- reducing shorthands ---

    fn psi_t(&self, i: usize) -> f64 {
        self.red.psi_t[i]
    }
    fn psi_rho(&self, i: usize) -> f64 {
        self.red.psi_rho[i]
    }
    fn d_psi_t(&self, i: usize, j: usize) -> f64 {
        self.red.d_psi_t[i * self.n + j]
    }
    fn d_psi_rho(&self, i: usize, j: usize) -> f64 {
        self.red.d_psi_rho[i * self.n + j]
    }
    fn d2_psi_t(&self, i: usize, j: usize, k: usize) -> f64 {
        self.red.d2_psi_t[(i * self.n + j) * self.n + k]
    }
    fn d2_psi_rho(&self, i: usize, j: usize, k: usize) -> f64 {
        self.red.d2_psi_rho[(i * self.n + j) * self.n + k]
    }

    // --- n d(tau)/dn_i, n d(delta)/dn_i and their derivatives ---

    fn ndtaudni(&self, i: usize) -> f64 {
        self.tau / self.red.tr * self.red.ndtr[i]
    }
    fn nddeltadni(&self, i: usize) -> f64 {
        self.delta - self.delta / self.red.rhor * self.red.ndrho[i]
    }
    fn d_ndtaudni_dTau(&self, i: usize) -> f64 {
        1.0 / self.red.tr * self.red.ndtr[i]
    }
    fn d_nddeltadni_dDelta(&self, i: usize) -> f64 {
        1.0 - 1.0 / self.red.rhor * self.red.ndrho[i]
    }
    fn d_ndtaudni_dxj(&self, i: usize, j: usize) -> f64 {
        let tr = self.red.tr;
        self.tau / tr
            * (self.red.d_ndtr[i * self.n + j] - 1.0 / tr * self.red.dtr[j] * self.red.ndtr[i])
    }
    fn d_nddeltadni_dxj(&self, i: usize, j: usize) -> f64 {
        let rhor = self.red.rhor;
        -self.delta / rhor
            * (self.red.d_ndrho[i * self.n + j] - 1.0 / rhor * self.red.drho[j] * self.red.ndrho[i])
    }
    fn d2_ndtaudni_dxj_dTau(&self, i: usize, j: usize) -> f64 {
        self.d_ndtaudni_dxj(i, j) / self.tau
    }
    fn d2_nddeltadni_dxj_dDelta(&self, i: usize, j: usize) -> f64 {
        self.d_nddeltadni_dxj(i, j) / self.delta
    }

    // --- n d(alphar)/dn_i and its tau/delta/x derivatives ---

    fn ndalphar_dni(&self, i: usize) -> f64 {
        let term1 = self.delta * self.a(1, 0) * (1.0 - 1.0 / self.red.rhor * self.red.ndrho[i]);
        let term2 = self.tau * self.a(0, 1) * (1.0 / self.red.tr) * self.red.ndtr[i];
        let mut s = 0.0;
        for k in 0..self.n {
            s += self.x[k] * self.rx(k, 0, 0);
        }
        let term3 = self.rx(i, 0, 0);
        term1 + term2 + term3 - s
    }
    fn d_ndalphardni_dTau(&self, i: usize) -> f64 {
        let term1 = self.delta * self.a(1, 1) * (1.0 - 1.0 / self.red.rhor * self.red.ndrho[i]);
        let term2 =
            (self.tau * self.a(0, 2) + self.a(0, 1)) * (1.0 / self.red.tr) * self.red.ndtr[i];
        let mut term3 = self.rx(i, 0, 1);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 0, 1);
        }
        term1 + term2 + term3
    }
    fn d_ndalphardni_dDelta(&self, i: usize) -> f64 {
        let term1 = (self.delta * self.a(2, 0) + self.a(1, 0))
            * (1.0 - 1.0 / self.red.rhor * self.red.ndrho[i]);
        let term2 = self.tau * self.a(1, 1) * (1.0 / self.red.tr) * self.red.ndtr[i];
        let mut term3 = self.rx(i, 1, 0);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 1, 0);
        }
        term1 + term2 + term3
    }
    fn d2_ndalphardni_dTau2(&self, i: usize) -> f64 {
        let term1 = self.delta * self.a(1, 2) * self.psi_rho(i);
        let term2 = (2.0 * self.a(0, 2) + self.tau * self.a(0, 3)) * self.psi_t(i);
        let mut term3 = self.rx(i, 0, 2);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 0, 2);
        }
        term1 + term2 + term3
    }
    fn d2_ndalphardni_dDelta2(&self, i: usize) -> f64 {
        let term1 = (2.0 * self.a(2, 0) + self.delta * self.a(3, 0)) * self.psi_rho(i);
        let term2 = self.tau * self.a(2, 1) * self.psi_t(i);
        let mut term3 = self.rx(i, 2, 0);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 2, 0);
        }
        term1 + term2 + term3
    }
    fn d2_ndalphardni_dDelta_dTau(&self, i: usize) -> f64 {
        let term1 = (self.a(1, 1) + self.delta * self.a(2, 1)) * self.psi_rho(i);
        let term2 = (self.tau * self.a(1, 2) + self.a(1, 1)) * self.psi_t(i);
        let mut term3 = self.rx(i, 1, 1);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 1, 1);
        }
        term1 + term2 + term3
    }
    fn d3_ndalphardni_dTau3(&self, i: usize) -> f64 {
        let term1 = self.delta * self.a(1, 3) * self.psi_rho(i);
        let term2 = (3.0 * self.a(0, 3) + self.tau * self.a(0, 4)) * self.psi_t(i);
        let mut term3 = self.rx(i, 0, 3);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 0, 3);
        }
        term1 + term2 + term3
    }
    fn d3_ndalphardni_dDelta3(&self, i: usize) -> f64 {
        let term1 = (3.0 * self.a(3, 0) + self.delta * self.a(4, 0)) * self.psi_rho(i);
        let term2 = self.tau * self.a(3, 1) * self.psi_t(i);
        let mut term3 = self.rx(i, 3, 0);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 3, 0);
        }
        term1 + term2 + term3
    }
    fn d3_ndalphardni_dDelta2_dTau(&self, i: usize) -> f64 {
        let term1 = (2.0 * self.a(2, 1) + self.delta * self.a(3, 1)) * self.psi_rho(i);
        let term2 = (self.tau * self.a(2, 2) + self.a(2, 1)) * self.psi_t(i);
        let mut term3 = self.rx(i, 2, 1);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 2, 1);
        }
        term1 + term2 + term3
    }
    fn d3_ndalphardni_dDelta_dTau2(&self, i: usize) -> f64 {
        let term1 = (self.a(1, 2) + self.delta * self.a(2, 2)) * self.psi_rho(i);
        let term2 = (self.tau * self.a(1, 3) + 2.0 * self.a(1, 2)) * self.psi_t(i);
        let mut term3 = self.rx(i, 1, 2);
        for k in 0..self.n {
            term3 -= self.x[k] * self.rx(k, 1, 2);
        }
        term1 + term2 + term3
    }

    fn d_ndalphardni_dxj(&self, i: usize, j: usize) -> f64 {
        let (tr, rhor) = (self.red.tr, self.red.rhor);
        let line1 = self.delta * self.rx(j, 1, 0) * (1.0 - 1.0 / rhor * self.red.ndrho[i]);
        let line3 = self.tau * self.rx(j, 0, 1) * (1.0 / tr) * self.red.ndtr[i];
        let line2 = -self.delta
            * self.a(1, 0)
            * (1.0 / rhor)
            * (self.red.d_ndrho[i * self.n + j]
                - 1.0 / rhor * self.red.drho[j] * self.red.ndrho[i]);
        let line4 = self.tau
            * self.a(0, 1)
            * (1.0 / tr)
            * (self.red.d_ndtr[i * self.n + j] - 1.0 / tr * self.red.dtr[j] * self.red.ndtr[i]);
        let mut s = 0.0;
        for k in 0..self.n {
            s += self.x[k] * self.rxx(j, k, 0, 0);
        }
        let line5 = self.rxx(i, j, 0, 0) - self.rx(j, 0, 0) - s;
        line1 + line2 + line3 + line4 + line5
    }
    fn d2_ndalphardni_dxj_dTau(&self, i: usize, j: usize) -> f64 {
        let term1 = self.delta * self.a(1, 1) * self.d_psi_rho(i, j);
        let term2 = self.delta * self.rx(j, 1, 1) * self.psi_rho(i);
        let term3 = (self.tau * self.a(0, 2) + self.a(0, 1)) * self.d_psi_t(i, j);
        let term4 = (self.tau * self.rx(j, 0, 2) + self.rx(j, 0, 1)) * self.psi_t(i);
        let mut term5 = self.rxx(i, j, 0, 1);
        for k in 0..self.n {
            term5 -= self.x[k] * self.rxx(k, j, 0, 1) + kd(k, j) * self.rx(k, 0, 1);
        }
        term1 + term2 + term3 + term4 + term5
    }
    fn d2_ndalphardni_dxj_dDelta(&self, i: usize, j: usize) -> f64 {
        let term1 = (self.a(1, 0) + self.delta * self.a(2, 0)) * self.d_psi_rho(i, j);
        let term2 = (self.rx(j, 1, 0) + self.delta * self.rx(j, 2, 0)) * self.psi_rho(i);
        let term3 = self.tau * self.a(1, 1) * self.d_psi_t(i, j);
        let term4 = self.tau * self.rx(j, 1, 1) * self.psi_t(i);
        let mut term5 = self.rxx(i, j, 1, 0);
        for k in 0..self.n {
            term5 -= self.x[k] * self.rxx(k, j, 1, 0) + kd(k, j) * self.rx(k, 1, 0);
        }
        term1 + term2 + term3 + term4 + term5
    }
    fn d2_ndalphardni_dxj_dxk(&self, i: usize, j: usize, k: usize) -> f64 {
        let term1 = self.delta
            * (self.rx(j, 1, 0) * self.d_psi_rho(i, k) + self.rx(k, 1, 0) * self.d_psi_rho(i, j));
        let term2 = self.delta * self.rxx(j, k, 1, 0) * self.psi_rho(i);
        let term3 = self.delta * self.a(1, 0) * self.d2_psi_rho(i, j, k);
        let term4 = self.tau
            * (self.rx(j, 0, 1) * self.d_psi_t(i, k) + self.rx(k, 0, 1) * self.d_psi_t(i, j));
        let term5 = self.tau * self.rxx(j, k, 0, 1) * self.psi_t(i);
        let term6 = self.tau * self.a(0, 1) * self.d2_psi_t(i, j, k);
        // d3alphardxidxjdxk is identically zero for the GERG-type models.
        let term7 = -2.0 * self.rxx(j, k, 0, 0);
        term1 + term2 + term3 + term4 + term5 + term6 + term7
    }
    fn d3_ndalphardni_dxj_dTau2(&self, i: usize, j: usize) -> f64 {
        let term1 = self.delta * self.a(1, 2) * self.d_psi_rho(i, j);
        let term2 = self.delta * self.rx(j, 1, 2) * self.psi_rho(i);
        let term3 = (self.tau * self.a(0, 3) + 2.0 * self.a(0, 2)) * self.d_psi_t(i, j);
        let term4 = (self.tau * self.rx(j, 0, 3) + 2.0 * self.rx(j, 0, 2)) * self.psi_t(i);
        let mut term5 = self.rxx(i, j, 0, 2);
        for k in 0..self.n {
            term5 -= self.x[k] * self.rxx(k, j, 0, 2) + kd(k, j) * self.rx(k, 0, 2);
        }
        term1 + term2 + term3 + term4 + term5
    }
    fn d3_ndalphardni_dxj_dDelta2(&self, i: usize, j: usize) -> f64 {
        let term1 = (2.0 * self.a(2, 0) + self.delta * self.a(3, 0)) * self.d_psi_rho(i, j);
        let term2 = (2.0 * self.rx(j, 2, 0) + self.delta * self.rx(j, 3, 0)) * self.psi_rho(i);
        let term3 = self.tau * self.a(2, 1) * self.d_psi_t(i, j);
        let term4 = self.tau * self.rx(j, 2, 1) * self.psi_t(i);
        let mut term5 = self.rxx(i, j, 2, 0);
        for k in 0..self.n {
            term5 -= self.x[k] * self.rxx(k, j, 2, 0) + kd(k, j) * self.rx(k, 2, 0);
        }
        term1 + term2 + term3 + term4 + term5
    }
    fn d3_ndalphardni_dxj_dDelta_dTau(&self, i: usize, j: usize) -> f64 {
        let term1 = (self.a(1, 1) + self.delta * self.a(2, 1)) * self.d_psi_rho(i, j);
        let term2 = (self.rx(j, 1, 1) + self.delta * self.rx(j, 2, 1)) * self.psi_rho(i);
        let term3 = (self.tau * self.a(1, 2) + self.a(1, 1)) * self.d_psi_t(i, j);
        let term4 = (self.tau * self.rx(j, 1, 2) + self.rx(j, 1, 1)) * self.psi_t(i);
        let mut term5 = self.rxx(i, j, 1, 1);
        for k in 0..self.n {
            term5 -= self.x[k] * self.rxx(k, j, 1, 1) + kd(k, j) * self.rx(k, 1, 1);
        }
        term1 + term2 + term3 + term4 + term5
    }
    fn d3_ndalphardni_dxj_dxk_dTau(&self, i: usize, j: usize, k: usize) -> f64 {
        let term1a = self.delta * self.rx(j, 1, 1) * self.d_psi_rho(i, k);
        let term1b = self.delta * self.rxx(j, k, 1, 1) * self.psi_rho(i);
        let term1c = self.delta * self.a(1, 1) * self.d2_psi_rho(i, j, k);
        let term1d = self.delta * self.rx(k, 1, 1) * self.d_psi_rho(i, j);
        let term1 = term1a + term1b + term1c + term1d;
        let term2a = (self.tau * self.rx(j, 0, 2) + self.rx(j, 0, 1)) * self.d_psi_t(i, k);
        let term2b = (self.tau * self.rxx(j, k, 0, 2) + self.rxx(j, k, 0, 1)) * self.psi_t(i);
        let term2c = (self.tau * self.a(0, 2) + self.a(0, 1)) * self.d2_psi_t(i, j, k);
        let term2d = (self.tau * self.rx(k, 0, 2) + self.rx(k, 0, 1)) * self.d_psi_t(i, j);
        let term2 = term2a + term2b + term2c + term2d;
        // d4alphar_dxi_dxj_dxk_dTau is identically zero.
        let term3 = -2.0 * self.rxx(j, k, 0, 1);
        term1 + term2 + term3
    }
    fn d3_ndalphardni_dxj_dxk_dDelta(&self, i: usize, j: usize, k: usize) -> f64 {
        let term1a = (self.delta * self.rx(j, 2, 0) + self.rx(j, 1, 0)) * self.d_psi_rho(i, k);
        let term1b = (self.delta * self.rxx(j, k, 2, 0) + self.rxx(j, k, 1, 0)) * self.psi_rho(i);
        let term1c = (self.delta * self.a(2, 0) + self.a(1, 0)) * self.d2_psi_rho(i, j, k);
        let term1d = (self.delta * self.rx(k, 2, 0) + self.rx(k, 1, 0)) * self.d_psi_rho(i, j);
        let term1 = term1a + term1b + term1c + term1d;
        let term2a = self.tau * self.rx(j, 1, 1) * self.d_psi_t(i, k);
        let term2b = self.tau * self.rxx(j, k, 1, 1) * self.psi_t(i);
        let term2c = self.tau * self.a(1, 1) * self.d2_psi_t(i, j, k);
        let term2d = self.tau * self.rx(k, 1, 1) * self.d_psi_t(i, j);
        let term2 = term2a + term2b + term2c + term2d;
        // d4alphar_dxi_dxj_dxk_dDelta is identically zero.
        let term3 = -2.0 * self.rxx(j, k, 1, 0);
        term1 + term2 + term3
    }

    // --- n d(n d alphar/dn_i)/dn_j and its derivatives ---

    fn nd_ndalphardni_dnj(&self, i: usize, j: usize) -> f64 {
        let line1 = self.d_ndalphardni_dDelta(i) * self.nddeltadni(j);
        let line2 = self.d_ndalphardni_dTau(i) * self.ndtaudni(j);
        let mut line3 = self.d_ndalphardni_dxj(i, j);
        for k in 0..self.n {
            line3 -= self.x[k] * self.d_ndalphardni_dxj(i, k);
        }
        line1 + line2 + line3
    }
    fn nd2nalphardnidnj(&self, i: usize, j: usize) -> f64 {
        self.ndalphar_dni(j) + self.nd_ndalphardni_dnj(i, j)
    }
    fn d_nd_ndalphardni_dnj_dTau(&self, i: usize, j: usize) -> f64 {
        let line1 = self.d2_ndalphardni_dDelta_dTau(i) * self.nddeltadni(j);
        let line2 = self.d2_ndalphardni_dTau2(i) * self.ndtaudni(j)
            + self.d_ndalphardni_dTau(i) * self.d_ndtaudni_dTau(j);
        let mut summer = 0.0;
        for k in 0..self.n {
            summer += self.x[k] * self.d2_ndalphardni_dxj_dTau(i, k);
        }
        let line3 = self.d2_ndalphardni_dxj_dTau(i, j) - summer;
        line1 + line2 + line3
    }
    fn d_nd_ndalphardni_dnj_dDelta(&self, i: usize, j: usize) -> f64 {
        let line1 = self.d2_ndalphardni_dDelta2(i) * self.nddeltadni(j)
            + self.d_ndalphardni_dDelta(i) * self.d_nddeltadni_dDelta(j);
        let line2 = self.d2_ndalphardni_dDelta_dTau(i) * self.ndtaudni(j);
        let mut summer = 0.0;
        for k in 0..self.n {
            summer += self.x[k] * self.d2_ndalphardni_dxj_dDelta(i, k);
        }
        let line3 = self.d2_ndalphardni_dxj_dDelta(i, j) - summer;
        line1 + line2 + line3
    }
    fn d_nd_ndalphardni_dnj_dxk(&self, i: usize, j: usize, k: usize) -> f64 {
        let line1 = self.d_ndalphardni_dDelta(i) * self.d_nddeltadni_dxj(j, k)
            + self.d2_ndalphardni_dxj_dDelta(i, k) * self.nddeltadni(j);
        let line2 = self.d_ndalphardni_dTau(i) * self.d_ndtaudni_dxj(j, k)
            + self.d2_ndalphardni_dxj_dTau(i, k) * self.ndtaudni(j);
        let mut line3 = self.d2_ndalphardni_dxj_dxk(i, j, k) - self.d_ndalphardni_dxj(i, k);
        for m in 0..self.n {
            line3 -= self.x[m] * self.d2_ndalphardni_dxj_dxk(i, m, k);
        }
        line1 + line2 + line3
    }
    fn d2_nd_ndalphardni_dnj_dTau2(&self, i: usize, j: usize) -> f64 {
        let line1 = self.d3_ndalphardni_dDelta_dTau2(i) * self.nddeltadni(j);
        let line2 = 2.0 * self.d2_ndalphardni_dTau2(i) * self.d_ndtaudni_dTau(j);
        let line3 = self.d3_ndalphardni_dTau3(i) * self.ndtaudni(j);
        let mut summer = 0.0;
        for k in 0..self.n {
            summer += self.x[k] * self.d3_ndalphardni_dxj_dTau2(i, k);
        }
        let line4 = self.d3_ndalphardni_dxj_dTau2(i, j) - summer;
        line1 + line2 + line3 + line4
    }
    fn d2_nd_ndalphardni_dnj_dDelta2(&self, i: usize, j: usize) -> f64 {
        let line1 = self.d3_ndalphardni_dDelta3(i) * self.nddeltadni(j);
        let line2 = 2.0 * self.d2_ndalphardni_dDelta2(i) * self.d_nddeltadni_dDelta(j);
        let line3 = self.d3_ndalphardni_dDelta2_dTau(i) * self.ndtaudni(j);
        let mut summer = 0.0;
        for k in 0..self.n {
            summer += self.x[k] * self.d3_ndalphardni_dxj_dDelta2(i, k);
        }
        let line4 = self.d3_ndalphardni_dxj_dDelta2(i, j) - summer;
        line1 + line2 + line3 + line4
    }
    fn d2_nd_ndalphardni_dnj_dDelta_dTau(&self, i: usize, j: usize) -> f64 {
        let line1 = self.d3_ndalphardni_dDelta2_dTau(i) * self.nddeltadni(j);
        let line2 = self.d2_ndalphardni_dDelta_dTau(i) * self.d_nddeltadni_dDelta(j);
        let line3 = self.d2_ndalphardni_dDelta_dTau(i) * self.d_ndtaudni_dTau(j);
        let line4 = self.d3_ndalphardni_dDelta_dTau2(i) * self.ndtaudni(j);
        let mut summer = 0.0;
        for k in 0..self.n {
            summer += self.x[k] * self.d3_ndalphardni_dxj_dDelta_dTau(i, k);
        }
        let line5 = self.d3_ndalphardni_dxj_dDelta_dTau(i, j) - summer;
        line1 + line2 + line3 + line4 + line5
    }
    fn d2_nd_ndalphardni_dnj_dxk_dTau(&self, i: usize, j: usize, k: usize) -> f64 {
        let line1 = self.d2_ndalphardni_dDelta_dTau(i) * self.d_nddeltadni_dxj(j, k)
            + self.d3_ndalphardni_dxj_dDelta_dTau(i, k) * self.nddeltadni(j);
        let line2 = self.d_ndalphardni_dTau(i) * self.d2_ndtaudni_dxj_dTau(j, k)
            + self.d2_ndalphardni_dTau2(i) * self.d_ndtaudni_dxj(j, k)
            + self.d2_ndalphardni_dxj_dTau(i, k) * self.d_ndtaudni_dTau(j)
            + self.d3_ndalphardni_dxj_dTau2(i, k) * self.ndtaudni(j);
        let mut line3 =
            self.d3_ndalphardni_dxj_dxk_dTau(i, j, k) - self.d2_ndalphardni_dxj_dTau(i, k);
        for m in 0..self.n {
            line3 -= self.x[m] * self.d3_ndalphardni_dxj_dxk_dTau(i, m, k);
        }
        line1 + line2 + line3
    }
    fn d2_nd_ndalphardni_dnj_dxk_dDelta(&self, i: usize, j: usize, k: usize) -> f64 {
        let line1 = self.d_ndalphardni_dDelta(i) * self.d2_nddeltadni_dxj_dDelta(j, k)
            + self.d2_ndalphardni_dDelta2(i) * self.d_nddeltadni_dxj(j, k)
            + self.d2_ndalphardni_dxj_dDelta(i, k) * self.d_nddeltadni_dDelta(j)
            + self.d3_ndalphardni_dxj_dDelta2(i, k) * self.nddeltadni(j);
        let line2 = self.d2_ndalphardni_dDelta_dTau(i) * self.d_ndtaudni_dxj(j, k)
            + self.d3_ndalphardni_dxj_dDelta_dTau(i, k) * self.ndtaudni(j);
        let mut line3 =
            self.d3_ndalphardni_dxj_dxk_dDelta(i, j, k) - self.d2_ndalphardni_dxj_dDelta(i, k);
        for m in 0..self.n {
            line3 -= self.x[m] * self.d3_ndalphardni_dxj_dxk_dDelta(i, m, k);
        }
        line1 + line2 + line3
    }

    // --- n d(ln f_i)/dn_j and its derivatives ---

    fn ndln_fugacity_i_dnj(&self, i: usize, j: usize) -> f64 {
        let s = if self.x[i] > f64::EPSILON {
            kd(i, j) / self.x[i]
        } else {
            0.0
        };
        s + self.nd2nalphardnidnj(i, j)
    }
    fn d_ndln_fugacity_i_dnj_dtau(&self, i: usize, j: usize) -> f64 {
        self.d_ndalphardni_dTau(j) + self.d_nd_ndalphardni_dnj_dTau(i, j)
    }
    fn d_ndln_fugacity_i_dnj_ddelta(&self, i: usize, j: usize) -> f64 {
        self.d_ndalphardni_dDelta(j) + self.d_nd_ndalphardni_dnj_dDelta(i, j)
    }
    fn d_ndln_fugacity_i_dnj_ddxk(&self, i: usize, j: usize, k: usize) -> f64 {
        let s = if self.x[i] > f64::EPSILON {
            -kd(i, j) * kd(i, k) / self.x[i].powi(2)
        } else {
            0.0
        };
        s + self.d_ndalphardni_dxj(j, k) + self.d_nd_ndalphardni_dnj_dxk(i, j, k)
    }
    fn d2_ndln_fugacity_i_dnj_dtau2(&self, i: usize, j: usize) -> f64 {
        self.d2_ndalphardni_dTau2(j) + self.d2_nd_ndalphardni_dnj_dTau2(i, j)
    }
    fn d2_ndln_fugacity_i_dnj_ddelta2(&self, i: usize, j: usize) -> f64 {
        self.d2_ndalphardni_dDelta2(j) + self.d2_nd_ndalphardni_dnj_dDelta2(i, j)
    }
    fn d2_ndln_fugacity_i_dnj_ddelta_dtau(&self, i: usize, j: usize) -> f64 {
        self.d2_ndalphardni_dDelta_dTau(j) + self.d2_nd_ndalphardni_dnj_dDelta_dTau(i, j)
    }
    fn d2_ndln_fugacity_i_dnj_dxk_dTau(&self, i: usize, j: usize, k: usize) -> f64 {
        self.d2_ndalphardni_dxj_dTau(j, k) + self.d2_nd_ndalphardni_dnj_dxk_dTau(i, j, k)
    }
    fn d2_ndln_fugacity_i_dnj_dxk_dDelta(&self, i: usize, j: usize, k: usize) -> f64 {
        self.d2_ndalphardni_dxj_dDelta(j, k) + self.d2_nd_ndalphardni_dnj_dxk_dDelta(i, j, k)
    }
    fn nd_ndln_fugacity_i_dnj_dnk(&self, i: usize, j: usize, k: usize) -> f64 {
        let mut sum = self.d_ndln_fugacity_i_dnj_dtau(i, j) * self.ndtaudni(k)
            + self.d_ndln_fugacity_i_dnj_ddelta(i, j) * self.nddeltadni(k)
            + self.d_ndln_fugacity_i_dnj_ddxk(i, j, k);
        for m in 0..self.n {
            sum -= self.x[m] * self.d_ndln_fugacity_i_dnj_ddxk(i, j, m);
        }
        sum
    }

    // --- the criticality matrices ---

    fn nAij(&self, i: usize, j: usize) -> f64 {
        self.ndln_fugacity_i_dnj(i, j)
    }
    fn n2Aijk(&self, i: usize, j: usize, k: usize) -> f64 {
        self.nd_ndln_fugacity_i_dnj_dnk(i, j, k) - self.nAij(i, j)
    }
    fn d_nAij_dX(&self, i: usize, j: usize, wrt: Wrt) -> f64 {
        match wrt {
            Wrt::Tau => self.d_ndln_fugacity_i_dnj_dtau(i, j),
            Wrt::Delta => self.d_ndln_fugacity_i_dnj_ddelta(i, j),
        }
    }
    fn d_n2Aijk_dX(&self, i: usize, j: usize, k: usize, wrt: Wrt) -> f64 {
        let mut summer = 0.0;
        match wrt {
            Wrt::Tau => {
                summer += self.d2_ndln_fugacity_i_dnj_dtau2(i, j) * self.ndtaudni(k);
                summer += self.d_ndln_fugacity_i_dnj_dtau(i, j) * self.d_ndtaudni_dTau(k);
                summer += self.d2_ndln_fugacity_i_dnj_ddelta_dtau(i, j) * self.nddeltadni(k);
                summer += self.d2_ndln_fugacity_i_dnj_dxk_dTau(i, j, k);
                for m in 0..self.n {
                    summer -= self.x[m] * self.d2_ndln_fugacity_i_dnj_dxk_dTau(i, j, m);
                }
            }
            Wrt::Delta => {
                summer += self.d2_ndln_fugacity_i_dnj_ddelta_dtau(i, j) * self.ndtaudni(k);
                summer += self.d2_ndln_fugacity_i_dnj_ddelta2(i, j) * self.nddeltadni(k);
                summer += self.d_ndln_fugacity_i_dnj_ddelta(i, j) * self.d_nddeltadni_dDelta(k);
                summer += self.d2_ndln_fugacity_i_dnj_dxk_dDelta(i, j, k);
                for m in 0..self.n {
                    summer -= self.x[m] * self.d2_ndln_fugacity_i_dnj_dxk_dDelta(i, j, m);
                }
            }
        }
        summer - self.d_nAij_dX(i, j, wrt)
    }

    fn lstar(&self) -> Mat {
        Mat::symmetric(self.n, |i, j| self.nAij(i, j))
    }
    fn dlstar_dx(&self, wrt: Wrt) -> Mat {
        Mat::symmetric(self.n, |i, j| self.d_nAij_dX(i, j, wrt))
    }
    fn d2lstar_dtau2(&self) -> Mat {
        Mat::symmetric(self.n, |i, j| self.d2_ndln_fugacity_i_dnj_dtau2(i, j))
    }
    /// `n2dLdni(j, k) = n2Aijk(j, k, i)`, symmetric in (j, k).
    fn n2dldni(&self, i: usize) -> Mat {
        Mat::symmetric(self.n, |j, k| self.n2Aijk(j, k, i))
    }
    fn mstar(&self, l: &Mat) -> Mat {
        let n = self.n;
        let mut m = l.clone();
        let adj_l = l.adjugate();
        for i in 0..n {
            m.set(n - 1, i, adj_l.trace_of_product(&self.n2dldni(i)));
        }
        m
    }
    fn dmstar_dx(&self, wrt: Wrt, l: &Mat, dl: &Mat) -> Mat {
        let n = self.n;
        let mut dm = dl.clone();
        let adj_l = l.adjugate();
        let d_adj_l = l.adjugate_derivative(dl);
        for i in 0..n {
            let n2 = self.n2dldni(i);
            let d_n2 = Mat::symmetric(n, |j, k| self.d_n2Aijk_dX(j, k, i, wrt));
            dm.set(
                n - 1,
                i,
                n2.trace_of_product(&d_adj_l) + adj_l.trace_of_product(&d_n2),
            );
        }
        dm
    }
}

#[derive(Clone, Copy)]
enum Wrt {
    Tau,
    Delta,
}

fn kd(i: usize, j: usize) -> f64 {
    if i == j { 1.0 } else { 0.0 }
}

/// Upstream `angle_difference`.
fn angle_difference(a1: f64, a2: f64) -> f64 {
    (a1 - a2 + PI) % (2.0 * PI) - PI
}

// ---------------------------------------------------------------------------
// The solvers' objective functions
// ---------------------------------------------------------------------------

/// Upstream `OneDimObjective`: `det(L*)` along tau at fixed delta.
struct OneDimObjective<'a, 'm> {
    st: &'a mut CritState<'m>,
    delta: f64,
}

impl Resid1D for OneDimObjective<'_, '_> {
    fn call(&mut self, tau: f64) -> f64 {
        let rhomolar = self.st.red.rhor * self.delta;
        let t = self.st.red.tr / tau;
        if self.st.update(rhomolar, t).is_err() {
            return f64::NAN;
        }
        self.st.lstar().determinant()
    }
    fn deriv(&mut self, _tau: f64) -> f64 {
        let adj = self.st.lstar().adjugate();
        adj.trace_of_product(&self.st.dlstar_dx(Wrt::Tau))
    }
    fn second_deriv(&mut self, _tau: f64) -> f64 {
        let l = self.st.lstar();
        let dl = self.st.dlstar_dx(Wrt::Tau);
        let d2l = self.st.d2lstar_dtau2();
        let adj = l.adjugate();
        let dadj = l.adjugate_derivative(&dl);
        dl.trace_of_product(&dadj) + adj.trace_of_product(&d2l)
    }
    fn third_deriv(&mut self, _tau: f64) -> f64 {
        unreachable!("Halley uses two derivatives")
    }
}

/// Upstream `Newton` (`FuncWrapper1DWithDeriv`): note the loop condition is
/// `iter < maxiter || |f| > ftol`, and convergence is on the step.
fn newton(
    mut call: impl FnMut(f64) -> f64,
    mut deriv: impl FnMut(f64) -> f64,
    x0: f64,
    ftol: f64,
    maxiter: i32,
) -> Result<f64> {
    let mut x = x0;
    let mut fval: f64 = 999.0;
    let mut iter = 0;
    while iter < maxiter || fval.abs() > ftol {
        fval = call(x);
        let dfdx = deriv(x);
        let dx = -fval / dfdx;
        if !fval.is_finite() {
            return Err(Error::Value(
                "Residual function in newton returned invalid number".into(),
            ));
        }
        x += dx;
        if (dx / x).abs() < 1e-11 {
            return Ok(x);
        }
        if iter > maxiter {
            return Err(Error::Solution(
                "Newton reached maximum number of iterations".into(),
            ));
        }
        iter += 1;
    }
    Ok(x)
}

/// Upstream `L0CurveTracer`: walk the `L1* = 0` contour and bracket sign
/// changes of `M1*`.
struct Tracer<'a, 'm> {
    st: &'a mut CritState<'m>,
    delta: f64,
    tau: f64,
    m1_last: f64,
    theta_last: f64,
    r_tau: f64,
    r_delta: f64,
    r_tau_tracer: f64,
    r_delta_tracer: f64,
    critical_points: Vec<CriticalPoint>,
    // State of the last `call`, as upstream keeps it in members.
    lstar: Option<Mat>,
    adj_lstar: Option<Mat>,
    dl_dtau: Option<Mat>,
    dl_ddelta: Option<Mat>,
}

impl Tracer<'_, '_> {
    fn tau_delta_at(&self, theta: f64) -> (f64, f64) {
        (
            self.tau + self.r_tau * theta.cos(),
            self.delta + self.r_delta * theta.sin(),
        )
    }
    /// `L1*` at the angle; upstream lets an invalid state propagate as an
    /// exception out of the solver, which here is a NaN the solvers refuse.
    fn call(&mut self, theta: f64) -> f64 {
        let (tau_new, delta_new) = self.tau_delta_at(theta);
        let rhomolar = self.st.red.rhor * delta_new;
        let t = self.st.red.tr / tau_new;
        if self.st.update(rhomolar, t).is_err() {
            return f64::NAN;
        }
        let l = self.st.lstar();
        self.adj_lstar = Some(l.adjugate());
        self.dl_dtau = Some(self.st.dlstar_dx(Wrt::Tau));
        self.dl_ddelta = Some(self.st.dlstar_dx(Wrt::Delta));
        let l1 = l.determinant();
        self.lstar = Some(l);
        l1
    }
    fn deriv(&self, theta: f64) -> f64 {
        let adj = self.adj_lstar.as_ref().expect("call precedes deriv");
        let dl1_dtau = adj.trace_of_product(self.dl_dtau.as_ref().expect("cached"));
        let dl1_ddelta = adj.trace_of_product(self.dl_ddelta.as_ref().expect("cached"));
        -self.r_tau * theta.sin() * dl1_dtau + self.r_delta * theta.cos() * dl1_ddelta
    }

    fn brent(&mut self, a: f64, b: f64) -> Result<f64> {
        brent(|th| self.call(th), a, b, f64::EPSILON, 1e-10, 100)
    }
    fn newton(&mut self, x0: f64) -> Result<f64> {
        // Newton calls call() then deriv() at each x; sharing &mut self
        // across two closures goes through a cell.
        let this = std::cell::RefCell::new(&mut *self);
        newton(
            |th| this.borrow_mut().call(th),
            |th| this.borrow().deriv(th),
            x0,
            1e-10,
            100,
        )
    }

    fn trace(&mut self) -> Result<()> {
        let mut theta = f64::NAN;
        for i in 0..300 {
            if i == 0 {
                self.r_tau = 0.001;
                self.r_delta = 0.001;
                theta = self.brent(0.0, PI)?;
            } else {
                self.r_tau = self.r_tau_tracer;
                self.r_delta = self.r_delta_tracer;
                match self.newton(self.theta_last) {
                    Ok(th) => theta = th,
                    Err(_) => {
                        if !self.critical_points.is_empty()
                            && self.delta
                                > 1.5 * self.critical_points[0].rhomolar / self.st.red.rhor
                        {
                            // A kink in the L1*=0 contour from poor
                            // low-temperature behaviour.
                            continue;
                        }
                    }
                }
                if angle_difference(theta, self.theta_last).abs() > PI / 2.0 {
                    if !self.critical_points.is_empty()
                        && self.delta > 1.2 * self.critical_points[0].rhomolar / self.st.red.rhor
                    {
                        continue;
                    }
                    theta = self.brent(self.theta_last - PI / 3.5, self.theta_last + PI / 3.5)?;
                }
            }

            // The second criticality condition, at the state the solver last
            // evaluated (upstream reads its members, not the returned angle).
            let lstar = self.lstar.as_ref().expect("a solver called call()");
            let m1 = self.st.mstar(lstar).determinant();
            let p_mpa = self.st.p() / 1e6;

            let (tau_new, delta_new) = self.tau_delta_at(theta);
            // Upstream `get_critical_is_terminated`.
            if p_mpa > 500.0 || delta_new > 5.0 || tau_new > 5.0 {
                break;
            }
            if i > 0 && m1 * self.m1_last < 0.0 {
                let rhomolar = self.st.red.rhor * (self.delta + delta_new) / 2.0;
                let t = self.st.red.tr / ((self.tau + tau_new) / 2.0);
                let crit = calc_critical_point(self.st, rhomolar, t)?;
                self.critical_points.push(crit);
            }
            self.tau = tau_new;
            self.delta = delta_new;
            self.m1_last = m1;
            self.theta_last = theta;
        }
        Ok(())
    }
}

/// Upstream `calc_critical_point`: `L1* = M1* = 0` by 2-D Newton in
/// (tau, delta) from (rho0, T0), tolerance 1e-10, 100 iterations.
fn calc_critical_point(st: &mut CritState<'_>, rho0: f64, t0: f64) -> Result<CriticalPoint> {
    let (tr, rhor) = (st.red.tr, st.red.rhor);
    let mut x = [tr / t0, rho0 / rhor];
    let tol = 1e-10;
    let maxiter = 100;
    let mut iter = 0;
    let mut error: f64 = 999.0;
    while iter == 0 || error.abs() > tol {
        // Resid::call
        st.update(x[1] * rhor, tr / x[0])?;
        let l = st.lstar();
        let m = st.mstar(&l);
        let f0 = [l.determinant(), m.determinant()];
        // Resid::Jacobian
        let adj_l = l.adjugate();
        let adj_m = m.adjugate();
        let dl_dtau = st.dlstar_dx(Wrt::Tau);
        let dl_ddelta = st.dlstar_dx(Wrt::Delta);
        let dm_dtau = st.dmstar_dx(Wrt::Tau, &l, &dl_dtau);
        let dm_ddelta = st.dmstar_dx(Wrt::Delta, &l, &dl_ddelta);
        let j = [
            [
                adj_l.trace_of_product(&dl_dtau),
                adj_l.trace_of_product(&dl_ddelta),
            ],
            [
                adj_m.trace_of_product(&dm_dtau),
                adj_m.trace_of_product(&dm_ddelta),
            ],
        ];
        // v = J^-1 (-r)
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let v = [
            (-f0[0] * j[1][1] + f0[1] * j[0][1]) / det,
            (-f0[1] * j[0][0] + f0[0] * j[1][0]) / det,
        ];
        let mut max_relchange: f64 = -1.0;
        for i in 0..2 {
            x[i] += v[i];
            let relchange = (v[i] / x[i]).abs();
            if x[i].abs() > 1e-16 && relchange > max_relchange {
                max_relchange = relchange;
            }
        }
        let max_abschange = v[0].abs().max(v[1].abs());
        if max_abschange < f64::EPSILON * 100.0 || max_relchange < 1e-12 {
            break;
        }
        error = (f0[0] * f0[0] + f0[1] * f0[1]).sqrt();
        if iter > maxiter {
            // Upstream poisons x[0] with _HUGE and keeps going; the next
            // update lands at T = 0 and refuses.
            x[0] = f64::INFINITY;
        }
        iter += 1;
    }
    let t = tr / x[0];
    let rhomolar = x[1] * rhor;
    let delta = rhomolar / rhor;
    let ar = st.model.alphar_all(&st.x, tr / t, delta);
    let p = rhomolar * st.model.gas_constant() * t * (1.0 + delta * ar.d10);
    Ok(CriticalPoint { t, p, rhomolar })
}

/// Upstream `calc_criticality_contour_values`: `(L1*, M1*)` = the
/// determinants of `L*` and `M*` at `(T, rhomolar)` — the two conditions a
/// critical point zeroes, exposed so they can be checked on their own.
pub fn criticality_contour_values(
    model: &MixtureModel,
    z: &[f64],
    t: f64,
    rhomolar: f64,
) -> Result<(f64, f64)> {
    let mut st = CritState::new(model, z);
    st.update(rhomolar, t)?;
    let l = st.lstar();
    let m = st.mstar(&l);
    Ok((l.determinant(), m.determinant()))
}

/// Upstream `_calc_all_critical_points(find_critical_points = true)`.
pub fn all_critical_points(model: &MixtureModel, z: &[f64]) -> Result<Vec<CriticalPoint>> {
    let mut st = CritState::new(model, z);
    // get_critical_point_starting_values: SPINODAL_MINIMUM_DELTA = 0.5.
    let delta0 = 0.5;
    let mut tau0 = 0.66;

    let tau_l0 = {
        let mut resid = OneDimObjective {
            st: &mut st,
            delta: delta0,
        };
        resid.call(tau0);
        // The deriv is taken at the state call() left behind, which the
        // bumps do not move — upstream's loop runs zero or three times.
        let mut bump_count = 0;
        while resid.deriv(tau0) > 0.0 && bump_count < 3 {
            tau0 *= 1.1;
            bump_count += 1;
        }
        halley(&mut resid, tau0, 1e-10, 100)?
    };

    let mut tracer = Tracer {
        st: &mut st,
        delta: delta0,
        tau: tau_l0,
        m1_last: f64::INFINITY,
        theta_last: f64::INFINITY,
        r_tau: 0.1,
        r_delta: 0.1,
        // get_critical_point_search_radii
        r_tau_tracer: 0.1,
        r_delta_tracer: 0.025,
        critical_points: Vec::new(),
        lstar: None,
        adj_lstar: None,
        dl_dtau: None,
        dl_ddelta: None,
    };
    tracer.trace()?;
    Ok(tracer.critical_points)
}

/// Upstream `calc_T_critical` / `calc_p_critical` / `calc_rhomolar_critical`
/// for a mixture: the critical point when the tracer finds exactly one.
pub fn critical_point(model: &MixtureModel, z: &[f64]) -> Result<CriticalPoint> {
    let pts = all_critical_points(model, z)?;
    if pts.len() == 1 {
        Ok(pts[0])
    } else {
        Err(Error::Value(format!(
            "critical point finding routine found {} critical points",
            pts.len()
        )))
    }
}
