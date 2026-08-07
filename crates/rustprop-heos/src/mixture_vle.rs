//! Mixture VLE, slice 10e: the blind QT/PQ flash chain (upstream
//! `FlashRoutines::QT_flash` / `PQ_flash` mixture branches) — Wilson/
//! preconditioner guesses, `successive_substitution`,
//! `newton_raphson_saturation`, and the `MixtureDerivatives` layer they
//! consume. The phase-envelope fast paths are PropsSI-dead upstream and
//! unported.

// Function names mirror upstream's `__const` suffix convention; loops keep
// upstream's index shapes for reviewability against the C++.
#![allow(non_snake_case)]
#![allow(clippy::needless_range_loop)]

use crate::alpha::HelmholtzDerivs;
use crate::mixture::MixtureModel;
use crate::mixture::XnFlag;
use crate::solvers::{Resid1D, brent, householder4, secant};
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};

// ---------------------------------------------------------------------------
// A live single-phase saturation state (upstream SatL / SatV instances)
// ---------------------------------------------------------------------------

/// One phase of the VLE problem: composition + (T, rhomolar) with every
/// derivative the fugacity layer needs cached at that state.
pub struct SatState<'m> {
    model: &'m MixtureModel,
    pub x: Vec<f64>,
    pub t: f64,
    pub rhomolar: f64,
    // Cached at (x, t, rhomolar):
    tr: f64,
    rhor: f64,
    tau: f64,
    delta: f64,
    /// Mixture residual derivatives.
    ar: HelmholtzDerivs,
    /// Per-component pure residual derivatives at the mixture (tau, delta).
    comp_ar: Vec<HelmholtzDerivs>,
    /// Per-excess-pair `(i, j, F, derivs)` at the mixture (tau, delta).
    pair_ar: Vec<(usize, usize, f64, HelmholtzDerivs)>,
}

impl<'m> SatState<'m> {
    pub fn new(model: &'m MixtureModel, x: Vec<f64>) -> Self {
        SatState {
            model,
            x,
            t: f64::NAN,
            rhomolar: f64::NAN,
            tr: f64::NAN,
            rhor: f64::NAN,
            tau: f64::NAN,
            delta: f64::NAN,
            ar: HelmholtzDerivs::default(),
            comp_ar: Vec::new(),
            pair_ar: Vec::new(),
        }
    }

    pub fn set_mole_fractions(&mut self, x: &[f64]) {
        self.x.clear();
        self.x.extend_from_slice(x);
    }

    /// Recompute every cached derivative at (T, rhomolar) — upstream
    /// `update_DmolarT_direct` plus the lazy derivative caches.
    fn set_state(&mut self, t: f64, rhomolar: f64) {
        self.t = t;
        self.rhomolar = rhomolar;
        self.tr = self.model.reducing.tr(&self.x);
        self.rhor = self.model.reducing.rhormolar(&self.x);
        self.tau = self.tr / t;
        self.delta = rhomolar / self.rhor;
        self.ar = self.model.alphar_all(&self.x, self.tau, self.delta);
        let n = self.model.n_components();
        self.comp_ar.clear();
        for i in 0..n {
            self.comp_ar
                .push(self.model.component_alphar_all(i, self.tau, self.delta));
        }
        self.pair_ar.clear();
        for (i, j, f, dep) in self.model.excess_pairs() {
            self.pair_ar.push((i, j, f, dep.all(self.tau, self.delta)));
        }
    }

    /// Upstream `update_TP_guessrho`: solve rho at (T, p) from the guess with
    /// the imposed phase's stability retries, then set the state.
    pub fn update_tp_guessrho(
        &mut self,
        t: f64,
        p: f64,
        rho_guess: f64,
        phase: Phase,
    ) -> Result<()> {
        let rho = self.solver_rho_tp_guessed(t, p, rho_guess, phase)?;
        self.set_state(t, rho);
        Ok(())
    }

    /// The main path of `solver_rho_Tp` with a provided guess: Householder4
    /// plus the imposed-phase mechanical-stability retries; the supercritical
    /// Brent catch (the guess<0 branches never run — a guess is always given).
    fn solver_rho_tp_guessed(
        &mut self,
        t: f64,
        p: f64,
        rho_guess: f64,
        phase: Phase,
    ) -> Result<f64> {
        struct GuessResid<'a, 'm> {
            state: &'a SatState<'m>,
            t: f64,
            p: f64,
            rhor: f64,
            delta: f64,
        }
        impl Resid1D for GuessResid<'_, '_> {
            fn call(&mut self, rhomolar: f64) -> f64 {
                self.delta = rhomolar / self.rhor;
                let peos = self.state.model.pressure(&self.state.x, self.t, rhomolar);
                (peos - self.p) / self.p
            }
            fn deriv(&mut self, _rhomolar: f64) -> f64 {
                let d = self.state.model.alphar_all(
                    &self.state.x,
                    self.state.model.reducing.tr(&self.state.x) / self.t,
                    self.delta,
                );
                self.state.model.gas_constant()
                    * self.t
                    * (1.0 + 2.0 * self.delta * d.d10 + self.delta * self.delta * d.d20)
                    / self.p
            }
            fn second_deriv(&mut self, _rhomolar: f64) -> f64 {
                let d = self.state.model.alphar_all(
                    &self.state.x,
                    self.state.model.reducing.tr(&self.state.x) / self.t,
                    self.delta,
                );
                self.state.model.gas_constant() * self.t / self.rhor
                    * (2.0 * d.d10 + 4.0 * self.delta * d.d20 + self.delta * self.delta * d.d30)
                    / self.p
            }
            fn third_deriv(&mut self, _rhomolar: f64) -> f64 {
                let d = self.state.model.alphar_all(
                    &self.state.x,
                    self.state.model.reducing.tr(&self.state.x) / self.t,
                    self.delta,
                );
                self.state.model.gas_constant() * self.t / (self.rhor * self.rhor)
                    * (6.0 * d.d20 + 6.0 * self.delta * d.d30 + self.delta * self.delta * d.d40)
                    / self.p
            }
        }
        let rhor = self.model.reducing.rhormolar(&self.x);
        let mut resid = GuessResid {
            state: self,
            t,
            p,
            rhor,
            delta: f64::NAN,
        };
        let attempt = householder4(&mut resid, rho_guess, 1e-8, 20).and_then(|rho| {
            if !rho.is_finite() || rho < 0.0 {
                return Err(Error::Value("invalid density".into()));
            }
            match phase {
                Phase::Liquid => {
                    let dpdrho = self.model.dpdrho_t(&self.x, t, rho);
                    let d2pdrho2 = self.model.d2pdrho2_t(&self.x, t, rho);
                    if dpdrho < 0.0 || d2pdrho2 < 0.0 {
                        let mut resid = GuessResid {
                            state: self,
                            t,
                            p,
                            rhor,
                            delta: f64::NAN,
                        };
                        return householder4(&mut resid, 3.0 * rhor, 1e-8, 100);
                    }
                    Ok(rho)
                }
                Phase::Gas => {
                    let dpdrho = self.model.dpdrho_t(&self.x, t, rho);
                    let d2pdrho2 = self.model.d2pdrho2_t(&self.x, t, rho);
                    if dpdrho < 0.0 || d2pdrho2 > 0.0 {
                        let mut resid = GuessResid {
                            state: self,
                            t,
                            p,
                            rhor,
                            delta: f64::NAN,
                        };
                        return householder4(&mut resid, 1e-6, 1e-8, 100);
                    }
                    Ok(rho)
                }
                _ => Ok(rho),
            }
        });
        match attempt {
            Ok(rho) => Ok(rho),
            Err(e) => {
                if matches!(phase, Phase::Supercritical | Phase::SupercriticalGas) {
                    let mut resid = GuessResid {
                        state: self,
                        t,
                        p,
                        rhor,
                        delta: f64::NAN,
                    };
                    return brent(
                        |v| resid.call(v),
                        1e-10,
                        3.0 * rhor,
                        f64::EPSILON,
                        1e-8,
                        100,
                    );
                }
                Err(Error::Value(format!(
                    "solver_rho_Tp was unable to find a solution for T={t:.10e}, p={p:.10e}, with guess value {rho_guess:.10e} with error: {e}"
                )))
            }
        }
    }

    /// EOS pressure of the current state.
    pub fn p(&self) -> f64 {
        self.rhomolar * self.model.gas_constant() * self.t * (1.0 + self.delta * self.ar.d10)
    }

    pub fn hmolar(&self) -> f64 {
        self.model.hmolar(&self.x, self.t, self.rhomolar)
    }

    pub fn smolar(&self) -> f64 {
        self.model.smolar(&self.x, self.t, self.rhomolar)
    }

    // -- residual_helmholtz composition derivatives (CS + Excess) ----------

    /// F[i][k] * (departure ik derivs); zero when no pair stores (i,k).
    fn f_pair(&self, i: usize, k: usize) -> Option<(f64, &HelmholtzDerivs)> {
        self.pair_ar
            .iter()
            .find(|(a, b, _, _)| (*a == i && *b == k) || (*a == k && *b == i))
            .map(|(_, _, f, d)| (*f, d))
    }

    /// `pick` selects which derivative of alphar_ik rides the formula
    /// (d00 for dalphar_dxi, d10 for d2alphar_dxi_dDelta, d01 for dTau).
    fn excess_dxi(&self, i: usize, flag: XnFlag, pick: fn(&HelmholtzDerivs) -> f64) -> f64 {
        let n = self.x.len();
        match flag {
            XnFlag::Independent => {
                let mut summer = 0.0;
                for k in 0..n {
                    if k != i {
                        if let Some((f, d)) = self.f_pair(i, k) {
                            summer += self.x[k] * f * pick(d);
                        }
                    }
                }
                summer
            }
            XnFlag::Dependent => {
                if i == n - 1 {
                    return 0.0;
                }
                let fin_arin = self
                    .f_pair(i, n - 1)
                    .map(|(f, d)| f * pick(d))
                    .unwrap_or(0.0);
                let mut dar_dxi = (1.0 - 2.0 * self.x[i]) * fin_arin;
                for k in 0..(n - 1) {
                    if i == k {
                        continue;
                    }
                    let fik_arik = self.f_pair(i, k).map(|(f, d)| f * pick(d)).unwrap_or(0.0);
                    let fkn_arkn = self
                        .f_pair(k, n - 1)
                        .map(|(f, d)| f * pick(d))
                        .unwrap_or(0.0);
                    dar_dxi += self.x[k] * (fik_arik - fin_arin - fkn_arkn);
                }
                dar_dxi
            }
        }
    }

    /// Upstream `ResidualHelmholtz::dalphar_dxi` = CS + Excess.
    fn dalphar_dxi(&self, i: usize, flag: XnFlag) -> f64 {
        let n = self.x.len();
        let cs = match flag {
            XnFlag::Independent => self.comp_ar[i].d00,
            XnFlag::Dependent => {
                if i == n - 1 {
                    0.0
                } else {
                    self.comp_ar[i].d00 - self.comp_ar[n - 1].d00
                }
            }
        };
        cs + self.excess_dxi(i, flag, |d| d.d00)
    }

    fn d2alphar_dxi_ddelta(&self, i: usize, flag: XnFlag) -> f64 {
        let n = self.x.len();
        let cs = match flag {
            XnFlag::Independent => self.comp_ar[i].d10,
            XnFlag::Dependent => {
                if i == n - 1 {
                    0.0
                } else {
                    self.comp_ar[i].d10 - self.comp_ar[n - 1].d10
                }
            }
        };
        cs + self.excess_dxi(i, flag, |d| d.d10)
    }

    fn d2alphar_dxi_dtau(&self, i: usize, flag: XnFlag) -> f64 {
        let n = self.x.len();
        let cs = match flag {
            XnFlag::Independent => self.comp_ar[i].d01,
            XnFlag::Dependent => {
                if i == n - 1 {
                    0.0
                } else {
                    self.comp_ar[i].d01 - self.comp_ar[n - 1].d01
                }
            }
        };
        cs + self.excess_dxi(i, flag, |d| d.d01)
    }

    /// Upstream `ResidualHelmholtz::d2alphardxidxj` — CS part is zero.
    fn d2alphardxidxj(&self, i: usize, j: usize, flag: XnFlag) -> f64 {
        let n = self.x.len();
        match flag {
            XnFlag::Independent => {
                if i != j {
                    self.f_pair(i, j).map(|(f, d)| f * d.d00).unwrap_or(0.0)
                } else {
                    0.0
                }
            }
            XnFlag::Dependent => {
                if i == n - 1 || j == n - 1 {
                    return 0.0;
                }
                let fin_arin = self.f_pair(i, n - 1).map(|(f, d)| f * d.d00).unwrap_or(0.0);
                if i == j {
                    return -2.0 * fin_arin;
                }
                let fij_arij = self.f_pair(i, j).map(|(f, d)| f * d.d00).unwrap_or(0.0);
                let fjn_arjn = self.f_pair(j, n - 1).map(|(f, d)| f * d.d00).unwrap_or(0.0);
                fij_arij - fin_arin - fjn_arjn
            }
        }
    }

    // -- MixtureDerivatives ------------------------------------------------

    /// GERG 7.47 sub-part: `ndalphar_dni__constT_V_nj`.
    fn ndalphar_dni__const_t_v_nj(&self, i: usize, flag: XnFlag) -> f64 {
        let red = &self.model.reducing;
        let term1 = self.delta
            * self.ar.d10
            * (1.0 - 1.0 / self.rhor * red.ndrhorbardni__constnj(&self.x, i, flag));
        let term2 =
            self.tau * self.ar.d01 * (1.0 / self.tr) * red.ndtrdni__constnj(&self.x, i, flag);
        let mut s = 0.0;
        let mut kmax = self.x.len();
        if flag == XnFlag::Dependent {
            kmax -= 1;
        }
        for k in 0..kmax {
            s += self.x[k] * self.dalphar_dxi(k, flag);
        }
        let term3 = self.dalphar_dxi(i, flag);
        term1 + term2 + term3 - s
    }

    /// `ln_fugacity_coefficient` = alphar + ndalphar_dni - ln(1 + delta ar_delta).
    pub fn ln_fugacity_coefficient(&self, i: usize, flag: XnFlag) -> f64 {
        self.ar.d00 + self.ndalphar_dni__const_t_v_nj(i, flag)
            - (1.0 + self.delta * self.ar.d10).ln()
    }

    /// `fugacity_i` = x_i rho R T exp(alphar + ndalphar_dni).
    pub fn fugacity_i(&self, i: usize, flag: XnFlag) -> f64 {
        self.x[i]
            * self.rhomolar
            * self.model.gas_constant()
            * self.t
            * (self.ar.d00 + self.ndalphar_dni__const_t_v_nj(i, flag)).exp()
    }

    /// `d_ndalphardni_dTau`.
    fn d_ndalphardni_dtau(&self, i: usize, flag: XnFlag) -> f64 {
        let red = &self.model.reducing;
        let term1 = self.delta
            * self.ar.d11
            * (1.0 - 1.0 / self.rhor * red.ndrhorbardni__constnj(&self.x, i, flag));
        let term2 = (self.tau * self.ar.d02 + self.ar.d01)
            * (1.0 / self.tr)
            * red.ndtrdni__constnj(&self.x, i, flag);
        let mut term3 = self.d2alphar_dxi_dtau(i, flag);
        let mut kmax = self.x.len();
        if flag == XnFlag::Dependent {
            kmax -= 1;
        }
        for k in 0..kmax {
            term3 -= self.x[k] * self.d2alphar_dxi_dtau(k, flag);
        }
        term1 + term2 + term3
    }

    /// `d_ndalphardni_dDelta`.
    fn d_ndalphardni_ddelta(&self, i: usize, flag: XnFlag) -> f64 {
        let red = &self.model.reducing;
        let term1 = (self.delta * self.ar.d20 + self.ar.d10)
            * (1.0 - 1.0 / self.rhor * red.ndrhorbardni__constnj(&self.x, i, flag));
        let term2 =
            self.tau * self.ar.d11 * (1.0 / self.tr) * red.ndtrdni__constnj(&self.x, i, flag);
        let mut term3 = self.d2alphar_dxi_ddelta(i, flag);
        let mut kmax = self.x.len();
        if flag == XnFlag::Dependent {
            kmax -= 1;
        }
        for k in 0..kmax {
            term3 -= self.x[k] * self.d2alphar_dxi_ddelta(k, flag);
        }
        term1 + term2 + term3
    }

    /// `d2nalphar_dni_dT` = -tau/T (dalphar_dTau + d_ndalphardni_dTau).
    fn d2nalphar_dni_dt(&self, i: usize, flag: XnFlag) -> f64 {
        -self.tau / self.t * (self.ar.d01 + self.d_ndalphardni_dtau(i, flag))
    }

    /// `dpdT__constV_n`.
    fn dpdt__const_v_n(&self) -> f64 {
        self.rhomolar
            * self.model.gas_constant()
            * (1.0 + self.delta * self.ar.d10 - self.delta * self.tau * self.ar.d11)
    }

    /// `ndpdV__constT_n`.
    fn ndpdv__const_t_n(&self) -> f64 {
        -self.rhomolar
            * self.rhomolar
            * self.model.gas_constant()
            * self.t
            * (1.0 + 2.0 * self.delta * self.ar.d10 + self.delta * self.delta * self.ar.d20)
    }

    /// `ndpdni__constT_V_nj` (GERG 7.63/7.64).
    fn ndpdni__const_t_v_nj(&self, i: usize, flag: XnFlag) -> f64 {
        let red = &self.model.reducing;
        let ndrhorbar = red.ndrhorbardni__constnj(&self.x, i, flag);
        let ndtr = red.ndtrdni__constnj(&self.x, i, flag);
        let mut summer = 0.0;
        let mut kmax = self.x.len();
        if flag == XnFlag::Dependent {
            kmax -= 1;
        }
        for k in 0..kmax {
            summer += self.x[k] * self.d2alphar_dxi_ddelta(k, flag);
        }
        let nd2alphar_dni_ddelta = self.delta * self.ar.d20 * (1.0 - 1.0 / self.rhor * ndrhorbar)
            + self.tau * self.ar.d11 / self.tr * ndtr
            + self.d2alphar_dxi_ddelta(i, flag)
            - summer;
        self.rhomolar
            * self.model.gas_constant()
            * self.t
            * (1.0
                + self.delta * self.ar.d10 * (2.0 - 1.0 / self.rhor * ndrhorbar)
                + self.delta * nd2alphar_dni_ddelta)
    }

    /// `partial_molar_volume`.
    fn partial_molar_volume(&self, i: usize, flag: XnFlag) -> f64 {
        -self.ndpdni__const_t_v_nj(i, flag) / self.ndpdv__const_t_n()
    }

    /// `dln_fugacity_coefficient_dT__constp_n`.
    pub fn dln_fugacity_coefficient_dt__constp_n(&self, i: usize, flag: XnFlag) -> f64 {
        let r_u = self.model.gas_constant();
        self.d2nalphar_dni_dt(i, flag) + 1.0 / self.t
            - self.partial_molar_volume(i, flag) / (r_u * self.t) * self.dpdt__const_v_n()
    }

    /// `dln_fugacity_coefficient_dp__constT_n` (GERG 7.30).
    pub fn dln_fugacity_coefficient_dp__const_t_n(&self, i: usize, flag: XnFlag) -> f64 {
        self.partial_molar_volume(i, flag) / (self.model.gas_constant() * self.t) - 1.0 / self.p()
    }

    /// `dln_fugacity_i_dp__constT_n` = above + 1/p.
    pub fn dln_fugacity_i_dp__const_t_n(&self, i: usize, flag: XnFlag) -> f64 {
        self.dln_fugacity_coefficient_dp__const_t_n(i, flag) + 1.0 / self.p()
    }

    /// Gernert 3.121: `ddelta_dxj__constT_V_xi`.
    fn ddelta_dxj__const_t_v_xi(&self, j: usize, flag: XnFlag) -> f64 {
        -self.delta / self.rhor * self.model.reducing.drhormolardxi__constxj(&self.x, j, flag)
    }

    /// Gernert 3.122: `dtau_dxj__constT_V_xi`.
    fn dtau_dxj__const_t_v_xi(&self, j: usize, flag: XnFlag) -> f64 {
        1.0 / self.t * self.model.reducing.dtrdxi__constxj(&self.x, j, flag)
    }

    /// Gernert 3.134: `d_dalpharddelta_dxj__constT_V_xi`.
    fn d_dalpharddelta_dxj__const_t_v_xi(&self, j: usize, flag: XnFlag) -> f64 {
        self.ar.d20 * self.ddelta_dxj__const_t_v_xi(j, flag)
            + self.ar.d11 * self.dtau_dxj__const_t_v_xi(j, flag)
            + self.d2alphar_dxi_ddelta(j, flag)
    }

    /// Gernert 3.119: `dalphar_dxj__constT_V_xi`.
    fn dalphar_dxj__const_t_v_xi(&self, j: usize, flag: XnFlag) -> f64 {
        self.ar.d10 * self.ddelta_dxj__const_t_v_xi(j, flag)
            + self.ar.d01 * self.dtau_dxj__const_t_v_xi(j, flag)
            + self.dalphar_dxi(j, flag)
    }

    /// Gernert 3.130: `dpdxj__constT_V_xi`.
    fn dpdxj__const_t_v_xi(&self, j: usize, flag: XnFlag) -> f64 {
        self.rhomolar
            * self.model.gas_constant()
            * self.t
            * (self.ddelta_dxj__const_t_v_xi(j, flag) * self.ar.d10
                + self.delta * self.d_dalpharddelta_dxj__const_t_v_xi(j, flag))
    }

    /// `d_ndalphardni_dxj__constdelta_tau_xi`.
    fn d_ndalphardni_dxj__constdelta_tau_xi(&self, i: usize, j: usize, flag: XnFlag) -> f64 {
        let red = &self.model.reducing;
        let line1 = self.delta
            * self.d2alphar_dxi_ddelta(j, flag)
            * (1.0 - 1.0 / self.rhor * red.ndrhorbardni__constnj(&self.x, i, flag));
        let line3 = self.tau
            * self.d2alphar_dxi_dtau(j, flag)
            * (1.0 / self.tr)
            * red.ndtrdni__constnj(&self.x, i, flag);
        let line2 = -self.delta
            * self.ar.d10
            * (1.0 / self.rhor)
            * (red.d_ndrhorbardni_dxj__constxi(&self.x, i, j, flag)
                - 1.0 / self.rhor
                    * red.drhormolardxi__constxj(&self.x, j, flag)
                    * red.ndrhorbardni__constnj(&self.x, i, flag));
        let line4 = self.tau
            * self.ar.d01
            * (1.0 / self.tr)
            * (red.d_ndtrdni_dxj__constxi(&self.x, i, j, flag)
                - 1.0 / self.tr
                    * red.dtrdxi__constxj(&self.x, j, flag)
                    * red.ndtrdni__constnj(&self.x, i, flag));
        let mut s = 0.0;
        let mut kmax = self.x.len();
        if flag == XnFlag::Dependent {
            kmax -= 1;
        }
        for k in 0..kmax {
            s += self.x[k] * self.d2alphardxidxj(j, k, flag);
        }
        let line5 = self.d2alphardxidxj(i, j, flag) - self.dalphar_dxi(j, flag) - s;
        line1 + line2 + line3 + line4 + line5
    }

    /// Gernert 3.118: `d_ndalphardni_dxj__constT_V_xi`.
    fn d_ndalphardni_dxj__const_t_v_xi(&self, i: usize, j: usize, flag: XnFlag) -> f64 {
        self.d_ndalphardni_dxj__constdelta_tau_xi(i, j, flag)
            + self.ddelta_dxj__const_t_v_xi(j, flag) * self.d_ndalphardni_ddelta(i, flag)
            + self.dtau_dxj__const_t_v_xi(j, flag) * self.d_ndalphardni_dtau(i, flag)
    }

    /// `d2nalphar_dxj_dni__constT_V`.
    fn d2nalphar_dxj_dni__const_t_v(&self, j: usize, i: usize, flag: XnFlag) -> f64 {
        self.d_ndalphardni_dxj__const_t_v_xi(i, j, flag) + self.dalphar_dxj__const_t_v_xi(j, flag)
    }

    /// Gernert 3.115: `dln_fugacity_coefficient_dxj__constT_p_xi`.
    fn dln_fugacity_coefficient_dxj__const_t_p_xi(&self, i: usize, j: usize, flag: XnFlag) -> f64 {
        let r_u = self.model.gas_constant();
        self.d2nalphar_dxj_dni__const_t_v(j, i, flag)
            - self.partial_molar_volume(i, flag) / (r_u * self.t)
                * self.dpdxj__const_t_v_xi(j, flag)
    }

    /// `dln_fugacity_dxj__constT_p_xi`.
    pub fn dln_fugacity_dxj__const_t_p_xi(&self, i: usize, j: usize, flag: XnFlag) -> f64 {
        let n = self.x.len();
        let mut val = self.dln_fugacity_coefficient_dxj__const_t_p_xi(i, j, flag);
        if i == n - 1 {
            val += -1.0 / self.x[n - 1];
        } else if i == j {
            val += 1.0 / self.x[j];
        }
        val
    }
}

// ---------------------------------------------------------------------------
// Guess chain (upstream VLERoutines.h inlines)
// ---------------------------------------------------------------------------

/// `Wilson_lnK_factor`: ln K_i = ln(pci/p) + 5.373 (1+omega_i)(1 - Tci/T).
fn wilson_ln_k_factor(model: &MixtureModel, t: f64, p: f64, i: usize) -> f64 {
    let pci = model.crit_p()[i];
    let tci = model.crit_t()[i];
    let omegai = model.acentric(i);
    (pci / p).ln() + 5.373 * (1.0 + omegai) * (1.0 - tci / t)
}

/// Which variable the successive-substitution chain iterates on
/// (upstream `sstype_enum`).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum SsType {
    ImposedT,
    ImposedP,
}

/// `saturation_preconditioner`: interpolate ln(p) vs T between the
/// mole-fraction-weighted triple and critical points.
pub fn saturation_preconditioner(
    model: &MixtureModel,
    input_value: f64,
    input_type: SsType,
    z: &[f64],
) -> f64 {
    let mut tcrit = 0.0;
    let mut pcrit = 0.0;
    let mut ttriple = 0.0;
    let mut ptriple = 0.0;
    for i in 0..z.len() {
        tcrit += model.crit_t()[i] * z[i];
        pcrit += model.crit_p()[i] * z[i];
        ttriple += model.triple_t()[i] * z[i];
        ptriple += model.triple_p()[i] * z[i];
    }
    match input_type {
        SsType::ImposedT => ((pcrit / ptriple).ln() / (tcrit - ttriple) * (input_value - ttriple)
            + ptriple.ln())
        .exp(),
        SsType::ImposedP => {
            1.0 / (1.0 / tcrit
                - (1.0 / ttriple - 1.0 / tcrit) / (pcrit / ptriple).ln()
                    * (input_value / pcrit).ln())
        }
    }
}

/// `saturation_Wilson`: explicit solution for beta in {0,1} with T imposed,
/// otherwise a bounded Brent (Secant fallback) on the Rachford-Rice residual.
/// Also fills K with the Wilson factors of the returned state.
pub fn saturation_wilson(
    model: &MixtureModel,
    beta: f64,
    input_value: f64,
    input_type: SsType,
    z: &[f64],
    guess: f64,
    k: &mut [f64],
) -> Result<f64> {
    if input_type == SsType::ImposedT && (beta.abs() < 1e-12 || (beta - 1.0).abs() < 1e-12) {
        let beta0 = beta.abs() < 1e-12;
        let mut out = 0.0;
        for i in 0..z.len() {
            let pci = model.crit_p()[i];
            let tci = model.crit_t()[i];
            let omegai = model.acentric(i);
            if beta0 {
                out += z[i] * pci * (5.373 * (1.0 + omegai) * (1.0 - tci / input_value)).exp();
            } else {
                out += z[i] / (pci * (5.373 * (1.0 + omegai) * (1.0 - tci / input_value)).exp());
            }
        }
        if !beta0 {
            out = 1.0 / out;
        }
        for i in 0..z.len() {
            let pci = model.crit_p()[i];
            let tci = model.crit_t()[i];
            let omegai = model.acentric(i);
            k[i] = pci / out * (5.373 * (1.0 + omegai) * (1.0 - tci / input_value)).exp();
        }
        Ok(out)
    } else {
        // Rachford-Rice residual over the Wilson K-factors.
        let resid = |iterate: f64, k: &mut [f64]| -> f64 {
            let (t, p) = match input_type {
                SsType::ImposedT => (input_value, iterate),
                SsType::ImposedP => (iterate, input_value),
            };
            let mut summer = 0.0;
            for i in 0..z.len() {
                k[i] = wilson_ln_k_factor(model, t, p, i).exp();
                summer += z[i] * (k[i] - 1.0) / (1.0 - beta + beta * k[i]);
            }
            summer
        };
        let (a, b) = match input_type {
            SsType::ImposedT => (1.0, 1e9),
            SsType::ImposedP => (50.0, 10000.0),
        };
        let first = brent(|v| resid(v, k), a, b, 1e-10, 1e-10, 100);
        let out = match first {
            Ok(v) => v,
            Err(_) => {
                if !guess.is_finite() || guess < 0.0 {
                    return Err(Error::Value(
                        "saturation_p_Wilson failed to get good output value".into(),
                    ));
                }
                secant(|v| resid(v, k), guess, 0.001, 1e-10, 100)?
            }
        };
        if !out.is_finite() {
            return Err(Error::Value(
                "saturation_p_Wilson failed to get good output value".into(),
            ));
        }
        Ok(out)
    }
}

/// `x_and_y_from_K`.
pub fn x_and_y_from_k(beta: f64, k: &[f64], z: &[f64], x: &mut [f64], y: &mut [f64]) {
    for i in 0..k.len() {
        let denominator = 1.0 - beta + beta * k[i];
        x[i] = z[i] / denominator;
        y[i] = k[i] * z[i] / denominator;
    }
}

fn normalize_vector(x: &mut [f64]) {
    let sum: f64 = x.iter().sum();
    for v in x.iter_mut() {
        *v /= sum;
    }
}

// ---------------------------------------------------------------------------
// successive_substitution
// ---------------------------------------------------------------------------

/// Outputs of the successive-substitution stage (upstream `mixture_VLE_IO`).
pub struct MixtureVleIo {
    pub t: f64,
    pub p: f64,
    pub rhomolar_liq: f64,
    pub rhomolar_vap: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// Upstream `SaturationSolvers::successive_substitution`.
#[allow(clippy::too_many_arguments)]
pub fn successive_substitution<'m>(
    model: &'m MixtureModel,
    beta: f64,
    mut t: f64,
    mut p: f64,
    z: &[f64],
    k: &mut [f64],
    sstype: SsType,
    nstep_max: i32,
    satl: &mut SatState<'m>,
    satv: &mut SatState<'m>,
) -> Result<MixtureVleIo> {
    let mut iter = 1;
    let n = z.len();
    let mut ln_phi_liq = vec![0.0; n];
    let mut ln_phi_vap = vec![0.0; n];

    let mut x = vec![0.0; n];
    let mut y = vec![0.0; n];
    x_and_y_from_k(beta, k, z, &mut x, &mut y);
    normalize_vector(&mut x);
    normalize_vector(&mut y);
    satl.set_mole_fractions(&x);
    satv.set_mole_fractions(&y);

    let mut rhomolar_liq = model.solver_rho_tp_srk(&x, t, p, Phase::Liquid)?;
    let rhomolar_vap = model.solver_rho_tp_srk(&y, t, p, Phase::Gas)?;

    // Peneloux volume translation for the liquid seed
    // (Horstmann doi:10.1016/j.fluid.2004.11.002; upstream hardcodes this R).
    let mut summer_c = 0.0;
    let v_srk = 1.0 / rhomolar_liq;
    for i in 0..n {
        let tc = model.crit_t()[i];
        let pc = model.crit_p()[i];
        let rhomolarc = model.crit_rhomolar()[i];
        let r = 8.3144598;
        summer_c += z[i] * (0.40768 * r * tc / pc * (0.29441 - pc / (rhomolarc * r * tc)));
    }
    rhomolar_liq = 1.0 / (v_srk - summer_c);
    satl.update_tp_guessrho(t, p, rhomolar_liq, Phase::Liquid)?;
    satv.update_tp_guessrho(t, p, rhomolar_vap, Phase::Gas)?;

    loop {
        satl.update_tp_guessrho(t, p, satl.rhomolar, Phase::Liquid)?;
        satv.update_tp_guessrho(t, p, satv.rhomolar, Phase::Gas)?;

        let mut f = 0.0;
        let mut df = 0.0;
        let flag = XnFlag::Independent;
        for i in 0..n {
            ln_phi_liq[i] = satl.ln_fugacity_coefficient(i, flag);
            ln_phi_vap[i] = satv.ln_fugacity_coefficient(i, flag);

            let (deriv_liq, deriv_vap) = match sstype {
                SsType::ImposedP => (
                    satl.dln_fugacity_coefficient_dt__constp_n(i, flag),
                    satv.dln_fugacity_coefficient_dt__constp_n(i, flag),
                ),
                SsType::ImposedT => (
                    satl.dln_fugacity_coefficient_dp__const_t_n(i, flag),
                    satv.dln_fugacity_coefficient_dp__const_t_n(i, flag),
                ),
            };

            k[i] = (ln_phi_liq[i] - ln_phi_vap[i]).exp();
            f += z[i] * (k[i] - 1.0) / (1.0 - beta + beta * k[i]);
            let dfdk = k[i] * z[i] / (1.0 - beta + beta * k[i]).powi(2);
            df += dfdk * (deriv_liq - deriv_vap);
        }

        let change = if df.abs() <= 1e-14 {
            if f.abs() <= 1e-12 {
                -f
            } else {
                return Err(Error::Value(format!(
                    "df very small (df = {df:e}) in successive_substitution but f is not converged (f = {f:e} > 1e-12)."
                )));
            }
        } else {
            -f / df
        };

        match sstype {
            SsType::ImposedP => {
                t += change;
            }
            SsType::ImposedT => {
                let omega = if change.abs() > 0.05 * p { 0.1 } else { 1.0 };
                p += omega * change;
            }
        }

        x_and_y_from_k(beta, k, z, &mut x, &mut y);
        normalize_vector(&mut x);
        normalize_vector(&mut y);
        satl.set_mole_fractions(&x);
        satv.set_mole_fractions(&y);

        iter += 1;
        if iter > 50 {
            return Err(Error::Value(
                "saturation_p was unable to reach a solution within 50 iterations".into(),
            ));
        }
        if !(f.abs() > 1e-12 && iter < nstep_max) {
            break;
        }
    }

    satl.update_tp_guessrho(t, p, satl.rhomolar, Phase::Liquid)?;
    satv.update_tp_guessrho(t, p, satv.rhomolar, Phase::Gas)?;

    Ok(MixtureVleIo {
        p: satl.p(),
        t: satl.t,
        rhomolar_liq: satl.rhomolar,
        rhomolar_vap: satv.rhomolar,
        x,
        y,
    })
}

// ---------------------------------------------------------------------------
// newton_raphson_saturation
// ---------------------------------------------------------------------------

/// Which variable stays fixed in the NR solve.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum NrImposed {
    TImposed,
    PImposed,
}

/// In/out block for `newton_raphson_saturation` (upstream
/// `newton_raphson_saturation_options`, RHOV_IMPOSED unused by QT/PQ).
pub struct NrSaturationOptions {
    pub imposed_variable: NrImposed,
    pub bubble_point: bool,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub rhomolar_liq: f64,
    pub rhomolar_vap: f64,
    pub t: f64,
    pub p: f64,
    pub nstep_max: i32,
}

/// Solve the dense linear system J v = -r in place (partial-pivot Gaussian
/// elimination; upstream uses Eigen's column-pivoted QR — same solution, and
/// the surrounding Newton converges both to the same fixed point).
fn solve_linear(j: &mut [Vec<f64>], r: &[f64]) -> Result<Vec<f64>> {
    let n = r.len();
    let mut aug: Vec<Vec<f64>> = (0..n)
        .map(|row| {
            let mut v = j[row].clone();
            v.push(-r[row]);
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
        if aug[piv][col] == 0.0 {
            return Err(Error::Value("singular Jacobian in NR saturation".into()));
        }
        aug.swap(col, piv);
        for row in (col + 1)..n {
            let factor = aug[row][col] / aug[col][col];
            for kcol in col..=n {
                aug[row][kcol] -= factor * aug[col][kcol];
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
    Ok(v)
}

/// Upstream `newton_raphson_saturation::call` + `build_arrays`
/// (T_IMPOSED / P_IMPOSED branches).
pub fn newton_raphson_saturation<'m>(
    model: &'m MixtureModel,
    io: &mut NrSaturationOptions,
    satl: &mut SatState<'m>,
    satv: &mut SatState<'m>,
) -> Result<()> {
    let n = io.x.len();
    let mut iter = 0;
    let mut t = io.t;
    let mut p = io.p;
    let mut rhomolar_liq = io.rhomolar_liq;
    let mut rhomolar_vap = io.rhomolar_vap;
    let mut x = io.x.clone();
    let mut y = io.y.clone();
    let mut error_rms;
    let mut min_rel_change;
    let _ = model;

    loop {
        // build_arrays
        satl.set_mole_fractions(&x);
        satv.set_mole_fractions(&y);
        satl.update_tp_guessrho(t, p, rhomolar_liq, Phase::Liquid)?;
        rhomolar_liq = satl.rhomolar;
        satv.update_tp_guessrho(t, p, rhomolar_vap, Phase::Gas)?;
        rhomolar_vap = satv.rhomolar;

        let p_liq = satl.p();
        let p_vap = satv.p();
        p = 0.5 * (p_liq + p_vap);

        let flag = XnFlag::Dependent;
        let mut r = vec![0.0; n];
        let mut j: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        for i in 0..n {
            let ln_f_liq = satl.fugacity_i(i, flag).ln();
            let ln_f_vap = satv.fugacity_i(i, flag).ln();
            r[i] = ln_f_liq - ln_f_vap;
            for jj in 0..(n - 1) {
                j[i][jj] = if io.bubble_point {
                    -satv.dln_fugacity_dxj__const_t_p_xi(i, jj, flag)
                } else {
                    satl.dln_fugacity_dxj__const_t_p_xi(i, jj, flag)
                };
            }
            j[i][n - 1] = match io.imposed_variable {
                NrImposed::PImposed => {
                    satl.dln_fugacity_coefficient_dt__constp_n(i, flag)
                        - satv.dln_fugacity_coefficient_dt__constp_n(i, flag)
                }
                NrImposed::TImposed => {
                    satl.dln_fugacity_i_dp__const_t_n(i, flag)
                        - satv.dln_fugacity_i_dp__const_t_n(i, flag)
                }
            };
        }
        error_rms = r.iter().map(|v| v * v).sum::<f64>().sqrt();

        // Newton step
        let v = solve_linear(&mut j, &r)?;
        let mut err_rel = vec![0.0; n];
        if io.bubble_point {
            for i in 0..(n - 1) {
                err_rel[i] = v[i] / y[i];
                y[i] += v[i];
            }
            y[n - 1] = 1.0 - y[..n - 1].iter().sum::<f64>();
        } else {
            for i in 0..(n - 1) {
                err_rel[i] = v[i] / x[i];
                x[i] += v[i];
            }
            x[n - 1] = 1.0 - x[..n - 1].iter().sum::<f64>();
        }
        match io.imposed_variable {
            NrImposed::PImposed => {
                t += v[n - 1];
                err_rel[n - 1] = v[n - 1] / t;
            }
            NrImposed::TImposed => {
                p += v[n - 1];
                err_rel[n - 1] = v[n - 1] / p;
            }
        }
        min_rel_change = err_rel.iter().fold(f64::INFINITY, |m, e| m.min(e.abs()));
        iter += 1;

        if iter == io.nstep_max {
            return Err(Error::Value(format!(
                "newton_raphson_saturation::call reached max number of iterations [{}]",
                io.nstep_max
            )));
        }
        if !(error_rms > 1e-7 && min_rel_change > 1000.0 * f64::EPSILON && iter < io.nstep_max) {
            break;
        }
    }

    io.p = p;
    io.x = x;
    io.y = y;
    io.t = t;
    io.rhomolar_liq = rhomolar_liq;
    io.rhomolar_vap = rhomolar_vap;
    Ok(())
}

// ---------------------------------------------------------------------------
// QT / PQ flash drivers
// ---------------------------------------------------------------------------

/// One converged two-phase mixture state: bulk quantities plus both phases.
pub struct MixtureTwoPhase {
    pub t: f64,
    pub p: f64,
    pub rhomolar: f64,
    pub q: f64,
    pub x_liq: Vec<f64>,
    pub y_vap: Vec<f64>,
    pub rhomolar_liq: f64,
    pub rhomolar_vap: f64,
    /// Phase enthalpies/entropies for the lever rule, evaluated at each
    /// phase's own composition and density.
    pub hmolar_liq: f64,
    pub hmolar_vap: f64,
    pub smolar_liq: f64,
    pub smolar_vap: f64,
}

impl MixtureModel {
    fn qx_flash(
        &self,
        q: f64,
        imposed_value: f64,
        sstype: SsType,
        z: &[f64],
    ) -> Result<MixtureTwoPhase> {
        let n = z.len();
        let mut k = vec![0.0; n];

        let (nstep_max_ss, guess) = match sstype {
            SsType::ImposedT => (
                20,
                saturation_preconditioner(self, imposed_value, SsType::ImposedT, z),
            ),
            SsType::ImposedP => (
                10,
                saturation_preconditioner(self, imposed_value, SsType::ImposedP, z),
            ),
        };
        let refined = saturation_wilson(self, q, imposed_value, sstype, z, guess, &mut k)?;

        let (t_ss, p_ss) = match sstype {
            SsType::ImposedT => (imposed_value, refined),
            SsType::ImposedP => (refined, imposed_value),
        };

        let mut satl = SatState::new(self, z.to_vec());
        let mut satv = SatState::new(self, z.to_vec());
        let options = successive_substitution(
            self,
            q,
            t_ss,
            p_ss,
            z,
            &mut k,
            sstype,
            nstep_max_ss,
            &mut satl,
            &mut satv,
        )?;

        let mut io = NrSaturationOptions {
            imposed_variable: match sstype {
                SsType::ImposedT => NrImposed::TImposed,
                SsType::ImposedP => NrImposed::PImposed,
            },
            bubble_point: q < 0.5,
            x: options.x,
            y: options.y,
            rhomolar_liq: options.rhomolar_liq,
            rhomolar_vap: options.rhomolar_vap,
            t: options.t,
            p: options.p,
            nstep_max: 30,
        };
        newton_raphson_saturation(self, &mut io, &mut satl, &mut satv)?;

        // Load the outputs exactly as upstream: p from SatV's state, T from
        // SatL's, lever-rule density from the phase states.
        Ok(MixtureTwoPhase {
            t: satl.t,
            p: satv.p(),
            rhomolar: 1.0 / (q / satv.rhomolar + (1.0 - q) / satl.rhomolar),
            q,
            hmolar_liq: satl.hmolar(),
            hmolar_vap: satv.hmolar(),
            smolar_liq: satl.smolar(),
            smolar_vap: satv.smolar(),
            x_liq: satl.x.clone(),
            y_vap: satv.x.clone(),
            rhomolar_liq: satl.rhomolar,
            rhomolar_vap: satv.rhomolar,
        })
    }

    /// Upstream `QT_flash` mixture branch (blind path).
    pub fn qt_flash(&self, q: f64, t: f64, z: &[f64]) -> Result<MixtureTwoPhase> {
        self.qx_flash(q, t, SsType::ImposedT, z)
    }

    /// Upstream `PQ_flash` mixture branch (blind path).
    pub fn pq_flash(&self, p: f64, q: f64, z: &[f64]) -> Result<MixtureTwoPhase> {
        self.qx_flash(q, p, SsType::ImposedP, z)
    }
}

impl MixtureTwoPhase {
    /// Lever-rule bulk enthalpy (upstream two-phase `calc_hmolar`).
    pub fn hmolar(&self) -> f64 {
        if self.q.abs() < f64::EPSILON {
            self.hmolar_liq
        } else if (self.q - 1.0).abs() < f64::EPSILON {
            self.hmolar_vap
        } else {
            self.q * self.hmolar_vap + (1.0 - self.q) * self.hmolar_liq
        }
    }

    /// Lever-rule bulk entropy.
    pub fn smolar(&self) -> f64 {
        if self.q.abs() < f64::EPSILON {
            self.smolar_liq
        } else if (self.q - 1.0).abs() < f64::EPSILON {
            self.smolar_vap
        } else {
            self.q * self.smolar_vap + (1.0 - self.q) * self.smolar_liq
        }
    }
}
