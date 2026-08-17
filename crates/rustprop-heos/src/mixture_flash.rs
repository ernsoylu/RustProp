//! Mixture PT flash, slice 10d (upstream `FlashRoutines::PT_flash_mixtures`'s
//! single-phase path + the mixture branches of `solver_rho_Tp` /
//! `solver_rho_Tp_SRK`) and the homogeneous mixture properties.
//!
//! Deliberate 10d scope: upstream runs the Michelsen/Gernert stability tester
//! and a Wilson-seeded split BEFORE the single-phase lambda; until slice 10f
//! lands, this port assumes the feed is stable — a PT query inside the
//! two-phase dome returns the metastable single-phase root instead of the
//! upstream two-phase split. Golden coverage stays away from the dome.

use crate::alpha::{DerivsMemo, HelmholtzDerivs};
use crate::mixture::MixtureModel;
use crate::solvers::{self, Resid1D, brent, householder4};
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};

/// One solved single-phase mixture state.
#[derive(Debug, Clone, Copy)]
pub struct MixtureState {
    pub t: f64,
    pub p: f64,
    pub rhomolar: f64,
    /// Upstream sets `_Q = -1` on the single-phase PT path.
    pub q: f64,
    pub phase: Phase,
}

/// Upstream `SolverTPResid` over the mixture model: relative pressure
/// residual with density derivatives.
struct MixSolverTpResid<'a> {
    model: &'a MixtureModel,
    x: &'a [f64],
    t: f64,
    p: f64,
    tr: f64,
    rhor: f64,
    delta: f64,
    memo: DerivsMemo,
}

impl MixSolverTpResid<'_> {
    /// The mixture alphar matrix at (tr/t, delta), computed once per point.
    fn ar(&mut self, delta: f64) -> HelmholtzDerivs {
        let (model, x) = (self.model, self.x);
        self.memo
            .get_or_compute(self.tr / self.t, delta, |tau, delta| {
                model.alphar_all(x, tau, delta)
            })
    }
}

impl Resid1D for MixSolverTpResid<'_> {
    fn call(&mut self, rhomolar: f64) -> f64 {
        self.delta = rhomolar / self.rhor;
        // `calc_pressure` off the shared matrix — `self.tr`/`self.rhor` are
        // the same `reducing.tr(x)`/`rhormolar(x)` values `pressure`
        // recomputes, so tau/delta and the arithmetic are identical.
        let d = self.ar(self.delta);
        let peos = rhomolar * self.model.gas_constant() * self.t * (1.0 + self.delta * d.d10);
        (peos - self.p) / self.p
    }
    fn deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self.ar(self.delta);
        self.model.gas_constant()
            * self.t
            * (1.0 + 2.0 * self.delta * d.d10 + self.delta * self.delta * d.d20)
            / self.p
    }
    fn second_deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self.ar(self.delta);
        self.model.gas_constant() * self.t / self.rhor
            * (2.0 * d.d10 + 4.0 * self.delta * d.d20 + self.delta * self.delta * d.d30)
            / self.p
    }
    fn third_deriv(&mut self, _rhomolar: f64) -> f64 {
        let d = self.ar(self.delta);
        self.model.gas_constant() * self.t / (self.rhor * self.rhor)
            * (6.0 * d.d20 + 6.0 * self.delta * d.d30 + self.delta * self.delta * d.d40)
            / self.p
    }
}

impl MixtureModel {
    fn tau_delta(&self, x: &[f64], t: f64, rhomolar: f64) -> (f64, f64) {
        (
            self.reducing.tr(x) / t,
            rhomolar / self.reducing.rhormolar(x),
        )
    }

    fn all_pair(
        &self,
        x: &[f64],
        t: f64,
        rhomolar: f64,
    ) -> (HelmholtzDerivs, HelmholtzDerivs, f64, f64) {
        let tr = self.reducing.tr(x);
        let rhor = self.reducing.rhormolar(x);
        let (tau, delta) = (tr / t, rhomolar / rhor);
        (
            self.alphar_all(x, tau, delta),
            self.alpha0_all(x, tau, delta, tr, rhor),
            tau,
            delta,
        )
    }

    /// `calc_pressure`: p = rho R T (1 + delta dalphar_ddelta).
    pub fn pressure(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(x, t, rhomolar);
        let d = self.alphar_all(x, tau, delta);
        rhomolar * self.gas_constant() * t * (1.0 + delta * d.d10)
    }

    /// `dp/drho|T` (upstream `first_partial_deriv(iP, iDmolar, iT)`).
    pub fn dpdrho_t(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(x, t, rhomolar);
        let d = self.alphar_all(x, tau, delta);
        self.gas_constant() * t * (1.0 + 2.0 * delta * d.d10 + delta * delta * d.d20)
    }

    /// `d2p/drho2|T`.
    pub fn d2pdrho2_t(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let rhor = self.reducing.rhormolar(x);
        let (tau, delta) = self.tau_delta(x, t, rhomolar);
        let d = self.alphar_all(x, tau, delta);
        self.gas_constant() * t / rhor * (2.0 * d.d10 + 4.0 * delta * d.d20 + delta * delta * d.d30)
    }

    /// `calc_gibbsmolar_nocache`: g = R T (1 + a0 + ar + delta dar_ddelta).
    pub fn gibbsmolar_nocache(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, _, delta) = self.all_pair(x, t, rhomolar);
        self.gas_constant() * t * (1.0 + a0.d00 + ar.d00 + delta * ar.d10)
    }

    /// `calc_molar_mass`: mole-fraction-weighted component masses.
    pub fn molar_mass(&self, x: &[f64]) -> f64 {
        self.molar_masses()
            .iter()
            .zip(x)
            .map(|(m, xi)| m * xi)
            .sum()
    }

    /// `calc_hmolar` (homogeneous branch).
    pub fn hmolar(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, tau, delta) = self.all_pair(x, t, rhomolar);
        self.gas_constant() * t * (1.0 + tau * (a0.d01 + ar.d01) + delta * ar.d10)
    }

    /// `calc_smolar` (homogeneous branch).
    pub fn smolar(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, tau, _) = self.all_pair(x, t, rhomolar);
        self.gas_constant() * (tau * (a0.d01 + ar.d01) - a0.d00 - ar.d00)
    }

    /// `calc_umolar` (homogeneous branch).
    pub fn umolar(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, tau, _) = self.all_pair(x, t, rhomolar);
        self.gas_constant() * t * tau * (a0.d01 + ar.d01)
    }

    /// `calc_cvmolar`.
    pub fn cvmolar(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, tau, _) = self.all_pair(x, t, rhomolar);
        -self.gas_constant() * tau * tau * (a0.d02 + ar.d02)
    }

    /// `calc_cpmolar`.
    pub fn cpmolar(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, tau, delta) = self.all_pair(x, t, rhomolar);
        self.gas_constant()
            * (-tau * tau * (a0.d02 + ar.d02)
                + (1.0 + delta * ar.d10 - delta * tau * ar.d11).powi(2)
                    / (1.0 + 2.0 * delta * ar.d10 + delta * delta * ar.d20))
    }

    /// `calc_speed_sound`.
    pub fn speed_sound(&self, x: &[f64], t: f64, rhomolar: f64) -> f64 {
        let (ar, a0, tau, delta) = self.all_pair(x, t, rhomolar);
        (self.gas_constant() * t / self.molar_mass(x)
            * (1.0 + 2.0 * delta * ar.d10 + delta * delta * ar.d20
                - (1.0 + delta * ar.d10 - delta * tau * ar.d11).powi(2)
                    / (tau * tau * (a0.d02 + ar.d02))))
            .sqrt()
    }

    /// `solver_rho_Tp_SRK` — the mixture-capable cubic seed: quadratic
    /// van-der-Waals mixing with `k_ij = 0`, component constants from
    /// `EOS().reduce` and `EOS().acentric`.
    // Omega literals verbatim from upstream (over-long digit strings included).
    #[allow(clippy::excessive_precision)]
    pub fn solver_rho_tp_srk(&self, x: &[f64], t: f64, p: f64, phase: Phase) -> Result<f64> {
        let r_u = self.gas_constant();
        let mut a = 0.0;
        let mut b = 0.0;
        let consts = self.srk_component_consts();
        for (i, &(tci, pci, acentric_i)) in consts.iter().enumerate() {
            let m_i = 0.480 + 1.574 * acentric_i - 0.176 * acentric_i.powf(2.0);
            let b_i = 0.08664034999649577215890158147700 * r_u * tci / pci;
            b += x[i] * b_i;
            let a_i = 0.42748023335403414043900347952220 * (r_u * tci).powf(2.0) / pci
                * (1.0 + m_i * (1.0 - (t / tci).sqrt())).powf(2.0);
            for (j, &(tcj, pcj, acentric_j)) in consts.iter().enumerate() {
                let m_j = 0.480 + 1.574 * acentric_j - 0.176 * acentric_j.powf(2.0);
                let a_j = 0.42748023335403414043900347952220 * (r_u * tcj).powf(2.0) / pcj
                    * (1.0 + m_j * (1.0 - (t / tcj).sqrt())).powf(2.0);
                let k_ij = 0.0;
                a += x[i] * x[j] * (a_i * a_j).sqrt() * (1.0 - k_ij);
            }
        }

        let big_a = a * p / (r_u * t).powf(2.0);
        let big_b = b * p / (r_u * t);

        let (nsolns, z0, z1, z2) =
            solvers::solve_cubic(1.0, -1.0, big_a - big_b - big_b * big_b, -big_a * big_b);
        if nsolns == 1 {
            return Ok(p / (z0 * r_u * t));
        }
        let rhomolar0 = p / (z0 * r_u * t);
        let rhomolar1 = p / (z1 * r_u * t);
        let rhomolar2 = p / (z2 * r_u * t);
        if rhomolar0 > 0.0 && rhomolar1 <= 0.0 && rhomolar2 <= 0.0 {
            return Ok(rhomolar0);
        }
        if rhomolar0 <= 0.0 && rhomolar1 > 0.0 && rhomolar2 <= 0.0 {
            return Ok(rhomolar1);
        }
        if rhomolar0 <= 0.0 && rhomolar1 <= 0.0 && rhomolar2 > 0.0 {
            return Ok(rhomolar2);
        }
        match phase {
            Phase::Liquid | Phase::SupercriticalLiquid => {
                Ok(rhomolar0.max(rhomolar1).max(rhomolar2))
            }
            Phase::Gas | Phase::SupercriticalGas | Phase::Supercritical => {
                Ok(rhomolar0.min(rhomolar1).min(rhomolar2))
            }
            _ => Err(Error::Value("Bad phase to solver_rho_Tp_SRK".into())),
        }
    }

    /// `solver_rho_Tp`, mixture branches (no guess supplied by any caller in
    /// the PT path; the pure-fluid ancillary branches never apply here).
    pub fn solver_rho_tp(&self, x: &[f64], t: f64, p: f64, phase: Phase) -> Result<f64> {
        let tr = self.reducing.tr(x);
        let rhor = self.reducing.rhormolar(x);
        let mut resid = MixSolverTpResid {
            model: self,
            x,
            t,
            p,
            tr,
            rhor,
            delta: f64::NAN,
            memo: DerivsMemo::default(),
        };

        let mut rhomolar_guess = self.solver_rho_tp_srk(x, t, p, phase)?;
        if matches!(
            phase,
            Phase::Gas | Phase::SupercriticalGas | Phase::Supercritical
        ) {
            if rhomolar_guess < 0.0 || !rhomolar_guess.is_finite() {
                rhomolar_guess = p / (self.gas_constant() * t);
            }
        } else if phase == Phase::Liquid {
            // Mixture liquid: 4th-order Householder from a very high density.
            return householder4(&mut resid, 3.0 * rhor, 1e-8, 100);
        }

        let attempt = householder4(&mut resid, rhomolar_guess, 1e-8, 20).and_then(|rho| {
            if !rho.is_finite() || rho < 0.0 {
                return Err(Error::Value("invalid density".into()));
            }
            if phase == Phase::Gas {
                let dpdrho = self.dpdrho_t(x, t, rho);
                let d2pdrho2 = self.d2pdrho2_t(x, t, rho);
                if dpdrho < 0.0 || d2pdrho2 > 0.0 {
                    // Retry with a tiny density to land on the gas branch
                    return householder4(&mut resid, 1e-6, 1e-8, 100);
                }
            }
            Ok(rho)
        });
        match attempt {
            Ok(rho) => Ok(rho),
            Err(e) => {
                if matches!(phase, Phase::Supercritical | Phase::SupercriticalGas) {
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
                    "solver_rho_Tp was unable to find a solution for T={t:.10e}, p={p:.10e}, with guess value {rhomolar_guess:.10e} with error: {e}"
                )))
            }
        }
    }

    /// The single-phase lambda of `PT_flash_mixtures`: SRK roots for both
    /// phases decide which HEOS branches to solve; both valid → lowest Gibbs.
    pub fn pt_flash(&self, x: &[f64], t: f64, p: f64) -> Result<MixtureState> {
        let rho_srk_gas = self.solver_rho_tp_srk(x, t, p, Phase::Gas);
        let rho_srk_liq = self.solver_rho_tp_srk(x, t, p, Phase::Liquid);
        let gas_ok = matches!(rho_srk_gas, Ok(r) if r > 0.0 && r.is_finite());
        let liq_ok = matches!(rho_srk_liq, Ok(r) if r > 0.0 && r.is_finite());

        let rho = if gas_ok && liq_ok {
            let rho_gas = self.solver_rho_tp(x, t, p, Phase::Gas).unwrap_or(-1.0);
            let rho_liq = self.solver_rho_tp(x, t, p, Phase::Liquid).unwrap_or(-1.0);
            if rho_gas > 0.0 && rho_liq > 0.0 {
                let g_gas = self.gibbsmolar_nocache(x, t, rho_gas);
                let g_liq = self.gibbsmolar_nocache(x, t, rho_liq);
                if g_liq <= g_gas { rho_liq } else { rho_gas }
            } else if rho_gas > 0.0 {
                rho_gas
            } else if rho_liq > 0.0 {
                rho_liq
            } else {
                return Err(Error::Value(
                    "Unable to obtain either HEOS density root in PT_flash_mixtures".into(),
                ));
            }
        } else {
            let primary = if gas_ok {
                Phase::Gas
            } else if liq_ok {
                Phase::Liquid
            } else {
                Phase::Gas
            };
            let fallback = if primary == Phase::Gas {
                Phase::Liquid
            } else {
                Phase::Gas
            };
            match self.solver_rho_tp(x, t, p, primary) {
                Ok(rho) => rho,
                Err(_) => self.solver_rho_tp(x, t, p, fallback)?,
            }
        };
        let phase = if rho < self.reducing.rhormolar(x) {
            Phase::Gas
        } else {
            Phase::Liquid
        };
        Ok(MixtureState {
            t,
            p,
            rhomolar: rho,
            q: -1.0,
            phase,
        })
    }
}
