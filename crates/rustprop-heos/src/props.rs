//! Single-phase properties at (T, rhomolar) from the Helmholtz derivative
//! matrix (PLAN.md 4.2) — formulas ported operation-for-operation from
//! `HelmholtzEOSMixtureBackend.cpp` @ v8.0.0 (`calc_pressure`,
//! `calc_hmolar_nocache`, `calc_smolar_nocache`, `calc_umolar_nocache`,
//! `calc_cvmolar`, `calc_cpmolar`, `calc_speed_sound`,
//! `calc_gibbsmolar_nocache`).
//!
//! For pure fluids `gas_constant()` is the fluid-specific EOS value
//! (`components[0].gas_constant()`), NOT the CODATA constant — verified
//! against the wheel (Water: 8.314371357587 J/mol/K).

use crate::alpha::HelmholtzEos;
use crate::derivs::DerivEos;
use rustprop_core::Result;

impl HelmholtzEos {
    /// Upstream `HelmholtzEOSMixtureBackend::calc_alpha0_deriv_nocache`'s
    /// closing validity gate at the state `(t, rhomolar)` — see
    /// [`crate::derivs::check_alpha0`] for the mechanism. Callers name the
    /// derivative upstream fetches for the output they are about to serve.
    ///
    /// The `pow(rhor/rhomolarc, nDelta)` / `pow(Tr/Tc, nTau)` rescale upstream
    /// applies before the check is identically 1 for a multiparameter EOS
    /// (upstream's own comment: it exists for cubics, whose `tau*`/`delta*`
    /// are shifted).
    pub fn check_alpha0_deriv(
        &self,
        n_tau: u32,
        n_delta: u32,
        t: f64,
        rhomolar: f64,
    ) -> Result<()> {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        crate::derivs::check_alpha0(&self.alpha0_all(tau, delta), n_tau, n_delta, tau, delta)
    }

    /// Pressure [Pa] (upstream `calc_pressure`).
    pub fn pressure(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let dar_ddelta = self.alphar_all(tau, delta).d10;
        rhomolar * self.gas_constant * t * (1.0 + delta * dar_ddelta)
    }

    /// Molar enthalpy [J/mol] (upstream `calc_hmolar_nocache`).
    pub fn hmolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant * t * (1.0 + tau * (ideal.d01 + residual.d01) + delta * residual.d10)
    }

    /// Molar entropy [J/mol/K] (upstream `calc_smolar_nocache`).
    pub fn smolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant * (tau * (ideal.d01 + residual.d01) - ideal.d00 - residual.d00)
    }

    /// Molar internal energy [J/mol] (upstream `calc_umolar_nocache`).
    pub fn umolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant * t * tau * (ideal.d01 + residual.d01)
    }

    /// Molar Gibbs energy [J/mol] (upstream `calc_gibbsmolar_nocache`).
    pub fn gibbsmolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant * t * (1.0 + ideal.d00 + residual.d00 + delta * residual.d10)
    }

    /// Molar isochoric heat capacity [J/mol/K] (upstream `calc_cvmolar`).
    pub fn cvmolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        -self.gas_constant * tau.powi(2) * (residual.d02 + ideal.d02)
    }

    /// Molar isobaric heat capacity [J/mol/K] (upstream `calc_cpmolar`).
    pub fn cpmolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant
            * (-tau.powi(2) * (residual.d02 + ideal.d02)
                + (1.0 + delta * residual.d10 - delta * tau * residual.d11).powi(2)
                    / (1.0 + 2.0 * delta * residual.d10 + delta.powi(2) * residual.d20))
    }

    /// Speed of sound [m/s] for a homogeneous phase (upstream
    /// `calc_speed_sound`, single-phase branch).
    pub fn speed_sound(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        (self.gas_constant * t / self.molar_mass
            * (1.0 + 2.0 * delta * residual.d10 + delta.powi(2) * residual.d20
                - (1.0 + delta * residual.d10 - delta * tau * residual.d11).powi(2)
                    / (tau.powi(2) * (residual.d02 + ideal.d02))))
            .sqrt()
    }
}
