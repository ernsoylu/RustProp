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
use rustprop_core::Result;

impl HelmholtzEos {
    fn tau_delta(&self, t: f64, rhomolar: f64) -> (f64, f64) {
        (self.t_reducing / t, rhomolar / self.rhomolar_reducing)
    }

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

    /// Compressibility factor [-] (upstream `calc_compressibility_factor`,
    /// `HelmholtzEOSMixtureBackend.h`): `Z = 1 + delta*dalphar_dDelta` — NOT
    /// `p/(rho*R*T)`, which differs in the dome where upstream still serves
    /// the raw formula at the bulk density.
    pub fn compressibility_factor(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        1.0 + delta * self.alphar_all(tau, delta).d10
    }

    /// Ideal-gas molar isobaric heat capacity [J/mol/K] (upstream
    /// `calc_cpmolar_idealgas`).
    pub fn cp0molar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        self.gas_constant * (1.0 - tau.powi(2) * self.alpha0_all(tau, delta).d02)
    }

    /// Molar Helmholtz energy [J/mol] for a homogeneous phase (upstream
    /// `calc_helmholtzmolar`, homogeneous branch).
    pub fn helmholtzmolar(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant * t * (ideal.d00 + residual.d00)
    }

    /// Residual molar enthalpy [J/mol] (upstream inline
    /// `calc_hmolar_residual` — raw at the bulk density, dome included).
    pub fn hmolar_residual(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        self.gas_constant * t * (tau * residual.d01 + delta * residual.d10)
    }

    /// Residual molar entropy [J/mol/K] (upstream inline
    /// `calc_smolar_residual` — raw at the bulk density, dome included).
    pub fn smolar_residual(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        self.gas_constant * (tau * residual.d01 - residual.d00)
    }

    /// Ideal-gas molar enthalpy [J/mol] (upstream `calc_hmolar_idealgas`
    /// wrapper chain, `AbstractState.cpp`): `R*T*(1 + tau*dalpha0_dtau)`.
    pub fn hmolar_idealgas(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        self.gas_constant * t * (1.0 + tau * self.alpha0_all(tau, delta).d01)
    }

    /// Ideal-gas molar entropy [J/mol/K]: `R*(tau*dalpha0_dtau - alpha0)`.
    /// The alpha0 VALUE carries the log-delta term, so this is
    /// delta-dependent and evaluates at the bulk density like everything
    /// else.
    pub fn smolar_idealgas(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let ideal = self.alpha0_all(tau, delta);
        self.gas_constant * (tau * ideal.d01 - ideal.d00)
    }

    /// Ideal-gas molar internal energy [J/mol]: `R*T*tau*dalpha0_dtau`.
    pub fn umolar_idealgas(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        self.gas_constant * t * tau * self.alpha0_all(tau, delta).d01
    }

    /// Residual molar Gibbs energy [J/mol] (upstream inline
    /// `calc_gibbsmolar_residual`): `R*T*(alphar + delta*dalphar_ddelta)`,
    /// raw at the bulk density, dome included.
    pub fn gmolar_residual(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        let residual = self.alphar_all(tau, delta);
        self.gas_constant * t * (residual.d00 + delta * residual.d10)
    }

    /// Second virial coefficient [m^3/mol] (upstream `calc_Bvirial`,
    /// `HelmholtzEOSMixtureBackend.cpp`): `dalphar_ddelta` at the CURRENT
    /// tau and `delta = 1e-12`, over the reducing density.
    pub fn bvirial(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, _) = self.tau_delta(t, rhomolar);
        self.alphar_all(tau, 1e-12).d10 / self.rhomolar_reducing
    }

    /// Third virial coefficient [m^6/mol^2] (upstream `calc_Cvirial`).
    pub fn cvirial(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, _) = self.tau_delta(t, rhomolar);
        self.alphar_all(tau, 1e-12).d20 / (self.rhomolar_reducing * self.rhomolar_reducing)
    }

    /// `dBvirial/dT` (upstream `calc_dBvirial_dT`): the mixed derivative at
    /// `delta = 1e-12` chained through `dtau/dT = -Tr/T^2`.
    pub fn dbvirial_dt(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, _) = self.tau_delta(t, rhomolar);
        self.alphar_all(tau, 1e-12).d11 / self.rhomolar_reducing * (-self.t_reducing / (t * t))
    }

    /// `dCvirial/dT` (upstream `calc_dCvirial_dT`).
    pub fn dcvirial_dt(&self, t: f64, rhomolar: f64) -> f64 {
        let (tau, _) = self.tau_delta(t, rhomolar);
        self.alphar_all(tau, 1e-12).d21 / (self.rhomolar_reducing * self.rhomolar_reducing)
            * (-self.t_reducing / (t * t))
    }

    /// The keyed `Tau`/`Delta` pair — upstream serves them literally as
    /// `_reducing.T/_T` and `_rhomolar/_reducing.rhomolar` (bulk density in
    /// the dome).
    pub fn tau_delta_keyed(&self, t: f64, rhomolar: f64) -> (f64, f64) {
        self.tau_delta(t, rhomolar)
    }

    /// Both derivative matrices at the bulk state, for the keyed
    /// alphar/alpha0 output strings.
    pub fn alpha_all(
        &self,
        t: f64,
        rhomolar: f64,
    ) -> (crate::alpha::HelmholtzDerivs, crate::alpha::HelmholtzDerivs) {
        let (tau, delta) = self.tau_delta(t, rhomolar);
        (self.alphar_all(tau, delta), self.alpha0_all(tau, delta))
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
