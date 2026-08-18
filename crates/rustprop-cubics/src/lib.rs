//! Cubic engine — SRK and Peng-Robinson (PLAN.md Phase 7), a port of
//! CoolProp 8's `src/Backends/Cubics/{GeneralizedCubic,CubicBackend}.{h,cpp}`
//! restricted to pure fluids.
//!
//! Upstream conventions this port reproduces exactly:
//! - The cubic backend's reducing state is `T_r = 1 K`, `rho_r = 1 mol/m^3`
//!   (never overwritten for the SRK/PR backends), so `tau = 1/T` and
//!   `delta = rho` throughout the residual part.
//! - alpha0 is evaluated at the conventional reduced variables
//!   `tau* = Tc/T`, `delta* = rho/rhomolarc` (the JSON `rhomolarc`, not the
//!   Kazakov critical density) and rescaled; products like `tau·dalpha0_dtau`
//!   are invariant, which the property formulas below exploit.
//! - The SRK Omega_a/Omega_b literals carry upstream's two corrupted digits
//!   (relative error ~4e-10 vs the exact `1/(9(2^(1/3)-1))` and
//!   `(2^(1/3)-1)/3`); reproduced verbatim per the fidelity mandate.
//! - `rhomolar_critical` (the OUTPUT) is upstream's Kazakov curve fit, not
//!   the JSON `rhomolarc` and not the cubic's own critical density.
//! - QT/PQ saturation is the classic equal-Gibbs secant with the
//!   acentric-factor (Pitzer) seed; PT with three roots below `pc` runs an
//!   inner PQ flash for `Tsat(p)` to pick the branch, exactly upstream's
//!   `solver_rho_Tp(T, p, rho_guess)`.

// Upstream constants are carried verbatim, including digits beyond f64
// precision (the SRK Omega literals) — the data is bit-faithful by mandate.
#![allow(clippy::excessive_precision)]

mod superanc_tables;

use rustprop_core::fluid::{ChebyshevInterval, CubicFluid};
use rustprop_core::params::Phase;
use rustprop_core::{Error, Result};
use rustprop_heos::alpha::{Alpha0Eval, HelmholtzDerivs};
use rustprop_heos::derivs::DerivEos;
use rustprop_heos::solvers::{bounded_secant, secant, solve_cubic};

// Re-exported so the facade's cubics-only feature can drive the generic
// (T, rho) Jacobian machinery without depending on rustprop-heos directly.
pub use rustprop_heos::derivs::StateDerivs;

/// Upstream `R_U_CODATA` (Configuration.h) — the gas constant handed to
/// `SRKBackend`/`PengRobinsonBackend` by default.
pub const R_U_CODATA: f64 = 8.314_462_618_153_24;

/// Which cubic equation of state.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CubicKind {
    Srk,
    PengRobinson,
}

/// The derivative bundle the property formulas consume, in the cubic's own
/// (tau = 1/T, delta = rho) variables. `dxy` = x delta-derivatives, y
/// tau-derivatives. Third order is carried for the generic (T, rho)
/// partial-derivative machinery (PIP, the compressibilities, virial
/// T-derivatives).
#[derive(Default, Clone, Copy)]
pub struct CubicDerivs {
    pub d00: f64,
    pub d10: f64,
    pub d01: f64,
    pub d20: f64,
    pub d11: f64,
    pub d02: f64,
    pub d30: f64,
    pub d21: f64,
    pub d12: f64,
    pub d03: f64,
}

/// One pure fluid under one cubic EOS.
pub struct CubicEos {
    pub kind: CubicKind,
    pub fluid: &'static CubicFluid,
    /// `a0 = Omega_a * R^2 * Tc^2 / pc`
    a0: f64,
    /// `b = Omega_b * R * Tc / pc`
    b: f64,
    /// `m(acentric)` of the Soave alpha function.
    m: f64,
    delta_1: f64,
    delta_2: f64,
    alpha0: Alpha0Eval,
}

impl CubicEos {
    pub fn new(kind: CubicKind, fluid: &'static CubicFluid) -> Self {
        let r = R_U_CODATA;
        let (tc, pc, omega) = (fluid.tc, fluid.pc, fluid.acentric);
        // Upstream literals verbatim (GeneralizedCubic.cpp:769-800); the SRK
        // pair differs from the exact constants in the 10th significant
        // digit — see the module docs.
        let (a0, b, m, delta_1, delta_2) = match kind {
            CubicKind::Srk => (
                0.42748023335403414043900347952220 * r * r * tc * tc / pc,
                0.08664034999649577215890158147700 * r * tc / pc,
                0.480 + 1.574 * omega - 0.176 * omega * omega,
                1.0,
                0.0,
            ),
            CubicKind::PengRobinson => (
                0.45723552892138218938000849856422 * r * r * tc * tc / pc,
                0.07779607390388455972148597969400 * r * tc / pc,
                0.37464 + 1.54226 * omega - 0.26992 * omega * omega,
                1.0 + std::f64::consts::SQRT_2,
                1.0 - std::f64::consts::SQRT_2,
            ),
        };
        CubicEos {
            kind,
            fluid,
            a0,
            b,
            m,
            delta_1,
            delta_2,
            alpha0: Alpha0Eval::new(fluid.alpha0),
        }
    }

    pub fn gas_constant(&self) -> f64 {
        R_U_CODATA
    }

    pub fn molar_mass(&self) -> f64 {
        self.fluid.molemass
    }

    /// The attractive parameter `b` [m^3/mol].
    pub fn b(&self) -> f64 {
        self.b
    }

    /// `a(T)` and its first four tau-derivatives at `tau = 1/T` — upstream
    /// `BasicMathiasCopemanAlphaFunction::calc_all_terms` (the batched hot
    /// path), the default alpha function for every bundled fluid.
    pub fn a_terms(&self, tau: f64) -> [f64; 5] {
        let tr_over_tci = 1.0 / self.fluid.tc; // T_r / Tc with T_r = 1
        let sqrt_tr_tci = tr_over_tci.sqrt();
        let (a0, m) = (self.a0, self.m);
        let sq = sqrt_tr_tci / tau.sqrt();
        let b_ = 1.0 + m * (1.0 - sq);
        let t15 = sqrt_tr_tci / (tau * tau.sqrt());
        let t25 = t15 / tau;
        let t35 = t25 / tau;
        let t45 = t35 / tau;
        let t3 = tr_over_tci / (tau * tau * tau);
        let t4 = t3 / tau;
        let t5 = t4 / tau;
        [
            a0 * b_ * b_,
            a0 * m * b_ * t15,
            a0 * m * 0.5 * (m * t3 - 3.0 * b_ * t25),
            (3.0 / 4.0) * a0 * m * (-3.0 * m * t4 + 5.0 * b_ * t35),
            (3.0 / 8.0) * a0 * m * (29.0 * m * t5 - 35.0 * b_ * t45),
        ]
    }

    /// `d^n/dtau^n (tau * a)` for n = 0..=3 (upstream `tau_times_a`,
    /// Leibniz: `tau*a^(n) + n*a^(n-1)`).
    fn tau_times_a(&self, tau: f64, a: &[f64; 5], n: usize) -> f64 {
        if n == 0 {
            tau * a[0]
        } else {
            tau * a[n] + (n as f64) * a[n - 1]
        }
    }

    /// `psi^- = -ln(1 - b*rho)` and its delta-derivatives (upstream
    /// `psi_minus`, pure fluid, rho_r = 1, c = 0).
    fn psi_minus(&self, delta: f64, idelta: usize) -> f64 {
        let bracket = 1.0 - self.b * delta;
        let r = self.b / bracket;
        match idelta {
            0 => -bracket.ln(),
            1 => r,
            2 => r * r,
            3 => 2.0 * r * r * r,
            _ => unreachable!(),
        }
    }

    /// `PI_12 = (1 + Delta_1*b*rho)(1 + Delta_2*b*rho)` and derivatives.
    fn pi_12(&self, delta: f64, idelta: usize) -> f64 {
        let b = self.b;
        match idelta {
            0 => (1.0 + self.delta_1 * b * delta) * (1.0 + self.delta_2 * b * delta),
            1 => {
                2.0 * (b * self.delta_1) * (b * self.delta_2) * delta
                    + (self.delta_1 + self.delta_2) * b
            }
            2 => 2.0 * (self.delta_1 * b) * (self.delta_2 * b),
            _ => unreachable!(),
        }
    }

    /// `psi^+` and its delta-derivatives (upstream `psi_plus`; `c_term = 1/b`).
    fn psi_plus(&self, delta: f64, idelta: usize) -> f64 {
        match idelta {
            0 => {
                // A_term groups delta * (Delta_i*bm) + 1 (GeneralizedCubic.h:619-623
                // with cm = 0, rho_r = 1).
                let a_term = ((delta * (self.delta_1 * self.b) + 1.0)
                    / (delta * (self.delta_2 * self.b) + 1.0))
                    .ln();
                // Upstream multiplies by c_term = 1/bm (GeneralizedCubic.cpp:553-556).
                // The ulp between these associations and direct division decides the
                // equal-Gibbs secant at sub-pascal PQ, where vapour-root cancellation
                // sits at the residual's tolerance floor.
                a_term * (1.0 / self.b) / (self.delta_1 - self.delta_2)
            }
            1 => 1.0 / self.pi_12(delta, 0),
            2 => {
                let p0 = self.pi_12(delta, 0);
                let p1 = self.pi_12(delta, 1);
                -p1 / (p0 * p0)
            }
            3 => {
                // Upstream: rho_r*(-p0*p2 + 2*p1*p1)/p0^3 with rho_r = 1
                // (PI_12(delta, 3) is identically zero, so no p3 term).
                let p0 = self.pi_12(delta, 0);
                let p1 = self.pi_12(delta, 1);
                let p2 = self.pi_12(delta, 2);
                (-p0 * p2 + 2.0 * p1 * p1) / (p0 * p0 * p0)
            }
            _ => unreachable!(),
        }
    }

    /// All residual derivatives through third order in the cubic's own
    /// (tau = 1/T, delta = rho) variables — upstream
    /// `CubicResidualHelmholtz::all`'s product structure: every mixed
    /// derivative is `psi_minus[n_d] - (tau*a)^(n_t) / R * psi_plus[n_d]`
    /// (T_r = 1 makes the denominator plain R), and `psi_minus` vanishes for
    /// any tau-derivative.
    pub fn alphar_all(&self, tau: f64, delta: f64) -> CubicDerivs {
        let a = self.a_terms(tau);
        let r = R_U_CODATA;
        let ta0 = self.tau_times_a(tau, &a, 0) / r;
        let ta1 = self.tau_times_a(tau, &a, 1) / r;
        let ta2 = self.tau_times_a(tau, &a, 2) / r;
        let ta3 = self.tau_times_a(tau, &a, 3) / r;
        let pm = [
            self.psi_minus(delta, 0),
            self.psi_minus(delta, 1),
            self.psi_minus(delta, 2),
            self.psi_minus(delta, 3),
        ];
        let pp = [
            self.psi_plus(delta, 0),
            self.psi_plus(delta, 1),
            self.psi_plus(delta, 2),
            self.psi_plus(delta, 3),
        ];
        CubicDerivs {
            d00: pm[0] - ta0 * pp[0],
            d10: pm[1] - ta0 * pp[1],
            d01: -ta1 * pp[0],
            d20: pm[2] - ta0 * pp[2],
            d11: -ta1 * pp[1],
            d02: -ta2 * pp[0],
            d30: pm[3] - ta0 * pp[3],
            d21: -ta1 * pp[2],
            d12: -ta2 * pp[1],
            d03: -ta3 * pp[0],
        }
    }

    /// alpha0 derivatives at the state, expressed as the INVARIANT products
    /// the property formulas need: `(a0, tau*da0_dtau, tau^2*d2a0_dtau2)` —
    /// evaluated at `tau* = Tc/T`, `delta* = rho/rhomolarc` where those
    /// products equal the cubic-variable ones.
    fn alpha0_products(&self, t: f64, rho: f64) -> (f64, f64, f64) {
        let tau_star = self.fluid.tc / t;
        let delta_star = rho / self.fluid.rhomolarc;
        let d = self.alpha0.all(tau_star, delta_star);
        (d.d00, tau_star * d.d01, tau_star * tau_star * d.d02)
    }

    /// Upstream `calc_alpha0_deriv_nocache`'s closing `ValidNumber` gate at
    /// the state `(t, rhomolar)`. `AbstractCubicBackend` overrides only the
    /// RESIDUAL `calc_alphar_deriv_nocache`, so it inherits the ideal-gas
    /// function and its throw whole; the message quotes the HEOS-convention
    /// `tau = _reducing.T/_T` and `delta = _rhomolar/_reducing.rhomolar`,
    /// which for a cubic are `1/T` and `rhomolar` (Tr = rhor = 1).
    pub fn check_alpha0_deriv(
        &self,
        n_tau: u32,
        n_delta: u32,
        t: f64,
        rhomolar: f64,
    ) -> Result<()> {
        let (tau, delta) = (1.0 / t, rhomolar);
        rustprop_heos::derivs::check_alpha0(
            &self.alpha0_all(tau, delta),
            n_tau,
            n_delta,
            tau,
            delta,
        )
    }

    /// The full alpha0 derivative matrix in the CUBIC (tau = 1/T,
    /// delta = rho) convention — upstream `calc_alpha0_deriv_nocache`
    /// (HelmholtzEOSMixtureBackend.cpp:3580): the fabricated fluid's alpha0
    /// container evaluated at `taustar = (Tc/Tr)*tau`, `deltastar =
    /// (rhor/rhomolarc)*delta`, then each derivative rescaled
    /// `* (rhor/rhomolarc)^nDelta / (Tr/Tc)^nTau` (Tr = rhor = 1 here, so
    /// `* (1/rhomolarc)^nDelta / (1/Tc)^nTau`), with `powf` mirroring
    /// upstream's `pow` calls. This is the rescale behind every keyed
    /// `dalpha0_*` output and the generic derivative machinery; only
    /// `calc_smolar`'s batched accessor skips it (see `smolar`).
    /// Fourth-order fields get the same chain-rule scale for consistency;
    /// upstream has no cubic call surface for them.
    pub fn alpha0_all(&self, tau: f64, delta: f64) -> HelmholtzDerivs {
        let tc = self.fluid.tc;
        let rc = self.fluid.rhomolarc;
        let tau_star = tc * tau;
        let delta_star = delta / rc;
        let d = self.alpha0.all(tau_star, delta_star);
        let sd = 1.0 / rc; // rhor / rhomolarc
        let st = 1.0 / tc; // Tr / Tc
        let s = |v: f64, nd: i32, nt: i32| v * sd.powf(nd as f64) / st.powf(nt as f64);
        HelmholtzDerivs {
            d00: s(d.d00, 0, 0),
            d10: s(d.d10, 1, 0),
            d01: s(d.d01, 0, 1),
            d20: s(d.d20, 2, 0),
            d11: s(d.d11, 1, 1),
            d02: s(d.d02, 0, 2),
            d30: s(d.d30, 3, 0),
            d21: s(d.d21, 2, 1),
            d12: s(d.d12, 1, 2),
            d03: s(d.d03, 0, 3),
            d40: s(d.d40, 4, 0),
            d31: s(d.d31, 3, 1),
            d22: s(d.d22, 2, 2),
            d13: s(d.d13, 1, 3),
            d04: s(d.d04, 0, 4),
        }
    }

    /// `p = rho R T (1 + delta * dalphar_ddelta)` — the corrected form of
    /// upstream `calc_pressure_nocache` (whose use of cached members is a
    /// latent bug that never fires because callers pass the live state).
    pub fn pressure(&self, t: f64, rho: f64) -> f64 {
        let d = self.alphar_all(1.0 / t, rho);
        rho * R_U_CODATA * t * (1.0 + rho * d.d10)
    }

    /// Molar enthalpy [J/mol].
    pub fn hmolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let (_, tau_da0, _) = self.alpha0_products(t, rho);
        let tau = 1.0 / t;
        R_U_CODATA * t * (1.0 + tau_da0 + tau * dr.d01 + rho * dr.d10)
    }

    /// Molar entropy [J/mol/K] — upstream DEFECT reproduced deliberately:
    /// `calc_smolar` alone reads the batched
    /// `calc_all_alpha0_derivs_nocache`, whose cubic-backend output is the
    /// tau*-derivative WITHOUT the rescale back to the cubic tau, and then
    /// multiplies it by the cubic `tau = 1/T` — i.e. `tau*da0/dtau*` where
    /// every sibling property forms the invariant `tau**da0/dtau*`. The
    /// entropy datum is therefore shifted by `R*(tau* - tau)*da0'` relative
    /// to the dimensionally consistent formula (h/u/cv/cp/w all use the
    /// rescaled accessors and are unaffected). Oracle-verified.
    pub fn smolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let tau_star = self.fluid.tc / t;
        let delta_star = rho / self.fluid.rhomolarc;
        let d0 = self.alpha0.all(tau_star, delta_star);
        let tau = 1.0 / t;
        R_U_CODATA * (tau * (d0.d01 + dr.d01) - d0.d00 - dr.d00)
    }

    /// Molar internal energy [J/mol].
    pub fn umolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let (_, tau_da0, _) = self.alpha0_products(t, rho);
        let tau = 1.0 / t;
        R_U_CODATA * t * (tau_da0 + tau * dr.d01)
    }

    /// Molar isochoric heat capacity [J/mol/K].
    pub fn cvmolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let (_, _, tau2_d2a0) = self.alpha0_products(t, rho);
        let tau = 1.0 / t;
        -R_U_CODATA * (tau2_d2a0 + tau * tau * dr.d02)
    }

    /// Molar isobaric heat capacity [J/mol/K].
    pub fn cpmolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let tau = 1.0 / t;
        let num = 1.0 + rho * dr.d10 - rho * tau * dr.d11;
        let den = 1.0 + 2.0 * rho * dr.d10 + rho * rho * dr.d20;
        self.cvmolar(t, rho) + R_U_CODATA * num * num / den
    }

    /// Speed of sound [m/s].
    pub fn speed_sound(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let (_, _, tau2_d2a0) = self.alpha0_products(t, rho);
        let tau = 1.0 / t;
        let num = 1.0 + rho * dr.d10 - rho * tau * dr.d11;
        let w2 = R_U_CODATA * t / self.fluid.molemass
            * (1.0 + 2.0 * rho * dr.d10 + rho * rho * dr.d20
                - num * num / (tau2_d2a0 + tau * tau * dr.d02));
        w2.sqrt()
    }

    /// Compressibility factor [-] — inherited `calc_compressibility_factor`:
    /// `Z = 1 + delta*dalphar_ddelta` at the bulk state (raw in the dome,
    /// like cp/cv — NOT `p/(rho*R*T)`, which differs two-phase).
    pub fn compressibility_factor(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        1.0 + rho * dr.d10
    }

    /// Ideal-gas molar isobaric heat capacity [J/mol/K] — inherited
    /// `calc_cpmolar_idealgas`: `R*(1 - tau^2*d2alpha0_dtau2)` (the invariant
    /// tau-squared product, convention-free).
    pub fn cp0molar(&self, t: f64, rho: f64) -> f64 {
        let (_, _, tau2_d2a0) = self.alpha0_products(t, rho);
        R_U_CODATA * (1.0 - tau2_d2a0)
    }

    /// Residual molar enthalpy [J/mol] — inherited `calc_hmolar_residual`:
    /// raw at the bulk state, dome included.
    pub fn hmolar_residual(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let tau = 1.0 / t;
        R_U_CODATA * t * (tau * dr.d01 + rho * dr.d10)
    }

    /// Residual molar entropy [J/mol/K] — inherited `calc_smolar_residual`:
    /// raw at the bulk state, dome included.
    pub fn smolar_residual(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let tau = 1.0 / t;
        R_U_CODATA * (tau * dr.d01 - dr.d00)
    }

    /// Molar Gibbs energy [J/mol] — inherited `calc_gibbsmolar_nocache`
    /// (`R*T*(1 + alpha0 + alphar + delta*dalphar_ddelta)`; the alpha0 VALUE
    /// needs no rescale, so the batched-accessor entropy defect does not
    /// touch it).
    pub fn gibbsmolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let (a0, _, _) = self.alpha0_products(t, rho);
        R_U_CODATA * t * (1.0 + a0 + dr.d00 + rho * dr.d10)
    }

    /// Molar Helmholtz energy [J/mol] — inherited `calc_helmholtzmolar`
    /// homogeneous branch: `A = R*T*(alpha0 + alphar)`.
    pub fn helmholtzmolar(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        let (a0, _, _) = self.alpha0_products(t, rho);
        R_U_CODATA * t * (a0 + dr.d00)
    }

    /// Residual molar Gibbs energy [J/mol] — inherited
    /// `calc_gibbsmolar_residual`: `R*T*(alphar + delta*dalphar_ddelta)`,
    /// raw at the bulk state.
    pub fn gibbsmolar_residual(&self, t: f64, rho: f64) -> f64 {
        let dr = self.alphar_all(1.0 / t, rho);
        R_U_CODATA * t * (dr.d00 + rho * dr.d10)
    }

    /// Ideal-gas molar enthalpy [J/mol] — `AbstractState::hmolar_idealgas`:
    /// `R*T*(1 + tau*dalpha0_dtau)` on the RESCALED individual accessor,
    /// whose `tau*da0_dtau` is the convention-invariant product (unlike
    /// `calc_smolar`'s batched-accessor defect).
    pub fn hmolar_idealgas(&self, t: f64, rho: f64) -> f64 {
        let (_, tau_da0, _) = self.alpha0_products(t, rho);
        R_U_CODATA * t * (1.0 + tau_da0)
    }

    /// Ideal-gas molar entropy [J/mol/K] — `AbstractState::smolar_idealgas`:
    /// `R*(tau*dalpha0_dtau - alpha0)`. Both operands come from the RESCALED
    /// individual accessors (`dalpha0_dTau()` / `alpha0()` →
    /// `calc_alpha0_deriv_nocache`), NOT the batched
    /// `calc_all_alpha0_derivs_nocache` that `calc_smolar` uses — so the
    /// tau-product is invariant and this shares no part of the smolar
    /// defect. Oracle: 67.93894818182274 at SRK::Water T=300, P=101325.
    pub fn smolar_idealgas(&self, t: f64, rho: f64) -> f64 {
        let (a0, tau_da0, _) = self.alpha0_products(t, rho);
        R_U_CODATA * (tau_da0 - a0)
    }

    /// Ideal-gas molar internal energy [J/mol] —
    /// `AbstractState::umolar_idealgas`: `R*T*tau*dalpha0_dtau`.
    pub fn umolar_idealgas(&self, t: f64, rho: f64) -> f64 {
        let (_, tau_da0, _) = self.alpha0_products(t, rho);
        R_U_CODATA * t * tau_da0
    }

    /// Second virial coefficient [m^3/mol] — inherited `calc_Bvirial`:
    /// `(1/rhor)*dalphar_ddelta(tau, 1e-12)` with rhor = 1.
    pub fn bvirial(&self, t: f64) -> f64 {
        self.alphar_all(1.0 / t, 1e-12).d10
    }

    /// Third virial coefficient [m^6/mol^2] — inherited `calc_Cvirial`:
    /// `(1/rhor^2)*d2alphar_ddelta2(tau, 1e-12)`.
    pub fn cvirial(&self, t: f64) -> f64 {
        self.alphar_all(1.0 / t, 1e-12).d20
    }

    /// `dB/dT` — inherited `calc_dBvirial_dT`:
    /// `(1/rhor)*d2alphar_ddelta_dtau(tau, 1e-12) * dtau_dT` with
    /// `dtau_dT = -Tr/T^2 = -1/T^2`.
    pub fn dbvirial_dt(&self, t: f64) -> f64 {
        self.alphar_all(1.0 / t, 1e-12).d11 * (-1.0 / (t * t))
    }

    /// `dC/dT` — inherited `calc_dCvirial_dT`:
    /// `(1/rhor^2)*d3alphar_ddelta2_dtau(tau, 1e-12) * dtau_dT`.
    pub fn dcvirial_dt(&self, t: f64) -> f64 {
        self.alphar_all(1.0 / t, 1e-12).d21 * (-1.0 / (t * t))
    }

    /// Upstream `calc_rhomolar_critical` — the Kazakov curve fit (NOT the
    /// JSON `rhomolarc`, and not the cubic's own critical density).
    pub fn rhomolar_critical(&self) -> f64 {
        let v_c_lmol = 2.14107171795 * (self.fluid.tc / self.fluid.pc * 1000.0) + 0.00773144012514;
        1.0 / (v_c_lmol / 1000.0)
    }

    /// The three-root cubic solve in rho (upstream `rho_Tp_cubic`):
    /// coefficients from `d1 = c - b`, `d2 = c + Delta_1 b`, `d3 = c +
    /// Delta_2 b` with `c = 0`; roots returned sorted ascending.
    pub fn rho_tp_cubic(&self, t: f64, p: f64) -> (i32, f64, f64, f64) {
        let r = R_U_CODATA;
        let am = self.a_terms(1.0 / t)[0];
        let bm = self.b;
        let d1 = -bm;
        let d2 = self.delta_1 * bm;
        let d3 = self.delta_2 * bm;
        let crho0 = -p;
        let crho1 = r * t - p * (d1 + d2 + d3);
        let crho2 = r * t * (d2 + d3) - p * (d1 * (d2 + d3) + d2 * d3) - am;
        let crho3 = r * t * d2 * d3 - p * d1 * d2 * d3 - d1 * am;
        let (n, r0, r1, r2) = solve_cubic(crho3, crho2, crho1, crho0);
        // Upstream `sort3` verbatim: a fixed three-comparison network on `>`.
        // The shape matters, not just the outcome — `solve_cubic` leaves NaN
        // in the slots it did not fill (the linear and quadratic branches
        // taken when the leading coefficient underflows, which happens for
        // real states: the Soave alpha function passes through zero at high
        // reduced temperature, and `crho3` carries a factor of `am`). Every
        // comparison against NaN is false, so upstream performs no swaps and
        // passes the NaN through. A `partial_cmp().unwrap()` sort panics
        // there instead — found by the seeded acceptance sweep.
        let (mut a, mut b, mut c) = (r0, r1, r2);
        if a > b {
            core::mem::swap(&mut a, &mut b);
        }
        if a > c {
            core::mem::swap(&mut a, &mut c);
        }
        if b > c {
            core::mem::swap(&mut b, &mut c);
        }
        (n, a, b, c)
    }

    /// Upstream `recalculate_singlephase_phase` with the cubic's critical
    /// parameters (JSON Tc/pc, Kazakov rho_c).
    fn singlephase_phase(&self, t: f64, p: f64, rho: f64) -> Phase {
        if p > self.fluid.pc {
            if t > self.fluid.tc {
                Phase::Supercritical
            } else {
                Phase::SupercriticalLiquid
            }
        } else if t > self.fluid.tc {
            Phase::SupercriticalGas
        } else if rho > self.rhomolar_critical() {
            Phase::Liquid
        } else {
            Phase::Gas
        }
    }

    /// Upstream `solver_rho_Tp(T, p, rho_guess)` — the PT path with no
    /// imposed phase: one root takes it; three roots below `pc` run a PQ
    /// saturation flash for `Tsat(p)` and pick gas above it (smallest
    /// positive root below the saturated-vapor density) or liquid (largest
    /// root); three roots at or above `pc` throw.
    pub fn pt_flash(&self, t: f64, p: f64) -> Result<(f64, Phase)> {
        let (n, r0, r1, r2) = self.rho_tp_cubic(t, p);
        let rho = if n == 1 {
            r0
        } else if n == 3 {
            if p < self.fluid.pc {
                let sat = self.pq_flash(p, 0.0)?;
                if t > sat.t {
                    let rho_v = sat.rho_v;
                    if r0 > 0.0 && r0 < rho_v {
                        r0
                    } else if r1 > 0.0 && r1 < rho_v {
                        r1
                    } else {
                        return Err(Error::Value(format!(
                            "Unable to find gaseous density for T: {t} K, p: {p} Pa"
                        )));
                    }
                } else {
                    r2
                }
            } else {
                return Err(Error::Value(
                    "Cubic has three roots, but phase not imposed and guess density not provided"
                        .into(),
                ));
            }
        } else {
            return Err(Error::Value("Obtained neither 1 nor three roots".into()));
        };
        // Upstream QUIRK, reproduced deliberately: on the PT path
        // `recalculate_singlephase_phase` runs INSIDE `solver_rho_Tp`, before
        // `update` assigns `_rhomolar` from the return value — so it reads
        // the `clear()`ed sentinel `-_HUGE`, and the subcritical
        // `rhomolar() > rhomolar_critical()` test is always false: every
        // subcritical PT state reports `iphase_gas`, liquid-like roots
        // included (wheel: SRK::Water at 300 K / 101325 Pa -> Phase 5.0 with
        // a 41892 mol/m^3 root). The stale density is modeled explicitly;
        // the DmolarT path keeps the real-density classification.
        Ok((rho, self.singlephase_phase(t, p, f64::NEG_INFINITY)))
    }

    /// The equal-Gibbs residual of upstream `SaturationResidual::call`:
    /// `ln dV - ln dL + ar(dV) - ar(dL) + dV*ar_d(dV) - dL*ar_d(dL)` with
    /// `dL` the largest and `dV` the smallest cubic root. When the cubic has
    /// one real root the three roots coincide and the residual is exactly
    /// zero — upstream's behavior, reproduced verbatim (the secant seeds
    /// start inside the dome for the bundled fluids).
    fn saturation_residual(&self, t: f64, p: f64) -> (f64, f64, f64) {
        let (_n, r0, _r1, r2) = self.rho_tp_cubic(t, p);
        let tau = 1.0 / t;
        let delta_l = r2;
        let delta_v = r0;
        let dl = self.alphar_all(tau, delta_l);
        let dv = self.alphar_all(tau, delta_v);
        let delta_gibbs =
            delta_v.ln() - delta_l.ln() + (dv.d00 - dl.d00) + (delta_v * dv.d10 - delta_l * dl.d10);
        (delta_gibbs, delta_l, delta_v)
    }

    /// Upstream `saturation(QT_INPUTS)`: Pitzer-seeded bounded secant on p.
    pub fn qt_flash(&self, t: f64, q: f64) -> Result<CubicSatState> {
        let (tc, pc, omega) = (self.fluid.tc, self.fluid.pc, self.fluid.acentric);
        let neg_log10_pr = (omega + 1.0) / (1.0 / 0.7 - 1.0) * (tc / t - 1.0);
        let ps_est = pc * 10.0f64.powf(-neg_log10_pr);
        let mut rho_l = f64::NAN;
        let mut rho_v = f64::NAN;
        let ps = bounded_secant(
            |p: f64| {
                let (resid, dl, dv) = self.saturation_residual(t, p);
                rho_l = dl;
                rho_v = dv;
                resid
            },
            ps_est,
            1e-10,
            pc,
            -0.01 * ps_est,
            1e-5,
            100,
        )?;
        Ok(CubicSatState {
            t,
            p: ps,
            rho_l,
            rho_v,
            rhomolar: 1.0 / (q / rho_v + (1.0 - q) / rho_l),
            q,
        })
    }

    /// Upstream `saturation(PQ_INPUTS)`: acentric-relation-seeded secant on T.
    pub fn pq_flash(&self, p: f64, q: f64) -> Result<CubicSatState> {
        let (tc, pc, omega) = (self.fluid.tc, self.fluid.pc, self.fluid.acentric);
        let theta = -(p / pc).log10() * (1.0 / 0.7 - 1.0) / (omega + 1.0);
        let ts_est = tc / (theta + 1.0);
        let mut rho_l = f64::NAN;
        let mut rho_v = f64::NAN;
        let ts = secant(
            |t: f64| {
                let (resid, dl, dv) = self.saturation_residual(t, p);
                rho_l = dl;
                rho_v = dv;
                resid
            },
            ts_est,
            -0.1,
            1e-10,
            100,
        )?;
        Ok(CubicSatState {
            t: ts,
            p,
            rho_l,
            rho_v,
            rhomolar: 1.0 / (q / rho_v + (1.0 - q) / rho_l),
            q,
        })
    }
}

impl CubicEos {
    fn superanc_tables(
        &self,
    ) -> (
        &'static [ChebyshevInterval],
        &'static [ChebyshevInterval],
        &'static [ChebyshevInterval],
    ) {
        match self.kind {
            CubicKind::Srk => (
                superanc_tables::SRK_P,
                superanc_tables::SRK_RHO_L,
                superanc_tables::SRK_RHO_V,
            ),
            CubicKind::PengRobinson => (
                superanc_tables::PR_P,
                superanc_tables::PR_RHO_L,
                superanc_tables::PR_RHO_V,
            ),
        }
    }

    /// `Ttilde = R * T * b / a(T)` — the generic-cubic reduction variable
    /// (the fluid enters only through `a(T)` and `b`).
    fn ttilde(&self, t: f64) -> f64 {
        let am = self.a_terms(1.0 / t)[0];
        R_U_CODATA * t * self.b / am
    }

    /// Upstream `calc_superanc_Tmax`: invert `Ttilde_max` analytically —
    /// `alpha(Tc) = 1`, so `a(Tc) = a0` exactly.
    pub fn superanc_tmax(&self) -> f64 {
        let (_, rho_l, _) = self.superanc_tables();
        let ttilde_max = rho_l[rho_l.len() - 1].xmax;
        ttilde_max * self.a0 / (R_U_CODATA * self.b)
    }

    /// Saturation properties off the cubic superancillary (upstream
    /// `calc_saturation_ancillary`): `p = ptilde * a / b^2`,
    /// `rho = rhotilde / b`.
    pub fn superanc_sat(&self, t: f64) -> Result<(f64, f64, f64)> {
        let (p_tab, rho_l_tab, rho_v_tab) = self.superanc_tables();
        let ttilde = self.ttilde(t);
        let am = self.a_terms(1.0 / t)[0];
        let p = superanc_eval(p_tab, ttilde)? * am / (self.b * self.b);
        let rho_l = superanc_eval(rho_l_tab, ttilde)? / self.b;
        let rho_v = superanc_eval(rho_v_tab, ttilde)? / self.b;
        Ok((p, rho_l, rho_v))
    }

    /// Upstream `update_DmolarT` — pure SRK/PR path: the superancillary
    /// brackets the dome at T; `rho >= rhoL` is liquid, `rho <= rhoV` gas,
    /// interior states take the saturation pressure and the lever-rule
    /// quality. At or above the superancillary Tmax the state is single
    /// phase directly.
    pub fn dmolar_t_flash(&self, rho: f64, t: f64) -> Result<CubicDtState> {
        if t >= self.superanc_tmax() {
            let p = self.pressure(t, rho);
            return Ok(CubicDtState::Single {
                t,
                rho,
                p,
                phase: self.singlephase_phase(t, p, rho),
            });
        }
        let (p_sat, rho_l_sat, rho_v_sat) = self.superanc_sat(t)?;
        if rho >= rho_l_sat {
            let p = self.pressure(t, rho);
            Ok(CubicDtState::Single {
                t,
                rho,
                p,
                phase: Phase::Liquid,
            })
        } else if rho <= rho_v_sat {
            let p = self.pressure(t, rho);
            Ok(CubicDtState::Single {
                t,
                rho,
                p,
                phase: Phase::Gas,
            })
        } else {
            // Inside the dome: lever rule 1/rho = (1-Q)/rhoL + Q/rhoV.
            let q = (rho_v_sat * (rho_l_sat - rho)) / (rho * (rho_l_sat - rho_v_sat));
            Ok(CubicDtState::TwoPhase(CubicSatState {
                t,
                p: p_sat,
                rho_l: rho_l_sat,
                rho_v: rho_v_sat,
                rhomolar: rho,
                q,
            }))
        }
    }
}

/// The generic (T, rho) partial-derivative machinery
/// (`rustprop_heos::derivs::StateDerivs`) over the cubic backend — upstream
/// runs `AbstractState::get_dT_drho` over the cubic's cached derivative
/// matrix exactly as over HEOS's. The cubic's own convention is Tr = 1,
/// rho_r = 1, R = R_U_CODATA; alphar comes from
/// `CubicResidualHelmholtz::all`, alpha0 from `calc_alpha0_deriv_nocache`'s
/// rescale (`alpha0_all` above). Fourth-order alphar fields are left zero:
/// `StateDerivs` reads at most third order, and upstream's cubic `all` fills
/// them only for its own 4th-order call surface, which is not ported.
impl DerivEos for CubicEos {
    fn alphar_all(&self, tau: f64, delta: f64) -> HelmholtzDerivs {
        let d = CubicEos::alphar_all(self, tau, delta);
        HelmholtzDerivs {
            d00: d.d00,
            d10: d.d10,
            d01: d.d01,
            d20: d.d20,
            d11: d.d11,
            d02: d.d02,
            d30: d.d30,
            d21: d.d21,
            d12: d.d12,
            d03: d.d03,
            ..Default::default()
        }
    }
    fn alpha0_all(&self, tau: f64, delta: f64) -> HelmholtzDerivs {
        CubicEos::alpha0_all(self, tau, delta)
    }
    fn t_reducing(&self) -> f64 {
        1.0
    }
    fn rhomolar_reducing(&self) -> f64 {
        1.0
    }
    fn gas_constant(&self) -> f64 {
        R_U_CODATA
    }
    fn molar_mass(&self) -> f64 {
        self.fluid.molemass
    }
}

/// A (Dmolar, T) cubic state.
pub enum CubicDtState {
    Single {
        t: f64,
        rho: f64,
        p: f64,
        phase: Phase,
    },
    TwoPhase(CubicSatState),
}

/// Clenshaw + Knuth-midpoint interval lookup over one superancillary table
/// (upstream `cubicsuperancillary.h`'s header-local implementations), with
/// upstream's out-of-domain errors.
fn superanc_eval(table: &[ChebyshevInterval], ttilde: f64) -> Result<f64> {
    let lo = table[0].xmin;
    let hi = table[table.len() - 1].xmax;
    if ttilde < lo {
        return Err(Error::Value(format!(
            "Ttilde ({ttilde}) is below the minimum of {lo}"
        )));
    }
    if ttilde > hi {
        return Err(Error::Value(format!(
            "Ttilde ({ttilde}) is above the maximum of {hi}"
        )));
    }
    // Knuth-midpoint bisection on xmin.
    let midpoint_knuth = |a: i32, b: i32| (a & b) + ((a ^ b) >> 1);
    let mut il = 0i32;
    let mut ir = table.len() as i32 - 1;
    while ir - il > 1 {
        let im = midpoint_knuth(il, ir);
        if ttilde >= table[im as usize].xmin {
            il = im;
        } else {
            ir = im;
        }
    }
    let e = if ttilde < table[il as usize].xmax {
        &table[il as usize]
    } else {
        &table[ir as usize]
    };
    // Clenshaw.
    let xscaled = (2.0 * ttilde - (e.xmax + e.xmin)) / (e.xmax - e.xmin);
    let norder = e.coef.len() - 1;
    let mut u_kp1 = e.coef[norder];
    let mut u_kp2 = 0.0;
    for k in (1..norder).rev() {
        let u_k = 2.0 * xscaled * u_kp1 - u_kp2 + e.coef[k];
        u_kp2 = u_kp1;
        u_kp1 = u_k;
    }
    Ok(e.coef[0] + xscaled * u_kp1 - u_kp2)
}

/// A converged cubic saturation state (both branch densities at one T, p).
pub struct CubicSatState {
    pub t: f64,
    pub p: f64,
    pub rho_l: f64,
    pub rho_v: f64,
    pub rhomolar: f64,
    pub q: f64,
}

pub use rustprop_core::UPSTREAM_VERSION;

#[cfg(test)]
mod tests {
    use super::*;
    use rustprop_core::params::Param;
    use rustprop_heos::derivs::StateDerivs;

    fn eos(kind: CubicKind, name: &str) -> CubicEos {
        let data = rustprop_data::cubics::CUBIC_FLUIDS
            .iter()
            .find(|f| f.name == name)
            .expect("fluid in cubic table");
        CubicEos::new(kind, data)
    }

    #[track_caller]
    fn assert_rel(actual: f64, expected: f64, rtol: f64, what: &str) {
        let err = ((actual - expected) / expected).abs();
        assert!(
            err <= rtol,
            "{what}: got {actual:.17e}, wheel {expected:.17e}, rel err {err:.3e} > {rtol:.1e}"
        );
    }

    /// The new third-order alphar fields against central differences of the
    /// golden-verified second-order fields.
    #[test]
    fn third_order_alphar_matches_finite_differences() {
        for (eos, tau, delta) in [
            (eos(CubicKind::Srk, "WATER"), 1.0 / 300.0, 41892.19240933652),
            (
                eos(CubicKind::PengRobinson, "N-PROPANE"),
                1.0 / 320.0,
                202.19964739479565,
            ),
        ] {
            let d = eos.alphar_all(tau, delta);
            let hd = 1e-6 * delta;
            let ht = 1e-6 * tau;
            let dp = eos.alphar_all(tau, delta + hd);
            let dm = eos.alphar_all(tau, delta - hd);
            let tp = eos.alphar_all(tau + ht, delta);
            let tm = eos.alphar_all(tau - ht, delta);
            let rtol = 1e-6;
            assert_rel(
                d.d30,
                (dp.d20 - dm.d20) / (2.0 * hd),
                rtol,
                "d30 vs FD(d20)",
            );
            assert_rel(
                d.d21,
                (dp.d11 - dm.d11) / (2.0 * hd),
                rtol,
                "d21 vs FD(d11)",
            );
            assert_rel(
                d.d12,
                (tp.d11 - tm.d11) / (2.0 * ht),
                rtol,
                "d12 vs FD(d11)",
            );
            assert_rel(
                d.d03,
                (tp.d02 - tm.d02) / (2.0 * ht),
                rtol,
                "d03 vs FD(d02)",
            );
            // Cross-checks through the other variable.
            assert_rel(
                d.d21,
                (tp.d20 - tm.d20) / (2.0 * ht),
                rtol,
                "d21 vs FD(d20)",
            );
            assert_rel(
                d.d12,
                (dp.d02 - dm.d02) / (2.0 * hd),
                rtol,
                "d12 vs FD(d02)",
            );
        }
    }

    /// The cubic-convention alpha0 matrix against the wheel's keyed
    /// `dalpha0_*` outputs (SRK::Water T=300, P=101325; rho is the wheel's
    /// density root).
    #[test]
    fn alpha0_matrix_matches_wheel() {
        let eos = eos(CubicKind::Srk, "WATER");
        let (t, rho) = (300.0, 41892.19240933652);
        let a0 = eos.alpha0_all(1.0 / t, rho);
        let rtol = 1e-13;
        assert_rel(a0.d01, 5228.34514185985, rtol, "dalpha0_dtau_constdelta");
        assert_rel(
            a0.d10,
            2.3870796501381716e-05,
            rtol,
            "dalpha0_ddelta_consttau",
        );
        assert_rel(
            a0.d20,
            -5.698149256103775e-10,
            rtol,
            "d2alpha0_ddelta2_consttau",
        );
        assert_rel(
            a0.d30,
            2.7203872265390565e-14,
            rtol,
            "d3alpha0_ddelta3_consttau",
        );
    }

    /// The new property methods against wheel values (SRK::Water T=300,
    /// P=101325 from the spec's Part 2 table, plus a fresh PR::n-Propane
    /// state probed 2026-08-17).
    #[test]
    fn new_property_methods_match_wheel() {
        let srk = eos(CubicKind::Srk, "WATER");
        let (t, rho) = (300.0, 41892.19240933652);
        let rtol = 1e-12;
        assert_rel(
            srk.hmolar_idealgas(t, rho),
            45965.219022242796,
            rtol,
            "Hmolar_idealgas",
        );
        assert_rel(
            srk.smolar_idealgas(t, rho),
            67.93894818182274,
            rtol,
            "Smolar_idealgas",
        );
        assert_rel(
            srk.umolar_idealgas(t, rho),
            43470.880236796824,
            rtol,
            "Umolar_idealgas",
        );
        assert_rel(
            srk.gibbsmolar_residual(t, rho),
            -26394.231653381616,
            rtol,
            "Gmolar_residual",
        );
        assert_rel(
            srk.helmholtzmolar(t, rho),
            -813.1157941411406,
            rtol,
            "Helmholtzmolar",
        );
        assert_rel(srk.bvirial(t), -0.00037031381518632956, rtol, "Bvirial");
        assert_rel(srk.cvirial(t), 8.716342737432009e-09, rtol, "Cvirial");
        assert_rel(
            srk.dbvirial_dt(t),
            1.9788439201291845e-06,
            rtol,
            "dBvirial_dT",
        );
        assert_rel(
            srk.dcvirial_dt(t),
            -4.1807133419866773e-11,
            rtol,
            "dCvirial_dT",
        );

        let pr = eos(CubicKind::PengRobinson, "N-PROPANE");
        let (t, rho) = (320.0, 202.19964739479565);
        assert_rel(
            pr.hmolar_idealgas(t, rho),
            29570.873486797595,
            rtol,
            "PR Hmolar_idealgas",
        );
        assert_rel(
            pr.smolar_idealgas(t, rho),
            117.17530240282125,
            rtol,
            "PR Smolar_idealgas",
        );
        assert_rel(
            pr.umolar_idealgas(t, rho),
            26910.245448988557,
            rtol,
            "PR Umolar_idealgas",
        );
        assert_rel(
            pr.gibbsmolar_residual(t, rho),
            -378.2725735622879,
            rtol,
            "PR Gmolar_residual",
        );
        assert_rel(
            pr.helmholtzmolar(t, rho),
            -10776.29937655468,
            rtol,
            "PR Helmholtzmolar",
        );
        assert_rel(pr.bvirial(t), -0.00035896737773585385, rtol, "PR Bvirial");
        assert_rel(pr.cvirial(t), 4.9907523151543885e-08, rtol, "PR Cvirial");
        assert_rel(
            pr.dbvirial_dt(t),
            1.9959910596557732e-06,
            rtol,
            "PR dBvirial_dT",
        );
        assert_rel(
            pr.dcvirial_dt(t),
            -2.2466814828132958e-10,
            rtol,
            "PR dCvirial_dT",
        );
    }

    /// End-to-end proof of the generic (T, rho) machinery over the cubic:
    /// PIP, both compressibilities, the isentropic expansion coefficient and
    /// the fundamental derivative of gas dynamics against the wheel, exactly
    /// the AbstractState formula chains props_api will wire up.
    #[test]
    fn derivative_machinery_matches_wheel() {
        struct Case {
            eos: CubicEos,
            t: f64,
            p: f64,
            rho: f64,
            pip: f64,
            kappa: f64,
            beta: f64,
            isentropic: f64,
            fd: f64,
        }
        for c in [
            Case {
                eos: eos(CubicKind::Srk, "WATER"),
                t: 300.0,
                p: 101325.0,
                rho: 41892.19240933652,
                pip: 14.350546351600986,
                kappa: 1.534208236534884e-10,
                beta: 0.0007050575617414767,
                isentropic: 87819.91820691254,
                fd: 10.175578115135734,
            },
            Case {
                eos: eos(CubicKind::PengRobinson, "N-PROPANE"),
                t: 320.0,
                p: 500000.0,
                rho: 202.19964739479565,
                pip: 0.7929063097130047,
                kappa: 2.159523584210258e-06,
                beta: 0.003832759413225277,
                isentropic: 1.069334363260059,
                fd: 1.0101341740122103,
            },
        ] {
            let (rho, _phase) = c.eos.pt_flash(c.t, c.p).expect("PT flash");
            assert_rel(rho, c.rho, 1e-9, "density root");
            let d = StateDerivs::new(&c.eos, c.t, rho);
            let rtol = 1e-9;

            // calc_PIP (HelmholtzEOSMixtureBackend.h:603)
            let pip = 2.0
                - rho
                    * (d.second_partial_deriv(
                        Param::P,
                        Param::Dmolar,
                        Param::T,
                        Param::T,
                        Param::Dmolar,
                    )
                    .unwrap()
                        / d.first_partial_deriv(Param::P, Param::T, Param::Dmolar)
                            .unwrap()
                        - d.second_partial_deriv(
                            Param::P,
                            Param::Dmolar,
                            Param::T,
                            Param::Dmolar,
                            Param::T,
                        )
                        .unwrap()
                            / d.first_partial_deriv(Param::P, Param::Dmolar, Param::T)
                                .unwrap());
            assert_rel(pip, c.pip, rtol, "PIP");

            // calc_isothermal_compressibility (AbstractState.cpp:950)
            let kappa = d
                .first_partial_deriv(Param::Dmolar, Param::P, Param::T)
                .unwrap()
                / rho;
            assert_rel(kappa, c.kappa, rtol, "isothermal_compressibility");

            // calc_isobaric_expansion_coefficient (AbstractState.cpp:953)
            let beta = -d
                .first_partial_deriv(Param::Dmolar, Param::T, Param::P)
                .unwrap()
                / rho;
            assert_rel(beta, c.beta, rtol, "isobaric_expansion_coefficient");

            // calc_isentropic_expansion_coefficient (AbstractState.cpp:956)
            let isentropic = rho / c.p
                * d.first_partial_deriv(Param::P, Param::Dmolar, Param::Smolar)
                    .unwrap();
            assert_rel(
                isentropic,
                c.isentropic,
                rtol,
                "isentropic_expansion_coefficient",
            );

            // fundamental_derivative_of_gas_dynamics (AbstractState.cpp:975)
            let w = c.eos.speed_sound(c.t, rho);
            let rhomass = rho * c.eos.molar_mass();
            let fd = 1.0
                + d.second_partial_deriv(
                    Param::P,
                    Param::Dmass,
                    Param::Smolar,
                    Param::Dmass,
                    Param::Smolar,
                )
                .unwrap()
                    * rhomass
                    / (2.0 * w * w);
            assert_rel(fd, c.fd, rtol, "fundamental_derivative_of_gas_dynamics");
        }
    }
}
