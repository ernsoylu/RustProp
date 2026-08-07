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
use rustprop_heos::alpha::Alpha0Eval;
use rustprop_heos::solvers::{bounded_secant, secant, solve_cubic};

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
/// (tau = 1/T, delta = rho) variables.
#[derive(Default, Clone, Copy)]
pub struct CubicDerivs {
    pub d00: f64,
    pub d10: f64,
    pub d01: f64,
    pub d20: f64,
    pub d11: f64,
    pub d02: f64,
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

    /// `d^n/dtau^n (tau * a)` for n = 0..=2 (upstream `tau_times_a`,
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
                let a_term = ((delta * self.delta_1 * self.b + 1.0)
                    / (delta * self.delta_2 * self.b + 1.0))
                    .ln();
                a_term / self.b / (self.delta_1 - self.delta_2)
            }
            1 => 1.0 / self.pi_12(delta, 0),
            2 => {
                let p0 = self.pi_12(delta, 0);
                let p1 = self.pi_12(delta, 1);
                -p1 / (p0 * p0)
            }
            _ => unreachable!(),
        }
    }

    /// All residual derivatives through second order in the cubic's own
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
        let pm = [
            self.psi_minus(delta, 0),
            self.psi_minus(delta, 1),
            self.psi_minus(delta, 2),
        ];
        let pp = [
            self.psi_plus(delta, 0),
            self.psi_plus(delta, 1),
            self.psi_plus(delta, 2),
        ];
        CubicDerivs {
            d00: pm[0] - ta0 * pp[0],
            d10: pm[1] - ta0 * pp[1],
            d01: -ta1 * pp[0],
            d20: pm[2] - ta0 * pp[2],
            d11: -ta1 * pp[1],
            d02: -ta2 * pp[0],
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
        Ok((rho, self.singlephase_phase(t, p, rho)))
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
