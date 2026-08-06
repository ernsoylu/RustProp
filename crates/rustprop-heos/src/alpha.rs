//! Helmholtz-energy term evaluation (PLAN.md 4.1/4.7) — port of the term
//! machinery in `include/CoolProp/fluids/Helmholtz.h` + `src/Helmholtz.cpp`
//! @ v8.0.0 for the families the ported fluids use:
//!
//! - residual: `ResidualHelmholtzGeneralizedExponential` (Power and Gaussian
//!   JSON terms convert into it exactly as upstream `add_Power`/`add_Gaussian`
//!   do), `ResidualHelmholtzNonAnalytic`, and `ResidualHelmholtzGaoB`;
//! - ideal-gas: `IdealHelmholtzLead`, `IdealHelmholtzEnthalpyEntropyOffset`
//!   (the document's `EnthalpyEntropyOffsetCore` slot), `IdealHelmholtzLogTau`,
//!   `IdealHelmholtzPower`, and `IdealHelmholtzPlanckEinsteinGeneralized`
//!   (JSON `PlanckEinstein` maps with `theta = -t`, `c = 1`, `d = -1`; JSON
//!   `PlanckEinsteinFunctionT` with `theta = -v/Tcrit`, `c = 1`, `d = -1`,
//!   as in upstream FluidLibrary.h; both extend one merged container).
//!
//! Evaluation order mirrors the upstream containers: GenExp writes into a
//! fresh accumulator and applies its trailing 1/delta, 1/tau scalings while
//! only its own sums are present; NonAnalytic then GaoB add their
//! already-scaled contributions (`ResidualHelmholtzContainer::all` order).
//! Ideal terms evaluate in the container's fixed member order — Lead,
//! EnthalpyEntropyOffsetCore, LogTau, Power, PlanckEinstein — NOT JSON
//! document order; the sums share accumulator fields, so term order is part
//! of the floating-point result. Floating-point operation order is kept
//! line-for-line, including upstream's sequential-multiply `powInt`.

use rustprop_core::fluid::{Alpha0Term, AlpharTerm};

/// All alpha derivatives at one (tau, delta), mirroring upstream
/// `HelmholtzDerivatives`. `dij` holds d^(i+j) alpha / d delta^i d tau^j.
#[derive(Debug, Clone, Copy, Default)]
pub struct HelmholtzDerivs {
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
    pub d40: f64,
    pub d31: f64,
    pub d22: f64,
    pub d13: f64,
    pub d04: f64,
}

/// Upstream `powInt` (CPnumerics.cpp): sequential multiplication.
fn pow_int(x: f64, y: i32) -> f64 {
    if y == 0 {
        return 1.0;
    }
    let (x_in, y_in) = if y < 0 { (1.0 / x, -y) } else { (x, y) };
    if y_in == 1 {
        return x_in;
    }
    let mut product = x_in;
    for _ in 1..y_in {
        product *= x_in;
    }
    product
}

fn pow2(x: f64) -> f64 {
    x * x
}
fn pow3(x: f64) -> f64 {
    pow2(x) * x
}
fn pow4(x: f64) -> f64 {
    pow2(x) * pow2(x)
}

// ---------------------------------------------------------------------------
// Residual: generalized exponential
// ---------------------------------------------------------------------------

/// One flattened element (upstream
/// `ResidualHelmholtzGeneralizedExponentialElement`); fields default to the
/// same zeros as the upstream constructor.
#[derive(Clone, Default)]
struct GenExpElement {
    n: f64,
    d: f64,
    t: f64,
    c: f64,
    l_double: f64,
    l_int: i32,
    l_is_int: bool,
    omega: f64,
    m_double: f64,
    eta2: f64,
    epsilon2: f64,
    beta2: f64,
    gamma2: f64,
}

/// Upstream `ResidualHelmholtzGeneralizedExponential`, restricted to the
/// containers Water's families set (`delta_li_in_u`, `eta2_in_u`,
/// `beta2_in_u`); the omega/eta1/beta1 channels join when a ported fluid
/// needs them.
pub(crate) struct GenExp {
    elements: Vec<GenExpElement>,
    delta_li_in_u: bool,
    tau_mi_in_u: bool,
    eta2_in_u: bool,
    beta2_in_u: bool,
}

impl GenExp {
    fn new() -> Self {
        GenExp {
            elements: Vec::new(),
            delta_li_in_u: false,
            tau_mi_in_u: false,
            eta2_in_u: false,
            beta2_in_u: false,
        }
    }

    /// Upstream `add_Power`.
    fn add_power(&mut self, n: &[f64], d: &[f64], t: &[f64], l: &[f64]) {
        for i in 0..n.len() {
            let mut el = GenExpElement {
                n: n[i],
                d: d[i],
                t: t[i],
                l_double: l[i],
                l_int: l[i] as i32,
                ..Default::default()
            };
            el.c = if el.l_double > 0.0 { 1.0 } else { 0.0 };
            self.elements.push(el);
        }
        self.delta_li_in_u = true;
    }

    /// Upstream `add_Exponential`: the g coefficient rides the c channel.
    fn add_exponential(&mut self, n: &[f64], d: &[f64], t: &[f64], g: &[f64], l: &[f64]) {
        for i in 0..n.len() {
            self.elements.push(GenExpElement {
                n: n[i],
                d: d[i],
                t: t[i],
                c: g[i],
                l_double: l[i],
                l_int: l[i] as i32,
                ..Default::default()
            });
        }
        self.delta_li_in_u = true;
    }

    /// Upstream `add_Lemmon2005`: c = omega = 1, l and m exponents.
    fn add_lemmon2005(&mut self, n: &[f64], d: &[f64], t: &[f64], l: &[f64], m: &[f64]) {
        for i in 0..n.len() {
            self.elements.push(GenExpElement {
                n: n[i],
                d: d[i],
                t: t[i],
                c: 1.0,
                omega: 1.0,
                l_double: l[i],
                m_double: m[i],
                l_int: l[i] as i32,
                ..Default::default()
            });
        }
        self.delta_li_in_u = true;
        self.tau_mi_in_u = true;
    }

    /// Upstream `add_DoubleExponential`: gd/ld ride the c/l channel,
    /// gt/lt the omega/m channel (argument list mirrors the upstream
    /// signature).
    #[allow(clippy::too_many_arguments)]
    fn add_double_exponential(
        &mut self,
        n: &[f64],
        d: &[f64],
        t: &[f64],
        gd: &[f64],
        ld: &[f64],
        gt: &[f64],
        lt: &[f64],
    ) {
        for i in 0..n.len() {
            self.elements.push(GenExpElement {
                n: n[i],
                d: d[i],
                t: t[i],
                c: gd[i],
                l_double: ld[i],
                l_int: ld[i] as i32,
                omega: gt[i],
                m_double: lt[i],
                ..Default::default()
            });
        }
        self.delta_li_in_u = true;
        self.tau_mi_in_u = true;
    }

    /// Upstream `add_Gaussian` (argument list mirrors the upstream signature).
    #[allow(clippy::too_many_arguments)]
    fn add_gaussian(
        &mut self,
        n: &[f64],
        d: &[f64],
        t: &[f64],
        eta: &[f64],
        epsilon: &[f64],
        beta: &[f64],
        gamma: &[f64],
    ) {
        for i in 0..n.len() {
            self.elements.push(GenExpElement {
                n: n[i],
                d: d[i],
                t: t[i],
                eta2: eta[i],
                epsilon2: epsilon[i],
                beta2: beta[i],
                gamma2: gamma[i],
                ..Default::default()
            });
        }
        self.eta2_in_u = true;
        self.beta2_in_u = true;
    }

    /// Upstream `finish()` — set the l-is-integer fast-path flags.
    fn finish(&mut self) {
        for el in &mut self.elements {
            el.l_is_int = ((el.l_double as i64) as f64 - el.l_double).abs() < 1e-14;
        }
    }

    /// Upstream `ResidualHelmholtzGeneralizedExponential::all`. Writes into
    /// `derivs`, which MUST be freshly zeroed: the trailing scaling
    /// multiplies the accumulated fields.
    fn all(&self, tau: f64, delta: f64, derivs: &mut HelmholtzDerivs) {
        let log_tau = tau.ln();
        let log_delta = delta.ln();
        let one_over_delta = 1.0 / delta;
        let one_over_tau = 1.0 / tau;

        for el in &self.elements {
            let (ni, di, ti) = (el.n, el.d, el.t);

            let mut u = 0.0;
            let mut du_ddelta = 0.0;
            let mut du_dtau = 0.0;
            let mut d2u_ddelta2 = 0.0;
            let mut d2u_dtau2 = 0.0;
            let mut d3u_ddelta3 = 0.0;
            let mut d3u_dtau3 = 0.0;
            let mut d4u_ddelta4 = 0.0;
            let mut d4u_dtau4 = 0.0;

            if self.delta_li_in_u {
                let (ci, l_double) = (el.c, el.l_double);
                if l_double.is_finite() && l_double > 0.0 && ci.abs() > f64::EPSILON {
                    let u_increment = if el.l_is_int {
                        -ci * pow_int(delta, el.l_int)
                    } else {
                        -ci * (l_double * log_delta).exp()
                    };
                    let du_ddelta_increment = l_double * u_increment * one_over_delta;
                    let d2u_ddelta2_increment =
                        (l_double - 1.0) * du_ddelta_increment * one_over_delta;
                    let d3u_ddelta3_increment =
                        (l_double - 2.0) * d2u_ddelta2_increment * one_over_delta;
                    let d4u_ddelta4_increment =
                        (l_double - 3.0) * d3u_ddelta3_increment * one_over_delta;
                    u += u_increment;
                    du_ddelta += du_ddelta_increment;
                    d2u_ddelta2 += d2u_ddelta2_increment;
                    d3u_ddelta3 += d3u_ddelta3_increment;
                    d4u_ddelta4 += d4u_ddelta4_increment;
                }
            }
            if self.tau_mi_in_u {
                let (omegai, m_double) = (el.omega, el.m_double);
                if m_double.abs() > 0.0 {
                    let u_increment = -omegai * (m_double * log_tau).exp();
                    let du_dtau_increment = m_double * u_increment * one_over_tau;
                    let d2u_dtau2_increment = (m_double - 1.0) * du_dtau_increment * one_over_tau;
                    let d3u_dtau3_increment = (m_double - 2.0) * d2u_dtau2_increment * one_over_tau;
                    let d4u_dtau4_increment = (m_double - 3.0) * d3u_dtau3_increment * one_over_tau;
                    u += u_increment;
                    du_dtau += du_dtau_increment;
                    d2u_dtau2 += d2u_dtau2_increment;
                    d3u_dtau3 += d3u_dtau3_increment;
                    d4u_dtau4 += d4u_dtau4_increment;
                }
            }
            if self.eta2_in_u {
                let (eta2, epsilon2) = (el.eta2, el.epsilon2);
                if eta2.is_finite() {
                    u += -eta2 * pow2(delta - epsilon2);
                    du_ddelta += -2.0 * eta2 * (delta - epsilon2);
                    d2u_ddelta2 += -2.0 * eta2;
                }
            }
            if self.beta2_in_u {
                let (beta2, gamma2) = (el.beta2, el.gamma2);
                if beta2.is_finite() {
                    u += -beta2 * pow2(tau - gamma2);
                    du_dtau += -2.0 * beta2 * (tau - gamma2);
                    d2u_dtau2 += -2.0 * beta2;
                }
            }

            let ndteu = ni * (ti * log_tau + di * log_delta + u).exp();

            let db_delta_ddelta = delta * d2u_ddelta2 + du_ddelta;
            let d2b_delta_ddelta2 = delta * d3u_ddelta3 + 2.0 * d2u_ddelta2;
            let d3b_delta_ddelta3 = delta * d4u_ddelta4 + 3.0 * d3u_ddelta3;

            let b_delta = delta * du_ddelta + di;
            let b_delta2 = delta * db_delta_ddelta + (b_delta - 1.0) * b_delta;
            let db_delta2_ddelta = delta * d2b_delta_ddelta2 + 2.0 * b_delta * db_delta_ddelta;
            let b_delta3 = delta * db_delta2_ddelta + (b_delta - 2.0) * b_delta2;
            let db_delta3_ddelta = delta * delta * d3b_delta_ddelta3
                + 3.0 * delta * b_delta * d2b_delta_ddelta2
                + 3.0 * delta * pow2(db_delta_ddelta)
                + 3.0 * b_delta * (b_delta - 1.0) * db_delta_ddelta;
            let b_delta4 = delta * db_delta3_ddelta + (b_delta - 3.0) * b_delta3;

            let db_tau_dtau = tau * d2u_dtau2 + du_dtau;
            let d2b_tau_dtau2 = tau * d3u_dtau3 + 2.0 * d2u_dtau2;
            let d3b_tau_dtau3 = tau * d4u_dtau4 + 3.0 * d3u_dtau3;

            let b_tau = tau * du_dtau + ti;
            let b_tau2 = tau * db_tau_dtau + (b_tau - 1.0) * b_tau;
            let db_tau2_dtau = tau * d2b_tau_dtau2 + 2.0 * b_tau * db_tau_dtau;
            let b_tau3 = tau * db_tau2_dtau + (b_tau - 2.0) * b_tau2;
            let db_tau3_dtau = tau * tau * d3b_tau_dtau3
                + 3.0 * tau * b_tau * d2b_tau_dtau2
                + 3.0 * tau * pow2(db_tau_dtau)
                + 3.0 * b_tau * (b_tau - 1.0) * db_tau_dtau;
            let b_tau4 = tau * db_tau3_dtau + (b_tau - 3.0) * b_tau3;

            derivs.d00 += ndteu;

            derivs.d10 += ndteu * b_delta;
            derivs.d01 += ndteu * b_tau;

            derivs.d20 += ndteu * b_delta2;
            derivs.d11 += ndteu * b_delta * b_tau;
            derivs.d02 += ndteu * b_tau2;

            derivs.d30 += ndteu * b_delta3;
            derivs.d21 += ndteu * b_delta2 * b_tau;
            derivs.d12 += ndteu * b_delta * b_tau2;
            derivs.d03 += ndteu * b_tau3;

            derivs.d40 += ndteu * b_delta4;
            derivs.d31 += ndteu * b_delta3 * b_tau;
            derivs.d22 += ndteu * b_delta2 * b_tau2;
            derivs.d13 += ndteu * b_delta * b_tau3;
            derivs.d04 += ndteu * b_tau4;
        }
        derivs.d10 *= one_over_delta;
        derivs.d01 *= one_over_tau;
        derivs.d20 *= pow2(one_over_delta);
        derivs.d02 *= pow2(one_over_tau);
        derivs.d11 *= one_over_delta * one_over_tau;

        derivs.d30 *= pow3(one_over_delta);
        derivs.d03 *= pow3(one_over_tau);
        derivs.d21 *= pow2(one_over_delta) * one_over_tau;
        derivs.d12 *= one_over_delta * pow2(one_over_tau);

        derivs.d40 *= pow4(one_over_delta);
        derivs.d04 *= pow4(one_over_tau);
        derivs.d31 *= pow3(one_over_delta) * one_over_tau;
        derivs.d22 *= pow2(one_over_delta) * pow2(one_over_tau);
        derivs.d13 *= one_over_delta * pow3(one_over_tau);
    }
}

// ---------------------------------------------------------------------------
// Residual: non-analytic (IAPWS-95 critical terms)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct NonAnalyticElement {
    n: f64,
    a: f64,
    b: f64,
    beta: f64,
    big_a: f64,
    big_b: f64,
    big_c: f64,
    big_d: f64,
}

/// Upstream `ResidualHelmholtzNonAnalytic::all`, ported line-for-line
/// including the near-critical nudge away from tau = 1 / delta = 1.
#[allow(clippy::too_many_lines)]
fn non_analytic_all(
    elements: &[NonAnalyticElement],
    tau_in: f64,
    delta_in: f64,
    derivs: &mut HelmholtzDerivs,
) {
    if elements.is_empty() {
        return;
    }

    let mut tau = tau_in;
    let mut delta = delta_in;
    if (tau_in - 1.0).abs() < 10.0 * f64::EPSILON {
        tau = 1.0 + 10.0 * f64::EPSILON;
    }
    if (delta_in - 1.0).abs() < 10.0 * f64::EPSILON {
        delta = 1.0 + 10.0 * f64::EPSILON;
    }

    for el in elements {
        let (ni, ai, bi, betai) = (el.n, el.a, el.b, el.beta);
        let (big_a, big_b, big_c, big_d) = (el.big_a, el.big_b, el.big_c, el.big_d);

        // Derivatives of theta (all others are zero)
        let theta = (1.0 - tau) + big_a * pow2(delta - 1.0).powf(1.0 / (2.0 * betai));
        let dtheta_dtau = -1.0;
        let dtheta_ddelta =
            big_a / betai * pow2(delta - 1.0).powf(1.0 / (2.0 * betai) - 1.0) * (delta - 1.0);
        let d2theta_ddelta2 =
            big_a / betai * (1.0 / betai - 1.0) * pow2(delta - 1.0).powf(1.0 / (2.0 * betai) - 1.0);
        let d3theta_ddelta3 = big_a / betai
            * (2.0 - 3.0 / betai + 1.0 / pow2(betai))
            * pow2(delta - 1.0).powf(1.0 / (2.0 * betai))
            / pow3(delta - 1.0);
        let d4theta_ddelta4 = big_a / betai
            * (-6.0 + 11.0 / betai - 6.0 / pow2(betai) + 1.0 / pow3(betai))
            * pow2(delta - 1.0).powf(1.0 / (2.0 * betai) - 2.0);

        // Derivatives of PSI
        let psi = (-big_c * pow2(delta - 1.0) - big_d * pow2(tau - 1.0)).exp();
        let dpsi_ddelta_over_psi = -2.0 * big_c * (delta - 1.0);
        let dpsi_ddelta = dpsi_ddelta_over_psi * psi;
        let dpsi_dtau_over_psi = -2.0 * big_d * (tau - 1.0);
        let dpsi_dtau = dpsi_dtau_over_psi * psi;
        let d2psi_ddelta2_over_psi = (2.0 * big_c * pow2(delta - 1.0) - 1.0) * 2.0 * big_c;
        let d2psi_ddelta2 = d2psi_ddelta2_over_psi * psi;
        let d3psi_ddelta3 = 2.0
            * big_c
            * psi
            * (-4.0 * big_c * big_c * pow3(delta - 1.0) + 6.0 * big_c * (delta - 1.0));
        let d4psi_ddelta4 = 4.0
            * big_c
            * big_c
            * psi
            * (4.0 * big_c * big_c * pow4(delta - 1.0) - 12.0 * big_c * pow2(delta - 1.0) + 3.0);
        let d2psi_dtau2 = (2.0 * big_d * pow2(tau - 1.0) - 1.0) * 2.0 * big_d * psi;
        let d3psi_dtau3 = 2.0
            * big_d
            * psi
            * (-4.0 * big_d * big_d * pow3(tau - 1.0) + 6.0 * big_d * (tau - 1.0));
        let d4psi_dtau4 = 4.0
            * big_d
            * big_d
            * psi
            * (4.0 * big_d * big_d * pow4(tau - 1.0) - 12.0 * big_d * pow2(tau - 1.0) + 3.0);
        let d2psi_ddelta_dtau = dpsi_ddelta * dpsi_dtau_over_psi;
        let d3psi_ddelta2_dtau = d2psi_ddelta2 * dpsi_dtau_over_psi;
        let d3psi_ddelta_dtau2 = d2psi_dtau2 * dpsi_ddelta_over_psi;
        let d4psi_ddelta_dtau3 = d3psi_dtau3 * dpsi_ddelta_over_psi;
        let d4psi_ddelta2_dtau2 = d2psi_dtau2 * d2psi_ddelta2_over_psi;
        let d4psi_ddelta3_dtau = d3psi_ddelta3 * dpsi_dtau_over_psi;

        // Derivatives of DELTA
        let cap_delta = pow2(theta) + big_b * pow2(delta - 1.0).powf(ai);
        let ddelta_dtau = 2.0 * theta * dtheta_dtau;
        let ddelta_ddelta = 2.0 * theta * dtheta_ddelta
            + 2.0 * big_b * ai * pow2(delta - 1.0).powf(ai - 1.0) * (delta - 1.0);
        let d2delta_dtau2 = 2.0;
        let d2delta_ddelta_dtau = 2.0 * dtheta_dtau * dtheta_ddelta;
        let d2delta_ddelta2 = 2.0
            * (theta * d2theta_ddelta2
                + pow2(dtheta_ddelta)
                + big_b * (2.0 * ai * ai - ai) * pow2(delta - 1.0).powf(ai - 1.0));
        let d3delta_dtau3 = 0.0;
        let d3delta_ddelta_dtau2 = 0.0;
        let d3delta_ddelta2_dtau = 2.0 * dtheta_dtau * d2theta_ddelta2;
        let d3delta_ddelta3 = 2.0
            * (theta * d3theta_ddelta3
                + 3.0 * dtheta_ddelta * d2theta_ddelta2
                + 2.0
                    * big_b
                    * ai
                    * (2.0 * ai * ai - 3.0 * ai + 1.0)
                    * pow2(delta - 1.0).powf(ai - 1.0)
                    / (delta - 1.0));
        let d4delta_dtau4 = 0.0;
        let d4delta_ddelta_dtau3 = 0.0;
        let d4delta_ddelta2_dtau2 = 0.0;
        let d4delta_ddelta3_dtau = 2.0 * dtheta_dtau * d3theta_ddelta3;
        let d4delta_ddelta4 = 2.0
            * (theta * d4theta_ddelta4
                + 4.0 * dtheta_ddelta * d3theta_ddelta3
                + 3.0 * pow2(d2theta_ddelta2)
                + 2.0
                    * big_b
                    * ai
                    * (4.0 * ai * ai * ai - 12.0 * ai * ai + 11.0 * ai - 3.0)
                    * pow2(delta - 1.0).powf(ai - 2.0));

        let ddeltabi_ddelta = bi * cap_delta.powf(bi - 1.0) * ddelta_ddelta;
        let ddeltabi_dtau = -2.0 * theta * bi * cap_delta.powf(bi - 1.0);
        let d2deltabi_ddelta2 = bi
            * (cap_delta.powf(bi - 1.0) * d2delta_ddelta2
                + (bi - 1.0) * cap_delta.powf(bi - 2.0) * ddelta_ddelta.powf(2.0));
        let d3deltabi_ddelta3 = bi
            * (cap_delta.powf(bi - 1.0) * d3delta_ddelta3
                + d2delta_ddelta2 * (bi - 1.0) * cap_delta.powf(bi - 2.0) * ddelta_ddelta
                + (bi - 1.0)
                    * (cap_delta.powf(bi - 2.0) * 2.0 * ddelta_ddelta * d2delta_ddelta2
                        + ddelta_ddelta.powf(2.0)
                            * (bi - 2.0)
                            * cap_delta.powf(bi - 3.0)
                            * ddelta_ddelta));
        let d2deltabi_ddelta_dtau = -big_a * bi * 2.0 / betai
            * cap_delta.powf(bi - 1.0)
            * (delta - 1.0)
            * (delta - 1.0).powf(2.0).powf(1.0 / (2.0 * betai) - 1.0)
            - 2.0 * theta * bi * (bi - 1.0) * cap_delta.powf(bi - 2.0) * ddelta_ddelta;
        let d2deltabi_dtau2 = 2.0 * bi * cap_delta.powf(bi - 1.0)
            + 4.0 * theta.powf(2.0) * bi * (bi - 1.0) * cap_delta.powf(bi - 2.0);
        let d3deltabi_dtau3 = -12.0 * theta * bi * (bi - 1.0) * cap_delta.powf(bi - 2.0)
            - 8.0 * theta.powf(3.0) * bi * (bi - 1.0) * (bi - 2.0) * cap_delta.powf(bi - 3.0);
        let d3deltabi_ddelta_dtau2 =
            2.0 * bi * (bi - 1.0) * cap_delta.powf(bi - 2.0) * ddelta_ddelta
                + 4.0
                    * theta.powf(2.0)
                    * bi
                    * (bi - 1.0)
                    * (bi - 2.0)
                    * cap_delta.powf(bi - 3.0)
                    * ddelta_ddelta
                + 8.0 * theta * bi * (bi - 1.0) * cap_delta.powf(bi - 2.0) * dtheta_ddelta;
        let d3deltabi_ddelta2_dtau = bi
            * ((bi - 1.0) * cap_delta.powf(bi - 2.0) * ddelta_dtau * d2delta_ddelta2
                + cap_delta.powf(bi - 1.0) * d3delta_ddelta2_dtau
                + (bi - 1.0)
                    * ((bi - 2.0)
                        * cap_delta.powf(bi - 3.0)
                        * ddelta_dtau
                        * ddelta_ddelta.powf(2.0)
                        + cap_delta.powf(bi - 2.0) * 2.0 * ddelta_ddelta * d2delta_ddelta_dtau));

        // Fourth partials
        let delta_bi = cap_delta.powf(bi);
        let d4deltabi_dtau4 = bi * delta_bi / cap_delta
            * ((pow3(bi) - 6.0 * pow2(bi) + 11.0 * bi - 6.0) * pow4(ddelta_dtau) / pow3(cap_delta)
                + 6.0 * (bi * bi - 3.0 * bi + 2.0) * pow2(ddelta_dtau / cap_delta) * d2delta_dtau2
                + 4.0 * (bi - 1.0) * ddelta_dtau / cap_delta * d3delta_dtau3
                + 3.0 * (bi - 1.0) * pow2(d2delta_dtau2) / cap_delta
                + d4delta_dtau4);
        let d4deltabi_ddelta4 = bi * delta_bi / cap_delta
            * ((pow3(bi) - 6.0 * pow2(bi) + 11.0 * bi - 6.0) * pow4(ddelta_ddelta)
                / pow3(cap_delta)
                + 6.0
                    * (bi * bi - 3.0 * bi + 2.0)
                    * pow2(ddelta_ddelta / cap_delta)
                    * d2delta_ddelta2
                + 4.0 * (bi - 1.0) * ddelta_ddelta / cap_delta * d3delta_ddelta3
                + 3.0 * (bi - 1.0) * pow2(d2delta_ddelta2) / cap_delta
                + d4delta_ddelta4);
        let d4deltabi_ddelta_dtau3 = bi * (bi - 1.0) * delta_bi / pow2(cap_delta)
            * ddelta_ddelta
            * ((bi - 1.0) * (bi - 2.0) * pow3(ddelta_dtau) / pow2(cap_delta)
                + 3.0 * (bi - 1.0) * ddelta_dtau / cap_delta * d2delta_dtau2
                + d3delta_dtau3)
            + bi * delta_bi / cap_delta
                * ((bi - 1.0)
                    * (bi - 2.0)
                    * (3.0 * pow2(ddelta_dtau) / pow2(cap_delta) * d2delta_ddelta_dtau
                        - 2.0 * pow3(ddelta_dtau) / pow3(cap_delta) * ddelta_ddelta)
                    + 3.0
                        * (bi - 1.0)
                        * (ddelta_dtau / cap_delta * d3delta_ddelta_dtau2
                            + d2delta_dtau2 / cap_delta * d2delta_ddelta_dtau
                            - ddelta_ddelta / pow2(cap_delta) * ddelta_dtau * d2delta_dtau2)
                    + d4delta_ddelta_dtau3);
        let d4deltabi_ddelta3_dtau = bi * (bi - 1.0) * delta_bi / pow2(cap_delta)
            * ddelta_dtau
            * ((bi - 1.0) * (bi - 2.0) * pow3(ddelta_ddelta) / pow2(cap_delta)
                + 3.0 * (bi - 1.0) * ddelta_ddelta / cap_delta * d2delta_ddelta2
                + d3delta_ddelta3)
            + bi * delta_bi / cap_delta
                * ((bi - 1.0)
                    * (bi - 2.0)
                    * (3.0 * pow2(ddelta_ddelta) / pow2(cap_delta) * d2delta_ddelta_dtau
                        - 2.0 * pow3(ddelta_ddelta) / pow3(cap_delta) * ddelta_dtau)
                    + 3.0
                        * (bi - 1.0)
                        * (ddelta_ddelta / cap_delta * d3delta_ddelta2_dtau
                            + d2delta_ddelta2 / cap_delta * d2delta_ddelta_dtau
                            - ddelta_dtau / pow2(cap_delta) * ddelta_ddelta * d2delta_ddelta2)
                    + d4delta_ddelta3_dtau);
        let d4deltabi_ddelta2_dtau2 = bi * delta_bi / pow4(cap_delta)
            * ((pow3(bi) - 6.0 * bi * bi + 11.0 * bi - 6.0)
                * pow2(ddelta_ddelta)
                * pow2(ddelta_dtau)
                + (bi - 1.0) * (bi - 2.0) * cap_delta * pow2(ddelta_ddelta) * d2delta_dtau2
                + (bi - 1.0) * (bi - 2.0) * cap_delta * pow2(ddelta_dtau) * d2delta_ddelta2
                + 4.0
                    * (bi - 1.0)
                    * (bi - 2.0)
                    * cap_delta
                    * ddelta_ddelta
                    * ddelta_dtau
                    * d2delta_ddelta_dtau
                + 2.0 * (bi - 1.0) * pow2(cap_delta * d2delta_ddelta_dtau)
                + 2.0 * (bi - 1.0) * pow2(cap_delta) * ddelta_dtau * d3delta_ddelta2_dtau
                + 2.0 * (bi - 1.0) * pow2(cap_delta) * ddelta_ddelta * d3delta_ddelta_dtau2
                + (bi - 1.0) * pow2(cap_delta) * d2delta_ddelta2 * d2delta_dtau2
                + pow3(cap_delta) * d4delta_ddelta2_dtau2);

        derivs.d00 += delta * ni * delta_bi * psi;

        // First partials
        derivs.d01 += ni * delta * (delta_bi * dpsi_dtau + ddeltabi_dtau * psi);
        derivs.d10 += ni * (delta_bi * (psi + delta * dpsi_ddelta) + ddeltabi_ddelta * delta * psi);

        // Second partials
        derivs.d02 += ni
            * delta
            * (d2deltabi_dtau2 * psi + 2.0 * ddeltabi_dtau * dpsi_dtau + delta_bi * d2psi_dtau2);
        derivs.d11 += ni
            * (delta_bi * (dpsi_dtau + delta * d2psi_ddelta_dtau)
                + delta * ddeltabi_ddelta * dpsi_dtau
                + ddeltabi_dtau * (psi + delta * dpsi_ddelta)
                + d2deltabi_ddelta_dtau * delta * psi);
        derivs.d20 += ni
            * (delta_bi * (2.0 * dpsi_ddelta + delta * d2psi_ddelta2)
                + 2.0 * ddeltabi_ddelta * (psi + delta * dpsi_ddelta)
                + d2deltabi_ddelta2 * delta * psi);

        // Third partials
        derivs.d03 += ni
            * delta
            * (d3deltabi_dtau3 * psi
                + 3.0 * d2deltabi_dtau2 * dpsi_dtau
                + 3.0 * ddeltabi_dtau * d2psi_dtau2
                + delta_bi * d3psi_dtau3);
        derivs.d12 += ni
            * delta
            * (d2deltabi_dtau2 * dpsi_ddelta
                + d3deltabi_ddelta_dtau2 * psi
                + 2.0 * ddeltabi_dtau * d2psi_ddelta_dtau
                + 2.0 * d2deltabi_ddelta_dtau * dpsi_dtau
                + delta_bi * d3psi_ddelta_dtau2
                + ddeltabi_ddelta * d2psi_dtau2)
            + ni * (d2deltabi_dtau2 * psi
                + 2.0 * ddeltabi_dtau * dpsi_dtau
                + delta_bi * d2psi_dtau2);
        derivs.d30 += ni
            * (delta_bi * (3.0 * d2psi_ddelta2 + delta * d3psi_ddelta3)
                + 3.0 * ddeltabi_ddelta * (2.0 * dpsi_ddelta + delta * d2psi_ddelta2)
                + 3.0 * d2deltabi_ddelta2 * (psi + delta * dpsi_ddelta)
                + d3deltabi_ddelta3 * psi * delta);
        let line1 = delta_bi * (2.0 * d2psi_ddelta_dtau + delta * d3psi_ddelta2_dtau)
            + ddeltabi_dtau * (2.0 * dpsi_ddelta + delta * d2psi_ddelta2);
        let line2 = 2.0 * ddeltabi_ddelta * (dpsi_dtau + delta * d2psi_ddelta_dtau)
            + 2.0 * d2deltabi_ddelta_dtau * (psi + delta * dpsi_ddelta);
        let line3 = d2deltabi_ddelta2 * delta * dpsi_dtau + d3deltabi_ddelta2_dtau * delta * psi;
        derivs.d21 += ni * (line1 + line2 + line3);

        // Fourth partials
        derivs.d04 += ni
            * delta
            * (delta_bi * d4psi_dtau4
                + 4.0 * ddeltabi_dtau * d3psi_dtau3
                + 6.0 * d2deltabi_dtau2 * d2psi_dtau2
                + 4.0 * d3deltabi_dtau3 * dpsi_dtau
                + psi * d4deltabi_dtau4);
        derivs.d40 += ni
            * (delta * delta_bi * d4psi_ddelta4
                + delta * psi * d4deltabi_ddelta4
                + 4.0 * delta * ddeltabi_ddelta * d3psi_ddelta3
                + 4.0 * delta * dpsi_ddelta * d3deltabi_ddelta3
                + 6.0 * delta * d2deltabi_ddelta2 * d2psi_ddelta2
                + 4.0 * delta_bi * d3psi_ddelta3
                + 4.0 * psi * d3deltabi_ddelta3
                + 12.0 * ddeltabi_ddelta * d2psi_ddelta2
                + 12.0 * dpsi_ddelta * d2deltabi_ddelta2);
        derivs.d13 += ni
            * (delta * delta_bi * d4psi_ddelta_dtau3
                + delta * psi * d4deltabi_ddelta_dtau3
                + delta * ddeltabi_ddelta * d3psi_dtau3
                + 3.0 * delta * ddeltabi_dtau * d3psi_ddelta_dtau2
                + delta * dpsi_ddelta * d3deltabi_dtau3
                + 3.0 * delta * dpsi_dtau * d3deltabi_ddelta_dtau2
                + 3.0 * delta * d2deltabi_ddelta_dtau * d2psi_dtau2
                + 3.0 * delta * d2deltabi_dtau2 * d2psi_ddelta_dtau
                + delta_bi * d3psi_dtau3
                + psi * d3deltabi_dtau3
                + 3.0 * ddeltabi_dtau * d2psi_dtau2
                + 3.0 * dpsi_dtau * d2deltabi_dtau2);
        derivs.d31 += ni
            * (delta * delta_bi * d4psi_ddelta3_dtau
                + delta * psi * d4deltabi_ddelta3_dtau
                + 3.0 * delta * ddeltabi_ddelta * d3psi_ddelta2_dtau
                + delta * ddeltabi_dtau * d3psi_ddelta3
                + 3.0 * delta * dpsi_ddelta * d3deltabi_ddelta2_dtau
                + delta * dpsi_dtau * d3deltabi_ddelta3
                + 3.0 * delta * d2deltabi_ddelta2 * d2psi_ddelta_dtau
                + 3.0 * delta * d2deltabi_ddelta_dtau * d2psi_ddelta2
                + 3.0 * delta_bi * d3psi_ddelta2_dtau
                + 3.0 * psi * d3deltabi_ddelta2_dtau
                + 6.0 * ddeltabi_ddelta * d2psi_ddelta_dtau
                + 3.0 * ddeltabi_dtau * d2psi_ddelta2
                + 6.0 * dpsi_ddelta * d2deltabi_ddelta_dtau
                + 3.0 * dpsi_dtau * d2deltabi_ddelta2);
        derivs.d22 += ni
            * (delta * delta_bi * d4psi_ddelta2_dtau2
                + delta * psi * d4deltabi_ddelta2_dtau2
                + 2.0 * delta * ddeltabi_ddelta * d3psi_ddelta_dtau2
                + 2.0 * delta * ddeltabi_dtau * d3psi_ddelta2_dtau
                + 2.0 * delta * dpsi_ddelta * d3deltabi_ddelta_dtau2
                + 2.0 * delta * dpsi_dtau * d3deltabi_ddelta2_dtau
                + delta * d2deltabi_ddelta2 * d2psi_dtau2
                + 4.0 * delta * d2deltabi_ddelta_dtau * d2psi_ddelta_dtau
                + delta * d2deltabi_dtau2 * d2psi_ddelta2
                + 2.0 * delta_bi * d3psi_ddelta_dtau2
                + 2.0 * psi * d3deltabi_ddelta_dtau2
                + 2.0 * ddeltabi_ddelta * d2psi_dtau2
                + 4.0 * ddeltabi_dtau * d2psi_ddelta_dtau
                + 2.0 * dpsi_ddelta * d2deltabi_dtau2
                + 4.0 * dpsi_dtau * d2deltabi_ddelta_dtau);
    }
}

// ---------------------------------------------------------------------------
// Residual: GaoB (Gao et al. modified-Gaussian terms)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct GaoBElement {
    n: f64,
    t: f64,
    d: f64,
    eta: f64,
    beta: f64,
    gamma: f64,
    epsilon: f64,
    b: f64,
}

/// Upstream `ResidualHelmholtzGaoB::all`, ported line-for-line (the tau/delta
/// factor derivatives `tau^k d^k F_tau/d tau^k` etc. keep upstream's exact
/// expression trees, `pow` maps to `powf`, `POWn` macros to the `pown`
/// helpers).
#[allow(clippy::too_many_lines)]
fn gaob_all(elements: &[GaoBElement], tau: f64, delta: f64, derivs: &mut HelmholtzDerivs) {
    for el in elements {
        let (n, t, d, eta, beta, gamma, epsilon, b) = (
            el.n, el.t, el.d, el.eta, el.beta, el.gamma, el.epsilon, el.b,
        );

        let ftau = tau.powf(t) * (1.0 / (b + beta * (-gamma + tau).powf(2.0))).exp();
        let fdelta = delta.powf(d) * (eta * (delta - epsilon).powf(2.0)).exp();
        let taudftaudtau = (2.0 * beta * tau.powf(t + 1.0) * (gamma - tau)
            + t * tau.powf(t) * (b + beta * (gamma - tau).powf(2.0)).powf(2.0))
            * (1.0 / (b + beta * (gamma - tau).powf(2.0))).exp()
            / (b + beta * (gamma - tau).powf(2.0)).powf(2.0);
        let tau2d2ftaudtau2 = tau.powf(t)
            * (4.0
                * beta
                * t
                * tau
                * (b + beta * (gamma - tau).powf(2.0)).powf(2.0)
                * (gamma - tau)
                + 2.0
                    * beta
                    * tau.powf(2.0)
                    * (4.0
                        * beta
                        * (b + beta * (gamma - tau).powf(2.0))
                        * (gamma - tau).powf(2.0)
                        + 2.0 * beta * (gamma - tau).powf(2.0)
                        - (b + beta * (gamma - tau).powf(2.0)).powf(2.0))
                + t * (b + beta * (gamma - tau).powf(2.0)).powf(4.0) * (t - 1.0))
            * (1.0 / (b + beta * (gamma - tau).powf(2.0))).exp()
            / (b + beta * (gamma - tau).powf(2.0)).powf(4.0);
        let tau3d3ftaudtau3 = tau.powf(t)
            * (4.0
                * beta.powf(2.0)
                * tau.powf(3.0)
                * (gamma - tau)
                * (12.0 * beta * (b + beta * (gamma - tau).powf(2.0)) * (gamma - tau).powf(2.0)
                    + 2.0 * beta * (gamma - tau).powf(2.0)
                    - 6.0 * (b + beta * (gamma - tau).powf(2.0)).powf(3.0)
                    + (b + beta * (gamma - tau).powf(2.0)).powf(2.0)
                        * (12.0 * beta * (gamma - tau).powf(2.0) - 3.0))
                + 6.0
                    * beta
                    * t
                    * tau.powf(2.0)
                    * (b + beta * (gamma - tau).powf(2.0)).powf(2.0)
                    * (4.0
                        * beta
                        * (b + beta * (gamma - tau).powf(2.0))
                        * (gamma - tau).powf(2.0)
                        + 2.0 * beta * (gamma - tau).powf(2.0)
                        - (b + beta * (gamma - tau).powf(2.0)).powf(2.0))
                + 6.0
                    * beta
                    * t
                    * tau
                    * (b + beta * (gamma - tau).powf(2.0)).powf(4.0)
                    * (gamma - tau)
                    * (t - 1.0)
                + t * (b + beta * (gamma - tau).powf(2.0)).powf(6.0)
                    * (t.powf(2.0) - 3.0 * t + 2.0))
            * (1.0 / (b + beta * (gamma - tau).powf(2.0))).exp()
            / (b + beta * (gamma - tau).powf(2.0)).powf(6.0);
        let tau4d4ftaudtau4 = tau.powf(t)
            * (16.0
                * beta.powf(2.0)
                * t
                * tau.powf(3.0)
                * (b + beta * (gamma - tau).powf(2.0)).powf(2.0)
                * (gamma - tau)
                * (12.0 * beta * (b + beta * (gamma - tau).powf(2.0)) * (gamma - tau).powf(2.0)
                    + 2.0 * beta * (gamma - tau).powf(2.0)
                    - 6.0 * (b + beta * (gamma - tau).powf(2.0)).powf(3.0)
                    + (b + beta * (gamma - tau).powf(2.0)).powf(2.0)
                        * (12.0 * beta * (gamma - tau).powf(2.0) - 3.0))
                + beta.powf(2.0)
                    * tau.powf(4.0)
                    * (beta.powf(2.0)
                        * (192.0 * b + 192.0 * beta * (gamma - tau).powf(2.0))
                        * (gamma - tau).powf(4.0)
                        + 16.0 * beta.powf(2.0) * (gamma - tau).powf(4.0)
                        + 96.0
                            * beta
                            * (b + beta * (gamma - tau).powf(2.0)).powf(3.0)
                            * (gamma - tau).powf(2.0)
                            * (4.0 * beta * (gamma - tau).powf(2.0) - 3.0)
                        + 48.0
                            * beta
                            * (b + beta * (gamma - tau).powf(2.0)).powf(2.0)
                            * (gamma - tau).powf(2.0)
                            * (12.0 * beta * (gamma - tau).powf(2.0) - 1.0)
                        + 24.0 * (b + beta * (gamma - tau).powf(2.0)).powf(5.0)
                        + (b + beta * (gamma - tau).powf(2.0)).powf(4.0)
                            * (-288.0 * beta * (gamma - tau).powf(2.0) + 12.0))
                + 12.0
                    * beta
                    * t
                    * tau.powf(2.0)
                    * (b + beta * (gamma - tau).powf(2.0)).powf(4.0)
                    * (t - 1.0)
                    * (4.0
                        * beta
                        * (b + beta * (gamma - tau).powf(2.0))
                        * (gamma - tau).powf(2.0)
                        + 2.0 * beta * (gamma - tau).powf(2.0)
                        - (b + beta * (gamma - tau).powf(2.0)).powf(2.0))
                + 8.0
                    * beta
                    * t
                    * tau
                    * (b + beta * (gamma - tau).powf(2.0)).powf(6.0)
                    * (gamma - tau)
                    * (t.powf(2.0) - 3.0 * t + 2.0)
                + t * (b + beta * (gamma - tau).powf(2.0)).powf(8.0)
                    * (t.powf(3.0) - 6.0 * t.powf(2.0) + 11.0 * t - 6.0))
            * (1.0 / (b + beta * (gamma - tau).powf(2.0))).exp()
            / (b + beta * (gamma - tau).powf(2.0)).powf(8.0);
        let deltadfdeltaddelta = (d * delta.powf(d)
            + 2.0 * delta.powf(d + 1.0) * eta * (delta - epsilon))
            * (eta * (delta - epsilon).powf(2.0)).exp();
        let delta2d2fdeltaddelta2 = delta.powf(d)
            * (4.0 * d * delta * eta * (delta - epsilon)
                + d * (d - 1.0)
                + 2.0 * delta.powf(2.0) * eta * (2.0 * eta * (delta - epsilon).powf(2.0) + 1.0))
            * (eta * (delta - epsilon).powf(2.0)).exp();
        let delta3d3fdeltaddelta3 = delta.powf(d)
            * (6.0 * d * delta.powf(2.0) * eta * (2.0 * eta * (delta - epsilon).powf(2.0) + 1.0)
                + 6.0 * d * delta * eta * (d - 1.0) * (delta - epsilon)
                + d * (d.powf(2.0) - 3.0 * d + 2.0)
                + 4.0
                    * delta.powf(3.0)
                    * eta.powf(2.0)
                    * (delta - epsilon)
                    * (2.0 * eta * (delta - epsilon).powf(2.0) + 3.0))
            * (eta * (delta - epsilon).powf(2.0)).exp();
        let delta4d4fdeltaddelta4 = delta.powf(d)
            * (16.0
                * d
                * delta.powf(3.0)
                * eta.powf(2.0)
                * (delta - epsilon)
                * (2.0 * eta * (delta - epsilon).powf(2.0) + 3.0)
                + 12.0
                    * d
                    * delta.powf(2.0)
                    * eta
                    * (d - 1.0)
                    * (2.0 * eta * (delta - epsilon).powf(2.0) + 1.0)
                + 8.0 * d * delta * eta * (delta - epsilon) * (d.powf(2.0) - 3.0 * d + 2.0)
                + d * (d.powf(3.0) - 6.0 * d.powf(2.0) + 11.0 * d - 6.0)
                + delta.powf(4.0)
                    * eta.powf(2.0)
                    * (16.0 * eta.powf(2.0) * (delta - epsilon).powf(4.0)
                        + 48.0 * eta * (delta - epsilon).powf(2.0)
                        + 12.0))
            * (eta * (delta - epsilon).powf(2.0)).exp();

        derivs.d00 += n * ftau * fdelta;

        derivs.d10 += n * ftau * deltadfdeltaddelta / delta;
        derivs.d01 += n * fdelta * taudftaudtau / tau;

        derivs.d20 += n * ftau * delta2d2fdeltaddelta2 / pow2(delta);
        derivs.d11 += n * taudftaudtau * deltadfdeltaddelta / tau / delta;
        derivs.d02 += n * fdelta * tau2d2ftaudtau2 / pow2(tau);

        derivs.d30 += n * ftau * delta3d3fdeltaddelta3 / pow3(delta);
        derivs.d21 += n * taudftaudtau * delta2d2fdeltaddelta2 / pow2(delta) / tau;
        derivs.d12 += n * tau2d2ftaudtau2 * deltadfdeltaddelta / pow2(tau) / delta;
        derivs.d03 += n * fdelta * tau3d3ftaudtau3 / pow3(tau);

        derivs.d40 += n * ftau * delta4d4fdeltaddelta4 / pow4(delta);
        derivs.d31 += n * taudftaudtau * delta3d3fdeltaddelta3 / pow3(delta) / tau;
        derivs.d22 += n * tau2d2ftaudtau2 * delta2d2fdeltaddelta2 / pow2(delta) / pow2(tau);
        derivs.d13 += n * tau3d3ftaudtau3 * deltadfdeltaddelta / pow3(tau) / delta;
        derivs.d04 += n * fdelta * tau4d4ftaudtau4 / pow4(tau);
    }
}

// ---------------------------------------------------------------------------
// Ideal-gas terms
// ---------------------------------------------------------------------------

/// Contents of the single merged `IdealHelmholtzPlanckEinsteinGeneralized`
/// slot; both `PlanckEinstein` and `PlanckEinsteinFunctionT` JSON terms land
/// here in document order, as upstream's parse-time `extend`.
#[derive(Default)]
struct PlanckEinsteinGen {
    n: Vec<f64>,
    theta: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

/// Upstream `IdealHelmholtzContainer`, restricted to the slots the ported
/// fluids use. `all()` evaluates in the container's fixed member order (see
/// module docs) — sums share accumulator fields, so this order is part of
/// the floating-point result.
#[derive(Default)]
struct IdealContainer {
    /// `IdealHelmholtzLead`: `a1 + a2*tau + ln(delta)`
    lead: Option<(f64, f64)>,
    /// `IdealHelmholtzEnthalpyEntropyOffset` (the Core slot): `a1 + a2*tau`
    offset_core: Option<(f64, f64)>,
    /// `IdealHelmholtzLogTau`: `a1*ln(tau)`
    log_tau: Option<f64>,
    /// `IdealHelmholtzPower`: `sum n_i*tau^t_i` — `(n, t)`
    power: Option<(Vec<f64>, Vec<f64>)>,
    /// `IdealHelmholtzPlanckEinsteinGeneralized`:
    /// `sum n_i*ln(c_i + d_i*exp(theta_i*tau))`
    planck_einstein: Option<PlanckEinsteinGen>,
    /// `IdealHelmholtzCP0Constant` — `(cp_over_r, tc, t0)`
    cp0_constant: Option<(f64, f64, f64)>,
    /// `IdealHelmholtzCP0PolyT` — `(c, t, tc, t0)`
    cp0_polyt: Option<(Vec<f64>, Vec<f64>, f64, f64)>,
}

impl IdealContainer {
    fn all(&self, tau: f64, delta: f64, derivs: &mut HelmholtzDerivs) {
        if let Some((a1, a2)) = self.lead {
            derivs.d00 += delta.ln() + a1 + a2 * tau;
            derivs.d10 += 1.0 / delta;
            derivs.d01 += a2;
            derivs.d20 += -1.0 / delta / delta;
            derivs.d30 += 2.0 / delta / delta / delta;
            derivs.d40 += -6.0 / pow4(delta);
        }
        if let Some((a1, a2)) = self.offset_core {
            derivs.d00 += a1 + a2 * tau;
            derivs.d01 += a2;
        }
        if let Some(a1) = self.log_tau {
            derivs.d00 += a1 * tau.ln();
            derivs.d01 += a1 / tau;
            derivs.d02 += -a1 / tau / tau;
            derivs.d03 += 2.0 * a1 / tau / tau / tau;
            derivs.d04 += -6.0 * a1 / pow4(tau);
        }
        if let Some((n, t)) = &self.power {
            // Upstream `IdealHelmholtzPower::all`: one fresh sum per
            // derivative order.
            let mut s = 0.0;
            for i in 0..n.len() {
                s += n[i] * tau.powf(t[i]);
            }
            derivs.d00 += s;
            let mut s = 0.0;
            for i in 0..n.len() {
                s += n[i] * t[i] * tau.powf(t[i] - 1.0);
            }
            derivs.d01 += s;
            let mut s = 0.0;
            for i in 0..n.len() {
                s += n[i] * t[i] * (t[i] - 1.0) * tau.powf(t[i] - 2.0);
            }
            derivs.d02 += s;
            let mut s = 0.0;
            for i in 0..n.len() {
                s += n[i] * t[i] * (t[i] - 1.0) * (t[i] - 2.0) * tau.powf(t[i] - 3.0);
            }
            derivs.d03 += s;
            let mut s = 0.0;
            for i in 0..n.len() {
                s +=
                    n[i] * t[i] * (t[i] - 1.0) * (t[i] - 2.0) * (t[i] - 3.0) * tau.powf(t[i] - 4.0);
            }
            derivs.d04 += s;
        }
        if let Some(pe) = &self.planck_einstein {
            let (n, theta, c, d) = (&pe.n, &pe.theta, &pe.c, &pe.d);
            let (mut s00, mut s01, mut s02, mut s03, mut s04) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for i in 0..n.len() {
                let expthetataui = (theta[i] * tau).exp();
                let para = c[i] + d[i] * expthetataui;
                s00 += n[i] * para.ln();
                s01 += n[i] * theta[i] * d[i] * expthetataui / para;
                s02 += n[i] * pow2(theta[i]) * c[i] * d[i] * expthetataui / pow2(para);
                s03 += n[i]
                    * pow3(theta[i])
                    * c[i]
                    * d[i]
                    * (c[i] - d[i] * expthetataui)
                    * expthetataui
                    / pow3(para);
                let bracket = 6.0 * pow3(d[i]) * pow3(expthetataui)
                    - 12.0 * d[i] * d[i] * para * pow2(expthetataui)
                    + 7.0 * d[i] * pow2(para) * expthetataui
                    - pow3(para);
                s04 += -n[i] * d[i] * pow4(theta[i]) * bracket * expthetataui / pow4(para);
            }
            derivs.d00 += s00;
            derivs.d01 += s01;
            derivs.d02 += s02;
            derivs.d03 += s03;
            derivs.d04 += s04;
        }
        if let Some((cp_over_r, tc, t0)) = self.cp0_constant {
            // Upstream `IdealHelmholtzCP0Constant::all` (tau0 = Tc/T0).
            let tau0 = tc / t0;
            derivs.d00 += cp_over_r - cp_over_r * tau / tau0 + cp_over_r * (tau / tau0).ln();
            derivs.d01 += cp_over_r / tau - cp_over_r / tau0;
            derivs.d02 += -cp_over_r / (tau * tau);
            derivs.d03 += 2.0 * cp_over_r / (tau * tau * tau);
            derivs.d04 += -6.0 * cp_over_r / pow4(tau);
        }
        if let Some((c, t, tc, t0)) = &self.cp0_polyt {
            // Upstream `IdealHelmholtzCP0PolyT::all` (tau0 = Tc/T0), with its
            // exact special cases at t = 0 and t = -1.
            let (tc, t0) = (*tc, *t0);
            let tau0 = tc / t0;
            let eps10 = 10.0 * f64::EPSILON;
            let mut sum = 0.0;
            for i in 0..c.len() {
                if t[i].abs() < eps10 {
                    sum += c[i] - c[i] * tau / tau0 + c[i] * (tau / tau0).ln();
                } else if (t[i] + 1.0).abs() < eps10 {
                    sum += c[i] * tau / tc * (tau0 / tau).ln() + c[i] / tc * (tau - tau0);
                } else {
                    sum += -c[i] * tc.powf(t[i]) * tau.powf(-t[i]) / (t[i] * (t[i] + 1.0))
                        - c[i] * t0.powf(t[i] + 1.0) * tau / (tc * (t[i] + 1.0))
                        + c[i] * t0.powf(t[i]) / t[i];
                }
            }
            derivs.d00 += sum;
            let mut sum = 0.0;
            for i in 0..c.len() {
                if t[i].abs() < eps10 {
                    sum += c[i] / tau - c[i] / tau0;
                } else if (t[i] + 1.0).abs() < eps10 {
                    sum += c[i] / tc * (tau0 / tau).ln();
                } else {
                    sum += c[i] * tc.powf(t[i]) * tau.powf(-t[i] - 1.0) / (t[i] + 1.0)
                        - c[i] * tc.powf(t[i]) / (tau0.powf(t[i] + 1.0) * (t[i] + 1.0));
                }
            }
            derivs.d01 += sum;
            let mut sum = 0.0;
            for i in 0..c.len() {
                if t[i].abs() < eps10 {
                    sum += -c[i] / (tau * tau);
                } else if (t[i] + 1.0).abs() < eps10 {
                    sum += -c[i] / (tau * tc);
                } else {
                    sum += -c[i] * (tc / tau).powf(t[i]) / (tau * tau);
                }
            }
            derivs.d02 += sum;
            let mut sum = 0.0;
            for i in 0..c.len() {
                if t[i].abs() < eps10 {
                    sum += 2.0 * c[i] / (tau * tau * tau);
                } else if (t[i] + 1.0).abs() < eps10 {
                    sum += c[i] / (tau * tau * tc);
                } else {
                    sum += c[i] * (tc / tau).powf(t[i]) * (t[i] + 2.0) / (tau * tau * tau);
                }
            }
            derivs.d03 += sum;
            let mut sum = 0.0;
            for i in 0..c.len() {
                if t[i].abs() < eps10 {
                    sum += -6.0 * c[i] / pow4(tau);
                } else if (t[i] + 1.0).abs() < eps10 {
                    sum += -3.0 * c[i] / (pow3(tau) * tc);
                } else {
                    sum += -c[i] * (t[i] + 2.0) * (t[i] + 3.0) * (tc / tau).powf(t[i])
                        / (tau * tau * tau * tau);
                }
            }
            derivs.d04 += sum;
        }
    }
}

// ---------------------------------------------------------------------------
// Assembled EOS for one pure fluid
// ---------------------------------------------------------------------------

/// Residual + ideal Helmholtz machinery for one pure fluid, built once from
/// its `FluidData` (upstream: FluidLibrary parsing into the term containers).
pub struct HelmholtzEos {
    genexp: GenExp,
    non_analytic: Vec<NonAnalyticElement>,
    gaob: Vec<GaoBElement>,
    ideal: IdealContainer,
    /// `EOS[0].STATES.reducing.T` [K]
    pub t_reducing: f64,
    /// `EOS[0].STATES.reducing.rhomolar` [mol/m^3]
    pub rhomolar_reducing: f64,
    /// `EOS[0].gas_constant` [J/mol/K]
    pub gas_constant: f64,
    /// `EOS[0].molar_mass` [kg/mol]
    pub molar_mass: f64,
}

impl HelmholtzEos {
    pub fn new(fluid: &rustprop_core::fluid::FluidData) -> Self {
        let mut genexp = GenExp::new();
        let mut non_analytic = Vec::new();
        let mut gaob = Vec::new();
        for term in fluid.eos.alphar {
            match term {
                AlpharTerm::Power { n, d, t, l } => genexp.add_power(n, d, t, l),
                AlpharTerm::Gaussian {
                    n,
                    d,
                    t,
                    eta,
                    beta,
                    gamma,
                    epsilon,
                } => {
                    genexp.add_gaussian(n, d, t, eta, epsilon, beta, gamma);
                }
                AlpharTerm::NonAnalytic {
                    n,
                    a,
                    b,
                    beta,
                    big_a,
                    big_b,
                    big_c,
                    big_d,
                } => {
                    for i in 0..n.len() {
                        non_analytic.push(NonAnalyticElement {
                            n: n[i],
                            a: a[i],
                            b: b[i],
                            beta: beta[i],
                            big_a: big_a[i],
                            big_b: big_b[i],
                            big_c: big_c[i],
                            big_d: big_d[i],
                        });
                    }
                }
                AlpharTerm::Exponential { n, d, t, g, l } => {
                    genexp.add_exponential(n, d, t, g, l);
                }
                AlpharTerm::DoubleExponential {
                    n,
                    d,
                    t,
                    gd,
                    ld,
                    gt,
                    lt,
                } => {
                    genexp.add_double_exponential(n, d, t, gd, ld, gt, lt);
                }
                AlpharTerm::Lemmon2005 { n, d, t, l, m } => {
                    genexp.add_lemmon2005(n, d, t, l, m);
                }
                AlpharTerm::GaoB {
                    n,
                    t,
                    d,
                    eta,
                    beta,
                    gamma,
                    epsilon,
                    b,
                } => {
                    for i in 0..n.len() {
                        gaob.push(GaoBElement {
                            n: n[i],
                            t: t[i],
                            d: d[i],
                            eta: eta[i],
                            beta: beta[i],
                            gamma: gamma[i],
                            epsilon: epsilon[i],
                            b: b[i],
                        });
                    }
                }
            }
        }
        genexp.finish();

        let mut ideal = IdealContainer::default();
        for term in fluid.eos.alpha0 {
            match term {
                Alpha0Term::Lead { a1, a2 } => ideal.lead = Some((*a1, *a2)),
                Alpha0Term::LogTau { a } => ideal.log_tau = Some(*a),
                Alpha0Term::Power { n, t } => {
                    // Upstream throws on a second Power term.
                    assert!(ideal.power.is_none(), "duplicate IdealGasHelmholtzPower");
                    ideal.power = Some((n.to_vec(), t.to_vec()));
                }
                Alpha0Term::EnthalpyEntropyOffset { a1, a2, .. } => {
                    // Upstream `set()`: first call stores, later calls increment.
                    match &mut ideal.offset_core {
                        None => ideal.offset_core = Some((*a1, *a2)),
                        Some((c1, c2)) => {
                            *c1 += *a1;
                            *c2 += *a2;
                        }
                    }
                }
                Alpha0Term::PlanckEinstein { n, t } => {
                    // Upstream FluidLibrary.h: flip the sign of theta; c = 1, d = -1.
                    let pe = ideal.planck_einstein.get_or_insert_with(Default::default);
                    pe.n.extend_from_slice(n);
                    pe.theta.extend(t.iter().map(|v| -v));
                    pe.c.extend(std::iter::repeat_n(1.0, n.len()));
                    pe.d.extend(std::iter::repeat_n(-1.0, n.len()));
                }
                Alpha0Term::PlanckEinsteinFunctionT { n, v, tcrit } => {
                    // Upstream FluidLibrary.h: theta = -v/Tcrit; c = 1, d = -1.
                    let pe = ideal.planck_einstein.get_or_insert_with(Default::default);
                    pe.n.extend_from_slice(n);
                    pe.theta.extend(v.iter().map(|x| -x / tcrit));
                    pe.c.extend(std::iter::repeat_n(1.0, n.len()));
                    pe.d.extend(std::iter::repeat_n(-1.0, n.len()));
                }
                Alpha0Term::PlanckEinsteinGeneralized { n, t, c, d } => {
                    // Direct generalized form: t is theta already.
                    let pe = ideal.planck_einstein.get_or_insert_with(Default::default);
                    pe.n.extend_from_slice(n);
                    pe.theta.extend_from_slice(t);
                    pe.c.extend_from_slice(c);
                    pe.d.extend_from_slice(d);
                }
                Alpha0Term::Cp0Constant { cp_over_r, tc, t0 } => {
                    // Upstream throws on a second CP0Constant term.
                    assert!(
                        ideal.cp0_constant.is_none(),
                        "duplicate IdealGasHelmholtzCP0Constant"
                    );
                    ideal.cp0_constant = Some((*cp_over_r, *tc, *t0));
                }
                Alpha0Term::Cp0PolyT { c, t, tc, t0 } => {
                    // Upstream throws on a second CP0PolyT term (AlyLee may
                    // extend an existing one below, but not vice versa).
                    assert!(
                        ideal.cp0_polyt.is_none(),
                        "duplicate IdealGasHelmholtzCP0PolyT"
                    );
                    ideal.cp0_polyt = Some((c.to_vec(), t.to_vec(), *tc, *t0));
                }
                Alpha0Term::Cp0AlyLee {
                    c: constants,
                    tc,
                    t0,
                } => {
                    // Upstream FluidLibrary.h converts Aly-Lee at parse time:
                    // nonzero constant -> CP0PolyT (extend if present),
                    // sinh -> PE (n=B, theta=-2C/Tc, c=1, d=-1),
                    // cosh -> PE (n=-D, theta=-2E/Tc, c=1, d=+1).
                    if constants[0].abs() > 1e-14 {
                        match &mut ideal.cp0_polyt {
                            Some((pc, pt, _, _)) => {
                                pc.push(constants[0]);
                                pt.push(0.0);
                            }
                            None => {
                                ideal.cp0_polyt = Some((vec![constants[0]], vec![0.0], *tc, *t0));
                            }
                        }
                    }
                    let mut n = Vec::new();
                    let mut t = Vec::new();
                    let mut cc = Vec::new();
                    let mut dd = Vec::new();
                    if constants[1].abs() > 1e-14 {
                        n.push(constants[1]);
                        t.push(-2.0 * constants[2] / tc);
                        cc.push(1.0);
                        dd.push(-1.0);
                    }
                    if constants[3].abs() > 1e-14 {
                        n.push(-constants[3]);
                        t.push(-2.0 * constants[4] / tc);
                        cc.push(1.0);
                        dd.push(1.0);
                    }
                    let pe = ideal.planck_einstein.get_or_insert_with(Default::default);
                    pe.n.extend_from_slice(&n);
                    pe.theta.extend_from_slice(&t);
                    pe.c.extend_from_slice(&cc);
                    pe.d.extend_from_slice(&dd);
                }
            }
        }

        HelmholtzEos {
            genexp,
            non_analytic,
            gaob,
            ideal,
            t_reducing: fluid.eos.reducing.t,
            rhomolar_reducing: fluid.eos.reducing.rhomolar,
            gas_constant: fluid.eos.gas_constant,
            molar_mass: fluid.eos.molar_mass,
        }
    }

    /// All residual derivatives (upstream `ResidualHelmholtzContainer::all`):
    /// GenExp first into a fresh accumulator, then NonAnalytic, then GaoB.
    pub fn alphar_all(&self, tau: f64, delta: f64) -> HelmholtzDerivs {
        let mut derivs = HelmholtzDerivs::default();
        self.genexp.all(tau, delta, &mut derivs);
        non_analytic_all(&self.non_analytic, tau, delta, &mut derivs);
        gaob_all(&self.gaob, tau, delta, &mut derivs);
        derivs
    }

    /// All ideal-gas derivatives (upstream `IdealHelmholtzContainer::all`).
    pub fn alpha0_all(&self, tau: f64, delta: f64) -> HelmholtzDerivs {
        let mut derivs = HelmholtzDerivs::default();
        self.ideal.all(tau, delta, &mut derivs);
        derivs
    }
}
