//! Viscosity models (PLAN.md 6.1, structured slice) — port of the
//! `TransportRoutines` viscosity families the structured fluids use, and the
//! assembly in `HelmholtzEOSMixtureBackend::calc_viscosity`:
//! `eta = dilute + initial_density + residual (+ critical, which no ported
//! fluid has)`.
//!
//! - dilute: kinetic theory (Neufeld Omega22), collision integral,
//!   powers of T / Tr, collision-integral powers of T*;
//! - initial density: Rainwater-Friend (`eta_dilute * B_eta * rhomolar`)
//!   and the empirical (Tariq) form;
//! - higher order: modified Batschinski-Hildebrand and friction theory
//!   (the latter needs the state's p and dp/dT|rho).
//!
//! Fluids whose blocks carry `hardcoded` section tags (or ECS/Chung/rhosr
//! top-level models, which datagen leaves as `transport: None`) error with
//! `NotImplemented` until their slices land.

use crate::alpha::HelmholtzEos;
use rustprop_core::fluid::{
    Conductivity, ConductivityCritical, ConductivityDilute, ConductivityResidual, Viscosity,
    ViscosityDilute, ViscosityHigherOrder, ViscosityInitialDensity,
};
use rustprop_core::{Error, Result};

/// Upstream `calc_viscosity` for a pure fluid at a fully-determined state:
/// the state's (T, rhomolar, p) — two-phase states evaluate at the mixture
/// density exactly as upstream does.
pub fn viscosity(eos: &HelmholtzEos, v: &Viscosity, t: f64, rhomolar: f64, p: f64) -> Result<f64> {
    let dilute = viscosity_dilute(eos, v, t)?;
    let initial_density = match &v.initial_density {
        None => 0.0,
        Some(ViscosityInitialDensity::RainwaterFriend { b, t: bt }) => {
            // B_eta* summed over powers of Tstar; B_eta = N_A * sigma^3 * B_eta*.
            let tstar = t / v.epsilon_over_k;
            let sigma = v.sigma_eta;
            let mut summer = 0.0;
            for i in 0..b.len() {
                summer += b[i] * tstar.powf(bt[i]);
            }
            let b_eta = 6.02214129e23 * sigma.powf(3.0) * summer; // [m^3/mol]
            dilute * b_eta * rhomolar
        }
        Some(ViscosityInitialDensity::Empirical {
            n,
            d,
            t: te,
            t_reducing,
            rhomolar_reducing,
        }) => {
            let tau = t_reducing / t;
            let delta = rhomolar / rhomolar_reducing;
            let mut summer = 0.0;
            for i in 0..n.len() {
                summer += n[i] * delta.powf(d[i]) * tau.powf(te[i]);
            }
            summer
        }
    };
    let residual = viscosity_higher_order(eos, v, t, rhomolar, p)?;
    Ok(dilute + initial_density + residual)
}

fn viscosity_dilute(eos: &HelmholtzEos, v: &Viscosity, t: f64) -> Result<f64> {
    Ok(match &v.dilute {
        ViscosityDilute::KineticTheory => {
            let tstar = t / v.epsilon_over_k;
            let sigma_nm = v.sigma_eta * 1e9;
            let molar_mass_kgkmol = eos.molar_mass * 1000.0;
            // Neufeld's empirical Omega22 collision integral.
            let omega22 = 1.16145 * tstar.powf(-0.14874)
                + 0.52487 * (-0.77320 * tstar).exp()
                + 2.16178 * (-2.43787 * tstar).exp();
            26.692e-9 * (molar_mass_kgkmol * t).sqrt() / (sigma_nm.powf(2.0) * omega22)
        }
        ViscosityDilute::CollisionIntegral {
            a,
            t: at,
            c,
            molar_mass,
        } => {
            let tstar = t / v.epsilon_over_k;
            let sigma_nm = v.sigma_eta * 1e9;
            let molar_mass_kgkmol = molar_mass * 1000.0;
            let ln_tstar = tstar.ln();
            let mut summer = 0.0;
            for i in 0..a.len() {
                summer += a[i] * ln_tstar.powf(at[i]);
            }
            let s = summer.exp();
            c * (molar_mass_kgkmol * t).sqrt() / (sigma_nm.powf(2.0) * s)
        }
        ViscosityDilute::PowersOfT { a, t: at } => {
            let mut summer = 0.0;
            for i in 0..a.len() {
                summer += a[i] * t.powf(at[i]);
            }
            summer
        }
        ViscosityDilute::PowersOfTr {
            a,
            t: at,
            t_reducing,
        } => {
            let tr = t / t_reducing;
            let mut summer = 0.0;
            for i in 0..a.len() {
                summer += a[i] * tr.powf(at[i]);
            }
            summer
        }
        ViscosityDilute::CollisionIntegralPowersOfTstar {
            a,
            t: at,
            c,
            t_reducing,
        } => {
            let tstar = t / t_reducing;
            let mut summer = 0.0;
            for i in 0..a.len() {
                summer += a[i] * tstar.powf(at[i]);
            }
            c * t.sqrt() / summer
        }
        ViscosityDilute::Hardcoded { name } => {
            return Err(Error::NotImplemented(format!(
                "hardcoded dilute viscosity [{name}] is not ported yet"
            )));
        }
    })
}

fn viscosity_higher_order(
    eos: &HelmholtzEos,
    v: &Viscosity,
    t: f64,
    rhomolar: f64,
    p: f64,
) -> Result<f64> {
    Ok(match &v.higher_order {
        ViscosityHigherOrder::ModifiedBatschinskiHildebrand {
            a,
            d1,
            t1,
            gamma,
            l,
            f,
            d2,
            t2,
            g,
            h,
            p: pp,
            q,
            t_reduce,
            rhomolar_reduce,
        } => {
            let delta = rhomolar / rhomolar_reduce;
            let tau = t_reduce / t;
            let mut s = 0.0;
            for i in 0..a.len() {
                s += a[i]
                    * delta.powf(d1[i])
                    * tau.powf(t1[i])
                    * (gamma[i] * delta.powf(l[i])).exp();
            }
            let mut big_f = 0.0;
            for i in 0..f.len() {
                big_f += f[i] * delta.powf(d2[i]) * tau.powf(t2[i]);
            }
            let mut summer_numer = 0.0;
            for i in 0..g.len() {
                summer_numer += g[i] * tau.powf(h[i]);
            }
            let mut summer_denom = 0.0;
            for i in 0..pp.len() {
                summer_denom += pp[i] * tau.powf(q[i]);
            }
            let delta0 = summer_numer / summer_denom;
            s + big_f * (1.0 / (delta0 - delta) - 1.0 / delta0)
        }
        ViscosityHigherOrder::FrictionTheory {
            ai,
            aa,
            ar,
            aaa,
            arr,
            adrdr,
            aii,
            arrr,
            aaaa,
            na,
            naa,
            nr,
            nrr,
            nii,
            nrrr,
            naaa,
            c1,
            c2,
            t_reduce,
        } => {
            let tau = t_reduce / t;
            let psi1 = tau.exp() - c1;
            let psi2 = (tau.powf(2.0)).exp() - c2;

            let ki = (ai[0] + ai[1] * psi1 + ai[2] * psi2) * tau;
            let ka = (aa[0] + aa[1] * psi1 + aa[2] * psi2) * tau.powf(*na);
            let kr = (ar[0] + ar[1] * psi1 + ar[2] * psi2) * tau.powf(*nr);
            let kaa = (aaa[0] + aaa[1] * psi1 + aaa[2] * psi2) * tau.powf(*naa);
            let (krr, kdrdr) = if arr.is_empty() {
                (
                    0.0,
                    (adrdr[0] + adrdr[1] * psi1 + adrdr[2] * psi2) * tau.powf(*nrr),
                )
            } else {
                (
                    (arr[0] + arr[1] * psi1 + arr[2] * psi2) * tau.powf(*nrr),
                    0.0,
                )
            };
            let kii = if aii.is_empty() {
                0.0
            } else {
                (aii[0] + aii[1] * psi1 + aii[2] * psi2) * tau.powf(*nii)
            };
            let (krrr, kaaa) = if !arrr.is_empty() && !aaaa.is_empty() {
                (
                    (arrr[0] + arrr[1] * psi1 + arrr[2] * psi2) * tau.powf(*nrrr),
                    (aaaa[0] + aaaa[1] * psi1 + aaaa[2] * psi2) * tau.powf(*naaa),
                )
            } else {
                (0.0, 0.0)
            };

            let p_bar = p / 1e5;
            // dp/dT|rho = rho*R*(1 + delta*d10 - delta*tau*d11)
            let tau_eos = eos.t_reducing / t;
            let delta = rhomolar / eos.rhomolar_reducing;
            let d = eos.alphar_all(tau_eos, delta);
            let dpdt =
                rhomolar * eos.gas_constant * (1.0 + delta * d.d10 - delta * tau_eos * d.d11);
            let pr = t * dpdt / 1e5;
            let pa = p_bar - pr;
            let pid = rhomolar * eos.gas_constant * t / 1e5;
            let deltapr = pr - pid;

            ka * pa
                + kr * deltapr
                + ki * pid
                + kaa * pa * pa
                + kdrdr * deltapr * deltapr
                + krr * pr * pr
                + kii * pid * pid
                + krrr * pr * pr * pr
                + kaaa * pa * pa * pa
        }
        ViscosityHigherOrder::Hardcoded { name } => {
            return Err(Error::NotImplemented(format!(
                "hardcoded higher-order viscosity [{name}] is not ported yet"
            )));
        }
    })
}

/// Upstream `calc_conductivity` (`lambda = dilute + residual + critical`)
/// for a pure fluid at a fully-determined state. The Olchowy-Sengers
/// critical enhancement needs the fluid's own viscosity and cp/cv — cp is
/// undefined in the two-phase region, where upstream throws.
#[allow(clippy::too_many_arguments)]
pub fn conductivity(
    eos: &HelmholtzEos,
    c: &Conductivity,
    viscosity_model: Option<&Viscosity>,
    reducing_p: f64,
    t: f64,
    rhomolar: f64,
    p: f64,
    two_phase: bool,
) -> Result<f64> {
    let dilute = match &c.dilute {
        ConductivityDilute::RatioOfPolynomials {
            a,
            n,
            b,
            m,
            t_reducing,
        } => {
            let tr = t / t_reducing;
            let mut summer1 = 0.0;
            for i in 0..a.len() {
                summer1 += a[i] * tr.powf(n[i]);
            }
            let mut summer2 = 0.0;
            for i in 0..b.len() {
                summer2 += b[i] * tr.powf(m[i]);
            }
            summer1 / summer2
        }
        ConductivityDilute::Eta0AndPoly { a, t: at } => {
            let v = viscosity_model.ok_or_else(|| {
                Error::NotImplemented(
                    "conductivity eta0_and_poly needs the fluid's (unported) viscosity model"
                        .into(),
                )
            })?;
            let eta0_upas = viscosity_dilute(eos, v, t)? * 1e6;
            let tau = eos.t_reducing / t;
            let mut summer = a[0] * eta0_upas;
            for i in 1..a.len() {
                summer += a[i] * tau.powf(at[i]);
            }
            summer
        }
        ConductivityDilute::Hardcoded { name } => {
            return Err(Error::NotImplemented(format!(
                "hardcoded dilute conductivity [{name}] is not ported yet"
            )));
        }
    };

    let residual = match &c.residual {
        ConductivityResidual::Polynomial {
            b,
            t: bt,
            d,
            t_reducing,
            rhomass_reducing,
        } => {
            let tau = t_reducing / t;
            let delta = rhomolar * eos.molar_mass / rhomass_reducing;
            let mut summer = 0.0;
            for i in 0..b.len() {
                summer += b[i] * tau.powf(bt[i]) * delta.powf(d[i]);
            }
            summer
        }
        ConductivityResidual::PolynomialAndExponential {
            a,
            t: at,
            d,
            gamma,
            l,
        } => {
            let tau = eos.t_reducing / t;
            let delta = rhomolar / eos.rhomolar_reducing;
            let mut summer = 0.0;
            for i in 0..a.len() {
                summer += a[i]
                    * tau.powf(at[i])
                    * delta.powf(d[i])
                    * (-gamma[i] * delta.powf(l[i])).exp();
            }
            summer
        }
    };

    let critical = match &c.critical {
        None => 0.0,
        Some(ConductivityCritical::SimplifiedOlchowySengers {
            k,
            r0,
            gamma,
            nu,
            big_gamma,
            zeta0,
            qd,
            t_ref,
        }) => {
            // cp/cv (and the enhancement) are undefined in the two-phase
            // region — upstream's cpmolar() throws there.
            if two_phase {
                return Err(Error::Value(
                    "Input is two-phase and the critical conductivity enhancement is not defined"
                        .into(),
                ));
            }
            let v = viscosity_model.ok_or_else(|| {
                Error::NotImplemented(
                    "the Olchowy-Sengers enhancement needs the fluid's (unported) viscosity model"
                        .into(),
                )
            })?;
            let tc = eos.t_reducing;
            let rhoc = eos.rhomolar_reducing;
            let pcrit = reducing_p;
            let tref = if t_ref.is_finite() { *t_ref } else { 1.5 * tc };

            let delta = rhomolar / rhoc;
            let tau = tc / t;
            let dd = eos.alphar_all(tau, delta);
            let dp_drho =
                eos.gas_constant * t * (1.0 + 2.0 * delta * dd.d10 + delta * delta * dd.d20);
            let x = pcrit / rhoc.powf(2.0) * rhomolar / dp_drho;

            let tau_ref = tc / tref;
            let dref = eos.alphar_all(tau_ref, delta);
            let dp_drho_ref =
                eos.gas_constant * tref * (1.0 + 2.0 * delta * dref.d10 + delta * delta * dref.d20);
            let xref = pcrit / rhoc.powf(2.0) * rhomolar / dp_drho_ref * tref / t;
            let num = x - xref;

            // No critical enhancement if the numerator is negative, zero, or
            // just a tiny bit positive due to roundoff (Lemmon, IJT, 2004).
            if num < f64::EPSILON * 10.0 {
                0.0
            } else {
                let zeta = zeta0 * (num / big_gamma).powf(nu / gamma);
                let cp = eos.cpmolar(t, rhomolar);
                let cv = eos.cvmolar(t, rhomolar);
                let mu = viscosity(eos, v, t, rhomolar, p)?;
                let pi = std::f64::consts::PI;
                let omega_tilde =
                    2.0 / pi * ((cp - cv) / cp * (zeta * qd).atan() + cv / cp * (zeta * qd));
                let omega_tilde0 = 2.0 / pi
                    * (1.0
                        - (-1.0
                            / (1.0 / (qd * zeta)
                                + 1.0 / 3.0 * (zeta * qd) * (zeta * qd) / delta / delta))
                            .exp());
                rhomolar * cp * r0 * k * t / (6.0 * pi * mu * zeta) * (omega_tilde - omega_tilde0)
            }
        }
        Some(ConductivityCritical::Hardcoded { name }) => {
            return Err(Error::NotImplemented(format!(
                "hardcoded critical conductivity [{name}] is not ported yet"
            )));
        }
    };

    Ok(dilute + residual + critical)
}
