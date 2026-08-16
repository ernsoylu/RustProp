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

// Hardcoded-model constants stay verbatim (excessive-precision literals
// included) and index loops mirror the upstream loops.
#![allow(clippy::excessive_precision, clippy::needless_range_loop)]

use crate::alpha::HelmholtzEos;
use rustprop_core::fluid::{
    Conductivity, ConductivityCritical, ConductivityDilute, ConductivityModel,
    ConductivityResidual, FluidData, Viscosity, ViscosityDilute, ViscosityHigherOrder,
    ViscosityInitialDensity, ViscosityModel,
};
use rustprop_core::{Error, Result};

fn pow2(x: f64) -> f64 {
    x * x
}
fn pow3(x: f64) -> f64 {
    x * x * x
}
fn pow4(x: f64) -> f64 {
    pow2(x) * pow2(x)
}
fn pow5(x: f64) -> f64 {
    pow4(x) * x
}

/// Upstream `powInt` (sequential multiplication), as in the alpha machinery.
fn pow_int(x: f64, y: i32) -> f64 {
    if y == 0 {
        return 1.0;
    }
    let (x_in, y_in) = if y < 0 { (1.0 / x, -y) } else { (x, y) };
    let mut product = x_in;
    for _ in 1..y_in {
        product *= x_in;
    }
    product
}

/// A resolved ECS reference fluid: its EOS, document, and transport models
/// (the caller — who owns the registry — resolves the name).
pub struct EcsRef<'a> {
    pub eos: &'a HelmholtzEos,
    pub fluid: &'a FluidData,
    pub viscosity: Option<&'a ViscosityModel>,
    pub conductivity: Option<&'a ConductivityModel>,
}

/// Resolver from a reference-fluid name to its pieces.
pub type EcsResolver<'a> = dyn Fn(&str) -> Result<EcsRef<'a>> + 'a;

/// Upstream `T_critical()`/`rhomolar_critical()`: superancillary numerical
/// values when present.
fn crit_of(fluid: &FluidData) -> (f64, f64) {
    match &fluid.eos.superancillary {
        Some(sa) => (sa.t_crit_num, sa.rho_crit_num),
        None => (fluid.states.critical.t, fluid.states.critical.rhomolar),
    }
}

/// Upstream `calc_viscosity` for a pure fluid at a fully-determined state:
/// the state's (T, rhomolar, p) — two-phase states evaluate at the mixture
/// density exactly as upstream does.
pub fn viscosity(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    model: &ViscosityModel,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
) -> Result<f64> {
    match model {
        ViscosityModel::Structured(v) => viscosity_structured(eos, v, t, rhomolar, p),
        ViscosityModel::Hardcoded { name } => match *name {
            "Water" => Ok(viscosity_water_hardcoded(eos, t, rhomolar)),
            "HeavyWater" => Ok(viscosity_heavywater_hardcoded(eos, t, rhomolar)),
            "Helium" => Ok(viscosity_helium_hardcoded(eos, t, rhomolar)),
            "R23" => Ok(viscosity_r23_hardcoded(t, rhomolar)),
            "Methanol" => Ok(viscosity_methanol_hardcoded(eos, t, rhomolar)),
            "m-Xylene" => Ok(viscosity_m_xylene_hardcoded(t, rhomolar)),
            "o-Xylene" => Ok(viscosity_o_xylene_hardcoded(t, rhomolar)),
            "p-Xylene" => Ok(viscosity_p_xylene_hardcoded(t, rhomolar)),
            other => Err(Error::NotImplemented(format!(
                "hardcoded viscosity [{other}] is not ported yet"
            ))),
        },
        ViscosityModel::Chung {
            rhomolar_critical,
            acentric,
            molar_mass,
            t_critical,
            dipole_moment_d,
            kappa: _,
        } => Ok(viscosity_chung(
            t,
            rhomolar,
            *rhomolar_critical,
            *acentric,
            *molar_mass,
            *t_critical,
            *dipole_moment_d,
        )),
        ViscosityModel::RhosrCs {
            c,
            c_liq,
            c_vap,
            rhosr_critical,
            x_crossover: _,
        } => Ok(viscosity_rhosr(
            eos,
            t,
            rhomolar,
            *c,
            c_liq,
            c_vap,
            *rhosr_critical,
        )),
        ViscosityModel::Ecs {
            reference_fluid,
            psi_a,
            psi_t,
            psi_rhomolar_reducing,
            sigma_eta,
            epsilon_over_k,
        } => {
            let resolver = ecs.ok_or_else(|| {
                Error::NotImplemented("ECS evaluation needs a reference-fluid resolver".into())
            })?;
            let reference = resolver(reference_fluid)?;
            viscosity_ecs(
                eos,
                fluid,
                &reference,
                t,
                rhomolar,
                psi_a,
                psi_t,
                *psi_rhomolar_reducing,
                *sigma_eta,
                *epsilon_over_k,
            )
        }
    }
}

/// Upstream `conformal_state_solver`: 2-D Newton matching the reference
/// fluid's (alphar, Z) to the fluid of interest's, with geometric step
/// halving, ftol 1e-9, 50-iteration cap.
fn conformal_state_solver(
    ref_eos: &HelmholtzEos,
    ref_tc: f64,
    ref_rhoc: f64,
    alphar_target: f64,
    z_target: f64,
    t0: &mut f64,
    rhomolar0: &mut f64,
) -> Result<()> {
    let eval = |t: f64, rho: f64| -> (f64, f64, crate::alpha::HelmholtzDerivs, f64) {
        let tau = ref_eos.t_reducing / t;
        let delta = rho / ref_eos.rhomolar_reducing;
        let d = ref_eos.alphar_all(tau, delta);
        (d.d00, 1.0 + delta * d.d10, d, delta)
    };
    let mut iter = 0;
    let mut resid = 9e30;
    let mut resid_old: f64;
    let (mut a0, mut z0, mut d, mut delta) = eval(*t0, *rhomolar0);
    loop {
        let dtau_dt = -ref_tc / (*t0 * *t0);
        let ddelta_drho = 1.0 / ref_rhoc;
        let r0 = a0 - alphar_target;
        let r1 = z0 - z_target;
        let j00 = d.d01 * dtau_dt;
        let j01 = d.d10 * ddelta_drho;
        let j10 = delta * d.d11 * dtau_dt;
        let j11 = (delta * d.d20 + d.d10) * ddelta_drho;
        // Direct 2x2 solve of J v = -r (upstream uses Eigen's QR; identical
        // solution to roundoff).
        let det = j00 * j11 - j01 * j10;
        if !det.is_finite() || det.abs() < 1e-300 {
            return Err(Error::Solution(
                "conformal state solver: singular Jacobian".into(),
            ));
        }
        let v0 = -(j11 * r0 - j01 * r1) / det;
        let v1 = -(-j10 * r0 + j00 * r1) / det;
        let mut good_solution = false;
        let (t0_init, rho0_init) = (*t0, *rhomolar0);
        resid_old = (r0 * r0 + r1 * r1).sqrt();
        let mut frac = 1.0;
        while frac > 0.001 {
            let t_new = t0_init + frac * v0;
            let rho_new = rho0_init + frac * v1;
            if t_new > 0.0 && rho_new > 0.0 {
                let (a_n, z_n, d_n, delta_n) = eval(t_new, rho_new);
                resid = ((a_n - alphar_target).powi(2) + (z_n - z_target).powi(2)).sqrt();
                if resid.is_finite() && resid <= resid_old {
                    good_solution = true;
                    *t0 = t_new;
                    *rhomolar0 = rho_new;
                    a0 = a_n;
                    z0 = z_n;
                    d = d_n;
                    delta = delta_n;
                    break;
                }
            }
            frac /= 2.0;
        }
        if !good_solution {
            return Err(Error::Value("Not able to get a solution".into()));
        }
        iter += 1;
        if iter > 50 {
            return Err(Error::Value(format!(
                "conformal_state_solver took too many iterations; residual is {resid}"
            )));
        }
        if resid.abs() <= 1e-9 {
            return Ok(());
        }
    }
}

/// The fluid's Lennard-Jones (sigma [m], epsilon/k [K]) pair as upstream
/// stores it on `transport`: from the viscosity block when given, else the
/// `default_transport` critical-point estimation.
fn fluid_lennard_jones(eos: &HelmholtzEos, viscosity_model: Option<&ViscosityModel>) -> (f64, f64) {
    match viscosity_model {
        Some(ViscosityModel::Ecs {
            sigma_eta,
            epsilon_over_k,
            ..
        }) if sigma_eta.is_finite() && epsilon_over_k.is_finite() => (*sigma_eta, *epsilon_over_k),
        Some(ViscosityModel::Structured(v))
            if v.sigma_eta.is_finite() && v.epsilon_over_k.is_finite() =>
        {
            (v.sigma_eta, v.epsilon_over_k)
        }
        _ => {
            let rho_crit_moll = eos.rhomolar_reducing / 1000.0;
            (
                0.809 / rho_crit_moll.powf(1.0 / 3.0) / 1e9,
                eos.t_reducing / 1.2593,
            )
        }
    }
}

/// Upstream `viscosity_dilute_kinetic_theory` with explicit L-J parameters
/// (Neufeld's Omega22), in Pa-s.
fn kinetic_theory_dilute(eos: &HelmholtzEos, t: f64, sigma: f64, epsilon_over_k: f64) -> f64 {
    let tstar = t / epsilon_over_k;
    let sigma_nm = sigma * 1e9;
    let molar_mass_kgkmol = eos.molar_mass * 1000.0;
    let omega22 = 1.16145 * tstar.powf(-0.14874)
        + 0.52487 * (-0.77320 * tstar).exp()
        + 2.16178 * (-2.43787 * tstar).exp();
    26.692e-9 * (molar_mass_kgkmol * t).sqrt() / (sigma_nm.powf(2.0) * omega22)
}

/// Upstream `viscosity_ECS`.
#[allow(clippy::too_many_arguments)]
fn viscosity_ecs(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    reference: &EcsRef,
    t: f64,
    rhomolar: f64,
    psi_a: &[f64],
    psi_t: &[f64],
    psi_rhomolar_reducing: f64,
    sigma_eta: f64,
    epsilon_over_k: f64,
) -> Result<f64> {
    let m = eos.molar_mass;
    let m0 = reference.eos.molar_mass;
    let (tc, rhocmolar) = crit_of(fluid);
    let (tc0, rhocmolar0) = crit_of(reference.fluid);

    let mut psi = 0.0;
    for i in 0..psi_a.len() {
        psi += psi_a[i] * (rhomolar / psi_rhomolar_reducing).powf(psi_t[i]);
    }

    // Dilute part: kinetic theory with the block's L-J parameters (or the
    // default_transport estimation when absent).
    let (sigma, eps) = if sigma_eta.is_finite() && epsilon_over_k.is_finite() {
        (sigma_eta, epsilon_over_k)
    } else {
        let rho_crit_moll = eos.rhomolar_reducing / 1000.0;
        (
            0.809 / rho_crit_moll.powf(1.0 / 3.0) / 1e9,
            eos.t_reducing / 1.2593,
        )
    };
    let eta_dilute = kinetic_theory_dilute(eos, t, sigma, eps);

    // Conformal state
    let tau = eos.t_reducing / t;
    let delta = rhomolar / eos.rhomolar_reducing;
    let dd = eos.alphar_all(tau, delta);
    let alphar_target = dd.d00;
    let z_target = 1.0 + delta * dd.d10;
    let mut t0 = t / (tc / tc0);
    let mut rhomolar0 = rhomolar * (rhocmolar0 / rhocmolar);
    conformal_state_solver(
        reference.eos,
        tc0,
        rhocmolar0,
        alphar_target,
        z_target,
        &mut t0,
        &mut rhomolar0,
    )?;

    // Reference background at (rho0*psi, T0)
    let rv = reference.viscosity.ok_or_else(|| {
        Error::NotImplemented("ECS reference fluid's viscosity model is unavailable".into())
    })?;
    let ViscosityModel::Structured(rvs) = rv else {
        return Err(Error::NotImplemented(
            "ECS reference fluid's viscosity model is not structured".into(),
        ));
    };
    let rho_ref = rhomolar0 * psi;
    let p_ref = reference.eos.pressure(t0, rho_ref);
    let eta_dilute_ref = viscosity_dilute(reference.eos, rvs, t0)?;
    let eta_resid = viscosity_background(reference.eos, rvs, eta_dilute_ref, t0, rho_ref, p_ref)?;

    let f = t / t0;
    let h = rhomolar0 / rhomolar;
    let f_eta = f.sqrt() * h.powf(-2.0 / 3.0) * (m / m0).sqrt();
    Ok(eta_dilute + eta_resid * f_eta)
}

/// The reference fluid's `calc_viscosity_background`:
/// initial-density + higher-order at the conformal state.
fn viscosity_background(
    eos: &HelmholtzEos,
    v: &Viscosity,
    eta_dilute: f64,
    t: f64,
    rhomolar: f64,
    p: f64,
) -> Result<f64> {
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
            eta_dilute * b_eta * rhomolar
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
    Ok(initial_density + residual)
}

/// Upstream `viscosity_Chung` (evaluates with kappa = 0 regardless of the
/// document's kappa, exactly as upstream's local `kappa = 0`).
#[allow(clippy::many_single_char_names)]
fn viscosity_chung(
    t: f64,
    rhomolar: f64,
    rhomolar_critical: f64,
    acentric: f64,
    molar_mass: f64,
    t_critical: f64,
    dipole_moment_d: f64,
) -> f64 {
    let a0 = [
        0.0, 6.32402, 0.12102e-2, 5.28346, 6.62263, 19.74540, -1.89992, 24.27450, 0.79716,
        -0.23816, 0.68629e-1,
    ];
    let a1 = [
        0.0,
        50.41190,
        -0.11536e-2,
        254.20900,
        38.09570,
        7.63034,
        -12.53670,
        3.44945,
        1.11764,
        0.67695e-1,
        0.34793,
    ];
    let a2 = [
        0.0,
        -51.68010,
        -0.62571e-2,
        -168.48100,
        -8.46414,
        -14.35440,
        4.98529,
        -11.29130,
        0.12348e-1,
        -0.81630,
        0.59256,
    ];
    let a3 = [
        0.0, 1189.02000, 0.37283e-1, 3898.27000, 31.41780, 31.52670, -18.15070, 69.34660, -4.11661,
        4.02528, -0.72663,
    ];
    let vc_cm3mol = 1.0 / (rhomolar_critical / 1e6);
    let m_gmol = molar_mass * 1000.0;
    let tc = t_critical;
    let mu_d = dipole_moment_d;
    let kappa = 0.0;
    let mu_r = 131.3 * mu_d / (vc_cm3mol * tc).sqrt();
    let mut a = [0.0f64; 11];
    for i in 1..=10 {
        a[i] = a0[i] + a1[i] * acentric + a2[i] * mu_r.powf(4.0) + a3[i] * kappa;
    }
    let f_c = 1.0 - 0.2756 * acentric + 0.059035 * mu_r.powf(4.0) + kappa;
    let epsilon_over_k = tc / 1.2593;
    let rho_molcm3 = rhomolar / 1e6;
    let tstar = t / epsilon_over_k;
    let omega_2_2 = 1.16145 * tstar.powf(-0.14874)
        + 0.52487 * (-0.77320 * tstar).exp()
        + 2.16178 * (-2.43787 * tstar).exp()
        - 6.435e-4 * tstar.powf(0.14874) * (18.0323 * tstar.powf(-0.76830) - 7.27371).sin();
    let eta0_p = 4.0785e-5 * (m_gmol * t).sqrt() / (vc_cm3mol.powf(2.0 / 3.0) * omega_2_2) * f_c;
    let y = rho_molcm3 * vc_cm3mol / 6.0;
    let g_1 = (1.0 - 0.5 * y) / (1.0 - y).powf(3.0);
    let g_2 = (a[1] * (1.0 - (-a[4] * y).exp()) / y + a[2] * g_1 * (a[5] * y).exp() + a[3] * g_1)
        / (a[1] * a[4] + a[2] + a[3]);
    let eta_k_p = eta0_p * (1.0 / g_2 + a[6] * y);
    let eta_p_p = (36.344e-6 * (m_gmol * tc).sqrt() / vc_cm3mol.powf(2.0 / 3.0))
        * a[7]
        * y.powf(2.0)
        * g_2
        * (a[8] + a[9] / tstar + a[10] / tstar.powf(2.0)).exp();
    (eta_k_p + eta_p_p) / 10.0
}

/// Upstream `viscosity_rhosr` — residual-entropy-scaled corresponding
/// states; the dilute part is kinetic theory with `default_transport`'s
/// Chung-estimated L-J parameters from the reducing state.
fn viscosity_rhosr(
    eos: &HelmholtzEos,
    t: f64,
    rhomolar: f64,
    c: f64,
    c_liq: &[f64],
    c_vap: &[f64],
    rhosr_critical: f64,
) -> f64 {
    // default_transport: sigma/epsilon estimated from the reducing state.
    let rho_crit_moll = eos.rhomolar_reducing / 1000.0;
    let sigma_eta = 0.809 / rho_crit_moll.powf(1.0 / 3.0) / 1e9;
    let epsilon_over_k = eos.t_reducing / 1.2593;
    let tstar = t / epsilon_over_k;
    let sigma_nm = sigma_eta * 1e9;
    let molar_mass_kgkmol = eos.molar_mass * 1000.0;
    let omega22 = 1.16145 * tstar.powf(-0.14874)
        + 0.52487 * (-0.77320 * tstar).exp()
        + 2.16178 * (-2.43787 * tstar).exp();
    let eta_dilute = 26.692e-9 * (molar_mass_kgkmol * t).sqrt() / (sigma_nm.powf(2.0) * omega22);

    let tau = eos.t_reducing / t;
    let delta = rhomolar / eos.rhomolar_reducing;
    let d = eos.alphar_all(tau, delta);
    let x = rhomolar * eos.gas_constant * (tau * d.d01 - d.d00) / rhosr_critical;
    let psi_liq = 1.0 / (1.0 + (-100.0 * (x - 2.0)).exp());
    let f_liq = c_liq[0] + x * (c_liq[1] + x * (c_liq[2] + x * c_liq[3]));
    let f_vap = c_vap[0] + x * (c_vap[1] + x * (c_vap[2] + x * c_vap[3]));
    let etastar_ref = (psi_liq * f_liq + (1.0 - psi_liq) * f_vap).exp();
    let etastar_fluid = 1.0 + c * (etastar_ref - 1.0);
    etastar_fluid * eta_dilute
}

fn viscosity_structured(
    eos: &HelmholtzEos,
    v: &Viscosity,
    t: f64,
    rhomolar: f64,
    p: f64,
) -> Result<f64> {
    let dilute = viscosity_dilute(eos, v, t)?;
    let background = viscosity_background(eos, v, dilute, t, rhomolar, p)?;
    Ok(dilute + background)
}

fn viscosity_dilute(eos: &HelmholtzEos, v: &Viscosity, t: f64) -> Result<f64> {
    Ok(match &v.dilute {
        ViscosityDilute::Hardcoded { name } => match *name {
            "Ethane" => viscosity_dilute_ethane(t),
            "Cyclohexane" => viscosity_dilute_cyclohexane(t),
            "CarbonDioxideLaeseckeJPCRD2017" => viscosity_dilute_co2_laesecke(t),
            other => {
                return Err(Error::NotImplemented(format!(
                    "hardcoded dilute viscosity [{other}] is not ported yet"
                )));
            }
        },
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
        ViscosityHigherOrder::Hardcoded { name } => match *name {
            "Ethane" => viscosity_ethane_higher_order(t, rhomolar),
            "Benzene" => viscosity_benzene_higher_order(eos, t, rhomolar),
            "Hydrogen" => viscosity_hydrogen_higher_order(eos, t, rhomolar),
            "Toluene" => viscosity_toluene_higher_order(eos, t, rhomolar),
            "n-Hexane" => viscosity_hexane_higher_order(eos, t, rhomolar),
            "n-Heptane" => viscosity_heptane_higher_order(eos, t, rhomolar),
            "CarbonDioxideLaeseckeJPCRD2017" => {
                viscosity_co2_higher_order_laesecke(eos, t, rhomolar)
            }
            other => {
                return Err(Error::NotImplemented(format!(
                    "hardcoded higher-order viscosity [{other}] is not ported yet"
                )));
            }
        },
    })
}

/// Upstream `calc_conductivity` (`lambda = dilute + residual + critical`)
/// for a pure fluid at a fully-determined state. Two-phase states evaluate
/// at the mixture density with the plain single-phase formulas — upstream's
/// `calc_cpmolar`/`calc_cvmolar` carry no two-phase guard, so the
/// Olchowy-Sengers enhancement (and everything else) computes verbatim.
#[allow(clippy::too_many_arguments)]
pub fn conductivity(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    model: &ConductivityModel,
    viscosity_model: Option<&ViscosityModel>,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
) -> Result<f64> {
    match model {
        ConductivityModel::Structured(c) => {
            conductivity_structured(eos, fluid, c, viscosity_model, t, rhomolar, p, ecs)
        }
        ConductivityModel::Hardcoded { name } => match *name {
            "Water" => {
                conductivity_water_hardcoded(eos, fluid, viscosity_model, t, rhomolar, p, ecs)
            }
            "HeavyWater" => Ok(conductivity_heavywater_hardcoded(eos, t, rhomolar)),
            "Helium" => {
                conductivity_helium_hardcoded(eos, fluid, viscosity_model, t, rhomolar, p, ecs)
            }
            "R23" => Ok(conductivity_r23_hardcoded(t, rhomolar)),
            "Methane" => Ok(conductivity_methane_hardcoded(eos, fluid, t, rhomolar)),
            other => Err(Error::NotImplemented(format!(
                "hardcoded conductivity [{other}] is not ported yet"
            ))),
        },
        ConductivityModel::Ecs {
            reference_fluid,
            psi_a,
            psi_t,
            psi_rhomolar_reducing,
            f_int_a,
            f_int_t,
            f_int_t_reducing,
        } => {
            let resolver = ecs.ok_or_else(|| {
                Error::NotImplemented("ECS evaluation needs a reference-fluid resolver".into())
            })?;
            let reference = resolver(reference_fluid)?;
            conductivity_ecs(
                eos,
                fluid,
                &reference,
                viscosity_model,
                t,
                rhomolar,
                p,
                ecs,
                psi_a,
                psi_t,
                *psi_rhomolar_reducing,
                f_int_a,
                f_int_t,
                *f_int_t_reducing,
            )
        }
    }
}

/// Upstream `conductivity_ECS`: `lambda = lambda_int + lambda_dilute +
/// lambda_resid·F_lambda + lambda_crit`, with the reference fluid evaluated
/// at the conformal state and the critical enhancement from the fluid of
/// interest with the Olchowy-Sengers struct DEFAULTS (`parse_ECS_conductivity`
/// never fills the critical block — the JSON `q_D` key is unread).
#[allow(clippy::too_many_arguments)]
fn conductivity_ecs(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    reference: &EcsRef,
    viscosity_model: Option<&ViscosityModel>,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
    psi_a: &[f64],
    psi_t: &[f64],
    psi_rhomolar_reducing: f64,
    f_int_a: &[f64],
    f_int_t: &[f64],
    f_int_t_reducing: f64,
) -> Result<f64> {
    let m = eos.molar_mass;
    let m_kmol = m * 1000.0;
    let m0 = reference.eos.molar_mass;
    let (tc, rhocmolar) = crit_of(fluid);
    let (tc0, rhocmolar0) = crit_of(reference.fluid);
    let r_u = eos.gas_constant;
    let r = r_u / m;
    let r_kjkgk = r_u / m_kmol;

    let mut psi = 0.0;
    for i in 0..psi_a.len() {
        psi += psi_a[i] * (rhomolar / psi_rhomolar_reducing).powf(psi_t[i]);
    }
    let mut fint = 0.0;
    for i in 0..f_int_a.len() {
        fint += f_int_a[i] * (t / f_int_t_reducing).powf(f_int_t[i]);
    }

    // Dilute viscosity of the fluid of interest [uPa-s]: kinetic theory with
    // the fluid's own L-J parameters (its ECS viscosity block, or the
    // default_transport estimation).
    let (sigma, eps) = fluid_lennard_jones(eos, viscosity_model);
    let eta_dilute_upas = kinetic_theory_dilute(eos, t, sigma, eps) * 1e6;

    // cp0 (ideal gas), mass-based: cp0/R = 1 - tau^2 * d2alpha0/dtau2.
    let tau = eos.t_reducing / t;
    let delta = rhomolar / eos.rhomolar_reducing;
    let a0 = eos.alpha0_all(tau, delta);
    let cp0 = eos.gas_constant * (1.0 - tau * tau * a0.d02) / m;

    let lambda_int = fint * eta_dilute_upas * (cp0 - 2.5 * r) / 1e3;
    let lambda_dilute = 15.0e-3 / 4.0 * r_kjkgk * eta_dilute_upas;

    // Conformal state of the reference fluid.
    let dd = eos.alphar_all(tau, delta);
    let alphar_target = dd.d00;
    let z_target = 1.0 + delta * dd.d10;
    let mut t0 = t / (tc / tc0);
    let mut rhomolar0 = rhomolar * (rhocmolar0 / rhocmolar);
    conformal_state_solver(
        reference.eos,
        tc0,
        rhocmolar0,
        alphar_target,
        z_target,
        &mut t0,
        &mut rhomolar0,
    )
    .map_err(|e| Error::Value(format!("Conformal state solver failed; error was: {e}")))?;

    // Reference residual (background) conductivity at (rho0*psi, T0) —
    // upstream `calc_conductivity_background` is the residual term only.
    let rc = reference.conductivity.ok_or_else(|| {
        Error::NotImplemented("ECS reference fluid's conductivity model is unavailable".into())
    })?;
    let ConductivityModel::Structured(rcs) = rc else {
        return Err(Error::NotImplemented(
            "ECS reference fluid's conductivity model is not structured".into(),
        ));
    };
    let rho_ref = rhomolar0 * psi;
    let lambda_resid = conductivity_residual(reference.eos, rcs, t0, rho_ref);

    let f = t / t0;
    let h = rhomolar0 / rhomolar;
    let f_lambda = f.sqrt() * h.powf(-2.0 / 3.0) * (m0 / m).sqrt();

    // Critical enhancement of the fluid of interest, pure struct defaults.
    let lambda_critical = olchowy_sengers(
        eos,
        fluid,
        viscosity_model,
        t,
        rhomolar,
        p,
        ecs,
        1.3806488e-23,
        1.03,
        1.239,
        0.63,
        0.0496,
        1.94e-10,
        2e9,
        f64::NAN,
    )?;

    Ok(lambda_int + lambda_dilute + lambda_resid * f_lambda + lambda_critical)
}

#[allow(clippy::too_many_arguments)]
fn conductivity_structured(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    c: &Conductivity,
    viscosity_model: Option<&ViscosityModel>,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
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
            let eta0_upas = fluid_dilute_viscosity(eos, viscosity_model, t)? * 1e6;
            let tau = eos.t_reducing / t;
            let mut summer = a[0] * eta0_upas;
            for i in 1..a.len() {
                summer += a[i] * tau.powf(at[i]);
            }
            summer
        }
        ConductivityDilute::Hardcoded { name } => match *name {
            "CarbonDioxideHuberJPCRD2016" => {
                let tau = eos.t_reducing / t;
                let l = [0.0151874307, 0.0280674040, 0.0228564190, -0.00741624210];
                // Huber 2016 Eq. (3), in mW/m/K.
                let lambda_0 = tau.powf(-0.5)
                    / (l[0] + l[1] * tau + l[2] * tau.powf(2.0) + l[3] * tau.powf(3.0));
                lambda_0 / 1000.0
            }
            "Ethane" => {
                let e_k = 245.0;
                let tau = 305.33 / t;
                let tstar = t / e_k;
                let fint = 1.7104147 - 0.6936482 / tstar;
                let eta0_upas = fluid_dilute_viscosity(eos, viscosity_model, t)? * 1e6;
                let a0 = eos.alpha0_all(eos.t_reducing / t, rhomolar / eos.rhomolar_reducing);
                0.276505e-3 * eta0_upas * (3.75 - fint * (tau * tau * a0.d02 + 1.5))
            }
            other => {
                return Err(Error::NotImplemented(format!(
                    "hardcoded dilute conductivity [{other}] is not ported yet"
                )));
            }
        },
    };

    let residual = conductivity_residual(eos, c, t, rhomolar);

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
        }) => olchowy_sengers(
            eos,
            fluid,
            viscosity_model,
            t,
            rhomolar,
            p,
            ecs,
            *k,
            *r0,
            *gamma,
            *nu,
            *big_gamma,
            *zeta0,
            *qd,
            *t_ref,
        )?,
        Some(ConductivityCritical::Hardcoded { name }) => match *name {
            // Upstream maps the JSON "None" tag to CONDUCTIVITY_CRITICAL_NONE
            // (FluidLibrary.h:902-903) and evaluates it as critical = 0.0
            // (HelmholtzEOSMixtureBackend.cpp:1023-1025) — the pseudo-pure
            // blends (R404A/R407C/R410A/R507A) carry it.
            "None" => 0.0,
            "Ammonia" => conductivity_critical_ammonia(t, rhomolar * eos.molar_mass),
            "R123" => {
                let tau = eos.t_reducing / t;
                let delta = rhomolar / eos.rhomolar_reducing;
                let (a13, a14, a15) = (0.486742e-2, -100.0, -7.08535);
                a13 * (a14 * (tau - 1.0).powf(4.0) + a15 * (delta - 1.0).powf(2.0)).exp()
            }
            other => {
                return Err(Error::NotImplemented(format!(
                    "hardcoded critical conductivity [{other}] is not ported yet"
                )));
            }
        },
    };

    Ok(dilute + residual + critical)
}

/// The residual conductivity term (upstream `calc_conductivity_background`
/// — the background is the residual contribution only), reused by the ECS
/// reference evaluation.
fn conductivity_residual(eos: &HelmholtzEos, c: &Conductivity, t: f64, rhomolar: f64) -> f64 {
    match &c.residual {
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
    }
}

// ---------------------------------------------------------------------------
// Hardcoded models (upstream TransportRoutines::*_hardcoded), ported
// line-for-line: constants stay verbatim (excessive-precision literals
// included) and index loops mirror the upstream loops.
// ---------------------------------------------------------------------------

/// Upstream `conductivity_critical_simplified_Olchowy_Sengers`, shared by
/// structured critical blocks (fluid parameters) and ECS conductivity
/// (upstream struct defaults).
#[allow(clippy::too_many_arguments)]
fn olchowy_sengers(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    viscosity_model: Option<&ViscosityModel>,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
    k: f64,
    r0: f64,
    gamma: f64,
    nu: f64,
    big_gamma: f64,
    zeta0: f64,
    qd: f64,
    t_ref: f64,
) -> Result<f64> {
    let tc = eos.t_reducing;
    let rhoc = eos.rhomolar_reducing;
    let pcrit = fluid.eos.reducing.p;
    let tref = if t_ref.is_finite() { t_ref } else { 1.5 * tc };

    let delta = rhomolar / rhoc;
    let tau = tc / t;
    let dd = eos.alphar_all(tau, delta);
    let dp_drho = eos.gas_constant * t * (1.0 + 2.0 * delta * dd.d10 + delta * delta * dd.d20);
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
        return Ok(0.0);
    }
    let v = viscosity_model.ok_or_else(|| {
        Error::NotImplemented(
            "the Olchowy-Sengers enhancement needs the fluid's (unported) viscosity model".into(),
        )
    })?;
    let zeta = zeta0 * (num / big_gamma).powf(nu / gamma);
    let cp = eos.cpmolar(t, rhomolar);
    let cv = eos.cvmolar(t, rhomolar);
    let mu = viscosity(eos, fluid, v, t, rhomolar, p, ecs)?;
    let pi = std::f64::consts::PI;
    let omega_tilde = 2.0 / pi * ((cp - cv) / cp * (zeta * qd).atan() + cv / cp * (zeta * qd));
    let omega_tilde0 = 2.0 / pi
        * (1.0
            - (-1.0 / (1.0 / (qd * zeta) + 1.0 / 3.0 * (zeta * qd) * (zeta * qd) / delta / delta))
                .exp());
    Ok(rhomolar * cp * r0 * k * t / (6.0 * pi * mu * zeta) * (omega_tilde - omega_tilde0))
}

/// The fluid's dilute viscosity (upstream `calc_viscosity_dilute` through
/// the model wrapper) — consumed by eta0_and_poly and the Ethane dilute
/// conductivity.
fn fluid_dilute_viscosity(
    eos: &HelmholtzEos,
    model: Option<&ViscosityModel>,
    t: f64,
) -> Result<f64> {
    match model {
        Some(ViscosityModel::Structured(v)) => viscosity_dilute(eos, v, t),
        Some(ViscosityModel::Hardcoded { name }) => Err(Error::NotImplemented(format!(
            "dilute viscosity of hardcoded model [{name}] is not separable"
        ))),
        Some(
            ViscosityModel::Chung { .. }
            | ViscosityModel::RhosrCs { .. }
            | ViscosityModel::Ecs { .. },
        ) => Err(Error::NotImplemented(
            "dilute viscosity of Chung/rhosr/ECS models is not separable".into(),
        )),
        None => Err(Error::NotImplemented(
            "this conductivity needs the fluid's (unported) viscosity model".into(),
        )),
    }
}

/// IAPWS viscosity helper (upstream `visc_Helper`): dilute + finite-density
/// parts in reduced units.
fn visc_helper(tbar: f64, rhobar: f64) -> (f64, f64) {
    let mubar_0 = 100.0 * tbar.sqrt()
        / (1.67752 + 2.20462 / tbar + 0.6366564 / pow_int(tbar, 2) - 0.241605 / pow_int(tbar, 3));
    let mut h = [[0.0f64; 7]; 6];
    h[0][0] = 5.20094e-1;
    h[1][0] = 8.50895e-2;
    h[2][0] = -1.08374;
    h[3][0] = -2.89555e-1;
    h[0][1] = 2.22531e-1;
    h[1][1] = 9.99115e-1;
    h[2][1] = 1.88797;
    h[3][1] = 1.26613;
    h[5][1] = 1.20573e-1;
    h[0][2] = -2.81378e-1;
    h[1][2] = -9.06851e-1;
    h[2][2] = -7.72479e-1;
    h[3][2] = -4.89837e-1;
    h[4][2] = -2.57040e-1;
    h[0][3] = 1.61913e-1;
    h[1][3] = 2.57399e-1;
    h[0][4] = -3.25372e-2;
    h[3][4] = 6.98452e-2;
    h[4][5] = 8.72102e-3;
    h[3][6] = -4.35673e-3;
    h[5][6] = -5.93264e-4;
    let mut sum = 0.0;
    for (i, row) in h.iter().enumerate() {
        for (j, hij) in row.iter().enumerate() {
            sum += pow_int(1.0 / tbar - 1.0, i as i32) * (hij * pow_int(rhobar - 1.0, j as i32));
        }
    }
    (mubar_0, (rhobar * sum).exp())
}

/// IAPWS 2008 water viscosity (upstream `viscosity_water_hardcoded`).
#[allow(clippy::many_single_char_names)]
fn viscosity_water_hardcoded(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let (x_mu, qc, qd, nu, gamma, zeta_0, lambda_0, tbar_r) =
        (0.068, 1.0 / 1.9, 1.0 / 1.1, 0.630, 1.239, 0.13, 0.06, 1.5);
    let pstar = 22.064e6;
    let tstar = 647.096;
    let rhostar = 322.0;
    let tbar = t / tstar;
    let rhobar = rhomolar * eos.molar_mass / rhostar;
    let r_water = eos.gas_constant / eos.molar_mass;

    let (mubar_0, mubar_1) = visc_helper(tbar, rhobar);

    // Critical enhancement. Upstream sets the local delta := rhobar (the
    // MASS-scaled value) while the state derivatives come from the cached
    // molar-delta state; both are mirrored exactly.
    let delta = rhobar;
    let state_tau = eos.t_reducing / t;
    let state_delta = rhomolar / eos.rhomolar_reducing;
    let d = eos.alphar_all(state_tau, state_delta);
    let drhodp = 1.0 / (r_water * t * (1.0 + 2.0 * delta * d.d10 + delta * delta * d.d20));
    let drhobar_dpbar = pstar / rhostar * drhodp;
    let tau = 1.0 / tbar_r;
    let dref = eos.alphar_all(tau, delta);
    let drhodp_r = 1.0
        / (r_water * tbar_r * tstar * (1.0 + 2.0 * rhobar * dref.d10 + delta * delta * dref.d20));
    let drhobar_dpbar_r = pstar / rhostar * drhodp_r;

    let mut delta_chibar = rhobar * (drhobar_dpbar - drhobar_dpbar_r * tbar_r / tbar);
    if delta_chibar < 0.0 {
        delta_chibar = 0.0;
    }
    let zeta = zeta_0 * (delta_chibar / lambda_0).powf(nu / gamma);
    let y = if zeta < 0.3817016416 {
        1.0 / 5.0
            * qc
            * zeta
            * pow_int(qd * zeta, 5)
            * (1.0 - qc * zeta + pow_int(qc * zeta, 2) - 765.0 / 504.0 * pow_int(qd * zeta, 2))
    } else {
        let psi_d = (1.0 + pow_int(qd * zeta, 2)).powf(-1.0 / 2.0).acos();
        let w = ((qc * zeta - 1.0) / (qc * zeta + 1.0)).abs().sqrt() * (psi_d / 2.0).tan();
        let l = if qc * zeta > 1.0 {
            ((1.0 + w) / (1.0 - w)).ln()
        } else {
            2.0 * w.abs().atan()
        };
        1.0 / 12.0 * (3.0 * psi_d).sin() - 1.0 / (4.0 * qc * zeta) * (2.0 * psi_d).sin()
            + 1.0 / pow_int(qc * zeta, 2) * (1.0 - 5.0 / 4.0 * pow_int(qc * zeta, 2)) * psi_d.sin()
            - 1.0 / pow_int(qc * zeta, 3)
                * ((1.0 - 3.0 / 2.0 * pow_int(qc * zeta, 2)) * psi_d
                    - (pow_int(qc * zeta, 2) - 1.0).abs().powf(3.0 / 2.0) * l)
    };
    let mubar_2 = (x_mu * y).exp();
    (mubar_0 * mubar_1 * mubar_2) / 1e6
}

/// IAPWS 2020 heavy-water viscosity.
fn viscosity_heavywater_hardcoded(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tbar = t / 643.847;
    let rhobar = rhomolar * eos.molar_mass / 358.0;
    let a = [1.000000, 0.940695, 0.578377, -0.202044];
    let i_idx = [
        0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2, 3, 0, 1, 3, 5, 0, 1, 5, 3,
    ];
    let j_idx = [
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6,
    ];
    let bij = [
        0.4864192,
        -0.2448372,
        -0.8702035,
        0.8716056,
        -1.051126,
        0.3458395,
        0.3509007,
        1.315436,
        1.297752,
        1.353448,
        -0.2847572,
        -1.037026,
        -1.287846,
        -0.02148229,
        0.07013759,
        0.4660127,
        0.2292075,
        -0.4857462,
        0.01641220,
        -0.02884911,
        0.1607171,
        -0.009603846,
        -0.01163815,
        -0.008239587,
        0.004559914,
        -0.003886659,
    ];
    let mu0 = tbar.sqrt() / (a[0] + a[1] / tbar + a[2] / pow2(tbar) + a[3] / pow3(tbar));
    let mut summer = 0.0;
    for i in 0..26 {
        summer += bij[i]
            * (1.0 / tbar - 1.0).powf(f64::from(i_idx[i]))
            * (rhobar - 1.0).powf(f64::from(j_idx[i]));
    }
    let mu1 = (rhobar * summer).exp();
    55.2651e-6 * (mu0 * mu1)
}

/// Arp/McCarty/Friend helium viscosity (NIST TN 1334).
fn viscosity_helium_hardcoded(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let rho = rhomolar * eos.molar_mass / 1000.0; // [g/cm^3]
    let x = if t <= 300.0 { t.ln() } else { 300.0f64.ln() };
    let b = -47.5295259 / x + 87.6799309 - 42.0741589 * x + 8.33128289 * x * x
        - 0.589252385 * x * x * x;
    let c =
        547.309267 / x - 904.870586 + 431.404928 * x - 81.4504854 * x * x + 5.37008433 * x * x * x;
    let d =
        -1684.39324 / x + 3331.08630 - 1632.19172 * x + 308.804413 * x * x - 20.2936367 * x * x * x;
    let eta_0_slash = -0.135311743 / x + 1.00347841 + 1.20654649 * x - 0.149564551 * x * x
        + 0.012520841 * x * x * x;
    let eta_e_slash = rho * b + rho * rho * c + rho * rho * rho * d;
    if t <= 100.0 {
        let ln_eta = eta_0_slash + eta_e_slash;
        ln_eta.exp() / 10.0 / 1e6
    } else {
        let ln_eta = eta_0_slash + eta_e_slash;
        let eta_0 = 196.0 * t.powf(0.71938) * (12.451 / t - 295.67 / t / t - 4.1249).exp();
        (ln_eta.exp() + eta_0 - eta_0_slash.exp()) / 10.0 / 1e6
    }
}

/// Shan R23 viscosity.
fn viscosity_r23_hardcoded(t: f64, rhomolar: f64) -> f64 {
    let (c1, c2, delta_gstar, rho_l, rhocbar, tc, delta_eta_max, ru, molar_mass) = (
        1.3163, 0.1832, 771.23, 32.174, 7.5114, 299.2793, 3.967, 8.31451, 70.014,
    );
    let a = [0.4425728, -0.5138403, 0.1547566, -0.02821844, 0.001578286];
    let (e_k, sigma) = (243.91, 0.4278);
    let tstar = t / e_k;
    let log_tstar = tstar.ln();
    let omega = (a[0]
        + a[1] * log_tstar
        + a[2] * log_tstar.powf(2.0)
        + a[3] * log_tstar.powf(3.0)
        + a[4] * log_tstar.powf(4.0))
    .exp();
    let eta_dg = 1.25 * 0.021357 * (molar_mass * t).sqrt() / (sigma * sigma * omega); // uPa-s
    let rhobar = rhomolar / 1000.0; // [mol/L]
    let eta_l = c2 * (rho_l * rho_l) / (rho_l - rhobar)
        * t.sqrt()
        * (rhobar / (rho_l - rhobar) * delta_gstar / (ru * t)).exp();
    let chi = rhobar - rhocbar;
    let tau = t - tc;
    let delta_eta_c =
        4.0 * delta_eta_max / ((chi.exp() + (-chi).exp()) * (tau.exp() + (-tau).exp()));
    (((rho_l - rhobar) / rho_l).powf(c1) * eta_dg + (rhobar / rho_l).powf(c1) * eta_l + delta_eta_c)
        / 1e6
}

/// Xiang/Laesecke methanol viscosity.
#[allow(clippy::many_single_char_names)]
fn viscosity_methanol_hardcoded(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let epsilon_over_k = 577.87;
    let sigma0: f64 = 0.3408e-9;
    let delta = 0.4575;
    let n_a = 6.02214129e23;
    let m = 32.04216; // kg/kmol
    let tstar = t / epsilon_over_k;
    let rhor = rhomolar * eos.molar_mass / 273.0;
    let tr = t / 512.6;

    let (b_eta, c_eta) = {
        let b = [
            -19.572881,
            219.73999,
            -1015.3226,
            2471.01251,
            -3375.1717,
            2491.6597,
            -787.26086,
            14.085455,
            -0.34664158,
        ];
        let bt = [0.0, -0.25, -0.5, -0.75, -1.0, -1.25, -1.5, -2.5, -5.5];
        let mut summer = 0.0;
        for i in 0..9 {
            summer += b[i] * tstar.powf(bt[i]);
        }
        let b_eta = n_a * sigma0.powf(3.0) * summer;
        let c = [1.86222085e-3, 9.990338];
        let c_eta_star = c[0] * tstar.powf(3.0) * (c[1] * tstar.powf(-0.5)).exp();
        let c_eta = (n_a * sigma0.powf(3.0)).powf(2.0) * c_eta_star;
        (b_eta, c_eta)
    };

    let eta_g = 1.0 + b_eta * rhomolar + c_eta * rhomolar * rhomolar;
    let a = [
        1.16145, -0.14874, 0.52487, -0.77320, 2.16178, -2.43787, 0.95976e-3, 0.10225, -0.97346,
        0.10657, -0.34528, -0.44557, -2.58055,
    ];
    let d = [
        -1.181909,
        0.5031030,
        -0.6268461,
        0.5169312,
        -0.2351349,
        5.3980235e-2,
        -4.9069617e-3,
    ];
    let e = [
        0.0,
        4.018368,
        -4.239180,
        2.245110,
        -0.5750698,
        2.3021026e-2,
        2.5696775e-2,
        -6.8372749e-3,
        7.2707189e-4,
        -2.9255711e-5,
    ];
    let omega_22_star_lj =
        a[0] * tstar.powf(a[1]) + a[2] * (a[3] * tstar).exp() + a[4] * (a[5] * tstar).exp();
    let omega_22_star_delta =
        a[7] * tstar.powf(a[8]) + a[9] * (a[10] * tstar).exp() + a[11] * (a[12] * tstar).exp();
    let omega_22_star_sm = omega_22_star_lj
        * (1.0 + delta * delta / (1.0 + a[6] * f64::powf(delta, 6.0)) * omega_22_star_delta);
    let eta_0 = 2.66957e-26 * (m * t).sqrt() / (sigma0.powf(2.0) * omega_22_star_sm);

    let mut summerd = 0.0;
    for (i, di) in d.iter().enumerate() {
        summerd += di / tr.powi(i as i32);
    }
    for (j, ej) in e.iter().enumerate().skip(1) {
        summerd += ej * rhor.powi(j as i32);
    }
    let sigmac = 0.7193422e-9;
    let sigma_hs = summerd * sigmac;
    let b = 2.0 * std::f64::consts::PI * n_a * sigma_hs.powf(3.0) / 3.0;
    let zeta = b * rhomolar / 4.0;
    let g_sigma_hs = (1.0 - 0.5 * zeta) / (1.0 - zeta).powf(3.0);
    let eta_e =
        1.0 / g_sigma_hs + 0.8 * b * rhomolar + 0.761 * g_sigma_hs * (b * rhomolar).powf(2.0);
    let f = 1.0 / (1.0 + (5.0 * (rhor - 1.0)).exp());
    eta_0 * (f * eta_g + (1.0 - f) * eta_e)
}

/// Cao (JPCRD 2016) m-xylene viscosity.
fn viscosity_m_xylene_hardcoded(t: f64, rhomolar: f64) -> f64 {
    let d = [-0.268950, -0.0290018, 0.0, 14.7728, 17.1128];
    let n = [6.8, 3.3, 22.0, 0.6, 0.4];
    let e = [0.320971, 0.0, 1.72866e-10, -18.9852, 0.0];
    let k = [0.3, 0.0, 3.2];
    let tr = t / 616.89;
    let rhor = rhomolar / 1000.0 / 2.665;
    let (a0, b0, c0) = (-1.4933, 473.2, -57033.0);
    let ln_seta = a0 + b0 / t + c0 / (t * t);
    let eta0 = 0.22115 * t.sqrt() / ln_seta.exp();
    let (a1, b1, c1) = (13.2814, -10862.4, 1664060.0);
    let rho_moll = rhomolar / 1000.0;
    let eta1 = (a1 + b1 / t + c1 / (t * t)) * rho_moll;
    let f = (d[0] + e[0] * tr.powf(-k[0])) * rhor.powf(n[0])
        + d[1] * rhor.powf(n[1])
        + e[2] * rhor.powf(n[2]) / tr.powf(k[2])
        + (d[3] * rhor + e[3] * tr) * rhor.powf(n[3])
        + d[4] * rhor.powf(n[4]);
    let delta_eta = rhor.powf(2.0 / 3.0) * tr.sqrt() * f;
    (eta0 + eta1 + delta_eta) / 1e6
}

/// Cao (JPCRD 2016) o-xylene viscosity.
fn viscosity_o_xylene_hardcoded(t: f64, rhomolar: f64) -> f64 {
    let d = [-2.05581e-3, 2.38762, 0.0, 10.4497, 15.9587];
    let n = [10.3, 3.3, 25.0, 0.7, 0.4];
    let e = [2.65651e-3, 0.0, 1.77616e-12, -18.2446, 0.0];
    let k = [0.8, 0.0, 4.4];
    let tr = t / 630.259;
    let rhor = rhomolar / 1000.0 / 2.6845;
    let (a0, b0, c0) = (-1.4933, 473.2, -57033.0);
    let ln_seta = a0 + b0 / t + c0 / (t * t);
    let eta0 = 0.22225 * t.sqrt() / ln_seta.exp();
    let (a1, b1, c1) = (13.2814, -10862.4, 1664060.0);
    let rho_moll = rhomolar / 1000.0;
    let eta1 = (a1 + b1 / t + c1 / (t * t)) * rho_moll;
    let f = (d[0] + e[0] * tr.powf(-k[0])) * rhor.powf(n[0])
        + d[1] * rhor.powf(n[1])
        + e[2] * rhor.powf(n[2]) / tr.powf(k[2])
        + (d[3] * rhor + e[3] * tr) * rhor.powf(n[3])
        + d[4] * rhor.powf(n[4]);
    let delta_eta = rhor.powf(2.0 / 3.0) * tr.sqrt() * f;
    (eta0 + eta1 + delta_eta) / 1e6
}

/// Balogun (JPCRD 2016) p-xylene viscosity.
fn viscosity_p_xylene_hardcoded(t: f64, rhomolar: f64) -> f64 {
    let tr = t / 616.168;
    let rhor = rhomolar / 1000.0 / 2.69392;
    let (a0, b0, c0) = (-1.4933, 473.2, -57033.0);
    let ln_seta = a0 + b0 / t + c0 / (t * t);
    let eta0 = 0.22005 * t.sqrt() / ln_seta.exp();
    let (a1, b1, c1) = (13.2814, -10862.4, 1664060.0);
    let rho_moll = rhomolar / 1000.0;
    let eta1 = (a1 + b1 / t + c1 / (t * t)) * rho_moll;
    let sum1 = 122.919 * rhor.powf(1.5) - 282.329 * rhor.powf(2.0) + 279.348 * rhor.powf(3.0)
        - 146.776 * rhor.powf(4.0)
        + 28.361 * rhor.powf(5.0)
        - 0.004585 * rhor.powf(11.0);
    let sum2 = 15.337 * rhor.powf(1.5) - 0.0004382 * rhor.powf(11.0) + 0.00002307 * rhor.powf(15.0);
    let delta_eta = rhor.powf(2.0 / 3.0) * (sum1 + 1.0 / tr.sqrt() * sum2);
    (eta0 + eta1 + delta_eta) / 1e6
}

/// Friend (JPCRD 1991) ethane dilute viscosity.
fn viscosity_dilute_ethane(t: f64) -> f64 {
    let c = [
        0.0,
        -3.0328138281,
        16.918880086,
        -37.189364917,
        41.288861858,
        -24.615921140,
        8.9488430959,
        -1.8739245042,
        0.20966101390,
        -9.6570437074e-3,
    ];
    let e_k = 245.0;
    let tstar = t / e_k;
    let mut omega_2_2 = 0.0;
    for (i, ci) in c.iter().enumerate().skip(1) {
        omega_2_2 += ci * tstar.powf((i as f64 - 1.0) / 3.0 - 1.0);
    }
    12.0085 * tstar.sqrt() * omega_2_2 / 1e6
}

/// Friend (JPCRD 1991) ethane higher-order viscosity.
fn viscosity_ethane_higher_order(t: f64, rhomolar: f64) -> f64 {
    let r = [0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 1.0, 1.0];
    let s = [0.0, 0.0, 1.0, 0.0, 1.0, 1.5, 0.0, 2.0, 0.0, 1.0, 0.0, 1.0];
    let g = [
        0.0,
        0.47177003,
        -0.23950311,
        0.39808301,
        -0.27343335,
        0.35192260,
        -0.21101308,
        -0.00478579,
        0.07378129,
        -0.030435255,
        -0.30435286,
        0.001215675,
    ];
    let tau = 305.33 / t;
    let delta = rhomolar / 6870.0;
    let mut sum1 = 0.0;
    for i in 1..=9 {
        sum1 += g[i] * delta.powf(r[i]) * tau.powf(s[i]);
    }
    let mut sum2 = 0.0;
    for i in 10..=11 {
        sum2 += g[i] * delta.powf(r[i]) * tau.powf(s[i]);
    }
    15.977 * sum1 / (1.0 + sum2) / 1e6
}

/// Tariq (JPCRD 2014) cyclohexane dilute viscosity.
fn viscosity_dilute_cyclohexane(t: f64) -> f64 {
    let s_eta = (-1.5093 + 364.87 / t - 39537.0 / t.powf(2.0)).exp();
    0.19592 * t.sqrt() / s_eta / 1e6
}

/// Laesecke (JPCRD 2017) CO2 dilute viscosity.
fn viscosity_dilute_co2_laesecke(t: f64) -> f64 {
    let a = [
        1749.354893188350,
        -369.069300007128,
        5423856.34887691,
        -2.21283852168356,
        -269503.247933569,
        73145.021531826,
        5.34368649509278,
    ];
    let den = a[0]
        + a[1] * t.powf(1.0 / 6.0)
        + a[2] * (a[3] * t.powf(1.0 / 3.0)).exp()
        + (a[4] + a[5] * t.powf(1.0 / 3.0)) / t.powf(1.0 / 3.0).exp()
        + a[6] * t.sqrt();
    0.0010055 * t.sqrt() / den
}

/// Laesecke (JPCRD 2017) CO2 residual viscosity.
fn viscosity_co2_higher_order_laesecke(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let (c1, c2, gamma) = (0.360603235428487, 0.121550806591497, 8.06282737481277);
    // Upstream `Ttriple()` = sat_min_liquid.T at runtime; equal to the JSON
    // Ttriple for CO2 (216.592 K).
    let tt = 216.592;
    let rho_tl = 1178.53;
    let tr = t / tt;
    let rhor = rhomolar * eos.molar_mass / rho_tl;
    let eta_tl = rho_tl.powf(2.0 / 3.0) * (eos.gas_constant * tt).sqrt()
        / (eos.molar_mass.powf(1.0 / 6.0) * 84446887.43579945);
    eta_tl * (c1 * tr * rhor.powf(3.0) + (rhor.powf(2.0) + rhor.powf(gamma)) / (tr - c2))
}

/// Avgeri (JPCRD 2014) benzene residual viscosity.
fn viscosity_benzene_higher_order(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tr = t / 562.02;
    let rhor = rhomolar * eos.molar_mass / 304.792;
    let c = [
        -9.98945, 86.06260, 2.74872, 1.11130, -1.0, -134.1330, -352.473, 6.60989, 88.4174,
    ];
    1e-6 * rhor.powf(2.0 / 3.0)
        * tr.sqrt()
        * (c[0] * rhor.powf(2.0)
            + c[1] * rhor / (c[2] + c[3] * tr + c[4] * rhor)
            + (c[5] * rhor + c[6] * rhor.powf(2.0)) / (c[7] + c[8] * rhor.powf(2.0)))
}

/// Muzny (JCED 2013) hydrogen residual viscosity.
fn viscosity_hydrogen_higher_order(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tr = t / 33.145;
    let rhor = rhomolar * eos.molar_mass * 0.011;
    let c = [
        0.0,
        6.43449673e-6,
        4.56334068e-2,
        2.32797868e-1,
        9.58326120e-1,
        1.27941189e-1,
        3.63576595e-1,
    ];
    c[1] * rhor.powf(2.0)
        * (c[2] * tr + c[3] / tr + c[4] * rhor.powf(2.0) / (c[5] + tr) + c[6] * rhor.powf(6.0))
            .exp()
}

/// Avgeri (JPCRD 2015) toluene residual viscosity.
fn viscosity_toluene_higher_order(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tr = t / 591.75;
    let rhor = rhomolar * eos.molar_mass / 291.987;
    let c = [
        19.919216,
        -2.6557905,
        -135.904211,
        -7.9962719,
        -11.014795,
        -10.113817,
    ];
    1e-6 * rhor.powf(2.0 / 3.0)
        * tr.sqrt()
        * ((c[0] * rhor + c[1] * rhor.powf(4.0)) / tr
            + c[2] * rhor * rhor * rhor / (rhor * rhor + c[3] + c[4] * tr)
            + c[5] * rhor)
}

/// Michailidou (JPCRD 2013) hexane residual viscosity.
fn viscosity_hexane_higher_order(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tr = t / 507.82;
    let rhor = rhomolar * eos.molar_mass / 233.182;
    let c = [
        2.53402335 / 1e6,
        -9.724061002 / 1e6,
        0.469437316,
        158.5571631,
        72.42916856 / 1e6,
        10.60751253,
        8.628373915,
        -6.61346441,
        -2.212724566,
    ];
    rhor.powf(2.0 / 3.0)
        * tr.sqrt()
        * (c[0] / tr
            + c[1] / (c[2] + tr + c[3] * rhor * rhor)
            + c[4] * (1.0 + rhor)
                / (c[5] + c[6] * tr + c[7] * rhor + rhor * rhor + c[8] * rhor * tr))
}

/// Michailidou (JPCRD 2014) heptane residual viscosity.
fn viscosity_heptane_higher_order(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tr = t / 540.13;
    let rhor = rhomolar * eos.molar_mass / 232.0;
    let c = [
        0.0,
        22.15000 / 1e6,
        -15.00870 / 1e6,
        3.71791 / 1e6,
        77.72818 / 1e6,
        9.73449,
        9.51900,
        -6.34076,
        -2.51909,
    ];
    rhor.powf(2.0 / 3.0)
        * tr.sqrt()
        * (c[1] * rhor
            + c[2] * rhor.powf(2.0)
            + c[3] * rhor.powf(3.0)
            + c[4] * rhor / (c[5] + c[6] * tr + c[7] * rhor + rhor * rhor + c[8] * rhor * tr))
}

/// IAPWS 2011 water conductivity.
#[allow(clippy::many_single_char_names)]
fn conductivity_water_hardcoded(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    viscosity_model: Option<&ViscosityModel>,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
) -> Result<f64> {
    let l: [[f64; 6]; 5] = [
        [
            1.60397357,
            -0.646013523,
            0.111443906,
            0.102997357,
            -0.0504123634,
            0.00609859258,
        ],
        [
            2.33771842,
            -2.78843778,
            1.53616167,
            -0.463045512,
            0.0832827019,
            -0.00719201245,
        ],
        [
            2.19650529,
            -4.54580785,
            3.55777244,
            -1.40944978,
            0.275418278,
            -0.0205938816,
        ],
        [
            -1.21051378,
            1.60812989,
            -0.621178141,
            0.0716373224,
            0.0,
            0.0,
        ],
        [
            -2.7203370,
            4.57586331,
            -3.18369245,
            1.1168348,
            -0.19268305,
            0.012913842,
        ],
    ];
    let tstar = 647.096;
    let rhostar = 322.0;
    let pstar = 22064000.0;
    let lambdastar = 1e-3;
    let mustar = 1e-6;
    let r = 461.51805; // [J/kg/K]

    let tbar = t / tstar;
    let rhobar = rhomolar * eos.molar_mass / rhostar;

    let lambdabar_0 = tbar.sqrt()
        / (2.443221e-3 + 1.323095e-2 / tbar + 6.770357e-3 / tbar.powf(2.0)
            - 3.454586e-3 / tbar.powf(3.0)
            + 4.096266e-4 / tbar.powf(4.0));

    let mut sum = 0.0;
    for (i, row) in l.iter().enumerate() {
        for (j, lij) in row.iter().enumerate() {
            sum += lij * pow_int(1.0 / tbar - 1.0, i as i32) * pow_int(rhobar - 1.0, j as i32);
        }
    }
    let lambdabar_1 = (rhobar * sum).exp();

    let (nu, big_gamma, gamma, xi_0, lambda_0, tr_bar, qd_bar) =
        (0.630, 177.8514, 1.239, 0.13, 0.06, 1.5, 1.0 / 0.4);
    #[allow(clippy::approx_constant)]
    let pi = 3.141592654; // upstream's literal
    let delta = rhomolar / eos.rhomolar_reducing;

    let state_tau = eos.t_reducing / t;
    let d = eos.alphar_all(state_tau, delta);
    let drhodp = 1.0 / (r * t * (1.0 + 2.0 * rhobar * d.d10 + rhobar * rhobar * d.d20));
    let drhobar_dpbar = pstar / rhostar * drhodp;
    let dref = eos.alphar_all(1.0 / tr_bar, delta);
    let drhodp_trbar =
        1.0 / (r * tr_bar * tstar * (1.0 + 2.0 * rhobar * dref.d10 + delta * delta * dref.d20));
    let drhobar_dpbar_trbar = pstar / rhostar * drhodp_trbar;
    let cp = eos.cpmolar(t, rhomolar) / eos.molar_mass; // [J/kg/K]
    let cv = eos.cvmolar(t, rhomolar) / eos.molar_mass;
    let cpbar = cp / r;
    let v = viscosity_model.ok_or_else(|| {
        Error::NotImplemented("water conductivity needs the water viscosity model".into())
    })?;
    let mubar = viscosity(eos, fluid, v, t, rhomolar, p, ecs)? / mustar;
    let delta_chibar_t = rhobar * (drhobar_dpbar - drhobar_dpbar_trbar * tr_bar / tbar);
    let xi = if delta_chibar_t < 0.0 {
        0.0
    } else {
        xi_0 * (delta_chibar_t / lambda_0).powf(nu / gamma)
    };
    let y = qd_bar * xi;
    let kappa = cp / cv;
    let z = if y < 1.2e-7 {
        0.0
    } else {
        2.0 / (pi * y)
            * (((1.0 - 1.0 / kappa) * y.atan() + y / kappa)
                - (1.0 - (-1.0 / (1.0 / y + y * y / 3.0 / rhobar / rhobar)).exp()))
    };
    let lambdabar_2 = big_gamma * rhobar * cpbar * tbar / mubar * z;
    Ok((lambdabar_0 * lambdabar_1 + lambdabar_2) * lambdastar)
}

/// IAPWS 2021 heavy-water conductivity.
fn conductivity_heavywater_hardcoded(eos: &HelmholtzEos, t: f64, rhomolar: f64) -> f64 {
    let tbar = t / 643.847;
    let rhobar = rhomolar * eos.molar_mass / 358.0;
    let a = [1.00000, 37.3223, 22.5485, 13.0465, 0.0, -2.60735];
    let lambda0 = a[0]
        + a[1] * tbar
        + a[2] * pow2(tbar)
        + a[3] * pow3(tbar)
        + a[4] * pow4(tbar)
        + a[5] * pow5(tbar);
    let be = -2.506;
    let b = [-167.310, 483.656, -191.039, 73.0358, -7.57467];
    let delta_lambda = b[0] * (1.0 - (be * rhobar).exp())
        + b[1] * rhobar
        + b[2] * pow2(rhobar)
        + b[3] * pow3(rhobar)
        + b[4] * pow4(rhobar);
    let f_1 = (0.144847 * tbar + -5.64493 * pow2(tbar)).exp();
    let f_2 = (-2.80000 * pow2(rhobar - 1.0)).exp()
        - 0.080738543 * (-17.9430 * pow2(rhobar - 0.125698)).exp();
    let tau = tbar / ((tbar - 1.1).abs() + 1.1);
    let f_3 = 1.0 + (60.0 * (tau - 1.0) + 20.0).exp();
    let f_4 = 1.0 + (100.0 * (tau - 1.0) + 15.0).exp();
    let delta_lambda_c =
        35429.6 * f_1 * f_2 * (1.0 + pow2(f_2) * (5000.0e6 * pow4(f_1) / f_3 + 3.5 * f_2 / f_4));
    let delta_lambda_l = -741.112 * f_1.powf(1.2) * (1.0 - (-((rhobar / 2.5).powf(10.0))).exp());
    (lambda0 + delta_lambda + delta_lambda_c + delta_lambda_l) * 0.742128e-3
}

/// Hands/Arp helium conductivity.
#[allow(clippy::many_single_char_names)]
fn conductivity_helium_hardcoded(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    viscosity_model: Option<&ViscosityModel>,
    t: f64,
    rhomolar: f64,
    p: f64,
    ecs: Option<&EcsResolver>,
) -> Result<f64> {
    let rhoc = 68.0;
    let rho = rhomolar * eos.molar_mass; // [kg/m^3]
    let summer = 3.739232544 / t - 2.620316969e1 / t / t + 5.982252246e1 / t / t / t
        - 4.926397634e1 / t / t / t / t;
    let lambda_0 = 2.7870034e-3 * t.powf(7.034007057e-1) * summer.exp();
    let c = [
        1.862970530e-4,
        -7.275964435e-7,
        -1.427549651e-4,
        3.290833592e-5,
        -5.213335363e-8,
        4.492659933e-8,
        -5.924416513e-9,
        7.087321137e-6,
        -6.013335678e-6,
        8.067145814e-7,
        3.995125013e-7,
    ];
    let lambda_e = (c[0] + c[1] * t + c[2] * t.powf(1.0 / 3.0) + c[3] * t.powf(2.0 / 3.0)) * rho
        + (c[4] + c[5] * t.powf(1.0 / 3.0) + c[6] * t.powf(2.0 / 3.0)) * rho * rho * rho
        + (c[7] + c[8] * t.powf(1.0 / 3.0) + c[9] * t.powf(2.0 / 3.0) + c[10] / t)
            * rho
            * rho
            * (rho / rhoc).ln();

    let mut lambda_c = 0.0;
    if 3.5 < t && t < 12.0 {
        let (x0, e1, e2, beta, gamma, delta, rhoc_crit, tc, pc) = (
            0.392, 2.8461, 0.27156, 0.3554, 1.1743, 4.304, 69.158, 5.18992, 2.2746e5,
        );
        let delta_t = (1.0 - t / tc).abs();
        let delta_rho = (1.0 - rho / rhoc_crit).abs();
        let v = viscosity_model.ok_or_else(|| {
            Error::NotImplemented("helium conductivity needs the helium viscosity model".into())
        })?;
        let eta = viscosity(eos, fluid, v, t, rhomolar, p, ecs)?;
        // Isothermal compressibility 1/(rho*dp/drho) and dp/dT|rho from the
        // alpha derivatives.
        let tau = eos.t_reducing / t;
        let del = rhomolar / eos.rhomolar_reducing;
        let d = eos.alphar_all(tau, del);
        let dpdrho = eos.gas_constant * t * (1.0 + 2.0 * del * d.d10 + del * del * d.d20);
        let k_t = 1.0 / (rhomolar * dpdrho);
        let dpdt = rhomolar * eos.gas_constant * (1.0 + del * d.d10 - del * tau * d.d11);

        let w = (delta_t / 0.2).powf(2.0) + (delta_rho / 0.25).powf(2.0);
        let k_tbar = if w > 1.0 {
            k_t
        } else {
            let x = (delta_t / delta_rho).powf(1.0 / beta);
            let h = e1
                * (1.0 + x / x0)
                * (1.0 + e2 * (1.0 + x / x0).powf(2.0 / beta)).powf((gamma - 1.0) / (2.0 * beta));
            let dhdx = e1
                * (e2
                    * ((x + x0) / x0).powf(2.0 / beta)
                    * (gamma - 1.0)
                    * (e2 * ((x + x0) / x0).powf(2.0 / beta) + 1.0)
                        .powf((1.0 / 2.0) * (gamma - 1.0) / beta)
                    + beta.powf(2.0)
                        * (e2 * ((x + x0) / x0).powf(2.0 / beta) + 1.0)
                            .powf((1.0 / 2.0) * (2.0 * beta + gamma - 1.0) / beta))
                / (beta.powf(2.0) * x0 * (e2 * ((x + x0) / x0).powf(2.0 / beta) + 1.0));
            let rhs = delta_rho.powf(delta - 1.0) * (delta * h - x / beta * dhdx);
            let k_tprime = 1.0 / (rhs * (rho / rhoc_crit).powf(2.0) * pc);
            w * k_t + (1.0 - w) * k_tprime
        };
        lambda_c = 3.4685233e-17 * 3.726229668 * k_tbar.sqrt() * t.powf(2.0) / rho / eta
            * dpdt.powf(2.0)
            * (-18.66 * delta_t.powf(2.0) - 4.25 * delta_rho.powf(4.0)).exp();
    }
    Ok(lambda_0 + lambda_e + lambda_c)
}

/// Shan R23 conductivity.
fn conductivity_r23_hardcoded(t: f64, rhomolar: f64) -> f64 {
    let (b1, b2, c1, c2, delta_gstar, rho_l, rhocbar, delta_lambda_max, ru, tc) = (
        -2.5370, 0.05366, 0.94215, 0.14914, 2508.58, 68.345, 7.5114, 25.0, 8.31451, 299.2793,
    );
    let lambda_dg = b1 + b2 * t;
    let rhobar = rhomolar / 1000.0;
    let lambda_l = c2 * (rho_l * rho_l) / (rho_l - rhobar)
        * t.sqrt()
        * (rhobar / (rho_l - rhobar) * delta_gstar / (ru * t)).exp();
    let chi = rhobar - rhocbar;
    let tau = t - tc;
    let delta_lambda_c =
        4.0 * delta_lambda_max / ((chi.exp() + (-chi).exp()) * (tau.exp() + (-tau).exp()));
    (((rho_l - rhobar) / rho_l).powf(c1) * lambda_dg
        + (rhobar / rho_l).powf(c1) * lambda_l
        + delta_lambda_c)
        / 1e3
}

/// Friend (JPCRD 1989) methane conductivity.
#[allow(clippy::many_single_char_names)]
fn conductivity_methane_hardcoded(
    eos: &HelmholtzEos,
    fluid: &FluidData,
    t: f64,
    rhomolar: f64,
) -> f64 {
    let delta = rhomolar / 10139.0;
    let tau = 190.55 / t;

    // Viscosity formulation from Friend, JPCRD, 1989 (self-contained).
    let c = [
        0.0,
        -3.0328138281,
        16.918880086,
        -37.189364917,
        41.288861858,
        -24.615921140,
        8.9488430959,
        -1.8739245042,
        0.20966101390,
        -9.6570437074e-3,
    ];
    let mut omega22_summer = 0.0;
    let tt = t / 174.0;
    for i in 1..=9 {
        omega22_summer += c[i] * tt.powf((i as f64 - 1.0) / 3.0 - 1.0);
    }
    let eta_dilute = 10.50 * tt.sqrt() * omega22_summer;
    let re = [0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 1.0, 1.0];
    let se = [0.0, 0.0, 1.0, 0.0, 1.0, 1.5, 0.0, 2.0, 0.0, 1.0, 0.0, 1.0];
    let ge = [
        0.0,
        0.41250137,
        -0.14390912,
        0.10366993,
        0.40287464,
        -0.24903524,
        -0.12953131,
        0.06575776,
        0.02566628,
        -0.03716526,
        -0.38798341,
        0.03533815,
    ];
    let mut summer1 = 0.0;
    let mut summer2 = 0.0;
    for i in 1..=9 {
        summer1 += ge[i] * delta.powf(re[i]) * tau.powf(se[i]);
    }
    for i in 10..=11 {
        summer2 += ge[i] * delta.powf(re[i]) * tau.powf(se[i]);
    }
    let eta_residual = 12.149 * summer1 / (1.0 + summer2);
    let eta = eta_residual + eta_dilute;

    // Dilute conductivity.
    let f_int = 1.458850 - 0.4377162 / tt;
    let state_tau = eos.t_reducing / t;
    let state_delta = rhomolar / eos.rhomolar_reducing;
    let a0 = eos.alpha0_all(state_tau, state_delta);
    let lambda_dilute =
        0.51828 * eta_dilute * (3.75 - f_int * (state_tau * state_tau * a0.d02 + 1.5));

    // Residual conductivity.
    let rl = [0.0, 1.0, 3.0, 4.0, 4.0, 5.0, 5.0, 2.0];
    let sl = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0];
    let jl = [
        0.0,
        2.4149207,
        0.55166331,
        -0.52837734,
        0.073809553,
        0.24465507,
        -0.047613626,
        1.5554612,
    ];
    let mut summer = 0.0;
    for i in 1..=6 {
        summer += jl[i] * delta.powf(rl[i]) * tau.powf(sl[i]);
    }
    // Upstream: T_critical()/rhomolar_critical() are the superancillary
    // NUMERICAL values; the saturated-vapor density comes from the CLASSIC
    // rhoV ancillary.
    let (t_crit, rho_crit) = match &fluid.eos.superancillary {
        Some(sa) => (sa.t_crit_num, sa.rho_crit_num),
        None => (fluid.states.critical.t, fluid.states.critical.rhomolar),
    };
    let mut delta_sigma_star = 1.0;
    if t < t_crit && rhomolar < rho_crit {
        delta_sigma_star = crate::ancillary::evaluate(&fluid.ancillaries.rho_v, t) / rho_crit;
    }
    let lambda_residual = 6.29638 * (summer + jl[7] * pow2(delta) / delta_sigma_star);

    // Critical region.
    let tstar = 1.0 - 1.0 / tau;
    let rhostar = 1.0 - delta;
    let (f_t, f_rho, f_a) = (2.646, 2.678, -0.637);
    let f = (-f_t * tstar.abs().sqrt() - f_rho * pow2(rhostar) - f_a * rhostar).exp();
    let d = eos.alphar_all(state_tau, state_delta);
    let chi_from_eq19a = 0.28631 * delta * tau / (1.0 + 2.0 * delta * d.d10 + pow2(delta) * d.d20);
    let chi_t_star = if tstar.abs() < 0.03 {
        if rhostar.abs() < 1e-16 {
            let (lambda_cap, gamma) = (0.0801, 1.190);
            lambda_cap * tstar.abs().powf(-gamma)
        } else if rhostar.abs() < 0.03 {
            let (beta, w, s_c, e, a, b, r_c, q) =
                (0.355, -1.401, -6.098, 0.287, 3.352, 0.732, 0.535, 0.1133);
            let omega = w * tstar * rhostar.abs().powf(-1.0 / beta);
            let mut theta = 1.0;
            if tstar < -rhostar.abs().powf(-1.0 / beta) / s_c {
                theta = 1.0
                    + e * (1.0 + s_c * tstar * rhostar.abs().powf(-1.0 / beta)).powf(2.0 * beta);
            }
            q * rhostar.abs().powf(-a) * theta.powf(b) / (theta + omega * (theta + r_c))
        } else {
            chi_from_eq19a
        }
    } else {
        chi_from_eq19a
    };
    let lambda_critical = 91.855 / (eta * pow2(tau))
        * pow2(1.0 + delta * d.d10 - delta * tau * d.d11)
        * chi_t_star.powf(0.4681)
        * f;
    (lambda_dilute + lambda_residual + lambda_critical) * 0.001
}

/// Tufeu (1984) ammonia critical conductivity enhancement.
#[allow(clippy::many_single_char_names)]
fn conductivity_critical_ammonia(t_in: f64, rhomass: f64) -> f64 {
    let (tc, rhoc) = (405.4, 235.0);
    let (lambda_cap, nu, gamma, delta_cap, zeta_0_plus, a_zeta, gamma_0_plus) =
        (1.2, 0.63, 1.24, 0.50, 1.34e-10, 1.0, 0.423e-8);
    #[allow(clippy::approx_constant)]
    let pi = 3.141592654; // upstream's literal
    let k_b = 1.3806504e-23;

    let rho = rhomass;
    let t = ((t_in - tc) / tc).abs();
    let a_chi = a_zeta / 0.7;
    let eta_b = (2.60 + 1.6 * t) * 1e-5;
    let dpdt = (2.18 - 0.12 / (17.8 * t).exp()) * 1e5;
    let x_t = 0.61 * rhoc + 16.5 * t.ln();
    let delta_lambda_i = lambda_cap * (k_b * t_in * t_in)
        / (6.0 * pi * eta_b * (zeta_0_plus * t.powf(-nu) * (1.0 + a_zeta * t.powf(delta_cap))))
        * dpdt
        * dpdt
        * gamma_0_plus
        * t.powf(-gamma)
        * (1.0 + a_chi * t.powf(delta_cap));
    let delta_lambda_id = delta_lambda_i * (-36.0 * t * t).exp();
    if rho < 0.6 * rhoc {
        delta_lambda_id * (x_t * x_t) / (x_t * x_t + pow_int(0.6 * rhoc - 0.96 * rhoc, 2))
            * pow_int(rho, 2)
            / pow_int(0.6 * rhoc, 2)
    } else {
        delta_lambda_id * (x_t * x_t) / (x_t * x_t + pow_int(rho - 0.96 * rhoc, 2))
    }
}
