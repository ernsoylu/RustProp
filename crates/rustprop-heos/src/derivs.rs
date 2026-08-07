//! Generic partial derivatives in the (T, rho) basis — port of upstream
//! `src/AbstractState.cpp`'s `get_dT_drho`, `get_dT_drho_second_derivatives`,
//! `calc_first_partial_deriv` and `calc_second_partial_deriv`.
//!
//! Every derivation follows Thorade 2013 (DOI 10.1007/s12665-013-2394-z), as
//! upstream's comment records. These are the machinery the Tabular backends
//! build their grids from (`first_partial_deriv(iT, xkey, ykey)` etc.) and
//! the same functions behind PropsSI's `d(X)/d(Y)|Z` output strings.
//!
//! Scope: the parameter set the tabular tables consume (T, P, Dmolar/Dmass,
//! Hmolar/Hmass, Smolar/Smass, Umolar/Umass, Gmolar/Gmass, Tau, Delta).
//! Upstream's `get_dT_drho` additionally handles Cvmolar/Cpmolar/speed_sound
//! — reachable only through derivative-output strings, which are not ported;
//! those error loudly here.

use crate::alpha::HelmholtzEos;
use rustprop_core::params::Param;
use rustprop_core::{Error, Result};

/// `(dX/dT|rho, dX/drho|T)` for parameter `index` at the live (T, rhomolar).
pub fn get_dt_drho(eos: &HelmholtzEos, index: Param, t: f64, rhomolar: f64) -> Result<(f64, f64)> {
    let rhor = eos.rhomolar_reducing;
    let tr = eos.t_reducing;
    let dt_dtau = -t.powi(2) / tr;
    let r = eos.gas_constant;
    let delta = rhomolar / rhor;
    let tau = tr / t;
    let ar = eos.alphar_all(tau, delta);
    let a0 = eos.alpha0_all(tau, delta);
    let mm = eos.molar_mass;

    let (mut dt, mut drho) = match index {
        Param::T => (1.0, 0.0),
        Param::Dmolar => (0.0, 1.0),
        Param::Dmass => (0.0, mm),
        Param::P => {
            // dp/drho|T and dp/dT|rho
            let drho = r * t * (1.0 + 2.0 * delta * ar.d10 + delta.powi(2) * ar.d20);
            let dt = rhomolar * r * (1.0 + delta * ar.d10 - tau * delta * ar.d11);
            (dt, drho)
        }
        Param::Hmolar | Param::Hmass => {
            let dt = r
                * (-tau.powi(2) * (a0.d02 + ar.d02)
                    + (1.0 + delta * ar.d10 - tau * delta * ar.d11));
            let drho =
                t * r / rhomolar * (tau * delta * ar.d11 + delta * ar.d10 + delta.powi(2) * ar.d20);
            (dt, drho)
        }
        Param::Smolar | Param::Smass => {
            let dt = r / t * (-tau.powi(2) * (a0.d02 + ar.d02));
            let drho = r / rhomolar * (-(1.0 + delta * ar.d10 - tau * delta * ar.d11));
            (dt, drho)
        }
        Param::Umolar | Param::Umass => {
            let dt = r * (-tau.powi(2) * (a0.d02 + ar.d02));
            let drho = t * r / rhomolar * (tau * delta * ar.d11);
            (dt, drho)
        }
        Param::Gmolar | Param::Gmass => {
            let dtau_dt = 1.0 / dt_dtau;
            let dt = r * t * (a0.d01 + ar.d01 + delta * ar.d11) * dtau_dt
                + r * (1.0 + a0.d00 + ar.d00 + delta * ar.d10);
            let ddelta_drho = 1.0 / rhor;
            // NOTE upstream adds dalphar_dDelta twice here (once inside the
            // bracket, once as the trailing term) — reproduced verbatim.
            let drho = t * r * (a0.d10 + ar.d10 + delta * ar.d20 + ar.d10) * ddelta_drho;
            (dt, drho)
        }
        Param::Tau => (1.0 / dt_dtau, 0.0),
        Param::Delta => (0.0, 1.0 / rhor),
        other => {
            return Err(Error::NotImplemented(format!(
                "get_dT_drho for [{}] is not ported yet",
                other.short_name()
            )));
        }
    };
    // Mass-basis siblings scale both partials by 1/M.
    if matches!(
        index,
        Param::Hmass | Param::Smass | Param::Umass | Param::Gmass
    ) {
        drho /= mm;
        dt /= mm;
    }
    Ok((dt, drho))
}

/// `(d2X/dT2|rho, d2X/drho/dT, d2X/drho2|T)` for parameter `index`.
pub fn get_dt_drho_second(
    eos: &HelmholtzEos,
    index: Param,
    t: f64,
    rhomolar: f64,
) -> Result<(f64, f64, f64)> {
    let rhor = eos.rhomolar_reducing;
    let tr = eos.t_reducing;
    let r = eos.gas_constant;
    let delta = rhomolar / rhor;
    let tau = tr / t;
    let ar = eos.alphar_all(tau, delta);
    let a0 = eos.alpha0_all(tau, delta);
    let mm = eos.molar_mass;

    let (mut dt2, mut drho_dt, mut drho2) = match index {
        Param::T | Param::Dmass | Param::Dmolar => (0.0, 0.0, 0.0),
        Param::Tau => (2.0 * tr / t.powi(3), 0.0, 0.0),
        Param::Delta => (0.0, 0.0, 0.0),
        Param::P => {
            let drho2 = t * r / rhomolar
                * (2.0 * delta * ar.d10 + 4.0 * delta.powi(2) * ar.d20 + delta.powi(3) * ar.d30);
            let dt2 = rhomolar * r / t * (tau.powi(2) * delta * ar.d12);
            let drho_dt = r
                * (1.0 + 2.0 * delta * ar.d10 + delta.powi(2) * ar.d20
                    - 2.0 * delta * tau * ar.d11
                    - tau * delta.powi(2) * ar.d21);
            (dt2, drho_dt, drho2)
        }
        Param::Hmolar | Param::Hmass => {
            let drho2 =
                r * t * (delta / rhomolar).powi(2) * (tau * ar.d21 + 2.0 * ar.d20 + delta * ar.d30);
            let dt2 = r / t
                * tau.powi(2)
                * (tau * (a0.d03 + ar.d03) + 2.0 * (a0.d02 + ar.d02) + delta * ar.d12);
            let drho_dt = r / rhomolar
                * delta
                * (delta * ar.d20 - tau.powi(2) * ar.d12 + ar.d10
                    - tau * delta * ar.d21
                    - tau * ar.d11);
            (dt2, drho_dt, drho2)
        }
        Param::Smolar | Param::Smass => {
            let drho2 = r / rhomolar.powi(2)
                * (1.0 - delta.powi(2) * ar.d20 + tau * delta.powi(2) * ar.d21);
            let dt2 = r * (tau / t).powi(2) * (tau * (a0.d03 + ar.d03) + 3.0 * (a0.d02 + ar.d02));
            let drho_dt = r / (t * rhomolar) * (-tau.powi(2) * delta * ar.d12);
            (dt2, drho_dt, drho2)
        }
        Param::Umolar | Param::Umass => {
            let drho2 = r * t * tau * (delta / rhomolar).powi(2) * ar.d21;
            let dt2 = r / t * tau.powi(2) * (tau * (a0.d03 + ar.d03) + 2.0 * (a0.d02 + ar.d02));
            let drho_dt = r / rhomolar * (-tau.powi(2) * delta * ar.d12);
            (dt2, drho_dt, drho2)
        }
        other => {
            return Err(Error::NotImplemented(format!(
                "get_dT_drho_second_derivatives for [{}] is not ported yet",
                other.short_name()
            )));
        }
    };
    if matches!(index, Param::Hmass | Param::Smass | Param::Umass) {
        drho2 /= mm;
        drho_dt /= mm;
        dt2 /= mm;
    }
    Ok((dt2, drho_dt, drho2))
}

/// `calc_first_partial_deriv(Of, Wrt, Constant)` — d(Of)/d(Wrt) at constant
/// (Constant), by the (T, rho) Jacobian ratio.
pub fn first_partial_deriv(
    eos: &HelmholtzEos,
    of: Param,
    wrt: Param,
    constant: Param,
    t: f64,
    rhomolar: f64,
) -> Result<f64> {
    let (dof_dt, dof_drho) = get_dt_drho(eos, of, t, rhomolar)?;
    let (dwrt_dt, dwrt_drho) = get_dt_drho(eos, wrt, t, rhomolar)?;
    let (dconstant_dt, dconstant_drho) = get_dt_drho(eos, constant, t, rhomolar)?;
    Ok((dof_dt * dconstant_drho - dof_drho * dconstant_dt)
        / (dwrt_dt * dconstant_drho - dwrt_drho * dconstant_dt))
}

/// `calc_second_partial_deriv(Of1, Wrt1, Constant1, Wrt2, Constant2)`.
#[allow(clippy::too_many_arguments)]
pub fn second_partial_deriv(
    eos: &HelmholtzEos,
    of1: Param,
    wrt1: Param,
    constant1: Param,
    wrt2: Param,
    constant2: Param,
    t: f64,
    rhomolar: f64,
) -> Result<f64> {
    // First and second partials of the terms in the first derivative
    let (dof1_dt, dof1_drho) = get_dt_drho(eos, of1, t, rhomolar)?;
    let (dwrt1_dt, dwrt1_drho) = get_dt_drho(eos, wrt1, t, rhomolar)?;
    let (dconstant1_dt, dconstant1_drho) = get_dt_drho(eos, constant1, t, rhomolar)?;
    let (d2of1_dt2, d2of1_drhodt, d2of1_drho2) = get_dt_drho_second(eos, of1, t, rhomolar)?;
    let (d2wrt1_dt2, d2wrt1_drhodt, d2wrt1_drho2) = get_dt_drho_second(eos, wrt1, t, rhomolar)?;
    let (d2constant1_dt2, d2constant1_drhodt, d2constant1_drho2) =
        get_dt_drho_second(eos, constant1, t, rhomolar)?;

    // First derivatives of the terms in the second derivative
    let (dwrt2_dt, dwrt2_drho) = get_dt_drho(eos, wrt2, t, rhomolar)?;
    let (dconstant2_dt, dconstant2_drho) = get_dt_drho(eos, constant2, t, rhomolar)?;

    // Numerator and denominator of the first partial derivative term
    let n = dof1_dt * dconstant1_drho - dof1_drho * dconstant1_dt;
    let d = dwrt1_dt * dconstant1_drho - dwrt1_drho * dconstant1_dt;

    // Their derivatives with respect to rho (T held constant)
    let dndrho_t = dof1_dt * d2constant1_drho2 + d2of1_drhodt * dconstant1_drho
        - dof1_drho * d2constant1_drhodt
        - d2of1_drho2 * dconstant1_dt;
    let dddrho_t = dwrt1_dt * d2constant1_drho2 + d2wrt1_drhodt * dconstant1_drho
        - dwrt1_drho * d2constant1_drhodt
        - d2wrt1_drho2 * dconstant1_dt;

    // ... and with respect to T (rho held constant)
    let dndt_rho = dof1_dt * d2constant1_drhodt + d2of1_dt2 * dconstant1_drho
        - dof1_drho * d2constant1_dt2
        - d2of1_drhodt * dconstant1_dt;
    let dddt_rho = dwrt1_dt * d2constant1_drhodt + d2wrt1_dt2 * dconstant1_drho
        - dwrt1_drho * d2constant1_dt2
        - d2wrt1_drhodt * dconstant1_dt;

    let dderiv1_drho = (d * dndrho_t - n * dddrho_t) / d.powi(2);
    let dderiv1_dt = (d * dndt_rho - n * dddt_rho) / d.powi(2);

    Ok(
        (dderiv1_dt * dconstant2_drho - dderiv1_drho * dconstant2_dt)
            / (dwrt2_dt * dconstant2_drho - dwrt2_drho * dconstant2_dt),
    )
}
