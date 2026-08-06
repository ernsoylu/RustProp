//! Typed `PropsSI`-style dispatch for the IF97 engine (PLAN.md 2.4).
//!
//! Resolves two (parameter, value) inputs through
//! [`rustprop_core::generate_update_pair`] and routes to the same engine
//! functions CoolProp's IF97 backend calls. Scope is deliberately the input
//! pairs and outputs exercised by the Phase 2 golden fixtures; everything else
//! returns [`Error::NotImplemented`] until Phase 5 brings full `PropsSI`
//! parity.

use rustprop_core::{Error, InputPair, Param, Result, generate_update_pair};
use rustprop_if97 as if97;

/// `PropsSI`-style single-output call for IF97 water.
///
/// `props(output, name1, value1, name2, value2)` mirrors
/// `PropsSI(output, name1, value1, name2, value2, "IF97::Water")`.
pub fn props(output: Param, name1: Param, value1: f64, name2: Param, value2: f64) -> Result<f64> {
    // Trivial outputs need no state.
    match output {
        Param::TCritical => return Ok(if97::get_tcrit()),
        Param::PCritical => return Ok(if97::get_pcrit()),
        Param::RhomassCritical => return Ok(if97::get_rhocrit()),
        Param::TTriple => return Ok(if97::get_ttrip()),
        Param::PTriple => return Ok(if97::get_ptrip()),
        Param::TMin => return Ok(if97::get_tmin()),
        Param::TMax => return Ok(if97::get_tmax()),
        // Golden-verified: the CoolProp IF97 backend reports the triple
        // pressure (611.657 Pa) as PMIN, not IF97's Pmin (611.213 Pa).
        Param::PMin => return Ok(if97::get_ptrip()),
        Param::PMax => return Ok(if97::get_pmax()),
        Param::MolarMass => return Ok(if97::get_mw()),
        Param::AcentricFactor => return if97::get_acentric(),
        _ => {}
    }

    let (pair, v1, v2) = generate_update_pair(name1, value1, name2, value2)
        .ok_or_else(|| Error::Value("This pair of inputs is not yet supported".into()))?;
    match pair {
        // canonical order: P, T
        InputPair::PT => forward_tp(output, v2, v1),
        // canonical order: P, Q
        InputPair::PQ => pq(output, v1, v2),
        // canonical order: Q, T
        InputPair::QT => {
            let (q, t) = (v1, v2);
            if output == Param::T {
                return Ok(t);
            }
            pq(output, if97::psat97(t)?, q)
        }
        // canonical order: Hmass, P
        InputPair::HmassP => hp(output, v2, v1),
        // canonical order: P, Smass
        InputPair::PSmass => ps(output, v1, v2),
        // canonical order: Hmass, Smass
        InputPair::HmassSmass => hs(output, v1, v2),
        _ => Err(Error::NotImplemented(format!(
            "Input pair {pair:?} not implemented for IF97"
        ))),
    }
}

fn forward_tp(output: Param, t: f64, p: f64) -> Result<f64> {
    match output {
        Param::T => Ok(t),
        Param::P => Ok(p),
        Param::Dmass => if97::rhomass_tp(t, p),
        Param::Hmass => if97::hmass_tp(t, p),
        Param::Smass => if97::smass_tp(t, p),
        Param::Umass => if97::umass_tp(t, p),
        Param::Cpmass => if97::cpmass_tp(t, p),
        Param::Cvmass => if97::cvmass_tp(t, p),
        Param::SpeedSound => if97::speed_sound_tp(t, p),
        Param::Viscosity => if97::visc_tp(t, p),
        Param::Conductivity => if97::tcond_tp(t, p),
        Param::Prandtl => if97::prandtl_tp(t, p),
        Param::Q => Err(Error::Input("Can't determine Q from T & P".into())),
        _ => Err(Error::NotImplemented(format!(
            "Output {output:?} not implemented for IF97 (T,p)"
        ))),
    }
}

fn pq(output: Param, p: f64, q: f64) -> Result<f64> {
    match output {
        Param::T => if97::tsat97(p),
        Param::P => Ok(p),
        Param::Q => Ok(q),
        Param::Dmass => if97::rhomass_pq(p, q),
        Param::Hmass => if97::hmass_pq(p, q),
        Param::Smass => if97::smass_pq(p, q),
        Param::Umass => if97::umass_pq(p, q),
        Param::SurfaceTension => if97::sigma97(if97::tsat97(p)?),
        Param::Viscosity if q == 0.0 => if97::viscliq_p(p),
        Param::Viscosity if q == 1.0 => if97::viscvap_p(p),
        Param::Conductivity if q == 0.0 => if97::tcondliq_p(p),
        Param::Conductivity if q == 1.0 => if97::tcondvap_p(p),
        Param::Cpmass if q == 0.0 => if97::cpliq_p(p),
        Param::Cpmass if q == 1.0 => if97::cpvap_p(p),
        Param::Cvmass if q == 0.0 => if97::cvliq_p(p),
        Param::Cvmass if q == 1.0 => if97::cvvap_p(p),
        Param::SpeedSound if q == 0.0 => if97::speed_soundliq_p(p),
        Param::SpeedSound if q == 1.0 => if97::speed_soundvap_p(p),
        _ => Err(Error::NotImplemented(format!(
            "Output {output:?} not implemented for IF97 (p,Q)"
        ))),
    }
}

fn hp(output: Param, p: f64, h: f64) -> Result<f64> {
    match output {
        Param::T => if97::t_phmass(p, h),
        Param::P => Ok(p),
        Param::Hmass => Ok(h),
        // Backend-faithful density (golden-verified): in the dome, inverse
        // phase-weighted mixing of the saturation-curve densities; outside,
        // forward density at the clipped backward temperature. Y_pX's
        // backward-T round trip differs by ~1e-5 here.
        Param::Dmass => {
            if if97::region_ph(p, h)? == 4 {
                let q = 1.0_f64.min(
                    0.0_f64.max((h - if97::hliq_p(p)?) / (if97::hvap_p(p)? - if97::hliq_p(p)?)),
                );
                Ok(1.0 / (q / if97::rhovap_p(p)? + (1.0 - q) / if97::rholiq_p(p)?))
            } else {
                if97::rhomass_tp(if97::t_phmass(p, h)?, p)
            }
        }
        Param::Smass => if97::smass_phmass(p, h),
        // CoolProp convention (golden-verified): Q = -1 outside the dome,
        // where raw IF97 Q_pX would clamp to 0/1.
        Param::Q => {
            if if97::region_ph(p, h)? == 4 {
                if97::q_phmass(p, h)
            } else {
                Ok(-1.0)
            }
        }
        _ => Err(Error::NotImplemented(format!(
            "Output {output:?} not implemented for IF97 (h,p)"
        ))),
    }
}

fn ps(output: Param, p: f64, s: f64) -> Result<f64> {
    match output {
        Param::T => if97::t_psmass(p, s),
        Param::P => Ok(p),
        Param::Smass => Ok(s),
        // Backend-faithful density — see the (h,p) case above.
        Param::Dmass => {
            if if97::region_ps(p, s)? == 4 {
                let q = 1.0_f64.min(
                    0.0_f64.max((s - if97::sliq_p(p)?) / (if97::svap_p(p)? - if97::sliq_p(p)?)),
                );
                Ok(1.0 / (q / if97::rhovap_p(p)? + (1.0 - q) / if97::rholiq_p(p)?))
            } else {
                if97::rhomass_tp(if97::t_psmass(p, s)?, p)
            }
        }
        Param::Hmass => if97::hmass_psmass(p, s),
        // CoolProp convention (golden-verified): Q = -1 outside the dome.
        Param::Q => {
            if if97::region_ps(p, s)? == 4 {
                if97::q_psmass(p, s)
            } else {
                Ok(-1.0)
            }
        }
        _ => Err(Error::NotImplemented(format!(
            "Output {output:?} not implemented for IF97 (p,s)"
        ))),
    }
}

fn hs(output: Param, h: f64, s: f64) -> Result<f64> {
    match output {
        Param::P => if97::p_hsmass(h, s),
        Param::T => if97::t_hsmass(h, s),
        Param::Hmass => Ok(h),
        Param::Smass => Ok(s),
        _ => Err(Error::NotImplemented(format!(
            "Output {output:?} not implemented for IF97 (h,s)"
        ))),
    }
}
