//! IF97 engine — IAPWS-IF97 industrial water/steam formulation.
//!
//! Pure-Rust port of the header-only CoolProp/IF97 library at commit
//! `7aaced02` (the exact revision CoolProp v8.0.0 pins via CPM; see PLAN.md
//! 2.1), in CoolProp's configuration: SI units (Pa, J), and the direct
//! (non-iterated) SR5-05 region-3 density — neither `IAPWS_UNITS` nor
//! `REGION3_ITERATE` is defined when CoolProp includes the header.
//!
//! Coefficient tables live in [`tables`] (generated verbatim by
//! `tools/if97-extract/`); the evaluation logic is a hand port that keeps
//! upstream's floating-point operation order.
//!
//! Public function names mirror upstream (`rhomass_Tp` → [`rhomass_tp`]).
//! Errors map upstream exception conditions: `std::out_of_range` →
//! [`Error::OutOfRange`], `std::invalid_argument` → [`Error::Input`],
//! convergence failure → [`Error::Solution`].

// Ported comparisons keep upstream's literal form rather than range syntax.
#![allow(clippy::manual_range_contains)]

use rustprop_core::Error;
pub use rustprop_core::Result;

mod backward;
mod gibbs;
pub mod region3;
mod region4;
pub mod tables;
mod transport;

// ---------------------------------------------------------------------------
// Constants (upstream IF97.h, SI configuration: p_fact = 1e6, R_fact = 1000).
// Multiplications are kept in source form so every value is bit-identical to
// the C++ compile-time results.
// ---------------------------------------------------------------------------

pub(crate) const P_FACT: f64 = 1e6;
pub(crate) const R_FACT: f64 = 1000.0;

pub const TCRIT: f64 = 647.096; // K
pub const PCRIT: f64 = 22.064 * P_FACT; // Pa
pub const RHOCRIT: f64 = 322.0; // kg/m^3
pub const SCRIT: f64 = 4.41202148223476 * R_FACT; // J/kg-K
pub const TTRIP: f64 = 273.16; // K
pub const PTRIP: f64 = 0.000611657 * P_FACT; // Pa
pub const TMIN: f64 = 273.15; // K
pub const TMAX: f64 = 1073.15; // K
pub const PMIN: f64 = 0.000611213 * P_FACT; // Pa
pub const PMAX: f64 = 100.0 * P_FACT; // Pa
pub const RGAS: f64 = 0.461526 * R_FACT; // J/kg-K (mass based!)
pub const MW: f64 = 0.018015268; // kg/mol
// Bounds for region determination
pub(crate) const TEXT: f64 = 2273.15; // Region 5 temperature limit [K]
pub(crate) const PEXT: f64 = 50.0 * P_FACT; // Region 5 pressure limit [Pa]
pub(crate) const P23MIN: f64 = 16.529164252605 * P_FACT; // min p on B23 curve [Pa]
pub(crate) const T23MIN: f64 = 623.15; // min T on B23 curve [K]
pub(crate) const P2AMAX: f64 = 4.0 * P_FACT; // max p on H2a2b boundary [Pa]
pub(crate) const P2BCMIN: f64 = 6.54670 * P_FACT; // min p on H2b2c boundary [Pa]
pub(crate) const S2BC: f64 = 5.85 * R_FACT;
// Bounds for backward p(h,s) / T(h,s) determination
pub(crate) const SMIN: f64 = 0.0;
pub(crate) const SMAX: f64 = 11.921054825051103 * R_FACT;
pub(crate) const STPMAX: f64 = 6.04048367171238 * R_FACT; // s(Tmax, Pmax)
pub(crate) const SGTRIP: f64 = 9.155492076509681 * R_FACT; // sat. vapor s at triple point
pub(crate) const SFTRIP: f64 = -4.09187776773977E-7 * R_FACT; // sat. liquid s at triple point
pub(crate) const HGTRIP: f64 = 2500.9109532932 * R_FACT; // sat. vapor h at triple point
pub(crate) const HFTRIP: f64 = 5.16837786577998E-4 * R_FACT; // sat. liquid h at triple point
pub(crate) const SFT23: f64 = 3.778281340 * R_FACT; // sat. liquid s at T23min
pub(crate) const SGT23: f64 = 5.210887825 * R_FACT; // sat. vapor s at T23min
pub(crate) const S13MIN: f64 = 3.397782955 * R_FACT; // s at (T13, Pmax)
pub(crate) const S23MIN: f64 = 5.048096828 * R_FACT; // B23 bounding box
pub(crate) const S23MAX: f64 = 5.260578707 * R_FACT; // B23 bounding box
pub(crate) const H23MIN: f64 = 2.563592004E3 * R_FACT; // B23 bounding box
pub(crate) const H23MAX: f64 = 2.812942061E3 * R_FACT; // B23 bounding box

/// Fast integer power (upstream `powi`) — kept operation-for-operation so
/// rounding matches the C++ implementation.
pub(crate) fn powi(mut x: f64, mut i: i32) -> f64 {
    let mut ans = 1.0;
    if i < 0 {
        x = 1.0 / x;
        i = -i;
    }
    while i > 0 {
        if i & 1 == 1 {
            ans *= x;
        }
        x *= x;
        i >>= 1;
    }
    ans
}

/// Property key, mirroring upstream `IF97parameters`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Prop {
    Dmass,
    Hmass,
    T,
    P,
    Smass,
    Umass,
    Cpmass,
    Cvmass,
    W,
    DrhoDp,
    Mu,
    K,
    Q,
}

/// Saturated-state request, mirroring upstream `IF97SatState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SatState {
    None,
    Liquid,
    Vapor,
}

/// IF97 region, mirroring upstream `IF97REGIONS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    R1,
    R2,
    R3,
    R4,
    R5,
}

// ---------------------------------------------------------------------------
// Region 2/3 boundary (upstream `Region23_T` / `Region23_p`)
// ---------------------------------------------------------------------------

/// B23 boundary pressure [Pa] at T [K] (upstream `Region23_T`).
pub(crate) fn b23_p_from_t(t: f64) -> f64 {
    let (p_star, t_star) = (1.0 * P_FACT, 1.0);
    let theta = t / t_star;
    let n = tables::REGION23_N;
    let pi = n[0] + n[1] * theta + n[2] * theta * theta;
    pi * p_star
}

/// B23 boundary temperature [K] at p [Pa] (upstream `Region23_p`).
pub(crate) fn b23_t_from_p(p: f64) -> f64 {
    let (p_star, t_star) = (1.0 * P_FACT, 1.0);
    let pi = p / p_star;
    let n = tables::REGION23_N;
    let theta = n[3] + ((pi - n[4]) / n[2]).sqrt();
    theta * t_star
}

// ---------------------------------------------------------------------------
// Region determination and dispatch (upstream `RegionDetermination_TP`,
// `RegionOutput`)
// ---------------------------------------------------------------------------

/// Determine the IF97 region for a (T, p) state.
pub fn region_determination_tp(t: f64, p: f64) -> Result<Region> {
    if p < PMIN || p > PMAX {
        return Err(Error::OutOfRange("Pressure out of range".into()));
    }
    if t > TEXT {
        Err(Error::OutOfRange("Temperature out of range".into()))
    } else if t > TMAX && t <= TEXT {
        if p <= PEXT {
            Ok(Region::R5)
        } else {
            Err(Error::OutOfRange("Pressure out of range".into()))
        }
    } else if t > T23MIN && t <= TMAX {
        if p < P23MIN {
            Ok(Region::R2)
        } else if p > b23_p_from_t(t) {
            Ok(Region::R3)
        } else {
            Ok(Region::R2)
        }
    } else if t >= TMIN && t <= T23MIN {
        let psat = region4::p_t(t)?;
        if p > psat {
            Ok(Region::R1)
        } else if p < psat {
            Ok(Region::R2)
        } else {
            Ok(Region::R4)
        }
    } else {
        Err(Error::OutOfRange("Temperature out of range".into()))
    }
}

/// Region-routed property evaluation (upstream `RegionOutput`).
pub(crate) fn region_output(outkey: Prop, t: f64, p: f64, state: SatState) -> Result<f64> {
    match region_determination_tp(t, p)? {
        Region::R1 => {
            if state == SatState::Vapor {
                gibbs::REGION2.output(outkey, t, p) // on sat. curve, vapor side
            } else {
                gibbs::REGION1.output(outkey, t, p)
            }
        }
        Region::R2 => {
            if state == SatState::Liquid {
                gibbs::REGION1.output(outkey, t, p) // on sat. curve, liquid side
            } else {
                gibbs::REGION2.output(outkey, t, p)
            }
        }
        Region::R3 => region3::output(outkey, t, p, state),
        Region::R4 => match state {
            SatState::Vapor => gibbs::REGION2.output(outkey, t, p),
            SatState::Liquid => gibbs::REGION1.output(outkey, t, p),
            SatState::None => Err(Error::OutOfRange(
                "Cannot use Region 4 with T and p as inputs".into(),
            )),
        },
        Region::R5 => gibbs::REGION5.output(outkey, t, p),
    }
}

// ---------------------------------------------------------------------------
// API — forward properties from (T, p)
// ---------------------------------------------------------------------------

/// Mass density [kg/m^3] from T [K] and p [Pa].
pub fn rhomass_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Dmass, t, p, SatState::None)
}
/// Mass enthalpy [J/kg] from T [K] and p [Pa].
pub fn hmass_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Hmass, t, p, SatState::None)
}
/// Mass entropy [J/kg/K] from T [K] and p [Pa].
pub fn smass_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Smass, t, p, SatState::None)
}
/// Mass internal energy [J/kg] from T [K] and p [Pa].
pub fn umass_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Umass, t, p, SatState::None)
}
/// Mass isobaric specific heat [J/kg/K] from T [K] and p [Pa].
pub fn cpmass_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Cpmass, t, p, SatState::None)
}
/// Mass isochoric specific heat [J/kg/K] from T [K] and p [Pa].
pub fn cvmass_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Cvmass, t, p, SatState::None)
}
/// Speed of sound [m/s] from T [K] and p [Pa].
pub fn speed_sound_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::W, t, p, SatState::None)
}
/// (d rho / d p) at constant T [kg/m^3/Pa] from T [K] and p [Pa].
pub fn drhodp_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::DrhoDp, t, p, SatState::None)
}

// ---------------------------------------------------------------------------
// API — transport properties
// ---------------------------------------------------------------------------

/// Viscosity [Pa-s] from T [K] and rho [kg/m^3] (region-free).
pub fn visc_trho(t: f64, rho: f64) -> f64 {
    transport::visc(t, rho)
}
/// Viscosity [Pa-s] from T [K] and p [Pa].
pub fn visc_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::Mu, t, p, SatState::None)
}
/// Thermal conductivity [W/m/K] from T [K] and p [Pa].
pub fn tcond_tp(t: f64, p: f64) -> Result<f64> {
    region_output(Prop::K, t, p, SatState::None)
}
/// Prandtl number [-] from T [K] and p [Pa].
pub fn prandtl_tp(t: f64, p: f64) -> Result<f64> {
    Ok(visc_tp(t, p)? * cpmass_tp(t, p)? * (1000.0 / R_FACT) / tcond_tp(t, p)?)
}

// ---------------------------------------------------------------------------
// API — saturated liquid/vapor properties from p
// ---------------------------------------------------------------------------

macro_rules! sat_fn {
    ($liq:ident, $vap:ident, $key:expr, $doc:literal) => {
        #[doc = concat!("Saturated liquid ", $doc, " from p [Pa].")]
        pub fn $liq(p: f64) -> Result<f64> {
            region_output($key, tsat97(p)?, p, SatState::Liquid)
        }
        #[doc = concat!("Saturated vapor ", $doc, " from p [Pa].")]
        pub fn $vap(p: f64) -> Result<f64> {
            region_output($key, tsat97(p)?, p, SatState::Vapor)
        }
    };
}

sat_fn!(rholiq_p, rhovap_p, Prop::Dmass, "mass density [kg/m^3]");
sat_fn!(hliq_p, hvap_p, Prop::Hmass, "mass enthalpy [J/kg]");
sat_fn!(sliq_p, svap_p, Prop::Smass, "mass entropy [J/kg/K]");
sat_fn!(uliq_p, uvap_p, Prop::Umass, "mass internal energy [J/kg]");
sat_fn!(
    cpliq_p,
    cpvap_p,
    Prop::Cpmass,
    "isobaric specific heat [J/kg/K]"
);
sat_fn!(
    cvliq_p,
    cvvap_p,
    Prop::Cvmass,
    "isochoric specific heat [J/kg/K]"
);
sat_fn!(
    speed_soundliq_p,
    speed_soundvap_p,
    Prop::W,
    "speed of sound [m/s]"
);
sat_fn!(viscliq_p, viscvap_p, Prop::Mu, "viscosity [Pa-s]");
sat_fn!(
    tcondliq_p,
    tcondvap_p,
    Prop::K,
    "thermal conductivity [W/m/K]"
);

/// Saturated liquid Prandtl number [-] from p [Pa].
pub fn prandtlliq_p(p: f64) -> Result<f64> {
    Ok(viscliq_p(p)? * cpliq_p(p)? * (1000.0 / R_FACT) / tcondliq_p(p)?)
}
/// Saturated vapor Prandtl number [-] from p [Pa].
pub fn prandtlvap_p(p: f64) -> Result<f64> {
    Ok(viscvap_p(p)? * cpvap_p(p)? * (1000.0 / R_FACT) / tcondvap_p(p)?)
}

// ---------------------------------------------------------------------------
// API — 2-phase and saturation-curve functions
// ---------------------------------------------------------------------------

/// Saturation temperature [K] from p [Pa].
pub fn tsat97(p: f64) -> Result<f64> {
    region4::t_p(p)
}
/// Saturation pressure [Pa] from T [K].
pub fn psat97(t: f64) -> Result<f64> {
    region4::p_t(t)
}
/// Surface tension [N/m] from T [K].
pub fn sigma97(t: f64) -> Result<f64> {
    region4::sigma_t(t)
}

// ---------------------------------------------------------------------------
// API — backward functions
// ---------------------------------------------------------------------------

/// T [K] from p [Pa] and mass enthalpy [J/kg].
pub fn t_phmass(p: f64, h: f64) -> Result<f64> {
    backward::region_output_backward(p, h, Prop::Hmass, true, SatState::None)
}
/// Mass density [kg/m^3] from p [Pa] and mass enthalpy [J/kg].
pub fn rhomass_phmass(p: f64, h: f64) -> Result<f64> {
    backward::y_px(Prop::Dmass, p, h, Prop::Hmass)
}
/// T [K] from p [Pa] and mass entropy [J/kg/K].
pub fn t_psmass(p: f64, s: f64) -> Result<f64> {
    backward::region_output_backward(p, s, Prop::Smass, true, SatState::None)
}
/// Mass density [kg/m^3] from p [Pa] and mass entropy [J/kg/K].
pub fn rhomass_psmass(p: f64, s: f64) -> Result<f64> {
    backward::y_px(Prop::Dmass, p, s, Prop::Smass)
}
/// p [Pa] from mass enthalpy [J/kg] and mass entropy [J/kg/K].
pub fn p_hsmass(h: f64, s: f64) -> Result<f64> {
    backward::backward_output_hs(Prop::P, h, s)
}
/// T [K] from mass enthalpy [J/kg] and mass entropy [J/kg/K].
pub fn t_hsmass(h: f64, s: f64) -> Result<f64> {
    backward::backward_output_hs(Prop::T, h, s)
}
/// Mass enthalpy [J/kg] from p [Pa] and mass entropy [J/kg/K].
pub fn hmass_psmass(p: f64, s: f64) -> Result<f64> {
    backward::y_px(Prop::Hmass, p, s, Prop::Smass)
}
/// Mass entropy [J/kg/K] from p [Pa] and mass enthalpy [J/kg].
pub fn smass_phmass(p: f64, h: f64) -> Result<f64> {
    backward::y_px(Prop::Smass, p, h, Prop::Hmass)
}
/// IF97 region (1..4) for a (p, h) state.
pub fn region_ph(p: f64, h: f64) -> Result<i32> {
    backward::backward_region(p, h, Prop::Hmass)
}
/// IF97 region (1..4) for a (p, s) state.
pub fn region_ps(p: f64, s: f64) -> Result<i32> {
    backward::backward_region(p, s, Prop::Smass)
}

// ---------------------------------------------------------------------------
// API — quality functions
// ---------------------------------------------------------------------------

/// Mass enthalpy [J/kg] from p [Pa] and quality Q [-].
pub fn hmass_pq(p: f64, q: f64) -> Result<f64> {
    backward::x_pq(Prop::Hmass, p, q)
}
/// Mass internal energy [J/kg] from p [Pa] and quality Q [-].
pub fn umass_pq(p: f64, q: f64) -> Result<f64> {
    backward::x_pq(Prop::Umass, p, q)
}
/// Mass entropy [J/kg/K] from p [Pa] and quality Q [-].
pub fn smass_pq(p: f64, q: f64) -> Result<f64> {
    backward::x_pq(Prop::Smass, p, q)
}
/// Specific volume [m^3/kg] from p [Pa] and quality Q [-].
pub fn v_pq(p: f64, q: f64) -> Result<f64> {
    Ok(1.0 / backward::x_pq(Prop::Dmass, p, q)?)
}
/// Mass density [kg/m^3] from p [Pa] and quality Q [-].
pub fn rhomass_pq(p: f64, q: f64) -> Result<f64> {
    backward::x_pq(Prop::Dmass, p, q)
}
/// Quality [-] from p [Pa] and mass enthalpy [J/kg].
pub fn q_phmass(p: f64, h: f64) -> Result<f64> {
    backward::q_px(p, h, Prop::Hmass)
}
/// Quality [-] from p [Pa] and mass internal energy [J/kg].
pub fn q_pumass(p: f64, u: f64) -> Result<f64> {
    backward::q_px(p, u, Prop::Umass)
}
/// Quality [-] from p [Pa] and mass entropy [J/kg/K].
pub fn q_psmass(p: f64, s: f64) -> Result<f64> {
    backward::q_px(p, s, Prop::Smass)
}
/// Quality [-] from p [Pa] and mass density [kg/m^3].
pub fn q_prhomass(p: f64, rho: f64) -> Result<f64> {
    backward::q_px(p, rho, Prop::Dmass)
}
/// Quality [-] from p [Pa] and specific volume [m^3/kg].
pub fn q_pv(p: f64, v: f64) -> Result<f64> {
    backward::q_px(p, 1.0 / v, Prop::Dmass)
}

// ---------------------------------------------------------------------------
// API — trivial functions
// ---------------------------------------------------------------------------

pub fn get_ttrip() -> f64 {
    TTRIP
}
pub fn get_ptrip() -> f64 {
    PTRIP
}
pub fn get_tcrit() -> f64 {
    TCRIT
}
pub fn get_pcrit() -> f64 {
    PCRIT
}
pub fn get_rhocrit() -> f64 {
    RHOCRIT
}
pub fn get_tmin() -> f64 {
    TMIN
}
pub fn get_pmin() -> f64 {
    PMIN
}
pub fn get_tmax() -> f64 {
    TMAX
}
pub fn get_pmax() -> f64 {
    PMAX
}
pub fn get_mw() -> f64 {
    MW
}
pub fn get_rgas() -> f64 {
    RGAS
}
/// Acentric factor [-] (upstream `get_Acentric`).
pub fn get_acentric() -> Result<f64> {
    Ok(-(psat97(0.7 * TCRIT)? / PCRIT).log10() - 1.0)
}
