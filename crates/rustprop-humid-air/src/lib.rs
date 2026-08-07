//! Humid-air engine (PLAN.md Phase 9) — port of CoolProp 8's
//! `src/HumidAirProp.cpp` (`HAPropsSI` semantics) on the ASHRAE RP-1485
//! virial mixture model, consuming this port's HEOS Water and (pseudo-pure)
//! Air plus IF97's region-4 saturation curve exactly where upstream does.
//!
//! Deliberate upstream reproductions:
//! - Three DIFFERENT gas constants in different formulas (global 8.314472;
//!   f_factor / water-IG-entropy 8.314371; the Lemmon air value 8.314510) —
//!   never unified.
//! - The dry-air molar mass is the HARDCODED 0.028966 kg/mol in `M_ha` and
//!   `epsilon = 0.621945`; `MM_Air()` (0.02896546…) appears only in the
//!   transport mixing rules.
//! - Water saturation pressure comes from the HEOS backend inside
//!   `f_factor` but from IF97 region 4 everywhere else; the ice/liquid split
//!   is `T > 273.16` in `f_factor` and the wet-bulb solver but
//!   `T >= 273.16` in `MoleFractionWater`/`RelativeHumidity`/dewpoint.
//! - `f_factor` is clamped to >= 1; its secant (like `MolarVolume`'s) exits
//!   silently at the iteration cap.
//! - `DewpointTemperature`'s dry-air guard tests `(1 - psi_w) < 1e-16` —
//!   the wrong side (fires for pure water) — and its initial guess uses the
//!   TOTAL pressure in `Tsat97`.
//!
//! One deliberate deviation, logged in PLAN.md: upstream's `HAPropsSI`
//! swallows every error and returns +inf (the message parked in a global);
//! this port's `ha_props_si` returns `Result` with the same message texts.

// Upstream magic anchor constants are carried verbatim, including digits
// beyond f64 precision — bit-faithful by mandate.
#![allow(clippy::excessive_precision)]

mod ice;

pub use ice::{h_ice, psub_ice};

use rustprop_core::fluid::FluidData;
use rustprop_core::{Error, Result};
use rustprop_heos::flash_pt::PtFlash;
use rustprop_heos::flash_px::HeosState;
use rustprop_heos::solvers::{brent, secant};

const EPSILON: f64 = 0.621945;
const R_BAR: f64 = 8.314472;
const R_BAR_WS: f64 = 8.314371;
const R_BAR_LEMMON: f64 = 8.314510;
const M_AIR_HARDCODED: f64 = 0.028966;

/// The accepted HAPropsSI parameter set (upstream `givens`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum HaParam {
    HumRat,
    PsiW,
    Tdp,
    Twb,
    Enthalpy,
    EnthalpyHa,
    InternalEnergy,
    InternalEnergyHa,
    Entropy,
    EntropyHa,
    Rh,
    T,
    P,
    Vda,
    Vha,
    Visc,
    Cond,
    Cp,
    CpHa,
    Cv,
    CvHa,
    PartialPressureWater,
    IsentropicExponent,
    SpeedOfSound,
    CompressibilityFactor,
}

/// Upstream `Name2Type` — exact, case-sensitive strings.
pub fn name_to_type(name: &str) -> Result<HaParam> {
    use HaParam::*;
    Ok(match name {
        "Omega" | "HumRat" | "W" => HumRat,
        "psi_w" | "Y" => PsiW,
        "Tdp" | "T_dp" | "DewPoint" | "D" => Tdp,
        "Twb" | "T_wb" | "WetBulb" | "B" => Twb,
        "Enthalpy" | "H" | "Hda" => Enthalpy,
        "Hha" => EnthalpyHa,
        "InternalEnergy" | "U" | "Uda" => InternalEnergy,
        "Uha" => InternalEnergyHa,
        "Entropy" | "S" | "Sda" => Entropy,
        "Sha" => EntropyHa,
        "RH" | "RelHum" | "R" => Rh,
        "Tdb" | "T_db" | "T" => T,
        "P" => P,
        "V" | "Vda" => Vda,
        "Vha" => Vha,
        "mu" | "Visc" | "M" => Visc,
        "k" | "Conductivity" | "K" => Cond,
        "C" | "cp" => Cp,
        "Cha" | "cp_ha" => CpHa,
        "CV" => Cv,
        "CVha" | "cv_ha" => CvHa,
        "P_w" => PartialPressureWater,
        "isentropic_exponent" => IsentropicExponent,
        "speed_of_sound" => SpeedOfSound,
        "Z" => CompressibilityFactor,
        _ => {
            return Err(Error::Value(format!(
                "Sorry, your input [{name}] was not understood to Name2Type. Acceptable values are T,P,R,W,D,B,H,S,M,K and aliases thereof\n"
            )));
        }
    })
}

/// Upstream `check_bounds`: returns (ok, min, max).
fn bounds_of(prop: HaParam) -> (f64, f64) {
    use HaParam::*;
    match prop {
        P => (0.00001e6, 10e6),
        T | Tdp | Twb => (-143.15 + 273.15, 350.0 + 273.15),
        HumRat => (0.0, 10.0),
        PsiW => (0.0, 0.94145),
        Rh => (0.0, 1.0),
        _ => (f64::NEG_INFINITY, f64::INFINITY),
    }
}

fn check_bounds(prop: HaParam, value: f64) -> bool {
    if !value.is_finite() {
        return false;
    }
    let (lo, hi) = bounds_of(prop);
    value >= lo && value <= hi
}

/// Everything the humid-air model needs: HEOS Water, HEOS pseudo-pure Air.
pub struct HumidAir {
    water: PtFlash,
    air: PtFlash,
    mm_water: f64,
    mm_air: f64,
    ref_offsets: std::sync::OnceLock<RefOffsets>,
}

struct RefOffsets {
    t_red_w: f64,
    rho_red_w: f64,
    hoffset_w: f64,
    hoffset_a: f64,
    soffset_w: f64,
    soffset_a: f64,
    ln_delta_air_s: f64,
}

struct Virials {
    b_aa: f64,
    db_aa: f64,
    c_aaa: f64,
    dc_aaa: f64,
    b_ww: f64,
    db_ww: f64,
    c_www: f64,
    dc_www: f64,
}

impl HumidAir {
    pub fn new(water: &'static FluidData, air: &'static FluidData) -> Self {
        let water = PtFlash::new(water);
        let air = PtFlash::new(air);
        let mm_water = water.eos.molar_mass;
        let mm_air = air.eos.molar_mass;
        HumidAir {
            water,
            air,
            mm_water,
            mm_air,
            ref_offsets: std::sync::OnceLock::new(),
        }
    }

    /// EOS virials in the delta -> 0 limit (upstream `calc_all_virials`).
    fn virials(&self, t: f64) -> Virials {
        let one = |fl: &PtFlash| -> (f64, f64, f64, f64) {
            let t_red = fl.eos.t_reducing;
            let rho_red = fl.eos.rhomolar_reducing;
            let tau = t_red / t;
            let dtau_dt = -t_red / (t * t);
            let d = fl.eos.alphar_all(tau, 1e-12);
            (
                d.d10 / rho_red,
                d.d11 / rho_red * dtau_dt,
                d.d20 / (rho_red * rho_red),
                d.d21 / (rho_red * rho_red) * dtau_dt,
            )
        };
        let (b_aa, db_aa, c_aaa, dc_aaa) = one(&self.air);
        let (b_ww, db_ww, c_www, dc_www) = one(&self.water);
        Virials {
            b_aa,
            db_aa,
            c_aaa,
            dc_aaa,
            b_ww,
            db_ww,
            c_www,
            dc_www,
        }
    }

    fn b_m(&self, t: f64, psi_w: f64) -> f64 {
        let v = self.virials(t);
        (1.0 - psi_w).powi(2) * v.b_aa
            + 2.0 * (1.0 - psi_w) * psi_w * b_aw(t)
            + psi_w * psi_w * v.b_ww
    }
    fn db_m_dt(&self, t: f64, psi_w: f64) -> f64 {
        let v = self.virials(t);
        (1.0 - psi_w).powi(2) * v.db_aa
            + 2.0 * (1.0 - psi_w) * psi_w * db_aw_dt(t)
            + psi_w * psi_w * v.db_ww
    }
    fn c_m(&self, t: f64, psi_w: f64) -> f64 {
        let v = self.virials(t);
        (1.0 - psi_w).powi(3) * v.c_aaa
            + 3.0 * (1.0 - psi_w).powi(2) * psi_w * c_aaw(t)
            + 3.0 * (1.0 - psi_w) * psi_w * psi_w * c_aww(t)
            + psi_w.powi(3) * v.c_www
    }
    fn dc_m_dt(&self, t: f64, psi_w: f64) -> f64 {
        let v = self.virials(t);
        (1.0 - psi_w).powi(3) * v.dc_aaa
            + 3.0 * (1.0 - psi_w).powi(2) * psi_w * dc_aaw_dt(t)
            + 3.0 * (1.0 - psi_w) * psi_w * psi_w * dc_aww_dt(t)
            + psi_w.powi(3) * v.dc_www
    }

    /// Reference offsets (upstream `ensure_ref_offsets`), computed once from
    /// the magic anchor states.
    fn refs(&self) -> &RefOffsets {
        self.ref_offsets.get_or_init(|| {
            let t_red_w = self.water.eos.t_reducing;
            let rho_red_w = self.water.eos.rhomolar_reducing;
            let rho_red_a = self.air.eos.rhomolar_reducing;
            let a0_at = |fl: &PtFlash, tau: f64, delta: f64| fl.eos.alpha0_all(tau, delta);

            // Water enthalpy offset (R_bar 8.314472)
            let hoffset_w = {
                let (tref, vref, href) = (473.15, 0.038837428192186184, 51885.582451893446);
                let tauref = t_red_w / tref;
                let d = a0_at(&self.water, tauref, (1.0 / vref) / rho_red_w);
                href - R_BAR * tref * (1.0 + tauref * d.d01)
            };
            // Water entropy offset (R_bar_ws 8.314371)
            let soffset_w = {
                let (tref, pref, sref) = (473.15, 101325.0, 141.18297895840303);
                let tauref = t_red_w / tref;
                let rho = pref / (R_BAR_WS * tref);
                let d = a0_at(&self.water, tauref, rho / rho_red_w);
                sref - R_BAR_WS * (tauref * d.d01 - d.d00)
            };
            // Air enthalpy offset (R_bar_Lemmon 8.314510)
            let hoffset_a = {
                let (tref, vref, href) = (473.15, 0.038837428192186184, 13782.240592933371);
                let tauref = 132.6312 / tref;
                let d = a0_at(&self.air, tauref, (1.0 / vref) / rho_red_a);
                href - R_BAR_LEMMON * tref * (1.0 + tauref * d.d01)
            };
            // Air entropy offset + the constant ln(delta) at 1/vmolar_a0
            let (soffset_a, ln_delta_air_s) = {
                let (t0, p0) = (273.15, 101325.0);
                let (tref, vref, sref) = (473.15, 0.038837605637863169, 212.22365283759311);
                let vmolar_a_0 = R_BAR_LEMMON * t0 / p0;
                let tauref = 132.6312 / tref;
                let d = a0_at(&self.air, tauref, (1.0 / vmolar_a_0) / rho_red_a);
                let sref_eos = R_BAR_LEMMON * (tauref * d.d01 - d.d00)
                    + R_BAR_LEMMON * (vref / vmolar_a_0).ln();
                (sref - sref_eos, (1.0 / (vmolar_a_0 * rho_red_a)).ln())
            };
            RefOffsets {
                t_red_w,
                rho_red_w,
                hoffset_w,
                hoffset_a,
                soffset_w,
                soffset_a,
                ln_delta_air_s,
            }
        })
    }

    /// alpha0 (value at delta = 1, tau-derivative) for both fluids.
    fn alpha0_pair(&self, t: f64) -> (f64, f64, f64, f64) {
        let dw = self
            .water
            .eos
            .alpha0_all(self.water.eos.t_reducing / t, 1.0);
        let da = self.air.eos.alpha0_all(132.6312 / t, 1.0);
        (da.d00, da.d01, dw.d00, dw.d01)
    }

    fn ideal_gas_h_water(&self, t: f64) -> f64 {
        let r = self.refs();
        let (_, _, _, da0_w) = self.alpha0_pair(t);
        let tau = r.t_red_w / t;
        -0.01102303806 + r.hoffset_w + R_BAR * t * (1.0 + tau * da0_w)
    }
    fn ideal_gas_s_water(&self, t: f64, p: f64) -> f64 {
        let r = self.refs();
        let (_, _, a0_w, da0_w) = self.alpha0_pair(t);
        let tau = r.t_red_w / t;
        let ln_delta = (p / (R_BAR_WS * t * r.rho_red_w)).ln();
        r.soffset_w + R_BAR_WS * (tau * da0_w - a0_w - ln_delta)
    }
    fn ideal_gas_h_air(&self, t: f64) -> f64 {
        let r = self.refs();
        let (_, da0_a, _, _) = self.alpha0_pair(t);
        let tau = 132.6312 / t;
        -7914.149298 + r.hoffset_a + R_BAR_LEMMON * t * (1.0 + tau * da0_a)
    }
    fn ideal_gas_s_air(&self, t: f64, vmolar_a: f64) -> f64 {
        let r = self.refs();
        let (a0_a, da0_a, _, _) = self.alpha0_pair(t);
        let (t0, p0) = (273.15, 101325.0);
        let vmolar_a_0 = R_BAR_LEMMON * t0 / p0;
        -196.1375815
            + r.soffset_a
            + R_BAR_LEMMON * (132.6312 / t * da0_a - a0_a - r.ln_delta_air_s)
            + R_BAR_LEMMON * (vmolar_a / vmolar_a_0).ln()
    }

    /// Upstream `HenryConstant` [1/Pa] (N2/O2/Ar at fixed mole fractions).
    fn henry_constant(&self, t: f64) -> Result<f64> {
        let tc = 647.096;
        let tr = t / tc;
        let tau = 1.0 - tr;
        let p_ws = rustprop_if97::psat97(t)?;
        let beta = |c0: f64, c1: f64, c2: f64| {
            p_ws * (c0 / tr + c1 * tau.powf(0.355) / tr + c2 * tr.powf(-0.41) * tau.exp()).exp()
        };
        let beta_n2 = beta(-9.67578, 4.72162, 11.70585);
        let beta_o2 = beta(-9.44833, 4.43822, 11.42005);
        let beta_ar = beta(-8.40954, 4.29587, 10.52779);
        let beta_a = 1.0 / (0.7812 / beta_n2 + 0.2095 / beta_o2 + 0.0093 / beta_ar);
        Ok(1.0 / (1.01325 * beta_a))
    }

    /// Upstream `isothermal_compressibility` [1/Pa].
    fn isothermal_compressibility(&self, t: f64, p: f64) -> Result<f64> {
        if t > 273.16 {
            let p_ws_t = rustprop_if97::psat97(t)?;
            if (p - p_ws_t).abs() <= p_ws_t * 3.3e-5 {
                // The saturation-adjacent polynomial (upstream's correlation
                // branch, forced near the IF97 curve).
                let tc = t - 273.15;
                return Ok((50.88496
                    + 0.6163813 * tc
                    + 1.459187e-3 * tc * tc
                    + 20.08438e-6 * tc.powi(3)
                    - 58.47727e-9 * tc.powi(4)
                    + 410.4110e-12 * tc.powi(5))
                    / (1.0 + 19.67348e-3 * tc)
                    * 1e-11);
            }
            let rho_mass = rustprop_if97::rhomass_tp(t, p)?;
            let rho_molar = rho_mass / self.mm_water;
            // k_T = 1/(rho * dp/drho|T) from the Water EOS.
            let tau = self.water.eos.t_reducing / t;
            let delta = rho_molar / self.water.eos.rhomolar_reducing;
            let d = self.water.eos.alphar_all(tau, delta);
            let dpdrho = self.water.eos.gas_constant
                * t
                * (1.0 + 2.0 * delta * d.d10 + delta * delta * d.d20);
            Ok(1.0 / (rho_molar * dpdrho))
        } else {
            Ok(ice::isotherm_compress_ice(t, p))
        }
    }

    /// Upstream `f_factor` — the enhancement factor, secant on ln f with
    /// the clamp at unity from below; NO convergence guard (returns the
    /// 100th iterate as upstream does).
    fn f_factor(&self, t: f64, p: f64) -> Result<f64> {
        let (p_ws, vbar_ws, mut beta_h);
        if t > 273.16 {
            let sat = self.water.sat().qt_flash(t, 0.0)?;
            p_ws = sat.p;
            vbar_ws = 1.0 / sat.rho_l;
            beta_h = self.henry_constant(t)?;
        } else {
            p_ws = ice::psub_ice(t);
            beta_h = 0.0;
            vbar_ws = ice::dg_dp_ice(t, p) * self.mm_water;
        }
        let mut k_t = self.isothermal_compressibility(t, p)?;
        if p_ws > p {
            k_t = 0.0;
            beta_h = 0.0;
        }
        let v = self.virials(t);
        let (b_aa, b_ww) = (v.b_aa, v.b_ww);
        let baw = b_aw(t);
        let (c_aaa, c_www) = (v.c_aaa, v.c_www);
        let caaw = c_aaw(t);
        let caww = c_aww(t);
        let rbar = R_BAR_WS; // f_factor's own 8.314371

        let mut x1 = 0.0f64;
        let mut x2 = 0.0f64;
        let mut y1 = 0.0;
        let mut f = 0.0f64;
        let mut change = f64::INFINITY;
        let mut iter = 1;
        while (iter <= 3 || change > 1e-8) && iter < 100 {
            if iter == 1 {
                x1 = 1.00;
                f = x1;
            }
            if iter == 2 {
                x2 = 1.00 + 0.000001;
                f = x2;
            }
            if iter > 2 {
                f = x2;
            }
            let lhs = f.ln();
            let psi_ws = f * p_ws / p;
            let line1 = ((1.0 + k_t * p_ws) * (p - p_ws) - k_t * (p * p - p_ws * p_ws) / 2.0)
                / (rbar * t)
                * vbar_ws
                + (1.0 - beta_h * (1.0 - psi_ws) * p).ln();
            let line2 = (1.0 - psi_ws).powi(2) * p / (rbar * t) * b_aa
                - 2.0 * (1.0 - psi_ws).powi(2) * p / (rbar * t) * baw
                - (p - p_ws - (1.0 - psi_ws).powi(2) * p) / (rbar * t) * b_ww;
            let line3 = (1.0 - psi_ws).powi(3) * p * p / (rbar * t).powi(2) * c_aaa
                + (3.0 * (1.0 - psi_ws).powi(2) * (1.0 - 2.0 * (1.0 - psi_ws)) * p * p)
                    / (2.0 * (rbar * t).powi(2))
                    * caaw;
            let line4 = -3.0 * (1.0 - psi_ws).powi(2) * psi_ws * p * p / (rbar * t).powi(2) * caww
                - ((3.0 - 2.0 * psi_ws) * psi_ws * psi_ws * p * p - p_ws * p_ws)
                    / (2.0 * (rbar * t).powi(2))
                    * c_www;
            let line5 = -((1.0 - psi_ws).powi(2) * (-2.0 + 3.0 * psi_ws) * psi_ws * p * p)
                / (rbar * t).powi(2)
                * b_aa
                * b_ww;
            let line6 = -(2.0 * (1.0 - psi_ws).powi(3) * (-1.0 + 3.0 * psi_ws) * p * p)
                / (rbar * t).powi(2)
                * b_aa
                * baw;
            let line7 = (6.0 * (1.0 - psi_ws).powi(2) * psi_ws * psi_ws * p * p)
                / (rbar * t).powi(2)
                * b_ww
                * baw
                - (3.0 * (1.0 - psi_ws).powi(4) * p * p) / (2.0 * (rbar * t).powi(2)) * b_aa * b_aa;
            let line8 = -(2.0 * (1.0 - psi_ws).powi(2) * psi_ws * (-2.0 + 3.0 * psi_ws) * p * p)
                / (rbar * t).powi(2)
                * baw
                * baw
                - (p_ws * p_ws - (4.0 - 3.0 * psi_ws) * psi_w_pow3(psi_ws) * p * p)
                    / (2.0 * (rbar * t).powi(2))
                    * b_ww
                    * b_ww;
            let resid = lhs - (line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8);
            if iter == 1 {
                y1 = resid;
            }
            if iter > 1 {
                let y2 = resid;
                let x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
                change = (y2 / (y2 - y1) * (x2 - x1)).abs();
                y1 = y2;
                x1 = x2;
                x2 = x3;
            }
            iter += 1;
        }
        Ok(if f >= 1.0 { f } else { 1.0 })
    }

    /// Upstream `MoleFractionWater`.
    fn mole_fraction_water(&self, t: f64, p: f64, given: HaParam, value: f64) -> Result<f64> {
        match given {
            HaParam::HumRat => Ok(value / (EPSILON + value)),
            HaParam::Rh => {
                let p_ws = if t >= 273.16 {
                    rustprop_if97::psat97(t)?
                } else {
                    ice::psub_ice(t)
                };
                let f = self.f_factor(t, p)?;
                Ok(value * f * p_ws / p)
            }
            HaParam::Tdp => {
                let tdp = value;
                let p_ws_dp = if tdp >= 273.16 {
                    rustprop_if97::psat97(tdp)?
                } else {
                    ice::psub_ice(tdp)
                };
                let f_dp = self.f_factor(tdp, p)?;
                Ok(f_dp * p_ws_dp / p)
            }
            _ => Ok(f64::NEG_INFINITY),
        }
    }

    fn relative_humidity(&self, t: f64, p: f64, psi_w: f64) -> Result<f64> {
        let p_ws = if t >= 273.16 {
            rustprop_if97::psat97(t)?
        } else {
            ice::psub_ice(t)
        };
        Ok(psi_w * p / (self.f_factor(t, p)? * p_ws))
    }

    /// Molar volume [m^3/mol_ha] — upstream's unguarded relative secant.
    fn molar_volume(&self, t: f64, p: f64, psi_w: f64) -> f64 {
        let bm = self.b_m(t, psi_w);
        let cm = self.c_m(t, psi_w);
        let f = |v: f64| (p - R_BAR * t / v * (1.0 + bm / v + cm / (v * v))) / p;
        let mut x1 = 0.0;
        let mut x2 = 0.0;
        let mut y1 = 0.0;
        let mut v_bar = 0.0;
        let mut resid = 999.0f64;
        let mut iter = 1;
        while (iter <= 3 || resid.abs() > 1e-11) && iter < 100 {
            if iter == 1 {
                x1 = R_BAR * t / p;
                v_bar = x1;
            }
            if iter == 2 {
                x2 = R_BAR * t / p + 1e-6;
                v_bar = x2;
            }
            if iter > 2 {
                v_bar = x2;
            }
            resid = f(v_bar);
            if iter == 1 {
                y1 = resid;
            }
            if iter > 1 {
                let y2 = resid;
                let x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
                y1 = y2;
                x1 = x2;
                x2 = x3;
            }
            iter += 1;
        }
        v_bar
    }

    fn pressure(&self, t: f64, v_bar: f64, psi_w: f64) -> f64 {
        let bm = self.b_m(t, psi_w);
        let cm = self.c_m(t, psi_w);
        R_BAR * t / v_bar * (1.0 + bm / v_bar + cm / (v_bar * v_bar))
    }

    /// Molar enthalpy [J/mol_ha].
    fn molar_enthalpy(&self, t: f64, _p: f64, psi_w: f64, vmolar: f64) -> f64 {
        let hbar_w = self.ideal_gas_h_water(t);
        let hbar_a = self.ideal_gas_h_air(t);
        (1.0 - psi_w) * hbar_a
            + psi_w * hbar_w
            + R_BAR
                * t
                * ((self.b_m(t, psi_w) - t * self.db_m_dt(t, psi_w)) / vmolar
                    + (self.c_m(t, psi_w) - t / 2.0 * self.dc_m_dt(t, psi_w)) / (vmolar * vmolar))
    }

    fn molar_internal_energy(&self, t: f64, p: f64, psi_w: f64, vmolar: f64) -> f64 {
        self.molar_enthalpy(t, p, psi_w, vmolar) - p * vmolar
    }

    /// Molar entropy [J/mol_ha/K] — including upstream's absolute-tolerance
    /// vbar_a secant (which regularly runs to its 100-iteration cap).
    fn molar_entropy(&self, t: f64, p: f64, psi_w: f64, v_bar: f64) -> f64 {
        let sbar_0 = 0.02366427495;
        let v = self.virials(t);
        let baa = v.b_aa;
        let caaa = v.c_aaa;
        let b = self.b_m(t, psi_w);
        let dbdt = self.db_m_dt(t, psi_w);
        let c = self.c_m(t, psi_w);
        let dcdt = self.dc_m_dt(t, psi_w);
        // vbar_a: dry-air molar volume at (T, p) from the dry virial EOS.
        let f = |va: f64| R_BAR_LEMMON * t / va * (1.0 + baa / va + caaa / (va * va)) - p;
        let mut x1 = 0.0;
        let mut x2 = 0.0;
        let mut y1 = 0.0;
        let mut vbar_a = 0.0;
        let mut resid = 999.0f64;
        let mut iter = 1;
        while (iter <= 3 || resid.abs() > 1e-8) && iter < 100 {
            if iter == 1 {
                x1 = R_BAR_LEMMON * t / p;
                vbar_a = x1;
            }
            if iter == 2 {
                x2 = R_BAR_LEMMON * t / p + 0.001;
                vbar_a = x2;
            }
            if iter > 2 {
                vbar_a = x2;
            }
            resid = f(vbar_a);
            if iter == 1 {
                y1 = resid;
            }
            if iter > 1 {
                let y2 = resid;
                let x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
                y1 = y2;
                x1 = x2;
                x2 = x3;
            }
            iter += 1;
        }
        let sbar_w = self.ideal_gas_s_water(t, p);
        let sbar_a = self.ideal_gas_s_air(t, vbar_a);
        if psi_w != 0.0 {
            sbar_0 + (1.0 - psi_w) * sbar_a + psi_w * sbar_w
                - R_BAR
                    * ((b + t * dbdt) / v_bar
                        + (c + t * dcdt) / (2.0 * v_bar * v_bar)
                        + (1.0 - psi_w) * (1.0 - psi_w).ln()
                        + psi_w * psi_w.ln())
        } else {
            sbar_0 + sbar_a
                - R_BAR * ((b + t * dbdt) / v_bar + (c + t * dcdt) / (2.0 * v_bar * v_bar))
        }
    }

    fn m_ha(&self, psi_w: f64) -> f64 {
        self.mm_water * psi_w + (1.0 - psi_w) * M_AIR_HARDCODED
    }

    /// Upstream `DewpointTemperature` — verbatim quirks included.
    fn dewpoint_temperature(&self, p: f64, psi_w: f64) -> Result<f64> {
        if (1.0 - psi_w) < 1e-16 {
            return Ok(-1.0);
        }
        let p_w = psi_w * p;
        let mut x1 = if p_w > 611.6547241637944 {
            rustprop_if97::tsat97(p)? - 1.0
        } else {
            268.0
        };
        let mut x2 = x1 + 0.1;
        let resid = |tdp: f64| -> Result<f64> {
            let p_ws_dp = if tdp >= 273.16 {
                rustprop_if97::psat97(tdp)?
            } else {
                ice::psub_ice(tdp)
            };
            Ok(p_w - self.f_factor(tdp, p)? * p_ws_dp)
        };
        let x0 = x1;
        let mut y1 = 0.0;
        let mut tdp = 0.0;
        let mut r = 999.0f64;
        let mut iter = 1;
        while (iter <= 3 || r.abs() > 1e-5) && iter < 100 {
            if iter == 1 {
                x1 = x0;
                tdp = x1;
            }
            if iter == 2 {
                x2 = x0 + 0.1;
                tdp = x2;
            }
            if iter > 2 {
                tdp = x2;
            }
            r = resid(tdp)?;
            if iter == 1 {
                y1 = r;
            }
            if iter > 1 {
                let y2 = r;
                let x3 = x2 - y2 / (y2 - y1) * (x2 - x1);
                y1 = y2;
                x1 = x2;
                x2 = x3;
            }
            iter += 1;
        }
        Ok(tdp)
    }

    /// The wet-bulb energy-balance residual (upstream `WetBulbSolver::call`).
    fn wetbulb_residual(&self, twb: f64, t: f64, p: f64, psi_w: f64) -> Result<f64> {
        let w = EPSILON * psi_w / (1.0 - psi_w);
        let v_bar_w = self.molar_volume(t, p, psi_w);
        let lhs = self.molar_enthalpy(t, p, psi_w, v_bar_w) * (1.0 + w) / self.m_ha(psi_w);
        let f_wb = self.f_factor(twb, p)?;
        let p_ws_wb = if twb > 273.16 {
            rustprop_if97::psat97(twb)?
        } else {
            ice::psub_ice(twb)
        };
        let p_s_wb = f_wb * p_ws_wb;
        let w_s_wb = EPSILON * p_s_wb / (p - p_s_wb);
        let psi_wb = w_s_wb / (EPSILON + w_s_wb);
        let h_w = if twb > 273.16 {
            let rho_mass = rustprop_if97::rhomass_tp(twb, p)?;
            let rho_molar = rho_mass / self.mm_water;
            self.water.eos.hmolar(twb, rho_molar) / self.mm_water
        } else {
            ice::h_ice(twb, p)
        };
        let v_bar_wb = self.molar_volume(twb, p, psi_wb);
        let rhs = self.molar_enthalpy(twb, p, psi_wb, v_bar_wb) * (1.0 + w_s_wb)
            / self.m_ha(psi_wb)
            + (w - w_s_wb) * h_w;
        Ok(lhs - rhs)
    }

    /// Upstream `WetbulbTemperature`: Brent on the energy balance with the
    /// staged fallbacks.
    fn wetbulb_temperature(&self, t: f64, p: f64, psi_w: f64) -> Result<f64> {
        let (_, tmax) = bounds_of(HaParam::T);
        let tsat = rustprop_if97::tsat97(p)?;
        let tupper = (tmax + 1.0).min(tsat - 0.1);
        let f = |twb: f64| self.wetbulb_residual(twb, t, p, psi_w).unwrap_or(f64::NAN);
        match brent(f, tupper, 100.0, f64::EPSILON, 1e-12, 50) {
            Ok(v) => Ok(v),
            Err(_) => match brent(f, 210.0, tsat - 1.0, 1e-12, 1e-12, 50) {
                Ok(v) => Ok(v),
                Err(_) => brent(f, 130.0 - 30.0, tmax - 1.0, 1e-12, 1e-12, 50),
            },
        }
    }

    /// Transport by the Tsilingiris mixing rule (uses `MM_Air()`, the real
    /// air molar mass, unlike everything else).
    fn transport(&self, t: f64, p: f64, psi_w: f64, conductivity: bool) -> Result<f64> {
        let (mw, ma) = (self.mm_water, self.mm_air);
        let (air_rho, _) = self.air.pt_flash(t, p)?;
        let air_visc = self.air_transport(t, air_rho, p, false)?;
        let sat = self.water.sat().pq_flash(p, 1.0)?;
        let wat_visc = self.water_transport(sat.t, sat.rho_v, p, false)?;
        let phi_av = 2.0f64.sqrt() / 4.0
            * (1.0 + ma / mw).powf(-0.5)
            * (1.0 + (air_visc / wat_visc).sqrt() * (mw / ma).powf(0.25)).powi(2);
        let phi_va = 2.0f64.sqrt() / 4.0
            * (1.0 + mw / ma).powf(-0.5)
            * (1.0 + (wat_visc / air_visc).sqrt() * (ma / mw).powf(0.25)).powi(2);
        let (qa, qw) = if conductivity {
            (
                self.air_transport(t, air_rho, p, true)?,
                self.water_transport(sat.t, sat.rho_v, p, true)?,
            )
        } else {
            (air_visc, wat_visc)
        };
        Ok((1.0 - psi_w) * qa / ((1.0 - psi_w) + psi_w * phi_av)
            + psi_w * qw / (psi_w + (1.0 - psi_w) * phi_va))
    }

    fn air_transport(&self, t: f64, rho: f64, p: f64, conductivity: bool) -> Result<f64> {
        transport_of(&self.air, t, rho, p, conductivity)
    }
    fn water_transport(&self, t: f64, rho: f64, p: f64, conductivity: bool) -> Result<f64> {
        transport_of(&self.water, t, rho, p, conductivity)
    }

    /// Upstream `_HAPropsSI_outputs`.
    #[allow(clippy::too_many_lines)]
    fn outputs(&self, output: HaParam, p: f64, t: f64, psi_w: f64) -> Result<f64> {
        use HaParam::*;
        let m_ha = self.m_ha(psi_w);
        let w = EPSILON * psi_w / (1.0 - psi_w);
        Ok(match output {
            T => t,
            P => p,
            PsiW => psi_w,
            HumRat => w,
            PartialPressureWater => psi_w * p,
            Rh => self.relative_humidity(t, p, psi_w)?,
            Tdp => self.dewpoint_temperature(p, psi_w)?,
            Twb => self.wetbulb_temperature(t, p, psi_w)?,
            Vda => self.molar_volume(t, p, psi_w) * (1.0 + w) / m_ha,
            Vha => self.molar_volume(t, p, psi_w) / m_ha,
            Enthalpy | EnthalpyHa => {
                let v_bar = self.molar_volume(t, p, psi_w);
                let h_bar = self.molar_enthalpy(t, p, psi_w, v_bar);
                if output == Enthalpy {
                    h_bar * (1.0 + w) / m_ha
                } else {
                    h_bar / m_ha
                }
            }
            InternalEnergy | InternalEnergyHa => {
                let v_bar = self.molar_volume(t, p, psi_w);
                let u_bar = self.molar_internal_energy(t, p, psi_w, v_bar);
                if output == InternalEnergy {
                    u_bar * (1.0 + w) / m_ha
                } else {
                    u_bar / m_ha
                }
            }
            Entropy | EntropyHa => {
                let v_bar = self.molar_volume(t, p, psi_w);
                let s_bar = self.molar_entropy(t, p, psi_w, v_bar);
                if output == Entropy {
                    s_bar * (1.0 + w) / m_ha
                } else {
                    s_bar / m_ha
                }
            }
            Visc => self.transport(t, p, psi_w, false)?,
            Cond => self.transport(t, p, psi_w, true)?,
            CpHa | Cp => {
                let dt = 1e-3;
                let v1 = self.molar_volume(t - dt, p, psi_w);
                let h1 = self.molar_enthalpy(t - dt, p, psi_w, v1);
                let v2 = self.molar_volume(t + dt, p, psi_w);
                let h2 = self.molar_enthalpy(t + dt, p, psi_w, v2);
                let cp_ha = (h2 - h1) / (2.0 * dt) / m_ha;
                if output == Cp {
                    cp_ha * (1.0 + w)
                } else {
                    cp_ha
                }
            }
            CvHa | Cv => {
                let dt = 1e-3;
                let v_bar = self.molar_volume(t, p, psi_w);
                let p1 = self.pressure(t - dt, v_bar, psi_w);
                let u1 = self.molar_internal_energy(t - dt, p1, psi_w, v_bar);
                let p2 = self.pressure(t + dt, v_bar, psi_w);
                let u2 = self.molar_internal_energy(t + dt, p2, psi_w, v_bar);
                let cv_ha = (u2 - u1) / (2.0 * dt) / m_ha;
                if output == Cv {
                    cv_ha * (1.0 + w)
                } else {
                    cv_ha
                }
            }
            IsentropicExponent => {
                let dv = 1e-8;
                let cp = self.outputs(CpHa, p, t, psi_w)?;
                let cv = self.outputs(CvHa, p, t, psi_w)?;
                let v_bar = self.molar_volume(t, p, psi_w);
                let p1 = self.pressure(t, v_bar - dv, psi_w);
                let p2 = self.pressure(t, v_bar + dv, psi_w);
                let dpdv = (p2 - p1) / (2.0 * dv);
                -cp / cv * dpdv * v_bar / p
            }
            SpeedOfSound => {
                let dv = 1e-8;
                let cp = self.outputs(CpHa, p, t, psi_w)?;
                let cv = self.outputs(CvHa, p, t, psi_w)?;
                let v_bar = self.molar_volume(t, p, psi_w);
                let p1 = self.pressure(t, v_bar - dv, psi_w);
                let p2 = self.pressure(t, v_bar + dv, psi_w);
                let dvdrho = -v_bar * v_bar / m_ha;
                let dpdrho = (p2 - p1) / (2.0 * dv) * dvdrho;
                (1.0 / m_ha * cp / cv * dpdrho * m_ha).sqrt()
            }
            CompressibilityFactor => {
                let v_bar = self.molar_volume(t, p, psi_w);
                p * v_bar / (R_BAR * t)
            }
        })
    }

    /// Upstream `_HAPropsSI_inputs` — resolve (T, psi_w) from (p, two
    /// non-pressure inputs).
    fn resolve_inputs(&self, p: f64, keys: [HaParam; 2], vals: [f64; 2]) -> Result<(f64, f64)> {
        use HaParam::*;
        // Regime (A): T given.
        if let Some(ti) = keys.iter().position(|k| *k == T) {
            let t = vals[ti];
            let (ok, ov) = (keys[1 - ti], vals[1 - ti]);
            let psi_w = match ok {
                Rh | HumRat | Tdp => self.mole_fraction_water(t, p, ok, ov)?,
                _ => {
                    // Iterate on W.
                    let resid = |w: f64| -> f64 {
                        (|| -> Result<f64> {
                            let psi = self.mole_fraction_water(t, p, HumRat, w)?;
                            Ok(self.outputs(ok, p, t, psi)? - ov)
                        })()
                        .unwrap_or(f64::NAN)
                    };
                    let w = match secant(resid, 0.0001, 0.00001, 1e-14, 100) {
                        Ok(w) if w.is_finite() => w,
                        _ => {
                            let psi_w_sat = self.mole_fraction_water(t, p, Rh, 1.0)?;
                            let w_max = psi_w_sat * EPSILON / (1.0 - psi_w_sat);
                            brent(resid, 0.0, w_max, 1e-7, 1e-7, 50)?
                        }
                    };
                    if !w.is_finite() {
                        return Err(Error::Value("Iterative value for W is invalid".into()));
                    }
                    self.mole_fraction_water(t, p, HumRat, w)?
                }
            };
            return Ok((t, psi_w));
        }

        // Regime (B): no T — iterate on the dry-bulb temperature.
        let main = keys
            .iter()
            .position(|k| *k == HumRat)
            .or_else(|| keys.iter().position(|k| *k == Tdp))
            .or_else(|| keys.iter().position(|k| *k == Rh))
            .ok_or_else(|| Error::Value(
                "Sorry, but currently at least one of the variables as an input to HAPropsSI() must be temperature, relative humidity, humidity ratio, or dewpoint\n  Eventually will add a 2-D NR solver to find T and psi_w simultaneously, but not included now"
                    .into(),
            ))?;
        let (main_key, main_val) = (keys[main], vals[main]);
        let (sec_key, sec_val) = (keys[1 - main], vals[1 - main]);
        let water_content = |k: HaParam| matches!(k, HumRat | Tdp | Rh);
        if water_content(sec_key) {
            let has_rh = keys.contains(&Rh);
            let has_tdp = keys.contains(&Tdp);
            let has_humrat = keys.contains(&HumRat);
            let valid_two_water = has_rh && (has_tdp || has_humrat) && !(has_tdp && has_humrat);
            if !valid_two_water {
                return Err(Error::Value(
                    "Sorry, but cannot provide two inputs that are both water-content (humidity ratio, relative humidity, absolute humidity"
                        .into(),
                ));
            }
        }

        let (mut t_min, mut t_max) = bounds_of(T);
        match main_key {
            Rh => {
                if main_val < 1e-10 {
                    t_max = 640.0;
                    if sec_key == Tdp {
                        return Err(Error::Value(
                            "For dry air, dewpoint is an invalid input variable\n".into(),
                        ));
                    }
                } else {
                    let flash = &self.water;
                    t_max = flash.sat().pq_flash(p, 0.0)?.t - 1.0;
                }
            }
            HumRat => {
                if main_val < 1e-10 {
                    t_min = 135.0;
                    t_max = 1000.0;
                } else {
                    let psi_w_sat = main_val / (EPSILON + main_val);
                    let pp_water_sat = psi_w_sat * p;
                    let mut t0 = if pp_water_sat > self.water.fluid().eos.sat_min_liquid.p {
                        rustprop_if97::tsat97(pp_water_sat)?
                    } else {
                        230.0
                    };
                    // Secant_Tdb_at_saturated_W: T where saturated psi
                    // equals psi_w_sat.
                    let resid = |t: f64| -> f64 {
                        self.mole_fraction_water(t, p, Rh, 1.0)
                            .map_or(f64::NAN, |ps| (ps - psi_w_sat) / psi_w_sat)
                    };
                    t0 = match secant(resid, t0, 0.1, 1e-7, 100) {
                        Ok(v) => v,
                        Err(_) => brent(resid, 100.0, 640.0, 1e-15, 1e-10, 100)?,
                    };
                    if !t0.is_finite() {
                        return Err(Error::Value("Intermediate value for Tdb is invalid".into()));
                    }
                    t_min = t0;
                }
            }
            Tdp => {
                let psi_w = self.mole_fraction_water(f64::NAN, p, Tdp, main_val)?;
                t_min = self.dewpoint_temperature(p, psi_w)?;
            }
            _ => unreachable!(),
        }

        // Brent_HAProps_T on the secondary output.
        let resid = |t: f64| -> f64 {
            (|| -> Result<f64> {
                let psi = self.mole_fraction_water(t, p, main_key, main_val)?;
                Ok(self.outputs(sec_key, p, t, psi)? - sec_val)
            })()
            .unwrap_or(f64::NAN)
        };
        let t = brent(resid, t_min, t_max, 1e-15, 1e-10, 50).or_else(|_| {
            // Endpoint-repaired secant fallback per upstream's helper.
            secant(resid, 0.5 * (t_min + t_max), 0.1, 1e-7, 50)
        })?;
        let r_check = resid(t);
        if !r_check.is_finite() || r_check.abs() > 1e-4 * sec_val.abs() + 1e-6 {
            return Err(Error::Solution(format!(
                "Brent_HAProps_T: no temperature in [{t_min}, {t_max}] K yields the requested output [{sec_val}] (closest residual {r_check}); the input is out of range"
            )));
        }
        let psi_w = self.mole_fraction_water(t, p, main_key, main_val)?;
        Ok((t, psi_w))
    }

    /// The `HAPropsSI` entry (upstream 2136-2267), minus the swallow-into-inf
    /// error transport: errors return as `Err` with upstream's messages.
    /// (Seven arguments mirror upstream's signature exactly.)
    #[allow(clippy::too_many_arguments)]
    pub fn ha_props_si(
        &self,
        output: &str,
        n1: &str,
        v1: f64,
        n2: &str,
        v2: f64,
        n3: &str,
        v3: f64,
    ) -> Result<f64> {
        let out = name_to_type(output)?;
        let k1 = name_to_type(n1)?;
        let k2 = name_to_type(n2)?;
        let k3 = name_to_type(n3)?;
        // Trivial echo BEFORE any validation (upstream's order).
        if out == k1 {
            return Ok(v1);
        }
        if out == k2 {
            return Ok(v2);
        }
        if out == k3 {
            return Ok(v3);
        }
        let (keys, vals, p) = if k1 == HaParam::P {
            ([k2, k3], [v2, v3], v1)
        } else if k2 == HaParam::P {
            ([k1, k3], [v1, v3], v2)
        } else if k3 == HaParam::P {
            ([k1, k2], [v1, v2], v3)
        } else {
            return Err(Error::Value(
                "Pressure must be one of the inputs to HAPropsSI".into(),
            ));
        };
        if keys[0] == keys[1] {
            return Err(Error::Value(
                "Other two inputs to HAPropsSI aside from pressure cannot be the same".into(),
            ));
        }
        for (k, v) in keys.iter().zip(vals.iter()) {
            if !check_bounds(*k, *v) {
                let (lo, hi) = bounds_of(*k);
                return Err(Error::Value(format!(
                    "The input for key ({k:?}) with value ({v}) is outside the range of validity: ({lo}) to ({hi})"
                )));
            }
        }
        let (t, psi_w) = self.resolve_inputs(p, keys, vals)?;
        if !check_bounds(HaParam::P, p) {
            let (lo, hi) = bounds_of(HaParam::P);
            return Err(Error::Value(format!(
                "The pressure value ({p}) is outside the range of validity: ({lo}) to ({hi})"
            )));
        }
        if !check_bounds(HaParam::T, t) {
            let (lo, hi) = bounds_of(HaParam::T);
            return Err(Error::Value(format!(
                "The temperature value ({t}) is outside the range of validity: ({lo}) to ({hi})"
            )));
        }
        if !check_bounds(HaParam::PsiW, psi_w) {
            let (lo, hi) = bounds_of(HaParam::PsiW);
            return Err(Error::Value(format!(
                "The water mole fraction value ({psi_w}) is outside the range of validity: ({lo}) to ({hi})"
            )));
        }
        let val = self.outputs(out, p, t, psi_w)?;
        if !check_bounds(out, val) {
            let (lo, hi) = bounds_of(out);
            return Err(Error::Value(format!(
                "The output for key ({out:?}) with value ({val}) is outside the range of validity: ({lo}) to ({hi})"
            )));
        }
        if !val.is_finite() {
            return Err(Error::Value("Invalid value about to be returned".into()));
        }
        Ok(val)
    }
}

/// Transport of one fluid at a resolved state through the ported models.
fn transport_of(flash: &PtFlash, t: f64, rho: f64, p: f64, conductivity: bool) -> Result<f64> {
    let data = flash.fluid();
    let tr = data
        .transport
        .as_ref()
        .ok_or_else(|| Error::Value("transport model missing".into()))?;
    use rustprop_core::fluid::TransportModel;
    if conductivity {
        let TransportModel::Model(c) = &tr.conductivity else {
            return Err(Error::Value("conductivity model missing".into()));
        };
        let v = match &tr.viscosity {
            TransportModel::Model(v) => Some(v),
            _ => None,
        };
        rustprop_heos::transport::conductivity(&flash.eos, data, c, v, t, rho, p, None)
    } else {
        let TransportModel::Model(v) = &tr.viscosity else {
            return Err(Error::Value("viscosity model missing".into()));
        };
        rustprop_heos::transport::viscosity(&flash.eos, data, v, t, rho, p, None)
    }
}

/// Hardcoded cross-virials (ASHRAE RP-1485).
fn b_aw(t: f64) -> f64 {
    let a = [0.0, 0.665687e2, -0.238834e3, -0.176755e3];
    let b = [0.0, -0.237, -1.048, -3.183];
    let (rhobarstar, tstar) = (1000.0, 100.0);
    1.0 / rhobarstar
        * (a[1] * (t / tstar).powf(b[1])
            + a[2] * (t / tstar).powf(b[2])
            + a[3] * (t / tstar).powf(b[3]))
        / 1000.0
}
fn db_aw_dt(t: f64) -> f64 {
    let a = [0.0, 0.665687e2, -0.238834e3, -0.176755e3];
    let b = [0.0, -0.237, -1.048, -3.183];
    let (rhobarstar, tstar) = (1000.0, 100.0);
    1.0 / rhobarstar / tstar
        * (a[1] * b[1] * (t / tstar).powf(b[1] - 1.0)
            + a[2] * b[2] * (t / tstar).powf(b[2] - 1.0)
            + a[3] * b[3] * (t / tstar).powf(b[3] - 1.0))
        / 1000.0
}
fn c_aaw(t: f64) -> f64 {
    let c = [
        0.0,
        0.482737e3,
        0.105678e6,
        -0.656394e8,
        0.294442e11,
        -0.319317e13,
    ];
    let rhobarstar = 1000.0;
    let mut summer = 0.0;
    for (i, ci) in c.iter().enumerate().skip(1) {
        summer += ci * t.powi(1 - i as i32);
    }
    1.0 / rhobarstar / rhobarstar * summer / 1e6
}
fn dc_aaw_dt(t: f64) -> f64 {
    let c = [
        0.0,
        0.482737e3,
        0.105678e6,
        -0.656394e8,
        0.294442e11,
        -0.319317e13,
    ];
    let rhobarstar = 1000.0;
    let mut summer = 0.0;
    for (i, ci) in c.iter().enumerate().skip(2) {
        summer += ci * (1.0 - i as f64) * t.powi(-(i as i32));
    }
    1.0 / rhobarstar / rhobarstar * summer / 1e6
}
fn c_aww(t: f64) -> f64 {
    let d = [0.0, -0.1072887e2, 0.347804e4, -0.383383e6, 0.334060e8];
    let mut summer = 0.0;
    for (i, di) in d.iter().enumerate().skip(1) {
        summer += di * t.powi(1 - i as i32);
    }
    -(summer.exp()) / 1e6
}
fn dc_aww_dt(t: f64) -> f64 {
    let d = [0.0, -0.1072887e2, 0.347804e4, -0.383383e6, 0.334060e8];
    let mut summer1 = 0.0;
    for (i, di) in d.iter().enumerate().skip(1) {
        summer1 += di * t.powi(1 - i as i32);
    }
    let mut summer2 = 0.0;
    for (i, di) in d.iter().enumerate().skip(2) {
        summer2 += di * (1.0 - i as f64) * t.powi(-(i as i32));
    }
    -(summer1.exp()) * summer2 / 1e6
}

fn psi_w_pow3(x: f64) -> f64 {
    x * x * x
}

// Keep the state type import used indirectly by the flashes.
#[allow(unused_imports)]
use HeosState as _HeosStateUsed;

pub use rustprop_core::UPSTREAM_VERSION;
