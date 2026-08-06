//! Backward equations — T(p,h) and T(p,s) for regions 1/2/3 (R7-97 and
//! SR3-03), p(h,s) for regions 1/2/3 (SR2-01, SR4-04), Tsat(h,s) for region 4,
//! the h/s boundary curves, and the (p,X) / (h,s) region determinations.
//! Ported from the `Backwards` namespace and the general machinery of IF97.h.

// Fidelity over style: 0.7853 is the upstream 2bS scale factor (not pi/4), and
// the branch chains mirror upstream's deliberate condition structure, some
// with identical arms for distinct physical conditions.
#![allow(
    clippy::approx_constant,
    clippy::collapsible_match,
    clippy::if_same_then_else
)]

use crate::tables::{
    BackResid, COEFF_1H, COEFF_1HS, COEFF_1S, COEFF_2AH, COEFF_2AHS, COEFF_2AS, COEFF_2BH,
    COEFF_2BHS, COEFF_2BS, COEFF_2CH, COEFF_2CHS, COEFF_2CS, COEFF_3AH, COEFF_3AHS, COEFF_3AS,
    COEFF_3BH, COEFF_3BHS, COEFF_3BS, COEFF_B2ABHS, COEFF_B2C3BHS, COEFF_B3A4HS, COEFF_B13HS,
    COEFF_B14HS, COEFF_T4HS, COEFF_TB23HS, HTMAX_N, REGION2AB_N, REGION2B2C_N, REGION3AB_N,
};
use crate::{
    H23MAX, H23MIN, HFTRIP, HGTRIP, P_FACT, P2AMAX, P2BCMIN, P23MIN, PCRIT, PMAX, PMIN, Prop,
    R_FACT, Region, S2BC, S13MIN, S23MAX, S23MIN, SCRIT, SFT23, SFTRIP, SGT23, SGTRIP, SMAX, SMIN,
    STPMAX, SatState, TCRIT, TMAX, TMIN, b23_p_from_t, gibbs, powi, region_output, region3, tsat97,
};
use rustprop_core::Error;

// ---------------------------------------------------------------------------
// Generic backward table evaluation (upstream `BackwardsRegion`)
// ---------------------------------------------------------------------------

struct BackTable {
    p_star: f64,
    t_star: f64,
    x_star: f64,
    h_star: f64,
    s_star: f64,
    s2_star: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    data: &'static [BackResid],
}

const UNSET: f64 = 0.0; // upstream leaves unused members uninitialized

const DEFAULTS: BackTable = BackTable {
    p_star: UNSET,
    t_star: UNSET,
    x_star: UNSET,
    h_star: UNSET,
    s_star: UNSET,
    s2_star: UNSET,
    a: UNSET,
    b: UNSET,
    c: UNSET,
    d: UNSET,
    e: UNSET,
    f: UNSET,
    data: &[],
};

impl BackTable {
    /// T(p,h) or T(p,s) — and v(p,h)/v(p,s) in region 3 (upstream `T_pX`).
    fn t_px(&self, p: f64, x: f64) -> f64 {
        let (pi, eta) = (p / self.p_star, x / self.x_star);
        let mut summer = 0.0;
        for el in self.data {
            summer +=
                el.n * (pi + self.a).powf(el.i) * (eta + self.b).powf(el.j) * self.f.powf(el.j);
        }
        summer * self.t_star
    }

    /// Boundary h'(s) / h"(s) fits (upstream `h_s`).
    /// c=1,e=0: straight sum; c>1,e=0: power fit; c=1,e=1: exp fit.
    fn h_s(&self, s: f64) -> f64 {
        let (sigma1, sigma2) = (s / self.s_star, s / self.s2_star);
        let mut summer = 0.0;
        for el in self.data {
            summer +=
                el.n * (sigma1.powf(self.d) + self.a).powf(el.i) * (sigma2 + self.b).powf(el.j);
        }
        ((1.0 - self.e) * summer.powf(self.c) + self.e * summer.exp()) * self.h_star
    }

    /// p(h,s) (upstream `p_hs`).
    fn p_hs(&self, h: f64, s: f64) -> f64 {
        let (eta, sigma) = (h / self.h_star, s / self.s_star);
        let mut summer = 0.0;
        for el in self.data {
            summer += el.n * (eta + self.a).powf(el.i) * (sigma + self.b).powf(el.j);
        }
        summer.powf(self.c) * self.p_star
    }

    /// Tb23(h,s) / Tsat(h,s) (upstream `t_hs`).
    fn t_hs(&self, h: f64, s: f64) -> f64 {
        let (eta, sigma) = (h / self.h_star, s / self.s_star);
        let mut summer = 0.0;
        for el in self.data {
            summer += el.n * (eta + self.a).powf(el.i) * (sigma + self.b).powf(el.j);
        }
        summer * self.t_star
    }
}

macro_rules! back_table {
    ($name:ident, $data:expr, { $($field:ident : $value:expr),* $(,)? }) => {
        const $name: BackTable = BackTable {
            $($field: $value,)*
            data: $data,
            ..DEFAULTS
        };
    };
}

back_table!(B1H, COEFF_1H, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 2500.0 * R_FACT, a: 0.0, b: 1.0, f: 1.0 });
back_table!(B1S, COEFF_1S, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 1.0 * R_FACT, a: 0.0, b: 2.0, f: 1.0 });
back_table!(B1HS, COEFF_1HS, { p_star: 100.0 * P_FACT, h_star: 3400.0 * R_FACT, s_star: 7.6 * R_FACT, a: 0.05, b: 0.05, c: 1.0 });
back_table!(B2AH, COEFF_2AH, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 2000.0 * R_FACT, a: 0.0, b: -2.1, f: 1.0 });
back_table!(B2BH, COEFF_2BH, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 2000.0 * R_FACT, a: -2.0, b: -2.6, f: 1.0 });
back_table!(B2CH, COEFF_2CH, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 2000.0 * R_FACT, a: 25.0, b: -1.8, f: 1.0 });
back_table!(B2AS, COEFF_2AS, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 2.0 * R_FACT, a: 0.0, b: -2.0, f: 1.0 });
back_table!(B2BS, COEFF_2BS, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 0.7853 * R_FACT, a: 0.0, b: -10.0, f: -1.0 });
back_table!(B2CS, COEFF_2CS, { p_star: 1.0 * P_FACT, t_star: 1.0, x_star: 2.9251 * R_FACT, a: 0.0, b: -2.0, f: -1.0 });
back_table!(B2AHS, COEFF_2AHS, { p_star: 4.0 * P_FACT, h_star: 4200.0 * R_FACT, s_star: 12.0 * R_FACT, a: -0.5, b: -1.2, c: 4.0 });
back_table!(B2BHS, COEFF_2BHS, { p_star: 100.0 * P_FACT, h_star: 4100.0 * R_FACT, s_star: 7.9 * R_FACT, a: -0.6, b: -1.01, c: 4.0 });
back_table!(B2CHS, COEFF_2CHS, { p_star: 100.0 * P_FACT, h_star: 3500.0 * R_FACT, s_star: 5.9 * R_FACT, a: -0.7, b: -1.1, c: 4.0 });
back_table!(B3AH, COEFF_3AH, { p_star: 100.0 * P_FACT, t_star: 760.0, x_star: 2300.0 * R_FACT, a: 0.240, b: -0.615, f: 1.0 });
back_table!(B3BH, COEFF_3BH, { p_star: 100.0 * P_FACT, t_star: 860.0, x_star: 2800.0 * R_FACT, a: 0.298, b: -0.720, f: 1.0 });
back_table!(B3AS, COEFF_3AS, { p_star: 100.0 * P_FACT, t_star: 760.0, x_star: 4.4 * R_FACT, a: 0.240, b: -0.703, f: 1.0 });
back_table!(B3BS, COEFF_3BS, { p_star: 100.0 * P_FACT, t_star: 860.0, x_star: 5.3 * R_FACT, a: 0.760, b: -0.818, f: 1.0 });
back_table!(B3AHS, COEFF_3AHS, { p_star: 99.0 * P_FACT, h_star: 2300.0 * R_FACT, s_star: 4.4 * R_FACT, a: -1.01, b: -0.750, c: 1.0 });
back_table!(B3BHS, COEFF_3BHS, { p_star: 16.6 * P_FACT, h_star: 2800.0 * R_FACT, s_star: 5.3 * R_FACT, a: -0.681, b: -0.792, c: -1.0 });
back_table!(B4HS, COEFF_T4HS, { h_star: 2800.0 * R_FACT, s_star: 9.2 * R_FACT, t_star: 550.0, a: -0.119, b: -1.07 });
back_table!(B14HS, COEFF_B14HS, { h_star: 1700.0 * R_FACT, s_star: 3.8 * R_FACT, s2_star: 3.8 * R_FACT, a: -1.09, b: 0.366E-4, c: 1.0, d: 1.0, e: 0.0 });
back_table!(B3A4HS, COEFF_B3A4HS, { h_star: 1700.0 * R_FACT, s_star: 3.8 * R_FACT, s2_star: 3.8 * R_FACT, a: -1.09, b: 0.366E-4, c: 1.0, d: 1.0, e: 0.0 });
back_table!(B2AB4HS, COEFF_B2ABHS, { h_star: 2800.0 * R_FACT, s_star: 5.21 * R_FACT, s2_star: 9.2 * R_FACT, a: -0.513, b: -0.524, c: 1.0, d: -1.0, e: 1.0 });
back_table!(B2C3B4HS, COEFF_B2C3BHS, { h_star: 2800.0 * R_FACT, s_star: 5.9 * R_FACT, s2_star: 5.9 * R_FACT, a: -1.02, b: -0.726, c: 4.0, d: 1.0, e: 0.0 });
back_table!(B13HS, COEFF_B13HS, { h_star: 1700.0 * R_FACT, s_star: 3.8 * R_FACT, s2_star: 3.8 * R_FACT, a: -0.884, b: -0.864, c: 1.0, d: 1.0, e: 0.0 });
back_table!(B23HS, COEFF_TB23HS, { h_star: 3000.0 * R_FACT, s_star: 5.3 * R_FACT, t_star: 900.0, a: -0.727, b: -0.864 });

// ---------------------------------------------------------------------------
// Simple boundary curves
// ---------------------------------------------------------------------------

/// Upstream `H2b2c_p`: h [J/kg] on the 2b/2c boundary at p [Pa].
fn h2b2c_p(p: f64) -> f64 {
    let (p_star, h_star) = (1.0 * P_FACT, 1.0 * R_FACT);
    let pi = p / p_star;
    let eta = REGION2B2C_N[3] + ((pi - REGION2B2C_N[4]) / REGION2B2C_N[2]).sqrt();
    eta * h_star
}

/// Upstream `H3ab_p`: h [J/kg] on the 3a/3b boundary at p [Pa].
fn h3ab_p(p: f64) -> f64 {
    let (p_star, h_star) = (1.0 * P_FACT, 1.0 * R_FACT);
    let pi = p / p_star;
    let eta = REGION3AB_N[0]
        + REGION3AB_N[1] * pi
        + REGION3AB_N[2] * pi * pi
        + REGION3AB_N[3] * pi * pi * pi;
    eta * h_star
}

/// Upstream `H2ab_s`: h [J/kg] on the 2a/2b boundary at s [J/kg/K].
fn h2ab_s(s: f64) -> f64 {
    let (s_star, h_star) = (1.0 * R_FACT, 1.0 * R_FACT);
    let sigma = s / s_star;
    let eta = REGION2AB_N[0]
        + REGION2AB_N[1] * sigma
        + REGION2AB_N[2] * powi(sigma, 2)
        + REGION2AB_N[3] * powi(sigma, 3);
    eta * h_star
}

/// Upstream `H13_s`.
fn h13_s(s: f64) -> f64 {
    B13HS.h_s(s)
}

/// Upstream `Hsat_s`: saturated-liquid/vapor enthalpy along s.
fn hsat_s(s: f64) -> Result<f64, Error> {
    if s < 0.0 {
        Err(Error::OutOfRange("Entropy out of range".into()))
    } else if s <= SFT23 {
        Ok(B14HS.h_s(s))
    } else if s <= SCRIT {
        Ok(B3A4HS.h_s(s))
    } else if s <= S2BC {
        Ok(B2C3B4HS.h_s(s))
    } else if s <= SGTRIP {
        Ok(B2AB4HS.h_s(s))
    } else {
        Err(Error::OutOfRange("Entropy out of range".into()))
    }
}

// ---------------------------------------------------------------------------
// Region determination from (p, h) or (p, s)
// ---------------------------------------------------------------------------

/// Upstream `RegionDetermination_pX`.
pub(crate) fn region_determination_px(p: f64, x: f64, inkey: Prop) -> Result<Region, Error> {
    if p < PMIN || p > PMAX {
        return Err(Error::OutOfRange("Pressure out of range".into()));
    }
    let xmin = gibbs::REGION1.output(inkey, TMIN, p)?;
    let xmax = gibbs::REGION2.output(inkey, TMAX, p)?;
    if x < xmin || x > (xmax + 1.0e-10) {
        return Err(if inkey == Prop::Hmass {
            Error::OutOfRange("Enthalpy out of range".into())
        } else {
            Error::OutOfRange("Entropy out of range".into())
        });
    }

    let mut xliq = 0.0;
    let mut xvap = 0.0;
    if p <= PCRIT {
        // Check saturation dome first
        let tsat = tsat97(p)?;
        xliq = region_output(inkey, tsat, p, SatState::Liquid)?; // Regions 1 & 3
        xvap = region_output(inkey, tsat, p, SatState::Vapor)?; // Regions 2 & 3
        if xliq <= x && x <= xvap {
            return Ok(Region::R4); // Within saturation dome (inclusive)
        }
    }

    if p <= P23MIN {
        if x < xliq {
            Ok(Region::R1)
        } else if x > xvap {
            Ok(Region::R2)
        } else {
            Ok(Region::R4) // already handled above
        }
    } else if x <= gibbs::REGION1.output(inkey, crate::T23MIN, p)? {
        Ok(Region::R1)
    } else if x >= gibbs::REGION2.output(inkey, crate::b23_t_from_p(p), p)? {
        Ok(Region::R2)
    } else {
        Ok(Region::R3) // R4 has already been accounted for above
    }
}

/// Upstream `BackwardRegion`: region as an integer, for `Region_ph`/`Region_ps`.
pub(crate) fn backward_region(p: f64, x: f64, inkey: Prop) -> Result<i32, Error> {
    if inkey != Prop::Hmass && inkey != Prop::Smass {
        return Err(Error::Input(
            "Backward Formulas take variable inputs of Enthalpy or Entropy only.".into(),
        ));
    }
    Ok(match region_determination_px(p, x, inkey)? {
        Region::R1 => 1,
        Region::R2 => 2,
        Region::R3 => 3,
        Region::R4 => 4,
        Region::R5 => 0,
    })
}

// ---------------------------------------------------------------------------
// Backward temperature: T(p,h) / T(p,s)
// ---------------------------------------------------------------------------

/// Upstream `RegionOutputBackward` — returns temperature only.
pub(crate) fn region_output_backward(
    p: f64,
    x: f64,
    inkey: Prop,
    clip: bool,
    state: SatState,
) -> Result<f64, Error> {
    if inkey != Prop::Hmass && inkey != Prop::Smass {
        return Err(Error::Input(
            "Backward Formulas take variable inputs of Enthalpy or Entropy only.".into(),
        ));
    }

    // The reverse functions carry ±25 mK of uncertainty; near saturation, clip
    // the result to the correct side of the curve (see the upstream comment).
    let mut tmin = TMIN;
    let mut tmax = TMAX;
    let eps = 1.0e-6; // saturation temperature offset of .001 mK
    if p < PCRIT && clip {
        let tsat = tsat97(p)?;
        tmin = tsat + eps;
        tmax = tsat - eps;
    } else if p == PCRIT {
        // Handle cases directly on the critical point
        match inkey {
            Prop::Hmass => {
                if x == h3ab_p(PCRIT) {
                    return Ok(TCRIT);
                }
            }
            Prop::Smass => {
                if x == SCRIT {
                    return Ok(TCRIT);
                }
            }
            _ => {}
        }
    }

    let mut region = region_determination_px(p, x, inkey)?;

    // Override region if a saturated state is requested
    if state == SatState::Liquid {
        region = if p <= P23MIN { Region::R1 } else { Region::R3 };
    } else if state == SatState::Vapor {
        region = if p <= P23MIN { Region::R2 } else { Region::R3 };
    }

    match region {
        Region::R1 => {
            if inkey == Prop::Hmass {
                Ok(tmax.min(B1H.t_px(p, x)))
            } else {
                Ok(tmax.min(B1S.t_px(p, x)))
            }
        }
        Region::R2 => {
            if inkey == Prop::Hmass {
                if p <= P2AMAX {
                    Ok(tmin.max(B2AH.t_px(p, x)))
                } else if p <= P2BCMIN {
                    Ok(tmin.max(B2BH.t_px(p, x)))
                } else if x >= h2b2c_p(p) {
                    Ok(tmin.max(B2BH.t_px(p, x)))
                } else {
                    Ok(tmin.max(B2CH.t_px(p, x)))
                }
            } else if p <= P2AMAX {
                Ok(tmin.max(B2AS.t_px(p, x)))
            } else if p <= P2BCMIN {
                Ok(tmin.max(B2BS.t_px(p, x)))
            } else if x >= S2BC {
                Ok(tmin.max(B2BS.t_px(p, x)))
            } else {
                Ok(tmin.max(B2CS.t_px(p, x)))
            }
        }
        Region::R3 => {
            if inkey == Prop::Hmass {
                if x <= h3ab_p(p) {
                    Ok(tmax.min(B3AH.t_px(p, x)))
                } else {
                    Ok(tmin.max(B3BH.t_px(p, x)))
                }
            } else if x <= SCRIT {
                Ok(tmax.min(B3AS.t_px(p, x)))
            } else {
                Ok(tmin.max(B3BS.t_px(p, x)))
            }
        }
        Region::R4 => tsat97(p), // just return Tsat in the 2-phase region
        Region::R5 => Err(Error::OutOfRange("Unable to match region".into())),
    }
}

// ---------------------------------------------------------------------------
// Generic backward property: Y(p, h) / Y(p, s), including the 2-phase dome
// ---------------------------------------------------------------------------

/// Upstream `Y_pX`: h(p,s), s(p,h), rho(p,h), rho(p,s), Q(p,X), ...
pub(crate) fn y_px(outkey: Prop, p: f64, x: f64, inkey: Prop) -> Result<f64, Error> {
    if inkey != Prop::Hmass && inkey != Prop::Smass {
        return Err(Error::Input(
            "Reverse state cannot be determined for these inputs.".into(),
        ));
    }

    let t = region_output_backward(p, x, inkey, false, SatState::None)?;

    if inkey == outkey {
        return Ok(x); // trivial result
    }

    match region_determination_px(p, x, inkey)? {
        Region::R4 => {
            // Saturation dome: get liquid/vapor values directly from region eqs.
            let tsat = tsat97(p)?;
            let (xliq, xvap) = if p > P23MIN {
                (
                    region3::output(inkey, tsat, p, SatState::Liquid)?,
                    region3::output(inkey, tsat, p, SatState::Vapor)?,
                )
            } else {
                (
                    gibbs::REGION1.output(inkey, tsat, p)?,
                    gibbs::REGION2.output(inkey, tsat, p)?,
                )
            };
            let q4 = 1.0_f64.min(0.0_f64.max((x - xliq) / (xvap - xliq)));
            match outkey {
                Prop::Dmass => {
                    let tl = region_output_backward(p, xliq, inkey, false, SatState::Liquid)?;
                    let tv = region_output_backward(p, xvap, inkey, false, SatState::Vapor)?;
                    let (yliq, yvap) = if p > P23MIN {
                        (
                            1.0 / region3::output(outkey, tl, p, SatState::Liquid)?,
                            1.0 / region3::output(outkey, tv, p, SatState::Vapor)?,
                        )
                    } else {
                        (
                            1.0 / gibbs::REGION1.output(outkey, tl, p)?,
                            1.0 / gibbs::REGION2.output(outkey, tv, p)?,
                        )
                    };
                    Ok(1.0 / (yliq * (1.0 - q4) + q4 * yvap)) // mixture density
                }
                Prop::T => Ok(tsat),
                Prop::Q => Ok(q4),
                Prop::Hmass | Prop::Smass => {
                    let tl = region_output_backward(p, xliq, inkey, false, SatState::Liquid)?;
                    let tv = region_output_backward(p, xvap, inkey, false, SatState::Vapor)?;
                    let (yliq, yvap) = if p > P23MIN {
                        (
                            region3::output(outkey, tl, p, SatState::Liquid)?,
                            region3::output(outkey, tv, p, SatState::Vapor)?,
                        )
                    } else {
                        (
                            gibbs::REGION1.output(outkey, tl, p)?,
                            gibbs::REGION2.output(outkey, tv, p)?,
                        )
                    };
                    Ok(yliq * (1.0 - q4) + q4 * yvap)
                }
                _ => Err(Error::Input(
                    "2-Phase: Requested output undefined in two-phase region.".into(),
                )),
            }
        }
        Region::R1 => {
            if outkey == Prop::Q {
                Ok(0.0)
            } else {
                gibbs::REGION1.output(outkey, t, p)
            }
        }
        Region::R2 => {
            if outkey == Prop::Q {
                Ok(1.0)
            } else {
                gibbs::REGION2.output(outkey, t, p)
            }
        }
        Region::R3 => {
            let liquid = if inkey == Prop::Hmass {
                x <= h3ab_p(p)
            } else {
                x <= SCRIT
            };
            if liquid {
                if outkey == Prop::Q {
                    Ok(0.0)
                } else {
                    region3::output(outkey, t, p, SatState::Liquid)
                }
            } else if outkey == Prop::Q {
                Ok(1.0)
            } else {
                region3::output(outkey, t, p, SatState::Vapor)
            }
        }
        Region::R5 => Err(Error::Input(
            "Reverse state functions not defined in REGION 5".into(),
        )),
    }
}

/// Upstream `Q_pX`: vapor quality from (p, h/s/u/rho).
pub(crate) fn q_px(p: f64, x: f64, inkey: Prop) -> Result<f64, Error> {
    if p < PMIN || p > PMAX {
        return Err(Error::OutOfRange("Pressure out of range".into()));
    }
    if p < crate::PTRIP {
        return Ok(0.0); // liquid, at all temperatures
    }
    if p > PCRIT {
        let t = match inkey {
            Prop::Hmass | Prop::Smass => {
                region_output_backward(p, x, inkey, false, SatState::None)?
            }
            _ => {
                return Err(Error::Input(
                    "Quality cannot be determined for these inputs.".into(),
                ));
            }
        };
        return if t < TCRIT {
            Ok(0.0) // liquid, at all pressures above critical point
        } else {
            Err(Error::Input(
                "Quality not defined in supercritical region.".into(),
            ))
        };
    }
    let tsat = tsat97(p)?;
    match inkey {
        Prop::Hmass | Prop::Smass | Prop::Umass => {
            let xliq = region_output(inkey, tsat, p, SatState::Liquid)?;
            let xvap = region_output(inkey, tsat, p, SatState::Vapor)?;
            Ok(1.0_f64.min(0.0_f64.max((x - xliq) / (xvap - xliq))))
        }
        Prop::Dmass => {
            let xliq = 1.0 / region_output(Prop::Dmass, tsat, p, SatState::Liquid)?;
            let xvap = 1.0 / region_output(Prop::Dmass, tsat, p, SatState::Vapor)?;
            let x = 1.0 / x;
            Ok(1.0_f64.min(0.0_f64.max((x - xliq) / (xvap - xliq))))
        }
        _ => Err(Error::Input(
            "Quality cannot be determined for these inputs.".into(),
        )),
    }
}

/// Upstream `X_pQ`: mixture h/s/u/rho from (p, Q).
pub(crate) fn x_pq(inkey: Prop, p: f64, q: f64) -> Result<f64, Error> {
    if p < crate::PTRIP || p > PCRIT {
        return Err(Error::OutOfRange("Pressure out of range".into()));
    }
    if q < 0.0 || q > 1.0 {
        return Err(Error::OutOfRange("Quality out of range".into()));
    }
    match inkey {
        Prop::Hmass | Prop::Smass | Prop::Umass => {
            let xliq = region_output(inkey, tsat97(p)?, p, SatState::Liquid)?;
            let xvap = region_output(inkey, tsat97(p)?, p, SatState::Vapor)?;
            Ok(q * xvap + (1.0 - q) * xliq)
        }
        Prop::Dmass => {
            let xliq = 1.0 / region_output(Prop::Dmass, tsat97(p)?, p, SatState::Liquid)?;
            let xvap = 1.0 / region_output(Prop::Dmass, tsat97(p)?, p, SatState::Vapor)?;
            Ok(1.0 / (q * xvap + (1.0 - q) * xliq))
        }
        _ => Err(Error::Input("Mixture property undefined".into())),
    }
}

// ---------------------------------------------------------------------------
// (h, s) domain boundaries and region determination
// ---------------------------------------------------------------------------

/// Upstream `Hmax`: enthalpy bound along Pmax / Tmax.
fn hmax(s: f64) -> Result<f64, Error> {
    let (s_star, h_star) = (1.0 * R_FACT, 1.0 * R_FACT);
    let sigma = s / s_star;
    if s < STPMAX {
        let t = region_output_backward(PMAX, s, Prop::Smass, false, SatState::None)?;
        region_output(Prop::Hmass, t, PMAX, SatState::None)
    } else {
        // Fitted h(s) = a*ln(s) + b/s + c/s^2 + d along the Tmax boundary
        let eta =
            HTMAX_N[0] * sigma.ln() + HTMAX_N[1] / sigma + HTMAX_N[2] / powi(sigma, 2) + HTMAX_N[3];
        Ok(eta * h_star)
    }
}

/// Upstream `Hmin`: enthalpy bound along Pmin / through the dome.
fn hmin(s: f64) -> Result<f64, Error> {
    if s < SGTRIP {
        Ok((s - SFTRIP) * (HGTRIP - HFTRIP) / (SGTRIP - SFTRIP) + HFTRIP)
    } else {
        let t = region_output_backward(PMIN, s, Prop::Smass, false, SatState::None)?;
        region_output(Prop::Hmass, t, PMIN, SatState::None)
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum BackRegion {
    B1,
    B2a,
    B2b,
    B2c,
    B3a,
    B3b,
    B4,
}

/// Upstream `RegionDetermination_HS`.
fn region_determination_hs(h: f64, s: f64) -> Result<BackRegion, Error> {
    if s < SMIN || s > SMAX {
        return Err(Error::OutOfRange("Entropy out of range".into()));
    }
    if h > hmax(s)? || h < hmin(s)? {
        return Err(Error::OutOfRange("Enthalpy out of range".into()));
    }

    if s <= SFT23 {
        // Region 1 check
        if h < hsat_s(s)? {
            Ok(BackRegion::B4)
        } else if s < S13MIN {
            Ok(BackRegion::B1)
        } else if h < h13_s(s) {
            Ok(BackRegion::B1)
        } else {
            Ok(BackRegion::B3a)
        }
    } else if s <= SCRIT {
        // Region 3a check (s < Scrit)
        if h < hsat_s(s)? {
            Ok(BackRegion::B4)
        } else {
            Ok(BackRegion::B3a)
        }
    } else if s <= S23MIN {
        // Region 3b check
        if h < hsat_s(s)? {
            Ok(BackRegion::B4)
        } else {
            Ok(BackRegion::B3b)
        }
    } else if s <= S23MAX {
        // Region 3b/2c check along the B23 curve
        if h < hsat_s(s)? {
            Ok(BackRegion::B4)
        } else if h < H23MIN {
            Ok(BackRegion::B3b)
        } else if h > H23MAX {
            Ok(BackRegion::B2c)
        } else {
            let tb23 = B23HS.t_hs(h, s);
            let pb23 = b23_p_from_t(tb23);
            let p = B2CHS.p_hs(h, s);
            if p > pb23 {
                Ok(BackRegion::B3b)
            } else {
                Ok(BackRegion::B2c)
            }
        }
    } else if s <= S2BC {
        if h < hsat_s(s)? {
            Ok(BackRegion::B4)
        } else {
            Ok(BackRegion::B2c)
        }
    } else if s < SGTRIP {
        // Region 2a/2b above the saturated-vapor curve
        if h < hsat_s(s)? {
            Ok(BackRegion::B4)
        } else if h > h2ab_s(s) {
            Ok(BackRegion::B2b)
        } else {
            Ok(BackRegion::B2a)
        }
    } else {
        Ok(BackRegion::B2a) // Region 2a fall-through
    }
}

/// Upstream `BackwardOutputHS`: p(h,s) or T(h,s).
pub(crate) fn backward_output_hs(outkey: Prop, h: f64, s: f64) -> Result<f64, Error> {
    if outkey != Prop::P && outkey != Prop::T {
        return Err(Error::Input(
            "Backward HS Formulas output Temperature or Pressure only.".into(),
        ));
    }

    let region = region_determination_hs(h, s)?;
    let mut pval = 0.0;
    let mut tval = 0.0;
    match region {
        BackRegion::B1 => pval = B1HS.p_hs(h, s),
        BackRegion::B2a => pval = B2AHS.p_hs(h, s),
        BackRegion::B2b => pval = B2BHS.p_hs(h, s),
        BackRegion::B2c => pval = B2CHS.p_hs(h, s),
        BackRegion::B3a => pval = B3AHS.p_hs(h, s),
        BackRegion::B3b => pval = B3BHS.p_hs(h, s),
        BackRegion::B4 => {
            // T(h,s) only defined over part of the 2-phase region
            if s >= SGT23 {
                tval = B4HS.t_hs(h, s);
            } else {
                return Err(Error::OutOfRange("Entropy out of range".into()));
            }
        }
    }
    if outkey == Prop::P {
        if region == BackRegion::B4 {
            crate::psat97(tval)
        } else {
            Ok(pval)
        }
    } else if region == BackRegion::B4 {
        Ok(tval)
    } else {
        region_output_backward(pval, h, Prop::Hmass, false, SatState::None)
    }
}
