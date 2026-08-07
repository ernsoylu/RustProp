//! Incompressible engine — brines and secondary working fluids (PLAN.md
//! Phase 8), a port of CoolProp 8's
//! `src/Backends/Incompressible/{IncompressibleFluid,IncompressibleBackend}.cpp`
//! and the `Polynomial2DFrac` machinery they consume.
//!
//! The thermodynamic model in full: `rho`/`c`/`conductivity`/`viscosity`/
//! `psat`/`Tfreeze` are direct evaluations of the document's coefficient
//! blocks in `(T - Tbase, x - xbase)` powers;
//! `h_raw = ∫c dT + p·(1/rho)(1 + (T/rho)·drho/dT)` and
//! `s_raw = ∫(c/T) dT + p·(1/rho^2)·drho/dT`, referenced against the
//! HARD-CODED state (T0 = 293.15 K, p0 = 101325 Pa, h0 = s0 = 0 — upstream's
//! per-fluid reference reader is commented out); `u = h - p/rho`;
//! `cp = cv = c`. Only five input pairs exist: PT, DmassP, PSmass, HmassP,
//! and QT with Q exactly 0 (saturated liquid).

use rustprop_core::fluid::{IncompData, IncompFluid};
use rustprop_core::{Error, Result};
use rustprop_heos::solvers::brent;

/// Upstream `INCOMP_EPSILON` / `INCOMP_DELTA` (IncompressibleFluid.cpp).
const INCOMP_EPSILON: f64 = f64::EPSILON * 100.0;
const INCOMP_DELTA: f64 = INCOMP_EPSILON * 10.0;

/// Upstream's hard-coded reference state (`set_reference_state` defaults).
const T_REF: f64 = 20.0 + 273.15;
const P_REF: f64 = 101325.0;

/// 2-D polynomial evaluation: `sum_j (T-Tb)^j * sum_i C[j][i] (x-xb)^i` —
/// upstream `Polynomial2DFrac::evaluate` (Horner over rows, then cols).
fn poly2d(coeffs: &[&[f64]], t: f64, x: f64, tbase: f64, xbase: f64) -> f64 {
    // Upstream initializes Horner FROM the highest row/coefficient (no
    // multiply on the first step) — with a single row the T argument is
    // never touched, which the trivial T_freeze path relies on (its p is
    // the cleared -HUGE sentinel; 0.0 * -inf would be NaN).
    let row = |r: &[f64]| -> f64 {
        let mut it = r.iter().rev();
        let mut acc = *it.next().unwrap_or(&0.0);
        for c in it {
            acc = acc * (x - xbase) + c;
        }
        acc
    };
    let mut it = coeffs.iter().rev();
    let mut acc = it.next().map_or(0.0, |r| row(r));
    for r in it {
        acc = acc * (t - tbase) + row(r);
    }
    acc
}

/// `d/dT` of the 2-D polynomial (upstream `derivative`, axis 0).
fn poly2d_dt(coeffs: &[&[f64]], t: f64, x: f64, tbase: f64, xbase: f64) -> f64 {
    let row = |r: &[f64]| -> f64 {
        let mut acc = 0.0;
        for c in r.iter().rev() {
            acc = acc * (x - xbase) + c;
        }
        acc
    };
    let mut acc = 0.0;
    for (j, r) in coeffs.iter().enumerate().skip(1).rev() {
        acc = acc * (t - tbase) + (j as f64) * row(r);
    }
    acc
}

/// `∫ f dT` of the 2-D polynomial (upstream `integral`, axis 0, x_exp 0):
/// `sum_j a_j(x) (T-Tb)^(j+1)/(j+1)`.
fn poly2d_int_dt(coeffs: &[&[f64]], t: f64, x: f64, tbase: f64, xbase: f64) -> f64 {
    let row = |r: &[f64]| -> f64 {
        let mut acc = 0.0;
        for c in r.iter().rev() {
            acc = acc * (x - xbase) + c;
        }
        acc
    };
    let mut acc = 0.0;
    for (j, r) in coeffs.iter().enumerate().rev() {
        acc = acc * (t - tbase) + row(r) / (j as f64 + 1.0);
    }
    acc * (t - tbase)
}

/// Binomial coefficient (upstream `binom` via factorials; the arguments
/// stay tiny — j < 6).
fn binom(n: usize, k: usize) -> f64 {
    let mut r = 1.0;
    for i in 0..k {
        r = r * ((n - i) as f64) / ((i + 1) as f64);
    }
    r
}

/// `∫ f/T dT` (upstream `integral`, axis 0, x_exp = -1): the `Tbase == 0`
/// closed form or `fracIntCentral` for centered fits.
fn poly2d_fracint_dt(coeffs: &[&[f64]], t: f64, x: f64, tbase: f64, xbase: f64) -> f64 {
    let row = |r: &[f64]| -> f64 {
        let mut acc = 0.0;
        for c in r.iter().rev() {
            acc = acc * (x - xbase) + c;
        }
        acc
    };
    if tbase.abs() < f64::EPSILON {
        // a_0 ln(T) + sum_{j>=1} a_j T^j / j
        let mut result = row(coeffs[0]) * t.ln();
        let mut acc = 0.0;
        for (j, r) in coeffs.iter().enumerate().skip(1).rev() {
            acc = acc * t + row(r) / (j as f64);
        }
        result += acc * t;
        result
    } else {
        // fracIntCentral: sum_j a_j * D_j(T, Tbase) with
        // D_j = (-1)^j ln(T) Tbase^j + sum_{k<j} binom(j,k)(-1)^k/(j-k) T^(j-k) Tbase^k
        let mut result = 0.0;
        for (j, r) in coeffs.iter().enumerate() {
            let mut d = (-1.0f64).powi(j as i32) * t.ln() * tbase.powi(j as i32);
            for k in 0..j {
                d += binom(j, k) * (-1.0f64).powi(k as i32) / ((j - k) as f64)
                    * t.powi((j - k) as i32)
                    * tbase.powi(k as i32);
            }
            result += row(r) * d;
        }
        result
    }
}

/// Upstream `baseExponential`: `exp(C0/((y-ybase)+C1) - C2)` with the
/// pole-linearization guard.
fn base_exponential(c: &[f64], y: f64, ybase: f64) -> f64 {
    let fnc = |x: f64| (c[0] / x - c[2]).exp();
    let x_den = (y - ybase) + c[1];
    if !(-INCOMP_EPSILON..=INCOMP_EPSILON).contains(&x_den) {
        return fnc(x_den);
    }
    let x_lo = -INCOMP_EPSILON - INCOMP_DELTA;
    let x_hi = INCOMP_EPSILON + INCOMP_DELTA;
    let f_lo = fnc(x_lo);
    let f_hi = fnc(x_hi);
    (f_hi - f_lo) / (x_hi - x_lo) * (x_den - x_lo) + f_lo
}

/// Upstream `baseLogexponential`: `exp(ln(1/u + 1/u^2)·C1 + C2)`,
/// `u = (y-ybase)+C0`, same guard.
fn base_logexponential(c: &[f64], y: f64, ybase: f64) -> f64 {
    let fnc = |x: f64| ((1.0 / x + 1.0 / (x * x)).ln() * c[1] + c[2]).exp();
    let x_den = (y - ybase) + c[0];
    if !(-INCOMP_EPSILON..=INCOMP_EPSILON).contains(&x_den) {
        return fnc(x_den);
    }
    let x_lo = -INCOMP_EPSILON - INCOMP_DELTA;
    let x_hi = INCOMP_EPSILON + INCOMP_DELTA;
    let f_lo = fnc(x_lo);
    let f_hi = fnc(x_hi);
    (f_hi - f_lo) / (x_hi - x_lo) * (x_den - x_lo) + f_lo
}

/// Upstream `basePolyOffset` for the JSON-loaded (column-vector) case:
/// `in = y`, centering value = coeffs[0], polynomial = coeffs[1..]. The
/// only data user is ExampleSecCool's `T_freeze`, whose `y` is the
/// PRESSURE — upstream's live wiring, reproduced (see the survey note).
fn base_polyoffset(c: &[f64], y: f64) -> f64 {
    let offset = c[0];
    let mut acc = 0.0;
    for coef in c[1..].iter().rev() {
        acc = acc * (y - offset) + coef;
    }
    acc
}

/// One property evaluation through the block's type dispatch (upstream
/// `rho`/`visc`/`cond`/`psat` share this shape; `ctx` names the property
/// for the error strings).
fn eval_block(
    block: &IncompData,
    t: f64,
    x: f64,
    tbase: f64,
    xbase: f64,
    ctx: &str,
) -> Result<f64> {
    Ok(match block {
        IncompData::Polynomial(c) => poly2d(c, t, x, tbase, xbase),
        IncompData::Exponential(c) => base_exponential(c, t, 0.0),
        IncompData::LogExponential(c) => base_logexponential(c, t, 0.0),
        IncompData::ExpPolynomial(c) => poly2d(c, t, x, tbase, xbase).exp(),
        IncompData::PolyOffset(c) => base_polyoffset(c, t),
        IncompData::NotSet => {
            return Err(Error::Value(format!(
                "The function type is not specified for {ctx}, are you sure the coefficients have been set?"
            )));
        }
    })
}

/// One incompressible fluid at one mass/volume fraction.
pub struct IncompEos {
    pub fluid: &'static IncompFluid,
}

impl IncompEos {
    pub fn new(fluid: &'static IncompFluid) -> Self {
        IncompEos { fluid }
    }

    /// Density [kg/m^3].
    pub fn rho(&self, t: f64, x: f64) -> Result<f64> {
        eval_block(
            &self.fluid.density,
            t,
            x,
            self.fluid.tbase,
            self.fluid.xbase,
            "density",
        )
    }

    /// Specific heat [J/kg/K] — POLYNOMIAL only upstream.
    pub fn c(&self, t: f64, x: f64) -> Result<f64> {
        match &self.fluid.specific_heat {
            IncompData::Polynomial(c) => {
                Ok(poly2d(c, t, x, self.fluid.tbase, self.fluid.xbase))
            }
            IncompData::NotSet => Err(Error::Value(
                "The function type is not specified for specific heat, are you sure the coefficients have been set?"
                    .into(),
            )),
            _ => Err(Error::Value(
                "There is no predefined way to use this function type for specific heat."
                    .into(),
            )),
        }
    }

    /// Dynamic viscosity [Pa-s].
    pub fn visc(&self, t: f64, x: f64) -> Result<f64> {
        eval_block(
            &self.fluid.viscosity,
            t,
            x,
            self.fluid.tbase,
            self.fluid.xbase,
            "viscosity",
        )
    }

    /// Thermal conductivity [W/m/K].
    pub fn cond(&self, t: f64, x: f64) -> Result<f64> {
        eval_block(
            &self.fluid.conductivity,
            t,
            x,
            self.fluid.tbase,
            self.fluid.xbase,
            "conductivity",
        )
    }

    /// Saturation pressure [Pa], with the v8 `TminPsat` gate (fires before
    /// the type dispatch).
    pub fn psat(&self, t: f64, x: f64) -> Result<f64> {
        if t <= self.fluid.tmin_psat {
            return Err(Error::Value(format!(
                "Saturation pressure is not available below TminPsat={} K (T={} K)",
                self.fluid.tmin_psat, t
            )));
        }
        eval_block(
            &self.fluid.saturation_pressure,
            t,
            x,
            self.fluid.tbase,
            self.fluid.xbase,
            "saturation pressure",
        )
    }

    /// Freezing temperature [K] — note upstream's per-type argument wiring:
    /// poly/exppoly see `(p, x)` with `Tbase = 0`; exponential forms see
    /// `x`; polyoffset sees `p` (the JSON column-vector quirk).
    pub fn t_freeze(&self, p: f64, x: f64) -> Result<f64> {
        let xb = self.fluid.xbase;
        Ok(match &self.fluid.t_freeze {
            IncompData::Polynomial(c) => poly2d(c, p, x, 0.0, xb),
            IncompData::ExpPolynomial(c) => poly2d(c, p, x, 0.0, xb).exp(),
            IncompData::Exponential(c) => base_exponential(c, x, 0.0),
            IncompData::LogExponential(c) => base_logexponential(c, x, 0.0),
            IncompData::PolyOffset(c) => base_polyoffset(c, p),
            IncompData::NotSet => {
                return Err(Error::Value(
                    "The function type is not specified for T_freeze, are you sure the coefficients have been set?"
                        .into(),
                ));
            }
        })
    }

    /// `drho/dT|p,x` — POLYNOMIAL only upstream.
    fn drhodt(&self, t: f64, x: f64) -> Result<f64> {
        match &self.fluid.density {
            IncompData::Polynomial(c) => Ok(poly2d_dt(c, t, x, self.fluid.tbase, self.fluid.xbase)),
            _ => Err(Error::Value(
                "There is no predefined way to use this function type for density.".into(),
            )),
        }
    }

    /// `∫c dT` and `∫c/T dT` — POLYNOMIAL only upstream.
    fn int_c_dt(&self, t: f64, x: f64) -> Result<f64> {
        match &self.fluid.specific_heat {
            IncompData::Polynomial(c) => {
                Ok(poly2d_int_dt(c, t, x, self.fluid.tbase, self.fluid.xbase))
            }
            _ => Err(Error::Value(
                "There is no predefined way to use this function type for entropy.".into(),
            )),
        }
    }
    fn fracint_c_dt(&self, t: f64, x: f64) -> Result<f64> {
        match &self.fluid.specific_heat {
            IncompData::Polynomial(c) => Ok(poly2d_fracint_dt(
                c,
                t,
                x,
                self.fluid.tbase,
                self.fluid.xbase,
            )),
            _ => Err(Error::Value(
                "There is no predefined way to use this function type for entropy.".into(),
            )),
        }
    }

    /// `h_raw = ∫c dT + p·(1/rho)(1 + (T/rho) drho/dT)`.
    pub fn raw_hmass(&self, t: f64, p: f64, x: f64) -> Result<f64> {
        let rho = self.rho(t, x)?;
        let drho = self.drhodt(t, x)?;
        Ok(self.int_c_dt(t, x)? + p * (1.0 / rho) * (1.0 + t / rho * drho))
    }

    /// `s_raw = ∫c/T dT + p·(1/rho^2) drho/dT`.
    pub fn raw_smass(&self, t: f64, p: f64, x: f64) -> Result<f64> {
        let rho = self.rho(t, x)?;
        let drho = self.drhodt(t, x)?;
        Ok(self.fracint_c_dt(t, x)? + p * (1.0 / (rho * rho)) * drho)
    }

    /// Reference-anchored enthalpy [J/kg] (h0 = 0 at 293.15 K, 101325 Pa).
    pub fn hmass(&self, t: f64, p: f64, x: f64) -> Result<f64> {
        Ok(self.raw_hmass(t, p, x)? - self.raw_hmass(T_REF, P_REF, x)?)
    }

    /// Reference-anchored entropy [J/kg/K].
    pub fn smass(&self, t: f64, p: f64, x: f64) -> Result<f64> {
        Ok(self.raw_smass(t, p, x)? - self.raw_smass(T_REF, P_REF, x)?)
    }

    /// Internal energy `u = h - p/rho` [J/kg].
    pub fn umass(&self, t: f64, p: f64, x: f64) -> Result<f64> {
        Ok(self.hmass(t, p, x)? - p / self.rho(t, x)?)
    }

    /// Upstream `checkTPX` — T, then P, then X, short-circuiting with
    /// verbatim messages.
    pub fn check_tpx(&self, t: f64, p: f64, x: f64) -> Result<()> {
        let fl = self.fluid;
        if fl.tmin <= 0.0 {
            return Err(Error::Value(
                "Please specify the minimum temperature.".into(),
            ));
        }
        if fl.tmax <= 0.0 {
            return Err(Error::Value(
                "Please specify the maximum temperature.".into(),
            ));
        }
        if fl.tmin > t || t > fl.tmax {
            return Err(Error::Value(format!(
                "Your temperature {t:.6} is not between {:.6} and {:.6}.",
                fl.tmin, fl.tmax
            )));
        }
        let mut tf = 0.0;
        if !matches!(fl.t_freeze, IncompData::NotSet) {
            tf = self.t_freeze(p, x)?;
        }
        if t < tf {
            return Err(Error::Value(format!(
                "Your temperature {t:.6} is below the freezing point of {tf:.6}."
            )));
        }
        let mut ps = 0.0;
        if !matches!(fl.saturation_pressure, IncompData::NotSet) {
            ps = self.psat(t, x).unwrap_or(0.0);
        }
        if p < 0.0 {
            return Err(Error::Value(format!(
                "You cannot use negative pressures: {p:.6} < {:.6}. ",
                0.0
            )));
        }
        if ps > 0.0 && p < ps {
            return Err(Error::Value(format!(
                "Equations are valid for liquid phase only: {p:.6} < {ps:.6} (psat). "
            )));
        }
        self.check_x(x)
    }

    /// Upstream `checkX`.
    pub fn check_x(&self, x: f64) -> Result<()> {
        let fl = self.fluid;
        if !(0.0..=1.0).contains(&fl.xmin) {
            return Err(Error::Value(
                "Please specify the minimum concentration between 0 and 1.".into(),
            ));
        }
        if !(0.0..=1.0).contains(&fl.xmax) {
            return Err(Error::Value(
                "Please specify the maximum concentration between 0 and 1.".into(),
            ));
        }
        if x < fl.xmin * (1.0 - INCOMP_EPSILON) || x > fl.xmax * (1.0 + INCOMP_EPSILON) {
            return Err(Error::Value(format!(
                "Your composition {x} is not between {} and {}.",
                fl.xmin, fl.xmax
            )));
        }
        Ok(())
    }

    /// The three T-solving flashes: Brent on [Tmin, Tmax], upstream's
    /// `macheps = DBL_EPSILON`, `tol = DBL_EPSILON*1e3`, `maxiter = 10`.
    pub fn t_from_rho(&self, rho: f64, x: f64) -> Result<f64> {
        let f = |t: f64| self.rho(t, x).map_or(f64::NAN, |r| r - rho);
        brent(
            f,
            self.fluid.tmin,
            self.fluid.tmax,
            f64::EPSILON,
            f64::EPSILON * 1e3,
            10,
        )
    }
    pub fn t_from_hmass(&self, h: f64, p: f64, x: f64) -> Result<f64> {
        let target = h + self.raw_hmass(T_REF, P_REF, x)?;
        let f = |t: f64| self.raw_hmass(t, p, x).map_or(f64::NAN, |v| v - target);
        brent(
            f,
            self.fluid.tmin,
            self.fluid.tmax,
            f64::EPSILON,
            f64::EPSILON * 1e3,
            10,
        )
    }
    pub fn t_from_smass(&self, s: f64, p: f64, x: f64) -> Result<f64> {
        let target = s + self.raw_smass(T_REF, P_REF, x)?;
        let f = |t: f64| self.raw_smass(t, p, x).map_or(f64::NAN, |v| v - target);
        brent(
            f,
            self.fluid.tmin,
            self.fluid.tmax,
            f64::EPSILON,
            f64::EPSILON * 1e3,
            10,
        )
    }
}

pub use rustprop_core::UPSTREAM_VERSION;
