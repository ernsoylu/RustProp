//! IAPWS-06 ice-Ih Gibbs formulation (port of upstream `src/Ice.cpp`) —
//! self-contained; consumed by the humid-air enhancement factor and the
//! wet-bulb energy balance below the triple point.

/// Minimal complex arithmetic for the two-term IAPWS-06 summation (Rust has
/// no std complex type; only the operations Ice.cpp uses are implemented).
#[derive(Clone, Copy)]
struct C {
    re: f64,
    im: f64,
}

impl C {
    const fn new(re: f64, im: f64) -> C {
        C { re, im }
    }
    fn add(self, o: C) -> C {
        C::new(self.re + o.re, self.im + o.im)
    }
    fn sub(self, o: C) -> C {
        C::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: C) -> C {
        C::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn scale(self, s: f64) -> C {
        C::new(self.re * s, self.im * s)
    }
    /// Principal-branch complex logarithm (matches `std::log`).
    fn ln(self) -> C {
        C::new(self.re.hypot(self.im).ln(), self.im.atan2(self.re))
    }
    /// Complex division.
    fn div(self, o: C) -> C {
        let d = o.re * o.re + o.im * o.im;
        C::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

const T1: C = C::new(0.368017112855051e-1, 0.510878114959572e-1);
const R1: C = C::new(0.447050716285388e2, 0.656876847463481e2);
const T2: C = C::new(0.337315741065416, 0.335449415919309);
const R20: C = C::new(-0.725974574329220e2, -0.781008427112870e2);
const R21: C = C::new(-0.557107698030123e-4, 0.464578634580806e-4);
const R22: C = C::new(0.234801409215913e-10, -0.285651142904972e-10);

const T_T: f64 = 273.16;
const P_T: f64 = 611.657;
const P_0: f64 = 101325.0;

const G00: f64 = -0.632020233449497e6;
const G01: f64 = 0.655022213658955;
const G02: f64 = -0.189369929326131e-7;
const G03: f64 = 0.339746123271053e-14;
const G04: f64 = -0.556464869058991e-21;
const S0: f64 = -0.332733756492168e4;

/// The shared bracket `(t-θ)ln(t-θ) + (t+θ)ln(t+θ) - 2t ln t - θ²/t`.
fn g_bracket(t: C, theta: f64) -> C {
    let th = C::new(theta, 0.0);
    t.sub(th)
        .mul(t.sub(th).ln())
        .add(t.add(th).mul(t.add(th).ln()))
        .sub(t.mul(t.ln()).scale(2.0))
        .add(C::new(-theta * theta, 0.0).div(t))
}

/// The T-derivative bracket `-ln(t-θ) + ln(t+θ) - 2θ/t`.
fn dg_bracket(t: C, theta: f64) -> C {
    let th = C::new(theta, 0.0);
    t.add(th)
        .ln()
        .sub(t.sub(th).ln())
        .add(C::new(-2.0 * theta, 0.0).div(t))
}

fn r2_of(pi: f64, pi_0: f64) -> C {
    R20.add(R21.scale(pi - pi_0))
        .add(R22.scale((pi - pi_0) * (pi - pi_0)))
}

/// Sublimation pressure [Pa] — the only saturation correlation the humid-air
/// module uses below 273.16 K.
pub fn psub_ice(t: f64) -> f64 {
    let a = [0.0, -0.212144006e2, 0.273203819e2, -0.610598130e1];
    let b = [0.0, 0.333333333e-2, 0.120666667e1, 0.170333333e1];
    let theta = t / T_T;
    let mut summer = 0.0;
    for i in 1..=3 {
        summer += a[i] * theta.powf(b[i]);
    }
    P_T * (summer / theta).exp()
}

/// Gibbs energy [J/kg].
pub fn g_ice(t: f64, p: f64) -> f64 {
    let theta = t / T_T;
    let pi = p / P_T;
    let pi_0 = P_0 / P_T;
    let d = pi - pi_0;
    let g0 = G00 + G01 * d + G02 * d * d + G03 * d * d * d + G04 * d * d * d * d;
    let r2 = r2_of(pi, pi_0);
    let term1 = R1.mul(g_bracket(T1, theta));
    let term2 = r2.mul(g_bracket(T2, theta));
    g0 - S0 * T_T * theta + T_T * term1.add(term2).re
}

/// `dg/dp` [m^3/kg].
pub fn dg_dp_ice(t: f64, p: f64) -> f64 {
    let theta = t / T_T;
    let pi = p / P_T;
    let pi_0 = P_0 / P_T;
    let d = pi - pi_0;
    let g0_p =
        G01 / P_T + G02 * 2.0 / P_T * d + G03 * 3.0 / P_T * d * d + G04 * 4.0 / P_T * d * d * d;
    let r2_p = R21.scale(1.0 / P_T).add(R22.scale(2.0 / P_T * d));
    g0_p + T_T * r2_p.mul(g_bracket(T2, theta)).re
}

/// `d²g/dp²`.
pub fn dg2_dp2_ice(t: f64, p: f64) -> f64 {
    let theta = t / T_T;
    let pi = p / P_T;
    let pi_0 = P_0 / P_T;
    let d = pi - pi_0;
    let g0_pp = G02 * 2.0 / P_T / P_T + G03 * 6.0 / P_T / P_T * d + G04 * 12.0 / P_T / P_T * d * d;
    let r2_pp = R22.scale(2.0 / P_T / P_T);
    g0_pp + T_T * r2_pp.mul(g_bracket(T2, theta)).re
}

/// `dg/dT` [J/kg/K].
pub fn dg_dt_ice(t: f64, p: f64) -> f64 {
    let theta = t / T_T;
    let pi = p / P_T;
    let pi_0 = P_0 / P_T;
    let r2 = r2_of(pi, pi_0);
    let term1 = R1.mul(dg_bracket(T1, theta));
    let term2 = r2.mul(dg_bracket(T2, theta));
    -S0 + term1.add(term2).re
}

/// Enthalpy [J/kg].
pub fn h_ice(t: f64, p: f64) -> f64 {
    g_ice(t, p) - t * dg_dt_ice(t, p)
}

/// Isothermal compressibility [1/Pa].
pub fn isotherm_compress_ice(t: f64, p: f64) -> f64 {
    -dg2_dp2_ice(t, p) / dg_dp_ice(t, p)
}
