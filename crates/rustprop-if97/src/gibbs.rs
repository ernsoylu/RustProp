//! Gibbs-form regions 1, 2, and 5 — port of `BaseRegion`/`Region1`/`Region2`/
//! `Region5` from IF97.h. Region-specific overrides (region 1 cv, speed of
//! sound, drhodp; region 5 conductivity without critical enhancement) are
//! dispatched on `Kind`, mirroring the upstream virtual methods.

use crate::tables::{
    Ideal, REGION1_RESID, REGION2_IDEAL, REGION2_RESID, REGION5_IDEAL, REGION5_RESID, Resid,
};
use crate::transport;
use crate::{P_FACT, Prop, R_FACT, RGAS, powi};
use rustprop_core::Error;

#[derive(Clone, Copy, PartialEq)]
pub(crate) enum Kind {
    Region1,
    Region2,
    Region5,
}

pub(crate) struct GibbsRegion {
    pub kind: Kind,
    t_star: f64,
    p_star: f64,
    /// `TAUrterm(T) = t_star/T - tau_shift`
    tau_shift: f64,
    /// `PIrterm(p) = p/p_star - pi_shift`
    pi_shift: f64,
    resid: &'static [Resid],
    /// Empty for region 1 (no ideal-gas part in its Gibbs formulation).
    ideal: &'static [Ideal],
}

pub(crate) const REGION1: GibbsRegion = GibbsRegion {
    kind: Kind::Region1,
    t_star: 1386.0,
    p_star: 16.53 * P_FACT,
    tau_shift: 1.222,
    pi_shift: 7.1,
    resid: REGION1_RESID,
    ideal: &[],
};

pub(crate) const REGION2: GibbsRegion = GibbsRegion {
    kind: Kind::Region2,
    t_star: 540.0,
    p_star: 1.0 * P_FACT,
    tau_shift: 0.5,
    pi_shift: 0.0,
    resid: REGION2_RESID,
    ideal: REGION2_IDEAL,
};

pub(crate) const REGION5: GibbsRegion = GibbsRegion {
    kind: Kind::Region5,
    t_star: 1000.0,
    p_star: 1.0 * P_FACT,
    tau_shift: 0.0,
    pi_shift: 0.0,
    resid: REGION5_RESID,
    ideal: REGION5_IDEAL,
};

impl GibbsRegion {
    fn pi_r(&self, p: f64) -> f64 {
        p / self.p_star - self.pi_shift
    }
    fn tau_r(&self, t: f64) -> f64 {
        self.t_star / t - self.tau_shift
    }
    /// `TAU0term` — only meaningful when the ideal part exists.
    fn tau_0(&self, t: f64) -> f64 {
        self.t_star / t
    }

    fn gammar(&self, t: f64, p: f64) -> f64 {
        let (pi, tau) = (self.pi_r(p), self.tau_r(t));
        self.resid
            .iter()
            .fold(0.0, |s, e| s + e.n * powi(pi, e.i) * powi(tau, e.j))
    }
    fn dgammar_dpi(&self, t: f64, p: f64) -> f64 {
        let (pi, tau) = (self.pi_r(p), self.tau_r(t));
        self.resid.iter().fold(0.0, |s, e| {
            s + e.n * f64::from(e.i) * powi(pi, e.i - 1) * powi(tau, e.j)
        })
    }
    fn d2gammar_dpi2(&self, t: f64, p: f64) -> f64 {
        let (pi, tau) = (self.pi_r(p), self.tau_r(t));
        self.resid.iter().fold(0.0, |s, e| {
            s + e.n * f64::from(e.i) * f64::from(e.i - 1) * powi(pi, e.i - 2) * powi(tau, e.j)
        })
    }
    fn dgammar_dtau(&self, t: f64, p: f64) -> f64 {
        let (pi, tau) = (self.pi_r(p), self.tau_r(t));
        self.resid.iter().fold(0.0, |s, e| {
            s + e.n * f64::from(e.j) * powi(pi, e.i) * powi(tau, e.j - 1)
        })
    }
    fn d2gammar_dpidtau(&self, t: f64, p: f64) -> f64 {
        let (pi, tau) = (self.pi_r(p), self.tau_r(t));
        self.resid.iter().fold(0.0, |s, e| {
            s + e.n * f64::from(e.j) * f64::from(e.i) * powi(pi, e.i - 1) * powi(tau, e.j - 1)
        })
    }
    fn d2gammar_dtau2(&self, t: f64, p: f64) -> f64 {
        let (pi, tau) = (self.pi_r(p), self.tau_r(t));
        self.resid.iter().fold(0.0, |s, e| {
            s + e.n * f64::from(e.j) * f64::from(e.j - 1) * powi(pi, e.i) * powi(tau, e.j - 2)
        })
    }

    fn gamma0(&self, t: f64, p: f64) -> f64 {
        if self.ideal.is_empty() {
            return 0.0; // Region 1 has no term
        }
        let pi = p / self.p_star;
        let tau = self.tau_0(t);
        self.ideal
            .iter()
            .fold(pi.ln(), |s, e| s + e.n * powi(tau, e.j))
    }
    fn dgamma0_dpi(&self, p: f64) -> f64 {
        if self.ideal.is_empty() {
            return 0.0; // Region 1 has no term
        }
        let pi = p / self.p_star;
        1.0 / pi
    }
    fn dgamma0_dtau(&self, t: f64) -> f64 {
        let tau = self.tau_0(t);
        self.ideal
            .iter()
            .fold(0.0, |s, e| s + e.n * f64::from(e.j) * powi(tau, e.j - 1))
    }
    fn d2gamma0_dtau2(&self, t: f64) -> f64 {
        let tau = self.tau_0(t);
        self.ideal.iter().fold(0.0, |s, e| {
            s + e.n * f64::from(e.j) * f64::from(e.j - 1) * powi(tau, e.j - 2)
        })
    }

    pub fn rhomass(&self, t: f64, p: f64) -> f64 {
        self.p_star
            / (RGAS * t)
            / (P_FACT / 1000.0 / R_FACT)
            / (self.dgamma0_dpi(p) + self.dgammar_dpi(t, p))
    }
    pub fn hmass(&self, t: f64, p: f64) -> f64 {
        RGAS * self.t_star * (self.dgamma0_dtau(t) + self.dgammar_dtau(t, p))
    }
    pub fn smass(&self, t: f64, p: f64) -> f64 {
        let tau = self.t_star / t;
        RGAS * (tau * (self.dgamma0_dtau(t) + self.dgammar_dtau(t, p))
            - (self.gammar(t, p) + self.gamma0(t, p)))
    }
    pub fn umass(&self, t: f64, p: f64) -> f64 {
        let (tau, pi) = (self.t_star / t, p / self.p_star);
        RGAS * t
            * (tau * (self.dgamma0_dtau(t) + self.dgammar_dtau(t, p))
                - pi * (self.dgamma0_dpi(p) + self.dgammar_dpi(t, p)))
    }
    pub fn cpmass(&self, t: f64, p: f64) -> f64 {
        let tau = self.t_star / t;
        -RGAS * tau * tau * (self.d2gammar_dtau2(t, p) + self.d2gamma0_dtau2(t))
    }
    pub fn cvmass(&self, t: f64, p: f64) -> f64 {
        let tau = self.t_star / t;
        if self.kind == Kind::Region1 {
            // Region 1 override (Table 3 of R7-97)
            return RGAS
                * (-tau * tau * self.d2gammar_dtau2(t, p)
                    + powi(
                        self.dgammar_dpi(t, p) - tau * self.d2gammar_dpidtau(t, p),
                        2,
                    ) / self.d2gammar_dpi2(t, p));
        }
        let pi = p / self.p_star;
        self.cpmass(t, p)
            - RGAS
                * powi(
                    1.0 + pi * self.dgammar_dpi(t, p) - tau * pi * self.d2gammar_dpidtau(t, p),
                    2,
                )
                / (1.0 - pi * pi * self.d2gammar_dpi2(t, p))
    }
    pub fn speed_sound(&self, t: f64, p: f64) -> f64 {
        let tau = self.t_star / t;
        if self.kind == Kind::Region1 {
            // Region 1 override (Table 3 of R7-97)
            let rhs = powi(self.dgammar_dpi(t, p), 2)
                / (powi(
                    self.dgammar_dpi(t, p) - tau * self.d2gammar_dpidtau(t, p),
                    2,
                ) / (tau * tau * self.d2gammar_dtau2(t, p))
                    - self.d2gammar_dpi2(t, p));
            return (RGAS * (1000.0 / R_FACT) * t * rhs).sqrt();
        }
        let pi = p / self.p_star;
        let rhs =
            (1.0 + 2.0 * pi * self.dgammar_dpi(t, p) + pi * pi * powi(self.dgammar_dpi(t, p), 2))
                / ((1.0 - pi * pi * self.d2gammar_dpi2(t, p))
                    + powi(
                        1.0 + pi * self.dgammar_dpi(t, p) - tau * pi * self.d2gammar_dpidtau(t, p),
                        2,
                    ) / (tau * tau * (self.d2gamma0_dtau2(t) + self.d2gammar_dtau2(t, p))));
        (RGAS * (1000.0 / R_FACT) * t * rhs).sqrt()
    }
    pub fn drhodp(&self, t: f64, p: f64) -> f64 {
        if self.kind == Kind::Region1 {
            // Region 1 override, from IAPWS Revised Advisory Note No. 3
            return -self.d2gammar_dpi2(t, p) / (powi(self.dgammar_dpi(t, p), 2) * RGAS * t)
                * (1000.0 * R_FACT / P_FACT);
        }
        // Regions 2 and 5, from IAPWS Revised Advisory Note No. 3
        let pi = p / self.p_star;
        (self.rhomass(t, p) / p)
            * ((1.0 - pi * pi * self.d2gammar_dpi2(t, p)) / (1.0 + pi * self.dgammar_dpi(t, p)))
    }

    fn tcond(&self, t: f64, p: f64, rho: f64) -> f64 {
        let lambda2 = if self.kind == Kind::Region5 {
            0.0 // No critical enhancement of thermal conductivity in Region 5
        } else {
            transport::lambda2_gibbs(
                t,
                rho,
                self.cpmass(t, p),
                self.cvmass(t, p),
                transport::visc(t, rho),
                self.drhodp(t, p),
            )
        };
        0.001 * (transport::lambda0(t) * transport::lambda1(t, rho) + lambda2)
    }

    pub fn output(&self, key: Prop, t: f64, p: f64) -> Result<f64, Error> {
        match key {
            Prop::T => Ok(t),
            Prop::P => Ok(p),
            Prop::Dmass => Ok(self.rhomass(t, p)),
            Prop::Hmass => Ok(self.hmass(t, p)),
            Prop::Smass => Ok(self.smass(t, p)),
            Prop::Umass => Ok(self.umass(t, p)),
            Prop::Cpmass => Ok(self.cpmass(t, p)),
            Prop::Cvmass => Ok(self.cvmass(t, p)),
            Prop::W => Ok(self.speed_sound(t, p)),
            Prop::Mu => Ok(transport::visc(t, self.rhomass(t, p))),
            Prop::K => Ok(self.tcond(t, p, self.rhomass(t, p))),
            Prop::DrhoDp => Ok(self.drhodp(t, p)),
            Prop::Q => Err(Error::Input("Can't determine Q from T & P".into())),
        }
    }
}
