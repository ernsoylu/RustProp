//! Fluid-data *types* mirroring the CoolProp v8.0.0 runtime JSON documents
//! (PLAN.md 3.1). Data *contents* are generated into `rustprop-data` by
//! `tools/rustprop-datagen` — this crate never holds values.
//!
//! Scope: exactly the structures Phase 4 (HEOS pure fluids) consumes, mapped
//! from `data/coolprop-json/Water.json`. Every field documents the JSON key it
//! mirrors. Known present-but-not-yet-ported document parts (SUPERANCILLARY,
//! critical_region_splines, hL/hLV/sL/sLV ancillaries, melting line, surface
//! tension, TRANSPORT, ENVIRONMENTAL) are listed in the data-fidelity test's
//! skip list and land with the phases that consume them.

/// One fluid document (top-level array element of the CoolProp JSON).
pub struct FluidData {
    /// `INFO.NAME`
    pub name: &'static str,
    /// `INFO.CAS`
    pub cas: &'static str,
    /// `INFO.ALIASES`
    pub aliases: &'static [&'static str],
    /// `EOS[0]` (CoolProp carries exactly one EOS per fluid document)
    pub eos: Eos,
    /// `ANCILLARIES`
    pub ancillaries: Ancillaries,
    /// `STATES`
    pub states: States,
}

/// `EOS[0]` — multiparameter Helmholtz EOS definition.
pub struct Eos {
    /// `EOS[0].gas_constant` [J/mol/K] — fluid-specific, not CODATA
    pub gas_constant: f64,
    /// `EOS[0].molar_mass` [kg/mol]
    pub molar_mass: f64,
    /// `EOS[0].p_max` [Pa]
    pub p_max: f64,
    /// `EOS[0].T_max` [K]
    pub t_max: f64,
    /// `EOS[0].Ttriple` [K]
    pub t_triple: f64,
    /// `EOS[0].acentric` [-]
    pub acentric: f64,
    /// `EOS[0].pseudo_pure`
    pub pseudo_pure: bool,
    /// `EOS[0].STATES.reducing`
    pub reducing: StatePoint,
    /// `EOS[0].STATES.sat_min_liquid`
    pub sat_min_liquid: StatePoint,
    /// `EOS[0].STATES.sat_min_vapor`
    pub sat_min_vapor: StatePoint,
    /// `EOS[0].STATES.hs_anchor`
    pub hs_anchor: StatePoint,
    /// `EOS[0].alpha0` — ideal-gas Helmholtz terms, in document order
    pub alpha0: &'static [Alpha0Term],
    /// `EOS[0].alphar` — residual Helmholtz terms, in document order
    pub alphar: &'static [AlpharTerm],
    /// `EOS[0].SUPERANCILLARY` — absent for fluids without one
    pub superancillary: Option<SuperAncillaryData>,
}

/// `EOS[0].SUPERANCILLARY` — the fields upstream's loader consumes
/// (`src/superancillary.cpp`); `crit_anc` and the remaining `meta` entries are
/// fitting-time artifacts and stay skip-listed.
pub struct SuperAncillaryData {
    /// `.jexpansions_p`
    pub p: &'static [ChebyshevInterval],
    /// `.jexpansions_rhoL`
    pub rho_l: &'static [ChebyshevInterval],
    /// `.jexpansions_rhoV`
    pub rho_v: &'static [ChebyshevInterval],
    /// `.meta."Tcrittrue / K"`
    pub t_crit_num: f64,
    /// `.meta."rhocrittrue / mol/m^3"`
    pub rho_crit_num: f64,
    /// `.check_points` — extended-precision verification states
    pub check_points: &'static [SuperAncCheckPoint],
}

/// One piecewise-Chebyshev interval (`jexpansions_*` array element).
pub struct ChebyshevInterval {
    /// `.xmin`
    pub xmin: f64,
    /// `.xmax`
    pub xmax: f64,
    /// `.coef`
    pub coef: &'static [f64],
}

/// One `check_points` element (keys carry units, e.g. `"T / K"`).
pub struct SuperAncCheckPoint {
    /// `."T / K"`
    pub t: f64,
    /// `."p(mp) / Pa"` — multiprecision reference
    pub p: f64,
    /// `."rho'(mp) / mol/m^3"`
    pub rho_l: f64,
    /// `."rho''(mp) / mol/m^3"`
    pub rho_v: f64,
    /// `."p(SA)/p(mp)"` — double-precision eval over multiprecision
    pub p_ratio: f64,
    /// `."rho'(SA)/rho'(mp)"`
    pub rho_l_ratio: f64,
    /// `."rho''(SA)/rho''(mp)"`
    pub rho_v_ratio: f64,
}

/// A state point (`T`/`p`/`rhomolar`/`hmolar`/`smolar` fields of the JSON
/// state objects; molar SI units as in the document).
pub struct StatePoint {
    /// `.T` [K]
    pub t: f64,
    /// `.p` [Pa]
    pub p: f64,
    /// `.rhomolar` [mol/m^3]
    pub rhomolar: f64,
    /// `.hmolar` [J/mol]
    pub hmolar: f64,
    /// `.smolar` [J/mol/K]
    pub smolar: f64,
}

/// `EOS[0].alpha0[i]`, tagged by `type`.
pub enum Alpha0Term {
    /// `IdealGasHelmholtzLead` — `a1 + a2*tau + ln(delta)`
    Lead {
        /// `.a1`
        a1: f64,
        /// `.a2`
        a2: f64,
    },
    /// `IdealGasHelmholtzLogTau` — `a*ln(tau)`
    LogTau {
        /// `.a`
        a: f64,
    },
    /// `IdealGasHelmholtzPlanckEinstein` — `sum n_k * ln(1 - exp(-t_k*tau))`
    PlanckEinstein {
        /// `.n`
        n: &'static [f64],
        /// `.t`
        t: &'static [f64],
    },
    /// `IdealGasHelmholtzPlanckEinsteinFunctionT` — Planck-Einstein terms
    /// with frequencies in temperature units; upstream maps
    /// `theta = -v/Tcrit` and merges into the generalized Planck-Einstein
    /// container at parse time.
    PlanckEinsteinFunctionT {
        /// `.n`
        n: &'static [f64],
        /// `.v`
        v: &'static [f64],
        /// `.Tcrit`
        tcrit: f64,
    },
    /// `IdealGasHelmholtzEnthalpyEntropyOffset` — `a1 + a2*tau`, the
    /// document's built-in reference-state offset (upstream slot
    /// `EnthalpyEntropyOffsetCore`).
    EnthalpyEntropyOffset {
        /// `.a1`
        a1: f64,
        /// `.a2`
        a2: f64,
        /// `.reference`
        reference: &'static str,
    },
    /// `IdealGasHelmholtzPower` — `sum n_k * tau^t_k`
    Power {
        /// `.n`
        n: &'static [f64],
        /// `.t`
        t: &'static [f64],
    },
}

/// `EOS[0].alphar[i]`, tagged by `type`.
pub enum AlpharTerm {
    /// `ResidualHelmholtzPower` —
    /// `sum n_k * delta^d_k * tau^t_k * exp(-delta^l_k)` (plain power when
    /// `l_k == 0`)
    Power {
        /// `.n`
        n: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.l`
        l: &'static [f64],
    },
    /// `ResidualHelmholtzGaussian` —
    /// `sum n_k * delta^d_k * tau^t_k * exp(-eta_k*(delta-epsilon_k)^2 - beta_k*(tau-gamma_k)^2)`
    Gaussian {
        /// `.n`
        n: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.eta`
        eta: &'static [f64],
        /// `.beta`
        beta: &'static [f64],
        /// `.gamma`
        gamma: &'static [f64],
        /// `.epsilon`
        epsilon: &'static [f64],
    },
    /// `ResidualHelmholtzNonAnalytic` (IAPWS-95 style critical terms)
    NonAnalytic {
        /// `.n`
        n: &'static [f64],
        /// `.a`
        a: &'static [f64],
        /// `.b`
        b: &'static [f64],
        /// `.beta`
        beta: &'static [f64],
        /// `.A`
        big_a: &'static [f64],
        /// `.B`
        big_b: &'static [f64],
        /// `.C`
        big_c: &'static [f64],
        /// `.D`
        big_d: &'static [f64],
    },
    /// `ResidualHelmholtzGaoB` (Gao et al. modified-Gaussian terms) —
    /// `sum n_k * tau^t_k * delta^d_k
    ///  * exp(eta_k*(delta-epsilon_k)^2 + 1/(beta_k*(tau-gamma_k)^2 + b_k))`
    GaoB {
        /// `.n`
        n: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.eta`
        eta: &'static [f64],
        /// `.beta`
        beta: &'static [f64],
        /// `.gamma`
        gamma: &'static [f64],
        /// `.epsilon`
        epsilon: &'static [f64],
        /// `.b`
        b: &'static [f64],
    },
}

/// `ANCILLARIES` — saturation ancillary equations used to seed solvers.
pub struct Ancillaries {
    /// `ANCILLARIES.pS`
    pub p_s: SaturationAncillary,
    /// `ANCILLARIES.rhoL`
    pub rho_l: SaturationAncillary,
    /// `ANCILLARIES.rhoV`
    pub rho_v: SaturationAncillary,
}

/// One saturation ancillary fit (`pS`/`rhoL`/`rhoV` objects).
pub struct SaturationAncillary {
    /// `.type` — upstream equation-form tag (`"pV"`, `"rhoLnoexp"`, ...)
    pub anc_type: &'static str,
    /// `.n`
    pub n: &'static [f64],
    /// `.t`
    pub t: &'static [f64],
    /// `.T_r` [K]
    pub t_r: f64,
    /// `.reducing_value` (Pa or mol/m^3)
    pub reducing_value: f64,
    /// `.using_tau_r`
    pub using_tau_r: bool,
    /// `.Tmin` [K]
    pub t_min: f64,
    /// `.Tmax` [K]
    pub t_max: f64,
}

/// `STATES` — tabulated characteristic states of the fluid document.
pub struct States {
    /// `STATES.critical`
    pub critical: StatePoint,
    /// `STATES.triple_liquid`
    pub triple_liquid: StatePoint,
    /// `STATES.triple_vapor`
    pub triple_vapor: StatePoint,
}
