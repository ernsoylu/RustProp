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
    /// `TRANSPORT` — None when the document has no TRANSPORT block at all
    pub transport: Option<Transport>,
}

/// `TRANSPORT` (ported subset).
pub struct Transport {
    /// `TRANSPORT.viscosity`
    pub viscosity: TransportModel<ViscosityModel>,
    /// `TRANSPORT.conductivity`
    pub conductivity: TransportModel<ConductivityModel>,
}

/// A viscosity model: the structured dilute/initial/higher form, or one of
/// upstream's fully-hardcoded per-fluid formulations.
// One static instance per fluid document; the variant size spread is
// irrelevant for rodata and boxing would break const construction.
#[allow(clippy::large_enum_variant)]
pub enum ViscosityModel {
    /// Structured sections
    Structured(Viscosity),
    /// top-level `.hardcoded` tag (Water, HeavyWater, Helium, R23,
    /// Methanol, the xylenes)
    Hardcoded {
        /// the `.hardcoded` tag
        name: &'static str,
    },
    /// `type: "Chung"` — the generalized Chung correlation from critical
    /// parameters (upstream evaluates with kappa = 0 regardless of the
    /// document's `kappa`)
    Chung {
        /// `.rhomolar_critical` [mol/m^3]
        rhomolar_critical: f64,
        /// `.acentric`
        acentric: f64,
        /// `.molar_mass` [kg/mol]
        molar_mass: f64,
        /// `.T_critical` [K]
        t_critical: f64,
        /// `.dipole_moment_D` [Debye]
        dipole_moment_d: f64,
        /// `.kappa` (stored for fidelity; unused by upstream's evaluation)
        kappa: f64,
    },
    /// `type: "ECS"` — extended corresponding states against a reference
    /// fluid (the L-J parameters may be absent; upstream then estimates
    /// them via `default_transport` — NaN encodes that here)
    Ecs {
        /// `.reference_fluid`
        reference_fluid: &'static str,
        /// `.psi.a`
        psi_a: &'static [f64],
        /// `.psi.t`
        psi_t: &'static [f64],
        /// `.psi.rhomolar_reducing` [mol/m^3]
        psi_rhomolar_reducing: f64,
        /// top-level `.sigma_eta` [m] (NaN when absent)
        sigma_eta: f64,
        /// top-level `.epsilon_over_k` [K] (NaN when absent)
        epsilon_over_k: f64,
    },
    /// `type: "rhosr-CS"` — residual-entropy-scaled corresponding states;
    /// the dilute part is kinetic theory with Chung-estimated L-J
    /// parameters from the reducing state (upstream `default_transport`)
    RhosrCs {
        /// `.C`
        c: f64,
        /// `.c_liq`
        c_liq: &'static [f64],
        /// `.c_vap`
        c_vap: &'static [f64],
        /// `.rhosr_critical`
        rhosr_critical: f64,
        /// `.x_crossover` (stored for fidelity; upstream's evaluation uses
        /// the literal 2)
        x_crossover: f64,
    },
}

/// A conductivity model: structured trio or fully-hardcoded.
#[allow(clippy::large_enum_variant)] // see ViscosityModel
pub enum ConductivityModel {
    /// Structured sections
    Structured(Conductivity),
    /// top-level `.hardcoded` tag (Water, HeavyWater, Helium, R23, Methane)
    Hardcoded {
        /// the `.hardcoded` tag
        name: &'static str,
    },
    /// `type: "ECS"` — extended corresponding states; the critical
    /// enhancement is Olchowy-Sengers with pure struct defaults (upstream
    /// does not read the block's informational `q_D`)
    Ecs {
        /// `.reference_fluid`
        reference_fluid: &'static str,
        /// `.psi.a`
        psi_a: &'static [f64],
        /// `.psi.t`
        psi_t: &'static [f64],
        /// `.psi.rhomolar_reducing` [mol/m^3]
        psi_rhomolar_reducing: f64,
        /// `.f_int.a`
        f_int_a: &'static [f64],
        /// `.f_int.t`
        f_int_t: &'static [f64],
        /// `.f_int.T_reducing` [K]
        f_int_t_reducing: f64,
    },
}

/// Per-property model slot: upstream distinguishes "no model provided"
/// (ValueError) from a model this port has not implemented yet.
pub enum TransportModel<T: 'static> {
    /// The property key is absent from the document
    Absent,
    /// The document carries a model class that is not ported yet
    /// (ECS/Chung/rhosr-CS lists/fully-hardcoded)
    Unported,
    /// A ported structured model
    Model(T),
}

/// `TRANSPORT.conductivity` (structured form); the assembly is
/// `lambda = dilute + residual + critical`.
pub struct Conductivity {
    /// `.dilute`
    pub dilute: ConductivityDilute,
    /// `.residual`
    pub residual: ConductivityResidual,
    /// `.critical` — absent means no critical enhancement
    pub critical: Option<ConductivityCritical>,
}

/// `TRANSPORT.conductivity.dilute`, tagged by `type` (or `hardcoded`).
pub enum ConductivityDilute {
    /// `ratio_of_polynomials` — `sum A_i Tr^n_i / sum B_i Tr^m_i`
    RatioOfPolynomials {
        /// `.A`
        a: &'static [f64],
        /// `.n`
        n: &'static [f64],
        /// `.B`
        b: &'static [f64],
        /// `.m`
        m: &'static [f64],
        /// `.T_reducing` [K]
        t_reducing: f64,
    },
    /// `eta0_and_poly` — `A_0*eta0[uPa-s] + sum_{i>=1} A_i*tau^t_i`
    Eta0AndPoly {
        /// `.A`
        a: &'static [f64],
        /// `.t`
        t: &'static [f64],
    },
    /// `.hardcoded` — ports with the hardcoded slice
    Hardcoded {
        /// the `.hardcoded` tag
        name: &'static str,
    },
}

/// `TRANSPORT.conductivity.residual`, tagged by `type`.
pub enum ConductivityResidual {
    /// `polynomial` — `sum B_i tau^t_i delta^d_i` with delta from the MASS
    /// density over `rhomass_reducing`
    Polynomial {
        /// `.B`
        b: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.T_reducing` [K]
        t_reducing: f64,
        /// `.rhomass_reducing` [kg/m^3]
        rhomass_reducing: f64,
    },
    /// `polynomial_and_exponential` — EOS tau/delta with
    /// `exp(-gamma_i*delta^l_i)` factors
    PolynomialAndExponential {
        /// `.A`
        a: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.gamma`
        gamma: &'static [f64],
        /// `.l`
        l: &'static [f64],
    },
}

/// `TRANSPORT.conductivity.critical`, tagged by `type` (or `hardcoded`).
pub enum ConductivityCritical {
    /// `simplified_Olchowy_Sengers`; absent JSON keys carry upstream's
    /// defaults (k = 1.3806488e-23, R0 = 1.03, gamma = 1.239, nu = 0.63,
    /// GAMMA = 0.0496, zeta0 = 1.94e-10, qD = 2e9); `t_ref` is NaN for
    /// upstream's 1.5*T_reducing default
    SimplifiedOlchowySengers {
        /// `.k` (default)
        k: f64,
        /// `.R0`
        r0: f64,
        /// `.gamma`
        gamma: f64,
        /// `.nu` (default)
        nu: f64,
        /// `.GAMMA`
        big_gamma: f64,
        /// `.zeta0` [m]
        zeta0: f64,
        /// `.qD` [1/m]
        qd: f64,
        /// `.T_ref` [K] (NaN -> 1.5*T_reducing)
        t_ref: f64,
    },
    /// `.hardcoded` — ports with the hardcoded slice
    Hardcoded {
        /// the `.hardcoded` tag
        name: &'static str,
    },
}

/// `TRANSPORT.viscosity` (structured form).
pub struct Viscosity {
    /// `.epsilon_over_k` [K] (NaN when the document omits it)
    pub epsilon_over_k: f64,
    /// `.sigma_eta` [m] (NaN when the document omits it)
    pub sigma_eta: f64,
    /// `.dilute`
    pub dilute: ViscosityDilute,
    /// `.initial_density` — absent for many fluids
    pub initial_density: Option<ViscosityInitialDensity>,
    /// `.higher_order`
    pub higher_order: ViscosityHigherOrder,
}

/// `TRANSPORT.viscosity.dilute`, tagged by `type` (or `hardcoded`).
pub enum ViscosityDilute {
    /// `kinetic_theory` — Neufeld Omega22 from the top-level
    /// epsilon_over_k/sigma_eta
    KineticTheory,
    /// `collision_integral` — note `.molar_mass` is the block's own value,
    /// not the EOS's
    CollisionIntegral {
        /// `.a`
        a: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.C`
        c: f64,
        /// `.molar_mass` [kg/mol]
        molar_mass: f64,
    },
    /// `powers_of_T`
    PowersOfT {
        /// `.a`
        a: &'static [f64],
        /// `.t`
        t: &'static [f64],
    },
    /// `powers_of_Tr`
    PowersOfTr {
        /// `.a`
        a: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.T_reducing` [K]
        t_reducing: f64,
    },
    /// `collision_integral_powers_of_Tstar`
    CollisionIntegralPowersOfTstar {
        /// `.a`
        a: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.C`
        c: f64,
        /// `.T_reducing` [K]
        t_reducing: f64,
    },
    /// `.hardcoded` — evaluation ports with the hardcoded slice
    Hardcoded {
        /// the `.hardcoded` tag
        name: &'static str,
    },
}

/// `TRANSPORT.viscosity.initial_density`, tagged by `type`.
pub enum ViscosityInitialDensity {
    /// `Rainwater-Friend` — returns B_eta [m^3/mol]; the contribution is
    /// `eta_dilute * B_eta * rhomolar`
    RainwaterFriend {
        /// `.b`
        b: &'static [f64],
        /// `.t`
        t: &'static [f64],
    },
    /// `empirical`
    Empirical {
        /// `.n`
        n: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.T_reducing` [K]
        t_reducing: f64,
        /// `.rhomolar_reducing` [mol/m^3]
        rhomolar_reducing: f64,
    },
}

/// `TRANSPORT.viscosity.higher_order`, tagged by `type` (or `hardcoded`).
pub enum ViscosityHigherOrder {
    /// `modified_Batschinski_Hildebrand`
    ModifiedBatschinskiHildebrand {
        /// `.a`
        a: &'static [f64],
        /// `.d1`
        d1: &'static [f64],
        /// `.t1`
        t1: &'static [f64],
        /// `.gamma`
        gamma: &'static [f64],
        /// `.l`
        l: &'static [f64],
        /// `.f`
        f: &'static [f64],
        /// `.d2`
        d2: &'static [f64],
        /// `.t2`
        t2: &'static [f64],
        /// `.g`
        g: &'static [f64],
        /// `.h`
        h: &'static [f64],
        /// `.p`
        p: &'static [f64],
        /// `.q`
        q: &'static [f64],
        /// `.T_reduce` [K]
        t_reduce: f64,
        /// `.rhomolar_reduce` [mol/m^3]
        rhomolar_reduce: f64,
    },
    /// `friction_theory` (optional channels carry empty slices / zero
    /// exponents exactly as upstream's empty vectors)
    FrictionTheory {
        /// `.Ai`
        ai: &'static [f64],
        /// `.Aa`
        aa: &'static [f64],
        /// `.Ar`
        ar: &'static [f64],
        /// `.Aaa`
        aaa: &'static [f64],
        /// `.Arr` (empty when `.Adrdr` is given)
        arr: &'static [f64],
        /// `.Adrdr` (empty when `.Arr` is given)
        adrdr: &'static [f64],
        /// `.Aii` (optional)
        aii: &'static [f64],
        /// `.Arrr` (optional, with `.Aaaa`)
        arrr: &'static [f64],
        /// `.Aaaa` (optional, with `.Arrr`)
        aaaa: &'static [f64],
        /// `.Na`
        na: f64,
        /// `.Naa`
        naa: f64,
        /// `.Nr`
        nr: f64,
        /// `.Nrr`
        nrr: f64,
        /// `.Nii` (0 when absent)
        nii: f64,
        /// `.Nrrr` (0 when absent)
        nrrr: f64,
        /// `.Naaa` (0 when absent)
        naaa: f64,
        /// `.c1`
        c1: f64,
        /// `.c2`
        c2: f64,
        /// `.T_reduce` [K]
        t_reduce: f64,
    },
    /// `.hardcoded` — evaluation ports with the hardcoded slice
    Hardcoded {
        /// the `.hardcoded` tag
        name: &'static str,
    },
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
    /// `EOS[0].STATES.temperature_max_sat` — the saturation-temperature
    /// maximum of a pseudo-pure fluid (upstream `max_sat_T`); `None` for
    /// pure fluids.
    pub max_sat_t: Option<StatePoint>,
    /// `EOS[0].STATES.pressure_max_sat` (upstream `max_sat_p`).
    pub max_sat_p: Option<StatePoint>,
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
    /// `IdealGasHelmholtzPlanckEinsteinGeneralized` —
    /// `sum n_k * ln(c_k + d_k*exp(t_k*tau))` (t is theta directly); merges
    /// into the generalized Planck-Einstein container at parse time.
    PlanckEinsteinGeneralized {
        /// `.n`
        n: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.c`
        c: &'static [f64],
        /// `.d`
        d: &'static [f64],
    },
    /// `IdealGasHelmholtzCP0Constant` — constant `cp0/R` over the range,
    /// anchored at `T0`.
    Cp0Constant {
        /// `.cp_over_R`
        cp_over_r: f64,
        /// `.Tc`
        tc: f64,
        /// `.T0`
        t0: f64,
    },
    /// `IdealGasHelmholtzCP0PolyT` — `cp0/R = sum c_k * T^t_k`.
    Cp0PolyT {
        /// `.c`
        c: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.Tc`
        tc: f64,
        /// `.T0`
        t0: f64,
    },
    /// `IdealGasHelmholtzCP0AlyLee` — Aly-Lee cp0 form (5 constants);
    /// upstream converts it at parse time into a CP0PolyT constant plus
    /// sinh/cosh Planck-Einstein-generalized entries.
    Cp0AlyLee {
        /// `.c` (A, B, C, D, E)
        c: &'static [f64],
        /// `.Tc`
        tc: f64,
        /// `.T0`
        t0: f64,
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
    /// `ResidualHelmholtzExponential` —
    /// `sum n_k * delta^d_k * tau^t_k * exp(-g_k*delta^l_k)`
    Exponential {
        /// `.n`
        n: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.g`
        g: &'static [f64],
        /// `.l`
        l: &'static [f64],
    },
    /// `ResidualHelmholtzDoubleExponential` —
    /// `sum n_k * delta^d_k * tau^t_k * exp(-gd_k*delta^ld_k - gt_k*tau^lt_k)`
    DoubleExponential {
        /// `.n`
        n: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.gd`
        gd: &'static [f64],
        /// `.ld`
        ld: &'static [f64],
        /// `.gt`
        gt: &'static [f64],
        /// `.lt`
        lt: &'static [f64],
    },
    /// `ResidualHelmholtzLemmon2005` —
    /// `sum n_k * delta^d_k * tau^t_k * exp(-delta^l_k - tau^m_k)`
    Lemmon2005 {
        /// `.n`
        n: &'static [f64],
        /// `.d`
        d: &'static [f64],
        /// `.t`
        t: &'static [f64],
        /// `.l`
        l: &'static [f64],
        /// `.m`
        m: &'static [f64],
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
    /// Upstream's `ancillaries.pL` slot: the `pS` curve for a pure fluid
    /// (upstream loads `pS` into BOTH pL and pV), the bubble-point `pL`
    /// curve for a pseudo-pure fluid.
    pub p_s: SaturationAncillary,
    /// Upstream's `ancillaries.pV` slot when it differs from pL — the
    /// dew-point curve of a pseudo-pure fluid. `None` for pure fluids (the
    /// slot aliases `p_s`, saving a duplicate static curve).
    pub p_v_split: Option<SaturationAncillary>,
    /// `ANCILLARIES.rhoL`
    pub rho_l: SaturationAncillary,
    /// `ANCILLARIES.rhoV`
    pub rho_v: SaturationAncillary,
    /// `ANCILLARIES.surface_tension` — absent for fluids without a curve
    pub surface_tension: Option<SurfaceTension>,
    /// `ANCILLARIES.melting_line` — absent for fluids without a curve
    pub melting_line: Option<MeltingLine>,
}

/// `ANCILLARIES.melting_line` (upstream `MeltingLineVariables`): the
/// pressure-temperature melting curve in one of three segment forms. The
/// per-part and aggregate p/T limits are computed at runtime (upstream
/// `set_limits`), not stored.
pub struct MeltingLine {
    /// `.T_m` [K] — the document's normal melting temperature (parsed and
    /// carried like upstream's field; -1 encodes "not provided").
    pub t_m: f64,
    /// `.type` + `.parts`
    pub kind: MeltingLineKind,
}

/// The three upstream melting-curve segment families.
pub enum MeltingLineKind {
    /// `Simon`: `p = p_0 + a*((T/T_0)^c - 1)`
    Simon { parts: &'static [SimonMeltPart] },
    /// `polynomial_in_Tr`: `p = p_0*(1 + sum a_i*((T/T_0)^t_i - 1))`
    PolynomialInTr { parts: &'static [PolyMeltPart] },
    /// `polynomial_in_Theta`: `p = p_0*(1 + sum a_i*(T/T_0 - 1)^t_i)`
    PolynomialInTheta { parts: &'static [PolyMeltPart] },
}

/// One Simon-type melting segment (upstream
/// `MeltingLinePiecewiseSimonSegment`).
pub struct SimonMeltPart {
    /// `.T_0` [K]
    pub t_0: f64,
    /// `.a` [Pa]
    pub a: f64,
    /// `.c`
    pub c: f64,
    /// `.p_0` [Pa]
    pub p_0: f64,
    /// `.T_min` [K]
    pub t_min: f64,
    /// `.T_max` [K]
    pub t_max: f64,
}

/// One polynomial melting segment (upstream
/// `MeltingLinePiecewisePolynomialIn{Tr,Theta}Segment`).
pub struct PolyMeltPart {
    /// `.T_0` [K]
    pub t_0: f64,
    /// `.p_0` [Pa]
    pub p_0: f64,
    /// `.T_min` [K]
    pub t_min: f64,
    /// `.T_max` [K]
    pub t_max: f64,
    /// `.a`
    pub a: &'static [f64],
    /// `.t`
    pub t: &'static [f64],
}

/// `ANCILLARIES.surface_tension` (upstream `SurfaceTensionCorrelation`):
/// `sigma = sum a_i * (1 - T/Tc)^n_i` [N/m].
pub struct SurfaceTension {
    /// `.a`
    pub a: &'static [f64],
    /// `.n`
    pub n: &'static [f64],
    /// `.Tc` [K]
    pub tc: f64,
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

/// One fluid of the cubic backend's library (upstream
/// `CubicsLibrary::CubicsValues`, loaded from `dev/cubics/all_cubic_fluids.json`).
/// SRK and Peng-Robinson share this data; only Tc, pc, acentric, and the
/// ideal-gas `alpha0` terms enter the equations of state.
pub struct CubicFluid {
    /// `.name`
    pub name: &'static str,
    /// `.CAS`
    pub cas: &'static str,
    /// `.aliases`
    pub aliases: &'static [&'static str],
    /// `.Tc` [K]
    pub tc: f64,
    /// `.pc` [Pa]
    pub pc: f64,
    /// `.rhomolarc` [mol/m^3] — the JSON critical density; used ONLY for the
    /// alpha0 reduced variables (the `rhomolar_critical` OUTPUT is upstream's
    /// Kazakov curve fit, a different number).
    pub rhomolarc: f64,
    /// `.acentric`
    pub acentric: f64,
    /// `.molemass` [kg/mol]
    pub molemass: f64,
    /// `.alpha0` — same term families and parse normalizations as the HEOS
    /// fluid library (upstream shares `parse_alpha0`).
    pub alpha0: &'static [Alpha0Term],
}

/// One property block of an incompressible fluid (upstream
/// `IncompressibleData`): the five coefficient forms plus the explicit
/// "notdefined" absence every document spells out.
pub enum IncompData {
    /// `polynomial` — 2-D matrix, rows = T powers, cols = x powers.
    Polynomial(&'static [&'static [f64]]),
    /// `exppolynomial` — `exp(polynomial)`.
    ExpPolynomial(&'static [&'static [f64]]),
    /// `exponential` — 3 coefficients.
    Exponential(&'static [f64]),
    /// `logexponential` — 3 coefficients.
    LogExponential(&'static [f64]),
    /// `polyoffset` — 1-D vector, first entry the centering value.
    PolyOffset(&'static [f64]),
    /// `notdefined` (or an unknown non-vital tag).
    NotSet,
}

/// The concentration basis of an incompressible fluid (upstream `xid`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IncompFrac {
    Pure,
    Mass,
    Volume,
}

/// One incompressible fluid (upstream `IncompressibleFluid`, loaded from
/// `dev/incompressible_liquids/json`). The dead-in-v8 `mass2input`/
/// `mole2input`/`volume2input` conversion polynomials are carried for
/// data fidelity but never evaluated (upstream's converters are stubs).
pub struct IncompFluid {
    /// `.name`
    pub name: &'static str,
    /// `.Tmin` [K]
    pub tmin: f64,
    /// `.Tmax` [K]
    pub tmax: f64,
    /// `.xmin`
    pub xmin: f64,
    /// `.xmax`
    pub xmax: f64,
    /// `.xid`
    pub xid: IncompFrac,
    /// `.TminPsat` [K]
    pub tmin_psat: f64,
    /// `.Tbase` [K]
    pub tbase: f64,
    /// `.xbase`
    pub xbase: f64,
    /// `.density` (vital)
    pub density: IncompData,
    /// `.specific_heat` (vital)
    pub specific_heat: IncompData,
    /// `.conductivity`
    pub conductivity: IncompData,
    /// `.viscosity`
    pub viscosity: IncompData,
    /// `.saturation_pressure`
    pub saturation_pressure: IncompData,
    /// `.T_freeze`
    pub t_freeze: IncompData,
    /// `.mass2input` (dead in v8)
    pub mass2input: IncompData,
    /// `.mole2input` (dead in v8)
    pub mole2input: IncompData,
    /// `.volume2input` (dead in v8)
    pub volume2input: IncompData,
}

/// One GERG-2008 binary interaction record (upstream
/// `mixture_binary_pairs.json`; the six Lemmon `xi`/`zeta` records are
/// converted to this form at datagen time as upstream converts them at load
/// time). Keys are CAS-sorted: `cas1 < cas2` lexicographically; `beta_*`
/// invert (`1/beta`) when the component order is swapped, `gamma_*` are
/// symmetric.
pub struct MixBinaryPair {
    pub cas1: &'static str,
    pub cas2: &'static str,
    pub beta_t: f64,
    pub gamma_t: f64,
    pub beta_v: f64,
    pub gamma_v: f64,
    /// The departure-function weight; 0 for the 848 reducing-only pairs.
    pub f: f64,
    /// Departure-function name (present only when `f != 0`).
    pub function: Option<&'static str>,
}

/// One mixture departure function (upstream
/// `mixture_departure_functions.json`): the three upstream types all map
/// onto the generalized-exponential term machinery.
/// One predefined mixture (upstream `predefined_mixtures.json` entry; the
/// registry key is `name + ".mix"` plus its uppercase form).
pub struct PredefinedMixture {
    /// Base name WITHOUT the `.mix` suffix (as shipped: "Air", "R410A", ...).
    pub name: &'static str,
    /// Component names as shipped (REFPROP-style, resolved through the
    /// fluid registry's alias map: "METHANE", "CO2", "ISOBUTAN", ...).
    pub fluids: &'static [&'static str],
    pub mole_fractions: &'static [f64],
}

pub struct MixDepartureFn {
    pub name: &'static str,
    pub kind: MixDepartureKind,
    /// Power-term count (the first `npower` terms of the arrays).
    pub npower: usize,
    pub n: &'static [f64],
    pub d: &'static [f64],
    pub t: &'static [f64],
    /// `l` exponents (Exponential and Gaussian+Exponential types).
    pub l: &'static [f64],
    /// Gaussian tail parameters (GERG-2008 and Gaussian+Exponential).
    pub eta: &'static [f64],
    pub epsilon: &'static [f64],
    pub beta: &'static [f64],
    pub gamma: &'static [f64],
}

/// The three upstream departure-function types.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MixDepartureKind {
    /// Power terms + the GERG gaussian `exp(-eta(d-e)^2 - beta(d-gamma))`
    /// (linear in delta).
    Gerg2008,
    /// Pure `n d^d t^t exp(-delta^l)` terms.
    Exponential,
    /// Power terms + standard gaussians.
    GaussianExponential,
}
