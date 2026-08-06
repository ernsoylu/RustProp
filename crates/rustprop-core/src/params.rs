//! Parameter system ported from upstream `include/CoolProp/DataStructures.h`
//! and `src/DataStructures.cpp` @ v8.0.0.
//!
//! Naming rule: upstream `iX_y` enumerators become CamelCase variants with the
//! `i` prefix stripped (`iT` → `T`, `igas_constant` → `GasConstant`,
//! `idalphar_dtau_constdelta` → `DalpharDtauConstdelta`). Discriminants equal
//! the upstream integer values; the golden test `params.rs` in
//! `rustprop-golden-tests` verifies every index, name, IO class, unit string,
//! description, and trivial flag against a dump from the CoolProp 8.0.0 wheel.
//!
//! Deviations from upstream (logged in PLAN.md):
//! - the `INVALID_PARAMETER`/`iundefined_parameter` and `INPUT_PAIR_INVALID`
//!   sentinels are not ported — Rust uses `Option`/`Result` instead;
//! - lookups return `Option` instead of throwing; callers produce the
//!   matching `Error` condition at the API boundary.

/// Input/output parameter, mirroring upstream `enum parameters`.
#[allow(clippy::upper_case_acronyms)] // names mirror upstream (PIP, GWP20, FH, HH, PH, ODP)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum Param {
    // General parameters
    GasConstant = 1,
    MolarMass,
    AcentricFactor,
    RhomolarReducing,
    RhomolarCritical,
    TReducing,
    TCritical,
    RhomassReducing,
    RhomassCritical,
    PCritical,
    PReducing,
    TTriple,
    PTriple,
    TMin,
    TMax,
    PMax,
    PMin,
    DipoleMoment,
    // Bulk properties
    T,
    P,
    Q,
    Qmass,
    Tau,
    Delta,
    // Molar specific thermodynamic properties
    Dmolar,
    Hmolar,
    Smolar,
    Cpmolar,
    Cp0molar,
    Cvmolar,
    Umolar,
    Gmolar,
    Helmholtzmolar,
    HmolarResidual,
    SmolarResidual,
    GmolarResidual,
    HmolarIdealgas,
    SmolarIdealgas,
    UmolarIdealgas,
    // Mass specific thermodynamic properties
    Dmass,
    Hmass,
    Smass,
    Cpmass,
    Cp0mass,
    Cvmass,
    Umass,
    Gmass,
    Helmholtzmass,
    HmassIdealgas,
    SmassIdealgas,
    UmassIdealgas,
    // Transport properties
    Viscosity,
    Conductivity,
    SurfaceTension,
    Prandtl,
    // Derivative-based terms
    SpeedSound,
    IsothermalCompressibility,
    IsobaricExpansionCoefficient,
    IsentropicExpansionCoefficient,
    // Fundamental derivative of gas dynamics
    FundamentalDerivativeOfGasDynamics,
    // Derivatives of the residual non-dimensionalized Helmholtz energy
    Alphar,
    DalpharDtauConstdelta,
    DalpharDdeltaConsttau,
    // Derivatives of the ideal-gas non-dimensionalized Helmholtz energy
    Alpha0,
    Dalpha0DtauConstdelta,
    Dalpha0DdeltaConsttau,
    D2alpha0Ddelta2Consttau,
    D3alpha0Ddelta3Consttau,
    // Other functions and derivatives
    Bvirial,
    Cvirial,
    DBvirialDT,
    DCvirialDT,
    Z,
    PIP,
    // Accessors for incompressibles
    FractionMin,
    FractionMax,
    TFreeze,
    // Environmental parameters
    GWP20,
    GWP100,
    GWP500,
    FH,
    HH,
    PH,
    ODP,
    Phase,
}

pub struct ParamInfo {
    pub param: Param,
    pub short: &'static str,
    /// `"IO"` if input/output, `"O"` if output only.
    pub io: &'static str,
    pub units: &'static str,
    pub long: &'static str,
    /// True if directly calculable (constants, critical parameters, ...).
    pub trivial: bool,
}

/// Transcribed from `parameter_info_list` in `src/DataStructures.cpp`
/// @ v8.0.0, in upstream table order.
#[rustfmt::skip]
pub const PARAM_INFO: &[ParamInfo] = &[
    ParamInfo { param: Param::T, short: "T", io: "IO", units: "K", long: "Temperature", trivial: false },
    ParamInfo { param: Param::P, short: "P", io: "IO", units: "Pa", long: "Pressure", trivial: false },
    ParamInfo { param: Param::Dmolar, short: "Dmolar", io: "IO", units: "mol/m^3", long: "Molar density", trivial: false },
    ParamInfo { param: Param::Hmolar, short: "Hmolar", io: "IO", units: "J/mol", long: "Molar specific enthalpy", trivial: false },
    ParamInfo { param: Param::Smolar, short: "Smolar", io: "IO", units: "J/mol/K", long: "Molar specific entropy", trivial: false },
    ParamInfo { param: Param::Umolar, short: "Umolar", io: "IO", units: "J/mol", long: "Molar specific internal energy", trivial: false },
    ParamInfo { param: Param::Gmolar, short: "Gmolar", io: "O", units: "J/mol", long: "Molar specific Gibbs energy", trivial: false },
    ParamInfo { param: Param::Helmholtzmolar, short: "Helmholtzmolar", io: "O", units: "J/mol", long: "Molar specific Helmholtz energy", trivial: false },
    ParamInfo { param: Param::Dmass, short: "Dmass", io: "IO", units: "kg/m^3", long: "Mass density", trivial: false },
    ParamInfo { param: Param::Hmass, short: "Hmass", io: "IO", units: "J/kg", long: "Mass specific enthalpy", trivial: false },
    ParamInfo { param: Param::Smass, short: "Smass", io: "IO", units: "J/kg/K", long: "Mass specific entropy", trivial: false },
    ParamInfo { param: Param::Umass, short: "Umass", io: "IO", units: "J/kg", long: "Mass specific internal energy", trivial: false },
    ParamInfo { param: Param::Gmass, short: "Gmass", io: "O", units: "J/kg", long: "Mass specific Gibbs energy", trivial: false },
    ParamInfo { param: Param::Helmholtzmass, short: "Helmholtzmass", io: "O", units: "J/kg", long: "Mass specific Helmholtz energy", trivial: false },
    ParamInfo { param: Param::Q, short: "Q", io: "IO", units: "mol/mol", long: "Molar vapor quality", trivial: false },
    ParamInfo { param: Param::Qmass, short: "Qmass", io: "IO", units: "kg/kg", long: "Mass-basis vapor quality", trivial: false },
    ParamInfo { param: Param::Delta, short: "Delta", io: "IO", units: "-", long: "Reduced density (rho/rhoc)", trivial: false },
    ParamInfo { param: Param::Tau, short: "Tau", io: "IO", units: "-", long: "Reciprocal reduced temperature (Tc/T)", trivial: false },
    ParamInfo { param: Param::Cpmolar, short: "Cpmolar", io: "O", units: "J/mol/K", long: "Molar specific constant pressure specific heat", trivial: false },
    ParamInfo { param: Param::Cpmass, short: "Cpmass", io: "O", units: "J/kg/K", long: "Mass specific constant pressure specific heat", trivial: false },
    ParamInfo { param: Param::Cvmolar, short: "Cvmolar", io: "O", units: "J/mol/K", long: "Molar specific constant volume specific heat", trivial: false },
    ParamInfo { param: Param::Cvmass, short: "Cvmass", io: "O", units: "J/kg/K", long: "Mass specific constant volume specific heat", trivial: false },
    ParamInfo { param: Param::Cp0molar, short: "Cp0molar", io: "O", units: "J/mol/K", long: "Ideal gas molar specific constant pressure specific heat", trivial: false },
    ParamInfo { param: Param::Cp0mass, short: "Cp0mass", io: "O", units: "J/kg/K", long: "Ideal gas mass specific constant pressure specific heat", trivial: false },
    ParamInfo { param: Param::HmolarResidual, short: "Hmolar_residual", io: "O", units: "J/mol", long: "Residual molar enthalpy", trivial: false },
    ParamInfo { param: Param::SmolarResidual, short: "Smolar_residual", io: "O", units: "J/mol/K", long: "Residual molar entropy (sr/R = s(T,rho) - s^0(T,rho))", trivial: false },
    ParamInfo { param: Param::GmolarResidual, short: "Gmolar_residual", io: "O", units: "J/mol", long: "Residual molar Gibbs energy", trivial: false },
    ParamInfo { param: Param::HmolarIdealgas, short: "Hmolar_idealgas", io: "O", units: "J/mol", long: "Ideal gas molar enthalpy", trivial: false },
    ParamInfo { param: Param::SmolarIdealgas, short: "Smolar_idealgas", io: "O", units: "J/mol/K", long: "Ideal gas molar entropy", trivial: false },
    ParamInfo { param: Param::UmolarIdealgas, short: "Umolar_idealgas", io: "O", units: "J/mol", long: "Ideal gas molar internal energy", trivial: false },
    ParamInfo { param: Param::HmassIdealgas, short: "Hmass_idealgas", io: "O", units: "J/kg", long: "Ideal gas specific enthalpy", trivial: false },
    ParamInfo { param: Param::SmassIdealgas, short: "Smass_idealgas", io: "O", units: "J/kg/K", long: "Ideal gas specific entropy", trivial: false },
    ParamInfo { param: Param::UmassIdealgas, short: "Umass_idealgas", io: "O", units: "J/kg", long: "Ideal gas specific internal energy", trivial: false },
    ParamInfo { param: Param::GWP20, short: "GWP20", io: "O", units: "-", long: "20-year global warming potential", trivial: true },
    ParamInfo { param: Param::GWP100, short: "GWP100", io: "O", units: "-", long: "100-year global warming potential", trivial: true },
    ParamInfo { param: Param::GWP500, short: "GWP500", io: "O", units: "-", long: "500-year global warming potential", trivial: true },
    ParamInfo { param: Param::FH, short: "FH", io: "O", units: "-", long: "Flammability hazard", trivial: true },
    ParamInfo { param: Param::HH, short: "HH", io: "O", units: "-", long: "Health hazard", trivial: true },
    ParamInfo { param: Param::PH, short: "PH", io: "O", units: "-", long: "Physical hazard", trivial: true },
    ParamInfo { param: Param::ODP, short: "ODP", io: "O", units: "-", long: "Ozone depletion potential", trivial: true },
    ParamInfo { param: Param::Bvirial, short: "Bvirial", io: "O", units: "-", long: "Second virial coefficient", trivial: false },
    ParamInfo { param: Param::Cvirial, short: "Cvirial", io: "O", units: "-", long: "Third virial coefficient", trivial: false },
    ParamInfo { param: Param::DBvirialDT, short: "dBvirial_dT", io: "O", units: "-", long: "Derivative of second virial coefficient with respect to T", trivial: false },
    ParamInfo { param: Param::DCvirialDT, short: "dCvirial_dT", io: "O", units: "-", long: "Derivative of third virial coefficient with respect to T", trivial: false },
    ParamInfo { param: Param::GasConstant, short: "gas_constant", io: "O", units: "J/mol/K", long: "Molar gas constant", trivial: true },
    ParamInfo { param: Param::MolarMass, short: "molar_mass", io: "O", units: "kg/mol", long: "Molar mass", trivial: true },
    ParamInfo { param: Param::AcentricFactor, short: "acentric", io: "O", units: "-", long: "Acentric factor", trivial: true },
    ParamInfo { param: Param::DipoleMoment, short: "dipole_moment", io: "O", units: "C-m", long: "Dipole moment", trivial: true },
    ParamInfo { param: Param::RhomassReducing, short: "rhomass_reducing", io: "O", units: "kg/m^3", long: "Mass density at reducing point", trivial: true },
    ParamInfo { param: Param::RhomolarReducing, short: "rhomolar_reducing", io: "O", units: "mol/m^3", long: "Molar density at reducing point", trivial: true },
    ParamInfo { param: Param::RhomolarCritical, short: "rhomolar_critical", io: "O", units: "mol/m^3", long: "Molar density at critical point", trivial: true },
    ParamInfo { param: Param::RhomassCritical, short: "rhomass_critical", io: "O", units: "kg/m^3", long: "Mass density at critical point", trivial: true },
    ParamInfo { param: Param::TReducing, short: "T_reducing", io: "O", units: "K", long: "Temperature at the reducing point", trivial: true },
    ParamInfo { param: Param::TCritical, short: "T_critical", io: "O", units: "K", long: "Temperature at the critical point", trivial: true },
    ParamInfo { param: Param::TTriple, short: "T_triple", io: "O", units: "K", long: "Temperature at the triple point", trivial: true },
    ParamInfo { param: Param::TMax, short: "T_max", io: "O", units: "K", long: "Maximum temperature limit", trivial: true },
    ParamInfo { param: Param::TMin, short: "T_min", io: "O", units: "K", long: "Minimum temperature limit", trivial: true },
    ParamInfo { param: Param::PMin, short: "P_min", io: "O", units: "Pa", long: "Minimum pressure limit", trivial: true },
    ParamInfo { param: Param::PMax, short: "P_max", io: "O", units: "Pa", long: "Maximum pressure limit", trivial: true },
    ParamInfo { param: Param::PCritical, short: "p_critical", io: "O", units: "Pa", long: "Pressure at the critical point", trivial: true },
    ParamInfo { param: Param::PReducing, short: "p_reducing", io: "O", units: "Pa", long: "Pressure at the reducing point", trivial: true },
    ParamInfo { param: Param::PTriple, short: "p_triple", io: "O", units: "Pa", long: "Pressure at the triple point (pure only)", trivial: true },
    ParamInfo { param: Param::FractionMin, short: "fraction_min", io: "O", units: "-", long: "Fraction (mole, mass, volume) minimum value for incompressible solutions", trivial: true },
    ParamInfo { param: Param::FractionMax, short: "fraction_max", io: "O", units: "-", long: "Fraction (mole, mass, volume) maximum value for incompressible solutions", trivial: true },
    ParamInfo { param: Param::TFreeze, short: "T_freeze", io: "O", units: "K", long: "Freezing temperature for incompressible solutions", trivial: true },
    ParamInfo { param: Param::SpeedSound, short: "speed_of_sound", io: "O", units: "m/s", long: "Speed of sound", trivial: false },
    ParamInfo { param: Param::Viscosity, short: "viscosity", io: "O", units: "Pa-s", long: "Viscosity", trivial: false },
    ParamInfo { param: Param::Conductivity, short: "conductivity", io: "O", units: "W/m/K", long: "Thermal conductivity", trivial: false },
    ParamInfo { param: Param::SurfaceTension, short: "surface_tension", io: "O", units: "N/m", long: "Surface tension", trivial: false },
    ParamInfo { param: Param::Prandtl, short: "Prandtl", io: "O", units: "-", long: "Prandtl number", trivial: false },
    ParamInfo { param: Param::IsothermalCompressibility, short: "isothermal_compressibility", io: "O", units: "1/Pa", long: "Isothermal compressibility", trivial: false },
    ParamInfo { param: Param::IsobaricExpansionCoefficient, short: "isobaric_expansion_coefficient", io: "O", units: "1/K", long: "Isobaric expansion coefficient", trivial: false },
    ParamInfo { param: Param::IsentropicExpansionCoefficient, short: "isentropic_expansion_coefficient", io: "O", units: "-", long: "Isentropic expansion coefficient", trivial: false },
    ParamInfo { param: Param::Z, short: "Z", io: "O", units: "-", long: "Compressibility factor", trivial: false },
    ParamInfo { param: Param::FundamentalDerivativeOfGasDynamics, short: "fundamental_derivative_of_gas_dynamics", io: "O", units: "-", long: "Fundamental derivative of gas dynamics", trivial: false },
    ParamInfo { param: Param::PIP, short: "PIP", io: "O", units: "-", long: "Phase identification parameter", trivial: false },
    ParamInfo { param: Param::Alphar, short: "alphar", io: "O", units: "-", long: "Residual Helmholtz energy", trivial: false },
    ParamInfo { param: Param::DalpharDtauConstdelta, short: "dalphar_dtau_constdelta", io: "O", units: "-", long: "Derivative of residual Helmholtz energy with tau", trivial: false },
    ParamInfo { param: Param::DalpharDdeltaConsttau, short: "dalphar_ddelta_consttau", io: "O", units: "-", long: "Derivative of residual Helmholtz energy with delta", trivial: false },
    ParamInfo { param: Param::Alpha0, short: "alpha0", io: "O", units: "-", long: "Ideal Helmholtz energy", trivial: false },
    ParamInfo { param: Param::Dalpha0DtauConstdelta, short: "dalpha0_dtau_constdelta", io: "O", units: "-", long: "Derivative of ideal Helmholtz energy with tau", trivial: false },
    ParamInfo { param: Param::Dalpha0DdeltaConsttau, short: "dalpha0_ddelta_consttau", io: "O", units: "-", long: "Derivative of ideal Helmholtz energy with delta", trivial: false },
    ParamInfo { param: Param::D2alpha0Ddelta2Consttau, short: "d2alpha0_ddelta2_consttau", io: "O", units: "-", long: "Second derivative of ideal Helmholtz energy with delta", trivial: false },
    ParamInfo { param: Param::D3alpha0Ddelta3Consttau, short: "d3alpha0_ddelta3_consttau", io: "O", units: "-", long: "Third derivative of ideal Helmholtz energy with delta", trivial: false },
    ParamInfo { param: Param::Phase, short: "Phase", io: "O", units: "-", long: "Phase index as a float", trivial: false },
];

/// Backward-compatibility aliases, transcribed from
/// `ParameterInformation::ParameterInformation()` in `src/DataStructures.cpp`.
#[rustfmt::skip]
pub const PARAM_ALIASES: &[(&str, Param)] = &[
    ("D", Param::Dmass),
    ("H", Param::Hmass),
    ("M", Param::MolarMass),
    ("S", Param::Smass),
    ("U", Param::Umass),
    ("C", Param::Cpmass),
    ("O", Param::Cvmass),
    ("G", Param::Gmass),
    ("V", Param::Viscosity),
    ("L", Param::Conductivity),
    ("pcrit", Param::PCritical),
    ("Pcrit", Param::PCritical),
    ("Tcrit", Param::TCritical),
    ("Ttriple", Param::TTriple),
    ("ptriple", Param::PTriple),
    ("rhocrit", Param::RhomassCritical),
    ("Tmin", Param::TMin),
    ("Tmax", Param::TMax),
    ("pmax", Param::PMax),
    ("pmin", Param::PMin),
    ("molemass", Param::MolarMass),
    ("molarmass", Param::MolarMass),
    ("A", Param::SpeedSound),
    ("I", Param::SurfaceTension),
];

/// Upstream `index_map_insert` stores each name plus its ASCII-uppercased
/// form, so a candidate matches a stored name iff it is byte-equal, or it
/// contains no lowercase ASCII and equals the name case-insensitively.
fn name_matches(candidate: &str, stored: &str) -> bool {
    candidate == stored
        || (!candidate.bytes().any(|b| b.is_ascii_lowercase())
            && candidate.eq_ignore_ascii_case(stored))
}

impl Param {
    /// Every parameter, in upstream enum order (indices 1..=85).
    #[rustfmt::skip]
    pub const ALL: [Param; 85] = [
        Param::GasConstant, Param::MolarMass, Param::AcentricFactor, Param::RhomolarReducing,
        Param::RhomolarCritical, Param::TReducing, Param::TCritical, Param::RhomassReducing,
        Param::RhomassCritical, Param::PCritical, Param::PReducing, Param::TTriple,
        Param::PTriple, Param::TMin, Param::TMax, Param::PMax, Param::PMin, Param::DipoleMoment,
        Param::T, Param::P, Param::Q, Param::Qmass, Param::Tau, Param::Delta,
        Param::Dmolar, Param::Hmolar, Param::Smolar, Param::Cpmolar, Param::Cp0molar,
        Param::Cvmolar, Param::Umolar, Param::Gmolar, Param::Helmholtzmolar,
        Param::HmolarResidual, Param::SmolarResidual, Param::GmolarResidual,
        Param::HmolarIdealgas, Param::SmolarIdealgas, Param::UmolarIdealgas,
        Param::Dmass, Param::Hmass, Param::Smass, Param::Cpmass, Param::Cp0mass, Param::Cvmass,
        Param::Umass, Param::Gmass, Param::Helmholtzmass, Param::HmassIdealgas,
        Param::SmassIdealgas, Param::UmassIdealgas,
        Param::Viscosity, Param::Conductivity, Param::SurfaceTension, Param::Prandtl,
        Param::SpeedSound, Param::IsothermalCompressibility, Param::IsobaricExpansionCoefficient,
        Param::IsentropicExpansionCoefficient, Param::FundamentalDerivativeOfGasDynamics,
        Param::Alphar, Param::DalpharDtauConstdelta, Param::DalpharDdeltaConsttau,
        Param::Alpha0, Param::Dalpha0DtauConstdelta, Param::Dalpha0DdeltaConsttau,
        Param::D2alpha0Ddelta2Consttau, Param::D3alpha0Ddelta3Consttau,
        Param::Bvirial, Param::Cvirial, Param::DBvirialDT, Param::DCvirialDT, Param::Z, Param::PIP,
        Param::FractionMin, Param::FractionMax, Param::TFreeze,
        Param::GWP20, Param::GWP100, Param::GWP500, Param::FH, Param::HH, Param::PH, Param::ODP,
        Param::Phase,
    ];

    /// Upstream integer value of this parameter.
    pub fn index(self) -> i32 {
        self as i32
    }

    /// Parameter for an upstream integer index.
    pub fn from_index(index: i32) -> Option<Param> {
        Param::ALL.iter().copied().find(|p| p.index() == index)
    }

    /// Upstream `is_valid_parameter` / `get_parameter_index`: resolves a short
    /// name or alias (exact case, or all-uppercase form).
    pub fn parse(name: &str) -> Option<Param> {
        PARAM_INFO
            .iter()
            .find(|pi| name_matches(name, pi.short))
            .map(|pi| pi.param)
            .or_else(|| {
                PARAM_ALIASES
                    .iter()
                    .find(|(alias, _)| name_matches(name, alias))
                    .map(|&(_, p)| p)
            })
    }

    fn info(self) -> &'static ParamInfo {
        PARAM_INFO
            .iter()
            .find(|pi| pi.param == self)
            .expect("every Param has a PARAM_INFO row (golden-verified)")
    }

    /// Upstream `get_parameter_information(key, "short")`.
    pub fn short_name(self) -> &'static str {
        self.info().short
    }

    /// Upstream `get_parameter_information(key, "IO")`: `"IO"` or `"O"`.
    pub fn io(self) -> &'static str {
        self.info().io
    }

    /// Upstream `get_parameter_information(key, "units")`.
    pub fn units(self) -> &'static str {
        self.info().units
    }

    /// Upstream `get_parameter_information(key, "long")`.
    pub fn long_desc(self) -> &'static str {
        self.info().long
    }

    /// Upstream `is_trivial_parameter`.
    pub fn is_trivial(self) -> bool {
        self.info().trivial
    }
}

/// Input pair for a state update, mirroring upstream `enum input_pairs`.
/// In each pair the input keys are sorted alphabetically.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum InputPair {
    QT = 1,
    QmassT,
    PQ,
    PQmass,
    QSmolar,
    QmassSmolar,
    QSmass,
    QmassSmass,
    HmolarQ,
    HmolarQmass,
    HmassQ,
    HmassQmass,
    DmolarQ,
    DmolarQmass,
    DmassQ,
    DmassQmass,
    PT,
    DmassT,
    DmolarT,
    HmolarT,
    HmassT,
    SmolarT,
    SmassT,
    TUmolar,
    TUmass,
    DmassP,
    DmolarP,
    HmassP,
    HmolarP,
    PSmass,
    PSmolar,
    PUmass,
    PUmolar,
    HmassSmass,
    HmolarSmolar,
    SmassUmass,
    SmolarUmolar,
    DmassHmass,
    DmolarHmolar,
    DmassSmass,
    DmolarSmolar,
    DmassUmass,
    DmolarUmolar,
}

pub struct InputPairInfo {
    pub pair: InputPair,
    pub short: &'static str,
    pub long: &'static str,
}

/// Transcribed from `input_pair_list` in `src/DataStructures.cpp` @ v8.0.0,
/// in upstream table order. Short descriptions are not unique ("QS_INPUTS",
/// "QmassS_INPUTS", "HQ_INPUTS", "HQmass_INPUTS" each appear twice);
/// upstream's `std::map::emplace` keeps the first insertion, which the
/// first-match scan in [`InputPair::parse`] reproduces.
#[rustfmt::skip]
pub const INPUT_PAIR_INFO: &[InputPairInfo] = &[
    InputPairInfo { pair: InputPair::QT, short: "QT_INPUTS", long: "Molar quality, Temperature in K" },
    InputPairInfo { pair: InputPair::QmassT, short: "QmassT_INPUTS", long: "Mass-basis quality, Temperature in K" },
    InputPairInfo { pair: InputPair::QSmolar, short: "QS_INPUTS", long: "Molar quality, Entropy in J/mol/K" },
    InputPairInfo { pair: InputPair::QmassSmolar, short: "QmassS_INPUTS", long: "Mass-basis quality, Entropy in J/mol/K" },
    InputPairInfo { pair: InputPair::QSmass, short: "QS_INPUTS", long: "Molar quality, Entropy in J/kg/K" },
    InputPairInfo { pair: InputPair::QmassSmass, short: "QmassS_INPUTS", long: "Mass-basis quality, Entropy in J/kg/K" },
    InputPairInfo { pair: InputPair::HmolarQ, short: "HQ_INPUTS", long: "Enthalpy in J/mol, Molar quality" },
    InputPairInfo { pair: InputPair::HmolarQmass, short: "HQmass_INPUTS", long: "Enthalpy in J/mol, Mass-basis quality" },
    InputPairInfo { pair: InputPair::HmassQ, short: "HQ_INPUTS", long: "Enthalpy in J/kg, Molar quality" },
    InputPairInfo { pair: InputPair::HmassQmass, short: "HQmass_INPUTS", long: "Enthalpy in J/kg, Mass-basis quality" },
    InputPairInfo { pair: InputPair::DmassQ, short: "DmassQ_INPUTS", long: "Molar density kg/m^3, Molar quality" },
    InputPairInfo { pair: InputPair::DmassQmass, short: "DmassQmass_INPUTS", long: "Mass density kg/m^3, Mass-basis quality" },
    InputPairInfo { pair: InputPair::DmolarQ, short: "DmolarQ_INPUTS", long: "Molar density in mol/m^3, Molar quality" },
    InputPairInfo { pair: InputPair::DmolarQmass, short: "DmolarQmass_INPUTS", long: "Molar density in mol/m^3, Mass-basis quality" },
    InputPairInfo { pair: InputPair::PQ, short: "PQ_INPUTS", long: "Pressure in Pa, Molar quality" },
    InputPairInfo { pair: InputPair::PQmass, short: "PQmass_INPUTS", long: "Pressure in Pa, Mass-basis quality" },
    InputPairInfo { pair: InputPair::PT, short: "PT_INPUTS", long: "Pressure in Pa, Temperature in K" },
    InputPairInfo { pair: InputPair::DmassT, short: "DmassT_INPUTS", long: "Mass density in kg/m^3, Temperature in K" },
    InputPairInfo { pair: InputPair::DmolarT, short: "DmolarT_INPUTS", long: "Molar density in mol/m^3, Temperature in K" },
    InputPairInfo { pair: InputPair::HmassT, short: "HmassT_INPUTS", long: "Enthalpy in J/kg, Temperature in K" },
    InputPairInfo { pair: InputPair::HmolarT, short: "HmolarT_INPUTS", long: "Enthalpy in J/mol, Temperature in K" },
    InputPairInfo { pair: InputPair::SmassT, short: "SmassT_INPUTS", long: "Entropy in J/kg/K, Temperature in K" },
    InputPairInfo { pair: InputPair::SmolarT, short: "SmolarT_INPUTS", long: "Entropy in J/mol/K, Temperature in K" },
    InputPairInfo { pair: InputPair::TUmass, short: "TUmass_INPUTS", long: "Temperature in K, Internal energy in J/kg" },
    InputPairInfo { pair: InputPair::TUmolar, short: "TUmolar_INPUTS", long: "Temperature in K, Internal energy in J/mol" },
    InputPairInfo { pair: InputPair::DmassP, short: "DmassP_INPUTS", long: "Mass density in kg/m^3, Pressure in Pa" },
    InputPairInfo { pair: InputPair::DmolarP, short: "DmolarP_INPUTS", long: "Molar density in mol/m^3, Pressure in Pa" },
    InputPairInfo { pair: InputPair::HmassP, short: "HmassP_INPUTS", long: "Enthalpy in J/kg, Pressure in Pa" },
    InputPairInfo { pair: InputPair::HmolarP, short: "HmolarP_INPUTS", long: "Enthalpy in J/mol, Pressure in Pa" },
    InputPairInfo { pair: InputPair::PSmass, short: "PSmass_INPUTS", long: "Pressure in Pa, Entropy in J/kg/K" },
    InputPairInfo { pair: InputPair::PSmolar, short: "PSmolar_INPUTS", long: "Pressure in Pa, Entropy in J/mol/K " },
    InputPairInfo { pair: InputPair::PUmass, short: "PUmass_INPUTS", long: "Pressure in Pa, Internal energy in J/kg" },
    InputPairInfo { pair: InputPair::PUmolar, short: "PUmolar_INPUTS", long: "Pressure in Pa, Internal energy in J/mol" },
    InputPairInfo { pair: InputPair::DmassHmass, short: "DmassHmass_INPUTS", long: "Mass density in kg/m^3, Enthalpy in J/kg" },
    InputPairInfo { pair: InputPair::DmolarHmolar, short: "DmolarHmolar_INPUTS", long: "Molar density in mol/m^3, Enthalpy in J/mol" },
    InputPairInfo { pair: InputPair::DmassSmass, short: "DmassSmass_INPUTS", long: "Mass density in kg/m^3, Entropy in J/kg/K" },
    InputPairInfo { pair: InputPair::DmolarSmolar, short: "DmolarSmolar_INPUTS", long: "Molar density in mol/m^3, Entropy in J/mol/K" },
    InputPairInfo { pair: InputPair::DmassUmass, short: "DmassUmass_INPUTS", long: "Mass density in kg/m^3, Internal energy in J/kg" },
    InputPairInfo { pair: InputPair::DmolarUmolar, short: "DmolarUmolar_INPUTS", long: "Molar density in mol/m^3, Internal energy in J/mol" },
    InputPairInfo { pair: InputPair::HmassSmass, short: "HmassSmass_INPUTS", long: "Enthalpy in J/kg, Entropy in J/kg/K" },
    InputPairInfo { pair: InputPair::HmolarSmolar, short: "HmolarSmolar_INPUTS", long: "Enthalpy in J/mol, Entropy in J/mol/K" },
    InputPairInfo { pair: InputPair::SmassUmass, short: "SmassUmass_INPUTS", long: "Entropy in J/kg/K, Internal energy in J/kg" },
    InputPairInfo { pair: InputPair::SmolarUmolar, short: "SmolarUmolar_INPUTS", long: "Entropy in J/mol/K, Internal energy in J/mol" },
];

impl InputPair {
    /// Upstream `get_input_pair_index` (case sensitive; first table match
    /// wins for the duplicated short descriptions).
    pub fn parse(name: &str) -> Option<InputPair> {
        INPUT_PAIR_INFO
            .iter()
            .find(|i| i.short == name)
            .map(|i| i.pair)
    }

    /// Upstream `get_input_pair_short_desc`.
    pub fn short_desc(self) -> &'static str {
        self.info().short
    }

    /// Upstream `get_input_pair_long_desc`.
    pub fn long_desc(self) -> &'static str {
        self.info().long
    }

    fn info(self) -> &'static InputPairInfo {
        INPUT_PAIR_INFO
            .iter()
            .find(|i| i.pair == self)
            .expect("every InputPair has an INPUT_PAIR_INFO row")
    }

    /// Upstream `split_input_pair`: the two parameters forming this pair.
    /// Total on the Rust enum (upstream throws only for its INVALID sentinel,
    /// which is not ported).
    #[rustfmt::skip]
    pub fn split(self) -> (Param, Param) {
        match self {
            InputPair::QT => (Param::Q, Param::T),
            InputPair::QmassT => (Param::Qmass, Param::T),
            InputPair::QSmolar => (Param::Q, Param::Smolar),
            InputPair::QmassSmolar => (Param::Qmass, Param::Smolar),
            InputPair::QSmass => (Param::Q, Param::Smass),
            InputPair::QmassSmass => (Param::Qmass, Param::Smass),
            InputPair::HmolarQ => (Param::Hmolar, Param::Q),
            InputPair::HmolarQmass => (Param::Hmolar, Param::Qmass),
            InputPair::HmassQ => (Param::Hmass, Param::Q),
            InputPair::HmassQmass => (Param::Hmass, Param::Qmass),
            InputPair::PQ => (Param::P, Param::Q),
            InputPair::PQmass => (Param::P, Param::Qmass),
            InputPair::PT => (Param::P, Param::T),
            InputPair::DmassT => (Param::Dmass, Param::T),
            InputPair::DmolarT => (Param::Dmolar, Param::T),
            InputPair::HmassT => (Param::Hmass, Param::T),
            InputPair::HmolarT => (Param::Hmolar, Param::T),
            InputPair::SmassT => (Param::Smass, Param::T),
            InputPair::SmolarT => (Param::Smolar, Param::T),
            InputPair::TUmass => (Param::T, Param::Umass),
            InputPair::TUmolar => (Param::T, Param::Umolar),
            InputPair::DmassP => (Param::Dmass, Param::P),
            InputPair::DmolarP => (Param::Dmolar, Param::P),
            InputPair::DmassQ => (Param::Dmass, Param::Q),
            InputPair::DmassQmass => (Param::Dmass, Param::Qmass),
            InputPair::DmolarQ => (Param::Dmolar, Param::Q),
            InputPair::DmolarQmass => (Param::Dmolar, Param::Qmass),
            InputPair::HmassP => (Param::Hmass, Param::P),
            InputPair::HmolarP => (Param::Hmolar, Param::P),
            InputPair::PSmass => (Param::P, Param::Smass),
            InputPair::PSmolar => (Param::P, Param::Smolar),
            InputPair::PUmass => (Param::P, Param::Umass),
            InputPair::PUmolar => (Param::P, Param::Umolar),
            InputPair::DmassHmass => (Param::Dmass, Param::Hmass),
            InputPair::DmolarHmolar => (Param::Dmolar, Param::Hmolar),
            InputPair::DmassSmass => (Param::Dmass, Param::Smass),
            InputPair::DmolarSmolar => (Param::Dmolar, Param::Smolar),
            InputPair::DmassUmass => (Param::Dmass, Param::Umass),
            InputPair::DmolarUmolar => (Param::Dmolar, Param::Umolar),
            InputPair::HmassSmass => (Param::Hmass, Param::Smass),
            InputPair::HmolarSmolar => (Param::Hmolar, Param::Smolar),
            InputPair::SmassUmass => (Param::Smass, Param::Umass),
            InputPair::SmolarUmolar => (Param::Smolar, Param::Umolar),
        }
    }

    /// Upstream `is_Qmass_pair`.
    pub fn is_qmass_pair(self) -> bool {
        matches!(
            self,
            InputPair::QmassT
                | InputPair::PQmass
                | InputPair::QmassSmolar
                | InputPair::QmassSmass
                | InputPair::HmolarQmass
                | InputPair::HmassQmass
                | InputPair::DmolarQmass
                | InputPair::DmassQmass
        )
    }
}

/// Upstream `generate_update_pair`: resolves two (parameter, value) inputs to
/// the canonical input pair, reordering the values to the pair's convention.
/// Returns `None` where upstream returns `INPUT_PAIR_INVALID`. Ported in the
/// exact upstream match order.
pub fn generate_update_pair(
    key1: Param,
    value1: f64,
    key2: Param,
    value2: f64,
) -> Option<(InputPair, f64, f64)> {
    #[rustfmt::skip]
    const MATCH_ORDER: &[(Param, Param, InputPair)] = &[
        (Param::Q, Param::T, InputPair::QT),
        (Param::Qmass, Param::T, InputPair::QmassT),
        (Param::P, Param::Q, InputPair::PQ),
        (Param::P, Param::Qmass, InputPair::PQmass),
        (Param::P, Param::T, InputPair::PT),
        (Param::Dmolar, Param::T, InputPair::DmolarT),
        (Param::Dmass, Param::T, InputPair::DmassT),
        (Param::Hmolar, Param::T, InputPair::HmolarT),
        (Param::Hmass, Param::T, InputPair::HmassT),
        (Param::Smolar, Param::T, InputPair::SmolarT),
        (Param::Smass, Param::T, InputPair::SmassT),
        (Param::T, Param::Umolar, InputPair::TUmolar),
        (Param::T, Param::Umass, InputPair::TUmass),
        (Param::Dmass, Param::Hmass, InputPair::DmassHmass),
        (Param::Dmolar, Param::Hmolar, InputPair::DmolarHmolar),
        (Param::Dmass, Param::Smass, InputPair::DmassSmass),
        (Param::Dmolar, Param::Smolar, InputPair::DmolarSmolar),
        (Param::Dmass, Param::Umass, InputPair::DmassUmass),
        (Param::Dmolar, Param::Umolar, InputPair::DmolarUmolar),
        (Param::Dmass, Param::P, InputPair::DmassP),
        (Param::Dmolar, Param::P, InputPair::DmolarP),
        (Param::Dmass, Param::Q, InputPair::DmassQ),
        (Param::Dmass, Param::Qmass, InputPair::DmassQmass),
        (Param::Dmolar, Param::Q, InputPair::DmolarQ),
        (Param::Dmolar, Param::Qmass, InputPair::DmolarQmass),
        (Param::Hmass, Param::P, InputPair::HmassP),
        (Param::Hmolar, Param::P, InputPair::HmolarP),
        (Param::P, Param::Smass, InputPair::PSmass),
        (Param::P, Param::Smolar, InputPair::PSmolar),
        (Param::P, Param::Umass, InputPair::PUmass),
        (Param::P, Param::Umolar, InputPair::PUmolar),
        (Param::Hmass, Param::Smass, InputPair::HmassSmass),
        (Param::Hmolar, Param::Smolar, InputPair::HmolarSmolar),
        (Param::Smass, Param::Umass, InputPair::SmassUmass),
        (Param::Smolar, Param::Umolar, InputPair::SmolarUmolar),
    ];
    for &(x1, x2, pair) in MATCH_ORDER {
        // Upstream match_pair: swap = !(key1 == x1)
        if (key1 == x1 && key2 == x2) || (key2 == x1 && key1 == x2) {
            let swap = key1 != x1;
            let (out1, out2) = if swap {
                (value2, value1)
            } else {
                (value1, value2)
            };
            return Some((pair, out1, out2));
        }
    }
    None
}

/// Fluid phase, mirroring upstream `enum phases`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum Phase {
    Liquid = 0,
    Supercritical,
    SupercriticalGas,
    SupercriticalLiquid,
    CriticalPoint,
    Gas,
    Twophase,
    Unknown,
    NotImposed,
}

/// Transcribed from `phase_info_list` in `src/DataStructures.cpp` @ v8.0.0.
#[rustfmt::skip]
pub const PHASE_INFO: &[(Phase, &str, &str)] = &[
    (Phase::Liquid, "phase_liquid", ""),
    (Phase::Gas, "phase_gas", ""),
    (Phase::Twophase, "phase_twophase", ""),
    (Phase::Supercritical, "phase_supercritical", ""),
    (Phase::SupercriticalGas, "phase_supercritical_gas", "p < pc, T > Tc"),
    (Phase::SupercriticalLiquid, "phase_supercritical_liquid", "p > pc, T < Tc"),
    (Phase::CriticalPoint, "phase_critical_point", "p = pc, T = Tc"),
    (Phase::Unknown, "phase_unknown", ""),
    (Phase::NotImposed, "phase_not_imposed", ""),
];

impl Phase {
    /// Every phase, in upstream enum order (indices 0..=8).
    pub const ALL: [Phase; 9] = [
        Phase::Liquid,
        Phase::Supercritical,
        Phase::SupercriticalGas,
        Phase::SupercriticalLiquid,
        Phase::CriticalPoint,
        Phase::Gas,
        Phase::Twophase,
        Phase::Unknown,
        Phase::NotImposed,
    ];

    /// Upstream integer value of this phase.
    pub fn index(self) -> i32 {
        self as i32
    }

    /// Upstream `is_valid_phase` / `get_phase_index` (case sensitive).
    pub fn parse(name: &str) -> Option<Phase> {
        PHASE_INFO
            .iter()
            .find(|(_, s, _)| *s == name)
            .map(|&(p, _, _)| p)
    }

    /// Upstream `get_phase_short_desc`.
    pub fn short_desc(self) -> &'static str {
        PHASE_INFO
            .iter()
            .find(|(p, _, _)| *p == self)
            .map(|&(_, s, _)| s)
            .expect("every Phase has a PHASE_INFO row")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Source-truth tests (no Python oracle exposes these): transcribed
    // semantics from src/DataStructures.cpp @ v8.0.0.

    #[test]
    fn param_case_rule_matches_upstream_index_map() {
        assert_eq!(Param::parse("Dmolar"), Some(Param::Dmolar));
        assert_eq!(Param::parse("DMOLAR"), Some(Param::Dmolar));
        assert_eq!(Param::parse("dmolar"), None); // only exact + uppercase stored
        assert_eq!(Param::parse("DMoLAR"), None);
        assert_eq!(Param::parse("viscosity"), Some(Param::Viscosity));
        assert_eq!(Param::parse("VISCOSITY"), Some(Param::Viscosity));
        assert_eq!(Param::parse("T"), Some(Param::T));
        assert_eq!(Param::parse("D"), Some(Param::Dmass)); // alias
        assert_eq!(Param::parse("Tcrit"), Some(Param::TCritical));
        assert_eq!(Param::parse("TCRIT"), Some(Param::TCritical));
        assert_eq!(Param::parse("nonsense"), None);
    }

    #[test]
    fn input_pair_parse_is_first_match_like_upstream_emplace() {
        // "QS_INPUTS" and "HQ_INPUTS" are duplicated in the upstream table;
        // std::map::emplace keeps the first insertion.
        assert_eq!(InputPair::parse("QS_INPUTS"), Some(InputPair::QSmolar));
        assert_eq!(
            InputPair::parse("QmassS_INPUTS"),
            Some(InputPair::QmassSmolar)
        );
        assert_eq!(InputPair::parse("HQ_INPUTS"), Some(InputPair::HmolarQ));
        assert_eq!(
            InputPair::parse("HQmass_INPUTS"),
            Some(InputPair::HmolarQmass)
        );
        assert_eq!(InputPair::parse("PT_INPUTS"), Some(InputPair::PT));
        assert_eq!(InputPair::parse("pt_inputs"), None); // case sensitive
    }

    #[test]
    fn split_covers_every_pair_and_matches_short_desc_ordering() {
        // The pair name encodes its split: e.g. DmassT -> (Dmass, T).
        for info in INPUT_PAIR_INFO {
            let (p1, p2) = info.pair.split();
            let expected_prefix = format!(
                "{}{}",
                p1.short_name().replace("_", ""),
                p2.short_name().replace("_", "")
            );
            let short = info.short.replace("_INPUTS", "").replace("_", "");
            // Upstream uses "QS"/"HQ" (molar dropped) for the molar-quality
            // entropy/enthalpy pairs; accept the documented exceptions.
            let exceptions = [
                "QSmolar",
                "QSmass",
                "QmassSmolar",
                "QmassSmass",
                "HmolarQ",
                "HmolarQmass",
                "HmassQ",
                "HmassQmass",
            ];
            if exceptions.contains(&format!("{:?}", info.pair).as_str()) {
                continue;
            }
            assert_eq!(short, expected_prefix, "pair {:?}", info.pair);
        }
    }

    #[test]
    fn generate_update_pair_swaps_values_like_upstream() {
        // match_pair(key1=T, x1=P) => swap=true => out1 gets the P value.
        assert_eq!(
            generate_update_pair(Param::T, 300.0, Param::P, 101325.0),
            Some((InputPair::PT, 101325.0, 300.0))
        );
        assert_eq!(
            generate_update_pair(Param::P, 101325.0, Param::T, 300.0),
            Some((InputPair::PT, 101325.0, 300.0))
        );
        assert_eq!(
            generate_update_pair(Param::T, 300.0, Param::Q, 0.5),
            Some((InputPair::QT, 0.5, 300.0))
        );
        assert_eq!(
            generate_update_pair(Param::T, 300.0, Param::Viscosity, 1.0),
            None
        );
        assert_eq!(generate_update_pair(Param::T, 300.0, Param::T, 400.0), None);
    }

    #[test]
    fn qmass_pairs_match_upstream_list() {
        let expected = [
            InputPair::QmassT,
            InputPair::PQmass,
            InputPair::QmassSmolar,
            InputPair::QmassSmass,
            InputPair::HmolarQmass,
            InputPair::HmassQmass,
            InputPair::DmolarQmass,
            InputPair::DmassQmass,
        ];
        for info in INPUT_PAIR_INFO {
            assert_eq!(
                info.pair.is_qmass_pair(),
                expected.contains(&info.pair),
                "pair {:?}",
                info.pair
            );
        }
    }

    #[test]
    fn phase_parse_round_trips() {
        for phase in Phase::ALL {
            assert_eq!(Phase::parse(phase.short_desc()), Some(phase));
        }
        assert_eq!(Phase::parse("PHASE_LIQUID"), None); // case sensitive
    }

    #[test]
    fn psmolar_long_desc_preserves_upstream_trailing_space() {
        assert_eq!(
            InputPair::PSmolar.long_desc(),
            "Pressure in Pa, Entropy in J/mol/K "
        );
    }
}
