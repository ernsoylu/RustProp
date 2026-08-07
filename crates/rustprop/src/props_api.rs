//! `PropsSI` string API (PLAN.md Phase 5) — port of upstream
//! `PropsSI`/`_PropsSImulti`/`_PropsSI_outputs` (src/CoolProp.cpp) restricted
//! to a single input point and a single output:
//!
//! - backend prefix parsing (`extract_backend`: `"HEOS::"`, `"IF97::"`, bare
//!   fluid -> `"?"` -> HEOS default);
//! - fluid resolution exactly as upstream `JSONFluidLibrary` registers keys:
//!   CAS, canonical name, every alias, and `upper(alias)` (EES compat), with
//!   later registrations replacing earlier ones;
//! - input-pair resolution via the golden-verified `generate_update_pair`,
//!   the trivial-output route (no state update; reachable with invalid or
//!   empty input names), the outputs-are-inputs echo route, and the
//!   mass-to-molar input conversion (`AbstractState::mass_to_molar_inputs`);
//! - `keyed_output` with upstream's mass-basis conversions and two-phase
//!   error conditions.
//!
//! Upstream returns `_HUGE` and sets a global error string; this port
//! returns the matching `Err` variant (error *conditions* are the fidelity
//! target, not message strings). Input pairs whose flashes are not yet
//! ported return `Error::NotImplemented` loudly.

use rustprop_core::params::{InputPair, Param, generate_update_pair};
use rustprop_core::{Error, Result};

#[cfg(feature = "heos")]
use rustprop_core::fluid::FluidData;
#[cfg(feature = "heos")]
use rustprop_heos::flash_pt::PtFlash;
#[cfg(feature = "heos")]
use rustprop_heos::flash_px::HeosState;

/// Upstream `extract_backend`: split `"BACKEND::FLUID"`; no prefix yields
/// backend `"?"` (the caller defaults it to HEOS).
fn extract_backend(fluid_string: &str) -> (&str, &str) {
    match fluid_string.find("::") {
        Some(i) => (&fluid_string[..i], &fluid_string[i + 2..]),
        None => ("?", fluid_string),
    }
}

/// `PropsSI(Output, Name1, Prop1, Name2, Prop2, FluidName)`.
pub fn props_si(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid_name: &str,
) -> Result<f64> {
    let (backend, fluid) = extract_backend(fluid_name);
    if fluid.contains('&') {
        return Err(Error::NotImplemented("mixtures are not ported yet".into()));
    }
    match backend {
        "IF97" => if97_route(output, name1, prop1, name2, prop2, fluid),
        "HEOS" | "?" => heos_route(output, name1, prop1, name2, prop2, fluid),
        "SRK" => cubic_route(
            CubicRouteKind::Srk,
            output,
            name1,
            prop1,
            name2,
            prop2,
            fluid,
        ),
        "PR" => cubic_route(
            CubicRouteKind::PengRobinson,
            output,
            name1,
            prop1,
            name2,
            prop2,
            fluid,
        ),
        other => Err(Error::Value(format!(
            "Invalid backend name [{other}] to factory function"
        ))),
    }
}

#[cfg(feature = "if97")]
fn if97_route(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid: &str,
) -> Result<f64> {
    if fluid != "Water" {
        return Err(Error::Value(format!(
            "IF97 backend supports only Water, got [{fluid}]"
        )));
    }
    let out = Param::parse(output).ok_or_else(|| {
        Error::Value(format!(
            "Output parameter parsing failed; error: Output string is invalid [{output}]"
        ))
    })?;
    let n1 = Param::parse(name1)
        .ok_or_else(|| Error::Value(format!("Input pair parsing failed for Name1: \"{name1}\"")))?;
    let n2 = Param::parse(name2)
        .ok_or_else(|| Error::Value(format!("Input pair parsing failed for Name2: \"{name2}\"")))?;
    crate::if97_api::props(out, n1, prop1, n2, prop2)
}

#[cfg(not(feature = "if97"))]
fn if97_route(_: &str, _: &str, _: f64, _: &str, _: f64, _: &str) -> Result<f64> {
    Err(Error::NotImplemented(
        "the `if97` feature is not enabled".into(),
    ))
}

#[cfg(not(feature = "heos"))]
fn heos_route(_: &str, _: &str, _: f64, _: &str, _: f64, _: &str) -> Result<f64> {
    Err(Error::NotImplemented(
        "the `heos` feature is not enabled".into(),
    ))
}

// ---------------------------------------------------------------------------
// Fluid registry (upstream JSONFluidLibrary string_to_index_map)
// ---------------------------------------------------------------------------

/// Resolve a fluid key exactly as upstream registers them: CAS, name, alias,
/// `upper(alias)`; later fluids replace earlier claims of the same key.
#[cfg(feature = "heos")]
pub fn resolve_fluid(key: &str) -> Result<&'static FluidData> {
    use std::collections::HashMap;
    use std::sync::OnceLock;
    static MAP: OnceLock<HashMap<String, &'static FluidData>> = OnceLock::new();
    let map = MAP.get_or_init(|| {
        let mut m: HashMap<String, &'static FluidData> = HashMap::new();
        for (_name, data) in rustprop_data::fluids::all() {
            m.insert(data.cas.to_string(), data);
            m.insert(data.name.to_string(), data);
            for alias in data.aliases {
                m.insert((*alias).to_string(), data);
                m.insert(alias.to_uppercase(), data);
            }
        }
        m
    });
    map.get(key).copied().ok_or_else(|| {
        Error::Value(format!(
            "key [{key}] was not found in string_to_index_map in JSONFluidLibrary"
        ))
    })
}

// ---------------------------------------------------------------------------
// HEOS route (upstream _PropsSImulti + _PropsSI_outputs, single point)
// ---------------------------------------------------------------------------

#[cfg(feature = "heos")]
fn heos_route(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid: &str,
) -> Result<f64> {
    // Multi-output strings ('&'-joined) fail upstream's output parsing with
    // "Output string is invalid" — they fall through to Param::parse below.
    let data = resolve_fluid(fluid)?;

    let out = Param::parse(output).ok_or_else(|| {
        Error::Value(format!(
            "Output parameter parsing failed; error: Output string is invalid [{output}]"
        ))
    })?;

    // Upstream: the pair stays INVALID when either NAME fails to parse (the
    // trivial route stays reachable); a pair of VALID params that is not a
    // valid combination throws immediately.
    // Upstream `generate_update_pair` RETURNS INPUT_PAIR_INVALID for an
    // unknown combination (it never throws) — PropsSImulti then errors only
    // if the outputs need a state update. Unknown parameter names land in
    // the same place.
    let keys = (Param::parse(name1), Param::parse(name2));
    let pair = match keys {
        (Some(k1), Some(k2)) => generate_update_pair(k1, prop1, k2, prop2),
        _ => None,
    };

    let flash = fluid_flash(data);
    let Some((pair, v1, v2)) = pair else {
        // If all outputs are trivial, never do a state update.
        if out.is_trivial() {
            return trivial_output(flash, data, out);
        }
        return Err(Error::Value(
            "Input pair variable is invalid and output(s) are non-trivial; cannot do state update"
                .into(),
        ));
    };
    if out.is_trivial() {
        return trivial_output(flash, data, out);
    }
    // If all outputs are also inputs, never do a state update.
    let (p1, p2) = pair.split();
    if out == p1 {
        return Ok(v1);
    }
    if out == p2 {
        return Ok(v2);
    }

    let state = update(flash, pair, v1, v2)?;
    keyed_output(flash, data, &state, out)
}

/// One cached `PtFlash` per fluid (upstream keeps backend instances alive
/// per factory call; the superancillary inverses and calorics are per-fluid
/// caches worth sharing across `props_si` calls).
#[cfg(feature = "heos")]
fn fluid_flash(data: &'static FluidData) -> &'static PtFlash {
    use std::collections::HashMap;
    use std::sync::{Mutex, OnceLock};
    static FLASHES: OnceLock<Mutex<HashMap<usize, &'static PtFlash>>> = OnceLock::new();
    let m = FLASHES.get_or_init(|| Mutex::new(HashMap::new()));
    let key = std::ptr::from_ref(data) as usize;
    let mut guard = m.lock().expect("fluid flash cache poisoned");
    guard
        .entry(key)
        .or_insert_with(|| Box::leak(Box::new(PtFlash::new(data))))
}

/// Upstream `HelmholtzEOSMixtureBackend::update`: mass-to-molar conversion,
/// Q-range validation, then the pair dispatch. Pairs whose flashes are not
/// ported yet error loudly.
#[cfg(feature = "heos")]
fn update(flash: &PtFlash, pair: InputPair, v1: f64, v2: f64) -> Result<HeosState> {
    // Upstream `AbstractState::mass_to_molar_inputs`.
    let mm = flash.eos.molar_mass;
    let (pair, v1, v2) = match pair {
        InputPair::QmassT => (InputPair::QT, v1, v2),
        InputPair::PQmass => (InputPair::PQ, v1, v2),
        InputPair::DmassT => (InputPair::DmolarT, v1 / mm, v2),
        InputPair::SmassT => (InputPair::SmolarT, v1 * mm, v2),
        InputPair::DmassP => (InputPair::DmolarP, v1 / mm, v2),
        InputPair::HmassP => (InputPair::HmolarP, v1 * mm, v2),
        InputPair::PSmass => (InputPair::PSmolar, v1, v2 * mm),
        InputPair::HmassSmass => (InputPair::HmolarSmolar, v1 * mm, v2 * mm),
        InputPair::PUmass => (InputPair::PUmolar, v1, v2 * mm),
        InputPair::SmassUmass => (InputPair::SmolarUmolar, v1 * mm, v2 * mm),
        // Pure fluids: Qmass == Qmolar exactly; rewrite to the molar sibling.
        InputPair::DmolarQmass => (InputPair::DmolarQ, v1, v2),
        InputPair::DmassQmass | InputPair::DmassQ => (InputPair::DmolarQ, v1 / mm, v2),
        InputPair::DmassHmass => (InputPair::DmolarHmolar, v1 / mm, v2 * mm),
        InputPair::DmassSmass => (InputPair::DmolarSmolar, v1 / mm, v2 * mm),
        InputPair::DmassUmass => (InputPair::DmolarUmolar, v1 / mm, v2 * mm),
        other => (other, v1, v2),
    };
    // Pseudo-pure fluids: the ported paths are the ones upstream's
    // ancillary machinery serves directly (QT at Q in {0,1}, PQ, PT); the
    // remaining pairs route into superancillary-only machinery and stay
    // loud errors until their classic legacy solvers are ported.
    if flash.fluid().eos.pseudo_pure
        && !matches!(pair, InputPair::PT | InputPair::QT | InputPair::PQ)
    {
        return Err(Error::NotImplemented(format!(
            "input pair {} is not ported yet for pseudo-pure fluids",
            pair.short_desc()
        )));
    }
    match pair {
        InputPair::PT => {
            let (rho, phase) = flash.pt_flash(v2, v1)?;
            Ok(HeosState::SinglePhase {
                t: v2,
                p: v1,
                rhomolar: rho,
                phase,
                q: -1.0,
            })
        }
        InputPair::QT => {
            if !(0.0..=1.0).contains(&v1) {
                return Err(Error::OutOfRange(
                    "Input vapor quality [Q] must be between 0 and 1".into(),
                ));
            }
            flash.qt_state(v2, v1)
        }
        InputPair::PQ => {
            if !(0.0..=1.0).contains(&v2) {
                return Err(Error::OutOfRange(
                    "Input vapor quality [Q] must be between 0 and 1".into(),
                ));
            }
            flash.pq_state(v1, v2)
        }
        InputPair::DmolarT => flash.dmolar_t_state(v1, v2),
        InputPair::SmolarT => flash.smolar_t_state(v1, v2),
        InputPair::DmolarP => flash.dmolar_p_state(v1, v2),
        InputPair::HmolarP => flash.hmolar_p_state(v1, v2),
        InputPair::PSmolar => flash.p_smolar_state(v1, v2),
        InputPair::HmolarSmolar => flash.hmolar_smolar_state(v1, v2),
        InputPair::HmolarT => flash.hmolar_t_state(v1, v2),
        InputPair::TUmolar => flash.umolar_t_state(v2, v1),
        InputPair::PUmolar => flash.p_umolar_state(v1, v2),
        InputPair::DmolarHmolar => flash.dmolar_hmolar_state(v1, v2),
        InputPair::DmolarSmolar => flash.dmolar_smolar_state(v1, v2),
        InputPair::DmolarUmolar => flash.dmolar_umolar_state(v1, v2),
        InputPair::DmolarQ => flash.dmolar_q_state(v1, v2),
        // Upstream's own dead ends: mass_to_molar_inputs never converts
        // these, and the backend update switch throws for them.
        InputPair::HmassT | InputPair::TUmass | InputPair::SmolarUmolar => {
            Err(Error::Value(format!(
                "This pair of inputs [{}] is not yet supported",
                pair.short_desc()
            )))
        }
        other => Err(Error::NotImplemented(format!(
            "input pair {} is not ported yet",
            other.short_desc()
        ))),
    }
}

/// Upstream `AbstractState::keyed_output` for the ported outputs, including
/// the mass-basis conversions and the two-phase error conditions.
#[cfg(feature = "heos")]
fn keyed_output(
    flash: &PtFlash,
    data: &'static FluidData,
    state: &HeosState,
    out: Param,
) -> Result<f64> {
    if out.is_trivial() {
        return trivial_output(flash, data, out);
    }
    let mm = flash.eos.molar_mass;
    // Second-derivative and single-phase-only properties are undefined in
    // the two-phase region (upstream throws).
    let single_phase_rho = |what: &str| -> Result<f64> {
        match state {
            HeosState::SinglePhase { rhomolar, .. } => Ok(*rhomolar),
            HeosState::TwoPhase { .. } => Err(Error::Value(format!(
                "Input is two-phase and {what} is not defined"
            ))),
        }
    };
    Ok(match out {
        Param::T => state.t(),
        Param::P => state.p(),
        Param::Q => state.q(),
        Param::Dmolar => state.rhomolar(),
        Param::Dmass => state.rhomolar() * mm,
        Param::Hmolar => flash.state_hmolar(state),
        Param::Hmass => flash.state_hmolar(state) / mm,
        Param::Smolar => flash.state_smolar(state),
        Param::Smass => flash.state_smolar(state) / mm,
        Param::Umolar => flash.state_umolar(state),
        Param::Umass => flash.state_umolar(state) / mm,
        Param::Cpmolar => flash.eos.cpmolar(state.t(), single_phase_rho("cpmolar")?),
        Param::Cpmass => flash.eos.cpmolar(state.t(), single_phase_rho("cpmass")?) / mm,
        Param::Cvmolar => flash.eos.cvmolar(state.t(), single_phase_rho("cvmolar")?),
        Param::Cvmass => flash.eos.cvmolar(state.t(), single_phase_rho("cvmass")?) / mm,
        Param::SpeedSound => flash
            .eos
            .speed_sound(state.t(), single_phase_rho("speed_sound")?),
        Param::Gmolar => flash
            .eos
            .gibbsmolar(state.t(), single_phase_rho("gibbsmolar")?),
        Param::Gmass => {
            flash
                .eos
                .gibbsmolar(state.t(), single_phase_rho("gibbsmass")?)
                / mm
        }
        Param::Viscosity => {
            // Upstream calc_viscosity evaluates at the state's (T, rhomolar)
            // regardless of phase (two-phase uses the mixture density).
            let v = viscosity_model(data)?;
            rustprop_heos::transport::viscosity(
                &flash.eos,
                data,
                v,
                state.t(),
                state.rhomolar(),
                state.p(),
                Some(&ecs_resolver),
            )?
        }
        Param::Conductivity => {
            let tr = data.transport.as_ref().ok_or_else(|| {
                Error::Value("Thermal conductivity model is not available for this fluid".into())
            })?;
            let c = match &tr.conductivity {
                rustprop_core::fluid::TransportModel::Absent => {
                    return Err(Error::Value(
                        "Thermal conductivity model is not available for this fluid".into(),
                    ));
                }
                rustprop_core::fluid::TransportModel::Unported => {
                    return Err(Error::NotImplemented(
                        "this fluid's conductivity model class is not ported yet".into(),
                    ));
                }
                rustprop_core::fluid::TransportModel::Model(c) => c,
            };
            // The Olchowy-Sengers term needs the fluid's viscosity; pass the
            // model when it is ported (its absence only errors if actually
            // needed). Two-phase states evaluate at the mixture density —
            // upstream has no two-phase guard anywhere in the conductivity
            // path (cp/cv are the raw single-phase formulas).
            let v = viscosity_model(data).ok();
            rustprop_heos::transport::conductivity(
                &flash.eos,
                data,
                c,
                v,
                state.t(),
                state.rhomolar(),
                state.p(),
                Some(&ecs_resolver),
            )?
        }
        Param::SurfaceTension => {
            // Upstream `calc_surface_tension`: two-phase (or critical-point)
            // states only; NotImplemented when the fluid has no curve.
            match state {
                HeosState::TwoPhase { t, .. } => {
                    let st = data.ancillaries.surface_tension.as_ref().ok_or_else(|| {
                        Error::NotImplemented("surface tension curve not provided".into())
                    })?;
                    rustprop_heos::ancillary::surface_tension(st, *t)?
                }
                HeosState::SinglePhase { .. } => {
                    return Err(Error::Value(
                        "surface tension is only defined within the two-phase region;                          Try PQ or QT inputs"
                            .into(),
                    ));
                }
            }
        }
        other => {
            return Err(Error::NotImplemented(format!(
                "output parameter {} is not ported yet",
                other.short_name()
            )));
        }
    })
}

/// ECS reference-fluid resolver: looks the name up in the compiled-in
/// registry (the reference must be feature-enabled alongside the fluid) and
/// hands transport its EOS, document, and transport models.
#[cfg(feature = "heos")]
fn ecs_resolver(name: &str) -> Result<rustprop_heos::transport::EcsRef<'static>> {
    let data = resolve_fluid(name)?;
    let flash = fluid_flash(data);
    let viscosity = match data.transport.as_ref().map(|tr| &tr.viscosity) {
        Some(rustprop_core::fluid::TransportModel::Model(v)) => Some(v),
        _ => None,
    };
    let conductivity = match data.transport.as_ref().map(|tr| &tr.conductivity) {
        Some(rustprop_core::fluid::TransportModel::Model(c)) => Some(c),
        _ => None,
    };
    Ok(rustprop_heos::transport::EcsRef {
        eos: &flash.eos,
        fluid: data,
        viscosity,
        conductivity,
    })
}

/// Resolve the fluid's ported viscosity model (upstream
/// `viscosity_model_provided` + the per-class dispatch).
#[cfg(feature = "heos")]
fn viscosity_model(
    data: &'static FluidData,
) -> Result<&'static rustprop_core::fluid::ViscosityModel> {
    let tr = data
        .transport
        .as_ref()
        .ok_or_else(|| Error::Value("Viscosity model is not available for this fluid".into()))?;
    match &tr.viscosity {
        rustprop_core::fluid::TransportModel::Absent => Err(Error::Value(
            "Viscosity model is not available for this fluid".into(),
        )),
        rustprop_core::fluid::TransportModel::Unported => Err(Error::NotImplemented(
            "this fluid's viscosity model class is not ported yet".into(),
        )),
        rustprop_core::fluid::TransportModel::Model(v) => Ok(v),
    }
}

/// Trivial outputs (no state update). `Ttriple`/`Tmin` follow upstream's
/// runtime semantics (`sat_min_liquid.T`), and the critical parameters are
/// the superancillary's NUMERICAL critical point (`calc_T_critical` et al.
/// prefer `get_Tcrit_num()`/`get_pmax()`/`get_rhocrit_num()`), not the
/// document's `STATES.critical`.
#[cfg(feature = "heos")]
fn trivial_output(flash: &PtFlash, data: &'static FluidData, out: Param) -> Result<f64> {
    Ok(match out {
        Param::TCritical => flash.t_critical(),
        Param::PCritical => flash.p_critical(),
        Param::RhomolarCritical => flash.rhomolar_critical(),
        Param::RhomassCritical => flash.rhomolar_critical() * data.eos.molar_mass,
        Param::TTriple => data.eos.sat_min_liquid.t,
        Param::PTriple => data.eos.sat_min_liquid.p,
        Param::TMin => data.eos.sat_min_liquid.t,
        Param::TMax => data.eos.t_max,
        Param::PMax => data.eos.p_max,
        Param::TReducing => data.eos.reducing.t,
        Param::RhomolarReducing => data.eos.reducing.rhomolar,
        Param::RhomassReducing => data.eos.reducing.rhomolar * data.eos.molar_mass,
        Param::MolarMass => data.eos.molar_mass,
        Param::AcentricFactor => data.eos.acentric,
        Param::GasConstant => data.eos.gas_constant,
        other => {
            return Err(Error::NotImplemented(format!(
                "trivial output parameter {} is not ported yet",
                other.short_name()
            )));
        }
    })
}

/// Backend selector for the cubic route (kept separate from the engine enum
/// so the non-cubics build still compiles the dispatch arm).
#[derive(Clone, Copy)]
enum CubicRouteKind {
    Srk,
    PengRobinson,
}

#[cfg(not(feature = "cubics"))]
#[allow(clippy::too_many_arguments)]
fn cubic_route(
    _: CubicRouteKind,
    _: &str,
    _: &str,
    _: f64,
    _: &str,
    _: f64,
    _: &str,
) -> Result<f64> {
    Err(Error::NotImplemented(
        "the `cubics` feature is not enabled".into(),
    ))
}

/// PropsSI route for `SRK::` / `PR::` (upstream `AbstractCubicBackend`
/// through the string API): PT / QT / PQ flashes, trivial outputs, and
/// upstream's error conditions for everything else.
#[cfg(feature = "cubics")]
#[allow(clippy::too_many_arguments)]
fn cubic_route(
    kind: CubicRouteKind,
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid: &str,
) -> Result<f64> {
    use rustprop_core::params::Phase;
    use rustprop_cubics::{CubicEos, CubicKind, CubicSatState};

    let kind = match kind {
        CubicRouteKind::Srk => CubicKind::Srk,
        CubicRouteKind::PengRobinson => CubicKind::PengRobinson,
    };
    // Upstream CubicsLibrary::get: uppercase name first, then alias (the
    // JSON aliases are stored uppercase already).
    let upper = fluid.to_uppercase();
    let data = rustprop_data::cubics::CUBIC_FLUIDS
        .iter()
        .find(|f| f.name == upper || f.aliases.contains(&upper.as_str()))
        .ok_or_else(|| {
            Error::Value(format!(
                "fluid name [{fluid}] is not the name of a cubic fluid"
            ))
        })?;

    // One engine per (kind, fluid), leaked once (same pattern as the HEOS
    // per-fluid flash cache).
    let eos: &'static CubicEos = {
        use std::collections::HashMap;
        use std::sync::{Mutex, OnceLock};
        static CACHE: OnceLock<Mutex<HashMap<(u8, &'static str), &'static CubicEos>>> =
            OnceLock::new();
        let key = (matches!(kind, CubicKind::PengRobinson) as u8, data.name);
        let mut map = CACHE
            .get_or_init(|| Mutex::new(HashMap::new()))
            .lock()
            .unwrap();
        map.entry(key)
            .or_insert_with(|| Box::leak(Box::new(CubicEos::new(kind, data))))
    };

    let out = Param::parse(output).ok_or_else(|| {
        Error::Value(format!(
            "Output parameter parsing failed; error: Output string is invalid [{output}]"
        ))
    })?;

    // Trivial outputs never need a state (upstream trivial_keyed_output).
    let trivial = |out: Param| -> Option<f64> {
        Some(match out {
            Param::TCritical => data.tc,
            Param::PCritical => data.pc,
            Param::AcentricFactor => data.acentric,
            Param::RhomolarCritical => eos.rhomolar_critical(),
            Param::RhomassCritical => eos.rhomolar_critical() * data.molemass,
            Param::MolarMass => data.molemass,
            Param::GasConstant => eos.gas_constant(),
            _ => return None,
        })
    };
    let keys = (Param::parse(name1), Param::parse(name2));
    let pair = match keys {
        (Some(k1), Some(k2)) => generate_update_pair(k1, prop1, k2, prop2),
        _ => None,
    };
    let Some((pair, v1, v2)) = pair else {
        if let Some(v) = trivial(out) {
            return Ok(v);
        }
        return Err(Error::Value(
            "Input pair variable is invalid and output(s) are non-trivial; cannot do state update"
                .into(),
        ));
    };
    if let Some(v) = trivial(out) {
        return Ok(v);
    }

    // Mass-basis inputs convert exactly as the HEOS route.
    let mm = data.molemass;
    let (pair, v1, v2) = match pair {
        InputPair::QmassT => (InputPair::QT, v1, v2),
        InputPair::PQmass => (InputPair::PQ, v1, v2),
        other => (other, v1, v2),
    };

    enum CubicState {
        Single {
            t: f64,
            rho: f64,
            p: f64,
            phase: Phase,
        },
        TwoPhase(CubicSatState),
    }
    let state = match pair {
        InputPair::PT => {
            let (rho, phase) = eos.pt_flash(v2, v1)?;
            CubicState::Single {
                t: v2,
                rho,
                p: v1,
                phase,
            }
        }
        // Upstream's cubic update applies NO quality-range validation.
        InputPair::QT => CubicState::TwoPhase(eos.qt_flash(v2, v1)?),
        InputPair::PQ => CubicState::TwoPhase(eos.pq_flash(v1, v2)?),
        InputPair::DmolarT | InputPair::DmassT => {
            return Err(Error::NotImplemented(
                "the (D,T) cubic route needs the cubic superancillary (PLAN 7.2)".into(),
            ));
        }
        // Upstream delegates every other pair to the HEOS flash routines,
        // which fail on the fabricated cubic fluid's unset ancillaries.
        _ => {
            return Err(Error::Value("type not set".into()));
        }
    };

    let (t, p, rho, q, phase) = match &state {
        CubicState::Single { t, rho, p, phase } => (*t, *p, *rho, -1.0, *phase),
        CubicState::TwoPhase(s) => (s.t, s.p, s.rhomolar, s.q, Phase::Twophase),
    };
    // Two-phase caloric outputs mix the saturated branches by quality.
    let mix = |f: &dyn Fn(f64, f64) -> f64| -> f64 {
        match &state {
            CubicState::Single { t, rho, .. } => f(*t, *rho),
            CubicState::TwoPhase(s) => q * f(s.t, s.rho_v) + (1.0 - q) * f(s.t, s.rho_l),
        }
    };
    Ok(match out {
        Param::T => t,
        Param::P => p,
        Param::Dmolar => rho,
        Param::Dmass => rho * mm,
        Param::Q => q,
        Param::Phase => f64::from(phase.index()),
        Param::Hmolar => mix(&|t, r| eos.hmolar(t, r)),
        Param::Hmass => mix(&|t, r| eos.hmolar(t, r)) / mm,
        Param::Smolar => mix(&|t, r| eos.smolar(t, r)),
        Param::Smass => mix(&|t, r| eos.smolar(t, r)) / mm,
        Param::Umolar => mix(&|t, r| eos.umolar(t, r)),
        Param::Umass => mix(&|t, r| eos.umolar(t, r)) / mm,
        Param::Cpmolar => match &state {
            CubicState::Single { t, rho, .. } => eos.cpmolar(*t, *rho),
            CubicState::TwoPhase(_) => {
                return Err(Error::Value(
                    "Input is two-phase and cp is not defined".into(),
                ));
            }
        },
        Param::Cvmolar => match &state {
            CubicState::Single { t, rho, .. } => eos.cvmolar(*t, *rho),
            CubicState::TwoPhase(_) => {
                return Err(Error::Value(
                    "Input is two-phase and cv is not defined".into(),
                ));
            }
        },
        Param::SpeedSound => match &state {
            CubicState::Single { t, rho, .. } => eos.speed_sound(*t, *rho),
            CubicState::TwoPhase(_) => {
                return Err(Error::Value(
                    "Input is two-phase and the speed of sound is not defined".into(),
                ));
            }
        },
        other => {
            return Err(Error::NotImplemented(format!(
                "output parameter {} is not ported for the cubic backends yet",
                other.short_name()
            )));
        }
    })
}
