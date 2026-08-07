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
    if fluid.contains('&') && backend != "HEOS" && backend != "?" && backend != "PCSAFT" {
        return Err(Error::NotImplemented("mixtures are not ported yet".into()));
    }
    match backend {
        "IF97" => if97_route(output, name1, prop1, name2, prop2, fluid),
        "HEOS" | "?" => {
            if (fluid.contains('[') && fluid.contains(']')) || fluid.contains('&') {
                heos_mixture_entry(output, name1, prop1, name2, prop2, fluid)
            } else {
                heos_route(output, name1, prop1, name2, prop2, fluid)
            }
        }
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
        "INCOMP" => incomp_route(output, name1, prop1, name2, prop2, fluid),
        "TTSE" | "TTSE&HEOS" => tabular_route(
            TabularScheme::Ttse,
            output,
            name1,
            prop1,
            name2,
            prop2,
            fluid,
        ),
        "BICUBIC" | "BICUBIC&HEOS" => tabular_route(
            TabularScheme::Bicubic,
            output,
            name1,
            prop1,
            name2,
            prop2,
            fluid,
        ),
        "PCSAFT" => pcsaft_route(output, name1, prop1, name2, prop2, fluid),
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
    // Upstream's HelmholtzEOSBackend constructor checks the predefined
    // mixture library BEFORE the pure-fluid library (keys are "<Name>.mix"
    // and its uppercase form, exact match).
    #[cfg(feature = "heos-mixtures")]
    if let Some(pm) = predefined_mixture(fluid) {
        return heos_mixture_route(
            output,
            name1,
            prop1,
            name2,
            prop2,
            pm.fluids,
            pm.mole_fractions,
        );
    }
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
        InputPair::DmassT => (InputPair::DmolarT, v1 / mm, v2),
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
        /// A two-phase (D,T) state: upstream's `update_DmolarT` dome branch
        /// fills SatL/SatV with `update_TDmolarP_unchecked`, which skips
        /// `pre_update` — the sub-states' reducing cache stays zeroed and
        /// every caloric read through them throws (tau = inf). P/T/D/Q
        /// remain readable. Oracle-verified.
        TwoPhaseDt(CubicSatState),
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
        InputPair::DmolarT => match eos.dmolar_t_flash(v1, v2)? {
            rustprop_cubics::CubicDtState::Single { t, rho, p, phase } => {
                CubicState::Single { t, rho, p, phase }
            }
            rustprop_cubics::CubicDtState::TwoPhase(s) => CubicState::TwoPhaseDt(s),
        },
        // Upstream delegates every other pair to the HEOS flash routines,
        // which fail on the fabricated cubic fluid's unset ancillaries.
        _ => {
            return Err(Error::Value("type not set".into()));
        }
    };

    let (t, p, rho, q, phase) = match &state {
        CubicState::Single { t, rho, p, phase } => (*t, *p, *rho, -1.0, *phase),
        CubicState::TwoPhase(s) | CubicState::TwoPhaseDt(s) => {
            (s.t, s.p, s.rhomolar, s.q, Phase::Twophase)
        }
    };
    // Two-phase caloric outputs mix the saturated branches by quality; the
    // (D,T) dome states reproduce upstream's broken-sub-state throw instead.
    let mix = |f: &dyn Fn(f64, f64) -> f64| -> Result<f64> {
        match &state {
            CubicState::Single { t, rho, .. } => Ok(f(*t, *rho)),
            CubicState::TwoPhase(s) => Ok(q * f(s.t, s.rho_v) + (1.0 - q) * f(s.t, s.rho_l)),
            CubicState::TwoPhaseDt(_) => Err(Error::Value(
                "calc_alpha0_deriv_nocache returned invalid number with inputs nTau: 1, nDelta: 0, tau: inf, delta: 0"
                    .into(),
            )),
        }
    };
    Ok(match out {
        Param::T => t,
        Param::P => p,
        Param::Dmolar => rho,
        Param::Dmass => rho * mm,
        Param::Q => q,
        Param::Phase => f64::from(phase.index()),
        Param::Hmolar => mix(&|t, r| eos.hmolar(t, r))?,
        Param::Hmass => mix(&|t, r| eos.hmolar(t, r))? / mm,
        Param::Smolar => mix(&|t, r| eos.smolar(t, r))?,
        Param::Smass => mix(&|t, r| eos.smolar(t, r))? / mm,
        Param::Umolar => mix(&|t, r| eos.umolar(t, r))?,
        Param::Umass => mix(&|t, r| eos.umolar(t, r))? / mm,
        Param::Cpmolar => match &state {
            CubicState::Single { t, rho, .. } => eos.cpmolar(*t, *rho),
            CubicState::TwoPhase(_) | CubicState::TwoPhaseDt(_) => {
                return Err(Error::Value(
                    "Input is two-phase and cp is not defined".into(),
                ));
            }
        },
        Param::Cvmolar => match &state {
            CubicState::Single { t, rho, .. } => eos.cvmolar(*t, *rho),
            CubicState::TwoPhase(_) | CubicState::TwoPhaseDt(_) => {
                return Err(Error::Value(
                    "Input is two-phase and cv is not defined".into(),
                ));
            }
        },
        Param::SpeedSound => match &state {
            CubicState::Single { t, rho, .. } => eos.speed_sound(*t, *rho),
            CubicState::TwoPhase(_) | CubicState::TwoPhaseDt(_) => {
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

#[cfg(not(feature = "incompressible"))]
fn incomp_route(_: &str, _: &str, _: f64, _: &str, _: f64, _: &str) -> Result<f64> {
    Err(Error::NotImplemented(
        "the `incompressible` feature is not enabled".into(),
    ))
}

/// PropsSI route for `INCOMP::` (upstream `IncompressibleBackend`):
/// `Name`, `Name[x]` (bracket fraction), and `Name-40%` (percent) forms;
/// five input pairs; mass-based outputs only.
#[cfg(feature = "incompressible")]
fn incomp_route(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid: &str,
) -> Result<f64> {
    use rustprop_incompressible::IncompEos;

    // Upstream `extract_fractions`: bracket syntax first, then the
    // percent/dash form (which requires BOTH '-' and '%').
    let (name, mut x, had_fraction) = if fluid.contains('[') && fluid.contains(']') {
        let inner = fluid.strip_suffix(']').ok_or_else(|| {
            Error::Value(format!("Fluid entry [{fluid}] must end with ']' character"))
        })?;
        let (name, fraction) = inner
            .split_once('[')
            .ok_or_else(|| Error::Value(format!("Could not break [{inner}] into name/fraction")))?;
        let f: f64 = fraction
            .parse()
            .map_err(|_| Error::Value(format!("fraction [{fraction}] was not converted fully")))?;
        if !(0.0..=1.0).contains(&f) {
            return Err(Error::Value(format!(
                "fraction [{fraction}] was not converted to a value between 0 and 1 inclusive"
            )));
        }
        (name.to_string(), f, true)
    } else if fluid.contains('-') && fluid.contains('%') {
        let parts: Vec<&str> = fluid.split('-').collect();
        if parts.len() != 2 {
            return Err(Error::Value(format!(
                "Format of incompressible solution string [{fluid}] is invalid, should be like \"EG-20%\" or \"EG-0.2\" "
            )));
        }
        let frac = parts[1];
        let (num, pct) = match frac.strip_suffix('%') {
            Some(n) => (n, true),
            None => (frac, false),
        };
        let mut f: f64 = num
            .parse()
            .map_err(|_| Error::Value(format!("fraction [{frac}] was not converted fully")))?;
        if pct {
            f *= 0.01;
        }
        (parts[0].to_string(), f, true)
    } else {
        (fluid.to_string(), 1.0, false)
    };

    let data = rustprop_data::incompressible::INCOMP_FLUIDS
        .iter()
        .find(|f| f.name == name)
        .ok_or_else(|| {
            Error::Value(format!(
                "key [{name}] was not found in string_to_index_map in JSONIncompressibleLibrary"
            ))
        })?;
    let eos = IncompEos::new(data);

    // Upstream fraction setters: pure fluids force x = 1.0 silently; the
    // PropsSI default (no bracket) is also 1.0.
    use rustprop_core::fluid::IncompFrac;
    if data.xid == IncompFrac::Pure {
        x = 1.0;
    }
    let _ = had_fraction;

    let out = Param::parse(output).ok_or_else(|| {
        Error::Value(format!(
            "Output parameter parsing failed; error: Output string is invalid [{output}]"
        ))
    })?;

    // Trivial outputs (no state update).
    match out {
        Param::TMin => return Ok(data.tmin),
        Param::TMax => return Ok(data.tmax),
        Param::FractionMin => return Ok(data.xmin),
        Param::FractionMax => return Ok(data.xmax),
        Param::TFreeze => {
            eos.check_x(x)?;
            // Trivial output: upstream never runs update(), so _p is the
            // cleared sentinel (-HUGE) — inert for every poly/exppoly fit
            // (single row) and for the exponential forms (x-based).
            return eos.t_freeze(f64::NEG_INFINITY, x);
        }
        _ => {}
    }

    // Molar outputs are not defined for the mass-based INCOMP backend.
    let molar_err = |what: &str, use_instead: &str| -> Error {
        Error::NotImplemented(format!(
            "{what} is not defined for the INCOMP backend; use {use_instead} instead."
        ))
    };
    match out {
        Param::MolarMass => {
            return Err(Error::NotImplemented(
                "Molar mass is not defined for the INCOMP (incompressible) backend; INCOMP fluids are mass-based."
                    .into(),
            ));
        }
        Param::Dmolar => return Err(molar_err("Dmolar / rhomolar", "Dmass / rhomass")),
        Param::Hmolar => return Err(molar_err("Hmolar / hmolar", "Hmass / hmass")),
        Param::Smolar => return Err(molar_err("Smolar / smolar", "Smass / smass")),
        Param::Umolar => return Err(molar_err("Umolar / umolar", "Umass / umass")),
        Param::Cpmolar => return Err(molar_err("Cpmolar / cpmolar", "Cpmass / cpmass")),
        Param::Cvmolar => return Err(molar_err("Cvmolar / cvmolar", "Cvmass / cvmass")),
        _ => {}
    }

    let keys = (Param::parse(name1), Param::parse(name2));
    let pair = match keys {
        (Some(k1), Some(k2)) => generate_update_pair(k1, prop1, k2, prop2),
        _ => None,
    };
    let Some((pair, v1, v2)) = pair else {
        return Err(Error::Value(
            "Input pair variable is invalid and output(s) are non-trivial; cannot do state update"
                .into(),
        ));
    };

    // Composition gates (upstream update()).
    if data.xid == IncompFrac::Pure {
        // x forced to 1.0 above; always valid.
    } else if !(0.0..=1.0).contains(&x) {
        return Err(Error::Value(format!(
            "{} is a solution or brine. Mass fractions must be set to a vector with one entry between 0 and 1. [{x}] is not valid.",
            data.name
        )));
    }

    // The five supported pairs.
    let (t, p) = match pair {
        InputPair::PT => (v2, v1),
        InputPair::DmassP => (eos.t_from_rho(v1, x)?, v2),
        InputPair::PSmass => (eos.t_from_smass(v2, v1, x)?, v1),
        InputPair::HmassP => (eos.t_from_hmass(v1, v2, x)?, v2),
        InputPair::QT => {
            if v1 != 0.0 {
                return Err(Error::Value(
                    "Incompressible fluids can only handle saturated liquid, Q=0.".into(),
                ));
            }
            (v2, eos.psat(v2, x)?)
        }
        other => {
            return Err(Error::Value(format!(
                "This pair of inputs [{}] is not yet supported",
                other.short_desc()
            )));
        }
    };
    if p < 0.0 {
        return Err(Error::Value("p is less than zero".into()));
    }
    if !p.is_finite() {
        return Err(Error::Value("p is not a valid number".into()));
    }
    if t < 0.0 {
        return Err(Error::Value("T is less than zero".into()));
    }
    if !t.is_finite() {
        return Err(Error::Value("T is not a valid number".into()));
    }
    eos.check_tpx(t, p, x)?;

    Ok(match out {
        Param::T => t,
        Param::P => p,
        Param::Dmass => eos.rho(t, x)?,
        Param::Hmass => eos.hmass(t, p, x)?,
        Param::Smass => eos.smass(t, p, x)?,
        Param::Umass => eos.umass(t, p, x)?,
        Param::Cpmass | Param::Cvmass => eos.c(t, x)?,
        Param::Viscosity => eos.visc(t, x)?,
        Param::Conductivity => eos.cond(t, x)?,
        Param::Prandtl => eos.c(t, x)? * eos.visc(t, x)? / eos.cond(t, x)?,
        Param::Phase => f64::from(rustprop_core::params::Phase::Liquid.index()),
        other => {
            return Err(Error::NotImplemented(format!(
                "output parameter {} is not ported for the INCOMP backend yet",
                other.short_name()
            )));
        }
    })
}

// ---------------------------------------------------------------------------
// HEOS mixture route (upstream extract_fractions + HelmholtzEOSMixtureBackend)
// ---------------------------------------------------------------------------

/// Upstream `extract_fractions`, mole-fraction branch: split on `&`, each
/// entry `Name[frac]`; fractions outside [0,1] and unparseable fractions
/// throw verbatim; components at or below `10 * DBL_EPSILON` are silently
/// dropped (unless the string names a single fluid).
#[cfg(feature = "heos-mixtures")]
fn extract_mole_fractions(fluid_string: &str) -> Result<(Vec<&str>, Vec<f64>)> {
    let mut names = Vec::new();
    let mut fractions = Vec::new();
    let pairs: Vec<&str> = fluid_string.split('&').collect();
    for entry in &pairs {
        if !entry.ends_with(']') {
            return Err(Error::Value(format!(
                "Fluid entry [{entry}] must end with ']' character"
            )));
        }
        let body = &entry[..entry.len() - 1];
        let mut parts = body.split('[');
        let (name, fraction) = match (parts.next(), parts.next(), parts.next()) {
            (Some(n), Some(f), None) => (n, f),
            _ => {
                return Err(Error::Value(format!(
                    "Could not break [{body}] into name/fraction"
                )));
            }
        };
        let f: f64 = fraction
            .parse()
            .map_err(|_| Error::Value(format!("fraction [{fraction}] was not converted fully")))?;
        if !(0.0..=1.0).contains(&f) {
            return Err(Error::Value(format!(
                "fraction [{fraction}] was not converted to a value between 0 and 1 inclusive"
            )));
        }
        if f > 10.0 * f64::EPSILON || pairs.len() == 1 {
            fractions.push(f);
            names.push(name);
        }
    }
    Ok((names, fractions))
}

/// Entry point for HEOS fluid strings containing `&` or `[...]`: parse the
/// composition, collapse to the pure route when a single component remains.
#[cfg(feature = "heos-mixtures")]
fn heos_mixture_entry(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid: &str,
) -> Result<f64> {
    let (names, fractions) = if fluid.contains('[') && fluid.contains(']') {
        extract_mole_fractions(fluid)?
    } else {
        // `A&B` without fractions: upstream's default fractions vector [1.0]
        // fails the set_mole_fractions size check inside factory init.
        let names: Vec<&str> = fluid.split('&').collect();
        if names.len() > 1 {
            return Err(Error::Value(format!(
                "Initialize failed for backend: \"HEOS\", fluid: \"{fluid}\" fractions \"[ 1.0000000000 ]\"; error: size of mole fraction vector [1] does not equal that of component vector [{}]",
                names.len()
            )));
        }
        (names, vec![1.0])
    };
    if names.len() == 1 {
        return heos_route(output, name1, prop1, name2, prop2, names[0]);
    }
    heos_mixture_route(output, name1, prop1, name2, prop2, &names, &fractions)
}

#[cfg(feature = "heos-mixtures")]
use rustprop_heos::mixture::MixtureModel;

/// Without the `heos-mixtures` feature, mixture strings stay a loud error
/// (the pair/departure tables are deliberately not linked in).
#[cfg(not(feature = "heos-mixtures"))]
fn heos_mixture_entry(
    _output: &str,
    _name1: &str,
    _prop1: f64,
    _name2: &str,
    _prop2: f64,
    _fluid: &str,
) -> Result<f64> {
    Err(Error::NotImplemented(
        "the `heos-mixtures` feature is not enabled".into(),
    ))
}

/// Upstream `is_predefined_mixture`: the library registers each blend as
/// `"<name>.mix"` plus the all-uppercase form (emplace: first key wins).
#[cfg(feature = "heos-mixtures")]
fn predefined_mixture(name: &str) -> Option<&'static rustprop_core::fluid::PredefinedMixture> {
    use std::collections::HashMap;
    use std::sync::OnceLock;
    static MAP: OnceLock<HashMap<String, &'static rustprop_core::fluid::PredefinedMixture>> =
        OnceLock::new();
    let map = MAP.get_or_init(|| {
        let mut m: HashMap<String, &'static rustprop_core::fluid::PredefinedMixture> =
            HashMap::new();
        for pm in rustprop_data::mixtures::MIX_PREDEFINED {
            let key = format!("{}.mix", pm.name);
            m.entry(key.to_uppercase()).or_insert(pm);
            m.entry(key).or_insert(pm);
        }
        m
    });
    map.get(name).copied()
}

/// One cached `MixtureModel` per component set (composition-independent).
#[cfg(feature = "heos-mixtures")]
fn mixture_model(components: &[&'static FluidData]) -> Result<&'static MixtureModel> {
    use std::collections::HashMap;
    use std::sync::{Mutex, OnceLock};
    static MODELS: OnceLock<Mutex<HashMap<Vec<usize>, &'static MixtureModel>>> = OnceLock::new();
    let m = MODELS.get_or_init(|| Mutex::new(HashMap::new()));
    let key: Vec<usize> = components
        .iter()
        .map(|c| std::ptr::from_ref(*c) as usize)
        .collect();
    let mut guard = m.lock().expect("mixture model cache poisoned");
    if let Some(model) = guard.get(&key) {
        return Ok(model);
    }
    let model = MixtureModel::new(
        components,
        rustprop_data::mixtures::MIX_BINARY_PAIRS,
        rustprop_data::mixtures::MIX_DEPARTURE_FNS,
    )?;
    let leaked: &'static MixtureModel = Box::leak(Box::new(model));
    guard.insert(key, leaked);
    Ok(leaked)
}

/// The state a mixture flash publishes (single- or two-phase).
#[cfg(feature = "heos-mixtures")]
enum MixState {
    Single(rustprop_heos::mixture_flash::MixtureState),
    Two(rustprop_heos::mixture_vle::MixtureTwoPhase),
}

#[cfg(feature = "heos-mixtures")]
fn heos_mixture_route(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    names: &[&str],
    fractions: &[f64],
) -> Result<f64> {
    let components: Vec<&'static FluidData> = names
        .iter()
        .map(|n| resolve_fluid(n))
        .collect::<Result<_>>()?;
    let model = mixture_model(&components)?;
    let z = fractions.to_vec();

    let out = Param::parse(output).ok_or_else(|| {
        Error::Value(format!(
            "Output parameter parsing failed; error: Output string is invalid [{output}]"
        ))
    })?;
    let keys = (Param::parse(name1), Param::parse(name2));
    let pair = match keys {
        (Some(k1), Some(k2)) => generate_update_pair(k1, prop1, k2, prop2),
        _ => None,
    };
    let Some((pair, v1, v2)) = pair else {
        if out.is_trivial() {
            return mixture_trivial_output(model, &components, &z, out);
        }
        return Err(Error::Value(
            "Input pair variable is invalid and output(s) are non-trivial; cannot do state update"
                .into(),
        ));
    };
    if out.is_trivial() {
        return mixture_trivial_output(model, &components, &z, out);
    }
    let (p1, p2) = pair.split();
    if out == p1 {
        return Ok(v1);
    }
    if out == p2 {
        return Ok(v2);
    }

    let state = mixture_update(model, &z, pair, v1, v2)?;
    mixture_keyed_output(model, &components, &z, &state, out)
}

/// Upstream trivial outputs for mixtures: mole-fraction-weighted limits,
/// reducing-state values, R and M; the critical/acentric family fails with
/// PropsSI's outer wrapper message (the underlying exceptions are
/// mixture-unsupported upstream).
#[cfg(feature = "heos-mixtures")]
fn mixture_trivial_output(
    model: &MixtureModel,
    components: &[&'static FluidData],
    z: &[f64],
    out: Param,
) -> Result<f64> {
    let weighted = |f: &dyn Fn(&'static FluidData) -> f64| -> f64 {
        components.iter().zip(z).map(|(c, x)| f(c) * x).sum()
    };
    Ok(match out {
        Param::MolarMass => model.molar_mass(z),
        Param::GasConstant => model.gas_constant(),
        Param::TReducing => model.reducing.tr(z),
        Param::RhomolarReducing => model.reducing.rhormolar(z),
        Param::TMax => weighted(&|c| c.eos.t_max),
        Param::TMin => weighted(&|c| c.eos.sat_min_liquid.t),
        Param::TTriple => weighted(&|c| c.eos.sat_min_liquid.t),
        Param::PMax => weighted(&|c| c.eos.p_max),
        Param::PTriple => weighted(&|c| c.eos.sat_min_liquid.p),
        Param::TCritical
        | Param::PCritical
        | Param::RhomolarCritical
        | Param::RhomassCritical
        | Param::AcentricFactor => {
            return Err(Error::Value("No outputs were able to be calculated".into()));
        }
        other => {
            return Err(Error::NotImplemented(format!(
                "trivial output {} is not ported for mixtures",
                other.short_name()
            )));
        }
    })
}

/// Upstream `HelmholtzEOSMixtureBackend::update` for mixtures: PT/QT/PQ are
/// ported; DP/DQ throw upstream's NotImplemented verbatim; the sweep-based
/// pairs (DmolarT, HmolarP, ...) need the 10f stability machinery and defer
/// loudly (upstream computes them — documented deviation until 10f).
#[cfg(feature = "heos-mixtures")]
fn mixture_update(
    model: &'static MixtureModel,
    z: &[f64],
    pair: InputPair,
    v1: f64,
    v2: f64,
) -> Result<MixState> {
    let mm = model.molar_mass(z);
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
        InputPair::DmolarQmass => (InputPair::DmolarQ, v1, v2),
        InputPair::DmassQmass | InputPair::DmassQ => (InputPair::DmolarQ, v1 / mm, v2),
        InputPair::DmassHmass => (InputPair::DmolarHmolar, v1 / mm, v2 * mm),
        InputPair::DmassSmass => (InputPair::DmolarSmolar, v1 / mm, v2 * mm),
        InputPair::DmassUmass => (InputPair::DmolarUmolar, v1 / mm, v2 * mm),
        other => (other, v1, v2),
    };
    match pair {
        InputPair::PT => {
            // Full PT_flash_mixtures: stability test + Wilson cross-check +
            // Michelsen split; single-phase fallback inside.
            match rustprop_heos::mixture_stability::pt_flash_mixtures(model, z, v2, v1)? {
                rustprop_heos::mixture_stability::PtFlashResult::Single(s) => {
                    Ok(MixState::Single(s))
                }
                rustprop_heos::mixture_stability::PtFlashResult::TwoPhase {
                    t,
                    p,
                    q,
                    rhomolar,
                    x,
                    y,
                    rhomolar_liq,
                    rhomolar_vap,
                } => Ok(MixState::Two(rustprop_heos::mixture_vle::MixtureTwoPhase {
                    t,
                    p,
                    rhomolar,
                    q,
                    hmolar_liq: model.hmolar(&x, t, rhomolar_liq),
                    hmolar_vap: model.hmolar(&y, t, rhomolar_vap),
                    smolar_liq: model.smolar(&x, t, rhomolar_liq),
                    smolar_vap: model.smolar(&y, t, rhomolar_vap),
                    umolar_liq: model.umolar(&x, t, rhomolar_liq),
                    umolar_vap: model.umolar(&y, t, rhomolar_vap),
                    x_liq: x,
                    y_vap: y,
                    rhomolar_liq,
                    rhomolar_vap,
                })),
            }
        }
        InputPair::QT => {
            if !(0.0..=1.0).contains(&v1) {
                return Err(Error::OutOfRange(
                    "Input vapor quality [Q] must be between 0 and 1".into(),
                ));
            }
            Ok(MixState::Two(model.qt_flash(v1, v2, z)?))
        }
        InputPair::PQ => {
            if !(0.0..=1.0).contains(&v2) {
                return Err(Error::OutOfRange(
                    "Input vapor quality [Q] must be between 0 and 1".into(),
                ));
            }
            Ok(MixState::Two(model.pq_flash(v1, v2, z)?))
        }
        InputPair::DmolarP => Err(Error::NotImplemented(
            "DP_flash not ready for mixtures".into(),
        )),
        InputPair::DmolarQ => Err(Error::NotImplemented(
            "DQ_flash not ready for mixtures".into(),
        )),
        InputPair::HmassT | InputPair::TUmass | InputPair::SmolarUmolar => {
            Err(Error::Value(format!(
                "This pair of inputs [{}] is not yet supported",
                pair.short_desc()
            )))
        }
        InputPair::DmolarT => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::dhsu_t_flash_mixtures(model, z, v2, ms::SweepVar::Dmolar, v1)
        }),
        InputPair::HmolarT => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::dhsu_t_flash_mixtures(model, z, v2, ms::SweepVar::Hmolar, v1)
        }),
        InputPair::SmolarT => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::dhsu_t_flash_mixtures(model, z, v2, ms::SweepVar::Smolar, v1)
        }),
        InputPair::TUmolar => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::dhsu_t_flash_mixtures(model, z, v1, ms::SweepVar::Umolar, v2)
        }),
        InputPair::HmolarP => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::hsu_p_flash_mixtures(model, z, v2, ms::SweepVar::Hmolar, v1)
        }),
        InputPair::PSmolar => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::hsu_p_flash_mixtures(model, z, v1, ms::SweepVar::Smolar, v2)
        }),
        InputPair::PUmolar => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::hsu_p_flash_mixtures(model, z, v1, ms::SweepVar::Umolar, v2)
        }),
        InputPair::DmolarHmolar => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::hsu_d_flash_mixtures(model, z, v1, ms::SweepVar::Hmolar, v2)
        }),
        InputPair::DmolarSmolar => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::hsu_d_flash_mixtures(model, z, v1, ms::SweepVar::Smolar, v2)
        }),
        InputPair::DmolarUmolar => mix_sweep_state(model, {
            use rustprop_heos::mixture_sweep as ms;
            ms::hsu_d_flash_mixtures(model, z, v1, ms::SweepVar::Umolar, v2)
        }),
        other => Err(Error::NotImplemented(format!(
            "mixture input pair {} is not ported yet",
            other.short_desc()
        ))),
    }
}

/// Converts a sweep-flash publish into the props MixState (two-phase results
/// carry per-phase compositions for the lever-rule caloric outputs).
#[cfg(feature = "heos-mixtures")]
fn mix_sweep_state(
    model: &'static MixtureModel,
    result: rustprop_core::Result<rustprop_heos::mixture_stability::PtFlashResult>,
) -> Result<MixState> {
    match result? {
        rustprop_heos::mixture_stability::PtFlashResult::Single(s) => Ok(MixState::Single(s)),
        rustprop_heos::mixture_stability::PtFlashResult::TwoPhase {
            t,
            p,
            q,
            rhomolar,
            x,
            y,
            rhomolar_liq,
            rhomolar_vap,
        } => Ok(MixState::Two(rustprop_heos::mixture_vle::MixtureTwoPhase {
            t,
            p,
            rhomolar,
            q,
            hmolar_liq: model.hmolar(&x, t, rhomolar_liq),
            hmolar_vap: model.hmolar(&y, t, rhomolar_vap),
            smolar_liq: model.smolar(&x, t, rhomolar_liq),
            smolar_vap: model.smolar(&y, t, rhomolar_vap),
            umolar_liq: model.umolar(&x, t, rhomolar_liq),
            umolar_vap: model.umolar(&y, t, rhomolar_vap),
            x_liq: x,
            y_vap: y,
            rhomolar_liq,
            rhomolar_vap,
        })),
    }
}

#[cfg(feature = "heos-mixtures")]
fn mixture_keyed_output(
    model: &MixtureModel,
    components: &[&'static FluidData],
    z: &[f64],
    state: &MixState,
    out: Param,
) -> Result<f64> {
    let mm = model.molar_mass(z);
    let (t, p, rhomolar, q) = match state {
        MixState::Single(s) => (s.t, s.p, s.rhomolar, s.q),
        MixState::Two(s) => (s.t, s.p, s.rhomolar, s.q),
    };
    let single_phase_rho = |what: &str| -> Result<f64> {
        match state {
            MixState::Single(s) => Ok(s.rhomolar),
            MixState::Two(_) => Err(Error::Value(format!(
                "Input is two-phase and {what} is not defined"
            ))),
        }
    };
    Ok(match out {
        Param::T => t,
        Param::P => p,
        Param::Q => q,
        Param::Dmolar => rhomolar,
        Param::Dmass => rhomolar * mm,
        Param::Hmolar | Param::Hmass => {
            let h = match state {
                MixState::Single(s) => model.hmolar(z, s.t, s.rhomolar),
                MixState::Two(s) => s.hmolar(),
            };
            if out == Param::Hmass { h / mm } else { h }
        }
        Param::Smolar | Param::Smass => {
            let sv = match state {
                MixState::Single(s) => model.smolar(z, s.t, s.rhomolar),
                MixState::Two(s) => s.smolar(),
            };
            if out == Param::Smass { sv / mm } else { sv }
        }
        Param::Umolar | Param::Umass => {
            let u = match state {
                MixState::Single(s) => model.umolar(z, s.t, s.rhomolar),
                MixState::Two(s) => s.umolar(),
            };
            if out == Param::Umass { u / mm } else { u }
        }
        Param::Cpmolar => model.cpmolar(z, t, single_phase_rho("cpmolar")?),
        Param::Cpmass => model.cpmolar(z, t, single_phase_rho("cpmass")?) / mm,
        Param::Cvmolar => model.cvmolar(z, t, single_phase_rho("cvmolar")?),
        Param::Cvmass => model.cvmolar(z, t, single_phase_rho("cvmass")?) / mm,
        Param::SpeedSound => model.speed_sound(z, t, single_phase_rho("speed_sound")?),
        Param::Gmolar => model.gibbsmolar_nocache(z, t, single_phase_rho("gibbsmolar")?),
        Param::Gmass => model.gibbsmolar_nocache(z, t, single_phase_rho("gibbsmass")?) / mm,
        Param::SurfaceTension => {
            return Err(Error::NotImplemented(
                "surface tension not implemented for mixtures".into(),
            ));
        }
        Param::Viscosity => {
            // Upstream's "highly approximate" mixture model: log-linear
            // mixing of PURE-component viscosities, each evaluated as a
            // pure fluid at the BULK (rhomolar, T).
            let mut summer = 0.0_f64;
            for (i, comp) in components.iter().enumerate() {
                let flash = fluid_flash(comp);
                let p_pure = flash.eos.pressure(t, rhomolar);
                let v = viscosity_model(comp)?;
                let eta = rustprop_heos::transport::viscosity(
                    &flash.eos,
                    comp,
                    v,
                    t,
                    rhomolar,
                    p_pure,
                    Some(&ecs_resolver),
                )?;
                summer += z[i] * eta.ln();
            }
            summer.exp()
        }
        Param::Conductivity => {
            // Linear mixing of pure-component conductivities at bulk state.
            let mut summer = 0.0;
            for (i, comp) in components.iter().enumerate() {
                let flash = fluid_flash(comp);
                let p_pure = flash.eos.pressure(t, rhomolar);
                let tr = comp.transport.as_ref().ok_or_else(|| {
                    Error::Value(
                        "Thermal conductivity model is not available for this fluid".into(),
                    )
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
                let v = viscosity_model(comp).ok();
                let lambda = rustprop_heos::transport::conductivity(
                    &flash.eos,
                    comp,
                    c,
                    v,
                    t,
                    rhomolar,
                    p_pure,
                    Some(&ecs_resolver),
                )?;
                summer += z[i] * lambda;
            }
            summer
        }
        other => {
            return Err(Error::NotImplemented(format!(
                "output parameter {} is not ported for mixtures",
                other.short_name()
            )));
        }
    })
}

// ---------------------------------------------------------------------------
// PC-SAFT route (upstream PCSAFTBackend + PropsSI surface, survey §4)
// ---------------------------------------------------------------------------

#[cfg(not(feature = "pcsaft"))]
fn pcsaft_route(_: &str, _: &str, _: f64, _: &str, _: f64, _: &str) -> Result<f64> {
    Err(Error::NotImplemented(
        "the `pcsaft` feature is not enabled".into(),
    ))
}

/// Upstream `PCSAFTLibrary` key registration: CAS, name, each alias, and
/// upper(alias) (the base name's uppercase is NOT registered).
#[cfg(feature = "pcsaft")]
fn resolve_pcsaft_fluid(key: &str) -> Result<&'static rustprop_core::fluid::PcsaftFluid> {
    use std::collections::HashMap;
    use std::sync::OnceLock;
    static MAP: OnceLock<HashMap<String, &'static rustprop_core::fluid::PcsaftFluid>> =
        OnceLock::new();
    let map = MAP.get_or_init(|| {
        let mut m: HashMap<String, &'static rustprop_core::fluid::PcsaftFluid> = HashMap::new();
        for f in rustprop_data::pcsaft::PCSAFT_FLUIDS {
            m.insert(f.cas.to_string(), f);
            m.insert(f.name.to_string(), f);
            for alias in f.aliases {
                m.insert((*alias).to_string(), f);
                m.insert(alias.to_uppercase(), f);
            }
        }
        m
    });
    map.get(key).copied().ok_or_else(|| {
        Error::Value(format!(
            "key [{key}] was not found in string_to_index_map in PCSAFTLibraryClass"
        ))
    })
}

#[cfg(feature = "pcsaft")]
fn pcsaft_route(
    output: &str,
    name1: &str,
    prop1: f64,
    name2: &str,
    prop2: f64,
    fluid: &str,
) -> Result<f64> {
    use rustprop_pcsaft::{PcsaftBackend, PcsaftInput};

    // Composition parsing shares upstream extract_fractions semantics.
    let (names, fractions) = if fluid.contains('[') && fluid.contains(']') {
        extract_mole_fractions_pcsaft(fluid)?
    } else {
        (fluid.split('&').collect::<Vec<_>>(), vec![1.0])
    };
    let fluids: Vec<_> = names
        .iter()
        .map(|n| resolve_pcsaft_fluid(n))
        .collect::<Result<_>>()?;
    let mut backend = PcsaftBackend::new(&fluids, rustprop_data::pcsaft::PCSAFT_BINARY_PAIRS)?;
    if names.len() > 1 {
        backend.set_mole_fractions(&fractions);
    }

    let out = Param::parse(output).ok_or_else(|| {
        Error::Value(format!(
            "Output parameter parsing failed; error: Output string is invalid [{output}]"
        ))
    })?;

    // Trivial outputs available without a state
    if out == Param::MolarMass {
        return Ok(backend.molar_mass());
    }

    let keys = (Param::parse(name1), Param::parse(name2));
    let pair = match keys {
        (Some(k1), Some(k2)) => generate_update_pair(k1, prop1, k2, prop2),
        _ => None,
    };
    let Some((pair, v1, v2)) = pair else {
        return Err(Error::Value(
            "Input pair variable is invalid and output(s) are non-trivial; cannot do state update"
                .into(),
        ));
    };
    let (p1, p2) = pair.split();
    if out == p1 {
        return Ok(v1);
    }
    if out == p2 {
        return Ok(v2);
    }

    // Mass->molar conversions with the composition's molar mass.
    let mm = backend.molar_mass();
    let (pcsaft_pair, u1, u2) = match pair {
        InputPair::PT => (PcsaftInput::Pt, v1, v2),
        InputPair::QT | InputPair::QmassT => (PcsaftInput::Qt, v1, v2),
        InputPair::PQ | InputPair::PQmass => (PcsaftInput::Pq, v1, v2),
        InputPair::DmolarT => (PcsaftInput::DmolarT, v1, v2),
        InputPair::DmassT => (PcsaftInput::DmolarT, v1 / mm, v2),
        other => {
            return Err(Error::Value(format!(
                "This pair of inputs [{}] is not yet supported",
                other.short_desc()
            )));
        }
    };
    backend.update(pcsaft_pair, u1, u2)?;

    Ok(match out {
        Param::T => backend.t,
        Param::P => backend.p,
        Param::Q => backend.q,
        Param::Phase => match backend.phase {
            rustprop_pcsaft::PcsaftPhase::Liquid => 0.0,
            rustprop_pcsaft::PcsaftPhase::Supercritical => 1.0,
            rustprop_pcsaft::PcsaftPhase::SupercriticalGas => 2.0,
            rustprop_pcsaft::PcsaftPhase::SupercriticalLiquid => 3.0,
            rustprop_pcsaft::PcsaftPhase::Gas => 5.0,
            rustprop_pcsaft::PcsaftPhase::TwoPhase => 6.0,
            rustprop_pcsaft::PcsaftPhase::Unknown => 7.0,
        },
        Param::Dmolar => backend.rhomolar,
        Param::Dmass => backend.rhomolar * mm,
        Param::Alphar => backend.calc_alphar(),
        Param::HmolarResidual => backend.calc_hmolar_residual(),
        Param::SmolarResidual => backend.calc_smolar_residual(),
        Param::GmolarResidual => backend.calc_gibbsmolar_residual(),
        // Every absolute caloric/transport output is a base-class
        // NotImplementedError upstream (PCSAFT overrides none of them).
        other => {
            return Err(Error::NotImplemented(format!(
                "Output [{}] is not implemented for this backend",
                other.short_name()
            )));
        }
    })
}

/// extract_fractions' mole-fraction branch for the PCSAFT route (same
/// semantics as the HEOS one; duplicated to stay independent of the
/// heos-mixtures feature gate).
#[cfg(feature = "pcsaft")]
fn extract_mole_fractions_pcsaft(fluid_string: &str) -> Result<(Vec<&str>, Vec<f64>)> {
    let mut names = Vec::new();
    let mut fractions = Vec::new();
    let pairs: Vec<&str> = fluid_string.split('&').collect();
    for entry in &pairs {
        if !entry.ends_with(']') {
            return Err(Error::Value(format!(
                "Fluid entry [{entry}] must end with ']' character"
            )));
        }
        let body = &entry[..entry.len() - 1];
        let mut parts = body.split('[');
        let (name, fraction) = match (parts.next(), parts.next(), parts.next()) {
            (Some(n), Some(f), None) => (n, f),
            _ => {
                return Err(Error::Value(format!(
                    "Could not break [{body}] into name/fraction"
                )));
            }
        };
        let f: f64 = fraction
            .parse()
            .map_err(|_| Error::Value(format!("fraction [{fraction}] was not converted fully")))?;
        if !(0.0..=1.0).contains(&f) {
            return Err(Error::Value(format!(
                "fraction [{fraction}] was not converted to a value between 0 and 1 inclusive"
            )));
        }
        if f > 10.0 * f64::EPSILON || pairs.len() == 1 {
            fractions.push(f);
            names.push(name);
        }
    }
    Ok((names, fractions))
}

// ---------------------------------------------------------------------------
// Tabular route: upstream REJECTS the tabular backends here
// ---------------------------------------------------------------------------

/// Which interpolation scheme the backend string selected.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum TabularScheme {
    Ttse,
    Bicubic,
}

/// Upstream `TabularBackend::available_in_high_level()` returns FALSE — "None
/// of the tabular methods are available from the high-level interface"
/// (TabularBackends.h:1077) — so `_PropsSImulti` rejects any TTSE/BICUBIC
/// backend string before it can update a state. The tables are a LOW-LEVEL
/// API only; use `rustprop_tabular::TabularState` for them, which is what
/// upstream's `AbstractState::factory("TTSE&HEOS", ...)` gives you.
fn tabular_route(
    _scheme: TabularScheme,
    _output: &str,
    _name1: &str,
    _prop1: f64,
    _name2: &str,
    _prop2: f64,
    _fluid: &str,
) -> Result<f64> {
    Err(Error::Value(
        "This AbstractState derived class cannot be used in the high-level interface; see www.coolprop.org/dev/coolprop/LowLevelAPI.html".into(),
    ))
}
