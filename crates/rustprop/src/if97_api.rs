//! Typed `PropsSI`-style dispatch for the IF97 engine (PLAN.md 2.4).
//!
//! Resolves two (parameter, value) inputs through
//! [`rustprop_core::generate_update_pair`] and routes to the same engine
//! functions CoolProp's IF97 backend calls, mirroring
//! `src/Backends/IF97/IF97Backend.h` structurally: `update()` resolves the
//! cached `(_T, _p, _Q, _phase)` state once, then every output is served
//! through the `calc_Flash` shape (saturated-branch shortcut within 1e-10 of
//! Q = 0/1, lever rule / loud `NotImplementedError` inside the dome, forward
//! (T, p) evaluators outside it), with the `_reverse` short-circuits for
//! s(p,h) / h(p,s) and AbstractState's mass<->molar conversions on top.

use rustprop_core::params::Phase;
use rustprop_core::{Error, InputPair, Param, Result, generate_update_pair};
use rustprop_if97 as if97;

/// `PropsSI`-style single-output call for IF97 water.
///
/// `props(output, name1, value1, name2, value2)` mirrors
/// `PropsSI(output, name1, value1, name2, value2, "IF97::Water")`.
pub fn props(output: Param, name1: Param, value1: f64, name2: Param, value2: f64) -> Result<f64> {
    // Trivial outputs need no state.
    match output {
        Param::TCritical => return Ok(if97::get_tcrit()),
        Param::PCritical => return Ok(if97::get_pcrit()),
        Param::RhomassCritical => return Ok(if97::get_rhocrit()),
        // calc_rhomolar_critical (IF97Backend.h): rhomass_critical / M.
        Param::RhomolarCritical => return Ok(if97::get_rhocrit() / if97::get_mw()),
        Param::TTriple => return Ok(if97::get_ttrip()),
        Param::PTriple => return Ok(if97::get_ptrip()),
        Param::TMin => return Ok(if97::get_tmin()),
        Param::TMax => return Ok(if97::get_tmax()),
        // Golden-verified: the CoolProp IF97 backend reports the triple
        // pressure (611.657 Pa) as PMIN, not IF97's Pmin (611.213 Pa).
        Param::PMin => return Ok(if97::get_ptrip()),
        Param::PMax => return Ok(if97::get_pmax()),
        Param::MolarMass => return Ok(if97::get_mw()),
        // calc_gas_constant (IF97Backend.h): mass-based Rgas put on the molar
        // basis per CoolProp convention (wheel: 8.314514578968002).
        Param::GasConstant => return Ok(if97::get_rgas() * if97::get_mw()),
        Param::AcentricFactor => return if97::get_acentric(),
        _ => {}
    }

    let (pair, v1, v2) = generate_update_pair(name1, value1, name2, value2)
        .ok_or_else(|| Error::Value("This pair of inputs is not yet supported".into()))?;

    // Upstream's generic outputs-are-inputs echo (src/CoolProp.cpp
    // `_PropsSI_outputs`, the `all_outputs_in_inputs` route at ~367-389 and
    // 440-453): once the input pair parses, an output whose parameter equals
    // either input's parameter returns that input's RAW value with NO state
    // update — even for pairs this backend cannot serve and for states no
    // flash could reach (wheel-confirmed: Dmass=5 from (Dmass, T), Q=0.5 from
    // (Q=0.5, T=1e6 K), Q=5 from (P, Q=5)). Trivial outputs never get here —
    // upstream a trivial output disables the echo — matching the early
    // returns above.
    let (p1, p2) = pair.split();
    if output == p1 {
        return Ok(v1);
    }
    if output == p2 {
        return Ok(v2);
    }

    match pair {
        // canonical order: P, T
        InputPair::PT => serve(output, &update_pt(v1, v2)?),
        // canonical order: P, Q
        InputPair::PQ => serve(output, &update_pq(v1, v2)?),
        // canonical order: Q, T
        InputPair::QT => serve(output, &update_qt(v1, v2)?),
        // canonical order: Hmass, P
        InputPair::HmassP => serve(output, &update_hp(v2, v1)?),
        // HmolarP_INPUTS: molar h converted to the mass basis at update, then
        // identical to HmassP (IF97Backend.h:213-217).
        InputPair::HmolarP => serve(output, &update_hp(v2, v1 / if97::MW)?),
        // canonical order: P, Smass
        InputPair::PSmass => serve(output, &update_ps(v1, v2)?),
        InputPair::PSmolar => serve(output, &update_ps(v1, v2 / if97::MW)?),
        // canonical order: Hmass, Smass
        InputPair::HmassSmass => serve(output, &update_hs(v1, v2)?),
        // Upstream DEFECT reproduced (IF97Backend.h:257-264): the
        // HmolarSmolar case converts to the mass basis, but the HmassSmass
        // fall-through then unconditionally reassigns `_hmass = value1;
        // _smass = value2;` — the RAW MOLAR values are used as mass values.
        // Wheel-confirmed: T(Hmolar=5e4, Smolar=120) = 281.154 K, i.e.
        // T(h = 5e4 J/kg, s = 120 J/kg/K).
        InputPair::HmolarSmolar => serve(output, &update_hs(v1, v2)?),
        // Upstream update()'s default arm, message verbatim (every pair the
        // backend serves is ported above, so anything else is this ValueError
        // in the wheel too).
        _ => Err(Error::Value(
            "This pair of inputs is not yet supported".into(),
        )),
    }
}

// ---------------------------------------------------------------------------
// update(): resolve the cached backend state once per call
// ---------------------------------------------------------------------------

/// Cached state after upstream `IF97Backend::update()`: `_T`, `_p`, `_Q`,
/// `_phase`, plus the `_reverse` h/s cache that redirects
/// `calc_smass`/`calc_hmass` to the dedicated reverse evaluators.
struct If97State {
    t: f64,
    p: f64,
    /// `_Q`: clamped lever quality inside the dome, -1 sentinel outside.
    q: f64,
    /// `_phase`, as `keyed_output(iPhase)` reports it.
    phase: Phase,
    /// `_reverse` with `_hmass` cached ((H,P) inputs): `calc_smass`
    /// short-circuits to `smass_phmass(p, h)`, in the dome too.
    reverse_h: Option<f64>,
    /// `_reverse` with `_smass` cached ((P,S) inputs): `calc_hmass`
    /// short-circuits to `hmass_psmass(p, s)`.
    reverse_s: Option<f64>,
}

impl If97State {
    fn two_phase(&self) -> bool {
        self.phase == Phase::Twophase
    }
}

/// Upstream `update(PT_INPUTS)` (IF97Backend.h:174-198 shape): cache `_p`,
/// `_T`, `_Q = -1`, and classify through `set_phase()`'s forward branch —
/// EAGERLY, so its `psat97(T)` throw poisons every output (wheel-confirmed:
/// at T = 200 K even Q, Phase, Z and Gmass error "Temperature out of range").
/// The tag source then throws a ValueError for an in-band twophase result;
/// the shipped wheel never does (see [`forward_phase`]), so no throw here.
fn update_pt(p: f64, t: f64) -> Result<If97State> {
    Ok(If97State {
        t,
        p,
        q: -1.0,
        phase: forward_phase(t, p)?,
        reverse_h: None,
        reverse_s: None,
    })
}

/// Upstream Q-range guard shared by PQ/QT (`OutOfRangeError`, verbatim).
// Ported comparison keeps upstream's literal `(_Q < 0) || (_Q > 1)` form.
#[allow(clippy::manual_range_contains)]
fn check_q_range(q: f64) -> Result<()> {
    if q < 0.0 || q > 1.0 {
        return Err(Error::OutOfRange(
            "Input vapor quality [Q] must be between 0 and 1".into(),
        ));
    }
    Ok(())
}

/// Upstream `update(PQ_INPUTS)` (IF97Backend.h:199-205): Q-range guard, then
/// `_T = Tsat97(p)` (throws off the saturation curve), phase hard twophase.
fn update_pq(p: f64, q: f64) -> Result<If97State> {
    check_q_range(q)?;
    Ok(If97State {
        t: if97::tsat97(p)?,
        p,
        q,
        phase: Phase::Twophase,
        reverse_h: None,
        reverse_s: None,
    })
}

/// Upstream `update(QT_INPUTS)` (IF97Backend.h:206-212). `_T` stays the RAW
/// input temperature — dome surface tension evaluates `sigma97` at it, not at
/// the `Tsat97(psat97(T))` round trip (wheel: 0.058916822384634644 at
/// T = 373.124295847684, vs 0.0589168215843171 through the PQ route).
fn update_qt(q: f64, t: f64) -> Result<If97State> {
    check_q_range(q)?;
    Ok(If97State {
        t,
        p: if97::psat97(t)?,
        q,
        phase: Phase::Twophase,
        reverse_h: None,
        reverse_s: None,
    })
}

/// Upstream `update(HmassP_INPUTS)` (IF97Backend.h:218-234): backward `_T`
/// first, then the BackwardRegion dome check; single phase classifies through
/// `set_phase()`'s reverse branch.
fn update_hp(p: f64, h: f64) -> Result<If97State> {
    let t = if97::t_phmass(p, h)?;
    let region = if97::region_ph(p, h)?;
    let (q, phase) = if region == 4 {
        // `_Q = min(1, max(0, (h - hL)/(hV - hL)))` — q_phmass builds the
        // identical clamped lever from the same saturation-curve evaluators.
        (if97::q_phmass(p, h)?, Phase::Twophase)
    } else {
        (-1.0, reverse_phase(region, t, p)?)
    };
    Ok(If97State {
        t,
        p,
        q,
        phase,
        reverse_h: Some(h),
        reverse_s: None,
    })
}

/// Upstream `update(PSmass_INPUTS)` (IF97Backend.h:240-256), mirror of the
/// (H,P) case on the entropy backward.
fn update_ps(p: f64, s: f64) -> Result<If97State> {
    let t = if97::t_psmass(p, s)?;
    let region = if97::region_ps(p, s)?;
    let (q, phase) = if region == 4 {
        (if97::q_psmass(p, s)?, Phase::Twophase)
    } else {
        (-1.0, reverse_phase(region, t, p)?)
    };
    Ok(If97State {
        t,
        p,
        q,
        phase,
        reverse_h: None,
        reverse_s: Some(s),
    })
}

/// Upstream `update(HmassSmass_INPUTS)` (IF97Backend.h:262-278): `_p` from the
/// (h,s) backward, `_T` from the (p,h) backward, dome check on (p,h).
/// `_reverse` is NOT set on this pair, so Hmolar/Smolar round-trip through the
/// forward evaluators — only the raw-echoed Hmass/Smass return the inputs —
/// and single phase classifies through `set_phase()`'s FORWARD branch.
fn update_hs(h: f64, s: f64) -> Result<If97State> {
    let p = if97::p_hsmass(h, s)?;
    let t = if97::t_phmass(p, h)?;
    let region = if97::region_ph(p, h)?;
    let (q, phase) = if region == 4 {
        (if97::q_phmass(p, h)?, Phase::Twophase)
    } else {
        (-1.0, forward_phase(t, p)?)
    };
    Ok(If97State {
        t,
        p,
        q,
        phase,
        reverse_h: None,
        reverse_s: None,
    })
}

/// `set_phase()`'s reverse branch (IF97Backend.h:110-155): CoolProp phase from
/// the IF97 backward region. Region 3 consults `Tsat97(p)`, which throws
/// "Pressure out of range" for p > pcrit — upstream defect reproduced
/// (wheel-confirmed: EVERY output of (H,P)/(P,S) errors so at, e.g.,
/// Hmass = 2.2e6, P = 30 MPa).
fn reverse_phase(region: i32, t: f64, p: f64) -> Result<Phase> {
    Ok(match region {
        1 => {
            if p <= if97::PCRIT {
                Phase::Liquid
            } else {
                Phase::SupercriticalLiquid
            }
        }
        2 => {
            if t <= if97::TCRIT {
                Phase::Gas
            } else {
                Phase::SupercriticalGas
            }
        }
        3 => {
            if t < if97::tsat97(p)? {
                if p <= if97::PCRIT {
                    Phase::Liquid
                } else {
                    Phase::SupercriticalLiquid
                }
            } else if t <= if97::TCRIT {
                Phase::Gas
            } else {
                Phase::SupercriticalGas
            }
        }
        // Region 5 / out of bounds (upstream `case 5: default:`) — in
        // practice dead, the backward solvers throw range errors first.
        _ => {
            return Err(Error::OutOfRange(
                "Outside of IF97 Reverse Function Bounds".into(),
            ));
        }
    })
}

/// `set_phase()`'s forward branch as the 8.0.0 WHEEL ships it (used by the
/// (P,T) and single-phase (H,S) updates). The tag source classifies the whole
/// ±3.3e-5 saturation band as twophase and makes the PT update throw there;
/// the wheel — gitrevision ae81610, same as the tag file — demonstrably does
/// neither. Probed 2026-08 against the wheel's Phase output on (T,P):
/// - critical-point band intact: |T-Tc| < eps/10 && |p-Pc| < eps -> 4
///   (Phase = 4 at dT = 3e-6 K, 2/3 just outside either band);
/// - supercritical splits as in the tag (1/2/3, boundaries non-inclusive:
///   Phase = 2 at T = Tc+1e-5, p = Pc exactly);
/// - subcritical HALF-band: gas only for p < psat97(T)*(1 - eps); on-curve
///   and the whole band (ladder: gas first at p = psat*(1 - 3.4e-5)) come
///   out LIQUID, never twophase, and nothing throws.
fn forward_phase(t: f64, p: f64) -> Result<Phase> {
    // IAPWS-IF97 RMS saturated-pressure inconsistency (upstream literal).
    const EPSILON: f64 = 3.3e-5;
    Ok(
        if (t - if97::TCRIT).abs() < EPSILON / 10.0 && (p - if97::PCRIT).abs() < EPSILON {
            Phase::CriticalPoint
        } else if t > if97::TCRIT {
            if p > if97::PCRIT {
                Phase::Supercritical
            } else {
                Phase::SupercriticalGas
            }
        } else if p > if97::PCRIT {
            Phase::SupercriticalLiquid
        } else if p < if97::psat97(t)? * (1.0 - EPSILON) {
            Phase::Gas
        } else {
            Phase::Liquid
        },
    )
}

// ---------------------------------------------------------------------------
// keyed_output(): serve one output from the resolved state
// ---------------------------------------------------------------------------

/// `calc_Flash` property key (IF97Backend.h:370-448).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FlashKey {
    Dmass,
    Hmass,
    Smass,
    Umass,
    Cpmass,
    Cvmass,
    SpeedSound,
    Viscosity,
    Conductivity,
    SurfaceTension,
    Prandtl,
}

/// Upstream `keyed_output` -> the backend `calc_*` overrides, plus
/// AbstractState's mass<->molar wrappers (specific molar = mass * M,
/// Dmolar = Dmass / M). Outputs the wheel also refuses (Z, Gmass, PIP,
/// Cp0mass, ...) keep erroring through the fallback arm — upstream raises
/// the base-class `calc_* is not implemented for this backend` there.
fn serve(output: Param, st: &If97State) -> Result<f64> {
    match output {
        Param::T => Ok(st.t),
        Param::P => Ok(st.p),
        Param::Q => Ok(st.q),
        Param::Qmass => qmass(st),
        Param::Phase => Ok(f64::from(st.phase.index())),
        Param::Dmass => flash(FlashKey::Dmass, st),
        Param::Dmolar => Ok(flash(FlashKey::Dmass, st)? / if97::MW),
        Param::Hmass => calc_hmass(st),
        Param::Hmolar => Ok(calc_hmass(st)? * if97::MW),
        Param::Smass => calc_smass(st),
        Param::Smolar => Ok(calc_smass(st)? * if97::MW),
        Param::Umass => flash(FlashKey::Umass, st),
        Param::Umolar => Ok(flash(FlashKey::Umass, st)? * if97::MW),
        Param::Cpmass => flash(FlashKey::Cpmass, st),
        Param::Cpmolar => Ok(flash(FlashKey::Cpmass, st)? * if97::MW),
        Param::Cvmass => flash(FlashKey::Cvmass, st),
        Param::Cvmolar => Ok(flash(FlashKey::Cvmass, st)? * if97::MW),
        Param::SpeedSound => flash(FlashKey::SpeedSound, st),
        Param::Viscosity => flash(FlashKey::Viscosity, st),
        Param::Conductivity => flash(FlashKey::Conductivity, st),
        Param::SurfaceTension => flash(FlashKey::SurfaceTension, st),
        Param::Prandtl => flash(FlashKey::Prandtl, st),
        _ => Err(Error::NotImplemented(format!(
            "Output {output:?} not implemented for IF97"
        ))),
    }
}

/// Upstream `AbstractState::calc_Qmass` (src/AbstractState.cpp:822-841):
/// out-of-[0,1] `_Q` (the -1 single-phase sentinel included) raises a
/// ValueError, EXACT endpoints return `_Q` (`== 0.0 || == 1.0` — NOT
/// calc_Flash's 1e-10 band: wheel errors at Q = 1e-11), and any other
/// quality reaches `calc_phase_molar_masses` -> `get_mole_fractions()`,
/// which this backend overrides to throw (messages verbatim).
fn qmass(st: &If97State) -> Result<f64> {
    if st.q < 0.0 || st.q > 1.0 {
        return Err(Error::Value(
            "Qmass requires a two-phase state (0 <= Q <= 1)".into(),
        ));
    }
    if st.q == 0.0 || st.q == 1.0 {
        return Ok(st.q);
    }
    Err(Error::NotImplemented(
        "get_mole_fractions composition has not been implemented.".into(),
    ))
}

/// `calc_hmass`: the `_reverse && _smass` short-circuit ((P,S) inputs) uses
/// the dedicated reverse h(p,s); everything else flows through `calc_Flash`.
fn calc_hmass(st: &If97State) -> Result<f64> {
    match st.reverse_s {
        Some(s) => if97::hmass_psmass(st.p, s),
        None => flash(FlashKey::Hmass, st),
    }
}

/// `calc_smass`: the `_reverse && _hmass` short-circuit ((H,P) inputs) uses
/// the dedicated reverse s(p,h).
fn calc_smass(st: &If97State) -> Result<f64> {
    match st.reverse_h {
        Some(h) => if97::smass_phmass(st.p, h),
        None => flash(FlashKey::Smass, st),
    }
}

/// Upstream `calc_Flash` (IF97Backend.h:370-448).
fn flash(key: FlashKey, st: &If97State) -> Result<f64> {
    if st.two_phase() {
        if st.q.abs() < 1e-10 {
            // bubble point (Q == 0) on the saturated-liquid curve
            sat_liquid(key, st)
        } else if (st.q - 1.0).abs() < 1e-10 {
            // dew point (Q == 1) on the saturated-vapor curve
            sat_vapor(key, st)
        } else {
            match key {
                // Density is inverse phase-weighted (reciprocal of v).
                FlashKey::Dmass => Ok(1.0
                    / (st.q / sat_vapor(FlashKey::Dmass, st)?
                        + (1.0 - st.q) / sat_liquid(FlashKey::Dmass, st)?)),
                // Messages verbatim from upstream `NotImplementedError`s.
                FlashKey::Cpmass => Err(Error::NotImplemented(
                    "Isobaric Specific Heat not valid in two phase region".into(),
                )),
                FlashKey::Cvmass => Err(Error::NotImplemented(
                    "Isochoric Specific Heat not valid in two phase region".into(),
                )),
                FlashKey::SpeedSound => Err(Error::NotImplemented(
                    "Speed of Sound not valid in two phase region".into(),
                )),
                FlashKey::Viscosity => Err(Error::NotImplemented(
                    "Viscosity not valid in two phase region".into(),
                )),
                FlashKey::Conductivity => Err(Error::NotImplemented(
                    "Conductivity not valid in two phase region".into(),
                )),
                // Surface tension is not a phase-weighted property.
                FlashKey::SurfaceTension => if97::sigma97(st.t),
                // keyed_output(iPrandtl) resolves through the NON-virtual base
                // AbstractState::Prandtl() = viscosity()*cpmass()/conductivity()
                // — the backend's calc_Flash(iPrandtl) shadow with its
                // "Prandtl number is not valid..." message is dead code for
                // the string API — and the cpmass throw surfaces first
                // (wheel-confirmed message).
                FlashKey::Prandtl => Err(Error::NotImplemented(
                    "Isobaric Specific Heat not valid in two phase region".into(),
                )),
                // Upstream default: phase-weighted combination.
                FlashKey::Hmass | FlashKey::Smass | FlashKey::Umass => {
                    Ok(st.q * sat_vapor(key, st)? + (1.0 - st.q) * sat_liquid(key, st)?)
                }
            }
        }
    } else {
        // Outside the saturation envelope: let IF97 determine the region.
        match key {
            FlashKey::Dmass => if97::rhomass_tp(st.t, st.p),
            FlashKey::Hmass => if97::hmass_tp(st.t, st.p),
            FlashKey::Smass => if97::smass_tp(st.t, st.p),
            FlashKey::Umass => if97::umass_tp(st.t, st.p),
            FlashKey::Cpmass => if97::cpmass_tp(st.t, st.p),
            FlashKey::Cvmass => if97::cvmass_tp(st.t, st.p),
            FlashKey::SpeedSound => if97::speed_sound_tp(st.t, st.p),
            FlashKey::Viscosity => if97::visc_tp(st.t, st.p),
            FlashKey::Conductivity => if97::tcond_tp(st.t, st.p),
            FlashKey::SurfaceTension => Err(Error::NotImplemented(
                "Surface Tension is only valid within the two phase region; Try PQ or QT inputs"
                    .into(),
            )),
            FlashKey::Prandtl => if97::prandtl_tp(st.t, st.p),
        }
    }
}

/// Upstream `calc_SatLiquid` (IF97Backend.h:292-330).
fn sat_liquid(key: FlashKey, st: &If97State) -> Result<f64> {
    match key {
        FlashKey::Dmass => if97::rholiq_p(st.p),
        FlashKey::Hmass => if97::hliq_p(st.p),
        FlashKey::Smass => if97::sliq_p(st.p),
        FlashKey::Cpmass => if97::cpliq_p(st.p),
        FlashKey::Cvmass => if97::cvliq_p(st.p),
        FlashKey::Umass => if97::uliq_p(st.p),
        FlashKey::SpeedSound => if97::speed_soundliq_p(st.p),
        FlashKey::Viscosity => if97::viscliq_p(st.p),
        FlashKey::Conductivity => if97::tcondliq_p(st.p),
        FlashKey::SurfaceTension => if97::sigma97(st.t),
        FlashKey::Prandtl => if97::prandtlliq_p(st.p),
    }
}

/// Upstream `calc_SatVapor` (IF97Backend.h:331-369).
fn sat_vapor(key: FlashKey, st: &If97State) -> Result<f64> {
    match key {
        FlashKey::Dmass => if97::rhovap_p(st.p),
        FlashKey::Hmass => if97::hvap_p(st.p),
        FlashKey::Smass => if97::svap_p(st.p),
        FlashKey::Cpmass => if97::cpvap_p(st.p),
        FlashKey::Cvmass => if97::cvvap_p(st.p),
        FlashKey::Umass => if97::uvap_p(st.p),
        FlashKey::SpeedSound => if97::speed_soundvap_p(st.p),
        FlashKey::Viscosity => if97::viscvap_p(st.p),
        FlashKey::Conductivity => if97::tcondvap_p(st.p),
        FlashKey::SurfaceTension => if97::sigma97(st.t),
        FlashKey::Prandtl => if97::prandtlvap_p(st.p),
    }
}

// ---------------------------------------------------------------------------
// Wheel-derived unit tests. Every expected value below is a verbatim
// CoolProp 8.0.0 wheel result (PropsSI(..., "IF97::Water"), probed 2026-08).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn p(output: Param, n1: Param, v1: f64, n2: Param, v2: f64) -> Result<f64> {
        props(output, n1, v1, n2, v2)
    }

    #[track_caller]
    fn assert_rel(actual: Result<f64>, expected: f64) {
        let actual = actual.expect("expected a value");
        let rel = ((actual - expected) / expected).abs();
        assert!(
            rel <= 1e-11,
            "got {actual:e}, wheel says {expected:e} (rel {rel:e})"
        );
    }

    #[track_caller]
    fn assert_err(actual: Result<f64>, expected: &Error) {
        match actual {
            Err(e) => assert_eq!(&e, expected),
            Ok(v) => panic!("expected {expected:?}, got Ok({v})"),
        }
    }

    fn dome_err(what: &str) -> Error {
        Error::NotImplemented(format!("{what} not valid in two phase region"))
    }

    // ---- Defect 1: the generic outputs-are-inputs echo --------------------

    #[test]
    fn echo_returns_raw_inputs_without_state_update() {
        // Pairs the backend cannot serve still echo.
        assert_eq!(p(Param::Dmass, Param::Dmass, 5.0, Param::T, 300.0), Ok(5.0));
        assert_eq!(
            p(Param::Umass, Param::Dmass, 5.0, Param::Umass, 1e5),
            Ok(1e5)
        );
        // Molar pair echo.
        assert_eq!(p(Param::Hmolar, Param::Hmolar, 5e4, Param::P, 1e5), Ok(5e4));
        // Echo of the second input.
        assert_eq!(p(Param::P, Param::Dmass, 5.0, Param::P, 1e5), Ok(1e5));
        // States no flash could reach: T = 1e6 K would throw in psat97.
        assert_eq!(p(Param::Q, Param::T, 1e6, Param::Q, 0.5), Ok(0.5));
        // Out-of-range Q echoes before the PQ range guard.
        assert_eq!(p(Param::Q, Param::P, 101325.0, Param::Q, 5.0), Ok(5.0));
        // In-dome raw h echoes, not the lever reconstruction.
        assert_eq!(
            p(Param::Hmass, Param::Hmass, 1.5e6, Param::P, 101325.0),
            Ok(1.5e6)
        );
    }

    // ---- Defect 2: Q on (T, P) is the -1 sentinel -------------------------

    #[test]
    fn q_from_tp_is_minus_one() {
        assert_eq!(p(Param::Q, Param::T, 300.0, Param::P, 101325.0), Ok(-1.0));
        // Exactly on the saturation curve: psat97(453.03565995709793) =
        // 1000000.6343433625; the wheel still serves Q = -1.0 (the other
        // outputs raise the Region-4 error).
        assert_eq!(
            p(
                Param::Q,
                Param::T,
                453.03565995709793,
                Param::P,
                1000000.6343433625
            ),
            Ok(-1.0)
        );
        assert_eq!(p(Param::Q, Param::T, 700.0, Param::P, 30e6), Ok(-1.0));
    }

    // ---- Defect 3: (H,P) ---------------------------------------------------

    #[test]
    fn hp_single_phase_liquid_outputs() {
        let (h, pr) = (100000.0, 101325.0);
        assert_rel(
            p(Param::Umass, Param::Hmass, h, Param::P, pr),
            99995.2890069095,
        );
        assert_rel(
            p(Param::Cpmass, Param::Hmass, h, Param::P, pr),
            4182.464480452732,
        );
        assert_rel(
            p(Param::Cvmass, Param::Hmass, h, Param::P, pr),
            4142.773686935153,
        );
        assert_rel(
            p(Param::Viscosity, Param::Hmass, h, Param::P, pr),
            0.0009139632319675285,
        );
        assert_rel(
            p(Param::Conductivity, Param::Hmass, h, Param::P, pr),
            0.6046091693267661,
        );
        assert_rel(
            p(Param::SpeedSound, Param::Hmass, h, Param::P, pr),
            1494.9390333682475,
        );
        assert_rel(
            p(Param::Prandtl, Param::Hmass, h, Param::P, pr),
            6.322462423784383,
        );
        // Molar wrappers: Hmolar is the FORWARD ROUND TRIP h(T(p,h), p) * M,
        // not raw h * M (raw would be 1803.2694... here).
        assert_rel(
            p(Param::Hmolar, Param::Hmass, h, Param::P, pr),
            1803.2722000981273,
        );
        assert_rel(
            p(Param::Smolar, Param::Hmass, h, Param::P, pr),
            6.323206728608124,
        );
        assert_rel(
            p(Param::Umolar, Param::Hmass, h, Param::P, pr),
            1801.4419301969285,
        );
        assert_rel(
            p(Param::Cpmolar, Param::Hmass, h, Param::P, pr),
            75.34821851583672,
        );
        assert_rel(
            p(Param::Cvmolar, Param::Hmass, h, Param::P, pr),
            74.63317823348488,
        );
        assert_rel(
            p(Param::Dmolar, Param::Hmass, h, Param::P, pr),
            55360.68747763173,
        );
        assert_eq!(p(Param::Phase, Param::Hmass, h, Param::P, pr), Ok(0.0));
        assert_err(
            p(Param::SurfaceTension, Param::Hmass, h, Param::P, pr),
            &Error::NotImplemented(
                "Surface Tension is only valid within the two phase region; Try PQ or QT inputs"
                    .into(),
            ),
        );
    }

    #[test]
    fn hp_dome_levers_and_refusals() {
        let (h, pr) = (1.5e6, 101325.0);
        assert_rel(
            p(Param::Umass, Param::Hmass, h, Param::P, pr),
            1418722.5881162882,
        );
        assert_rel(
            p(Param::Dmass, Param::Hmass, h, Param::P, pr),
            1.2466563298665507,
        );
        assert_rel(
            p(Param::Q, Param::Hmass, h, Param::P, pr),
            0.4790559545800459,
        );
        // Hmolar levers back through calc_Flash: exactly h * M.
        assert_rel(
            p(Param::Hmolar, Param::Hmass, h, Param::P, pr),
            27022.902000000002,
        );
        // Surface tension serves in the dome (sigma97 at the backward T).
        assert_rel(
            p(Param::SurfaceTension, Param::Hmass, h, Param::P, pr),
            0.0589168215843171,
        );
        assert_eq!(p(Param::Phase, Param::Hmass, h, Param::P, pr), Ok(6.0));
        assert_err(
            p(Param::Cpmass, Param::Hmass, h, Param::P, pr),
            &dome_err("Isobaric Specific Heat"),
        );
        assert_err(
            p(Param::Cvmass, Param::Hmass, h, Param::P, pr),
            &dome_err("Isochoric Specific Heat"),
        );
        assert_err(
            p(Param::SpeedSound, Param::Hmass, h, Param::P, pr),
            &dome_err("Speed of Sound"),
        );
        assert_err(
            p(Param::Viscosity, Param::Hmass, h, Param::P, pr),
            &dome_err("Viscosity"),
        );
        assert_err(
            p(Param::Conductivity, Param::Hmass, h, Param::P, pr),
            &dome_err("Conductivity"),
        );
        assert_err(
            p(Param::Prandtl, Param::Hmass, h, Param::P, pr),
            &Error::NotImplemented("Isobaric Specific Heat not valid in two phase region".into()),
        );
    }

    #[test]
    fn hp_dome_edges_use_the_saturated_branch() {
        // h exactly at hL(101325): Q clamps to 0.0 -> saturated-liquid branch.
        let h_l = 418990.71780418983;
        assert_eq!(p(Param::Q, Param::Hmass, h_l, Param::P, 101325.0), Ok(0.0));
        assert_rel(
            p(Param::Cpmass, Param::Hmass, h_l, Param::P, 101325.0),
            4216.612690426814,
        );
        assert_rel(
            p(Param::Prandtl, Param::Hmass, h_l, Param::P, 101325.0),
            1.7537547048901514,
        );
        // One J/kg inside: Q = 4.43e-7 is beyond the 1e-10 tolerance.
        assert_err(
            p(Param::Cpmass, Param::Hmass, h_l + 1.0, Param::P, 101325.0),
            &dome_err("Isobaric Specific Heat"),
        );
    }

    #[test]
    fn hp_reverse_phase_quirks() {
        // Region 1 above pcrit -> supercritical liquid.
        assert_eq!(p(Param::Phase, Param::Hmass, 1e5, Param::P, 30e6), Ok(3.0));
        // Region 2 above pcrit, T > Tcrit -> supercritical gas.
        assert_eq!(
            p(Param::Phase, Param::Hmass, 3.5e6, Param::P, 30e6),
            Ok(2.0)
        );
        // Region 3 above pcrit: set_phase consults Tsat97(30 MPa), which
        // throws — upstream quirk, every output errors identically.
        assert_err(
            p(Param::Umass, Param::Hmass, 2.2e6, Param::P, 30e6),
            &Error::OutOfRange("Pressure out of range".into()),
        );
    }

    // ---- Defect 3: (P,S) ---------------------------------------------------

    #[test]
    fn ps_outputs() {
        let (pr, s) = (101325.0, 300.0);
        assert_rel(
            p(Param::Umass, Param::P, pr, Param::Smass, s),
            84949.4311057666,
        );
        // calc_hmass reverse short-circuit: hmass_psmass, dome or not.
        assert_rel(
            p(Param::Hmass, Param::P, pr, Param::Smass, s),
            85050.94343470575,
        );
        assert_rel(
            p(Param::Cpmass, Param::P, pr, Param::Smass, s),
            4184.62129999503,
        );
        assert_rel(
            p(Param::Viscosity, Param::P, pr, Param::Smass, s),
            0.000995541431129626,
        );
        // Steam side.
        assert_rel(
            p(Param::SpeedSound, Param::P, pr, Param::Smass, 7500.0),
            490.7564321489067,
        );
        // Dome.
        assert_rel(
            p(Param::Umass, Param::P, pr, Param::Smass, 4000.0),
            1348357.6875674292,
        );
        assert_err(
            p(Param::Cpmass, Param::P, pr, Param::Smass, 4000.0),
            &dome_err("Isobaric Specific Heat"),
        );
        assert_eq!(p(Param::Phase, Param::P, pr, Param::Smass, 4000.0), Ok(6.0));
    }

    // ---- Defect 3: (H,S) ---------------------------------------------------

    #[test]
    fn hs_outputs() {
        // Dome state from the task brief.
        let (h, s) = (2e6, 6000.0);
        assert_rel(
            p(Param::Dmass, Param::Hmass, h, Param::Smass, s),
            0.2213157144039711,
        );
        assert_rel(
            p(Param::T, Param::Hmass, h, Param::Smass, s),
            338.37860224064286,
        );
        assert_rel(
            p(Param::P, Param::Hmass, h, Param::Smass, s),
            25298.265163108077,
        );
        assert_rel(
            p(Param::Q, Param::Hmass, h, Param::Smass, s),
            0.7364863113496587,
        );
        assert_rel(
            p(Param::Umass, Param::Hmass, h, Param::Smass, s),
            1885691.5098358956,
        );
        assert_rel(
            p(Param::SurfaceTension, Param::Hmass, h, Param::Smass, s),
            0.06532571554898685,
        );
        assert_err(
            p(Param::Viscosity, Param::Hmass, h, Param::Smass, s),
            &dome_err("Viscosity"),
        );
        // Single-phase liquid: Q sentinel and the forward round trips.
        assert_eq!(
            p(Param::Q, Param::Hmass, 1e5, Param::Smass, 350.0),
            Ok(-1.0)
        );
        assert_rel(
            p(Param::Umass, Param::Hmass, 1e5, Param::Smass, 350.0),
            99796.36373551821,
        );
        assert_rel(
            p(Param::Viscosity, Param::Hmass, 1e5, Param::Smass, 350.0),
            0.0009148637092052269,
        );
        // No _reverse on (H,S): Hmolar is the calc_Flash round trip.
        assert_rel(
            p(Param::Hmolar, Param::Hmass, 1e5, Param::Smass, 350.0),
            1803.2480547058053,
        );
    }

    #[test]
    fn hs_forward_phase_classification() {
        // (h, s) generated from wheel PT states.
        // (650 K, 25 MPa): supercritical.
        assert_eq!(
            p(
                Param::Phase,
                Param::Hmass,
                1876359.1164090303,
                Param::Smass,
                4075.9789901397203
            ),
            Ok(1.0)
        );
        // (700 K, 0.1 MPa): supercritical gas.
        assert_eq!(
            p(
                Param::Phase,
                Param::Hmass,
                3334336.1008003005,
                Param::Smass,
                8626.349167486484
            ),
            Ok(2.0)
        );
        // (300 K, 30 MPa): supercritical liquid.
        assert_eq!(
            p(
                Param::Phase,
                Param::Hmass,
                139885.036558159,
                Param::Smass,
                384.4768076233224
            ),
            Ok(3.0)
        );
        // Subcritical liquid / steam.
        assert_eq!(
            p(Param::Phase, Param::Hmass, 1e5, Param::Smass, 350.0),
            Ok(0.0)
        );
        assert_eq!(
            p(Param::Phase, Param::Hmass, 2.8e6, Param::Smass, 7500.0),
            Ok(5.0)
        );
    }

    // ---- Defect 3: (P,Q) and the QT route ---------------------------------

    #[test]
    fn pq_saturated_branches_and_dome() {
        let pr = 101325.0;
        assert_rel(
            p(Param::Prandtl, Param::P, pr, Param::Q, 0.0),
            1.7537547048901514,
        );
        assert_rel(
            p(Param::Prandtl, Param::P, pr, Param::Q, 1.0),
            1.0342483422571493,
        );
        assert_rel(
            p(Param::Umass, Param::P, pr, Param::Q, 0.0),
            418884.99171568773,
        );
        // Q within 1e-10 of an endpoint takes the saturated branch, not the
        // lever (wheel: identical to Q = 0 for every output).
        assert_rel(
            p(Param::Cpmass, Param::P, pr, Param::Q, 1e-11),
            4216.612690426814,
        );
        assert_rel(
            p(Param::Dmass, Param::P, pr, Param::Q, 1e-11),
            958.3727293380055,
        );
        // Mid-dome.
        assert_rel(
            p(Param::Umass, Param::P, pr, Param::Q, 0.5),
            1462434.9015347678,
        );
        assert_rel(
            p(Param::Dmass, Param::P, pr, Param::Q, 0.5),
            1.1945013625647232,
        );
        assert_err(
            p(Param::Prandtl, Param::P, pr, Param::Q, 0.5),
            &Error::NotImplemented("Isobaric Specific Heat not valid in two phase region".into()),
        );
        assert_eq!(p(Param::Phase, Param::P, pr, Param::Q, 0.5), Ok(6.0));
        // Range guard (message verbatim; the echo above bypasses it only for
        // outputs that ARE inputs).
        assert_err(
            p(Param::T, Param::P, pr, Param::Q, 5.0),
            &Error::OutOfRange("Input vapor quality [Q] must be between 0 and 1".into()),
        );
        assert_err(
            p(Param::Dmass, Param::Q, 5.0, Param::T, 300.0),
            &Error::OutOfRange("Input vapor quality [Q] must be between 0 and 1".into()),
        );
    }

    #[test]
    fn qt_uses_the_raw_input_temperature() {
        // sigma97 evaluates at the INPUT T, not the Tsat97(psat97(T)) round
        // trip (which the PQ route would give: 0.0589168215843171).
        assert_rel(
            p(
                Param::SurfaceTension,
                Param::Q,
                0.5,
                Param::T,
                373.124295847684,
            ),
            0.058916822384634644,
        );
        assert_rel(
            p(Param::Prandtl, Param::Q, 0.0, Param::T, 300.0),
            5.857530906793436,
        );
        assert_rel(
            p(Param::Umass, Param::Q, 0.5, Param::T, 300.0),
            1262123.6282264711,
        );
        assert_rel(
            p(Param::Hmolar, Param::Q, 0.5, Param::T, 300.0),
            23982.537272782884,
        );
    }

    // ---- (T,P): full output family, phase, update-time errors --------------

    #[test]
    fn pt_serves_molar_family_and_verbatim_errors() {
        let (t, pr) = (300.0, 101325.0);
        assert_rel(
            p(Param::Dmolar, Param::T, t, Param::P, pr),
            55317.41609929838,
        );
        assert_rel(
            p(Param::Hmolar, Param::T, t, Param::P, pr),
            2029.6909514166302,
        );
        assert_rel(
            p(Param::Smolar, Param::T, t, Param::P, pr),
            7.081742084427582,
        );
        assert_rel(
            p(Param::Umolar, Param::T, t, Param::P, pr),
            2027.859249808983,
        );
        assert_rel(
            p(Param::Cpmolar, Param::T, t, Param::P, pr),
            75.32358887591938,
        );
        assert_rel(
            p(Param::Cvmolar, Param::T, t, Param::P, pr),
            74.41539910700445,
        );
        assert_err(
            p(Param::SurfaceTension, Param::T, t, Param::P, pr),
            &Error::NotImplemented(
                "Surface Tension is only valid within the two phase region; Try PQ or QT inputs"
                    .into(),
            ),
        );
    }

    #[test]
    fn pt_phase_classification_matches_the_shipped_wheel() {
        let ph = |t, pr| p(Param::Phase, Param::T, t, Param::P, pr);
        assert_eq!(ph(300.0, 101325.0), Ok(0.0)); // liquid
        assert_eq!(ph(400.0, 1e5), Ok(5.0)); // gas
        assert_eq!(ph(650.0, 25e6), Ok(1.0)); // supercritical (region 3)
        assert_eq!(ph(700.0, 30e6), Ok(1.0)); // supercritical
        assert_eq!(ph(300.0, 30e6), Ok(3.0)); // supercritical liquid
        assert_eq!(ph(1500.0, 1e6), Ok(2.0)); // supercritical gas (region 5)
        assert_eq!(ph(640.0, 20e6), Ok(5.0)); // region-3 vapor side, T < Tc
        // Phase classifies even where the forward evaluators throw.
        assert_eq!(ph(12345.0, 1e6), Ok(2.0));
        assert_eq!(ph(300.0, 200e6), Ok(3.0));
        // Saturation HALF-band: on-curve and in-band are LIQUID (never
        // twophase, no throw); gas starts below psat*(1 - 3.3e-5).
        let (t, psat) = (453.03565995709793, 1000000.6343433625);
        assert_eq!(ph(t, psat), Ok(0.0));
        assert_eq!(ph(t, psat * (1.0 + 1e-5)), Ok(0.0));
        assert_eq!(ph(t, psat * (1.0 - 1e-5)), Ok(0.0));
        assert_eq!(ph(t, psat * (1.0 - 3.4e-5)), Ok(5.0));
        // Critical-point band: |T-Tc| < 3.3e-6 && |p-Pc| < 3.3e-5 (absolute).
        assert_eq!(ph(647.096, 22.064e6), Ok(4.0));
        assert_eq!(ph(647.096 + 3e-6, 22.064e6), Ok(4.0));
        assert_eq!(ph(647.096 + 3.4e-6, 22.064e6), Ok(2.0)); // T > Tc, p == Pc
        assert_eq!(ph(647.096, 22.064e6 + 1.0), Ok(3.0)); // T <= Tc, p > Pc
        assert_eq!(ph(647.096, 22.064e6 - 1e5), Ok(5.0)); // gas side
    }

    #[test]
    fn pt_update_time_psat_throw_poisons_every_output() {
        // T = 200 K < Tmin: set_phase's psat97 throws at update, so even Q
        // and Phase error (wheel-confirmed for all outputs incl. Z/Gmass).
        let expected = Error::OutOfRange("Temperature out of range".into());
        assert_err(p(Param::Q, Param::T, 200.0, Param::P, 1e5), &expected);
        assert_err(p(Param::Phase, Param::T, 200.0, Param::P, 1e5), &expected);
        assert_err(p(Param::Dmass, Param::T, 200.0, Param::P, 1e5), &expected);
        // T > Tcrit skips the psat call: Q/Phase serve while Dmass throws.
        assert_eq!(p(Param::Q, Param::T, 12345.0, Param::P, 1e6), Ok(-1.0));
        assert_err(p(Param::Dmass, Param::T, 12345.0, Param::P, 1e6), &expected);
        // Exactly on the saturation curve the update succeeds (liquid) and
        // the forward evaluators raise the Region-4 error instead.
        assert_err(
            p(
                Param::Dmass,
                Param::T,
                453.03565995709793,
                Param::P,
                1000000.6343433625,
            ),
            &Error::OutOfRange("Cannot use Region 4 with T and p as inputs".into()),
        );
    }

    #[test]
    fn qmass_semantics() {
        // Single-phase sentinel _Q = -1 -> range ValueError (verbatim).
        assert_err(
            p(Param::Qmass, Param::T, 300.0, Param::P, 101325.0),
            &Error::Value("Qmass requires a two-phase state (0 <= Q <= 1)".into()),
        );
        assert_err(
            p(Param::Qmass, Param::Hmass, 1e5, Param::P, 101325.0),
            &Error::Value("Qmass requires a two-phase state (0 <= Q <= 1)".into()),
        );
        // EXACT endpoints return _Q; the clamped HP dome edge counts.
        assert_eq!(p(Param::Qmass, Param::P, 101325.0, Param::Q, 0.0), Ok(0.0));
        assert_eq!(p(Param::Qmass, Param::P, 101325.0, Param::Q, 1.0), Ok(1.0));
        assert_eq!(
            p(
                Param::Qmass,
                Param::Hmass,
                418990.71780418983,
                Param::P,
                101325.0
            ),
            Ok(0.0)
        );
        // Any other quality (1e-11 included — no 1e-10 band here) reaches the
        // backend's get_mole_fractions override, which throws.
        let mole_err = Error::NotImplemented(
            "get_mole_fractions composition has not been implemented.".into(),
        );
        assert_err(
            p(Param::Qmass, Param::P, 101325.0, Param::Q, 0.5),
            &mole_err,
        );
        assert_err(
            p(Param::Qmass, Param::P, 101325.0, Param::Q, 1e-11),
            &mole_err,
        );
    }

    // ---- Molar input pairs and pair-error parity ---------------------------

    #[test]
    fn molar_input_pairs_convert_to_mass() {
        assert_rel(
            p(Param::T, Param::Hmolar, 5e4, Param::P, 1e5),
            422.5633348137442,
        );
        // Non-echoed Hmass on (Hmolar, P) is the forward round trip, not the
        // raw conversion (2775423.6018026485).
        assert_rel(
            p(Param::Hmass, Param::Hmolar, 5e4, Param::P, 1e5),
            2775426.789249523,
        );
        assert_rel(
            p(Param::T, Param::P, 1e5, Param::Smolar, 100.0),
            372.75591861133773,
        );
        // (Hmolar, Smolar) does NOT convert — upstream's HmassSmass
        // fall-through clobbers the converted values with the raw molar ones
        // (defect reproduced): 281.154 K is T(h = 5e4 J/kg, s = 120 J/kg/K).
        assert_rel(
            p(Param::T, Param::Hmolar, 5e4, Param::Smolar, 120.0),
            281.15430258287824,
        );
    }

    #[test]
    fn unsupported_pairs_match_upstreams_value_error() {
        let expected = Error::Value("This pair of inputs is not yet supported".into());
        assert_err(p(Param::T, Param::Dmass, 900.0, Param::P, 1e5), &expected);
        assert_err(
            p(Param::Dmass, Param::Qmass, 0.5, Param::T, 300.0),
            &expected,
        );
        assert_err(
            p(Param::T, Param::Smass, 300.0, Param::Umass, 5e4),
            &expected,
        );
    }

    // ---- Trivial additions -------------------------------------------------

    #[test]
    fn trivial_gas_constant_and_rhomolar_critical() {
        assert_rel(
            p(Param::GasConstant, Param::T, 300.0, Param::P, 101325.0),
            8.314514578968002,
        );
        assert_rel(
            p(Param::RhomolarCritical, Param::T, 300.0, Param::P, 101325.0),
            17873.727995609057,
        );
    }
}
