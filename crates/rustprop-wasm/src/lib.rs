//! JS bindings for rustprop.
//!
//! Which engines and which fluids exist in the built `.wasm` is decided at
//! COMPILE time by Cargo features, not at runtime — that is the whole point
//! of this project, and it is why there is no "load a fluid" call here. A
//! bundle built with `--features heos,water` carries the Helmholtz engine and
//! Water's coefficients and nothing else.
//!
//! ```bash
//! wasm-pack build crates/rustprop-wasm --target web --features if97
//! ```
//!
//! ```js
//! import init, { props_si } from "./pkg/rustprop_wasm.js";
//! await init();
//! props_si("D", "T", 400, "P", 101325, "IF97::Water");  // 0.5548...
//! ```
//!
//! Errors surface as thrown JS exceptions carrying rustprop's message, so
//! `try`/`catch` works the way a JS caller expects.

use wasm_bindgen::prelude::*;

/// The CoolProp release this port tracks.
#[wasm_bindgen]
pub fn upstream_version() -> String {
    rustprop::UPSTREAM_VERSION.to_string()
}

/// `PropsSI(output, name1, val1, name2, val2, fluid)`.
// Mirrors the facade's own gate on `props_si`: the string API exists when at
// least one engine that answers it is compiled in. A humid-air-only bundle
// legitimately exports `ha_props_si` and nothing else.
#[cfg(any(
    feature = "heos",
    feature = "if97",
    feature = "cubics",
    feature = "incompressible",
    feature = "pcsaft"
))]
#[wasm_bindgen]
pub fn props_si(
    output: &str,
    name1: &str,
    val1: f64,
    name2: &str,
    val2: f64,
    fluid: &str,
) -> Result<f64, JsError> {
    rustprop::props_si(output, name1, val1, name2, val2, fluid).map_err(to_js)
}

/// The typed fast path: one output over many states, with both input vectors
/// crossing the boundary as `Float64Array`s.
///
/// A JS caller paying per-call string marshalling dominates its own runtime
/// long before the thermodynamics does; this keeps the strings on one side
/// and the numbers on the other. `vals1` and `vals2` must be the same length.
/// A state that fails yields NaN in its slot rather than aborting the batch —
/// a sweep over a grid routinely clips the phase envelope, and one bad cell
/// should not cost the caller the other ten thousand.
// Mirrors the facade's own gate on `props_si`: the string API exists when at
// least one engine that answers it is compiled in. A humid-air-only bundle
// legitimately exports `ha_props_si` and nothing else.
#[cfg(any(
    feature = "heos",
    feature = "if97",
    feature = "cubics",
    feature = "incompressible",
    feature = "pcsaft"
))]
#[wasm_bindgen]
pub fn props_si_many(
    output: &str,
    name1: &str,
    vals1: &[f64],
    name2: &str,
    vals2: &[f64],
    fluid: &str,
) -> Result<Vec<f64>, JsError> {
    if vals1.len() != vals2.len() {
        return Err(JsError::new(&format!(
            "props_si_many: vals1 has {} entries, vals2 has {}",
            vals1.len(),
            vals2.len()
        )));
    }
    Ok(vals1
        .iter()
        .zip(vals2.iter())
        .map(|(&v1, &v2)| {
            rustprop::props_si(output, name1, v1, name2, v2, fluid).unwrap_or(f64::NAN)
        })
        .collect())
}

/// `HAPropsSI(output, n1, v1, n2, v2, n3, v3)` — humid air.
#[cfg(feature = "humid-air")]
#[wasm_bindgen]
pub fn ha_props_si(
    output: &str,
    name1: &str,
    val1: f64,
    name2: &str,
    val2: f64,
    name3: &str,
    val3: f64,
) -> Result<f64, JsError> {
    rustprop::ha_props_si(output, name1, val1, name2, val2, name3, val3).map_err(to_js)
}

fn to_js(e: rustprop::Error) -> JsError {
    JsError::new(&e.to_string())
}
