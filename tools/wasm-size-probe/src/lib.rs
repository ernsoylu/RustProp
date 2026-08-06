//! Minimal cdylib whose compiled `.wasm` size CI reports per feature set.
//!
//! The linker strips unreachable code, so the numbers only mean something if
//! `probe()` actually calls into every enabled engine's public entry points —
//! extend it under `#[cfg(feature = ...)]` as facade API lands (first: the
//! IF97 baseline, PLAN.md step 2.6).

// The exported symbol needs `#[unsafe(no_mangle)]` (edition 2024), which the
// workspace-wide `deny(unsafe_code)` would reject; this probe is the one
// deliberate exception and ships nowhere.
#![allow(unsafe_code)]

#[unsafe(no_mangle)]
pub extern "C" fn probe() -> usize {
    rustprop::UPSTREAM_VERSION.len()
}
