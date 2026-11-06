//! A build with the mixture engine but without `mixture-data`, before any
//! tables are installed: no predefined blend exists, and an explicit mixture
//! is refused as unavailable rather than computed from missing tables.
//! Runs only in that configuration (CI: `-p rustprop --features
//! heos-mixture-engine,rustprop-data/water`); the workspace build links the
//! tables. The installed path is covered end to end by frees-wasm.
#![cfg(all(feature = "heos-mixture-engine", not(feature = "mixture-data")))]

use rustprop::props_api::install_mixture_tables;
use rustprop::{Error, props_si};

#[test]
fn without_tables_there_are_no_blends_and_mixtures_are_unavailable() {
    assert!(props_si("T", "P", 101_325.0, "Q", 0.0, "R410A.mix").is_err());
    let mixed = props_si("T", "P", 101_325.0, "Q", 0.0, "Water[0.6]&Water[0.4]");
    assert!(matches!(mixed, Err(Error::NotAvailable(_))), "{mixed:?}");
    assert!(install_mixture_tables(b"RPMIXTB1 truncated").is_err());
    // A failed install leaves the build without tables.
    let mixed = props_si("T", "P", 101_325.0, "Q", 0.0, "Water[0.6]&Water[0.4]");
    assert!(matches!(mixed, Err(Error::NotAvailable(_))));
}
