//! The packed fluid format (`rustprop_core::pack`): every fluid and the
//! mixture tables decode bit for bit to the linked data, decoded fluids
//! compute the same states, and malformed input is refused, never a panic.

use rustprop_core::fluid::FluidData;
use rustprop_core::pack::{self, MixtureTables};
use rustprop_heos::flash_pt::PtFlash;

/// `encode` writes every field (it destructures without `..`), so a decoded
/// value that re-encodes to the same bytes equals the original field for
/// field, NaN payloads and infinities included.
#[test]
fn every_fluid_round_trips_bit_for_bit() {
    let all = rustprop_data::fluids::all();
    assert_eq!(all.len(), 136);
    for (name, data) in all {
        let bytes = pack::encode_fluid(data);
        assert_eq!(pack::fluid_name(&bytes).unwrap(), name);
        let decoded = pack::decode_fluid(&bytes).unwrap();
        assert_eq!(pack::encode_fluid(decoded), bytes, "{name}");
    }
}

#[test]
fn the_mixture_tables_round_trip_bit_for_bit() {
    let linked = MixtureTables {
        binary_pairs: rustprop_data::mixtures::MIX_BINARY_PAIRS,
        departure_fns: rustprop_data::mixtures::MIX_DEPARTURE_FNS,
        predefined: rustprop_data::mixtures::MIX_PREDEFINED,
    };
    let bytes = pack::encode_mixtures(&linked);
    let decoded = pack::decode_mixtures(&bytes).unwrap();
    assert_eq!(decoded.binary_pairs.len(), 888);
    assert_eq!(decoded.departure_fns.len(), 28);
    assert_eq!(decoded.predefined.len(), 154);
    assert_eq!(pack::encode_mixtures(decoded), bytes);
}

/// The same code on the same data is deterministic, so a decoded fluid's
/// states equal the linked fluid's bitwise (refusals included).
#[test]
fn decoded_fluids_compute_the_same_states() {
    let probe = |data: &'static FluidData| {
        let flash = PtFlash::new(data);
        let (tc, pc) = (data.states.critical.t, data.states.critical.p);
        let sat = flash
            .qt_state(0.8 * tc, 0.0)
            .map(|s| {
                (
                    s.p().to_bits(),
                    s.rhomolar().to_bits(),
                    flash.state_hmolar(&s).to_bits(),
                )
            })
            .map_err(|e| e.to_string());
        let pt = flash
            .pt_flash(1.1 * tc, 1.5 * pc)
            .map(|(rho, _)| rho.to_bits())
            .map_err(|e| e.to_string());
        (sat, pt)
    };
    for (name, data) in rustprop_data::fluids::all() {
        let decoded = pack::decode_fluid(&pack::encode_fluid(data)).unwrap();
        assert_eq!(probe(decoded), probe(data), "{name}");
    }
}

#[test]
fn malformed_input_is_refused() {
    let bytes = pack::encode_fluid(&rustprop_data::fluids::air::AIR);
    // Every truncation, including an empty input.
    for n in 0..bytes.len() {
        assert!(pack::decode_fluid(&bytes[..n]).is_err(), "prefix {n}");
    }
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(pack::decode_fluid(&trailing).is_err());
    // The other record kind, and a format version this build does not read.
    assert!(pack::decode_mixtures(&bytes).is_err());
    let mut version = bytes.clone();
    version[7] = b'9';
    assert!(pack::decode_fluid(&version).is_err());
    // A name length far past the data is refused before any allocation.
    let mut long = bytes.clone();
    long[8..12].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(pack::fluid_name(&long).is_err());
    assert!(pack::decode_fluid(&long).is_err());
}
