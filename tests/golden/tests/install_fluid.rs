//! `props_api::install_fluid`: a fluid installed from packed bytes is served
//! by `props_si` exactly as the linked one is. Water is packed under another
//! name (and no CAS or aliases), so the installed copy is what answers.
//! Its own test binary, because installing changes the process's registry.

use rustprop::props_api::{install_fluid, resolve_fluid};
use rustprop::props_si;
use rustprop_core::pack;

fn str_field(name: &str) -> Vec<u8> {
    let mut v = (name.len() as u32).to_le_bytes().to_vec();
    v.extend_from_slice(name.as_bytes());
    v
}

/// Packed Water renamed: the name, CAS and aliases lead the record.
fn packed_water_as(name: &str) -> Vec<u8> {
    let water = &rustprop_data::fluids::water::WATER;
    let bytes = pack::encode_fluid(water);
    let head = 8
        + 4
        + water.name.len()
        + 4
        + water.cas.len()
        + 4
        + water.aliases.iter().map(|a| 4 + a.len()).sum::<usize>();
    let mut out = bytes[..8].to_vec();
    out.extend(str_field(name));
    out.extend(str_field("0-00-0-packed"));
    out.extend(0u32.to_le_bytes());
    out.extend_from_slice(&bytes[head..]);
    out
}

#[test]
fn an_installed_fluid_answers_as_the_linked_one() {
    assert!(resolve_fluid("PackedWater").is_err());
    assert_eq!(
        install_fluid(&packed_water_as("PackedWater")).unwrap(),
        "PackedWater"
    );
    let installed = resolve_fluid("PackedWater").unwrap();
    assert!(!std::ptr::eq(
        installed,
        &rustprop_data::fluids::water::WATER
    ));

    // A second install keeps the first: the caches are keyed by address.
    assert_eq!(
        install_fluid(&packed_water_as("PackedWater")).unwrap(),
        "PackedWater"
    );
    assert!(std::ptr::eq(
        resolve_fluid("PackedWater").unwrap(),
        installed
    ));

    for (out, n1, v1, n2, v2) in [
        ("Hmass", "T", 300.0, "P", 101_325.0),
        ("Dmass", "T", 600.0, "P", 30e6),
        ("viscosity", "T", 350.0, "P", 1e6),
        ("conductivity", "T", 500.0, "P", 1e5),
        ("P", "T", 373.15, "Q", 0.0),
        ("Smass", "P", 2e6, "Q", 1.0),
    ] {
        let linked = props_si(out, n1, v1, n2, v2, "Water").unwrap();
        let packed = props_si(out, n1, v1, n2, v2, "PackedWater").unwrap();
        assert_eq!(
            packed.to_bits(),
            linked.to_bits(),
            "{out}({n1}={v1}, {n2}={v2})"
        );
    }
}

#[test]
fn installing_a_linked_fluid_keeps_the_linked_data() {
    let water = &rustprop_data::fluids::water::WATER;
    assert_eq!(install_fluid(&pack::encode_fluid(water)).unwrap(), "Water");
    assert!(std::ptr::eq(resolve_fluid("Water").unwrap(), water));
    assert!(install_fluid(b"not a fluid").is_err());
}
