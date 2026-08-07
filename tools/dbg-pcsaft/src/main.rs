use rustprop_data::pcsaft::{PCSAFT_BINARY_PAIRS, PCSAFT_FLUIDS};
use rustprop_pcsaft::{PcsaftBackend, PcsaftInput};

fn main() {
    let w = PCSAFT_FLUIDS.iter().find(|f| f.name == "WATER").unwrap();
    let mut b = PcsaftBackend::new(&[w], PCSAFT_BINARY_PAIRS).unwrap();
    match b.update(PcsaftInput::Pt, 82825.16480101226, 350.0) {
        Ok(()) => println!("WATER PT: rho={} phase={:?}", b.rhomolar, b.phase),
        Err(e) => println!("WATER PT FAIL {e:?}"),
    }
}
