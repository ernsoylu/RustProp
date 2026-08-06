//! Example CLI that makes the rustprop libraries and calculations available
//! over stdout. Grows a `PropsSI`-style command surface as engines land.

fn main() {
    println!(
        "rustprop-cli {} — pure-Rust port of CoolProp {}",
        env!("CARGO_PKG_VERSION"),
        rustprop::UPSTREAM_VERSION,
    );
    println!("No engines are ported yet. Planned PropsSI-style usage:");
    println!("  rustprop-cli props T P 101325 Q 0 Water");
}
