//! Data pipeline: reads fluid/mixture JSON from an upstream CoolProp v8.0.0
//! checkout and emits per-fluid, feature-gated Rust modules into
//! `crates/rustprop-data`. Runs at development time only — JSON parsing never
//! reaches the shipped library or wasm binaries.

fn main() {
    eprintln!("rustprop-datagen: not implemented yet.");
    eprintln!("Planned: <coolprop-checkout>/dev/fluids/*.json -> crates/rustprop-data/src/");
    std::process::exit(1);
}
