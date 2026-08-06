# rustprop

Pure-Rust port of [CoolProp 8](https://github.com/CoolProp/CoolProp) — thermophysical property calculations — built for WebAssembly with **modular, opt-in calculation engines** so every application compiles in only the parts it actually needs.

Fidelity is the rule: the same algorithms and the same fluid data as upstream CoolProp v8.0.0, validated by tests against upstream results.

## Workspace layout

| Crate | Role |
|---|---|
| `crates/rustprop-core` | Shared foundations: parameter/state/error **types**, engine traits. Fluid-data *types* live here so engines never depend on the data crate. |
| `crates/rustprop-data` | Generated fluid/mixture data **contents**, one Cargo feature per fluid (data dominates WASM size — always opt-in). |
| `crates/rustprop-heos` | HEOS engine — multiparameter Helmholtz EOS (upstream `src/Backends/Helmholtz`). |
| `crates/rustprop-cubics` | SRK / Peng-Robinson with Chebyshev superancillaries (upstream `src/Backends/Cubics`). |
| `crates/rustprop-if97` | IAPWS-IF97 industrial water/steam formulation (upstream `src/Backends/IF97`). |
| `crates/rustprop-incompressible` | Brines and secondary working fluids (upstream `src/Backends/Incompressible`). |
| `crates/rustprop-pcsaft` | PC-SAFT EOS (upstream `src/Backends/PCSAFT`). |
| `crates/rustprop-tabular` | TTSE / bicubic table interpolation (upstream `src/Backends/Tabular`). |
| `crates/rustprop-svdsbtl` | SVD-compressed tabular lookup, new in v8 (upstream `src/Backends/SVDSBTL`). |
| `crates/rustprop-humid-air` | Humid air / psychrometrics (`HAPropsSI`). |
| `crates/rustprop` | Facade: `PropsSI`-style API; every engine behind a Cargo feature, `default = []`. |
| `apps/rustprop-cli` | Example CLI exposing the libraries and calculations over stdout. |
| `tools/rustprop-datagen` | Codegen: upstream CoolProp JSON → Rust data modules in `rustprop-data`. |

A `rustprop-wasm` bindings crate (wasm-bindgen, prebuilt `.wasm` release artifacts) is planned once the first engine works.

## Building

`~/.cargo/bin` may not be on `PATH` in this environment:

```bash
export PATH="$HOME/.cargo/bin:$PATH"

cargo build                                           # native
cargo test                                            # all tests
cargo run -p rustprop-cli                             # example CLI
cargo build -p rustprop --features all-backends \
    --target wasm32-unknown-unknown --release         # wasm facade
```

## Status

Scaffold only — no algorithms ported yet. The porting plan, fidelity rules, and upstream mapping live in `CLAUDE.md`.

## License

MIT. Derivative work of CoolProp (MIT, © 2012–2018 Ian H. Bell and other CoolProp developers) — see `LICENSE`.
