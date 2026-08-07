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
| `crates/rustprop-wasm` | wasm-bindgen bindings: JS-facing `props_si` / `ha_props_si` plus a `Float64Array` batch path. |
| `apps/rustprop-cli` | Example CLI exposing the libraries and calculations over stdout. |
| `tools/rustprop-datagen` | Codegen: upstream CoolProp JSON → Rust data modules in `rustprop-data`. |
| `tools/rustprop-svdgen` | Converts upstream `.svd.bin.z` surfaces into the flat `.svds` blobs `rustprop-svdsbtl` reads. |

### WASM bundle sizes

Engine and fluid selection is a compile-time choice, so a bundle carries only what
the application asked for. Measured bytes per feature set are in
[WASM-SIZES.md](WASM-SIZES.md) — IF97 alone is 124 KB, HEOS with Water 304 KB, and
everything at once 4.2 MB. Regenerate with `tools/wasm-size-table.sh`.

```bash
wasm-pack build crates/rustprop-wasm --target web --features heos,water
```

## Project status

All fifteen phases of `PLAN.md` are complete — every engine ported, ~33,900
committed oracle records, CI green. See **[NEXT-STEPS.md](NEXT-STEPS.md)** for
current status, known divergences from upstream, and what to work on next.

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

## Quickstart

`PropsSI` semantics, from the CLI (any of the 130 ported pure fluids, or `IF97::Water`):

```bash
$ cargo run -p rustprop-cli -- props Dmolar T 300 P 101325 Water
55317.35277350119
$ cargo run -p rustprop-cli -- props H T 300 P 101325 IF97::Water
112665.04341853978
```

or from Rust (features select the engines and, per fluid, the data your binary carries —
`heos` + `rustprop-data/water` compiles to ~136 KB of wasm):

```rust
// Equivalent to PropsSI("Dmolar", "T", 300, "P", 101325, "Water")
let d = rustprop::props_si("Dmolar", "T", 300.0, "P", 101325.0, "Water")?;
```

## Status

Golden-verified against the CoolProp 8.0.0 oracle wheel (~12,000 committed records):

- **IF97** steam engine (356 records at rtol 1e-11; all IAPWS published check tables pass).
- **HEOS** pure fluids: all 130 superancillary fluids generated, bitwise data-fidelity-walked,
  and smoke-tested; 12 fluids carry full per-suite batteries (Helmholtz terms through every
  flash pair, including (H,S)); all eight input pairs; `PropsSI`-style string API with
  mass-basis aliases, trivial outputs, and upstream error conditions.

The porting plan, fidelity rules, and upstream mapping live in `CLAUDE.md` and `PLAN.md`.

## License

MIT. Derivative work of CoolProp (MIT, © 2012–2018 Ian H. Bell and other CoolProp developers) — see `LICENSE`.
