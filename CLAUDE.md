# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

rustprop is a from-scratch port of **CoolProp 8** — the C++ thermophysical property library (upstream: https://github.com/CoolProp/CoolProp, current release v8.0.0, June 2026) — to **pure Rust**, targeting **WebAssembly** for the owner's Rust-based web apps. Repository: https://github.com/ernsoylu/RustProp.

The defining constraint is **modularity for WASM binary size**: all CoolProp data and algorithms get ported, but as independently selectable calculation engines (workspace crates / Cargo features), so application logic compiles in only the specific parts a calculation needs — never a monolithic all-fluids, all-backends binary. Fluid *data* selection matters as much as algorithm selection here: per-fluid JSON data dominates size and must be opt-in too.

- Pure Rust only — no C/C++ FFI (would defeat modular WASM compilation).
- Primary compile target is `wasm32-unknown-unknown`; keep every dependency wasm-compatible (no threads, filesystem, or clock assumptions in core crates). Native targets remain useful for tests.
- **Fidelity is a hard requirement**: implement exactly the algorithms CoolProp 8 implements — no reformulations or "improvements" — and carry over all fluid/mixture data unchanged. Upstream is the correctness oracle; comprehensive tests must validate the ported data and computed property values against upstream results.

## Deliverables

1. The modular engine crates (the core of the project).
2. Comprehensive test suites validating ported data and calculation results against upstream CoolProp values.
3. An example CLI application that makes the libraries and calculations available over stdout.
4. Releases published both as source code and as prebuilt binaries for consumption by other WASM applications.
5. CI/CD pipelines covering build, lint, the full test suite, and release packaging — GitHub Actions. `.github/workflows/ci.yml` exists; the release pipeline (prebuilt wasm artifacts) is still to come.

## Status (as of 2026-08)

Workspace scaffolded — **no algorithms ported yet**. Everything compiles with zero external dependencies; `cargo test` runs a placeholder version test; the wasm facade builds. The crate map is in README.md. Architecture rules the scaffold encodes:

- **Types/contents split**: fluid-data *types* live in `rustprop-core`, generated data *contents* in `rustprop-data` (one Cargo feature per fluid, `default = []`). Engines depend only on core — apps link data solely for the fluids they opt into.
- The facade crate `rustprop` puts every engine behind a Cargo feature with `default = []`; `all-backends` turns everything on (used by the CLI and CI).
- `rustprop-data` contents come only from `tools/rustprop-datagen` codegen, never hand edits; JSON parsing must stay out of shipped binaries.
- Workspace lints deny `unsafe_code`; release profile uses fat LTO, `panic = "abort"`, symbol stripping.

Natural next steps: implement `rustprop-datagen` against an upstream v8.0.0 checkout; port the first engine (IF97 is self-contained — no per-fluid JSON data needed); add a `rustprop-wasm` bindings crate and the release pipeline.

## Toolchain

**Gotcha:** `~/.cargo/bin` is NOT on PATH in this environment, while a stale distro `rustc 1.75` sits at `/usr/bin/rustc` (and plain `cargo` is "not found"). Start shell work with:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
```

Installed under `~/.cargo/bin`: rustup (stable toolchain, cargo/rustc 1.97.1), the `wasm32-unknown-unknown` target, and `wasm-pack` 0.13.1.

Common commands (CI runs exactly these):

```bash
cargo build                                     # native build
cargo test --workspace --all-features           # all tests
cargo test -p rustprop <test_name>              # single test by name
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo fmt --all
cargo run -p rustprop-cli                       # example CLI
cargo build -p rustprop --features all-backends \
    --target wasm32-unknown-unknown --release   # wasm facade build
```

## Porting reference: upstream CoolProp 8 layout

Upstream core lives in `include/` + `src/` (C++17). The calculation backends under `src/Backends/` are the natural boundaries for this project's modular engines:

| Upstream backend | What it is |
|---|---|
| `Helmholtz` | HEOS — core multiparameter Helmholtz EOS (pure fluids + mixtures) |
| `Cubics` | SRK / Peng-Robinson, with Chebyshev superancillaries (new in v8) |
| `IF97` | IAPWS-IF97 industrial water/steam formulation |
| `Incompressible` | Brines, secondary working fluids |
| `PCSAFT` | PC-SAFT EOS |
| `Tabular` | TTSE / bicubic table interpolation |
| `SVDSBTL` | SVD-compressed tabular lookup (new in v8) |
| `REFPROP` | Shim to proprietary NIST REFPROP — not portable; out of scope for a pure-Rust port |

Other porting-relevant v8 facts:

- Fluid/mixture JSON data lives under upstream `dev/`. v8 replaced RapidJSON with nlohmann/json; serde is the natural Rust-side equivalent.
- v8 removed the legacy non-SI `Props`/`HAProps` API — only `PropsSI`/`HAPropsSI` semantics need porting.
- Upstream already ships an Emscripten-based JS/WASM wrapper (v8: native JS arrays). This project deliberately replaces that approach — idiomatic Rust crates instead of one emscripten blob — to get the modular, minimal binaries the wrapper cannot provide.
