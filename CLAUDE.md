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

## Implementation rules

- The roadmap is **`PLAN.md`** — phased, checkbox-tracked, every step carrying a `→ verify:` clause. Work the phases in order; tick boxes and append to its Decisions log in the same commit as the work; update the Status section below at phase gates.
- Follow `.claude/skills/karpathy-guidelines/SKILL.md` for all code: surface assumptions before coding; write the minimum code that solves the step (nothing speculative); make surgical diffs; frame every task as a verifiable goal and loop until its check passes.

## Status (as of 2026-08)

PLAN.md **Phases 0–3 are complete**: verification infrastructure; core parameter system + error types; and the **IF97 steam engine, fully ported and golden-verified** — 356 oracle records match at rtol 1e-11 (most at 1e-12), all IAPWS published check tables pass, steam properties work end-to-end from the CLI (`cargo run -p rustprop-cli -- props T P 101325 Q 0 IF97::Water`) and compile to ~113 KB of wasm (`if97` feature). Live now: the pinned upstream checkout (sibling `~/homecloud/dev/CoolProp`, tag `v8.0.0`), the golden-fixture oracle (CoolProp 8.0.0 wheel in `tools/golden-gen/.venv`), the comparison harness (`tests/golden`), and a CI wasm-size report (70-byte baseline). Shipped crates have zero external dependencies — serde/serde_json are confined to the unpublished test harness. The crate map is in README.md. Architecture rules the scaffold encodes:

- **Types/contents split**: fluid-data *types* live in `rustprop-core`, generated data *contents* in `rustprop-data` (one Cargo feature per fluid, `default = []`). Engines depend only on core — apps link data solely for the fluids they opt into.
- The facade crate `rustprop` puts every engine behind a Cargo feature with `default = []`; `all-backends` turns everything on (used by the CLI and CI).
- `rustprop-data` contents come only from `tools/rustprop-datagen` codegen, never hand edits; JSON parsing must stay out of shipped binaries.
- Workspace lints deny `unsafe_code`; release profile uses fat LTO, `panic = "abort"`, symbol stripping.

Phase 3 added the data pipeline: fluid JSON dumped verbatim from the oracle wheel into `data/coolprop-json/` (pinned, attributed), `rustprop-datagen` emitting feature-gated modules into `rustprop-data` (Water first, bitwise fidelity-tested, CI regeneration guard).

Phase 4 (HEOS pure fluids) is complete through 4.7 except the (H,S) pair: Helmholtz term
machinery (GenExp with B-recursion, NonAnalytic, GaoB; ideal container in upstream's fixed member
order), single-phase properties, classic + super-ancillaries (piecewise Chebyshev, dyadic-split
inverse-on-ln(p)), QT/PQ/PT flashes with upstream's solver strategy tree (SRK seed, Halley,
Householder4, Brent), and the (D,T)/(H,P)/(P,S)/(D,P) pairs. **Six fluids ported and
golden-verified end to end** — Water, Nitrogen, CarbonDioxide, R134a, n-Propane, Ammonia — each
running the full suite battery (terms 1e-13/1e-12, props 1e-9, ancillaries 1e-12, saturation 1e-8
policy with observed ≤4e-12 for the new fluids, PT 1e-9 with Cp/A at the solver-dependent 1e-8
tier, flash pairs 1e-9/1e-8) against ~5,900 committed oracle records, plus bitwise data-fidelity
walks of every fluid document. Ammonia's document carries two EOS blocks — upstream evaluates
`EOSVector[0]` only (Gao-2020), the Tillner-Roth alternate is not ported.

Next work: `PLAN.md` 4.6 (H,S) — the superancillary hs cascade (own session; structure scoped in
the Decisions log) — then 4.8 all-fluids sweep.

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
