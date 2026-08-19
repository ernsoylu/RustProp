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
5. CI/CD pipelines covering build, lint, the full test suite, and release packaging — GitHub Actions. Both exist: `.github/workflows/ci.yml` (every push, plus a weekly `sweep` job for the `#[ignore]`d heavy suites) and `.github/workflows/release.yml` (fires on a `vX.Y.Z` tag).

## Implementation rules

- **`PLAN.md` is no longer a roadmap** — all fifteen phases are done. It is now the append-only **Decisions log**: every non-obvious choice gets an entry, in the same commit as the work, and no existing entry is ever rewritten.
- Follow `.claude/skills/karpathy-guidelines/SKILL.md` for all code: surface assumptions before coding; write the minimum code that solves the step (nothing speculative); make surgical diffs; frame every task as a verifiable goal and loop until its check passes.
- Nothing is committed until the full gate is green. The commands are in `NEXT-STEPS.md`; `cargo fmt --all` runs last, because formatting invalidates exact-string edits.

## Picking the work back up

**Read `NEXT-STEPS.md` first.** It carries the current status, the two items
blocked on the owner (crates.io names + token), the full list of known
divergences and by-design exclusions, the gate commands, the pitfalls learned,
and ranked candidates for the next piece of work. `PLAN.md`'s Decisions log
remains the authoritative record of *why* each choice was made.

## Status (2026-08-18)

**All fifteen PLAN.md phases are complete** — every checkbox ticked, every gate
passed. rustprop is a feature-complete port of CoolProp 8.0.0's `PropsSI` /
`HAPropsSI` surface. Work no longer advances phases; it lands in **waves**
driven by an external consumer (see below).

| | |
|---|---|
| Engines | HEOS (pure, pseudo-pure, mixtures), IF97, cubics (SRK/PR), incompressible, PC-SAFT, tabular (TTSE/bicubic), SVDSBTL, humid air, transport, surface tension |
| Fluids | 136 HEOS (130 pure + 6 pseudo-pure), 154 predefined mixtures, 116 cubic, 126 incompressible, 180 PC-SAFT |
| Oracle records | 41,629 in 123 committed fixtures, read by 35 suites (`cat tests/golden/fixtures/*.jsonl \| wc -l`) |
| Deliverables | engine crates, `rustprop-cli`, `rustprop-wasm`, `release.yml`, CI |
| Bundle sizes | measured in `WASM-SIZES.md` — 128 KB (IF97) to 4.2 MB (all-backends) |

**Blocked on the owner, and only on the owner**: claim the crates.io names, add
the `CARGO_REGISTRY_TOKEN` secret, then tag `v0.1.0`. Publication is
irreversible; nothing else stands between the tree and a release.

Nothing from `NEXT-STEPS.md` or `PLAN.md` is restated below, on purpose: this
file goes stale the moment it duplicates them.

### Architecture rules the code encodes

- **Types/contents split**: fluid-data *types* live in `rustprop-core`,
  generated *contents* in `rustprop-data` (one Cargo feature per fluid,
  `default = []`). Engines depend only on core, so an app links data solely for
  the fluids it opts into.
- The facade crate `rustprop` puts every engine behind a Cargo feature with
  `default = []`; `all-backends` turns everything on (CLI and CI use it).
- `rustprop-data` contents come only from `tools/rustprop-datagen` codegen,
  never hand edits. JSON parsing must stay out of shipped binaries — serde is
  confined to dev tooling and the unpublished test harness. `wasm-bindgen` (in
  `rustprop-wasm`) is the only external dependency in any shipped crate.
- Workspace lints deny `unsafe_code`; the release profile uses fat LTO,
  `panic = "abort"`, symbol stripping.
- Golden fixtures are generated ONLY by `tools/golden-gen/gen_fixtures.py` and
  committed. Never hand-edit one; regenerate through its generator.

### Fidelity rules, learned the expensive way

These override any instinct about what the code "should" do:

- **A guard upstream does not have is a DEFECT**, however defensive it looks.
  Roughly twenty invented guards have been found and removed this way; every
  hunt for them has found more.
- Reproduce upstream's bugs, quirks and failure windows as **error parity** —
  same refusal, same state, verbatim message where the port can reach it.
- **When the shipped 8.0.0 wheel disagrees with the v8.0.0 tag source, port the
  WHEEL** (two such discoveries so far: IF97 `set_phase`, HEOS DmolarT phase
  labels). Every golden came from the wheel.
- Where upstream's answer is an artifact of its own mutable backend state
  (stale caches, cross-call corruption), the port does not reproduce it — a
  stateless port cannot express history-dependence, and would be serving
  knowingly worse numbers if it tried. Such cases are pinned as divergences
  with evidence that the wheel contradicts itself; see the NEXT-STEPS table.
- **Randomized and unread coverage finds what hand-chosen goldens do not.** The
  seeded acceptance sweep has produced a real defect nearly every time it
  widened, and the six fixture batteries nobody had wired up produced another
  the day they were first run. Prefer widening coverage over adding more
  hand-picked states.

### Waves since phase completion

Driven by integrating rustprop into the sibling `frees-wasm` project as its
property backend (that repo's decision D8). Full accounts live in PLAN.md's
Decisions log.

- **Wave 1** — Helmholtz derivative-matrix memoization (bit-identical by
  construction; HP flash 963 → 380 µs, LogPH build 36.2 → 12.4 s); the
  rational-polynomial caloric ancillary family; the cubic sub-pascal
  `psi_plus(0)` fix; pre-tag hygiene (MSRV 1.88, tag-gated release, crates.io
  metadata).
- **Wave 2** — the last pseudo-pure input pairs `(H,P)`/`(P,S)`/`(P,U)`/`(D,P)`
  with a 665-record suite, and the closure of the single-phase HSU_P bisection
  stand-in by upstream's real TOMS748 plus a warm-density carry (median
  displacement over 1,433 goldens 1.77e-10 → 2.04e-16; HP liquid Water 283.6 →
  80.2 µs). Three port bugs fixed to bitwise along the way.
- **Wave 3** — the `Ok(non-finite)` validity gap closed by porting upstream's
  two real gates (`calc_alpha0_deriv_nocache`'s `ValidNumber` throw and the
  binding-level `_raise_if_invalid`), pinned by a 1,626-record `validity` suite;
  the acceptance sweep widened to 6,380 with pseudo-pure caloric draws; and the
  six Phase-4.8 fixture batteries that no test had ever read wired in, which
  exposed upstream's PT stale-cache quirk (`heos_pt.rs`).

## Toolchain

**Gotcha:** `~/.cargo/bin` is NOT on PATH in this environment, while a stale distro `rustc 1.75` sits at `/usr/bin/rustc` (and plain `cargo` is "not found"). Start shell work with:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
```

Installed under `~/.cargo/bin`: rustup (stable toolchain, cargo/rustc 1.97.1), the `wasm32-unknown-unknown` target, and `wasm-pack` 0.13.1.

Common commands (CI runs exactly these):

```bash
cargo build                                     # native build
cargo test --workspace --all-features           # all tests, DEBUG (ci.yml)
cargo test --workspace --all-features --release # all tests, RELEASE (release.yml verify)
cargo test -p rustprop <test_name>              # single test by name
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo fmt --all
cargo run -p rustprop-cli                       # example CLI
cargo build -p rustprop --features all-backends \
    --target wasm32-unknown-unknown --release   # wasm facade build
```

**Gate the test suite in BOTH profiles.** `ci.yml` runs debug, `release.yml`'s
verify job runs `--release`, so a one-profile gate can pass here and fail on
the tag. It has happened: `heos_pt` asserted bitwise on a PT density that the
two profiles round to doubles 42 ulp apart — both reproducing the requested
pressure — and the assertion had to be relaxed to 1e-13 (`d5a7331`). Never
assert bitwise on anything an iterative solver produced.

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
