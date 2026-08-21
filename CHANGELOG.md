# Changelog

Format loosely follows [Keep a Changelog](https://keepachangelog.com/1.1.0/).
This project is a fidelity port, so the sections that matter most to a consumer
are not "Added" and "Fixed" but **Divergences from upstream** and **Not
ported** — a difference from CoolProp is the only kind of surprise this library
can hand you.

## [0.1.0] — 2026-08-21

Tagged and released on GitHub with prebuilt wasm bundles and CLI binaries.
crates.io publication is deferred: the workspace creates twelve new crates
against crates.io's new-crate rate limit (burst 5, one per 10 minutes), so
`release.yml`'s preflight refuses to publish until a limit override is in
place — see `RELEASE-CHECKLIST.md` §0. Until then, depend on the tag:
`rustprop = { git = "https://github.com/ernsoylu/RustProp", tag = "v0.1.0" }`.

First release. A from-scratch pure-Rust port of
[CoolProp 8.0.0](https://github.com/CoolProp/CoolProp)'s `PropsSI` /
`HAPropsSI` surface, built so a WebAssembly bundle carries only the engines and
the fluids it actually uses.

### Added

**Engines**, each an independently selectable crate behind a Cargo feature
(`default = []` everywhere):

| Feature | Engine | Fluids |
|---|---|---|
| `heos` | Multiparameter Helmholtz (upstream `Backends/Helmholtz`) | 130 pure + 6 pseudo-pure, opt-in per fluid |
| `heos-mixtures` | HEOS mixtures: binary pairs, departure functions, VLE, flashes | 154 predefined blends |
| `if97` | IAPWS-IF97 industrial water/steam | self-contained |
| `cubics` | SRK / Peng-Robinson with v8's Chebyshev superancillaries | 116 |
| `incompressible` | Brines and secondary working fluids | 126 |
| `pcsaft` | PC-SAFT | 180 |
| `humid-air` | `HAPropsSI` psychrometrics | — |
| `tabular` | TTSE / bicubic tables (low-level API, as upstream) | — |
| `svdsbtl` | SVD-compressed tabular lookup, new in v8 (low-level API) | — |

Transport properties (viscosity, conductivity), surface tension, melting and
saturation lines, and partial derivatives come with the engines that own them.

**Crates** (twelve published): `rustprop` (facade), `rustprop-core`,
`rustprop-data`, `rustprop-heos`, `rustprop-if97`, `rustprop-cubics`,
`rustprop-incompressible`, `rustprop-pcsaft`, `rustprop-tabular`,
`rustprop-svdsbtl`, `rustprop-humid-air`, `rustprop-wasm`.

**Modularity.** Fluid *data* is opt-in per fluid, not just per engine, because
data dominates bundle size. Measured `wasm-pack --target web --release` output
(full table in `WASM-SIZES.md`): IF97 alone **127.8 KB**, HEOS + Water
**339.3 KB**, humid air **252.2 KB**, everything at once **4.2 MB**.

**Deliverables.** `rustprop-cli` (an example CLI over stdout), `rustprop-wasm`
(wasm-bindgen bindings with a `Float64Array` batch path), prebuilt wasm bundles
in five presets × three wasm-pack targets, and CLI binaries for linux-x64 /
macos-arm64 / windows-x64 — the last two attached to the GitHub release.

**Engineering.** Rust 1.88 MSRV, edition 2024, `unsafe_code = "deny"`
workspace-wide, no C/C++ FFI, and `wasm-bindgen` as the only external
dependency in any shipped crate. `serde` is confined to dev tooling and the
test harness — no JSON parser ships in a binary.

### Verification

41,629 oracle records in 123 committed fixtures, read by 35 test suites. Every
record is the answer of one specific CoolProp 8.0.0 binary — see
`tools/golden-gen/ORACLE.md`, which pins it by sha256 and archives it. Most
suites assert relative agreement at 1e-12 or tighter and a large fraction match
bitwise; the tolerances that are looser are looser for a reason and each one is
justified where it is written.

### Divergences from upstream

Every one of these is asserted by a test, so it can neither widen nor silently
heal. The full table with mechanisms lives in `NEXT-STEPS.md`; this section is
the subset a consumer can actually reach.

**API-shaped — you will notice these immediately:**

- **`HAPropsSI` errors return `Err`**, not upstream's `+inf` with a global
  error slot. Same for the rest of the API: this port returns `Result`. A
  global error slot is not a thing a WASM library should have.
- **`rustprop-cli` is not published** to crates.io — `cargo install
  rustprop-cli` does not work. It is the example app; take the prebuilt binary
  from the GitHub release, or build it from the repository.
- **Derivative output strings are not parsed.** `PropsSI("d(Hmolar)/d(T)|P",
  ...)` and its siblings are rejected where upstream answers. The machinery
  exists (`rustprop_heos::derivs`, 207 goldens); only the string parser is
  missing. This is the one item in this list that is a plain missing feature
  rather than a decision.
- **Pseudo-pure fluids serve a subset of input pairs.** `PT`, `PQ`, `QT` and
  the four classic-ancillary caloric pairs `(H,P)` / `(P,S)` / `(P,U)` /
  `(D,P)` work; `DmolarT`, `HS`, `DQ`, `HSU_D` and the rest raise a loud
  `NotImplemented`. Upstream routes those through legacy solvers that are dead
  code for the 130 superancillary fluids.
- **`(Dmolar, P)` below the triple-point pressure** raises `NotImplemented`
  where upstream answers (148 states in a systematic scan). Not reachable at
  physical inputs.

**Numerical — only visible if you compare against CoolProp digit by digit:**

- **HEOS `PT` flash: `h`/`s`/`cp`/`w` differ by up to 5.0e-8 near the critical
  point.** Density matches bitwise. Upstream serves these properties off its
  density solver's *last trial iterate* rather than the root it returns, so a
  CoolProp `PT` state is internally inconsistent; this port evaluates at the
  root. The port's answer is the self-consistent one.
- **Three mixture states** where the port answers and the recorded upstream
  value is provably not upstream's own equilibrium (a shared-backend corruption
  in its `HSU_P` residual). A fresh CoolProp `PT` flash at the port's converged
  temperature reproduces the port bitwise.
- **SVDSBTL evaluator agrees to a few ulp, not bitwise** (700 of 745 records
  bitwise, worst 1.8e-15). Upstream's reference build uses
  `-ffp-contract=fast`; matching it would mean matching a compiler flag.
- **One state in the whole 136-fluid registry** where upstream answers and this
  port refuses: MethylLinoleate `(P, H/S/U)` at `p = 1.001 · p_triple`. Found
  by a 157,374-state scan that agrees everywhere else. Upstream's own validity
  gate refuses the same state 614 ulp lower.
- **R507A** gas-classified caloric `(P,X)` states at exactly
  `p = 0.995 · p_sat,max` raise where upstream converges; 0.9925 and 0.9975 of
  the same maximum agree to ≤2e-9. A chaotic retry trajectory, not a
  systematic error.
- **Refusal messages** are verbatim upstream text wherever the port reaches
  it, but a few are the port's own diagnostic — pseudo-pure `PY`-flash
  refusals, and one `post_update` message that says "rhomolar is not a valid
  number" where upstream says "less than zero". Refusal-vs-answer always
  agrees; only the text differs.

**Behavioural:**

- **No disk cache for tabular tables.** Upstream memoizes TTSE/bicubic tables
  under `~/.CoolProp/Tabular`; WASM has no home directory. Cost: a LogPH table
  build runs ~100 s per process. This is the one gap a consumer actually feels.
- **PC-SAFT `WATER` `PT`/`DT` errors loudly**, matching upstream's failure but
  not its output — upstream computes on children whose `sigma` is still a −1
  sentinel and returns garbage densities.

**Where upstream's shipped wheel disagrees with its own v8.0.0 tag source, this
port follows the wheel** — that is what every golden was generated from. Two
such places are known: IF97 `set_phase` (the wheel ships pre-refactor logic)
and HEOS `DmolarT` phase labels (the wheel reclassifies by final pressure).

### Not ported, by design

- **The REFPROP backend** — a shim to proprietary NIST code, out of scope for a
  pure-Rust port.
- **The SVDSBTL *builder*** (Eigen BDCSVD). An SVD is unique only up to sign
  and rotation within degenerate singular subspaces, so no independent
  implementation reproduces upstream's factors bitwise. Coefficients are
  *ingested* from upstream artifacts, never recomputed. Also unported on that
  engine: the dome blend, `fast_evaluate`, the REFPROP/IF97 sources, and
  `PiecewiseChebyshevCurve` — the loader refuses that kind loudly rather than
  guessing at the basis.
- **The SVDSBTL critical patch**, which upstream resolves by returning the
  *source backend's* value inside a calibrated bounding box — for this port
  that would just re-test its own HEOS.
- **Mixture phase-envelope machinery and mixture `HS_flash`.** Both are
  unreachable from `PropsSI` upstream.
- **Third-order / PSI derivative machinery for mixtures**, deferred with the
  phase envelope.
- **Ammonia's Tillner-Roth alternate EOS.** The fluid document carries two EOS
  blocks; upstream evaluates only `EOSVector[0]` (Gao 2020), so this port does
  too.
- **Upstream's Emscripten JS/WASM wrapper**, deliberately replaced by idiomatic
  Rust crates — the whole point of the modular bundle sizes above.

Two upstream behaviours are *not* reproduced on principle: answers that are
artifacts of upstream's own mutable backend state (stale caches, cross-call
corruption). A stateless port cannot express history-dependence, and
reproducing it would mean serving knowingly worse numbers. Each such case is
pinned as a divergence with evidence that upstream contradicts itself.

[0.1.0]: https://github.com/ernsoylu/RustProp/releases/tag/v0.1.0
