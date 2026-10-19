# Changelog

Format loosely follows [Keep a Changelog](https://keepachangelog.com/1.1.0/).
This project is a fidelity port, so the sections that matter most to a consumer
are not "Added" and "Fixed" but **Divergences from upstream** and **Not
ported** — a difference from CoolProp is the only kind of surprise this library
can hand you.

## [0.2.0] — 2026-10-19

### Added — every mixture keyed output upstream serves

Mixtures (`A[x]&B[y]` strings and the predefined `<Name>.mix` blends) now
answer the 33 keyed outputs that 0.1.x refused with "not ported for
mixtures": `Z`, `Phase`, `Prandtl`, `Cp0molar`/`Cp0mass`, `Helmholtzmolar`/
`Helmholtzmass`, the residual `Hmolar`/`Smolar`/`Gmolar`, the ideal-gas
`H`/`S`/`U` in molar and mass bases, `Bvirial`/`Cvirial` and their
temperature derivatives, `alphar`, `alpha0` and the five `dalpha*` strings,
`Tau`, `Delta`, `PIP`, and the isothermal-compressibility, isobaric- and
isentropic-expansion coefficients. They match the CoolProp 8.0.0 wheel to
1e-8 (1e-7 for the four Jacobian-ratio outputs) on gas, compressed-liquid and
two-phase states, pinned by the new `mixture_outputs.jsonl` (1,567 records).

With this, every output the wheel answers for a mixture at a gas, liquid or
two-phase state, the port answers identically — measured over binaries, a
ternary, the nine predefined low-GWP refrigerant blends and the classic
R4xx blends.

### Changed

- The single-state formula family (`compressibility_factor`, `cp0molar`,
  `helmholtzmolar`, the residual and ideal-gas parts, the virials,
  `tau_delta`, `alpha_all`) moved from `HelmholtzEos` inherent methods to
  default methods on `rustprop_heos::derivs::DerivEos`, which
  `HelmholtzEos`, `CubicEos` and the new `mixture_flash::MixtureAt`
  implement. Same expressions, same results; callers need `DerivEos` in
  scope.

### Added — mixture critical points

`Tcrit`, `pcrit`, `rhomolar_critical` and `rhomass_critical` now answer for
mixtures, through a port of upstream's `calc_all_critical_points`: the
criticality-contour tracer that walks `L1* = 0` and solves `L1* = M1* = 0`
wherever `M1*` changes sign, built on the third-order composition
derivatives of the GERG reducing functions and the `MixtureDerivatives` tree
(`rustprop_heos::mixture_critical`). As upstream, an answer comes only when
the tracer finds exactly one critical point; otherwise the call refuses with
"No outputs were able to be calculated".

`L1*`/`M1*` match the wheel's `criticality_contour_values()` bit for bit on
binaries and to ~1e-14 on ternaries; the critical point agrees with the wheel
for 146 of its 147 predefined blends (answer and refusal alike), pinned by
`mixture_critical.jsonl`.

### Divergences from upstream

- **HighN2.mix critical point** — the wheel answers, the port refuses. After
  its first critical point the `L1* = 0` contour folds; `L1*(theta)` has no
  root near the tracer's last angle, so the Newton step wanders, and where it
  stops depends on the derivative's last bits. Pinned with heal detection.

### Test infrastructure

- `mixture_critical.jsonl` replays in two parts: blends of up to six
  components on every run, the natural gases (7–21 components, minutes per
  blend in a debug build) in the weekly `sweep` job via
  `mixture_critical_points_heavy_match_oracle`.

- The four libm-knife-edge validity states (cubic Propane at T = 1e20 and
  1e30 K) are skipped on every platform, not only off Linux: a newer glibc on
  Linux/x86-64 lands them two ulps from the wheel too.

## [0.1.1] — 2026-08-31

### Added — a native SDK for C, C++ and desktop/cloud Rust

rustprop was built WebAssembly-first. It now ships for native targets too,
without giving up the modularity that motivated the project.

- **`rustprop-capi`: a C ABI.** `librustprop.so` / `.dylib` / `.a` /
  `rustprop.dll` with a hand-written `rustprop.h`, so C, C++, Python (ctypes),
  C#, Go, Julia, MATLAB and Fortran can call `PropsSI` and `HAPropsSI`
  directly. Twelve functions: the two calculations, a batch form, per-thread
  error reporting, and introspection.

  Unlike the JS bindings, **every function is exported by every build**. A
  call into an engine your copy was not compiled with returns
  `RUSTPROP_UNAVAILABLE` rather than failing to link, so one header works
  against any build — and `rustprop_backends()`, `rustprop_has_backend()`,
  `rustprop_fluid_count()` and `rustprop_fluid_name()` tell you what you have,
  which a prebuilt binary otherwise cannot say.

  Every entry point is safe to call concurrently from any number of threads,
  with no initialisation call and nothing to free.

- **Prebuilt SDKs for ten targets**, each with the shared and static library,
  the header, pkg-config and CMake package files, worked examples and the CLI:
  Linux x86-64 (four instruction-set baselines), Linux arm64, Linux armv7,
  macOS arm64 and x86-64, Windows x86-64 and arm64. Every target ran the full
  golden suite on its own hardware before shipping, except Linux armv7 and
  macOS x86-64, which are cross-built because no runner can execute them —
  their `BUILD-INFO.txt` says `BUILT ONLY`.

  No musl artifact. It was in the matrix and was removed after the first
  release rehearsal: musl's libm disagrees with glibc's in the `validity`
  golden suite — seven parity failures on `PR::Propane` at T = 1e30 K, by
  factors of exactly 2, 4 and 8. Every other suite passed, so the divergence
  is narrow and lives in the extreme-value tail, but this library does not
  ship numbers it has not checked. Build from source if you need musl.

- **Instruction-set variants for x86-64**: `x86-64-v2`, `-v3` and `-v4`
  alongside the portable baseline. These change only which processors the
  binary runs on. The numbers are identical — verified, not assumed: the full
  suite passes at all four baselines, and a 29,848-value sweep across every
  engine returns byte-identical results from each.

- **`USAGE.md`** — how to call rustprop from Rust, C, C++, Python, Go, Java,
  Fortran and JavaScript, which prebuilt artifact to take, containers, and
  troubleshooting. The bindings it documents are working programs that CI runs.

- **A `rust-sources` bundle** for air-gapped and vendored Rust builds,
  carrying both the published `.crate` files and the workspace source tree.
  crates.io remains the normal route. No prebuilt `.rlib` ships: Rust has no
  stable ABI, so one would link only against the exact compiler that built it
  — use the C ABI if you need a prebuilt binary.

### Changed

- The CLI is now built for all ten targets instead of three, and ships both
  standalone and inside each SDK.
- The workspace publishes **thirteen** crates rather than twelve, which raises
  the crates.io new-crate count that `release.yml`'s preflight checks against.
  See `RELEASE-CHECKLIST.md` §0.

### Fixed

- `rustprop-wasm`'s doc comment quoted `0.5548` for
  `props_si("D","T",400,"P",101325,"IF97::Water")`, which returns
  `0.5549215811909716`.
- `artifact_hardening`'s allocation check measured process-wide `VmSize` and
  so failed whenever glibc mapped a 64 MB thread arena during the measurement.
  It now records the largest single allocation request instead, which is what
  it always meant and is immune to that noise.

### Unchanged

- The WebAssembly bundles, the engines and every computed value. `rustprop-capi`
  is a leaf crate with no external dependencies that nothing else depends on,
  so no existing consumer's build is affected by its arrival.

## [0.1.0] — 2026-08-21

Tagged and released on GitHub with prebuilt wasm bundles and CLI binaries.
crates.io publication is deferred: the workspace creates twelve new crates
(thirteen as of 0.1.1, above)
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
  `NotImplemented`. Upstream routes those through legacy solvers, which are
  dead code for the 130 superancillary fluids only where the superancillary
  cascade succeeds — the ported `HS` legacy leg is reachable for Water when it
  does not (see the HS divergence in NEXT-STEPS.md).
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
