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

Phase 4 (HEOS pure fluids) is complete through 4.7 — **all seven flash pairs**: Helmholtz term
machinery (GenExp with B-recursion, NonAnalytic, GaoB; ideal container in upstream's fixed member
order), single-phase properties, classic + super-ancillaries (piecewise Chebyshev, dyadic-split
inverse-on-ln(p)), QT/PQ/PT flashes with upstream's solver strategy tree (SRK seed, Halley,
Householder4, Brent), the (D,T)/(H,P)/(P,S)/(D,P) pairs, and the (H,S) pair: runtime-built
caloric superancillaries (degree-12 L/U refit of h/s along both saturation branches),
colleague-matrix extrema (EISPACK `hqr` — logged stand-in for Eigen's RealSchur), the
two-phase Qh==Qs screen + Brent, the three-leg single-phase cascade with the (T, ln rho)
homotopy corrector, and the legacy TS-scan sad path with its (Smolar,T) inner flash — which
low-quality two-phase inputs genuinely require (upstream's endpoint Brent cannot bracket
them). Single-phase HS results carry upstream's `_Q = 10000` sentinel (oracle-confirmed).
**Six fluids golden-verified end to end** — Water, Nitrogen, CarbonDioxide, R134a, n-Propane,
Ammonia — full suite battery (terms 1e-13/1e-12, props 1e-9, ancillaries 1e-12, saturation
1e-8 policy observed ≤4e-12, PT 1e-9 with Cp/A at 1e-8, flash pairs 1e-9/1e-8, hs 1e-8 with
documented Dmolar/P scale guards) against ~6,200 committed oracle records, plus bitwise
data-fidelity walks of every fluid document. Ammonia's document carries two EOS blocks —
upstream evaluates `EOSVector[0]` only (Gao-2020), the Tillner-Roth alternate is not ported.
The melting-line caloric cascade leg is deferred with the melting line (no committed golden
needs it).

**4.8 all-fluids sweep is done**: all 130 pure fluids (every one has a superancillary) are
datagen-generated behind per-fluid features (`all-fluids` aggregate; datagen-emitted registry
`rustprop_data::fluids::all()`), bitwise fidelity-walked, and smoke-tested against the wheel
(1,922-record `#[ignore]`d suite, weekly/dispatch CI job). Six new term families landed with
representative-fluid full batteries: CP0Constant/CP0PolyT/CP0AlyLee, direct
PlanckEinsteinGeneralized, Exponential/DoubleExponential/Lemmon2005 (GenExp tau_mi channel).
Runtime `Ttriple()`/`Tmin()` = `sat_min_liquid.T` everywhere (differs from the JSON `Ttriple`
key for 27 fluids). Wasm cost: HEOS+Water 136 KB; all 130 fluids 3.31 MB (~26 KB/fluid).
The 6 pseudo-pure fluids (Air, R404A, ...) are deferred with the Maxwell fallback and pL/pV
ancillary shape.

**Phase 5 (PropsSI string API) is done**: `rustprop::props_si` with upstream's backend-prefix
parsing, registry resolution (CAS/name/alias/upper-alias, 639 queries wheel-verified), trivial
and echo routes, mass<->molar conversions, Q validation, and error conditions (196 goldens +
variant-asserted errors). Critical-parameter DISCOVERY encoded: superancillary fluids report the
NUMERICAL critical point (Tcrit_num/pmax/rhocrit_num) through every consumer, not
STATES.critical. CLI: `rustprop-cli props Dmolar T 300 P 101325 Water`; README quickstart is
real (e2e + doc-tests).

Phase 10 (mixtures) in progress: slices 10a (888 binary pairs + 28 departure fns
datagen'd behind `mixture-data`, 6 Lemmon pairs converted at datagen), 10b
(`Gerg2008Reducing`: Yr + five f_Y blocks + first/second composition derivs, both
XN conventions), and 10c (`MixtureModel`: CS sum + excess term via the new GenExp
eta1 delta-linear channel; GERG Table B5 alpha0 on STATES.critical scales with
R_mix = R_U_CODATA) are done — 864 goldens at 1e-12 across six pairs/all three
departure kinds/F=0, and 10d (PT single-phase flash: SRK-seeded lowest-Gibbs
root selection + homogeneous props, 696 goldens at 1e-9 against the wheel's real
PT update). Wheel discovery: mixture DmolarT updates re-solve density in phase
determination (delta = 1+6e-13), so Helmholtz-assembly fixtures impose
iphase_supercritical. In-dome PT returns the metastable single-phase root until
10f (stability + Michelsen split). 10e (QT/PQ VLE) is done: Wilson/preconditioner
seeds, successive_substitution (Peneloux SRK seed, hardcoded R=8.3144598),
newton_raphson_saturation (XN_DEPENDENT Jacobian, Gaussian solve for Eigen QR),
full MixtureDerivatives fugacity layer — 720 goldens at 1e-8, Q in {0,0.3,0.5,1}.
Fixed latent 10b bug: d2Yrdxidxj lacked the Gernert Table S1 XN_DEPENDENT
branches (FD regression test added). PropsSI mixture routing shipped behind the
opt-in `heos-mixtures` facade feature (~358 KB wasm): extract_fractions verbatim,
PT/QT/PQ + weighted trivials + mass basis + error parity, 284 goldens. Deviations
until 10f: sweep-based pairs (DmolarT/HmolarP/...) and mixture transport error
loudly where upstream computes them. Predefined mixtures done: 154 blends
datagen'd (`MIX_PREDEFINED`), "<Name>.mix"/uppercase registry checked before the
pure library, 175 goldens incl. Air ternary VLE and 10-component Amarillo PT.
10f part 1 done: Michelsen TPD
stability (SS+GDEM, minimize_tpd trust-region), Wilson cross-check split,
solver_rho_Tp_global (spinodal finder, omega-Halley ladder), PTflash_twophase
solve_michelsen (log-K RR, scaled Gibbs Newton, Jacobi min-eigenvalue) — full
PT_flash_mixtures glue; 192 in-dome goldens at 1e-6 (Q conditioning documented).
10f part 2 done: all ten sweep
pairs (DHSU_T/HSU_P/HSU_D with upstream's fast paths + verify gates) — 138
goldens (weekly CI job, #[ignore]d). Boost TOMS748 ported verbatim after
bisection walked into a wrong-root pocket the wheel reproduces bitwise (its
TOMS748 interpolates past). Mixture transport done (log-linear
eta / linear lambda over pure components at bulk state, 24 goldens bitwise).
**PHASE 10 COMPLETE** (phase-envelope machinery is PropsSI-dead upstream and
unported by design). Phase 11 started: 11a data done (180
fluids + 140 CAS-sorted kij behind `pcsaft-fluids`, python-verified bitwise).
Workspace serde_json now enables `float_roundtrip` — its DEFAULT parse is
best-effort and put 1-ulp errors into ~140 generated files (now corrected).
11b EOS kernels done: alphar/dadt/Z +
residual h/s/g over all five term families (shared Prep, per-kernel XA
tolerances 1e-15/1e-14, quirks documented) — 80 goldens, calorics bitwise.
11c part 1 done: fugacity coefficients
+ solver_rho_Tp (two-grid bracket scan + Brent + min-Gibbs root pick) —
TOLUENE 1 ulp / PROPANE bitwise. **PHASE 11 COMPLETE**: 11c part 2
shipped the inside-out QT/PQ flashes (kb-shadowing quirk, _HUGE=+inf clear
semantics), PT/DT phase determination, and the PropsSI PCSAFT route — 164
flash goldens at 1e-7 (TOLUENE PQ 6e-15, NaCl(aq) VLE 1.2e-9). Documented
deviation: WATER PT/DT errors loudly where upstream returns sigma-sentinel
garbage. Phase 12 started: 12a done — the
generic (T,rho) partial-derivative machinery (`rustprop_heos::derivs`:
get_dT_drho + second derivatives + the Jacobian-ratio first/second
partial_deriv) that the Tabular grid build calls at every node; 207 goldens
vs the wheel's own d(X)/d(Y)|Z strings (firsts 1e-9, seconds 1e-8).
Next: 12b table construction (LogPH/LogPT 200x200 grids + saturation table),
12c TTSE eval, 12d bicubic, 12e inversion + PropsSI routing.

Phase 6.2 (surface tension) is done: 104 curves ported and bitwise-walked, 518 goldens at
1e-12 through `props_si("I", ...)` with upstream's two-phase gating and error conditions.

**Phase 6.1 (transport) is DONE — every transport class ported**: structured families,
fully-/section-hardcoded models, Chung, rhosr-CS, and ECS (conformal-state 2-D Newton;
references Propane/R134a/Nitrogen resolved through the registry via an `EcsRef` resolver
seam; conductivity's OS critical term uses pure struct defaults — upstream never reads the
JSON `q_D`). TRANSPORT slots are per-property tri-state (Absent/Unported/Model); zero
Unported remain. Viscosity: 61 fluids/482 goldens; conductivity: 58 fluids/571 goldens, 1e-8.
Fidelity discovery: upstream v8 has NO two-phase conductivity guard (cp/cv are raw
single-phase formulas at mixture density) — two-phase states evaluate; the only two-phase
errors are conformal-solver failures (R32 et al., error parity asserted).

**Tier-2 deferrals are closed** except pseudo-pure fluids: every PropsSI-reachable input
pair is ported — (H,T)/(T,U) via generalized DHSU_T, (P,U) via HSU_P, (D,H)/(D,S)/(D,U)
via HSU_D's superancillary happy path, (D,Q) with strict-mode root enumeration; HQ/QS are
upstream string-API dead ends (generate_update_pair has no rows) and (Hmass,T)/(T,Umass)/
(Smolar,Umolar) are upstream "not yet supported" — all with exact error-message parity.
Melting lines fully ported (3 segment families, 29 fluids, PT below-Tmelt check, (P,X)
Tmelt bracket floor, HS cascade leg 4 via MeltingCaloric). Sub-triple (P,X) gas states
and multi-output-'&' parity done. ~15,800 committed oracle records.

**Pseudo-pure fluids are ported** (Air, R404A/407C/410A/507A, SES36 — 136 registry
fluids total): datagen pL/pV split (`p_s` + `p_v_split`, pure fluids alias one curve),
max_sat_T/max_sat_p state points, bitwise walker. Flashes per upstream: QT strictly
Q∈{0,1} (ancillary p + guessed PT solve), PQ with per-branch temperatures across the
glide (`HeosState::TwoPhase` carries `t_l`/`t_v`), PT via the 1.02·pL/0.98·pV ancillary
arbiter (in-band throws, as upstream). Discovery: upstream's pseudo-pure QT/PQ never
call `saturation_T_pure` — the ancillaries are used explicitly. Remaining pairs are loud
NotImplemented for pseudo-pure (upstream serves them through legacy solvers that are
dead code for the 130 superancillary fluids). 330 goldens + verbatim error parity.
~16,500 committed oracle records overall. Tier-2 deferral list: EMPTY.

**Phase 7 (cubics) is DONE**: `rustprop-cubics` ports SRK/Peng-Robinson per upstream
(T_r = 1/rho_r = 1 so tau = 1/T, delta = rho; SRK Omega literals verbatim incl. upstream's
corrupted digits; Kazakov rhomolar_critical; the smolar tau*-rescale DEFECT reproduced;
QT/PQ equal-Gibbs secant with Pitzer seeds; PT cubic-root selection with the inner PQ
branch pick; DmolarT via the six extracted Chebyshev superancillary tables with the
broken-sub-state two-phase caloric throw). `SRK::`/`PR::` PropsSI routes; 116-fluid table
behind one `cubic-fluids` feature; cubics-only wasm = 146 KB total. 2,328 cubic goldens.

**Phase 8 (incompressible) is DONE**: 126 fluids (74 pure + solutions/brines) behind
`incompressible-fluids`; Polynomial2DFrac machinery (Horner-from-top, fracIntCentral),
the five block forms, the hard-coded reference state, five input pairs, INCOMP:: with
Name/Name[x]/Name-40% parsing. 935 goldens, direct evaluations bit-identical.

**Phase 9 (humid air) is DONE**: HAPropsSI on RP-1485 virials — IAPWS-06 ice, EOS +
hardcoded virials, enhancement factor, three distinct gas constants, upstream's solver
loop shapes and quirks all reproduced; `ha_props_si` facade + CLI `ha` subcommand.
897 goldens. Deviation (logged): errors return as Result instead of +inf-with-global.
~21,000 committed oracle records overall. Next: Phase 10 (HEOS mixtures).

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
