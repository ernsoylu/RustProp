# Where this project stands, and what to do next

Read this first when picking the work back up. `PLAN.md` is the phase-by-phase
roadmap and its Decisions log is the authoritative record of *why* things are
the way they are; this file is the short version plus the open ends.

Last updated: 2026-08-17, after the output-tail + latents round (see the
2026-08-16/17 blocks in PLAN.md's Decisions log for the full record).

---

## Status

**All fifteen PLAN.md phases are complete.** Every checkbox is ticked, every
phase gate has passed, and CI is green. What exists:

| | |
|---|---|
| Engines ported | HEOS (pure + mixtures), IF97, cubics (SRK/PR), incompressible, PC-SAFT, tabular (TTSE/bicubic), SVDSBTL, humid air, transport, surface tension |
| Fluids | 136 HEOS (130 pure + 6 pseudo-pure), 154 predefined mixtures, 116 cubic, 126 incompressible, 180 PC-SAFT |
| Oracle records | 39,468 committed, generated from the CoolProp 8.0.0 wheel |
| Deliverables | library crates, `rustprop-cli`, `rustprop-wasm` (wasm-bindgen), `release.yml`, CI |
| Smallest useful bundle | 128 KB (IF97) — see `WASM-SIZES.md` |

The last work done (2026-08-16) was the post-completion audit that the
2026-08-07 handoff ranked as candidates #1 and #2 — and it vindicated the
premise twice over. The audit found and fixed TEN more invented guards, two
echo-route defects and one invented answer; the widened acceptance sweep
(3,720 → 5,485 records, plus 1,950 tabular and 1,000 SVDSBTL low-level
records) then exposed five genuinely new mixture defects, of which two were
port bugs (both fixed to bitwise agreement) and three were the WHEEL failing
against its own equilibrium (pinned in `acceptance.rs`). The lesson stands,
stronger: **randomized coverage finds what hand-chosen goldens do not** —
and a new output is only real once the sweep draws it (the Phase output
exposed a labeling divergence within minutes of existing).

The 2026-08-17 round then closed the output tail (every wheel-served output
now ported on both routes, environmental data included; sweep 6,020 records,
zero failures on first run) and both mixture latents — one pinned unported
with the wheel contradicting itself, one proved not-constructible by a
~26,000-state exact-replica scan — and fixed the hunt's bonus find, a
swallowed density-solve failure that let a degenerate Wilson split through.

**Since 2026-08-17 the work has had an external driver**: integrating rustprop
into the sibling `frees-wasm` project as its property backend (that repo's
decision D8). Two waves have landed, both merged to `main` with the full gate
green:

- **Wave 1** — the Helmholtz derivative-matrix memoization (bit-identical by
  construction; HP flash 963/839 → 380/305 µs, LogPH build 36.2 → 12.4 s),
  the rational-polynomial caloric ancillary family for pseudo-pure fluids,
  the cubic sub-pascal `psi_plus(0)` fix, and the pre-tag hygiene bundle
  (real MSRV 1.88, tag-gated release, crates.io metadata in all 12 packages).
- **Wave 2** — the last pseudo-pure input pairs, ported and goldened: `(H,P)`,
  `(P,S)`, `(P,U)`, `(D,P)` for all six pseudo-pure fluids (665-record suite,
  654 value at 1e-8 plus 11 error-parity pins, including upstream's own Air
  1-bar failure window), and the **closure of the single-phase HSU_P bisection
  stand-in** — upstream's real TOMS748 with its midpoint re-evaluation, plus
  the warm-density carry. That closure is a fidelity win first and a speed win
  second: median displacement over 1,433 (P, caloric) goldens fell from
  1.77e-10 to 2.04e-16, bitwise-exact records went 262 → 608 of 1,433, no
  fixture moved, and HP liquid Water dropped 283.6 → 80.2 µs (3.54x, quiet
  box), taking `acceptance_tabular` from 204 s to 12.6 s with it.

Wave 2 also turned up three port bugs the goldens had never reached, all
fixed to bitwise agreement: `solver_rho_tp_guessed` was missing the
phase-imposed stability retries (R410A's PQ liquid branch had been converging
to the wrong root), `ancillary::invert` used upstream's `Secant` where
`SaturationAncillaryFunction::invert` calls `ExtrapolatingSecant` (the two
differ precisely in the non-finite-residual handling the band states need),
and `update()` had no `post_update` validity gate, so some states answered
NaN where the wheel refuses.

---

## Blocked on the owner

Two things cannot be done from a coding session. Both are irreversible, which
is why they were left alone.

1. **Claim the crates.io names.** Every `rustprop*` name was verified free
   (control-checked against `serde`, which returns 200 through the same
   request). `rustprop` and `rustprop-wasm` are also free on npm.
2. **Add the `CARGO_REGISTRY_TOKEN` repository secret**, which
   `.github/workflows/release.yml` reads.

Then tagging `v0.1.0` runs the release pipeline: verify job → crates.io
publication → five wasm presets across three targets → CLI binaries for
linux-x64 / macos-arm64 / windows-x64 → GitHub release.

```bash
git tag v0.1.0 && git push origin v0.1.0
```

A crates.io version can be yanked but never replaced, so the `verify` job
re-runs fmt, clippy, the full test suite and `cargo publish --dry-run
--workspace` before anything is uploaded.

---

## Known divergences from upstream

Each of these is asserted somewhere, so it can neither widen nor silently
heal. Do not "fix" one without checking the assertion that pins it.

| Divergence | Where | Why it stands |
|---|---|---|
| THREE pinned mixture records where the port answers and the wheel's recorded value is provably not the wheel's own equilibrium (mixture HSU_P shared-state corruption ×2; shallow-TPD metastable root ×1) | `acceptance.rs` `mixture_divergences`, each pinned to the PORT's value with heal detection | Upstream's HSU_P residual mutates the shared backend (a Tmax-endpoint PT evaluation corrupts SatL/SatV and disables the two-phase split for the rest of the solve); a fresh wheel PT flash at the port's converged T reproduces the port BITWISE. The corruption is history-dependence the port's stateless flashes deliberately cannot express |
| The imposition-clear channel of the same corruption (upstream's stability feed fallback permanently clears SatL's constructor liquid imposition for the backend instance): the wheel's Nitrogen[0.97]&Water[0.03] sweep-pair inversions contradict its own forward flashes by up to 1.7 K or error outright | `hsu_p_imposition_clear_divergence_pinned` in `tests/golden/tests/mixtures.rs`, asserting the port's self-consistent inversions | Only observable through sweep-pair flashes (scalar PT builds a fresh backend per call); reproducing it would make PS/HP flashes history-dependent — same ruling as the row above. Full mechanism in the 2026-08-17 Decisions entries |
| `HAPropsSI` errors return `Result` instead of upstream's `+inf`-with-a-global | humid-air suite | A global error slot is not a thing a WASM library should have |
| PC-SAFT `WATER` PT/DT errors loudly | PC-SAFT suite, error parity asserted | Upstream computes on children whose sigma is still the −1 sentinel and returns garbage densities |
| Tabular msgpack+zlib disk cache under `~/.CoolProp/Tabular` not ported | documented in `PLAN.md` | No home directory in WASM. Cost: a LogPH table build runs ~100 s per process (40k HP flashes) — exactly the cost upstream's cache exists to avoid |
| Pseudo-pure fluids serve PT/PQ/QT plus the classic-ancillary (H,P)/(P,S)/(P,U)/(D,P) flashes (Wave-2 R6, goldened over six regions by R7); the remaining pairs (DmolarT, HS, DQ, HSU_D, ...) are loud `NotImplemented` | pseudo-pure suite (665 records: 654 value at 1e-8, 11 error-parity), verbatim error parity | Upstream routes the rest through legacy solvers that are dead code for the 130 superancillary fluids |
| R507A gas-classified caloric (P,X) states at p = 0.995·max_sat_p error loudly where the wheel converges | PLAN.md 2026-08-17 R6 entry | Both implementations fire the same gas stability retry from 1e-6 at the Tmin probe (TVanc−0.01, 0.24 K below Tcrit); the retry trajectory is chaotic through the vdW loop and only bitwise alphar arithmetic would reproduce the wheel's lucky convergence (EOS sums agree at 1e-13). A needle: 0.9925·pmax and 0.9975·pmax both agree at ≤2e-9 |
| Pseudo-pure PY-flash refusal MESSAGES are the port's own bracket diagnostic, not upstream's "unable to solve 1phase PY flash with Tmin=…, Tmax=… due to error: …" wrapper | `pseudo_pure_error_parity_matches_upstream` (refusal-vs-answer asserted for all 11 error records; oracle text carried in each record's `error` field) | The wrapper is `HSU_P_flash`'s catch around the single-phase solve and what it quotes is the INNER diagnostic. R8 closed the bisection stand-in (the inner solve is upstream's own TOMS748 now), but the no-bracket derivative path and the 2-D Newton fallback a refusal falls through to are still unported, so the text a refusal carries is still the port's bracket diagnostic. `post_update` and `solver_rho_Tp` refusals in the same suite ARE verbatim |
| Upstream's `post_update` validity gate is ported for the HEOS pure/pseudo-pure arm only, not for `mixture_update` | `props_api.rs::post_update` | No reachable NaN mixture state has been observed, and the 10f divergence pins would need re-validating; a gate without evidence is an invented guard |
| SVDSBTL evaluator agrees to a few ulp, not bitwise (700 of 745 records bitwise, worst 1.8e-15) | `tests/golden/tests/svdsbtl.rs` | GCC compiles the reference build with `-ffp-contract=fast`. Fusing the obvious candidate makes agreement *worse*, so the contraction sits elsewhere; chasing it would match a compiler flag, not port an algorithm |
| Upstream's `PXFLASH_DIRECT_EOS` cache-bypass for warm (P,X) probes is not ported (R8); the port takes upstream's own `catch (...)` fallback, the cached `update_TP_guessrho` path | `PxResid` doc comment in `flash_px.rs`; the (P,X) golden suites at their 1e-8 policy | It exists to dodge a `CachedElement` layer this port does not have (each residual owns a `DerivsMemo`), and upstream calls it "bit-equivalent within ULP" to the path we do run. Measured cost: 2.0e-16 median displacement over 1 433 (P, caloric) goldens, 608 of them bitwise |

**Wheel-vs-tag discoveries** (not divergences — the port follows the SHIPPED
wheel, which is what every golden was generated from): the 8.0.0 wheel does
not match the v8.0.0 tag source in two places even though its gitrevision is
the tag commit. (1) IF97 `set_phase`: the wheel ships pre-refactor logic —
critical-point band present, the subcritical saturation band collapses to
LIQUID, no PT two-phase throw (documented at `forward_phase` in
`if97_api.rs`). (2) HEOS DmolarT phase labels: the wheel reclassifies by the
final pressure (compressed liquid with p > pc → supercritical_liquid) where
the tag's `T_phase_determination` alone would say liquid (documented at
`dmolar_t_state`). When source and wheel disagree, this project ports the
wheel.

---

## Unported by design

Not gaps — decisions, each with its reasoning in the Decisions log.

- **The SVDSBTL builder** (Eigen BDCSVD). An SVD is unique only up to sign and
  rotation within degenerate singular subspaces, so no independent
  implementation reproduces upstream's U and V bitwise. Coefficients are
  *ingested* by `tools/rustprop-svdgen`, never recomputed. Also unported on
  that engine: the dome blend, `fast_evaluate`, REFPROP/IF97 sources, and
  `PiecewiseChebyshevCurve` (the loader refuses kind 2 loudly rather than
  guessing at the basis).
- **The SVDSBTL critical patch** is a caller seam. Upstream defaults it to
  `"auto"` and inside a calibrated bbox returns the *source backend's* value
  rather than evaluating the SVD at all — which for this port would just be
  re-testing rustprop's HEOS.
- **Mixture phase-envelope machinery** and mixture `HS_flash`. Both are
  PropsSI-dead upstream.
- **Third-order / PSI derivative machinery** for mixtures, deferred with the
  phase envelope.
- **Ammonia's Tillner-Roth alternate EOS.** The document carries two EOS
  blocks; upstream evaluates `EOSVector[0]` (Gao-2020) only.
- **REFPROP backend.** Proprietary shim, out of scope for a pure-Rust port.

---

## How to work on this

### Toolchain

`~/.cargo/bin` is not on `PATH` and a stale distro `rustc 1.75` sits at
`/usr/bin/rustc`. Start every shell with:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
```

### The gate

Nothing gets committed until all of these are zero. CI runs the same set.

```bash
cargo fmt --all
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features --release
cargo run -q -p rustprop-datagen && git diff --exit-code   # datagen determinism
cargo build -p rustprop --features all-backends --target wasm32-unknown-unknown --release
cargo publish --dry-run --workspace
node tests/wasm-smoke/smoke.mjs        # after a wasm-pack --target nodejs build
```

### Heavy suites

Six suites are `#[ignore]`d for runtime and run in the weekly CI job
(`sweep`, Mondays 04:17 UTC) or on manual dispatch:

```bash
cargo test -p rustprop-golden-tests --test heos_all_smoke -- --ignored --nocapture
cargo test --release -p rustprop-golden-tests --test mixtures mixture_sweep -- --ignored
cargo test --release -p rustprop-golden-tests --test tabular_state tabular_pairs -- --ignored
cargo test --release -p rustprop-golden-tests --test acceptance -- --ignored --nocapture
cargo test --release -p rustprop-golden-tests --test acceptance_tabular -- --ignored
cargo test --release -p rustprop-golden-tests --test acceptance_svdsbtl -- --ignored
```

### Regenerating fixtures

The oracle is the pinned CoolProp 8.0.0 wheel in `tools/golden-gen/.venv`.
Fixtures are committed; regenerate only through the generator.

```bash
cd tools/golden-gen && ./.venv/bin/python gen_fixtures.py     # all
./.venv/bin/python -c "import gen_fixtures as g; g.write_jsonl('x.jsonl', g.gen_x())"  # one
```

The acceptance sweep is seeded (`random.Random(20260807)`), so raising
`N_PER` widens coverage without invalidating existing records — the first N
draws are unchanged.

### Pitfalls learned the hard way

- **Python edit scripts write only at the end.** A mid-script `assert` failure
  silently discards every earlier replacement. This cost three debugging
  sessions. Use the `Edit` tool for targeted changes.
- `cargo fmt` reformats code, so exact-string replacements stop matching
  after it runs. Format last.
- `cargo publish --dry-run -p <crate>` **cannot** work before the first
  publish — cargo resolves each crate's rustprop deps against the crates.io
  index, and `--no-verify` does not help because packaging resolves too.
  `--workspace` works: it verifies against the local crates in dependency
  order.
- The `--workspace` dry-run's overlay registry is IMMUTABLE at name+version,
  like any registry: the first dry-run caches each `rustprop-*-0.1.0`
  tarball into `~/.cargo/registry/{cache,src}` and its compiled rlib into
  `target/`, and no later dry-run re-reads changed sources. A cross-crate
  API added AFTER the first dry-run fails verification with a baffling
  E0432 while every local file and freshly packaged tarball is correct.
  Purge all three layers before re-running after cross-crate API changes:
  `rm -rf target/package ~/.cargo/registry/{src,cache}/*/rustprop-*` plus
  `cargo clean -p` on the touched crates. (Three stacked cache layers of
  debugging, 2026-08-17.)
- wasm-pack's bundled wasm-opt 117 rejects the bulk-memory instructions rustc
  now emits. The flags are already in `crates/rustprop-wasm/Cargo.toml`
  metadata; do not remove them.
- Feature forwarding must name **both** this crate's features and the
  facade's. An `all-backends` bundle once built to 13 KB with no exports
  because only the facade's were enabled.

---

## Candidate next work

Ranked by value per unit of effort. Nothing here is required for a release.

*(2026-08-16: the previous #1 — widen the sweep — and #2 — audit the defect
classes — are DONE; see PLAN.md's Decisions log. The sweep now covers
mixtures, blends, wide outputs, IF97 pairs, pseudo-pure transport, and both
low-level backends: HEOS 2580 + 865 mix, SRK/PR 450 each, INCOMP 360,
PCSAFT 280, HA 240, IF97 260, plus tabular 1,950 and SVDSBTL 1,000 in their
own fixtures. The partial_cmp sites were proven NaN-unreachable and closed.)*

*(2026-08-17: the previous #1 — the output tail — and #2 — the two mixture
latents — are DONE; see PLAN.md's Decisions log. Every output the wheel
serves is ported on both routes; the imposition-clear latent is pinned
unported by `hsu_p_imposition_clear_divergence_pinned`, the Brent-throw
latent proved not-constructible, and the hunt's bonus find — swallowed
density-solve failures in `successive_substitution_guessrho` — is fixed.)*

*(2026-08-17, later: the previous #1 — the cubic sub-pascal secant — is
FIXED, not excused: two association deviations in `psi_plus(0)` put ulp-level
noise into the equal-Gibbs residual exactly where near-vacuum cancellation
left no headroom. All 58 gate-band flashes now converge, the nine formerly
erroring records match the wheel to ≤2 ulp, and the acceptance allowance is
deleted. See the Decisions log.)*

### 1. Close the `Ok(non-finite)` validity gap (small-to-medium, evidence in hand)

**rustprop returns `Ok(NaN)` in places where the wheel raises.** A ~1.8-million
combination scan over (output, pair, value, fluid), run during the Wave-2
integration, found **1,754 such calls**. Witness, re-confirmed against the
wheel while writing this:

```
PropsSI("L","T",1e30,"P",101325,"Water")
  rustprop -> Ok(NaN)
  wheel    -> raises "calc_alpha0_deriv_nocache returned invalid number
              with inputs nTau: 2, nDelta: 0, tau: 6.47096e-28"
```

It looks like a missing `ValidNumber` check on the ideal-gas derivative — the
same *class* of defect as the invented guards the 2026-08-16 audit removed,
but inverted: a guard upstream HAS that the port does not. Note this is
exactly what the sweep is for, and the sweep did not draw it — the scan that
found it varied *inputs* far outside the fluid's range, which the seeded
acceptance sweep deliberately does not do. Worth deciding whether that
input-abuse dimension belongs in the sweep permanently.

Not urgent for the frees integration (that repo's backend rejects non-finite
values on both sides of every call, so this cannot reach its engine), but it
is a real divergence and currently pinned nowhere on this side.

### 2. Decide the `post_update` refusal text (small)

R7's `post_update` port turned `PropsSI("Hmass","T",0,"Smass",101325,"Water")`
from `Ok(NaN)` into a refusal, matching the wheel — but the message differs:
rustprop says *"rhomolar is not a valid number"*, the wheel *"rhomolar is less
than zero"*, because rustprop's density is NaN there where upstream's is
negative. The `< 0` arm IS carried, in upstream's order, so this is a
value-level difference rather than a missing branch. Either pin it as a
divergence row or chase the density difference; do not leave it unrecorded.

### 3. Performance (measured; the obvious targets are gone)

This is no longer unprofiled ground — `tools/perf-bench` and `PERF.md` exist,
and the two biggest costs have been paid: the derivative-matrix memoization
(Wave 1) and the TOMS748 + warm-carry closure (Wave 2) together took the HP
flash from ~963 µs to ~80 µs and the LogPH-backed suites from ~200 s to
~13 s. What remains is smaller and should be measured before it is touched:
the mixture-side stability pairs in `mixture_vle.rs` still evaluate
`dpdrho`/`d2pdrho2` separately and would take the same bit-identical memo
pattern.

### 4. Documentation for consumers (small, medium)

The README quickstart is real and doc-tested. What does not exist: per-engine
guidance on *which* backend to choose, and a worked WASM example beyond the
size table.

---

## Orientation for a new session

- `PLAN.md` — phases, checkboxes, and the Decisions log. The log is
  append-only and is where every non-obvious choice is justified. Read the
  entries for whatever you are about to touch.
- `CLAUDE.md` — the working status summary, kept current at phase gates.
- `WASM-SIZES.md` — measured bundle sizes, regenerate with
  `tools/wasm-size-table.sh`.
- Upstream lives at `~/homecloud/dev/CoolProp`, pinned at tag `v8.0.0`. It is
  the correctness oracle; when in doubt, read its source rather than
  reasoning about what it probably does. Several findings in the Decisions
  log exist because the code disagreed with the documentation.
