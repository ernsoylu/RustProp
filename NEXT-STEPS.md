# Where this project stands, and what to do next

Read this first when picking the work back up. `PLAN.md` is the phase-by-phase
roadmap and its Decisions log is the authoritative record of *why* things are
the way they are; this file is the short version plus the open ends.

Last updated: 2026-08-07, at commit `011b31a`.

---

## Status

**All fifteen PLAN.md phases are complete.** Every checkbox is ticked, every
phase gate has passed, and CI is green. What exists:

| | |
|---|---|
| Engines ported | HEOS (pure + mixtures), IF97, cubics (SRK/PR), incompressible, PC-SAFT, tabular (TTSE/bicubic), SVDSBTL, humid air, transport, surface tension |
| Fluids | 136 HEOS (130 pure + 6 pseudo-pure), 154 predefined mixtures, 116 cubic, 126 incompressible, 180 PC-SAFT |
| Oracle records | ~33,900 committed, generated from the CoolProp 8.0.0 wheel |
| Deliverables | library crates, `rustprop-cli`, `rustprop-wasm` (wasm-bindgen), `release.yml`, CI |
| Smallest useful bundle | 124 KB (IF97) — see `WASM-SIZES.md` |

The last work done was Phase 15's seeded acceptance sweep, which found five
real defects on its first run (a panic, three wrongly-refused two-phase
outputs, and a missing echo route on two backends). All are fixed. That is
the strongest single argument for where to spend effort next: **randomized
coverage finds what hand-chosen goldens do not.**

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
| Cubic PQ flashes below **1 Pa**: upstream's equal-Gibbs secant converges, this port's gives up (4 records of 3,720) | `tests/golden/tests/acceptance.rs`, asserted to STILL reproduce | Seed, step, tolerance and iteration cap are upstream's verbatim; the difference is root conditioning inside the residual, at the extreme cold end of the cubic's own saturation range (SRK CO2 bottoms out at 91 K / 0.18 Pa — the real triple is 217 K) |
| `HAPropsSI` errors return `Result` instead of upstream's `+inf`-with-a-global | humid-air suite | A global error slot is not a thing a WASM library should have |
| PC-SAFT `WATER` PT/DT errors loudly | PC-SAFT suite, error parity asserted | Upstream computes on children whose sigma is still the −1 sentinel and returns garbage densities |
| Tabular msgpack+zlib disk cache under `~/.CoolProp/Tabular` not ported | documented in `PLAN.md` | No home directory in WASM. Cost: a LogPH table build runs ~100 s per process (40k HP flashes) — exactly the cost upstream's cache exists to avoid |
| Pseudo-pure fluids serve only PT/PQ/QT; other pairs are loud `NotImplemented` | pseudo-pure suite, verbatim error parity | Upstream routes them through legacy solvers that are dead code for the 130 superancillary fluids |
| SVDSBTL evaluator agrees to a few ulp, not bitwise (700 of 745 records bitwise, worst 1.8e-15) | `tests/golden/tests/svdsbtl.rs` | GCC compiles the reference build with `-ffp-contract=fast`. Fusing the obvious candidate makes agreement *worse*, so the contraction sits elsewhere; chasing it would match a compiler flag, not port an algorithm |

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

Four suites are `#[ignore]`d for runtime and run in the weekly CI job
(`sweep`, Mondays 04:17 UTC) or on manual dispatch:

```bash
cargo test -p rustprop-golden-tests --test heos_all_smoke -- --ignored --nocapture
cargo test --release -p rustprop-golden-tests --test mixtures mixture_sweep -- --ignored
cargo test --release -p rustprop-golden-tests --test tabular_state tabular_pairs -- --ignored
cargo test --release -p rustprop-golden-tests --test acceptance -- --ignored --nocapture
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
- wasm-pack's bundled wasm-opt 117 rejects the bulk-memory instructions rustc
  now emits. The flags are already in `crates/rustprop-wasm/Cargo.toml`
  metadata; do not remove them.
- Feature forwarding must name **both** this crate's features and the
  facade's. An `all-backends` bundle once built to 13 KB with no exports
  because only the facade's were enabled.

---

## Candidate next work

Ranked by value per unit of effort. Nothing here is required for a release.

### 1. Widen the acceptance sweep (small effort, high value)

It found five real defects at 3,720 records. Raise `N_PER` in
`gen_acceptance_sweep`, add the fluids and pairs it does not yet reach
(mixtures, predefined blends, tabular low-level states, SVDSBTL), and let the
weekly job run it. The seed makes this additive.

The sweep currently covers: HEOS 2040, SRK 360, PR 360, INCOMP 360, PCSAFT
240, HA 240, IF97 120. It does **not** yet cover mixtures or the two
low-level backends.

### 2. Audit for the defect classes the sweep exposed (small, high)

Two of the five defects were the same mistake in different places:

- **Invented guards.** Grep for errors this port raises that upstream does
  not. The pattern: a guard that "obviously" should exist (cp is undefined in
  the dome, right?) but upstream just evaluates the formula.
- **Missing echo routes.** `props_si` must return an output that *is* an
  input without updating the state. HEOS, mixtures, cubics and incompressible
  have it now; check any route added later.
- **`partial_cmp().unwrap()`.** Six more sites exist in `rustprop-heos`
  (`chebappr.rs`, `flash_hs.rs`, `flash_px.rs`). Each is a latent panic if a
  NaN ever reaches it. Check whether upstream's counterpart is a comparison
  network that passes NaN through.

### 3. The cubic sub-pascal secant (medium, low)

The one open divergence. The call parameters are verbatim; the difference is
inside `saturation_residual`'s density root selection at near-vacuum. Worth
doing only if someone actually needs cubic VLE below a pascal.

### 4. Performance (medium, unmeasured)

Nothing in this project has been profiled. The obvious candidate is the
LogPH table build at ~100 s per process. Before optimising anything, measure
— and note that `StateDerivs` caching already took LogPT from 50 s to 3.2 s
without changing a single result.

### 5. Documentation for consumers (small, medium)

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
