# Where this project stands, and what to do next

Read this first when picking the work back up. `PLAN.md` is the phase-by-phase
roadmap and its Decisions log is the authoritative record of *why* things are
the way they are; this file is the short version plus the open ends.

Last updated: 2026-08-19, after R10 (release engineering) — off-Linux platform
tests, the cargo-deny gate, the oracle's archival, `RELEASE-CHECKLIST.md`,
`CHANGELOG.md`, and a crates.io rate limit that would have half-published the
project. Before it, the Wave-4b integration merged R12 (registry-wide (P,X)
refusal parity, the superancillary inverse's TOMS748, the assertion audit) and
re-measured its claims. See the dated blocks in PLAN.md's Decisions log for the
full record of each round.

**Before doing anything release-shaped, read `RELEASE-CHECKLIST.md` §0.**

---

## Status

**All fifteen PLAN.md phases are complete.** Every checkbox is ticked, every
phase gate has passed, and CI is green. What exists:

| | |
|---|---|
| Engines ported | HEOS (pure + mixtures), IF97, cubics (SRK/PR), incompressible, PC-SAFT, tabular (TTSE/bicubic), SVDSBTL, humid air, transport, surface tension |
| Fluids | 136 HEOS (130 pure + 6 pseudo-pure), 154 predefined mixtures, 116 cubic, 126 incompressible, 180 PC-SAFT |
| Oracle records | 41,629 in 123 committed fixtures, read by 35 suites — every fixture is now consumed by a test (Wave-3 R9) |
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
decision D8). Three waves have landed, each with the full gate green:

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
- **Wave 3** — the `Ok(non-finite)` validity gap, closed by porting upstream's
  two real gates (`calc_alpha0_deriv_nocache`'s `ValidNumber` throw and the
  scalar bindings' `_raise_if_invalid`) at their own sites and pinning them
  with `tests/golden/tests/validity.rs` (1,626 records); a 58,926-combination
  abuse grid took "port non-finite / wheel raises" from 2,898 to 0 with "port
  raises / wheel answers" set-identical at 314. Then R9: the seeded acceptance
  sweep widened 6,020 → 6,380 with pseudo-pure `(H,P)`/`(P,S)`/`(P,U)`/`(D,P)`
  draws (clean on first run — the first widening that found nothing, which is
  the honest verdict on Wave 2's pseudo-pure work), and the **six Phase-4.8
  fixture batteries that no test had ever read** were wired into
  `heos_fluids()`. Those 4,662 records exposed an upstream quirk the six
  original fluids had only grazed: `PT_flash` serves properties off the density
  solver's LAST TRIAL (see the divergence table).
- **Wave 4** — R12 measured `(P, caloric)` refusal parity across the whole
  registry instead of assuming it: 136 fluids × 157,374 states classified
  answer-vs-refuse on both sides. The wheel refuses at 297 (Air 177, SES36 89,
  R407C 31) and the port at the same 297 but three. Air's low-quality `(P,h)`
  band — Wave-3 F8's open question — IS upstream parity: its upper edge is
  bitwise identical on both sides at nine pressures and its lower edge agrees
  to ≤7.1e-9 of the dome width, because upstream's unported Halley / 2-D Newton
  continuations fail everywhere in the band too. The three exceptions are one
  MethylLinoleate knife-edge state, now a pinned divergence (see the table).
  R12 also closed the superancillary inverse's bisection stand-in (upstream
  calls boost TOMS748 there) and audited every exact/≤1e-12 assertion in the
  suite, relaxing the one that pinned a tabular density bitwise.
- **Wave-4b integration** — R12 merged, and its three headline claims
  re-measured from scratch rather than taken on report. All 175 fixture records
  replay from the wheel exactly (128 answers bitwise, 47 refusal messages
  verbatim); both Air band edges were re-bisected on both sides at nine
  pressures and reproduce R12's table digit for digit (upper edge bitwise
  everywhere, lower edge worst 2.24e-6 relative / 7.08e-9 absolute in quality);
  and the MethylLinoleate knife-edge was re-derived, which sharpened it — the
  wheel's own refusal edge is **614 ulp** below the Q = 0 enthalpy it serves,
  the port's is 9,435 ulp above it, so the disagreement window is 10,048 ulp of
  h ≈ 2.9e-12 of quality, and the wheel refuses with the port's exact message
  just outside it. The size table and the oracle-record counters had drifted
  and were re-measured (see the Decisions log).

---

## Release readiness

**Restated 2026-08-19 on the R10 merge (`5c93a15`), superseding the Wave-4b
version.** Nothing below is carried forward: every claim was re-derived on the
merged tree, and each one that R10 reported was re-measured independently
rather than trusted. What changed since the Wave-4b statement is listed at the
end of this section.

**What is green, measured on this tree:** `fmt --check`; `clippy --workspace
--all-targets --all-features -D warnings`; `cargo test --workspace
--all-features` in **debug** and in **`--release`** — 134 passed, 0 failed,
7 ignored, across 71 test binaries, identical in both profiles; the datagen
determinism check (regenerate → empty diff, no untracked files); the wasm32
`all-backends` release build; all six `#[ignore]`d weekly suites run exactly as
`ci.yml`'s `sweep` job runs them (130 fluids; mixture sweep; tabular_pairs
1,632 records, worst 5.9e-9; acceptance 6,380 records over 8 engines with the 3
pinned mixture divergences; acceptance_tabular 1,950 records; acceptance_svdsbtl
995/1000 bitwise, worst 8.2e-16); `cargo publish --dry-run --workspace` from
genuinely cold caches (12 packaged, 12 verified, 12 uploads aborted, no errors);
the MSRV job's `cargo +1.88 check` on both shipped facades; `wasm-pack --target
nodejs` + `node tests/wasm-smoke/smoke.mjs` (61 checks); `cargo deny --workspace
--all-features check` (advisories ok, bans ok, licenses ok, sources ok); and the
CLI smoke `props T P 101325 Q 0 IF97::Water` → `373.12430000048056`.

The **seventh** ignored test is not a coverage gap and not part of the weekly
sweep: `crates/rustprop-svdsbtl/tests/roundtrip.rs::eval_from_env_artifact` is a
scratch harness that requires `RUSTPROP_SVDS` to point at a locally built
upstream table.

What could **not** be run here, and is therefore only ever green on a runner:
the macOS and Windows CLI builds in `release.yml`'s `cli-binaries` matrix, the
five-preset × three-target wasm bundle matrix, the crates.io upload itself, and
— new with R10 — the `schedule-keepalive` job *as a job* (it needs a `schedule`
event). Note the distinction for `supply-chain`: its **content** was already
fully proven here (the pinned tarball was re-downloaded, its sha256 matched the
digest `ci.yml` hard-codes, and that exact binary produced the verdict above);
only the YAML wrapping it was unexercised, and that ran green on 2026-08-21.

**`platform` is no longer unknown — it has run, and it found things.** See
*Platform-scoped assertions* below for what off-Linux/x86-64 execution actually
established. As of run `32470882030` (commit `19003b8`) the whole workflow is
green: `ci`, `msrv`, `supply-chain`, and `platform` on **macos-latest** and
**windows-latest** in **both** profiles plus the CLI smoke.

Read **`RELEASE-CHECKLIST.md`** before doing anything in (a) — it is the
execution order, and it opens with a hazard that would otherwise half-publish
the project:

> **The crates.io new-crate rate limit blocks a one-shot release.** crates.io
> allows a burst of **5** new crates per user, refilling one per 10 minutes
> (`rust-lang/crates.io`, `rate_limiter.rs`). This workspace publishes
> **twelve**, all new, and `cargo publish` does not retry a 429 — so the first
> tag would upload about five and stop, permanently, with the `rustprop` facade
> (11th in the upload order) among the casualties. `release.yml`'s `crates-io`
> job now refuses to start in that state, before any upload. The remedy is a
> rate-limit override requested from the crates.io team *before* tagging;
> `RELEASE-CHECKLIST.md` §0 has the procedure and the alternative.

That hazard was re-derived independently on the merge, not taken on R10's word.
The preflight's crate set (`cargo metadata --no-deps`, `select(.publish==null)`)
resolves to exactly the twelve crates the dry-run uploads and excludes the six
`publish = false` members; the dry-run from cold caches confirms the upload order
is **core, data, heos, if97, pcsaft, svdsbtl, cubics, humid-air, incompressible,
tabular, `rustprop`, `rustprop-wasm`** — so the facade really is 11th, past a
burst of 5, and the packaging order cargo prints first is *not* the upload order.
All twelve names still returned 404 (control `serde` 200), and
`gh variable list` is empty, so `CRATES_IO_NEW_CRATE_BURST` is unset and the
preflight uses its default of 5: **as the repository stands, the `crates-io` job
would refuse to publish.** That is the intended behaviour, not a defect — it is
the guard doing its job, and clearing it is step 0 below.

### (a) Owner-only, and irreversible

**`RELEASE-CHECKLIST.md` is the executable version of this list** — same steps,
in order, with the checks and the resume story. The summary:

0. **Clear the crates.io new-crate rate limit** (see the box above). New since
   R10, and it comes first because everything after it is irreversible.
1. **Re-check the twelve crates.io names.** crates.io has no reservation
   mechanism — a name is claimed by the first successful publish, so this is a
   check, not an action. All twelve returned **404** on 2026-08-19 (re-checked
   on the R10 merge), with `serde` at 200 through the same request as the
   control: `rustprop`,
   `rustprop-core`, `rustprop-data`, `rustprop-heos`, `rustprop-if97`,
   `rustprop-cubics`, `rustprop-incompressible`, `rustprop-pcsaft`,
   `rustprop-tabular`, `rustprop-svdsbtl`, `rustprop-humid-air`,
   `rustprop-wasm`. (`rustprop-cli` is `publish = false` — it is the example
   app, not a deliverable crate, which also means `cargo install rustprop-cli`
   will not work.)
2. **Add the `CARGO_REGISTRY_TOKEN` repository secret** that
   `.github/workflows/release.yml` reads. `gh secret list` is still **empty**.
3. **Push `main`. DONE, 2026-08-21** — and it earned its place on this list
   twice over. The first push went red on `ci` (a floating toolchain: the
   runner's clippy 1.98 against this box's 1.97.1) and on `platform` (three
   separate off-Linux findings across two more rounds; see *Platform-scoped
   assertions*). The tree that ships is now one CI has actually built: run
   `32470882030`, commit `19003b8`, every job green. The rule stands for any
   future tag — never tag a tree no CI run has built.
4. **Tag `v0.1.0` and push the tag**, which runs verify → crates.io publish →
   five wasm presets × three targets → CLI binaries for linux-x64 /
   macos-arm64 / windows-x64 → GitHub release.

```bash
git push origin main            # let CI go green on the real tree first
git tag v0.1.0 && git push origin v0.1.0
```

Nothing has been published and nothing is half-published: **no tags exist**
locally or on `origin`, there are **no GitHub releases**, and the `release`
workflow has **never run** — the Actions history is CI runs only. A crates.io
version can be yanked but never replaced, which is why the `verify` job re-runs
fmt, clippy, the release-profile suite and the dry-run before anything uploads.

One de-risking step is available to the owner and costs nothing irreversible:
after pushing `main`, **`workflow_dispatch` the release workflow**. The
`crates-io` job is gated on `startsWith(github.ref, 'refs/tags/v')`, so a
dispatch run exercises verify + the five wasm bundles + all three CLI targets
and publishes nothing. It has never been run even once; the first real tag
would otherwise be the first execution of that pipeline.

### (b) Still doable in a session, and nobody has done it

**Essentially nothing — and that is the finding.** The local verification
surface is exhausted: every command in "The gate" above, plus all six heavy
suites, plus the MSRV checks, the node smoke and the cold-cache dry-run, are
green on this tree. What is left cannot be reached from a workstation:

- **`schedule-keepalive` still needs a `schedule` event**, which cannot be
  simulated here. `platform` and `supply-chain` no longer belong on this list:
  both have run, on the pushed tree, and both are green (run `32470882030`).
- **`release.yml` needs a dispatch**, which is a workflow run — an owner action
  by policy, written out in `RELEASE-CHECKLIST.md` §3 rather than performed.
- **Cross-compiling the off-Linux targets locally is not a substitute** even if
  the toolchains were installed: the risk `platform` exists to cover is a
  *libm* difference in the last bits of `exp`/`log`/`pow`, which only shows up
  when the suite actually executes on that platform.

The one open judgement call that needs a person rather than a machine: the
**archived oracle wheel** (10.3 MB, `.git` 64 → 75 MB). It sits alone in commit
`6f50785` specifically so it can be dropped before the first push if the owner
would rather not carry a binary in history; the PLAN.md entry names the exact
sentences that would then need deleting.

### (b′) The R10 round — landed 2026-08-19, merged as `5c93a15`

The Wave-4b R10 round is **done** and merged to `main`. What it produced,
and what is left in each item:

- **Supply chain — closed.** `deny.toml` is an exact licence allowlist over the
  29-crate external graph (all permissive: `0BSD`, `Apache-2.0`, `MIT`,
  `Unicode-3.0`, `Unlicense`, `Zlib`), zero advisories, `sources` denying git
  and non-crates.io registries. `ci.yml`'s `supply-chain` job installs
  cargo-deny 0.20.2 from a sha256-pinned tarball rather than an unpinned Docker
  action. The first run found two real problems in `tools/dbg-pcsaft` (no
  licence field, three wildcard path deps), fixed at source.
- **Platform coverage — jobs added; HAVE NOW RUN, and found three things.**
  `ci.yml`'s `platform` job runs the golden suites on `macos-latest` and
  `windows-latest` — debug on every event, release on pushes to `main`, plus a
  CLI smoke. `.gitattributes` pins LF so the Windows checkout cannot mangle the
  fixtures. It went red on its first two runs and is green on the third; what
  it established is in *Platform-scoped assertions* below. It also gained
  `--no-fail-fast`, because `cargo test` stops at the first failing test BINARY
  and that hid a whole round's worth of findings behind one failure per host.
- **Oracle archival — closed.** `requirements.txt` pins CoolProp by sha256 for
  `--require-hashes`; `tools/golden-gen/ORACLE.md` records the full identity
  and a cold-start walkthrough; `verify-oracle.sh` answers "is this THE oracle"
  and reads its expected digest from `manifest.json` so it cannot go stale; and
  the wheel itself is archived at `tools/golden-gen/vendor/` (10.3 MB, MIT,
  dev-only — see the PLAN.md entry for the size-vs-reproducibility judgement,
  and note it is a separate commit so the judgement stays reversible).
- **A release checklist — closed, and it found the release's real blocker.**
  `RELEASE-CHECKLIST.md` is the owner's execution order. Writing it surfaced
  the crates.io new-crate rate limit quoted above, which no amount of testing
  the tree could have found.

Also closed:

- **`manifest.json`'s drift**, fixed at the generator. Its `files` key came
  from the in-process `WRITTEN` list, so regenerating one fixture *shrank* the
  manifest; `write_manifest()` now scans the fixture directory (and records the
  oracle's sha256). It lists all 123 files, and refreshing it no longer costs a
  ~100 s tabular rebuild. Verified by counting rather than by reading the claim:
  the manifest's `files` key and `ls tests/golden/fixtures/*.jsonl` are the same
  123 names, set-for-set, with nothing missing on either side.
- **`.gitignore`** now covers `release.yml`'s `pkg-<preset>-<target>` fan-out
  and the `rustprop-{wasm,cli}-*.tar.gz` archives it packages at the repo root.
- **`CHANGELOG.md`** exists, with the divergences a consumer can actually hit
  separated from the ones only a digit-by-digit comparison would show.

The ranked candidate work further down this file (the answer-vs-refuse residue
at unphysical inputs, the mixture-side memo, consumer documentation, a fifth
sweep stream) is unchanged and none of it is required for a release either.

**What changed since the Wave-4b statement.** Only the release *apparatus*
changed; **the shipped numbers did not move**. No crate source, fixture or
golden value differs from the Wave-4b tree — R10 touched workflows, docs, the
generator's manifest writer, `deny.toml`, `.gitattributes`, `.gitignore` and
`tools/dbg-pcsaft`'s manifest, and the record/fixture/test counts are unchanged
at 41,629 / 123 / 134. Concretely:

1. **A release blocker was found that testing could never have surfaced** — the
   crates.io new-crate rate limit. Wave-4b's statement said the tree was ready
   to tag; on the evidence it was ready to *half-publish*. This is the single
   most important delta.
2. **The supply-chain gate went from absent to real**, and is now part of the
   gate command list rather than an aspiration.
3. **The oracle stopped being un-reproducible.** Wave-4b identified it as
   "CoolProp 8.0.0"; it is now pinned by the sha256 of the actual compiled
   extension, archived, and checkable with one script.
4. **Off-Linux went from untested to *scheduled to be tested***. That is a
   smaller step than it sounds: the coverage does not exist until CI runs it.
5. **The manifest stopped lying** (111 of 123 → 123 of 123), fixed at the
   generator so the documented one-file regeneration recipe can no longer
   shrink it.
6. **Two counts corrected in the owner's favour**: "claim the crates.io names"
   was wrong — crates.io has no reservation mechanism, so it is a re-check
   immediately before tagging, not an action that can be completed in advance;
   and `cargo install rustprop-cli` will not work, because the CLI is
   `publish = false` and reaches users only as a release binary.

### (c) Known and accepted gaps

Not restated here — they have their own tables, and each is asserted somewhere
so it can neither widen nor silently heal:

- **Known divergences from upstream** (the table below), now including R12's
  MethylLinoleate row: one state in 157,374 where the wheel answers and the
  port refuses.
- **Unported by design** (the list below): the SVDSBTL builder and its critical
  patch, mixture phase-envelope machinery and mixture `HS_flash`, third-order /
  PSI mixture derivatives, Ammonia's Tillner-Roth alternate EOS, the REFPROP
  backend.
- **The tabular msgpack+zlib disk cache**, which costs ~100 s per process for a
  LogPH build. WASM has no home directory; this is the one gap a consumer
  actually feels.
- **`HAPropsSI` returns `Result` rather than upstream's `+inf`-with-a-global**,
  and pseudo-pure fluids serve PT/PQ/QT plus the four classic-ancillary caloric
  pairs while the rest refuse loudly.

---

## Known divergences from upstream

Each of these is asserted somewhere, so it can neither widen nor silently
heal. Do not "fix" one without checking the assertion that pins it.

| Divergence | Where | Why it stands |
|---|---|---|
| THREE pinned mixture records where the port answers and the wheel's recorded value is provably not the wheel's own equilibrium (mixture HSU_P shared-state corruption ×2; shallow-TPD metastable root ×1) | `acceptance.rs` `mixture_divergences`, each pinned to the PORT's value with heal detection | Upstream's HSU_P residual mutates the shared backend (a Tmax-endpoint PT evaluation corrupts SatL/SatV and disables the two-phase split for the rest of the solve); a fresh wheel PT flash at the port's converged T reproduces the port BITWISE. The corruption is history-dependence the port's stateless flashes deliberately cannot express |
| Upstream's `PT_flash` serves every property off the density solver's LAST TRIAL, not the root it returns; the port evaluates at the root, so its density matches BITWISE while h/s/cp/w differ by up to 5.0e-8 near the critical point | `heos_pt.rs::stale_cache_allowance` — a per-record bound, `\|dX/dln ρ\| · ftol/\|dln p/dln ρ\|`, replacing the old blanket 1e-8 tier for Cp/A; the 218 density records are asserted to 1e-13 (`DENSITY_RTOL`), the premise the bound is derived from. Not bitwise: a PT root is not a single double — Ammonia at T=425.838 K, P=17.045 MPa lands 42 ulp away in a debug build and neither density is wrong, so the bitwise COUNT is printed (217/218 debug, 218/218 release) rather than asserted (`d5a7331`) | `FlashRoutines.cpp:336` assigns `_rhomolar` and nothing else, while `SolverTPResid::call` (`HelmholtzEOSMixtureBackend.cpp:2803`) mutated the same backend at every trial; the calculators then recompute `_delta` fresh but read `dalphar_*` from those caches, so a wheel PT state is internally inconsistent. Reproducing it would mean threading a stale-iterate cache through a stateless port to serve knowingly worse numbers — same ruling as the mixture rows above |
| The imposition-clear channel of the same corruption (upstream's stability feed fallback permanently clears SatL's constructor liquid imposition for the backend instance): the wheel's Nitrogen[0.97]&Water[0.03] sweep-pair inversions contradict its own forward flashes by up to 1.7 K or error outright | `hsu_p_imposition_clear_divergence_pinned` in `tests/golden/tests/mixtures.rs`, asserting the port's self-consistent inversions | Only observable through sweep-pair flashes (scalar PT builds a fresh backend per call); reproducing it would make PS/HP flashes history-dependent — same ruling as the row above. Full mechanism in the 2026-08-17 Decisions entries |
| `HAPropsSI` errors return `Result` instead of upstream's `+inf`-with-a-global | humid-air suite | A global error slot is not a thing a WASM library should have |
| PC-SAFT `WATER` PT/DT errors loudly | PC-SAFT suite, error parity asserted | Upstream computes on children whose sigma is still the −1 sentinel and returns garbage densities |
| Tabular msgpack+zlib disk cache under `~/.CoolProp/Tabular` not ported | documented in `PLAN.md` | No home directory in WASM. Cost: a LogPH table build runs ~100 s per process (40k HP flashes) — exactly the cost upstream's cache exists to avoid |
| Pseudo-pure fluids serve PT/PQ/QT plus the classic-ancillary (H,P)/(P,S)/(P,U)/(D,P) flashes (Wave-2 R6, goldened over six regions by R7); the remaining pairs (DmolarT, HS, DQ, HSU_D, ...) are loud `NotImplemented` | pseudo-pure suite (665 records: 654 value at 1e-8, 11 error-parity), verbatim error parity | Upstream routes the rest through legacy solvers that are dead code for the 130 superancillary fluids |
| R507A gas-classified caloric (P,X) states at p = 0.995·max_sat_p error loudly where the wheel converges | PLAN.md 2026-08-17 R6 entry | Both implementations fire the same gas stability retry from 1e-6 at the Tmin probe (TVanc−0.01, 0.24 K below Tcrit); the retry trajectory is chaotic through the vdW loop and only bitwise alphar arithmetic would reproduce the wheel's lucky convergence (EOS sums agree at 1e-13). A needle: 0.9925·pmax and 0.9975·pmax both agree at ≤2e-9 |
| Pseudo-pure PY-flash refusal MESSAGES are the port's own bracket diagnostic, not upstream's "unable to solve 1phase PY flash with Tmin=…, Tmax=… due to error: …" wrapper | `pseudo_pure_error_parity_matches_upstream` (refusal-vs-answer asserted for all 11 error records; oracle text carried in each record's `error` field) | The wrapper is `HSU_P_flash`'s catch around the single-phase solve and what it quotes is the INNER diagnostic. R8 closed the bisection stand-in (the inner solve is upstream's own TOMS748 now), but the no-bracket derivative path and the 2-D Newton fallback a refusal falls through to are still unported, so the text a refusal carries is still the port's bracket diagnostic. `post_update` and `solver_rho_Tp` refusals in the same suite ARE verbatim |
| Upstream's `post_update` validity gate is ported for the HEOS pure/pseudo-pure arm only, not for `mixture_update` | `props_api.rs::post_update` | No reachable NaN mixture state has been observed, and the 10f divergence pins would need re-validating; a gate without evidence is an invented guard |
| `post_update`'s refusal TEXT at `PropsSI("Hmass","T",0,"Smass",101325,"Water")`: rustprop "rhomolar is not a valid number", the wheel "rhomolar is less than zero" | `post_update_refusal_text_divergence_pinned` in `tests/golden/tests/validity.rs` | Both arms are carried, in upstream's order, and refusal-vs-answer agrees; the two trip DIFFERENT arms because the port's `(Smass,T)` flash leaves `rhomolar` NaN where upstream's leaves it negative. Matching the text would mean reproducing a divergent iteration's garbage bit for bit, on a path where neither implementation has an answer |
| At inputs far outside the fluid's range, refusal-vs-answer agrees but the alpha0 message's `tau`/`delta` can differ (9 rows in the R11 scan, all sub-triple `(P,Q)` or negative-`T` `(Q,T)`) | `tests/golden/tests/validity.rs` asserts the 216 records where the states agree; the 9 are described in the 2026-08-18 Decisions entry | Same cause as the row above: the `nTau`/`nDelta` always match — it is the garbage state underneath that differs |
| SVDSBTL evaluator agrees to a few ulp, not bitwise (700 of 745 records bitwise, worst 1.8e-15) | `tests/golden/tests/svdsbtl.rs` | GCC compiles the reference build with `-ffp-contract=fast`. Fusing the obvious candidate makes agreement *worse*, so the contraction sits elsewhere; chasing it would match a compiler flag, not port an algorithm |
| Upstream's `PXFLASH_DIRECT_EOS` cache-bypass for warm (P,X) probes is not ported (R8); the port takes upstream's own `catch (...)` fallback, the cached `update_TP_guessrho` path | `PxResid` doc comment in `flash_px.rs`; the (P,X) golden suites at their 1e-8 policy | It exists to dodge a `CachedElement` layer this port does not have (each residual owns a `DerivsMemo`), and upstream calls it "bit-equivalent within ULP" to the path we do run. Measured cost: 2.0e-16 median displacement over 1 433 (P, caloric) goldens, 608 of them bitwise |
| ONE state in the whole registry where the wheel answers and the port refuses: MethylLinoleate `(P, H/S/U)` at p = 1.001·ptriple = 1.3126e-06 Pa with the caloric input taken from the wheel's own Q=0 value (R12; found by a 157,374-state answer-vs-refuse scan over all 136 fluids, which agrees everywhere else) | `methyl_linoleate_superanc_extrapolation_divergence_pinned` and the pinned count in `px_refusal_parity_matches_upstream` (`tests/golden/tests/px_refusal_parity.rs`), plus the fixture's ±10,000-ulp control rows | This fluid's JSON `satminL` (T=260 K, p=1.3113e-06 Pa) and its superancillary p-curve (p(260 K)=1.4986e-06 Pa) disagree by 14%, and the phase-determination band gates on the JSON value — so `get_T_from_p` EXTRAPOLATES the inverted-pressure expansion below its own domain, on both sides. The residual coefficient noise (Eigen's `L*f` reduction order vs a sequential sum; ~2e-15 inside the domain, 2.0e-12 extrapolated here) then decides the SIGN of `q = (u−u_L)/(u_V−u_L)` for an input that sits exactly on Q=0: wheel 0, port −3.4e-12. With rho_V = 1.6e-10 mol/m³ that makes the lever-rule density negative, and upstream's own `post_update` refuses it — 10,000 ulp lower the WHEEL refuses identically. Healing it means matching Eigen bit for bit, same ruling as the SVDSBTL row |

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

## Platform-scoped assertions

These are NOT divergences from upstream — the port and the wheel agree on
Linux/x86-64 at every state below. They are the places where a fixture cannot
speak for a platform it was never generated on, and the suites say so
explicitly rather than pretending.

**The premise.** Every one of the 41,629 oracle records came from a
**Linux/x86-64** CoolProp 8.0.0 wheel. Apple's, MSVC's and glibc's libm are
each free to return different last bits for `exp`/`log`/`pow`/`sinh`/`asinh`,
and nothing in this repo can say what CoolProp built against a different one
would have answered. Where that decides the result, the honest scope of the
assertion is the platform the oracle ran on.

**What off-Linux execution has established so far** (`ci.yml`'s `platform`
job, first run 2026-08-21):

| Finding | Hosts | Resolution |
|---|---|---|
| `heos_pt`'s PT-density assertion was too tight at 1e-13; the host libm moves an iterative density root (worst: R22 `Dmolar(T=367.44852500798726, P=9980000.000913477)`, 9.82e-13 relative) | macOS-arm64 **and** Windows-x64 | `DENSITY_RTOL` relaxed to 1e-11 on every platform — 10x the observed worst, still three orders below the 1e-8 scale `stale_cache_allowance` derives from (`3a94971`) |
| `PR::Propane` PT at T=1e30 and T=1e20, P=101325: the returned density is the rounding residue of a `-713.06...` cancellation, an integer number of ulps. Linux +4 ulp / Windows +2 ulp at 1e30; Linux exactly 0 (so `delta = 0`, so the alpha0 `ValidNumber` gate fires) / Windows +1 ulp (so it answers) at 1e20 | Windows-x64 | Per-STATE allowlist in `tests/golden/tests/validity.rs`, consulted only off Linux/x86-64, skip counts ASSERTED (132/48/6). The two `SRK::Propane` states at the same temperatures are listed too — Windows agreed there, but the measured margin is two ulps, and at T=1e20 the SIGN of the residue is what produces the refusal 31 records assert (`19003b8`) |
| SES36 `(P, Smolar)` at P=2846151 (0.1% under `max_sat_p`): `t_l == t_v`, `rho_l`/`rho_v` 13 ulp apart, and the two entropies BITWISE EQUAL — so the lever rule is 0/0 = NaN and `post_update` refuses "T is not a valid number", exactly as the wheel does. macOS puts one ulp between them and answers the fixture's own neighbouring `(P, Umolar)` values | macOS-arm64 | Per-STATE allowlist in `tests/golden/tests/pseudo_pure.rs`, same mechanism, skip count ASSERTED (2) (`19003b8`) |

**The rule the allowlists follow**, and the reason they are not a blanket
skip: *keep asserting values everywhere; stop asserting knife-edge
refusal-versus-answer classification off Linux/x86-64.* Each entry is one
enumerated state carrying its measured evidence, `cfg!`-gated so Linux/x86-64
skips nothing, and each test asserts its exact skip count — so a set that
grows turns the suite red instead of quietly shrinking its coverage. Every
other record in those fixtures (the HEOS alpha0 gates at T=1e30, the
melting-line refusals, the negative-pressure `post_update` arm, the cubic
states at 300 K and 500 K) keeps asserting on all three platforms.

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
cargo test --workspace --all-features             # DEBUG — what ci.yml runs
cargo test --workspace --all-features --release   # RELEASE — release.yml's verify job
cargo run -q -p rustprop-datagen && git diff --exit-code   # datagen determinism
cargo build -p rustprop --features all-backends --target wasm32-unknown-unknown --release
cargo publish --dry-run --workspace
cargo deny --workspace --all-features check    # licences + advisories (deny.toml)
node tests/wasm-smoke/smoke.mjs        # after a wasm-pack --target nodejs build
```

> **Run the test suite in BOTH profiles, every time.** `ci.yml` runs debug and
> `release.yml`'s verify job runs `--release`, so a one-profile gate can pass
> here and fail on the tag. It has already happened: `heos_pt` asserted bitwise
> on a PT density that debug and release round to doubles 42 ulp apart — both
> reproducing the requested pressure — and the assertion had to be relaxed to
> 1e-13 (`d5a7331`). Anything an iterative solver produces can split this way.

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

The acceptance sweep is seeded, so widening it never invalidates what is
already there. It now runs **four** independent streams, each frozen at the
moment it landed: `20260807` (the original 3,720), `20260816` (+1,765),
`20260817` (+535) and `20260818` (+360 pseudo-pure caloric pairs). **Add a new
section with a fresh `random.Random(...)` rather than touching an existing
stream** — every earlier record must stay bitwise identical, and the check is
that regeneration produces a diff of pure insertions. Raising `N_PER` inside
the original stream is also safe (the first N draws are unchanged).

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

Everything ranked #1 or #2 here since 2026-08-16 has been done or deliberately
pinned, in this order: widen the sweep; audit the invented-guard classes; the
output tail; the two mixture latents; the cubic sub-pascal secant; the
`Ok(non-finite)` validity gap; and (Wave-3 R9) the unread fixture batteries.
Each has a dated entry in PLAN.md's Decisions log with its evidence — that log,
not this list, is the record.

Two standing answers worth not re-deriving:

- **Input abuse does not belong in the seeded acceptance sweep.**
  `sweep_propssi` skips every state the wheel rejects, by design, so it is
  structurally blind to answer-vs-refuse defects. `validity.rs` (1,626 records,
  0.02 s, runs in the ordinary `cargo test`) is that dimension's home.
- **The sweep's value is in widening it, not in re-running it.** Its 2026-08-18
  widening (pseudo-pure caloric pairs) was the first that found nothing — a
  clean verdict on Wave 2's work, not a reason to stop widening.

### 1. Close the answer-vs-refuse residue at unphysical inputs (medium)

The R11 abuse scan left 136 rows where the port answers and the wheel raises,
and 314 where the port raises and the wheel answers. None is a validity-gate
question, all are pre-existing, and none is reachable at physical inputs. By
count:

| n | class | what happens |
|---|---|---|
| 148 | `(Dmolar,P)` below the triple-point pressure | port: loud `NotImplemented`; wheel: answers |
| 130 | pseudo-pure `DmolarT` / `SmolarT` / `HmolarSmolar` | port: loud `NotImplemented` (by design — see the divergence table); wheel: answers |
| 103 | sub-triple `(P,Q)`, e.g. Water at 1e-30 Pa | upstream's ancillary extrapolation yields a NEGATIVE density and `post_update` refuses; the port converges to a positive one (8.4e204 mol/m³) and answers garbage |
| 36 | `d(Hmolar)/d(T)\|P` and the other derivative OUTPUT strings | `Param::parse` does not accept them, so `props_si` rejects outright while the wheel answers. The machinery exists (`rustprop_heos::derivs`, 207 goldens) — only the string parser is missing, which makes this the cheapest item here and the only one that is a real missing feature |
| 22 | `(T,Dmolar)` at T ≤ 0 | the port's flash reports a TWO-PHASE state and serves the ancillary pressure; upstream's pressure is NaN and `post_update` refuses |
| 10 | `solver_rho_Tp` supercritical-liquid failures at T = 1e-30 | upstream gives up, the port finds a root |

Everything except the derivative-string row means reproducing a divergent
iteration's garbage bit for bit, which is only worth doing if a real caller
lands there.

### 2. Performance (measured; the obvious targets are gone)

This is no longer unprofiled ground — `tools/perf-bench` and `PERF.md` exist,
and the two biggest costs have been paid: the derivative-matrix memoization
(Wave 1) and the TOMS748 + warm-carry closure (Wave 2) together took the HP
flash from ~963 µs to ~80 µs and the LogPH-backed suites from ~200 s to
~13 s. What remains is smaller and should be measured before it is touched:
the mixture-side stability pairs in `mixture_vle.rs` still evaluate
`dpdrho`/`d2pdrho2` separately and would take the same bit-identical memo
pattern.

### 3. Documentation for consumers (small)

The README quickstart is real and doc-tested. What does not exist: per-engine
guidance on *which* backend to choose, and a worked WASM example beyond the
size table.

### 4. Widen the sweep again, at the pairs no stream draws (small)

The four seeded streams cover PT/PQ/QT/DT/HP/PS broadly and the pseudo-pure
caloric pairs, but several PropsSI-reachable pure-fluid pairs have only
hand-chosen goldens: `(D,H)`/`(D,S)`/`(D,U)` (HSU_D), `(H,T)`/`(T,U)`
(generalized DHSU_T), `(P,U)` for the 130 pure fluids, and `(D,Q)`. Add a
fifth stream rather than touching an existing one.

---

## Orientation for a new session

- `PLAN.md` — phases, checkboxes, and the Decisions log. The log is
  append-only and is where every non-obvious choice is justified. Read the
  entries for whatever you are about to touch.
- `CLAUDE.md` — the working status summary, kept current at phase gates.
- `RELEASE-CHECKLIST.md` — the owner's pre-tag execution order, the resume
  story for a partial crates.io publish, and everything about `release.yml`
  that has never been observed running.
- `CHANGELOG.md` — what v0.1.0 ships, what it deliberately does not, and which
  divergences a consumer can actually reach.
- `tools/golden-gen/ORACLE.md` — which exact CoolProp binary every golden came
  from, why a version number does not identify it, and how to reproduce it from
  cold. Read before touching a fixture.
- `WASM-SIZES.md` — measured bundle sizes, regenerate with
  `tools/wasm-size-table.sh`.
- Upstream lives at `~/homecloud/dev/CoolProp`, pinned at tag `v8.0.0`. It is
  the correctness oracle; when in doubt, read its source rather than
  reasoning about what it probably does. Several findings in the Decisions
  log exist because the code disagreed with the documentation.
