# Release checklist — v0.1.0

For the owner, in execution order. Everything before §4 is reversible;
everything from §4 on is not.

**What "irreversible" means here.** A crates.io version can be *yanked* but
never deleted and never replaced: the name, the version number and the exact
tarball are permanent. Yanking only stops *new* dependency resolution — anyone
with a lockfile keeps downloading it. There is no second first release.

Read `NEXT-STEPS.md` §"Release readiness" alongside this: it says what is
green, what is open, and what is an accepted gap. This file says what to *do*.

---

## 0. Before anything: the rate limit that will bite

**This is the one thing most likely to break the first release, and it is not
a bug in the tree.**

crates.io rate-limits the creation of *new* crates per user: a burst of **5**,
refilling **one per 10 minutes**
([`rate_limiter.rs`](https://github.com/rust-lang/crates.io/blob/main/src/rate_limiter.rs),
`LimitedAction::PublishNew`). This workspace publishes **twelve** crates, all
twelve of which are new. `cargo publish` does not retry an HTTP 429.

So a naive `cargo publish --workspace` on the tag would upload roughly five
crates, be refused for the sixth, and stop — leaving a **partial, permanent**
publish, with the facade `rustprop` (the crate anyone would actually depend on)
among the ones that did not make it.

`release.yml`'s `crates-io` job now refuses to start in that state: a preflight
step counts how many of the twelve are new and fails *before the first upload*
if that exceeds the allowance. It can only ever refuse to publish; it cannot
cause one.

**Do one of these before tagging:**

- **Preferred — get an override.** Ask the crates.io team to raise your
  `PublishNew` limit, saying you are releasing a 12-crate workspace in one
  shot. The current route is an issue on
  [`rust-lang/crates.io`](https://github.com/rust-lang/crates.io/issues) using
  their rate-limit-increase template, or `help@crates.io`. Allow days, not
  hours. When it is granted, set the repository variable to the allowance you
  were given, which disarms the preflight:

  ```bash
  gh variable set CRATES_IO_NEW_CRATE_BURST --body 12 -R ernsoylu/RustProp
  ```

- **Or publish by hand, paced, then tag.** Publish the twelve yourself from a
  clean checkout in the order in §6, waiting ~10 minutes between each after the
  first five. Then the tag's preflight sees zero new crates and passes, and
  `cargo publish --workspace` has nothing left to do — *except* that it will
  fail on "crate version already uploaded". If you take this route, expect the
  `crates-io` job to fail and to need the §6 resume treatment, and understand
  that `github-release` will not run until it succeeds.

The preferred route is preferred precisely because it keeps the pipeline the
single thing that publishes.

*Not verified from the development machine: the burst/refill numbers are the
defaults in crates.io's source as of 2026-08-19. The deployed service may run
different values, and per-user overrides exist. The preflight measures the
input (how many crates are new), not the live limit.*

---

## 1. The gate must be green, locally, in both profiles

Run all of it. `ci.yml` runs the debug suite and `release.yml`'s `verify` job
runs the release suite, so a one-profile pass can go green here and red on the
tag — it has happened (`d5a7331`).

```bash
export PATH="$HOME/.cargo/bin:$PATH"

cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features             # DEBUG
cargo test --workspace --all-features --release   # RELEASE
cargo run -q -p rustprop-datagen && git diff --exit-code   # datagen determinism
cargo build -p rustprop --features all-backends \
    --target wasm32-unknown-unknown --release
cargo publish --dry-run --workspace
cargo deny --workspace --all-features check       # licences + advisories
```

Then the wasm deliverable end to end:

```bash
wasm-pack build crates/rustprop-wasm --target nodejs --out-dir pkg-node \
    --features if97,heos,water,humid-air
node tests/wasm-smoke/smoke.mjs
```

And the six `#[ignore]`d heavy suites, which no push-triggered job runs:

```bash
cargo test -p rustprop-golden-tests --test heos_all_smoke -- --ignored --nocapture
cargo test --release -p rustprop-golden-tests --test mixtures mixture_sweep -- --ignored
cargo test --release -p rustprop-golden-tests --test tabular_state tabular_pairs -- --ignored
cargo test --release -p rustprop-golden-tests --test acceptance -- --ignored --nocapture
cargo test --release -p rustprop-golden-tests --test acceptance_tabular -- --ignored
cargo test --release -p rustprop-golden-tests --test acceptance_svdsbtl -- --ignored
```

> If `cargo publish --dry-run --workspace` fails with a baffling `E0432` after
> any cross-crate API change, it is the overlay registry's immutable cache, not
> your code. Purge and re-run:
> `rm -rf target/package ~/.cargo/registry/{src,cache}/*/rustprop-*`

Checks:

- [ ] fmt, clippy, both test profiles, datagen determinism, wasm build, dry-run, cargo-deny — all zero
- [ ] Node smoke passes (61 checks)
- [ ] all six heavy suites pass
- [ ] `git status` clean

## 2. Push `main` — this is a precondition, not a nicety

`origin/main` is **seven commits behind** local `main` as of 2026-08-19, so CI
has never built the tree that would ship. Tagging first would make the release
pipeline the first thing that ever compiled it.

```bash
git push origin main
```

Then watch CI, and note that **three of its jobs have never executed anywhere**
— `platform` (macOS + Windows), `supply-chain` (cargo-deny), and
`schedule-keepalive`. They were written on a Linux box that cannot run a macOS
or Windows runner.

```bash
gh run watch -R ernsoylu/RustProp
```

- [ ] `ci` green
- [ ] `msrv` green
- [ ] `supply-chain` green
- [ ] **`platform` green on BOTH macOS and Windows** — if it is red, read the
      failing assertion before assuming the workflow is wrong. A golden
      mismatch here is a real finding: it means a property value differs off
      Linux/x86-64, which nothing in this project has ever measured. Do not tag
      on a red `platform` job; nothing forces you to stop, because `release.yml`
      gates only on its own ubuntu `verify` job, so this one is on you.

## 3. Rehearse the release pipeline — recommended, publishes nothing

`release.yml` has **never executed in any mode**. The first tag should not also
be the first run of the pipeline.

```bash
gh workflow run release.yml -R ernsoylu/RustProp --ref main
gh run watch -R ernsoylu/RustProp
```

**Why this is safe.** The publishing job is gated on the ref being a tag:

```yaml
  crates-io:
    needs: verify
    # Publication must never fire from a workflow_dispatch (which exists for
    # verify-only runs) — only from a pushed vX.Y.Z tag.
    if: startsWith(github.ref, 'refs/tags/v')
```

A `workflow_dispatch --ref main` gives `github.ref = refs/heads/main`, so
`crates-io` is **skipped**. `github-release` declares
`needs: [crates-io, wasm-bundles, cli-binaries]`, and a job whose dependency was
skipped is skipped too — so no GitHub release is created either. What the
rehearsal *does* exercise is everything expensive and untested: `verify` (fmt,
clippy, the release-profile suite, the dry-run), all five wasm presets × three
targets, and CLI builds for linux-x64 / macos-arm64 / windows-x64. The
tag-version check also skips itself on a non-tag ref, by design.

- [ ] rehearsal run: `verify`, `wasm-bundles` (5 presets), `cli-binaries` (3 targets) all green
- [ ] `crates-io` and `github-release` both show as **skipped**

## 4. Owner-only, and from here nothing can be undone

### 4a. The token

`release.yml` reads `secrets.CARGO_REGISTRY_TOKEN`. It is **absent** from the
repository secrets as of 2026-08-19.

Create a scoped token at <https://crates.io/settings/tokens> — scopes
`publish-new` and `publish-update`, no crate-name restriction (the twelve
crates do not exist yet, so they cannot be named in a scope), then:

```bash
gh secret set CARGO_REGISTRY_TOKEN -R ernsoylu/RustProp
```

- [ ] `gh secret list -R ernsoylu/RustProp` shows `CARGO_REGISTRY_TOKEN`

### 4b. The names

crates.io has **no reservation mechanism** — a name is claimed by the first
successful publish. "Claiming" is therefore a check, not an action. All twelve
returned 404 on 2026-08-19 (with `serde` at 200 as the control). Re-check
immediately before tagging, because someone else claiming `rustprop` in the
interim is the one failure this checklist cannot recover from:

```bash
for c in rustprop rustprop-core rustprop-data rustprop-heos rustprop-if97 \
         rustprop-cubics rustprop-incompressible rustprop-pcsaft \
         rustprop-tabular rustprop-svdsbtl rustprop-humid-air rustprop-wasm; do
  printf '%s  %s\n' \
    "$(curl -sS -o /dev/null -w '%{http_code}' -A "rustprop-release-check" \
       "https://crates.io/api/v1/crates/$c")" "$c"
done
```

- [ ] all twelve return 404

(`rustprop-cli` is `publish = false`. It is the example app, not a deliverable
crate — which also means **`cargo install rustprop-cli` will not work**. The
CLI reaches users only as a prebuilt binary attached to the GitHub release.)

### 4c. Tag

The `verify` job checks that the tag matches the workspace version, which is
`0.1.0` in the root `Cargo.toml`.

```bash
git tag v0.1.0
git push origin v0.1.0
```

- [ ] tag pushed; `release` workflow started

## 5. What the tag runs

```
verify ──┬── crates-io ─────┐
         ├── wasm-bundles ──┼── github-release
         └── cli-binaries ──┘
```

| Job | Does | Irreversible? |
|---|---|---|
| `verify` | tag-vs-version check, fmt, clippy, `--release` suite, `publish --dry-run` | no |
| `crates-io` | rate-limit preflight, then `cargo publish --workspace` | **YES** |
| `wasm-bundles` | 5 presets × 3 targets, `all-backends` node smoke, uploads artifacts | no |
| `cli-binaries` | linux-x64 / macos-arm64 / windows-x64, uploads artifacts | no |
| `github-release` | attaches every `*.tar.gz`, `generate_release_notes: true` | reversible (delete the release) |

`concurrency: release-${{ github.ref }}` keeps a re-run queued behind an
in-flight one rather than interleaving uploads.

Note `github-release` runs `softprops/action-gh-release@v2` with
`generate_release_notes: true`, which writes GitHub's auto-generated commit
summary — **not** `CHANGELOG.md`. Paste the changelog section into the release
body afterwards if you want it there.

## 6. If `crates-io` fails partway — the resume story

This is the only failure that cannot be fixed by re-running.

**Upload order.** Taken from the actual `cargo publish --dry-run --workspace`
output on 2026-08-19 (packaging order differs; this is the *upload* order):

| # | crate | if the run dies after this one, it is live forever |
|---|---|---|
| 1 | `rustprop-core` | types only; harmless alone |
| 2 | `rustprop-data` | generated data; useless without an engine but harmless |
| 3 | `rustprop-heos` | |
| 4 | `rustprop-if97` | |
| 5 | `rustprop-pcsaft` | ← **a default 5-burst runs out about here** |
| 6 | `rustprop-svdsbtl` | |
| 7 | `rustprop-cubics` | |
| 8 | `rustprop-humid-air` | |
| 9 | `rustprop-incompressible` | |
| 10 | `rustprop-tabular` | |
| 11 | `rustprop` | the facade — until this is up, `cargo add rustprop` fails |
| 12 | `rustprop-wasm` | |

**Diagnose first.** Read the job log for the last successful `Uploading` line,
then confirm against the registry rather than the log:

```bash
for c in rustprop-core rustprop-data rustprop-heos rustprop-if97 rustprop-pcsaft \
         rustprop-svdsbtl rustprop-cubics rustprop-humid-air rustprop-incompressible \
         rustprop-tabular rustprop rustprop-wasm; do
  printf '%s  %s\n' \
    "$(curl -sS -o /dev/null -w '%{http_code}' -A "rustprop-release-check" \
       "https://crates.io/api/v1/crates/$c")" "$c"
done
```

**Do not count on simply re-running the job.** `cargo publish --workspace`
re-attempts every member, and crates.io rejects a version that already exists
("crate version `0.1.0` is already uploaded"). Whether cargo filters those out
before uploading or surfaces the rejection **could not be determined from the
development machine** — neither `cargo publish --help` nor the cargo book says,
and testing it requires a real partial publish. Plan for the failure; if a
plain re-run happens to work, nothing is lost.

**Resume**, from a clean checkout of the tagged commit, excluding what is
already live:

```bash
git fetch --tags && git checkout v0.1.0
export CARGO_REGISTRY_TOKEN=...        # same token
cargo publish --workspace \
  --exclude rustprop-core --exclude rustprop-data --exclude rustprop-heos \
  --exclude rustprop-if97 --exclude rustprop-pcsaft
```

If it was the rate limit, wait ~10 minutes per remaining crate, or publish one
at a time in the order above:

```bash
cargo publish -p rustprop-svdsbtl        # then --cubics, --humid-air, ...
```

Single-crate publishes work *after* the first upload — the pitfall in
`NEXT-STEPS.md` about `cargo publish -p` failing is specific to the state where
nothing is on the index yet, because packaging has to resolve each crate's
`rustprop-*` dependencies against crates.io.

**Then get `github-release` to run.** It needs `crates-io` to succeed. Once the
registry is complete, re-run the failed jobs:

```bash
gh run rerun <run-id> --failed -R ernsoylu/RustProp
```

The preflight will now see zero new crates and pass. If `cargo publish
--workspace` still refuses the already-uploaded versions, the pragmatic finish
is to create the release from the artifacts by hand:

```bash
gh run download <run-id> -R ernsoylu/RustProp -D artifacts
gh release create v0.1.0 -R ernsoylu/RustProp --generate-notes artifacts/**/*.tar.gz
```

**Never** re-tag or force-push a tag to "retry". The version on crates.io is
already fixed to the tarball that went up; a moved tag makes the git history
disagree with the registry permanently.

**Yanking.** If what went up is wrong rather than incomplete:

```bash
cargo yank --version 0.1.0 rustprop      # repeat per crate, deepest last
```

Yank does not free the version number. The next release is `0.1.1`, always.

## 7. After a successful release

- [ ] all twelve crates resolve: `cargo search rustprop`
- [ ] a scratch project builds against the registry, not the path deps:
      `cargo new /tmp/rp && cd /tmp/rp && cargo add rustprop --features if97 && cargo run`
- [ ] `docs.rs` built each crate (check the badge on each crate page; a
      `docs.rs` failure is fixable without a republish only via a new version,
      so read the build log if one is red)
- [ ] the GitHub release carries 5 wasm tarballs + 3 CLI tarballs
- [ ] paste the `CHANGELOG.md` v0.1.0 section into the release body
- [ ] bump the workspace version to `0.1.1` on `main` so the next tag cannot
      collide with a published version

## 8. What this checklist could not verify

Written and gated on Linux/x86-64. Everything below is reasoned from source or
documentation, not observed:

- **`release.yml` has never executed**, in any mode. Every claim about its job
  graph comes from reading the YAML.
- **The `platform` job has never executed.** Whether the golden suites actually
  pass on macOS-arm64 and Windows-x64 is unknown — that is the point of adding
  it, and a first-run failure is a finding about the port, not necessarily
  about the workflow.
- **`supply-chain` and `schedule-keepalive` have never executed.** cargo-deny's
  verdict *was* reproduced locally with the same pinned 0.20.2 binary
  (advisories ok, bans ok, licenses ok, sources ok); the GitHub job wrapping it
  was not.
- **The crates.io rate limit** is quoted from crates.io's source defaults. The
  deployed values are not observable from outside.
- **The crates.io upload itself** — obviously — and therefore the resume
  procedure in §6, which is derived from `cargo publish --help` and the
  observed dry-run ordering, not from a real partial publish. In particular,
  **what `cargo publish --workspace` does when some members are already
  published is unknown**; §6 assumes the pessimistic answer.
- **The macOS and Windows CLI builds and the 5×3 wasm bundle matrix** in
  `release.yml` have never been built anywhere.
