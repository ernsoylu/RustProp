# PERF.md — hot-path measurements

One section per landed perf change, newest last. All numbers come from
`tools/perf-bench` (standalone harness; release profile: fat LTO,
`codegen-units = 1`) on the same machine — 11th Gen Intel Core i7-1185G7 @
3.00 GHz, 8 threads, **powersave governor**, Linux 7.0.0-28-generic,
rustc 1.97.1. Each section states its own load conditions: the box is shared,
so a measurement taken under concurrent `cargo` builds says so.

## 1. Helmholtz derivative-matrix memoization

Measurements from `tools/perf-bench` (standalone harness, commit 8906b61; release
profile: fat LTO, `codegen-units = 1`). Machine: 11th Gen Intel Core i7-1185G7 @
3.00 GHz (8 threads, **powersave governor** — the harness spins a 2 s global
warmup and a per-bench 0.3 s warmup to ramp it), Linux 7.0.0-28-generic,
rustc 1.97.1 (2026-07-14). Single process, no other load.

BASELINE = worktree at commit 8906b61 with the memoization edits stashed.
AFTER = the same tree with the memoization applied (this commit).

| bench | BASELINE | AFTER | ratio |
|---|---:|---:|---:|
| warm PT Water (Dmolar) | 17 693 ns/op | 7 092 ns/op | 2.49x |
| HP liquid Water (T) | 962 929 ns/op | 379 570 ns/op | 2.54x |
| HP gas Water (T) | 839 126 ns/op | 305 400 ns/op | 2.75x |
| HS gas Water (T) | 54 374 ns/op | 46 400 ns/op | 1.17x |
| QT Water (P) | 244.4 ns/op | 226.0 ns/op | 1.08x |
| mixture PT CH4[0.6]/C2H6[0.4] (Dmolar) | 16.89 ms/op | 6.75 ms/op | 2.50x |
| alphar_all Water (raw kernel) | 2 119 ns | 1 908 ns | 1.11x |
| alpha0_all Water (raw kernel) | 131.3 ns | 114.8 ns | 1.14x |
| LogPT 200x200 grid build | 1.301 s | 0.597 s | 2.18x |
| LogPH 60x60 grid build | 3.139 s | 1.031 s | 3.05x |
| LogPH 200x200 grid build | 36.157 s | 12.402 s | 2.92x |

Reading the table: the raw-kernel rows do not touch the memo at all, so their
~10-14% movement is this box's run-to-run variance band (powersave governor);
the 2.2-3.1x flash and grid-build improvements are far outside it. QT is
ancillary-only and was expected flat. HS gas spends most of its time in the
caloric-superancillary cascade rather than the density solvers, hence the
smaller win.

### What changed

The five `Resid1D` residuals (`SolverTpResid` in flash_pt.rs, `CaloricTResid`
in flash_px.rs, `MixSolverTpResid`, `DpdrhoResid`, and the local `GuessResid`
in the mixture solvers) each recomputed the identical full
`alphar_all(tau, delta)` derivative matrix in `call`/`deriv`/`second_deriv`/
`third_deriv` at the SAME point every solver iteration — Householder4 makes
four same-point calls per iteration, Halley three. Each residual now owns a
`DerivsMemo` (alpha.rs): a single-slot memo keyed on the exact bit patterns
`(tau.to_bits(), delta.to_bits())` holding the last computed matrix, so each
solver point computes the matrix once. The post-Halley/post-Householder
stability checks in flash_pt.rs likewise evaluated `dpdrho_t` and
`d2pdrho2_t` as two separate full-matrix evaluations at one point; they are
merged into `dpdrho_d2pdrho2_t`, which reads both derivatives off one matrix.

### Bit-identity witness

`alphar_all` is a deterministic pure function of `(tau, delta)`, and the memo
key is the exact bit pattern of both inputs, so a memo hit returns f64s
bit-identical to the recompute it elides; the caller-side formulas were
inlined operation-for-operation from the `pressure`/`smolar`/`hmolar`/`umolar`
methods they replace. Witness: `perf-bench --dump-grids` writes every f64 of
the LogPT 200x200 and LogPH 60x60 grid builds (values, all four derivative
matrices per property, transport, nearest-node indices) as `to_bits()` hex
lines — 24.9 MB + 2.2 MB of dumps. The BASELINE and AFTER dumps are
**byte-identical** (`cmp` clean; sha256
`9ba7cb05d12c5a17…` for LogPT and `03221efdc4d01097…` for LogPH on both
sides), and the fold-every-f64 checksums of all three grid builds, including
LogPH 200x200 (`f3e773ea7beeb681`), are unchanged.

## 2. Single-phase HSU_P: TOMS748 + warm-density carry (Wave-2 R8)

This one is a FIDELITY change that happens to be fast: it removed the
documented 30-bit bisection stand-in in `px_solve_single_phase` in favour of
upstream's own TOMS748 plus the warm-density carry across probes. The speedup
is a side effect of running upstream's ~9 probes instead of bisection's ~34.

BASELINE = the pre-R8 tree (the `Merge Wave-2 R7` commit). AFTER = the same
tree with the swap applied. The two binaries were run **back to back**,
alternating, four times each, after a discarded warmup pass of each.

**These numbers were re-measured on a QUIET box at Wave-2 integration**
(load average 1.3-1.4, one `cargo`/`rustc` process, nothing else building).
They supersede the figures taken during R8's own run, when 8 concurrent
`cargo`/`rustc` processes from parallel agents were live (load average 6.7).
The ratio survived that re-measurement essentially unchanged, but the
absolute numbers did not, and the correction is large enough to record:
R8's run estimated the contention penalty at ~6-8% by comparing against
section 1's quiet-box figure, and that estimate was WRONG by roughly an
order of magnitude. The same BASELINE bench measures 276.7-288.7 μs quiet
against the 399.7-409.4 μs recorded under load — a ~40-48% inflation, not
6-8%. Absolute μs/op on this box is strongly load- and frequency-dependent
and should never be compared across sections measured in different sessions;
only the within-session A/B ratio is meaningful.

| bench | BASELINE (median of 4) | AFTER (median of 4) | ratio |
|---|---:|---:|---:|
| HP liquid Water (T) | 283.6 μs/op | 80.2 μs/op | 3.54x |
| HP gas Water (T) | 257.1 μs/op | 68.5 μs/op | 3.75x |
| HS gas Water (T) | 38.5 μs/op | 40.9 μs/op | flat (in the noise) |
| warm PT Water (Dmolar) | 6.2 μs/op | 6.1 μs/op | flat (untouched path) |
| QT Water (P) | 217 ns/op | 210 ns/op | flat (untouched path) |
| alphar_all Water | 1.7 μs/op | 1.6 μs/op | flat (untouched path) |
| mixture PT CH4/C2H6 (Dmolar) | 4942 μs/op | 5136 μs/op | flat (untouched path) |

Per-pass ranges, so the noise is visible: HP liquid BASELINE
276.7 / 280.0 / 287.2 / 288.7 μs against AFTER 75.6 / 77.4 / 83.0 / 88.1 μs;
HP gas BASELINE 244.9 / 249.5 / 264.7 / 273.0 μs against AFTER
67.2 / 67.9 / 69.2 / 69.4 μs. Every AFTER pass is well under the 150 μs target,
and the two distributions do not overlap on either bench.

One correction of substance beyond the scale: `HS gas Water (T)` was recorded
under load as 61.6 -> 57.2 μs and read as a ~1.1x gain. Quiet, it is
38.5 -> 40.9 μs with the AFTER passes spread 37.9-46.2 μs — i.e. flat, inside
the variance band, with no gain to claim. `(H,S)` reaches the (P,X) solve only
through its single-phase cascade, so flat is the expected result.

The remaining rows were left flat on purpose — `warm PT`, `QT`, `alphar_all`
and `mixture PT` do not touch the (P,X) solve at all, so they are the control
rows: their movement is this box's variance band.

### Downstream wall-clock

Everything that builds a LogPH table rides on HP flashes, so the tabular
suites fall with the solve (same `cargo test --release` invocations, same box):

| suite | before | after (under load) | after (quiet box) |
|---|---:|---:|---:|
| `acceptance` (6 020 records) | 64.6 s | 19.7 s | 9.7 s |
| `acceptance_tabular` (1 950 records) | 204.0 s | 17.7 s | 12.6 s |
| `tabular_state` (1 632 records) | 183.5 s | 20.8 s | 8.4 s |

The `before` column was measured under load and has no quiet-box counterpart —
re-measuring it would mean rebuilding the pre-R8 tree for a number the ratio
does not need. Read the third column as the wall-clock a CI runner actually
sees, and the first two only against each other.

### Where the time went

The stand-in bisected `[t_min, t_max]` — often a ~1 500 K supercritical
bracket — to a 2^-29 relative width, which is ~31 halvings, and every one of
them ran a cold SRK-seeded `solver_rho_Tp`. TOMS748 at the same tolerance
converges in a mean of 8.79 probes over the 1 128 (P, caloric) golden solves
(mode 8; 7-9 for the representative subcritical liquid/gas HP and PS solves),
and 2.25 of those probes on average skip the cold solve entirely because the
warm carry seeds Householder4 from the previous probe's density. Supercritical-
pressure solves keep every probe cold by upstream's own
`force_robust_density = (p > p_crit)` rule, which is why the widest brackets
(13-16 probes, 22 solves of 1 128) improve least.

### Accuracy witness

The point of the change: replaying all 1 433 committed (P, caloric) goldens
across 16 fixtures, relative displacement from the wheel moved
median 1.77e-10 -> 2.04e-16, p90 7.28e-10 -> 1.40e-13, p99 2.42e-9 -> 8.72e-10,
max 1.52e-8 -> 1.07e-8, with bitwise agreement rising from 262/1433 to
608/1433. No golden moved past its tolerance, so no fixture was regenerated.
