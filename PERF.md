# PERF.md — Helmholtz derivative-matrix memoization

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

## What changed

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

## Bit-identity witness

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
