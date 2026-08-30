//! Port of upstream `include/CoolProp/svd/SVDDecomposition.h` and
//! `SVDEvaluator.h` — the rank-r reconstruction hot path.
//!
//! The builder (`src/SVD/SVDBuilder.cpp`) is deliberately NOT ported: it is
//! Eigen's BDCSVD, and an SVD is unique only up to sign and rotation within
//! degenerate singular subspaces, so no independent implementation
//! reproduces upstream's U and V bitwise. The coefficients are ingested from
//! upstream's own artifacts instead — see the `artifact` module and the
//! PLAN.md Phase 13.1 design note.

// Upstream writes these guards as negated comparisons (`!(a_hi > a_lo)`,
// `!(w[0] < w[1])`) so a NaN bound fails the check instead of passing it.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use crate::hermite::{HermiteBasis, hermite_basis};
use rustprop_core::{Error, Result};

/// Upstream `OutputTransform`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum OutputTransform {
    /// `value = sum_k S_k u_k v_k`
    Identity = 0,
    /// `value = exp(sum_k S_k u_k v_k)` — the matrix was the log of the
    /// property (density, pressure: strictly positive, wide dynamic range).
    Exp = 1,
}

/// Upstream `SlopeSource`. Provenance only — the eval kernel reads slopes
/// from arrays and is agnostic to how they were computed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SlopeSource {
    NaturalCubicSpline = 0,
    HermiteFd = 1,
    Pchip = 2,
}

impl OutputTransform {
    pub fn from_index(i: u8) -> Result<Self> {
        Ok(match i {
            0 => OutputTransform::Identity,
            1 => OutputTransform::Exp,
            other => return Err(Error::Value(format!("unknown OutputTransform {other}"))),
        })
    }
}

impl SlopeSource {
    pub fn from_index(i: u8) -> Result<Self> {
        Ok(match i {
            0 => SlopeSource::NaturalCubicSpline,
            1 => SlopeSource::HermiteFd,
            2 => SlopeSource::Pchip,
            other => return Err(Error::Value(format!("unknown SlopeSource {other}"))),
        })
    }
}

/// Upstream `SVDDecomposition`. The singular values are pre-folded into
/// `v_s` (`v_s[j*rank + k] = sigma_k * V[j, k]`) so the kernel is a plain
/// dot product; `s` is kept for diagnostics only.
#[derive(Clone, Debug)]
pub struct SvdDecomposition {
    pub nx: usize,
    pub ny: usize,
    pub rank: usize,
    pub out_transform: OutputTransform,
    pub slope_source: SlopeSource,
    pub x_grid: Vec<f64>,
    pub y_grid: Vec<f64>,
    /// NX rows x rank columns, row-major.
    pub u: Vec<f64>,
    pub du_dx: Vec<f64>,
    /// NY rows x rank columns, row-major, singular values folded in.
    pub v_s: Vec<f64>,
    pub dv_s_dy: Vec<f64>,
    pub s: Vec<f64>,
}

impl SvdDecomposition {
    /// The invariant check upstream's `SVDEvaluator` constructor runs, so
    /// the eval path can assume shapes and monotone grids.
    pub fn validate(&self) -> Result<()> {
        let (nx, ny, r) = (self.nx, self.ny, self.rank);
        // Checked, not bare `nx * r`: on a 32-bit target the bare product
        // wraps, and a wrapped 0 would match an empty coefficient vector and
        // let an unusable decomposition pass this very check.
        let (nx_r, ny_r) = (
            crate::artifact::grid_len(nx, r)?,
            crate::artifact::grid_len(ny, r)?,
        );
        if nx < 2
            || ny < 2
            || r == 0
            || self.x_grid.len() != nx
            || self.y_grid.len() != ny
            || self.u.len() != nx_r
            || self.du_dx.len() != nx_r
            || self.v_s.len() != ny_r
            || self.dv_s_dy.len() != ny_r
        {
            return Err(Error::Value(
                "SVDEvaluator: SVDDecomposition has inconsistent dimensions".into(),
            ));
        }
        if self.x_grid.windows(2).any(|w| !(w[0] < w[1])) {
            return Err(Error::Value(
                "SVDEvaluator: x_grid must be strictly increasing".into(),
            ));
        }
        if self.y_grid.windows(2).any(|w| !(w[0] < w[1])) {
            return Err(Error::Value(
                "SVDEvaluator: y_grid must be strictly increasing".into(),
            ));
        }
        Ok(())
    }

    /// Upstream `SVDEvaluator::make_context`: locate both cells and build
    /// the Hermite weights once, so K property evaluations at the same
    /// (x, y) pay for it once.
    #[inline]
    pub fn make_context(&self, x: f64, y: f64) -> EvalContext {
        let i = locate(&self.x_grid, x);
        let j = locate(&self.y_grid, y);
        let hx = self.x_grid[i + 1] - self.x_grid[i];
        let hy = self.y_grid[j + 1] - self.y_grid[j];
        let tx = (x - self.x_grid[i]) / hx;
        let ty = (y - self.y_grid[j]) / hy;
        EvalContext {
            i,
            j,
            hx,
            hy,
            bx: hermite_basis(tx),
            by: hermite_basis(ty),
        }
    }

    /// Upstream `SVDEvaluator::eval_with_context`: the rank-r dot product of
    /// the two Hermite-interpolated mode vectors, optionally exponentiated.
    #[inline]
    pub fn eval_with_context(&self, c: &EvalContext) -> f64 {
        let r = self.rank;
        let (u0, u1) = (&self.u[c.i * r..], &self.u[(c.i + 1) * r..]);
        let (du0, du1) = (&self.du_dx[c.i * r..], &self.du_dx[(c.i + 1) * r..]);
        let (v0, v1) = (&self.v_s[c.j * r..], &self.v_s[(c.j + 1) * r..]);
        let (dv0, dv1) = (&self.dv_s_dy[c.j * r..], &self.dv_s_dy[(c.j + 1) * r..]);
        let bx10_hx = c.bx.h10 * c.hx;
        let bx11_hx = c.bx.h11 * c.hx;
        let by10_hy = c.by.h10 * c.hy;
        let by11_hy = c.by.h11 * c.hy;
        let mut acc = 0.0;
        for k in 0..r {
            let u_k = c.bx.h00 * u0[k] + bx10_hx * du0[k] + c.bx.h01 * u1[k] + bx11_hx * du1[k];
            let v_k = c.by.h00 * v0[k] + by10_hy * dv0[k] + c.by.h01 * v1[k] + by11_hy * dv1[k];
            acc += u_k * v_k;
        }
        match self.out_transform {
            OutputTransform::Exp => acc.exp(),
            OutputTransform::Identity => acc,
        }
    }

    /// `eval(x, y)`. No clamping: outside the grid the Hermite kernel
    /// extrapolates from the boundary cell, so callers gate on the region.
    #[inline]
    pub fn eval(&self, x: f64, y: f64) -> f64 {
        self.eval_with_context(&self.make_context(x, y))
    }
}

/// Upstream `SVDEvalContext`, shared across the property evaluators of one
/// region (they all share that region's grids).
#[derive(Clone, Copy, Debug)]
pub struct EvalContext {
    pub i: usize,
    pub j: usize,
    pub hx: f64,
    pub hy: f64,
    pub bx: HermiteBasis,
    pub by: HermiteBasis,
}

/// Upstream `SVDEvaluator::locate`: `i` with `grid[i] <= x <= grid[i+1]`,
/// clamped to `[0, n-2]`. Closed at both ends, so `x == grid.last()` gives
/// `n-2` and the kernel sees `t = 1` on the last cell.
#[inline]
fn locate(g: &[f64], x: f64) -> usize {
    let n = g.len();
    if x <= g[0] {
        return 0;
    }
    if x >= g[n - 1] {
        return n - 2;
    }
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if g[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}
