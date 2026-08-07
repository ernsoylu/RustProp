//! Port of upstream `include/CoolProp/region/` — the (a, b) plane geometry
//! an `SVDSurface` normalises through: `AxisTransform`, the boundary-curve
//! family, `Region` and `RegionAtlas`.
//!
//! Only the curve kinds that upstream's serializer can write are ported:
//! `ConstantCurve` (kind 0) and `CubicSplineCurve` (kind 1). A
//! `PiecewiseChebyshevCurve` (kind 2) is rejected loudly at load time — no
//! artifact this port has seen uses one, and guessing at its evaluation
//! would be worse than saying so.

// Upstream writes these guards as negated comparisons (`!(a_hi > a_lo)`,
// `!(w[0] < w[1])`) so a NaN bound fails the check instead of passing it.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use rustprop_core::{Error, Result};

/// Upstream `region::AxisScale`, in its enum order (the serialized value).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AxisScale {
    Linear = 0,
    Log = 1,
    Power = 2,
    PowerLo = 3,
}

impl AxisScale {
    pub fn from_index(i: u8) -> Result<Self> {
        Ok(match i {
            0 => AxisScale::Linear,
            1 => AxisScale::Log,
            2 => AxisScale::Power,
            3 => AxisScale::PowerLo,
            other => {
                return Err(Error::Value(format!("unknown AxisScale index {other}")));
            }
        })
    }
}

/// Upstream `region::AxisTransform`: the primary-axis map to `xi` in [0, 1].
#[derive(Clone, Copy, Debug)]
pub struct AxisTransform {
    pub scale: AxisScale,
    pub a_lo: f64,
    pub a_hi: f64,
    pub a_lo_t: f64,
    pub a_hi_t: f64,
    /// `1 / (a_hi_t - a_lo_t)`; for POWER/POWER_LO the LINEAR span
    /// `1 / (a_hi - a_lo)` instead, so the chain rule has one factor.
    pub inv_span_t: f64,
}

impl AxisTransform {
    /// `AxisTransform::make`.
    pub fn make(scale: AxisScale, a_lo: f64, a_hi: f64) -> Result<Self> {
        if !(a_hi > a_lo) {
            return Err(Error::Value("AxisTransform: a_hi must exceed a_lo".into()));
        }
        if scale == AxisScale::Log && !(a_lo > 0.0) {
            return Err(Error::Value("AxisTransform: LOG requires a_lo > 0".into()));
        }
        if matches!(scale, AxisScale::Power | AxisScale::PowerLo) {
            return Ok(AxisTransform {
                scale,
                a_lo,
                a_hi,
                a_lo_t: 0.0,
                a_hi_t: 0.0,
                inv_span_t: 1.0 / (a_hi - a_lo),
            });
        }
        let (lo_t, hi_t) = if scale == AxisScale::Log {
            (a_lo.ln(), a_hi.ln())
        } else {
            (a_lo, a_hi)
        };
        Ok(AxisTransform {
            scale,
            a_lo,
            a_hi,
            a_lo_t: lo_t,
            a_hi_t: hi_t,
            inv_span_t: 1.0 / (hi_t - lo_t),
        })
    }

    /// `a -> xi`. No clamping; out-of-range inputs map outside [0, 1].
    #[inline]
    pub fn forward(&self, a: f64) -> f64 {
        match self.scale {
            AxisScale::Power => 1.0 - ((self.a_hi - a) * self.inv_span_t).cbrt(),
            AxisScale::PowerLo => ((a - self.a_lo) * self.inv_span_t).cbrt(),
            AxisScale::Log => (a.ln() - self.a_lo_t) * self.inv_span_t,
            AxisScale::Linear => (a - self.a_lo_t) * self.inv_span_t,
        }
    }

    /// `xi -> a`.
    #[inline]
    pub fn inverse(&self, xi: f64) -> f64 {
        match self.scale {
            AxisScale::Power => {
                let t = 1.0 - xi;
                self.a_hi - (self.a_hi - self.a_lo) * t * t * t
            }
            AxisScale::PowerLo => self.a_lo + (self.a_hi - self.a_lo) * xi * xi * xi,
            AxisScale::Log => (self.a_lo_t + xi * (self.a_hi_t - self.a_lo_t)).exp(),
            AxisScale::Linear => self.a_lo_t + xi * (self.a_hi_t - self.a_lo_t),
        }
    }
}

/// The boundary curves upstream's serializer can write.
#[derive(Clone, Debug)]
pub enum BoundaryCurve {
    /// `ConstantCurve`: `b(a) = b` on `[a_lo, a_hi]`.
    Constant { a_lo: f64, a_hi: f64, b: f64 },
    /// `CubicSplineCurve`: a natural cubic spline stored as knots plus the
    /// second-derivative table `m`, with the tight bounds `build()` found
    /// by cubic-extremum root finding.
    CubicSpline(CubicSpline),
}

/// A natural cubic spline reconstructed from upstream's `State` snapshot —
/// no tridiagonal re-solve, exactly like `CubicSplineCurve::from_state`.
#[derive(Clone, Debug)]
pub struct CubicSpline {
    a: Vec<f64>,
    b: Vec<f64>,
    m: Vec<f64>,
    b_min: f64,
    b_max: f64,
    /// `kBuckets / (a.last - a.first)` for the O(1) indexed search.
    inv_step: f64,
    bucket_to_knot: Vec<u16>,
}

impl CubicSpline {
    const BUCKETS: usize = 256;

    /// `CubicSplineCurve::from_state`.
    pub fn from_state(
        a: Vec<f64>,
        b: Vec<f64>,
        m: Vec<f64>,
        b_min: f64,
        b_max: f64,
    ) -> Result<Self> {
        let n = a.len();
        if b.len() != n || m.len() != n {
            return Err(Error::Value(
                "CubicSplineCurve::from_state: a/b/M size mismatch".into(),
            ));
        }
        if n < 2 {
            return Err(Error::Value(
                "CubicSplineCurve::from_state: need at least 2 knots".into(),
            ));
        }
        if a.windows(2).any(|w| !(w[0] < w[1])) {
            return Err(Error::Value(
                "CubicSplineCurve::from_state: a must be strictly increasing".into(),
            ));
        }
        let mut s = CubicSpline {
            a,
            b,
            m,
            b_min,
            b_max,
            inv_step: 0.0,
            bucket_to_knot: vec![0; Self::BUCKETS],
        };
        s.build_bucket_table();
        Ok(s)
    }

    /// `build_bucket_table_`.
    fn build_bucket_table(&mut self) {
        let n = self.a.len();
        let a_lo = self.a[0];
        let a_hi = self.a[n - 1];
        let span = a_hi - a_lo;
        self.inv_step = if span > 0.0 {
            Self::BUCKETS as f64 / span
        } else {
            0.0
        };
        let mut i = 0usize;
        let max_seg = n - 2;
        for k in 0..Self::BUCKETS {
            let a_bucket = a_lo + (k as f64) * span / (Self::BUCKETS as f64);
            while i + 1 < n - 1 && self.a[i + 1] <= a_bucket {
                i += 1;
            }
            self.bucket_to_knot[k] = i.min(max_seg) as u16;
        }
    }

    /// `locate`: bucket hash then a short forward walk. Non-finite input
    /// returns the leftmost segment so the NaN propagates through eval.
    #[inline]
    fn locate(&self, a: f64) -> usize {
        let n = self.a.len();
        let a_lo = self.a[0];
        if !a.is_finite() || a <= a_lo {
            return 0;
        }
        if a >= self.a[n - 1] {
            return n - 2;
        }
        let k = ((a - a_lo) * self.inv_step) as isize;
        let k = k.clamp(0, Self::BUCKETS as isize - 1) as usize;
        let mut i = self.bucket_to_knot[k] as usize;
        while i + 1 < n - 1 && self.a[i + 1] <= a {
            i += 1;
        }
        i
    }

    #[inline]
    pub fn eval(&self, a: f64) -> f64 {
        let i = self.locate(a);
        let h = self.a[i + 1] - self.a[i];
        let t = (a - self.a[i]) / h;
        segment_eval(t, h, self.b[i], self.b[i + 1], self.m[i], self.m[i + 1])
    }

    #[inline]
    pub fn eval_da(&self, a: f64) -> f64 {
        let i = self.locate(a);
        let h = self.a[i + 1] - self.a[i];
        let t = (a - self.a[i]) / h;
        segment_eval_da(t, h, self.b[i], self.b[i + 1], self.m[i], self.m[i + 1])
    }
}

/// The natural-spline segment form upstream evaluates.
#[inline]
fn segment_eval(t: f64, h: f64, b_i: f64, b_ip1: f64, m_i: f64, m_ip1: f64) -> f64 {
    let s = 1.0 - t;
    let h2_6 = h * h / 6.0;
    s * b_i + t * b_ip1 + h2_6 * ((s * s * s - s) * m_i + (t * t * t - t) * m_ip1)
}

#[inline]
fn segment_eval_da(t: f64, h: f64, b_i: f64, b_ip1: f64, m_i: f64, m_ip1: f64) -> f64 {
    let s = 1.0 - t;
    (b_ip1 - b_i) / h + (h / 6.0) * ((1.0 - 3.0 * s * s) * m_i + (3.0 * t * t - 1.0) * m_ip1)
}

impl BoundaryCurve {
    #[inline]
    pub fn eval(&self, a: f64) -> f64 {
        match self {
            BoundaryCurve::Constant { b, .. } => *b,
            BoundaryCurve::CubicSpline(s) => s.eval(a),
        }
    }

    #[inline]
    pub fn eval_da(&self, a: f64) -> f64 {
        match self {
            BoundaryCurve::Constant { .. } => 0.0,
            BoundaryCurve::CubicSpline(s) => s.eval_da(a),
        }
    }

    /// Upstream's `eval_fast` is a surrogate only for the SuperAncillary
    /// curves, which never reach a serialized artifact; for these two kinds
    /// the default forwards to `eval`.
    #[inline]
    pub fn eval_fast(&self, a: f64) -> f64 {
        self.eval(a)
    }

    /// `(b_min, b_max)` over the build interval.
    pub fn bounds(&self) -> (f64, f64) {
        match self {
            BoundaryCurve::Constant { b, .. } => (*b, *b),
            BoundaryCurve::CubicSpline(s) => (s.b_min, s.b_max),
        }
    }

    pub fn a_range(&self) -> (f64, f64) {
        match self {
            BoundaryCurve::Constant { a_lo, a_hi, .. } => (*a_lo, *a_hi),
            BoundaryCurve::CubicSpline(s) => (s.a[0], s.a[s.a.len() - 1]),
        }
    }
}

/// Upstream `region::Region`.
#[derive(Clone, Debug)]
pub struct Region {
    pub primary: AxisTransform,
    pub secondary_scale: AxisScale,
    pub b_lo: BoundaryCurve,
    pub b_hi: BoundaryCurve,
    bbox: BBox,
}

#[derive(Clone, Copy, Debug)]
struct BBox {
    a_lo: f64,
    a_hi: f64,
    b_min: f64,
    b_max: f64,
}

impl Region {
    pub fn new(
        primary: AxisTransform,
        b_lo: BoundaryCurve,
        b_hi: BoundaryCurve,
        secondary_scale: AxisScale,
    ) -> Result<Self> {
        if !matches!(secondary_scale, AxisScale::Linear | AxisScale::Log) {
            return Err(Error::Value(
                "Region: secondary axis scale must be LINEAR or LOG".into(),
            ));
        }
        // `compute_bbox`: the primary axis bounds, plus b_lo's MIN and
        // b_hi's MAX (not the union of both curves' ranges).
        let bbox = BBox {
            a_lo: primary.a_lo,
            a_hi: primary.a_hi,
            b_min: b_lo.bounds().0,
            b_max: b_hi.bounds().1,
        };
        Ok(Region {
            primary,
            secondary_scale,
            b_lo,
            b_hi,
            bbox,
        })
    }

    #[inline]
    pub fn aabb_contains(&self, a: f64, b: f64) -> bool {
        a >= self.bbox.a_lo && a <= self.bbox.a_hi && b >= self.bbox.b_min && b <= self.bbox.b_max
    }

    /// Sign-only bracketing test; only call after `aabb_contains`.
    #[inline]
    pub fn curve_contains(&self, a: f64, b: f64) -> bool {
        b >= self.b_lo.eval_fast(a) && b <= self.b_hi.eval_fast(a)
    }

    /// `(a, b) -> (xi, eta)`. A degenerate boundary span (the curves
    /// pinching together) yields NaN eta rather than a silently-wrong zero.
    #[inline]
    pub fn to_normalized(&self, a: f64, b: f64) -> (f64, f64) {
        let xi = self.primary.forward(a);
        let b_lo_val = self.b_lo.eval(a);
        let b_hi_val = self.b_hi.eval(a);
        let log = self.secondary_scale == AxisScale::Log;
        let gb = if log { b.ln() } else { b };
        let g_lo = if log { b_lo_val.ln() } else { b_lo_val };
        let g_hi = if log { b_hi_val.ln() } else { b_hi_val };
        let span = g_hi - g_lo;
        let tol = f64::EPSILON * (1.0 + g_lo.abs() + g_hi.abs());
        let eta = if span.abs() <= tol {
            f64::NAN
        } else {
            (gb - g_lo) / span
        };
        (xi, eta)
    }

    /// `(xi, eta) -> (a, b)`.
    #[inline]
    pub fn from_normalized(&self, xi: f64, eta: f64) -> (f64, f64) {
        let a = self.primary.inverse(xi);
        let b_lo_val = self.b_lo.eval(a);
        let b_hi_val = self.b_hi.eval(a);
        let b = if self.secondary_scale == AxisScale::Log {
            let g_lo = b_lo_val.ln();
            let g_hi = b_hi_val.ln();
            (g_lo + eta * (g_hi - g_lo)).exp()
        } else {
            b_lo_val + eta * (b_hi_val - b_lo_val)
        };
        (a, b)
    }
}

/// Upstream `region::RegionAtlas`: bounding-box-first dispatch in
/// REGISTRATION ORDER — first match wins, which is part of the contract
/// because disjoint regions routinely have overlapping AABBs.
#[derive(Clone, Debug, Default)]
pub struct RegionAtlas {
    regions: Vec<Region>,
}

impl RegionAtlas {
    pub fn add(&mut self, region: Region) -> usize {
        self.regions.push(region);
        self.regions.len() - 1
    }

    pub fn len(&self) -> usize {
        self.regions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    pub fn region(&self, i: usize) -> &Region {
        &self.regions[i]
    }

    /// The first region whose AABB *and* curve envelope contain (a, b).
    pub fn find_region(&self, a: f64, b: f64) -> Option<usize> {
        self.regions
            .iter()
            .position(|r| r.aabb_contains(a, b) && r.curve_contains(a, b))
    }

    /// Every region whose curve envelope contains (a, b) — debug helper for
    /// asserting disjointness, not the hot path.
    pub fn find_all_curve_hits(&self, a: f64, b: f64) -> Vec<usize> {
        (0..self.regions.len())
            .filter(|&i| self.regions[i].curve_contains(a, b))
            .collect()
    }
}
