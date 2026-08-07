//! The `.svds` blob format: how ingested SVD coefficients reach a shipped
//! binary.
//!
//! Upstream caches its surfaces as zlib-compressed msgpack. This port does
//! NOT parse that at runtime — a msgpack decoder in a wasm binary is exactly
//! the kind of weight this project avoids, and the same rule already keeps
//! JSON out of `rustprop-data`. Instead `tools/rustprop-svdgen` converts an
//! upstream artifact once, offline, into the flat little-endian layout below,
//! which a consumer embeds with `include_bytes!` and reads with plain
//! `f64::from_le_bytes`.
//!
//! Every field is little-endian; `f64` fields are IEEE-754 bit patterns, so
//! the coefficients are bitwise upstream's.
//!
//! ```text
//! header   "RPSVDS01" (8 bytes) | input_pair i32 | n_props u32
//!          | n_regions u32 | reserved u32
//! props    n_props   x i32                (upstream parameter keys)
//! regions  n_regions x region
//! decomps  n_regions*n_props x decomp     (region-major)
//!
//! region   primary_scale u8 | secondary_scale u8 | pad u16
//!          | a_lo f64 | a_hi f64 | a_lo_t f64 | a_hi_t f64 | inv_span_t f64
//!          | curve b_lo | curve b_hi
//!
//! curve    kind u32
//!          kind 0 (constant): a_lo f64 | a_hi f64 | b f64
//!          kind 1 (spline):   n u32 | pad u32 | b_min f64 | b_max f64
//!                             | a[n] | b[n] | m[n]
//!
//! decomp   nx u32 | ny u32 | rank u32 | out_transform u8 | slope_source u8
//!          | pad u16
//!          | x_grid[nx] | y_grid[ny] | u[nx*rank] | du_dx[nx*rank]
//!          | v_s[ny*rank] | dv_s_dy[ny*rank] | s[rank]
//! ```

use crate::region::{AxisScale, AxisTransform, BoundaryCurve, CubicSpline, Region, RegionAtlas};
use crate::surface::SvdSurface;
use crate::svd::{OutputTransform, SlopeSource, SvdDecomposition};
use rustprop_core::params::Param;
use rustprop_core::{Error, Result};

/// The magic every `.svds` blob starts with.
pub const MAGIC: &[u8; 8] = b"RPSVDS01";

struct Reader<'a> {
    b: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    fn new(b: &'a [u8]) -> Self {
        Reader { b, at: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.at.checked_add(n).ok_or_else(|| short(self.at, n))?;
        if end > self.b.len() {
            return Err(short(self.at, n));
        }
        let s = &self.b[self.at..end];
        self.at = end;
        Ok(s)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }

    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(
            self.take(8)?.try_into().expect("8 bytes"),
        ))
    }

    fn f64s(&mut self, n: usize) -> Result<Vec<f64>> {
        let bytes = self.take(n.checked_mul(8).ok_or_else(|| short(self.at, n))?)?;
        Ok(bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().expect("8 bytes")))
            .collect())
    }
}

fn short(at: usize, want: usize) -> Error {
    Error::Value(format!(
        "svds artifact truncated: wanted {want} bytes at offset {at}"
    ))
}

fn read_curve(r: &mut Reader) -> Result<BoundaryCurve> {
    match r.u32()? {
        0 => Ok(BoundaryCurve::Constant {
            a_lo: r.f64()?,
            a_hi: r.f64()?,
            b: r.f64()?,
        }),
        1 => {
            let n = r.u32()? as usize;
            let _pad = r.u32()?;
            let b_min = r.f64()?;
            let b_max = r.f64()?;
            let a = r.f64s(n)?;
            let b = r.f64s(n)?;
            let m = r.f64s(n)?;
            Ok(BoundaryCurve::CubicSpline(CubicSpline::from_state(
                a, b, m, b_min, b_max,
            )?))
        }
        // Upstream's kind 2 is PiecewiseChebyshevCurve. No artifact this
        // port has seen uses one; refusing beats guessing at its basis.
        other => Err(Error::NotImplemented(format!(
            "svds artifact uses boundary-curve kind {other}, which is not ported"
        ))),
    }
}

fn read_region(r: &mut Reader) -> Result<Region> {
    let primary_scale = AxisScale::from_index(r.u8()?)?;
    let secondary_scale = AxisScale::from_index(r.u8()?)?;
    let _pad = r.take(2)?;
    let a_lo = r.f64()?;
    let a_hi = r.f64()?;
    let a_lo_t = r.f64()?;
    let a_hi_t = r.f64()?;
    let inv_span_t = r.f64()?;
    let b_lo = read_curve(r)?;
    let b_hi = read_curve(r)?;
    // The transformed bounds are carried rather than recomputed so the
    // normalisation uses upstream's exact bits.
    let primary = AxisTransform {
        scale: primary_scale,
        a_lo,
        a_hi,
        a_lo_t,
        a_hi_t,
        inv_span_t,
    };
    Region::new(primary, b_lo, b_hi, secondary_scale)
}

fn read_decomp(r: &mut Reader) -> Result<SvdDecomposition> {
    let nx = r.u32()? as usize;
    let ny = r.u32()? as usize;
    let rank = r.u32()? as usize;
    let out_transform = OutputTransform::from_index(r.u8()?)?;
    let slope_source = SlopeSource::from_index(r.u8()?)?;
    let _pad = r.take(2)?;
    let x_grid = r.f64s(nx)?;
    let y_grid = r.f64s(ny)?;
    let u = r.f64s(nx * rank)?;
    let du_dx = r.f64s(nx * rank)?;
    let v_s = r.f64s(ny * rank)?;
    let dv_s_dy = r.f64s(ny * rank)?;
    let s = r.f64s(rank)?;
    Ok(SvdDecomposition {
        nx,
        ny,
        rank,
        out_transform,
        slope_source,
        x_grid,
        y_grid,
        u,
        du_dx,
        v_s,
        dv_s_dy,
        s,
    })
}

/// Load a surface from a `.svds` blob.
pub fn load(fluid_name: &str, bytes: &[u8]) -> Result<SvdSurface> {
    let mut r = Reader::new(bytes);
    if r.take(8)? != MAGIC {
        return Err(Error::Value(
            "not an svds artifact (bad magic; expected RPSVDS01)".into(),
        ));
    }
    let input_pair = r.i32()?;
    let n_props = r.u32()? as usize;
    let n_regions = r.u32()? as usize;
    let _reserved = r.u32()?;

    let mut properties = Vec::with_capacity(n_props);
    for _ in 0..n_props {
        let key = r.i32()?;
        properties.push(param_from_index(key)?);
    }

    let mut atlas = RegionAtlas::default();
    for _ in 0..n_regions {
        atlas.add(read_region(&mut r)?);
    }

    let mut decomps = Vec::with_capacity(n_regions);
    for _ in 0..n_regions {
        let mut row = Vec::with_capacity(n_props);
        for _ in 0..n_props {
            row.push(read_decomp(&mut r)?);
        }
        decomps.push(row);
    }

    SvdSurface::new(
        fluid_name.to_string(),
        input_pair,
        properties,
        atlas,
        decomps,
    )
}

/// Upstream `parameters` values, for the keys its surfaces tabulate.
fn param_from_index(key: i32) -> Result<Param> {
    Ok(match key {
        19 => Param::T,
        20 => Param::P,
        40 => Param::Dmass,
        41 => Param::Hmass,
        42 => Param::Smass,
        46 => Param::Umass,
        52 => Param::Viscosity,
        53 => Param::Conductivity,
        56 => Param::SpeedSound,
        other => {
            return Err(Error::NotImplemented(format!(
                "svds artifact tabulates upstream parameter {other}, which this port does not map"
            )));
        }
    })
}
