//! Port of upstream `include/CoolProp/sbtl/SVDSurface.h`: one input pair's
//! property lookup for one pure fluid.
//!
//! ```text
//! eval(prop, a, b)
//!   -> atlas.find_region(a, b)      region index, or a miss
//!   -> region.to_normalized(a, b)   (xi, eta) in the unit square
//!   -> decomposition.eval(eta, xi)  rank-r Hermite reconstruction
//! ```
//!
//! Note the axis convention, which is upstream's and easy to get backwards:
//! the SVD's x axis is the SECONDARY normalised coordinate (eta) and its y
//! axis is the PRIMARY one (xi).

use crate::region::RegionAtlas;
use crate::svd::{EvalContext, SvdDecomposition};
use rustprop_core::params::Param;
use rustprop_core::{Error, Result};

/// Upstream `SVDSurface`.
pub struct SvdSurface {
    pub fluid_name: String,
    /// Upstream's `input_pairs` enum value, carried verbatim from the
    /// artifact (this port has no equivalent enum of its own).
    pub input_pair: i32,
    properties: Vec<Param>,
    atlas: RegionAtlas,
    /// `decomps[region][property]`.
    decomps: Vec<Vec<SvdDecomposition>>,
}

/// `SVDSurface::ResolvedPoint`.
#[derive(Clone, Copy, Debug)]
pub struct ResolvedPoint {
    pub region_idx: usize,
    /// The value for the SVD's x axis — the normalised SECONDARY coordinate.
    pub svd_x: f64,
    /// The value for the SVD's y axis — the normalised PRIMARY coordinate.
    pub svd_y: f64,
}

impl SvdSurface {
    pub fn new(
        fluid_name: String,
        input_pair: i32,
        properties: Vec<Param>,
        atlas: RegionAtlas,
        decomps: Vec<Vec<SvdDecomposition>>,
    ) -> Result<Self> {
        // `seal()`: every region must carry every property.
        if decomps.len() != atlas.len() {
            return Err(Error::Value(format!(
                "SVDSurface::seal: {} regions but {} decomposition rows",
                atlas.len(),
                decomps.len()
            )));
        }
        for (i, row) in decomps.iter().enumerate() {
            if row.len() != properties.len() {
                return Err(Error::Value(format!(
                    "SVDSurface::seal: region {i} has {} decompositions, expected {}",
                    row.len(),
                    properties.len()
                )));
            }
            for d in row {
                d.validate()?;
            }
            // `eval_with_region_multi` builds ONE context from `row[0]` and
            // reuses it across the row, which is only sound if every property
            // in a region really does share that region's grids. The artifact
            // stores a grid per region x property, so nothing upstream of
            // here enforces it — check it once at load rather than indexing
            // one decomposition's coefficients with another's cell index.
            if let Some((first, rest)) = row.split_first() {
                for (p, d) in rest.iter().enumerate() {
                    if d.nx != first.nx
                        || d.ny != first.ny
                        || d.x_grid != first.x_grid
                        || d.y_grid != first.y_grid
                    {
                        return Err(Error::Value(format!(
                            "SVDSurface::seal: region {i} property {} is on a {}x{} grid but \
                             property 0 is on {}x{}; every property in a region must share \
                             that region's grids",
                            p + 1,
                            d.nx,
                            d.ny,
                            first.nx,
                            first.ny
                        )));
                    }
                }
            }
        }
        Ok(SvdSurface {
            fluid_name,
            input_pair,
            properties,
            atlas,
            decomps,
        })
    }

    pub fn properties(&self) -> &[Param] {
        &self.properties
    }

    pub fn atlas(&self) -> &RegionAtlas {
        &self.atlas
    }

    pub fn region_count(&self) -> usize {
        self.atlas.len()
    }

    pub fn contains_property(&self, prop: Param) -> bool {
        self.properties.contains(&prop)
    }

    fn property_index(&self, prop: Param) -> Result<usize> {
        self.properties
            .iter()
            .position(|p| *p == prop)
            .ok_or_else(|| {
                Error::Value(format!(
                    "SVDSurface: property {} is not tabulated on this surface",
                    prop.short_name()
                ))
            })
    }

    /// The decomposition row for a region. `region_idx` reaches the public
    /// `eval_with_region*` entry points from the caller, so it is checked
    /// rather than indexed — these return `Result` and should use it.
    fn region_row(&self, region_idx: usize) -> Result<&[SvdDecomposition]> {
        self.decomps
            .get(region_idx)
            .map(Vec::as_slice)
            .ok_or_else(|| {
                Error::Value(format!(
                    "SVDSurface: region index {region_idx} is out of range; this surface has {} \
                 regions",
                    self.decomps.len()
                ))
            })
    }

    /// `resolve(a, b)` — region dispatch plus normalisation, with the
    /// eta/xi -> x/y swap upstream applies.
    pub fn resolve(&self, a: f64, b: f64) -> Option<ResolvedPoint> {
        let region_idx = self.atlas.find_region(a, b)?;
        let (xi, eta) = self.atlas.region(region_idx).to_normalized(a, b);
        Some(ResolvedPoint {
            region_idx,
            svd_x: eta,
            svd_y: xi,
        })
    }

    /// `eval(prop, a, b)`. A point outside every region yields NaN, as
    /// upstream does.
    pub fn eval(&self, prop: Param, a: f64, b: f64) -> Result<f64> {
        match self.resolve(a, b) {
            None => Ok(f64::NAN),
            Some(r) => self.eval_with_region(prop, r.region_idx, r.svd_x, r.svd_y),
        }
    }

    /// `eval_with_region`: the fast path after `resolve`.
    pub fn eval_with_region(
        &self,
        prop: Param,
        region_idx: usize,
        svd_x: f64,
        svd_y: f64,
    ) -> Result<f64> {
        let p = self.property_index(prop)?;
        Ok(self.region_row(region_idx)?[p].eval(svd_x, svd_y))
    }

    /// `eval_with_region_multi`: every property in one region shares that
    /// region's grids, so locate + basis setup is done once and amortised.
    pub fn eval_with_region_multi(
        &self,
        region_idx: usize,
        svd_x: f64,
        svd_y: f64,
        props: &[Param],
        out: &mut [f64],
    ) -> Result<()> {
        if out.len() < props.len() {
            return Err(Error::Value(
                "SVDSurface::eval_with_region_multi: output slice too short".into(),
            ));
        }
        let row = self.region_row(region_idx)?;
        let ctx: Option<EvalContext> = row.first().map(|d| d.make_context(svd_x, svd_y));
        for (n, prop) in props.iter().enumerate() {
            let p = self.property_index(*prop)?;
            let d = &row[p];
            out[n] = match ctx {
                Some(ref c) => d.eval_with_context(c),
                None => d.eval(svd_x, svd_y),
            };
        }
        Ok(())
    }
}
