//! TTSE evaluation — port of upstream `TTSEBackend` (`src/Backends/Tabular/
//! TTSEBackend.cpp`) plus the cell-search helpers it calls in
//! `SinglePhaseGriddedTableData` (`bisect_vector`,
//! `find_native_nearest_neighbor`, `find_native_nearest_good_neighbor`).
//!
//! TTSE is a second-order Taylor expansion about the NEAREST grid node:
//!
//! ```text
//! z(x, y) = z + dx*z_x + dy*z_y + dx^2/2*z_xx + dy^2/2*z_yy + dx*dy*z_xy
//! ```
//!
//! Transport properties carry no stored derivatives, so they interpolate
//! bilinearly over the surrounding cell instead — with upstream's guard that
//! all four corners must be valid (tested on `smolar`, not on the property
//! being interpolated: upstream's own choice, reproduced).

// Upstream shapes kept verbatim: the root-selection tests are written as
// negated comparisons (`!(|dx2| < spacing)`), and the derivative evaluator
// carries upstream's full argument list.
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::too_many_arguments)]

use crate::tables::{GriddedTable, PropGrid, valid_number};
use rustprop_core::params::Param;
use rustprop_core::{Error, Result};

/// Upstream `bisect_vector`: index `i` such that `val` lies in
/// `[vec[i], vec[i+1]]`, skipping invalid entries at the ends and in the
/// middle exactly as upstream does.
pub fn bisect_vector(vec: &[f64], val: f64) -> Result<usize> {
    let n = vec.len();
    let (mut l, mut r) = (0usize, n - 1);
    while !valid_number(vec[r]) {
        if r == 1 {
            return Err(Error::Value(
                "All the values in bisection vector are invalid".into(),
            ));
        }
        r -= 1;
    }
    while !valid_number(vec[l]) {
        if l == n - 1 {
            return Err(Error::Value(
                "All the values in bisection vector are invalid".into(),
            ));
        }
        l += 1;
    }
    let mut m = (l + r) / 2;
    let mut r_l = vec[l] - val;
    let mut r_r = vec[r] - val;
    while r - l > 1 {
        if !valid_number(vec[m]) {
            let (mut mr, mut ml) = (m, m);
            while !valid_number(vec[mr]) {
                if mr == n - 1 {
                    return Err(Error::Value(
                        "All the values in bisection vector are invalid".into(),
                    ));
                }
                mr += 1;
            }
            while !valid_number(vec[ml]) {
                if ml == 1 {
                    return Err(Error::Value(
                        "All the values in bisection vector are invalid".into(),
                    ));
                }
                ml -= 1;
            }
            let r_ml = vec[ml] - val;
            let r_mr = vec[mr] - val;
            if r_r * r_ml > 0.0 && r_l * r_ml < 0.0 {
                r = ml;
                r_r = vec[ml] - val;
            } else if r_r * r_mr < 0.0 && r_l * r_mr > 0.0 {
                l = mr;
                r_l = vec[mr] - val;
            } else {
                return Err(Error::Value(format!(
                    "Unable to bisect segmented vector; neither chunk contains the solution val:{val} left:({},{}) right:({},{})",
                    vec[l], vec[ml], vec[mr], vec[r]
                )));
            }
            m = (l + r) / 2;
        } else {
            let r_m = vec[m] - val;
            if r_r * r_m > 0.0 && r_l * r_m < 0.0 {
                r = m;
                r_r = vec[r] - val;
            } else {
                l = m;
                r_l = vec[l] - val;
            }
            m = (l + r) / 2;
        }
    }
    Ok(l)
}

impl GriddedTable {
    /// Upstream `is_in_bounds` (with its `e` slack on every side).
    pub fn in_bounds(&self, x: f64, y: f64) -> bool {
        let e = 0.0;
        x >= self.xmin - e && x <= self.xmax + e && y >= self.ymin - e && y <= self.ymax + e
    }

    /// `find_native_nearest_neighbor`: bisect, then round to whichever
    /// bracketing node is nearer — arithmetic midpoint on a linear axis,
    /// GEOMETRIC midpoint on a log axis.
    pub fn find_native_nearest_neighbor(&self, x: f64, y: f64) -> Result<(usize, usize)> {
        let mut i = bisect_vector(&self.xvec, x)?;
        if i != self.nx - 1 {
            // x is linear for both ported layouts (logx = false upstream).
            if x > (self.xvec[i] + self.xvec[i + 1]) / 2.0 {
                i += 1;
            }
        }
        let mut j = bisect_vector(&self.yvec, y)?;
        if j != self.ny - 1 {
            // y is log spaced for both layouts.
            if y > (self.yvec[j] * self.yvec[j + 1]).sqrt() {
                j += 1;
            }
        }
        Ok((i, j))
    }

    /// `find_native_nearest_good_neighbor`: the nearest node, replaced by its
    /// cached good neighbour when it is a hole.
    pub fn find_native_nearest_good_neighbor(&self, x: f64, y: f64) -> Result<(usize, usize)> {
        let (mut i, mut j) = self.find_native_nearest_neighbor(x, y)?;
        if !valid_number(self.t.val[i][j]) {
            let (inew, jnew) = (self.nearest_i[i][j], self.nearest_j[i][j]);
            i = inew;
            j = jnew;
        }
        Ok((i, j))
    }

    /// `find_native_nearest_good_cell`: the lower-left corner of the cell
    /// containing (x, y) — plain bisection, no rounding (bicubic needs the
    /// cell, not the nearest node).
    pub fn find_native_nearest_good_cell(&self, x: f64, y: f64) -> Result<(usize, usize)> {
        Ok((bisect_vector(&self.xvec, x)?, bisect_vector(&self.yvec, y)?))
    }

    /// The stored grid for one output (upstream `get(key)`).
    pub fn grid_for(&self, key: Param) -> Result<&PropGrid> {
        Ok(match key {
            Param::Dmolar => &self.rhomolar,
            Param::T => &self.t,
            Param::Umolar => &self.umolar,
            Param::Hmolar => &self.hmolar,
            Param::Smolar => &self.smolar,
            Param::P => &self.p,
            other => {
                return Err(Error::Value(format!(
                    "invalid key [{}]",
                    other.short_name()
                )));
            }
        })
    }

    /// Transport grids, which have no derivatives.
    fn transport_for(&self, key: Param) -> Result<&Vec<Vec<f64>>> {
        Ok(match key {
            Param::Viscosity => &self.visc,
            Param::Conductivity => &self.cond,
            other => {
                return Err(Error::Value(format!(
                    "invalid key [{}]",
                    other.short_name()
                )));
            }
        })
    }
}

/// `TTSEBackend::evaluate_single_phase`: the second-order Taylor expansion
/// about node (i, j).
pub fn evaluate_single_phase(
    table: &GriddedTable,
    output: Param,
    x: f64,
    y: f64,
    i: usize,
    j: usize,
) -> Result<f64> {
    let g = table.grid_for(output)?;
    let deltax = x - table.xvec[i];
    let deltay = y - table.yvec[j];
    Ok(g.val[i][j]
        + deltax * g.dx[i][j]
        + deltay * g.dy[i][j]
        + 0.5 * deltax * deltax * g.dxx[i][j]
        + 0.5 * deltay * deltay * g.dyy[i][j]
        + deltay * deltax * g.dxy[i][j])
}

/// `TTSEBackend::evaluate_single_phase_derivative` for the first derivatives
/// (upstream throws `NotImplementedError` beyond first order). `nx`/`ny`
/// select d/dx or d/dy; the native-key short circuits come first, exactly as
/// upstream has them (twice — before and after `connect_pointers`).
pub fn evaluate_single_phase_derivative(
    table: &GriddedTable,
    output: Param,
    x: f64,
    y: f64,
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
) -> Result<f64> {
    let (xkey, ykey) = (table.kind.xkey(), table.kind.ykey());
    if nx == 1 && ny == 0 {
        if output == xkey {
            return Ok(1.0);
        }
        if output == ykey {
            return Ok(0.0);
        }
    } else if ny == 1 && nx == 0 {
        if output == ykey {
            return Ok(1.0);
        }
        if output == xkey {
            return Ok(0.0);
        }
    }
    let g = table.grid_for(output)?;
    let deltax = x - table.xvec[i];
    let deltay = y - table.yvec[j];
    if nx == 1 && ny == 0 {
        Ok(g.dx[i][j] + deltax * g.dxx[i][j] + deltay * g.dxy[i][j])
    } else if ny == 1 && nx == 0 {
        Ok(g.dy[i][j] + deltay * g.dyy[i][j] + deltax * g.dxy[i][j])
    } else {
        Err(Error::NotImplemented(
            "only first derivatives currently supported".into(),
        ))
    }
}

/// `TTSEBackend::evaluate_single_phase_transport`: bilinear interpolation
/// over the cell, with upstream's four-valid-corner guard (checked on
/// `smolar`, whatever the requested output).
pub fn evaluate_single_phase_transport(
    table: &GriddedTable,
    output: Param,
    x: f64,
    y: f64,
    i: usize,
    j: usize,
) -> Result<f64> {
    let in_bounds = i < table.xvec.len() - 1 && j < table.yvec.len() - 1;
    if !in_bounds {
        return Err(Error::Value(
            "Cell to TTSEBackend::evaluate_single_phase_transport is not valid".into(),
        ));
    }
    let is_valid = valid_number(table.smolar.val[i][j])
        && valid_number(table.smolar.val[i + 1][j])
        && valid_number(table.smolar.val[i][j + 1])
        && valid_number(table.smolar.val[i + 1][j + 1]);
    if !is_valid {
        return Err(Error::Value(
            "Cell to TTSEBackend::evaluate_single_phase_transport must have four valid corners for now".into(),
        ));
    }
    let f = table.transport_for(output)?;
    let (x1, x2) = (table.xvec[i], table.xvec[i + 1]);
    let (y1, y2) = (table.yvec[j], table.yvec[j + 1]);
    let (f11, f12) = (f[i][j], f[i][j + 1]);
    let (f21, f22) = (f[i + 1][j], f[i + 1][j + 1]);
    Ok(1.0 / ((x2 - x1) * (y2 - y1))
        * (f11 * (x2 - x) * (y2 - y)
            + f21 * (x - x1) * (y2 - y)
            + f12 * (x2 - x) * (y - y1)
            + f22 * (x - x1) * (y - y1)))
}

/// `TTSEBackend::invert_single_phase_x`: solve the quadratic in `deltax` for
/// the x that produces the requested `output` value, with upstream's
/// spacing-based root selection.
pub fn invert_single_phase_x(
    table: &GriddedTable,
    output: Param,
    x: f64,
    y: f64,
    i: usize,
    j: usize,
) -> Result<f64> {
    let g = table.grid_for(output)?;
    let deltay = y - table.yvec[j];
    let a = 0.5 * g.dxx[i][j];
    let b = g.dx[i][j] + deltay * g.dxy[i][j];
    let c = g.val[i][j] - x + deltay * g.dy[i][j] + 0.5 * deltay * deltay * g.dyy[i][j];

    let deltax1 = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
    let deltax2 = (-b - (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);

    // x is linear for both ported layouts.
    let xspacing = table.xvec[1] - table.xvec[0];
    if deltax1.abs() < xspacing && !(deltax2.abs() < xspacing) {
        Ok(deltax1 + table.xvec[i])
    } else if deltax2.abs() < xspacing && !(deltax1.abs() < xspacing) {
        Ok(deltax2 + table.xvec[i])
    } else if deltax1.abs() < deltax2.abs() && deltax1.abs() < 10.0 * xspacing {
        Ok(deltax1 + table.xvec[i])
    } else {
        Err(Error::Value(format!(
            "Cannot find the x solution; xspacing: {xspacing} dx1: {deltax1} dx2: {deltax2}"
        )))
    }
}

/// `TTSEBackend::invert_single_phase_y`: the same quadratic in `deltay`.
/// Upstream's y axis is log spaced, so the root test compares RATIOS against
/// the grid ratio (and its `xvec[j]`/`yvec[j]` mix-up in the x-variant is
/// not reachable here).
pub fn invert_single_phase_y(
    table: &GriddedTable,
    output: Param,
    y: f64,
    x: f64,
    i: usize,
    j: usize,
) -> Result<f64> {
    let g = table.grid_for(output)?;
    let deltax = x - table.xvec[i];
    let a = 0.5 * g.dyy[i][j];
    let b = g.dy[i][j] + deltax * g.dxy[i][j];
    let c = g.val[i][j] - y + deltax * g.dx[i][j] + 0.5 * deltax * deltax * g.dxx[i][j];

    let deltay1 = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
    let deltay2 = (-b - (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);

    let yratio = table.yvec[1] / table.yvec[0];
    let yj = table.yvec[j];
    let yratio1 = (yj + deltay1) / yj;
    let yratio2 = (yj + deltay2) / yj;
    if yratio1 < yratio && yratio1 > 1.0 / yratio {
        Ok(deltay1 + table.yvec[j])
    } else if yratio2 < yratio && yratio2 > 1.0 / yratio {
        Ok(deltay2 + table.yvec[j])
    } else if yratio1 < yratio * 5.0 && yratio1 > 1.0 / yratio / 5.0 {
        Ok(deltay1 + table.yvec[j])
    } else {
        Err(Error::Value(format!(
            "Cannot find the y solution; yj: {yj} yratio: {yratio} yratio1: {yratio1} yratio2: {yratio2} a: {a} b^2-4*a*c {}",
            b * b - 4.0 * a * c
        )))
    }
}
