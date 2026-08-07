//! Bicubic evaluation — port of upstream `TabularDataSet::build_coeffs`
//! (`src/Backends/Tabular/TabularBackends.cpp`) and
//! `BicubicBackend::evaluate_single_phase`
//! (`src/Backends/Tabular/BicubicBackend.cpp`).
//!
//! Each cell carries 16 alpha coefficients per tabulated property, obtained
//! from the classic bicubic inverse matrix applied to the cell's corner
//! values, first derivatives and cross derivatives — all scaled into the
//! unit square (xhat, yhat). Evaluation is Horner in yhat then xhat.
//!
//! Upstream stores `Ainv` as a flat array loaded into an Eigen matrix
//! (COLUMN-major) and then multiplies by `Ainv.transpose()`, which is the
//! same as reading the flat array ROW-major — that is what `AINV` below is.

// Cell loops keep upstream's explicit (i, j) index form for reviewability.
#![allow(clippy::needless_range_loop)]

use crate::tables::{GriddedTable, valid_number};
use rustprop_core::params::Param;
use rustprop_core::{Error, Result};

/// Upstream `Ainv_data`, indexed row-major (see the module note).
#[rustfmt::skip]
const AINV: [f64; 256] = [
    1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    -3., 3., 0., 0., -2., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    2., -2., 0., 0., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
    0., 0., 0., 0., 0., 0., 0., 0., -3., 3., 0., 0., -2., -1., 0., 0.,
    0., 0., 0., 0., 0., 0., 0., 0., 2., -2., 0., 0., 1., 1., 0., 0.,
    -3., 0., 3., 0., 0., 0., 0., 0., -2., 0., -1., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., -3., 0., 3., 0., 0., 0., 0., 0., -2., 0., -1., 0.,
    9., -9., -9., 9., 6., 3., -6., -3., 6., -6., 3., -3., 4., 2., 2., 1.,
    -6., 6., 6., -6., -3., -3., 3., 3., -4., 4., -2., 2., -2., -2., -1., -1.,
    2., 0., -2., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 2., 0., -2., 0., 0., 0., 0., 0., 1., 0., 1., 0.,
    -6., 6., 6., -6., -4., -2., 4., 2., -3., 3., -3., 3., -2., -1., -2., -1.,
    4., -4., -4., 4., 2., 2., -2., -2., 2., -2., 2., -2., 1., 1., 1., 1.,
];

/// The six properties upstream builds coefficients for, in its order.
const PARAM_LIST: [Param; 6] = [
    Param::Dmolar,
    Param::T,
    Param::Smolar,
    Param::Hmolar,
    Param::P,
    Param::Umolar,
];

fn slot(param: Param) -> Option<usize> {
    PARAM_LIST.iter().position(|p| *p == param)
}

/// One cell's bicubic coefficients (upstream `CellCoeffs`).
#[derive(Clone, Default)]
pub struct CellCoeffs {
    /// 16 alpha coefficients per property slot (empty when not built).
    alpha: [Vec<f64>; 6],
    pub dx_dxhat: f64,
    pub dy_dyhat: f64,
    valid: bool,
    /// Upstream `set_alternate`: an invalid cell points at a good neighbour.
    alternate: Option<(usize, usize)>,
}

impl CellCoeffs {
    pub fn valid(&self) -> bool {
        self.valid
    }
    pub fn alternate(&self) -> Option<(usize, usize)> {
        self.alternate
    }
    pub fn get(&self, param: Param) -> Result<&[f64]> {
        let s = slot(param)
            .ok_or_else(|| Error::Value(format!("invalid key [{}]", param.short_name())))?;
        Ok(&self.alpha[s])
    }
}

/// The coefficient grid for one table (upstream's
/// `std::vector<std::vector<CellCoeffs>>`), sized `(Nx-1) x (Ny-1)`.
pub struct CellCoeffGrid {
    pub nx_cells: usize,
    pub ny_cells: usize,
    cells: Vec<Vec<CellCoeffs>>,
}

impl CellCoeffGrid {
    pub fn cell(&self, i: usize, j: usize) -> &CellCoeffs {
        &self.cells[i][j]
    }

    /// `TabularDataSet::build_coeffs`: build the 16 coefficients per property
    /// for every cell whose four corners are valid, then remap invalid cells
    /// onto a good neighbour in upstream's fixed offset order.
    pub fn build(table: &GriddedTable) -> Self {
        let (nxc, nyc) = (table.nx - 1, table.ny - 1);
        let mut cells = vec![vec![CellCoeffs::default(); nyc]; nxc];

        for param in PARAM_LIST {
            // Tables never carry coefficients for their own axes.
            if param == table.kind.xkey() || param == table.kind.ykey() {
                continue;
            }
            let g = table.grid_for(param).expect("param in list is tabulated");
            let s = slot(param).expect("param in list");
            for i in 0..nxc {
                for j in 0..nyc {
                    if valid_number(g.val[i][j])
                        && valid_number(g.val[i + 1][j])
                        && valid_number(g.val[i][j + 1])
                        && valid_number(g.val[i + 1][j + 1])
                    {
                        let dx_dxhat = table.xvec[i + 1] - table.xvec[i];
                        let dy_dyhat = table.yvec[j + 1] - table.yvec[j];
                        let f = [
                            g.val[i][j],
                            g.val[i + 1][j],
                            g.val[i][j + 1],
                            g.val[i + 1][j + 1],
                            g.dx[i][j] * dx_dxhat,
                            g.dx[i + 1][j] * dx_dxhat,
                            g.dx[i][j + 1] * dx_dxhat,
                            g.dx[i + 1][j + 1] * dx_dxhat,
                            g.dy[i][j] * dy_dyhat,
                            g.dy[i + 1][j] * dy_dyhat,
                            g.dy[i][j + 1] * dy_dyhat,
                            g.dy[i + 1][j + 1] * dy_dyhat,
                            g.dxy[i][j] * dy_dyhat * dx_dxhat,
                            g.dxy[i + 1][j] * dy_dyhat * dx_dxhat,
                            g.dxy[i][j + 1] * dy_dyhat * dx_dxhat,
                            g.dxy[i + 1][j + 1] * dy_dyhat * dx_dxhat,
                        ];
                        // alpha = Ainv(row-major) * F
                        let mut alpha = vec![0.0; 16];
                        for (r, a) in alpha.iter_mut().enumerate() {
                            let mut acc = 0.0;
                            for c in 0..16 {
                                acc += AINV[r * 16 + c] * f[c];
                            }
                            *a = acc;
                        }
                        let cell = &mut cells[i][j];
                        cell.alpha[s] = alpha;
                        cell.dx_dxhat = dx_dxhat;
                        cell.dy_dyhat = dy_dyhat;
                        cell.valid = true;
                    } else {
                        cells[i][j].valid = false;
                    }
                }
            }
        }

        // Remap invalid cells onto a good neighbour — same offsets and the
        // same strict interior bound as `make_good_neighbors`.
        const XOFF: [isize; 8] = [-1, 1, 0, 0, -1, 1, 1, -1];
        const YOFF: [isize; 8] = [0, 0, 1, -1, -1, -1, 1, 1];
        for i in 0..nxc {
            for j in 0..nyc {
                if !cells[i][j].valid {
                    for k in 0..8 {
                        let iplus = i as isize + XOFF[k];
                        let jplus = j as isize + YOFF[k];
                        if iplus > 0
                            && (iplus as usize) < nxc
                            && jplus > 0
                            && (jplus as usize) < nyc
                            && cells[iplus as usize][jplus as usize].valid
                        {
                            cells[i][j].alternate = Some((iplus as usize, jplus as usize));
                            break;
                        }
                    }
                }
            }
        }

        CellCoeffGrid {
            nx_cells: nxc,
            ny_cells: nyc,
            cells,
        }
    }
}

/// `BicubicBackend::evaluate_single_phase`: Horner in yhat, then in xhat,
/// over the cell's 16 coefficients. An invalid cell raises upstream's
/// verbatim error rather than dereferencing empty coefficients (#1950).
pub fn evaluate_single_phase(
    table: &GriddedTable,
    coeffs: &CellCoeffGrid,
    output: Param,
    x: f64,
    y: f64,
    i: usize,
    j: usize,
) -> Result<f64> {
    let cell = coeffs.cell(i, j);
    if !cell.valid() {
        return Err(Error::Value(format!(
            "evaluate_single_phase called on cell ({i}, {j}) with no bicubic coefficients (x={x}, y={y})"
        )));
    }
    let alpha = cell.get(output)?;

    // Normalized position in the unit square
    let xhat = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    let yhat = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);

    let b0 = ((alpha[3 * 4] * yhat + alpha[2 * 4]) * yhat + alpha[4]) * yhat + alpha[0];
    let b1 = ((alpha[3 * 4 + 1] * yhat + alpha[2 * 4 + 1]) * yhat + alpha[5]) * yhat + alpha[1];
    let b2 = ((alpha[3 * 4 + 2] * yhat + alpha[2 * 4 + 2]) * yhat + alpha[6]) * yhat + alpha[2];
    let b3 = ((alpha[3 * 4 + 3] * yhat + alpha[2 * 4 + 3]) * yhat + alpha[7]) * yhat + alpha[3];

    Ok(((b3 * xhat + b2) * xhat + b1) * xhat + b0)
}

/// `BicubicBackend::evaluate_single_phase_transport`: the same bilinear
/// interpolation TTSE uses, but WITHOUT its in-range and four-valid-corner
/// guards — upstream notes "By definition i,i+1,j,j+1 are all in range and
/// valid" because a bicubic cell only exists when its corners are valid.
pub fn evaluate_single_phase_transport(
    table: &GriddedTable,
    output: Param,
    x: f64,
    y: f64,
    i: usize,
    j: usize,
) -> Result<f64> {
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
