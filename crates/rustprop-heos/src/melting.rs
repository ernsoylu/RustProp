//! Melting-line evaluation (upstream `MeltingLineVariables::evaluate` in
//! `src/Backends/Helmholtz/Fluids/Ancillaries.cpp`): `p(T)` by segment
//! lookup on the closed T-range, `T(p)` by the Simon closed form or a Brent
//! solve per segment on its p-range, and the aggregate limits computed on
//! demand (upstream `set_limits` fills them post-parse; nothing is stored).

use rustprop_core::fluid::{MeltingLine, MeltingLineKind, PolyMeltPart};
use rustprop_core::{Error, Result};

fn in_closed(a: f64, b: f64, x: f64) -> bool {
    x >= a.min(b) && x <= a.max(b)
}

/// One polynomial segment's `p(T)` — `evaluate` of the InTr/InTheta segment
/// classes; `theta` selects the `(T/T_0 - 1)^t` form over `((T/T_0)^t - 1)`.
fn poly_part_eval(part: &PolyMeltPart, t: f64, theta: bool) -> f64 {
    let mut summer = 0.0;
    for i in 0..part.a.len() {
        summer += if theta {
            part.a[i] * (t / part.t_0 - 1.0).powf(part.t[i])
        } else {
            part.a[i] * ((t / part.t_0).powf(part.t[i]) - 1.0)
        };
    }
    part.p_0 * (1.0 + summer)
}

/// A segment's melting pressure at its own T-limit (upstream `set_limits`).
fn part_p_limits(ml: &MeltingLine, idx: usize) -> (f64, f64) {
    match &ml.kind {
        MeltingLineKind::Simon { parts } => {
            let p = &parts[idx];
            (
                p.p_0 + p.a * ((p.t_min / p.t_0).powf(p.c) - 1.0),
                p.p_0 + p.a * ((p.t_max / p.t_0).powf(p.c) - 1.0),
            )
        }
        MeltingLineKind::PolynomialInTr { parts } => {
            let p = &parts[idx];
            (
                poly_part_eval(p, p.t_min, false),
                poly_part_eval(p, p.t_max, false),
            )
        }
        MeltingLineKind::PolynomialInTheta { parts } => {
            let p = &parts[idx];
            (
                poly_part_eval(p, p.t_min, true),
                poly_part_eval(p, p.t_max, true),
            )
        }
    }
}

fn part_count(ml: &MeltingLine) -> usize {
    match &ml.kind {
        MeltingLineKind::Simon { parts } => parts.len(),
        MeltingLineKind::PolynomialInTr { parts } => parts.len(),
        MeltingLineKind::PolynomialInTheta { parts } => parts.len(),
    }
}

/// Upstream `get_parts_pranges`: each part's (p_min, p_max).
pub(crate) fn part_pranges(ml: &MeltingLine) -> Vec<(f64, f64)> {
    (0..part_count(ml)).map(|i| part_p_limits(ml, i)).collect()
}

/// Upstream `pmin` after `set_limits`: the first part's `p(T_min)`.
pub fn p_min(ml: &MeltingLine) -> f64 {
    part_p_limits(ml, 0).0
}

/// Upstream `pmax` after `set_limits`: the last part's `p(T_max)`.
pub fn p_max(ml: &MeltingLine) -> f64 {
    part_p_limits(ml, part_count(ml) - 1).1
}

/// Upstream `evaluate(iP, iT, T)`.
pub fn p_of_t(ml: &MeltingLine, t: f64) -> Result<f64> {
    match &ml.kind {
        MeltingLineKind::Simon { parts } => {
            for part in *parts {
                if in_closed(part.t_min, part.t_max, t) {
                    return Ok(part.p_0 + part.a * ((t / part.t_0).powf(part.c) - 1.0));
                }
            }
            Err(Error::Value(
                "unable to calculate melting line (p,T) for Simon curve".into(),
            ))
        }
        MeltingLineKind::PolynomialInTr { parts } => {
            for part in *parts {
                if in_closed(part.t_min, part.t_max, t) {
                    return Ok(poly_part_eval(part, t, false));
                }
            }
            Err(Error::Value(
                "unable to calculate melting line (p,T) for polynomial_in_Tr curve".into(),
            ))
        }
        MeltingLineKind::PolynomialInTheta { parts } => {
            for part in *parts {
                if in_closed(part.t_min, part.t_max, t) {
                    return Ok(poly_part_eval(part, t, true));
                }
            }
            Err(Error::Value(
                "unable to calculate melting line (p,T) for polynomial_in_Theta curve".into(),
            ))
        }
    }
}

/// Upstream `evaluate(iT, iP, p)`.
pub fn t_of_p(ml: &MeltingLine, p: f64) -> Result<f64> {
    match &ml.kind {
        MeltingLineKind::Simon { parts } => {
            for part in *parts {
                let t = ((p - part.p_0) / part.a + 1.0).powf(1.0 / part.c) * part.t_0;
                if t >= part.t_0 && t <= part.t_max {
                    return Ok(t);
                }
            }
            Err(Error::Value(format!(
                "unable to calculate melting line T(p) for Simon curve for p={p}; bounds are {},{} Pa",
                p_min(ml),
                p_max(ml)
            )))
        }
        MeltingLineKind::PolynomialInTr { parts } => {
            for (i, part) in parts.iter().enumerate() {
                let (pmin_i, pmax_i) = part_p_limits(ml, i);
                if in_closed(pmin_i, pmax_i, p) {
                    return crate::solvers::brent(
                        |t| p - poly_part_eval(part, t, false),
                        part.t_min,
                        part.t_max,
                        f64::EPSILON,
                        1e-12,
                        100,
                    );
                }
            }
            Err(Error::Value(format!(
                "unable to calculate melting line T(p) for polynomial_in_Theta curve for p={p}; bounds are {},{} Pa",
                p_min(ml),
                p_max(ml)
            )))
        }
        MeltingLineKind::PolynomialInTheta { parts } => {
            for (i, part) in parts.iter().enumerate() {
                let (pmin_i, pmax_i) = part_p_limits(ml, i);
                if in_closed(pmin_i, pmax_i, p) {
                    return crate::solvers::brent(
                        |t| p - poly_part_eval(part, t, true),
                        part.t_min,
                        part.t_max,
                        f64::EPSILON,
                        1e-12,
                        100,
                    );
                }
            }
            Err(Error::Value(format!(
                "unable to calculate melting line T(p) for polynomial_in_Theta curve for p={p}; bounds are {},{} Pa",
                p_min(ml),
                p_max(ml)
            )))
        }
    }
}
