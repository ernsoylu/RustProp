//! Regression tests for `bisect_vector` / `bisect_segmented_vector_slice` on
//! degenerate inputs.
//!
//! Both are `pub` in a `pub mod` and both compute `n - 1` on a length taken
//! from their argument. On an empty slice that underflows: debug builds trap
//! on the subtraction, and release builds — where `overflow-checks` is off —
//! wrap to `usize::MAX` and panic one line later on the index instead.
//!
//! No well-formed table reaches these states, so refusing them cannot change
//! a result that currently computes; it only replaces a panic with the `Err`
//! the signature already promises.

use rustprop_tabular::ttse::{bisect_segmented_vector_slice, bisect_vector};

#[test]
fn empty_vector_does_not_panic() {
    let e = bisect_vector(&[], 0.5).expect_err("an empty vector cannot bracket anything");
    assert!(
        e.to_string().contains("at least 2 nodes"),
        "expected the two-node refusal, got: {e}"
    );
}

#[test]
fn single_node_vector_does_not_return_a_bogus_bracket() {
    // Before the fix this returned Ok(0) — a bracket index whose `i + 1` node
    // does not exist, leaving the panic to whichever caller dereferenced it.
    let e = bisect_vector(&[1.0], 0.5).expect_err("one node is not a bracket");
    assert!(
        e.to_string().contains("at least 2 nodes"),
        "expected the two-node refusal, got: {e}"
    );
}

#[test]
fn a_well_formed_vector_still_bisects() {
    let v = [0.0, 1.0, 2.0, 3.0, 4.0];
    let i = bisect_vector(&v, 2.5).expect("a 5-node vector brackets 2.5");
    assert!(i < v.len() - 1, "index {i} must leave an i+1 node");
    assert!(
        v[i] <= 2.5 && 2.5 <= v[i + 1],
        "2.5 must lie in [{}, {}]",
        v[i],
        v[i + 1]
    );
}

#[test]
fn matrix_bisect_refuses_degenerate_shapes() {
    let empty: Vec<Vec<f64>> = vec![];
    assert!(
        bisect_segmented_vector_slice(&empty, 0, 0.5).is_err(),
        "an empty matrix has no row 0"
    );

    let one_short_row = vec![vec![1.0f64]];
    assert!(
        bisect_segmented_vector_slice(&one_short_row, 0, 0.5).is_err(),
        "a 1x1 matrix cannot bracket"
    );

    let rows = vec![vec![0.0, 1.0], vec![2.0, 3.0]];
    assert!(
        bisect_segmented_vector_slice(&rows, 5, 0.5).is_err(),
        "row 5 is out of range"
    );
}
