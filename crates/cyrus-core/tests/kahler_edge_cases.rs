#![allow(missing_docs)]
use cyrus_core::{Point, Triangulation, compute_mori_generators};

#[test]
fn test_mori_empty_points() {
    let tri = Triangulation::new(vec![]);
    let points: Vec<Point> = Vec::new();
    let res = compute_mori_generators(&tri, &points);
    assert!(res.is_err());
}
