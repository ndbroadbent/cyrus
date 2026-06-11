//! Regression tests for intersection numbers on very small polytopes
//! (h11 <= 3, i.e. at most 8 lattice points in 4D).
//!
//! These exercise the landscape-smoke failure class: with n <= dim + 4
//! points, `build_linear_relations_from_rref` used to cap the relation rows
//! at h11 = n - (dim + 1), dropping a genuine relation and leaving the
//! intersection-number least-squares system rank-deficient (loud
//! SingularMatrix failures, or NaN solutions before the residual gate).
//!
//! Fixtures generated with CYTools (vendored cytools_latest) on
//! Kreuzer-Skarke polytopes: the expected values are
//! `cy.intersection_numbers(in_basis=True)` in the recorded divisor basis.

use cyrus_core::{
    Point, Polytope, compute_glsm_and_linrels, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_regular_triangulation, intersection_in_basis,
};

struct Fixture {
    name: &'static str,
    points: &'static [[i64; 4]],
    heights: &'static [f64],
    cytools_basis: &'static [usize],
    /// Entries (i, j, k, value) of the in-basis tensor, sorted indices;
    /// every combination not listed must be zero.
    kappa: &'static [(usize, usize, usize, i64)],
    n_simplices: usize,
    expected_linrel_rows: usize,
}

const H11_2: Fixture = Fixture {
    name: "KS h11=2 (7 points)",
    points: &[
        [0, 0, 0, 0],
        [-2, -2, -2, -1],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [-1, -1, -1, 0],
    ],
    heights: &[0.0, 13.0, 1.0, 1.0, 1.0, 1.0, 3.0],
    cytools_basis: &[2, 6],
    kappa: &[(1, 1, 1, -16), (0, 1, 1, 4)],
    n_simplices: 8,
    expected_linrel_rows: 4,
};

const H11_3: Fixture = Fixture {
    name: "KS h11=3 (8 points)",
    points: &[
        [0, 0, 0, 0],
        [-1, -1, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [-1, 0, -1, -1],
        [0, -1, -1, -1],
        [0, 1, 0, 0],
        [1, 0, 0, 0],
    ],
    heights: &[
        -0.8333333333333334,
        0.0,
        0.0,
        -0.3333333333333334,
        1.0,
        0.0,
        0.0,
        0.0,
    ],
    cytools_basis: &[2, 5, 7],
    kappa: &[
        (0, 1, 2, 3),
        (0, 0, 0, 2),
        (0, 0, 1, 1),
        (0, 0, 2, 3),
        (0, 1, 1, -1),
        (0, 2, 2, 1),
        (1, 1, 1, -4),
        (1, 1, 2, 3),
        (1, 2, 2, 1),
        (2, 2, 2, -1),
    ],
    n_simplices: 13,
    expected_linrel_rows: 4,
};

fn check_fixture(fixture: &Fixture) {
    let points: Vec<Point> = fixture
        .points
        .iter()
        .map(|p| Point::new(p.to_vec()))
        .collect();
    let polytope = Polytope::from_vertices(points).expect("polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("triangulation points");
    assert_eq!(
        triangulation_points.len(),
        fixture.points.len(),
        "{}: all fixture points should be triangulation points",
        fixture.name
    );

    let points_i64: Vec<Vec<i64>> = triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linrels = compute_linear_relations_no_origin(&points_i64);
    assert_eq!(
        linrels.len(),
        fixture.expected_linrel_rows,
        "{}: small polytopes need every non-origin RREF pivot relation \
         (the old h11 cap dropped one and made the system rank-deficient)",
        fixture.name
    );

    let triangulation = compute_regular_triangulation(&triangulation_points, fixture.heights)
        .expect("regular triangulation");
    assert_eq!(
        triangulation.simplices().len(),
        fixture.n_simplices,
        "{}: simplex count must match CYTools",
        fixture.name
    );

    let (_, _, _) = compute_glsm_and_linrels(&triangulation_points).expect("glsm");
    let linrels_i64: Vec<Vec<i64>> = linrels
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("linrel fits i64"))
                .collect()
        })
        .collect();
    let kappa_full =
        compute_intersection_cytools(&triangulation, &triangulation_points, &linrels_i64)
            .expect("intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, fixture.cytools_basis);

    let dim = fixture.cytools_basis.len();
    for i in 0..dim {
        for j in i..dim {
            for k in j..dim {
                let expected = fixture
                    .kappa
                    .iter()
                    .find(|&&(a, b, c, _)| (a, b, c) == (i, j, k))
                    .map_or(0, |&(_, _, _, value)| value);
                let actual = kappa_basis.get(i, j, k).to_string();
                assert_eq!(
                    actual,
                    expected.to_string(),
                    "{}: kappa({i},{j},{k}) must match CYTools",
                    fixture.name
                );
            }
        }
    }
}

#[test]
fn small_h11_2_intersection_numbers_match_cytools() {
    check_fixture(&H11_2);
}

#[test]
fn small_h11_3_intersection_numbers_match_cytools() {
    check_fixture(&H11_3);
}
