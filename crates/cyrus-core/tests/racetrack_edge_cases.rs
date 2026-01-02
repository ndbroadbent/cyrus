#![allow(missing_docs)]
use cyrus_core::racetrack::{
    GvInvariant, RacetrackResult, RacetrackTerm, ZETA, build_racetrack_terms, compute_w0,
    solve_racetrack,
};

#[test]
fn test_zeta_value() {
    assert!((ZETA - 1.202_056_903_159_594).abs() < 1e-15);
}

#[test]
fn test_solve_racetrack_same_sign_coefficients() {
    // Both coefficients positive - no solution
    let terms = vec![
        RacetrackTerm {
            curve: vec![1],
            coefficient: 1.0,
            exponent: 1.0,
        },
        RacetrackTerm {
            curve: vec![2],
            coefficient: 2.0, // Same sign as first
            exponent: 2.0,
        },
    ];
    assert!(solve_racetrack(&terms).is_none());
}

#[test]
fn test_solve_racetrack_both_negative() {
    // Both coefficients negative - no solution
    let terms = vec![
        RacetrackTerm {
            curve: vec![1],
            coefficient: -1.0,
            exponent: 1.0,
        },
        RacetrackTerm {
            curve: vec![2],
            coefficient: -2.0,
            exponent: 2.0,
        },
    ];
    assert!(solve_racetrack(&terms).is_none());
}

#[test]
fn test_solve_racetrack_ratio_negative() {
    // Ratio is negative due to exponent signs
    let terms = vec![
        RacetrackTerm {
            curve: vec![1],
            coefficient: 1.0,
            exponent: -1.0, // Negative exponent
        },
        RacetrackTerm {
            curve: vec![2],
            coefficient: -1.0,
            exponent: 1.0,
        },
    ];
    // ratio = -((-1) * 1) / ((1) * (-1)) = -1 / -1 = 1 (positive)
    // But let's check if it leads to valid g_s
    let res = solve_racetrack(&terms);
    // Could be valid or invalid depending on g_s bounds
    if let Some(r) = res {
        assert!(r.g_s > 0.0 && r.g_s <= 1.0);
    }
}

#[test]
fn test_solve_racetrack_gs_out_of_bounds() {
    // Create terms where g_s > 1
    let terms = vec![
        RacetrackTerm {
            curve: vec![1],
            coefficient: 1.0,
            exponent: 10.0,
        },
        RacetrackTerm {
            curve: vec![2],
            coefficient: -1.01, // Almost same magnitude
            exponent: 10.1,
        },
    ];
    // This should produce g_s out of [0,1] range
    let res = solve_racetrack(&terms);
    if let Some(r) = res {
        assert!(r.g_s > 0.0 && r.g_s <= 1.0);
    }
}

#[test]
fn test_build_racetrack_terms() {
    let gv = vec![
        GvInvariant {
            curve: vec![1, 0],
            value: 10.0,
        },
        GvInvariant {
            curve: vec![0, 1],
            value: 20.0,
        },
    ];
    let m = vec![2, 3];
    let p = vec![1.0, 1.0];

    let terms = build_racetrack_terms(&gv, &m, &p);
    assert!(!terms.is_empty());
}

#[test]
fn test_build_racetrack_terms_grouping() {
    // Two GV invariants with same q·p should be grouped
    let gv = vec![
        GvInvariant {
            curve: vec![1],
            value: 10.0,
        },
        GvInvariant {
            curve: vec![1], // Same curve
            value: 20.0,
        },
    ];
    let m = vec![1];
    let p = vec![1.0];

    let terms = build_racetrack_terms(&gv, &m, &p);
    // Should be grouped into 1 term
    assert_eq!(terms.len(), 1);
    // coefficient = 1*10 + 1*20 = 30
    assert!((terms[0].coefficient - 30.0).abs() < 1e-10);
}

#[test]
fn test_compute_w0() {
    let result = RacetrackResult {
        g_s: 0.1,
        im_tau: 10.0,
        delta: 0.0,
        epsilon: 0.0,
    };

    let term1 = RacetrackTerm {
        curve: vec![1],
        coefficient: 1.0,
        exponent: 1.0,
    };
    let term2 = RacetrackTerm {
        curve: vec![2],
        coefficient: -1.0,
        exponent: 2.0,
    };

    let w0 = compute_w0(&result, &term1, &term2);
    // Just verify it computes something finite
    assert!(w0.is_finite());
}

#[test]
fn test_racetrack_result_default() {
    let res = RacetrackResult::default();
    assert!(res.g_s.abs() < f64::EPSILON);
    assert!(res.im_tau.abs() < f64::EPSILON);
    assert!(res.delta.abs() < f64::EPSILON);
    assert!(res.epsilon.abs() < f64::EPSILON);
}

#[test]
fn test_solve_racetrack_empty() {
    let terms: Vec<RacetrackTerm> = vec![];
    assert!(solve_racetrack(&terms).is_none());
}
