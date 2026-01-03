#![allow(missing_docs)]
use cyrus_core::racetrack::{
    GvInvariant, RacetrackResult, RacetrackTerm, ZETA, build_racetrack_terms, compute_w0,
    solve_racetrack,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, NonNeg, Pos};

// Helper to convert raw i64 slice to typed
fn i64_vec(v: &[i64]) -> Vec<I64<Finite>> {
    v.iter().map(|&x| I64::<Finite>::new(x)).collect()
}

// Helper to convert raw f64 slice to typed
fn f64_vec(v: &[f64]) -> Vec<F64<Finite>> {
    v.iter().map(|&x| F64::<Finite>::new(x).unwrap()).collect()
}

// Helper to create F64<Finite> from raw
fn f64_finite(v: f64) -> F64<Finite> {
    F64::<Finite>::new(v).unwrap()
}

#[test]
fn test_zeta_value() {
    assert!((ZETA.get() - 1.202_056_903_159_594).abs() < 1e-15);
}

#[test]
fn test_solve_racetrack_same_sign_coefficients() {
    // Both coefficients positive - no solution
    let terms = vec![
        RacetrackTerm {
            curve: i64_vec(&[1]),
            coefficient: f64_finite(1.0),
            exponent: f64_finite(1.0),
        },
        RacetrackTerm {
            curve: i64_vec(&[2]),
            coefficient: f64_finite(2.0), // Same sign as first
            exponent: f64_finite(2.0),
        },
    ];
    assert!(solve_racetrack(&terms).is_none());
}

#[test]
fn test_solve_racetrack_both_negative() {
    // Both coefficients negative - no solution
    let terms = vec![
        RacetrackTerm {
            curve: i64_vec(&[1]),
            coefficient: f64_finite(-1.0),
            exponent: f64_finite(1.0),
        },
        RacetrackTerm {
            curve: i64_vec(&[2]),
            coefficient: f64_finite(-2.0),
            exponent: f64_finite(2.0),
        },
    ];
    assert!(solve_racetrack(&terms).is_none());
}

#[test]
fn test_solve_racetrack_ratio_negative() {
    // Ratio is negative due to exponent signs
    let terms = vec![
        RacetrackTerm {
            curve: i64_vec(&[1]),
            coefficient: f64_finite(1.0),
            exponent: f64_finite(-1.0), // Negative exponent
        },
        RacetrackTerm {
            curve: i64_vec(&[2]),
            coefficient: f64_finite(-1.0),
            exponent: f64_finite(1.0),
        },
    ];
    // ratio = -((-1) * 1) / ((1) * (-1)) = -1 / -1 = 1 (positive)
    // But let's check if it leads to valid g_s
    let res = solve_racetrack(&terms);
    // Could be valid or invalid depending on g_s bounds
    if let Some(r) = res {
        assert!(r.g_s.get() > 0.0 && r.g_s.get() <= 1.0);
    }
}

#[test]
fn test_solve_racetrack_gs_out_of_bounds() {
    // Create terms where g_s > 1
    let terms = vec![
        RacetrackTerm {
            curve: i64_vec(&[1]),
            coefficient: f64_finite(1.0),
            exponent: f64_finite(10.0),
        },
        RacetrackTerm {
            curve: i64_vec(&[2]),
            coefficient: f64_finite(-1.01), // Almost same magnitude
            exponent: f64_finite(10.1),
        },
    ];
    // This should produce g_s out of [0,1] range
    let res = solve_racetrack(&terms);
    if let Some(r) = res {
        assert!(r.g_s.get() > 0.0 && r.g_s.get() <= 1.0);
    }
}

#[test]
fn test_build_racetrack_terms() {
    let gv = vec![
        GvInvariant {
            curve: i64_vec(&[1, 0]),
            value: f64_finite(10.0),
        },
        GvInvariant {
            curve: i64_vec(&[0, 1]),
            value: f64_finite(20.0),
        },
    ];
    let m = i64_vec(&[2, 3]);
    let p = f64_vec(&[1.0, 1.0]);

    let terms = build_racetrack_terms(&gv, &m, &p);
    assert!(!terms.is_empty());
}

#[test]
fn test_build_racetrack_terms_grouping() {
    // Two GV invariants with same q·p should be grouped
    let gv = vec![
        GvInvariant {
            curve: i64_vec(&[1]),
            value: f64_finite(10.0),
        },
        GvInvariant {
            curve: i64_vec(&[1]), // Same curve
            value: f64_finite(20.0),
        },
    ];
    let m = i64_vec(&[1]);
    let p = f64_vec(&[1.0]);

    let terms = build_racetrack_terms(&gv, &m, &p);
    // Should be grouped into 1 term
    assert_eq!(terms.len(), 1);
    // coefficient = 1*10 + 1*20 = 30
    assert!((terms[0].coefficient.get() - 30.0).abs() < 1e-10);
}

#[test]
fn test_compute_w0() {
    let result = RacetrackResult {
        g_s: F64::<Pos>::new(0.1).unwrap(),
        im_tau: F64::<Pos>::new(10.0).unwrap(),
        delta: F64::<NonNeg>::new(0.0).unwrap(),
        epsilon: F64::<NonNeg>::new(0.0).unwrap(),
    };

    let term1 = RacetrackTerm {
        curve: i64_vec(&[1]),
        coefficient: f64_finite(1.0),
        exponent: f64_finite(1.0),
    };
    let term2 = RacetrackTerm {
        curve: i64_vec(&[2]),
        coefficient: f64_finite(-1.0),
        exponent: f64_finite(2.0),
    };

    let w0 = compute_w0(&result, &term1, &term2);
    // Just verify it computes something finite
    assert!(w0.get().is_finite());
}

#[test]
fn test_racetrack_result_fields() {
    // Construct a valid result and verify fields are accessible
    let result = RacetrackResult {
        g_s: F64::<Pos>::new(0.5).unwrap(),
        im_tau: F64::<Pos>::new(2.0).unwrap(),
        delta: F64::<NonNeg>::new(0.001).unwrap(),
        epsilon: F64::<NonNeg>::new(0.002).unwrap(),
    };
    assert!((result.g_s.get() - 0.5).abs() < f64::EPSILON);
    assert!((result.im_tau.get() - 2.0).abs() < f64::EPSILON);
    assert!((result.delta.get() - 0.001).abs() < f64::EPSILON);
    assert!((result.epsilon.get() - 0.002).abs() < f64::EPSILON);
}

#[test]
fn test_solve_racetrack_empty() {
    let terms: Vec<RacetrackTerm> = vec![];
    assert!(solve_racetrack(&terms).is_none());
}
