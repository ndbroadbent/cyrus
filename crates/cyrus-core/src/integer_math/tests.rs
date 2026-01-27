//! Tests for integer math module.

use super::*;
use malachite::Integer;
use malachite::Rational;

fn i(n: i64) -> Integer {
    Integer::from(n)
}

fn r(n: i64) -> Rational {
    Rational::from(n)
}

#[test]
fn test_solve_linear_system_rational_simple() {
    // 2x + y = 5
    // x + 3y = 5
    // Sol: x=2, y=1
    let m = vec![vec![r(2), r(1)], vec![r(1), r(3)]];
    let c = vec![r(5), r(5)];

    let x = solve_linear_system_rational(&m, &c).unwrap();
    assert_eq!(x[0], r(2));
    assert_eq!(x[1], r(1));
}

#[test]
fn test_integer_kernel_simple() {
    // A = [1, 1, 1] (1x3 matrix)
    let a = vec![vec![i(1), i(1), i(1)]];
    let kernel = integer_kernel(&a);

    if std::env::var_os("CYRUS_DEBUG_KERNEL").is_some() {
        eprintln!("kernel: {kernel:?}");
    }

    assert_eq!(kernel.len(), 2);
    for v in &kernel {
        let mut sum = i(0);
        for val in v {
            sum += val;
        }
        assert_eq!(sum, i(0));
    }
}

#[test]
fn test_integer_kernel_independent() {
    // A = [1, 0, 0]
    //     [0, 1, 0]
    let a = vec![vec![i(1), i(0), i(0)], vec![i(0), i(1), i(0)]];
    let kernel = integer_kernel(&a);
    assert_eq!(kernel.len(), 1);
    assert_eq!(kernel[0], vec![i(0), i(0), i(1)]);
}

#[test]
fn test_compute_linear_relations_from_points() {
    // McAllister dual polytope points (0-8)
    let points: Vec<Vec<i64>> = vec![
        vec![0, 0, 0, 0],    // 0: origin
        vec![-1, 2, -1, -1], // 1
        vec![1, -1, 0, 0],   // 2
        vec![-1, -1, 1, 1],  // 3
        vec![-1, -1, 1, 2],  // 4
        vec![-1, -1, 2, 1],  // 5
        vec![-1, -1, 2, 3],  // 6
        vec![-1, -1, 3, 2],  // 7
        vec![-1, -1, 2, 2],  // 8
    ];

    let lr = compute_linear_relations_no_origin(&points);

    // Expected CYTools linear_relations (from fixture, origin excluded)
    // Shape: 4 x 8 (h11=4, points 1-8)
    // CYTools uses +1 on pivot columns (identity matrix structure)
    let expected = [
        vec![i(1), i(0), i(0), i(0), i(-2), i(-2), i(-4), i(-2)],
        vec![i(0), i(1), i(0), i(0), i(-3), i(-3), i(-6), i(-3)],
        vec![i(0), i(0), i(1), i(0), i(1), i(-1), i(0), i(0)],
        vec![i(0), i(0), i(0), i(1), i(-1), i(1), i(-1), i(0)],
    ];

    println!("Computed linear_relations:");
    for (idx, row) in lr.iter().enumerate() {
        println!("  Row {idx}: {row:?}");
    }
    println!("Expected:");
    for (idx, row) in expected.iter().enumerate() {
        println!("  Row {idx}: {row:?}");
    }

    assert_eq!(lr.len(), expected.len(), "Wrong number of rows");
    for (row_idx, (actual, exp)) in lr.iter().zip(expected.iter()).enumerate() {
        assert_eq!(
            actual, exp,
            "Row {row_idx} mismatch:\nActual: {actual:?}\nExpected: {exp:?}"
        );
    }
}
