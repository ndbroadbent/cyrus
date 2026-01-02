#![allow(missing_docs)]
use cyrus_core::compute_v0;
use cyrus_core::integer_math::integer_kernel;
use malachite::Integer;

fn i(n: i64) -> Integer {
    Integer::from(n)
}

#[test]
fn test_kernel_sign_flip() {
    // A = [1, 1]
    // Kernel is [1, -1].
    // If we force the first element to be negative, does it flip?
    // integer_kernel implementation is deterministic based on HNF.
    // Let's check what it returns.

    let a = vec![vec![i(1), i(1)]];
    let kernel = integer_kernel(&a);
    // Standard HNF might return [1, -1] or [-1, 1].
    // My extract_mori_relation flips if [0] < 0.
    // But integer_kernel itself just returns basis.

    assert_eq!(kernel.len(), 1);
    // Just verifying it works.
}

#[test]
fn test_vacuum_error() {
    // compute_v0 inputs: ek0, g_s, v_string, w0
    // V0 = -3 * ek0 * (g_s^7 / (4*v_string)^2) * w0^2
    // If v_string is 0, division by zero -> infinity.

    let v0 = compute_v0(1.0, 0.1, 0.0, 1.0);
    assert!(v0.is_infinite());
}

#[test]
fn test_polytope_duplicate_vertex() {
    use cyrus_core::{Point, Polytope};
    // Create polytope with duplicate vertices
    let vertices = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![0, 0]), // Duplicate of first
    ];
    let result = Polytope::from_vertices(vertices);
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(err.to_string().contains("duplicate vertex"));
}

#[test]
fn test_polytope_vertices_accessor() {
    use cyrus_core::{Point, Polytope};
    let vertices = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![0, 1]),
    ];
    let poly = Polytope::from_vertices(vertices).unwrap();
    assert_eq!(poly.vertices().len(), 3);
    assert_eq!(poly.dim(), 2);
    assert_eq!(poly.num_vertices(), 3);
}

#[test]
fn test_kernel_complex_reduction() {
    // Matrix that requires multiple row swaps and reductions
    // This exercises the reduce_column loop more thoroughly
    let a = vec![
        vec![i(2), i(4), i(6)], // All multiples of 2
        vec![i(3), i(6), i(9)], // All multiples of 3
    ];
    let kernel = integer_kernel(&a);
    // Should have a 1-dimensional kernel
    assert!(!kernel.is_empty());
    // Verify kernel vectors satisfy A*v = 0
    for v in &kernel {
        for row in &a {
            let mut sum = i(0);
            for (j, coef) in row.iter().enumerate() {
                sum += coef * &v[j];
            }
            assert_eq!(sum, i(0));
        }
    }
}

#[test]
fn test_kernel_larger_matrix() {
    // Larger matrix to exercise more of the reduction paths
    let a = vec![
        vec![i(1), i(2), i(3), i(4)],
        vec![i(5), i(6), i(7), i(8)],
        vec![i(9), i(10), i(11), i(12)],
    ];
    let kernel = integer_kernel(&a);
    // Rank 2 matrix (row 3 is linear combination), so 2D kernel
    assert!(!kernel.is_empty());
}

#[test]
fn test_kernel_with_zeros() {
    // Matrix with zero elements requiring careful pivot selection
    let a = vec![vec![i(0), i(1), i(2)], vec![i(1), i(0), i(1)]];
    let kernel = integer_kernel(&a);
    // Should still compute kernel correctly
    assert!(!kernel.is_empty());
}
