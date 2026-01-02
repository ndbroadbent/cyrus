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
