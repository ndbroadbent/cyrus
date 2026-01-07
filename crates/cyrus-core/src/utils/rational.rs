//! Rational number conversion utilities.
//!
//! Direct port of CYTools `utils.py` float/fmpq conversion functions.
//! Uses malachite instead of flint.

use malachite::num::arithmetic::traits::Abs;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use malachite::{Integer, Rational};

/// Convert a float to a rational number.
///
/// Uses continued fraction expansion to find a "reasonable" rational approximation.
/// This is equivalent to Python's `fractions.Fraction(c).limit_denominator()`.
///
/// Direct port of CYTools `utils.py` `float_to_fmpq` function.
///
/// # Arguments
///
/// * `c` - The input float
/// * `max_denominator` - Maximum denominator (default: 10^6)
///
/// # Returns
///
/// A rational approximation of the input.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::float_to_rational;
/// use malachite::Rational;
///
/// let r = float_to_rational(0.333333333333, None);
/// assert_eq!(r, Rational::from_signeds(1, 3));
///
/// let r = float_to_rational(0.1, None);
/// assert_eq!(r, Rational::from_signeds(1, 10));
/// ```
pub fn float_to_rational(c: f64, max_denominator: Option<i64>) -> Rational {
    let max_denom = max_denominator.unwrap_or(1_000_000);

    if c == 0.0 {
        return Rational::from(0);
    }

    // Use continued fraction expansion (like Python's Fraction.limit_denominator)
    // This finds the rational p/q with q <= max_denom that best approximates c
    let neg = c < 0.0;
    let c_abs = c.abs();

    // Initial bounds: p0/q0 = 0/1, p1/q1 = 1/0 (represents infinity)
    let mut p0: i64 = 0;
    let mut q0: i64 = 1;
    let mut p1: i64 = 1;
    let mut q1: i64 = 0;

    let mut x = c_abs;

    loop {
        let a = x.floor() as i64;

        // Compute next convergent
        let p2 = a * p1 + p0;
        let q2 = a * q1 + q0;

        if q2 > max_denom {
            break;
        }

        p0 = p1;
        q0 = q1;
        p1 = p2;
        q1 = q2;

        let frac = x - a as f64;
        if frac < 1e-15 {
            break;
        }
        x = 1.0 / frac;

        // Prevent infinite loops from floating point errors
        if x > 1e15 {
            break;
        }
    }

    let (p, q) = if neg { (-p1, q1) } else { (p1, q1) };
    Rational::from_signeds(p, q)
}

/// Convert a rational number to a float.
///
/// Direct port of CYTools `utils.py` `fmpq_to_float` function.
///
/// # Arguments
///
/// * `r` - The input rational number
///
/// # Returns
///
/// The float approximation.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::rational_to_float;
/// use malachite::Rational;
///
/// let f = rational_to_float(&Rational::from_signeds(1, 2));
/// assert!((f - 0.5).abs() < 1e-15);
///
/// let f = rational_to_float(&Rational::from_signeds(1, 3));
/// assert!((f - 0.333333333333333).abs() < 1e-10);
/// ```
pub fn rational_to_float(r: &Rational) -> f64 {
    f64::rounding_from(r, RoundingMode::Nearest).0
}

/// Convert an array of floats to rationals.
///
/// Direct port of CYTools `utils.py` `array_to_flint` function (for float -> rational).
///
/// # Arguments
///
/// * `arr` - The input array of floats
/// * `max_denominator` - Maximum denominator for conversion
///
/// # Returns
///
/// Array of rational approximations.
pub fn array_to_rational(arr: &[f64], max_denominator: Option<i64>) -> Vec<Rational> {
    arr.iter()
        .map(|&x| float_to_rational(x, max_denominator))
        .collect()
}

/// Convert an array of rationals to floats.
///
/// Direct port of CYTools `utils.py` `array_from_flint` function (for rational -> float).
///
/// # Arguments
///
/// * `arr` - The input array of rationals
///
/// # Returns
///
/// Array of float approximations.
pub fn array_from_rational(arr: &[Rational]) -> Vec<f64> {
    arr.iter().map(rational_to_float).collect()
}

/// Compute the GCD of a list of integers.
///
/// Returns 1 if the list is empty or all zeros.
pub fn gcd_integers(values: &[Integer]) -> Integer {
    let mut result = Integer::from(0);
    for v in values {
        result = gcd_two(&result, v);
    }
    if result == 0 {
        Integer::from(1)
    } else {
        result.abs()
    }
}

/// Compute GCD of two integers using Euclidean algorithm.
fn gcd_two(a: &Integer, b: &Integer) -> Integer {
    let mut a = a.clone().abs();
    let mut b = b.clone().abs();
    while b != 0 {
        let t = b.clone();
        b = &a % &b;
        a = t;
    }
    a
}

/// Compute the integer nullspace (kernel) of a matrix, optionally reduced by GCD.
///
/// This is the direct port of CYTools `utils.py` `integral_nullspace` function.
///
/// # Arguments
///
/// * `matrix` - The input matrix (row-major)
/// * `reduce_by_gcd` - Whether to divide each basis vector by its GCD (default: true)
///
/// # Returns
///
/// The nullspace basis vectors as rows.
///
/// # Note
///
/// CYTools returns column vectors; we return row vectors. The mathematical content
/// is the same - each returned Vec is one basis element of the kernel.
pub fn integral_nullspace(matrix: &[Vec<Integer>], reduce_by_gcd: bool) -> Vec<Vec<Integer>> {
    use crate::integer_math::integer_kernel;

    let mut kernel = integer_kernel(matrix);

    if reduce_by_gcd {
        for vec in &mut kernel {
            let g = gcd_integers(vec);
            if g > 1 {
                for x in vec.iter_mut() {
                    *x = &*x / &g;
                }
            }
        }
    }

    kernel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_float_to_rational_simple() {
        let r = float_to_rational(0.5, None);
        assert_eq!(r, Rational::from_signeds(1, 2));
    }

    #[test]
    fn test_float_to_rational_third() {
        let r = float_to_rational(0.333333333333, None);
        assert_eq!(r, Rational::from_signeds(1, 3));
    }

    #[test]
    fn test_float_to_rational_tenth() {
        let r = float_to_rational(0.1, None);
        assert_eq!(r, Rational::from_signeds(1, 10));
    }

    #[test]
    fn test_float_to_rational_negative() {
        let r = float_to_rational(-0.25, None);
        assert_eq!(r, Rational::from_signeds(-1, 4));
    }

    #[test]
    fn test_float_to_rational_integer() {
        let r = float_to_rational(5.0, None);
        assert_eq!(r, Rational::from(5));
    }

    #[test]
    fn test_float_to_rational_zero() {
        let r = float_to_rational(0.0, None);
        assert_eq!(r, Rational::from(0));
    }

    #[test]
    fn test_rational_to_float_half() {
        let f = rational_to_float(&Rational::from_signeds(1, 2));
        assert!((f - 0.5).abs() < 1e-15);
    }

    #[test]
    fn test_rational_to_float_third() {
        let f = rational_to_float(&Rational::from_signeds(1, 3));
        assert!((f - 0.333333333333333).abs() < 1e-10);
    }

    #[test]
    fn test_array_roundtrip() {
        let original = vec![0.5, 0.25, 0.125, 0.1];
        let rationals = array_to_rational(&original, None);
        let back = array_from_rational(&rationals);

        for (a, b) in original.iter().zip(back.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_gcd_integers() {
        let v = vec![Integer::from(12), Integer::from(18), Integer::from(24)];
        assert_eq!(gcd_integers(&v), Integer::from(6));
    }

    #[test]
    fn test_gcd_integers_with_negative() {
        let v = vec![Integer::from(-12), Integer::from(18), Integer::from(-24)];
        assert_eq!(gcd_integers(&v), Integer::from(6));
    }

    #[test]
    fn test_gcd_integers_coprime() {
        let v = vec![Integer::from(7), Integer::from(11), Integer::from(13)];
        assert_eq!(gcd_integers(&v), Integer::from(1));
    }

    #[test]
    fn test_integral_nullspace_simple() {
        // Matrix [[1, 0, 0], [0, 1, 0]] has a 1D nullspace: [0, 0, 1]
        let m = vec![
            vec![Integer::from(1), Integer::from(0), Integer::from(0)],
            vec![Integer::from(0), Integer::from(1), Integer::from(0)],
        ];
        let kernel = integral_nullspace(&m, true);
        assert_eq!(kernel.len(), 1);
        assert_eq!(
            kernel[0],
            vec![Integer::from(0), Integer::from(0), Integer::from(1)]
        );
    }

    #[test]
    fn test_integral_nullspace_2d() {
        // Matrix [[1, 1, 1]] has a 2D nullspace
        let m = vec![vec![Integer::from(1), Integer::from(1), Integer::from(1)]];
        let kernel = integral_nullspace(&m, true);
        assert_eq!(kernel.len(), 2);

        // Check all vectors are in the kernel
        for v in &kernel {
            let sum: Integer = m[0].iter().zip(v.iter()).map(|(a, b)| a * b).sum();
            assert_eq!(sum, Integer::from(0));
        }
    }

    #[test]
    fn test_integral_nullspace_gcd_reduction() {
        // Create a matrix whose kernel has a vector that needs GCD reduction
        let m = vec![vec![Integer::from(2), Integer::from(4)]];
        let kernel = integral_nullspace(&m, true);
        assert_eq!(kernel.len(), 1);

        // Should be reduced to [-2, 1] or [2, -1] (coprime)
        let v = &kernel[0];
        let g = gcd_integers(v);
        assert_eq!(g, Integer::from(1));
    }
}
