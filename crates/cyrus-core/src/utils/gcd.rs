//! GCD (Greatest Common Divisor) utilities.
//!
//! Direct port of CYTools `utils.py` gcd functions.

/// Compute the greatest common divisor of two floats using Euclidean algorithm.
///
/// This finds the largest floating-point number that divides both `a` and `b`.
///
/// # Warning
///
/// Unexpected behavior can occur if `b` starts below `tol`.
/// E.g., `gcd_float(100.0, 0.1, 0.2)` returns `100.0`.
///
/// # Arguments
///
/// * `a` - First number
/// * `b` - Second number
/// * `tol` - Tolerance for considering a value as zero (default: 1e-5)
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::gcd_float;
///
/// let g = gcd_float(0.2, 0.5, 1e-5);
/// assert!((g - 0.1).abs() < 1e-10);
/// ```
#[must_use]
pub fn gcd_float(a: f64, b: f64, tol: f64) -> f64 {
    if b.abs() < tol {
        a.abs()
    } else {
        gcd_float(b, a % b, tol)
    }
}

/// Compute the GCD of two floats with default tolerance (1e-5).
#[must_use]
pub fn gcd_float_default(a: f64, b: f64) -> f64 {
    gcd_float(a, b, 1e-5)
}

/// Compute the GCD of a slice of floats.
///
/// # Arguments
///
/// * `values` - Slice of floating-point values
/// * `tol` - Tolerance for the GCD computation
///
/// # Returns
///
/// The GCD of all values, or 0.0 if the slice is empty.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::gcd_list;
///
/// let values = [0.2, 0.4, 0.6];
/// let g = gcd_list(&values, 1e-5);
/// assert!((g - 0.2).abs() < 1e-10);
/// ```
#[must_use]
pub fn gcd_list(values: &[f64], tol: f64) -> f64 {
    values
        .iter()
        .copied()
        .reduce(|acc, x| gcd_float(acc, x, tol))
        .unwrap_or(0.0)
}

/// Compute the GCD of a slice of floats with default tolerance (1e-5).
#[must_use]
pub fn gcd_list_default(values: &[f64]) -> f64 {
    gcd_list(values, 1e-5)
}

/// Compute the GCD of a slice of integers.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::gcd_list_int;
///
/// assert_eq!(gcd_list_int(&[12, 18, 24]), 6);
/// assert_eq!(gcd_list_int(&[7, 11, 13]), 1);
/// assert_eq!(gcd_list_int(&[]), 0);
/// ```
#[must_use]
pub fn gcd_list_int(values: &[i64]) -> i64 {
    values.iter().copied().reduce(gcd_int).unwrap_or(0)
}

/// Compute the GCD of two integers using Euclidean algorithm.
#[must_use]
pub const fn gcd_int(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd_float_basic() {
        let g = gcd_float(0.2, 0.5, 1e-5);
        assert!((g - 0.1).abs() < 1e-9, "Expected ~0.1, got {g}");
    }

    #[test]
    fn test_gcd_float_same_value() {
        let g = gcd_float(0.3, 0.3, 1e-5);
        assert!((g - 0.3).abs() < 1e-9);
    }

    #[test]
    fn test_gcd_float_one_zero() {
        let g = gcd_float(0.5, 0.0, 1e-5);
        assert!((g - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_gcd_float_integers() {
        let g = gcd_float(12.0, 18.0, 1e-5);
        assert!((g - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_gcd_list_basic() {
        let values = [0.2, 0.4, 0.6];
        let g = gcd_list(&values, 1e-5);
        assert!((g - 0.2).abs() < 1e-9);
    }

    #[test]
    fn test_gcd_list_single() {
        let g = gcd_list(&[0.7], 1e-5);
        assert!((g - 0.7).abs() < 1e-9);
    }

    #[test]
    fn test_gcd_list_empty() {
        let g = gcd_list(&[], 1e-5);
        assert!((g - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_gcd_int_basic() {
        assert_eq!(gcd_int(12, 18), 6);
        assert_eq!(gcd_int(48, 18), 6);
        assert_eq!(gcd_int(7, 11), 1);
    }

    #[test]
    fn test_gcd_int_negative() {
        assert_eq!(gcd_int(-12, 18), 6);
        assert_eq!(gcd_int(12, -18), 6);
        assert_eq!(gcd_int(-12, -18), 6);
    }

    #[test]
    fn test_gcd_int_zero() {
        assert_eq!(gcd_int(12, 0), 12);
        assert_eq!(gcd_int(0, 18), 18);
        assert_eq!(gcd_int(0, 0), 0);
    }

    #[test]
    fn test_gcd_list_int() {
        assert_eq!(gcd_list_int(&[12, 18, 24]), 6);
        assert_eq!(gcd_list_int(&[100, 50, 25]), 25);
        assert_eq!(gcd_list_int(&[7]), 7);
        assert_eq!(gcd_list_int(&[]), 0);
    }
}
