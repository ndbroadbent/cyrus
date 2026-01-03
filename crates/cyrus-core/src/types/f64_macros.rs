//! Compile-time verified macros for creating F64 values.
//!
//! These macros verify constraints at compile time, so invalid values
//! cause compilation errors rather than runtime panics.
//!
//! # Example
//!
//! ```
//! use cyrus_core::{f64_pos, f64_neg, f64_finite};
//!
//! let x = f64_pos!(3.14);       // F64<Pos>
//! let y = f64_neg!(-2.5);       // F64<Neg>
//! let z = f64_finite!(0.0);     // F64 (F64<Finite>)
//!
//! // These would fail at compile time:
//! // let bad = f64_pos!(-1.0);        // negative value
//! // let bad = f64_pos!(f64::NAN);    // NaN
//! // let bad = f64_neg!(1.0);         // positive value
//! ```

/// Check if a const f64 is finite (not NaN, not ±∞).
///
/// Used by compile-time macros since `is_finite()` is not const.
#[doc(hidden)]
#[must_use]
pub const fn const_is_finite(v: f64) -> bool {
    // NaN != NaN, so v == v is false for NaN
    // Infinity and -Infinity are equal to themselves but fail the bound check
    v == v && v != f64::INFINITY && v != f64::NEG_INFINITY
}

/// Create an `F64<Pos>` from a positive literal. Compile-time verified.
///
/// # Example
///
/// ```
/// use cyrus_core::f64_pos;
///
/// let x = f64_pos!(3.14);
/// assert!(x.get() > 0.0);
///
/// // These would fail at compile time:
/// // let bad = f64_pos!(-1.0);
/// // let bad = f64_pos!(f64::NAN);
/// // let bad = f64_pos!(f64::INFINITY);
/// ```
#[macro_export]
macro_rules! f64_pos {
    ($val:expr) => {{
        const V: f64 = $val;
        const _: () = assert!(
            $crate::types::f64_macros::const_is_finite(V) && V > 0.0,
            "f64_pos! requires finite positive value"
        );
        $crate::types::f64::F64::<$crate::types::tags::Pos>::from_raw(V)
    }};
}

/// Create an `F64<Neg>` from a negative literal. Compile-time verified.
///
/// # Example
///
/// ```
/// use cyrus_core::f64_neg;
///
/// let x = f64_neg!(-2.5);
/// assert!(x.get() < 0.0);
///
/// // These would fail at compile time:
/// // let bad = f64_neg!(1.0);
/// // let bad = f64_neg!(f64::NEG_INFINITY);
/// ```
#[macro_export]
macro_rules! f64_neg {
    ($val:expr) => {{
        const V: f64 = $val;
        const _: () = assert!(
            $crate::types::f64_macros::const_is_finite(V) && V < 0.0,
            "f64_neg! requires finite negative value"
        );
        $crate::types::f64::F64::<$crate::types::tags::Neg>::from_raw(V)
    }};
}

/// Create an `F64<Finite>` (or just `F64`) from a finite literal. Compile-time verified.
///
/// # Example
///
/// ```
/// use cyrus_core::f64_finite;
///
/// let x = f64_finite!(0.0);
/// let y = f64_finite!(-3.5);
///
/// // These would fail at compile time:
/// // let bad = f64_finite!(f64::NAN);
/// // let bad = f64_finite!(f64::INFINITY);
/// ```
#[macro_export]
macro_rules! f64_finite {
    ($val:expr) => {{
        const V: f64 = $val;
        const _: () = assert!(
            $crate::types::f64_macros::const_is_finite(V),
            "f64_finite! requires finite value (not NaN, not infinity)"
        );
        $crate::types::f64::F64::<$crate::types::tags::Finite>::from_raw(V)
    }};
}

/// Create an `F64<NonNeg>` from a non-negative literal. Compile-time verified.
///
/// # Example
///
/// ```
/// use cyrus_core::f64_nonneg;
///
/// let x = f64_nonneg!(0.0);
/// let y = f64_nonneg!(5.0);
/// ```
#[macro_export]
macro_rules! f64_nonneg {
    ($val:expr) => {{
        const V: f64 = $val;
        const _: () = assert!(
            $crate::types::f64_macros::const_is_finite(V) && V >= 0.0,
            "f64_nonneg! requires finite non-negative value"
        );
        $crate::types::f64::F64::<$crate::types::tags::NonNeg>::from_raw(V)
    }};
}

/// Create an `F64<NonPos>` from a non-positive literal. Compile-time verified.
///
/// # Example
///
/// ```
/// use cyrus_core::f64_nonpos;
///
/// let x = f64_nonpos!(0.0);
/// let y = f64_nonpos!(-5.0);
/// ```
#[macro_export]
macro_rules! f64_nonpos {
    ($val:expr) => {{
        const V: f64 = $val;
        const _: () = assert!(
            $crate::types::f64_macros::const_is_finite(V) && V <= 0.0,
            "f64_nonpos! requires finite non-positive value"
        );
        $crate::types::f64::F64::<$crate::types::tags::NonPos>::from_raw(V)
    }};
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
#[allow(clippy::float_cmp, clippy::approx_constant)]
mod tests {
    #[test]
    fn test_f64_pos_macro() {
        let x = f64_pos!(3.17); // Use non-pi value
        assert!(x.get() > 0.0);
    }

    #[test]
    fn test_f64_neg_macro() {
        let x = f64_neg!(-2.5);
        assert!(x.get() < 0.0);
    }

    #[test]
    fn test_f64_finite_macro() {
        let x = f64_finite!(0.0);
        assert_eq!(x.get(), 0.0);

        let y = f64_finite!(-5.5);
        assert_eq!(y.get(), -5.5);
    }

    #[test]
    fn test_f64_nonneg_macro() {
        let x = f64_nonneg!(0.0);
        assert_eq!(x.get(), 0.0);

        let y = f64_nonneg!(5.0);
        assert!(y.get() > 0.0);
    }

    #[test]
    fn test_f64_nonpos_macro() {
        let x = f64_nonpos!(0.0);
        assert_eq!(x.get(), 0.0);

        let y = f64_nonpos!(-5.0);
        assert!(y.get() < 0.0);
    }
}
