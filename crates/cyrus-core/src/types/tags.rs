//! Type tags for compile-time numeric invariants.
//!
//! Tags are zero-sized types that encode what we know about a value:
//! - [`Finite`] - No NaN or infinity (always true for Rational/Integer)
//! - [`NonZero`] - Not zero (≠ 0) - critical for safe division
//! - [`Pos`] - Strictly positive (> 0)
//! - [`Neg`] - Strictly negative (< 0)
//! - [`Zero`] - Exactly zero
//! - [`One`] - Exactly 1
//! - [`MinusOne`] - Exactly -1
//! - [`NonNeg`] - Non-negative (≥ 0)
//! - [`NonPos`] - Non-positive (≤ 0)
//!
//! ## Trait Hierarchy
//!
//! ```text
//! IsFinite
//!    ├── IsNonZero
//!    │      ├── IsPositive (→ IsOne)
//!    │      └── IsNegative (→ IsMinusOne)
//!    ├── IsNonNeg
//!    ├── IsNonPos
//!    └── IsZero
//! ```
//!
//! Key insight: if `x: NonZero`, then `x.abs(): Pos` (not just `NonNeg`).

use std::marker::PhantomData;

// ============================================================================
// Tag Types (Zero-Sized)
// ============================================================================

/// Tag: finite value (no NaN, no ±∞). Default for all numeric wrappers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Finite;

/// Tag: strictly positive value (> 0)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pos;

/// Tag: strictly negative value (< 0)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Neg;

/// Tag: exactly zero
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Zero;

/// Tag: exactly 1
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct One;

/// Tag: exactly 2
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Two;

/// Tag: greater than or equal to 1 (≥ 1)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GTEOne;

/// Tag: exactly -1
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinusOne;

/// Tag: non-zero value (≠ 0) - critical for safe division
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonZero;

/// Tag: non-negative value (≥ 0)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonNeg;

/// Tag: non-positive value (≤ 0)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonPos;

// ============================================================================
// Trait Lattice
// ============================================================================

/// Marker trait for types that wrap a numeric value with a tag.
///
/// This is implemented by `F64<Tag>`, `Rational<Tag>`, `Integer<Tag>`.
/// Note: No `Copy` bound since `Rational` and `Integer` are heap-allocated.
pub trait Tagged<Tag> {
    /// The underlying numeric type (f64, Rational, Integer).
    type Inner;

    /// Get the underlying value.
    fn get(&self) -> &Self::Inner;

    /// Create from inner value without validation.
    ///
    /// # Safety (logical)
    /// Caller must ensure the value satisfies the tag's invariant.
    fn from_inner_unchecked(inner: Self::Inner) -> Self;
}

/// Guaranteed finite (no NaN, no ±∞). Base marker for all typed numerics.
pub trait IsFinite {}

/// Guaranteed non-zero (≠ 0). Safe to divide by.
pub trait IsNonZero: IsFinite {}

/// Guaranteed strictly positive (> 0).
pub trait IsPositive: IsNonZero {}

/// Guaranteed strictly negative (< 0).
pub trait IsNegative: IsNonZero {}

/// Guaranteed exactly zero.
pub trait IsZero: IsFinite {}

/// Guaranteed non-negative (≥ 0).
pub trait IsNonNeg: IsFinite {}

/// Guaranteed non-positive (≤ 0).
pub trait IsNonPos: IsFinite {}

/// Guaranteed exactly 1.
pub trait IsOne: IsPositive {}

/// Guaranteed exactly 2.
pub trait IsTwo: IsPositive {}

/// Guaranteed ≥ 1.
pub trait IsGTEOne: IsPositive {}

/// Guaranteed exactly -1.
pub trait IsMinusOne: IsNegative {}

// ============================================================================
// PhantomData helper
// ============================================================================

/// Helper to create PhantomData for tags.
#[inline]
pub const fn tag<T>() -> PhantomData<T> {
    PhantomData
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tags_are_zst() {
        assert_eq!(std::mem::size_of::<Finite>(), 0);
        assert_eq!(std::mem::size_of::<NonZero>(), 0);
        assert_eq!(std::mem::size_of::<Pos>(), 0);
        assert_eq!(std::mem::size_of::<Neg>(), 0);
        assert_eq!(std::mem::size_of::<Zero>(), 0);
        assert_eq!(std::mem::size_of::<One>(), 0);
        assert_eq!(std::mem::size_of::<Two>(), 0);
        assert_eq!(std::mem::size_of::<MinusOne>(), 0);
        assert_eq!(std::mem::size_of::<NonNeg>(), 0);
        assert_eq!(std::mem::size_of::<NonPos>(), 0);
        assert_eq!(std::mem::size_of::<GTEOne>(), 0);
    }

    #[test]
    fn test_finite_is_default() {
        let _: Finite = Finite;
    }

    #[test]
    fn test_tag_helper() {
        let _pd: PhantomData<Pos> = tag();
        let _pd2: PhantomData<NonZero> = tag();
    }

    #[test]
    fn test_tags_debug() {
        // Test Debug derives - use inline format to satisfy clippy
        let finite = Finite;
        let pos = Pos;
        let neg = Neg;
        let zero = Zero;
        let one = One;
        let two = Two;
        let gte_one = GTEOne;
        let minus_one = MinusOne;
        let non_zero = NonZero;
        let non_neg = NonNeg;
        let non_pos = NonPos;

        assert!(format!("{finite:?}").contains("Finite"));
        assert!(format!("{pos:?}").contains("Pos"));
        assert!(format!("{neg:?}").contains("Neg"));
        assert!(format!("{zero:?}").contains("Zero"));
        assert!(format!("{one:?}").contains("One"));
        assert!(format!("{two:?}").contains("Two"));
        assert!(format!("{gte_one:?}").contains("GTEOne"));
        assert!(format!("{minus_one:?}").contains("MinusOne"));
        assert!(format!("{non_zero:?}").contains("NonZero"));
        assert!(format!("{non_neg:?}").contains("NonNeg"));
        assert!(format!("{non_pos:?}").contains("NonPos"));
    }

    #[test]
    fn test_tags_copy_eq() {
        // Test Copy and PartialEq, Eq
        let p1 = Pos;
        let p2 = p1; // Copy
        assert_eq!(p1, p2);

        let n1 = Neg;
        let n2 = n1;
        assert_eq!(n1, n2);

        let z1 = Zero;
        let z2 = z1;
        assert_eq!(z1, z2);

        let o1 = One;
        let o2 = o1;
        assert_eq!(o1, o2);

        let t1 = Two;
        let t2 = t1;
        assert_eq!(t1, t2);

        let g1 = GTEOne;
        let g2 = g1;
        assert_eq!(g1, g2);

        let m1 = MinusOne;
        let m2 = m1;
        assert_eq!(m1, m2);

        let nz1 = NonZero;
        let nz2 = nz1;
        assert_eq!(nz1, nz2);

        let nn1 = NonNeg;
        let nn2 = nn1;
        assert_eq!(nn1, nn2);

        let np1 = NonPos;
        let np2 = np1;
        assert_eq!(np1, np2);

        let f1 = Finite;
        let f2 = f1;
        assert_eq!(f1, f2);
    }
}
