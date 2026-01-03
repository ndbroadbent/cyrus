//! Type tags for compile-time numeric invariants.
//!
//! Tags are zero-sized types that encode what we know about a value:
//! - [`Finite`]: No NaN or infinity (always true for Rational/Integer)
//! - [`NonZero`]: Not zero (≠ 0) - critical for safe division
//! - [`Pos`]: Strictly positive (> 0)
//! - [`Neg`]: Strictly negative (< 0)
//! - [`Zero`]: Exactly zero
//! - [`One`]: Exactly 1
//! - [`MinusOne`]: Exactly -1
//! - [`NonNeg`]: Non-negative (≥ 0)
//! - [`NonPos`]: Non-positive (≤ 0)
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

/// Guaranteed finite (no NaN, no ±∞).
///
/// This is the base constraint - all tagged numeric types implement it.
/// For `Rational` and `Integer`, this is always true.
/// Note: No `Copy` bound since `Rational` and `Integer` are heap-allocated.
pub trait IsFinite {
    /// The numeric wrapper type (e.g., `F64<Finite>`).
    type Finite: IsFinite;

    /// Widen to the `Finite` variant, forgetting more specific constraints.
    fn to_finite(self) -> Self::Finite;
}

/// Guaranteed non-zero (≠ 0).
///
/// Critical for safe division. If `x: IsNonZero`, then:
/// - `1/x` is valid (no division by zero)
/// - `x.abs()` returns `Pos` (not just `NonNeg`)
pub trait IsNonZero: IsFinite {
    /// The non-zero variant type.
    type NonZero: IsNonZero;

    /// Widen to `NonZero`, forgetting `Pos`/`Neg` constraint if present.
    fn to_non_zero(self) -> Self::NonZero;
}

/// Guaranteed strictly positive (> 0).
///
/// Implies `IsNonZero` - positive values are never zero.
pub trait IsPositive: IsNonZero {
    /// The positive variant type.
    type Positive: IsPositive;

    /// Widen to `Pos`, forgetting `One` constraint if present.
    fn to_positive(self) -> Self::Positive;
}

/// Guaranteed strictly negative (< 0).
///
/// Implies `IsNonZero` - negative values are never zero.
pub trait IsNegative: IsNonZero {
    /// The negative variant type.
    type Negative: IsNegative;

    /// Widen to `Neg`, forgetting `MinusOne` constraint if present.
    fn to_negative(self) -> Self::Negative;
}

/// Guaranteed exactly zero.
pub trait IsZero: IsFinite {}

/// Guaranteed non-negative (≥ 0).
pub trait IsNonNeg: IsFinite {
    /// The non-negative variant type.
    type NonNeg: IsNonNeg;

    /// Widen to `NonNeg`.
    fn to_non_neg(self) -> Self::NonNeg;
}

/// Guaranteed non-positive (≤ 0).
pub trait IsNonPos: IsFinite {
    /// The non-positive variant type.
    type NonPos: IsNonPos;

    /// Widen to `NonPos`.
    fn to_non_pos(self) -> Self::NonPos;
}

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
#[inline(always)]
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
        assert_eq!(std::mem::size_of::<MinusOne>(), 0);
        assert_eq!(std::mem::size_of::<NonNeg>(), 0);
        assert_eq!(std::mem::size_of::<NonPos>(), 0);
    }

    #[test]
    fn test_finite_is_default() {
        let _: Finite = Default::default();
    }
}
