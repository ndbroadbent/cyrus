//! Type-safe f64 wrapper with compile-time invariant tracking.
//!
//! `F64<Tag>` wraps an `f64` with a phantom type that encodes constraints:
//!
//! - `F64` or `F64<Finite>`: Any finite value (no NaN, no ±∞)
//! - `F64<Pos>`: Strictly positive (> 0)
//! - `F64<Neg>`: Strictly negative (< 0)
//! - `F64<Zero>`: Exactly zero
//! - `F64<One>`: Exactly 1.0
//! - `F64<MinusOne>`: Exactly -1.0
//! - `F64<NonNeg>`: Non-negative (≥ 0)
//! - `F64<NonPos>`: Non-positive (≤ 0)
//!
//! # Example
//!
//! ```
//! use cyrus_core::{pos, finite};
//! use cyrus_core::types::f64::F64;
//!
//! let x = pos!(4.0);        // Compile-time verified F64<Pos>
//! let y = x.sqrt();         // sqrt of Pos is Pos
//! assert!((y.get() - 2.0).abs() < 1e-10);
//!
//! let z = finite!(-3.5);    // Compile-time verified F64 (F64<Finite>)
//! ```

use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

pub use super::tags::{Finite, GTEOne, IsFinite, IsGTEOne, IsNegative, IsNonNeg, IsNonPos, IsNonZero, IsOne, IsPositive, IsMinusOne, IsTwo, IsZero, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero, tag};

// ============================================================================
// F64<Tag> Wrapper
// ============================================================================

/// A type-safe f64 wrapper with compile-time invariant tracking.
///
/// The default tag is `Finite`, so `F64` means `F64<Finite>`.
///
/// # Representation
///
/// `#[repr(transparent)]` ensures this has the same memory layout as `f64`.
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct F64<Tag = Finite>(pub(crate) f64, pub(crate) PhantomData<Tag>);

// ============================================================================
// Core Impls
// ============================================================================

impl<Tag> F64<Tag> {
    /// Get the underlying f64 value.
    #[inline(always)]
    #[must_use]
    pub fn get(self) -> f64 {
        self.0
    }

    /// Create from a raw f64 without validation.
    ///
    /// Used internally and by macros after compile-time verification.
    #[doc(hidden)]
    #[inline(always)]
    #[must_use]
    pub const fn from_raw(x: f64) -> Self {
        Self(x, PhantomData)
    }
}

impl<Tag> std::fmt::Debug for F64<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let full_name = std::any::type_name::<Tag>();
        let tag_name = full_name.rsplit("::").next().unwrap_or(full_name);
        write!(f, "F64<{tag_name}>({:?})", self.0)
    }
}

impl<Tag> std::fmt::Display for F64<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<Tag> PartialEq for F64<Tag> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

// Eq is safe because all tags guarantee finite values (no NaN).
impl<Tag> Eq for F64<Tag> {}

impl<Tag> Hash for F64<Tag> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Normalize -0.0 to 0.0 for consistent hashing
        let v = if self.0 == 0.0 { 0.0 } else { self.0 };
        v.to_bits().hash(state);
    }
}

impl<Tag> PartialOrd for F64<Tag> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

// ============================================================================
// Constructors
// ============================================================================

impl F64<Finite> {
    /// Create a finite f64. Returns `None` if NaN or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        x.is_finite().then(|| Self(x, PhantomData))
    }

    /// Zero constant.
    pub const ZERO: Self = Self(0.0, PhantomData);

    /// Compute absolute value. Finite.abs() is NonNeg (could be zero).
    #[must_use]
    pub fn abs(self) -> F64<NonNeg> {
        F64(self.0.abs(), PhantomData)
    }

    /// Try to narrow to positive. Returns `None` if ≤ 0.
    #[must_use]
    pub fn try_to_pos(self) -> Option<F64<Pos>> {
        (self.0 > 0.0).then(|| F64(self.0, PhantomData))
    }

    /// Try to narrow to negative. Returns `None` if ≥ 0.
    #[must_use]
    pub fn try_to_neg(self) -> Option<F64<Neg>> {
        (self.0 < 0.0).then(|| F64(self.0, PhantomData))
    }

    /// Try to narrow to non-zero. Returns `None` if = 0.
    #[must_use]
    pub fn try_to_non_zero(self) -> Option<F64<NonZero>> {
        (self.0 != 0.0).then(|| F64(self.0, PhantomData))
    }
}

impl F64<Pos> {
    /// Create a positive f64. Returns `None` if ≤ 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x > 0.0).then(|| Self(x, PhantomData))
    }

    /// Compute absolute value. Pos.abs() is Pos (already positive).
    #[must_use]
    pub fn abs(self) -> F64<Pos> {
        self
    }
}

impl F64<Neg> {
    /// Create a negative f64. Returns `None` if ≥ 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x < 0.0).then(|| Self(x, PhantomData))
    }

    /// Compute absolute value. Neg.abs() is Pos.
    #[must_use]
    pub fn abs(self) -> F64<Pos> {
        F64(self.0.abs(), PhantomData)
    }
}

impl F64<NonZero> {
    /// Create a non-zero f64. Returns `None` if = 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x != 0.0).then(|| Self(x, PhantomData))
    }

    /// Compute absolute value. NonZero.abs() is Pos (not just NonNeg).
    #[must_use]
    pub fn abs(self) -> F64<Pos> {
        F64(self.0.abs(), PhantomData)
    }
}

impl F64<NonNeg> {
    /// Create a non-negative f64. Returns `None` if < 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x >= 0.0).then(|| Self(x, PhantomData))
    }
}

impl F64<NonPos> {
    /// Create a non-positive f64. Returns `None` if > 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x <= 0.0).then(|| Self(x, PhantomData))
    }
}

impl F64<Zero> {
    /// The zero constant.
    pub const ZERO: Self = Self(0.0, PhantomData);

    /// Create a zero. Returns `None` if not exactly 0.0.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x == 0.0).then_some(Self::ZERO)
    }
}

impl F64<One> {
    /// The one constant.
    pub const ONE: Self = Self(1.0, PhantomData);

    /// Create a one. Returns `None` if not exactly 1.0.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x == 1.0).then_some(Self::ONE)
    }
}

impl F64<MinusOne> {
    /// The minus-one constant.
    pub const MINUS_ONE: Self = Self(-1.0, PhantomData);

    /// Create a minus-one. Returns `None` if not exactly -1.0.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x == -1.0).then_some(Self::MINUS_ONE)
    }
}

impl F64<Two> {
    /// The two constant.
    pub const TWO: Self = Self(2.0, PhantomData);

    /// Create a two. Returns `None` if not exactly 2.0.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x == 2.0).then_some(Self::TWO)
    }
}

impl F64<GTEOne> {
    /// Create a value ≥ 1. Returns `None` if < 1, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x >= 1.0).then(|| Self(x, PhantomData))
    }
}

impl Default for F64<Zero> {
    fn default() -> Self {
        Self::ZERO
    }
}

// ============================================================================
// Trait Implementations
// ============================================================================

impl IsFinite for F64<Finite> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { self }
}

impl IsFinite for F64<Pos> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<Neg> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<Zero> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<One> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<MinusOne> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<NonZero> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<NonNeg> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<NonPos> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<Two> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

impl IsFinite for F64<GTEOne> {
    type Finite = F64<Finite>;
    fn to_finite(self) -> Self::Finite { F64(self.0, PhantomData) }
}

// --- IsNonZero ---

impl IsNonZero for F64<NonZero> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { self }
}

impl IsNonZero for F64<Pos> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { F64(self.0, PhantomData) }
}

impl IsNonZero for F64<Neg> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { F64(self.0, PhantomData) }
}

impl IsNonZero for F64<One> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { F64(self.0, PhantomData) }
}

impl IsNonZero for F64<MinusOne> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { F64(self.0, PhantomData) }
}

impl IsNonZero for F64<Two> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { F64(self.0, PhantomData) }
}

impl IsNonZero for F64<GTEOne> {
    type NonZero = F64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero { F64(self.0, PhantomData) }
}

// --- IsPositive ---

impl IsPositive for F64<Pos> {
    type Positive = F64<Pos>;
    fn to_positive(self) -> Self::Positive { self }
}

impl IsPositive for F64<One> {
    type Positive = F64<Pos>;
    fn to_positive(self) -> Self::Positive { F64(self.0, PhantomData) }
}

impl IsPositive for F64<Two> {
    type Positive = F64<Pos>;
    fn to_positive(self) -> Self::Positive { F64(self.0, PhantomData) }
}

impl IsPositive for F64<GTEOne> {
    type Positive = F64<Pos>;
    fn to_positive(self) -> Self::Positive { F64(self.0, PhantomData) }
}

// --- IsNegative ---

impl IsNegative for F64<Neg> {
    type Negative = F64<Neg>;
    fn to_negative(self) -> Self::Negative { self }
}

impl IsNegative for F64<MinusOne> {
    type Negative = F64<Neg>;
    fn to_negative(self) -> Self::Negative { F64(self.0, PhantomData) }
}

// --- IsZero ---

impl IsZero for F64<Zero> {}

// --- IsNonNeg ---

impl IsNonNeg for F64<NonNeg> {
    type NonNeg = F64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg { self }
}

impl IsNonNeg for F64<Pos> {
    type NonNeg = F64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg { F64(self.0, PhantomData) }
}

impl IsNonNeg for F64<One> {
    type NonNeg = F64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg { F64(self.0, PhantomData) }
}

impl IsNonNeg for F64<Zero> {
    type NonNeg = F64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg { F64(self.0, PhantomData) }
}

impl IsNonNeg for F64<Two> {
    type NonNeg = F64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg { F64(self.0, PhantomData) }
}

impl IsNonNeg for F64<GTEOne> {
    type NonNeg = F64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg { F64(self.0, PhantomData) }
}

// --- IsNonPos ---

impl IsNonPos for F64<NonPos> {
    type NonPos = F64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos { self }
}

impl IsNonPos for F64<Neg> {
    type NonPos = F64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos { F64(self.0, PhantomData) }
}

impl IsNonPos for F64<MinusOne> {
    type NonPos = F64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos { F64(self.0, PhantomData) }
}

impl IsNonPos for F64<Zero> {
    type NonPos = F64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos { F64(self.0, PhantomData) }
}

// --- IsOne / IsMinusOne ---

impl IsOne for F64<One> {}
impl IsTwo for F64<Two> {}
impl IsMinusOne for F64<MinusOne> {}
impl IsGTEOne for F64<GTEOne> {}
impl IsGTEOne for F64<One> {}
impl IsGTEOne for F64<Two> {}

// ============================================================================
// Arithmetic Operations (generated from algebra rules)
// ============================================================================

crate::impl_ops!(F64, f64);

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f64_size() {
        assert_eq!(std::mem::size_of::<F64>(), std::mem::size_of::<f64>());
        assert_eq!(std::mem::size_of::<F64<Pos>>(), std::mem::size_of::<f64>());
    }

    #[test]
    fn test_finite_new() {
        assert!(F64::<Finite>::new(0.0).is_some());
        assert!(F64::<Finite>::new(1.0).is_some());
        assert!(F64::<Finite>::new(-1.0).is_some());
        assert!(F64::<Finite>::new(f64::NAN).is_none());
        assert!(F64::<Finite>::new(f64::INFINITY).is_none());
    }

    #[test]
    fn test_pos_new() {
        assert!(F64::<Pos>::new(1.0).is_some());
        assert!(F64::<Pos>::new(0.0).is_none());
        assert!(F64::<Pos>::new(-1.0).is_none());
    }

    #[test]
    fn test_neg_new() {
        assert!(F64::<Neg>::new(-1.0).is_some());
        assert!(F64::<Neg>::new(0.0).is_none());
        assert!(F64::<Neg>::new(1.0).is_none());
    }

    #[test]
    fn test_constants() {
        assert_eq!(F64::<Zero>::ZERO.get(), 0.0);
        assert_eq!(F64::<One>::ONE.get(), 1.0);
        assert_eq!(F64::<MinusOne>::MINUS_ONE.get(), -1.0);
    }

    #[test]
    fn test_display() {
        let x = F64::<Pos>::new(3.14).unwrap();
        assert_eq!(format!("{x}"), "3.14");
    }

    #[test]
    fn test_debug() {
        let x = F64::<Pos>::new(3.14).unwrap();
        let s = format!("{x:?}");
        assert!(s.contains("Pos"));
        assert!(s.contains("3.14"));
    }

    #[test]
    fn test_hash() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(F64::<Pos>::new(1.0).unwrap());
        set.insert(F64::<Pos>::new(2.0).unwrap());
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_widening() {
        let p = F64::<Pos>::new(1.0).unwrap();
        let o = F64::<One>::ONE;

        let _: F64<Finite> = p.to_finite();
        let _: F64<Pos> = o.to_positive();
        let _: F64<NonNeg> = p.to_non_neg();
    }
}
