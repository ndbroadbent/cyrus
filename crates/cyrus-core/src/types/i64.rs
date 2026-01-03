//! Type-safe i64 wrapper with compile-time invariant tracking.
//!
//! `I64<Tag>` wraps an `i64` with a phantom type that encodes constraints:
//!
//! - `I64` or `I64<Finite>`: Any value (always finite for integers)
//! - `I64<Pos>`: Strictly positive (> 0)
//! - `I64<Neg>`: Strictly negative (< 0)
//! - `I64<Zero>`: Exactly zero
//! - `I64<One>`: Exactly 1
//! - `I64<MinusOne>`: Exactly -1
//! - `I64<NonNeg>`: Non-negative (≥ 0)
//! - `I64<NonPos>`: Non-positive (≤ 0)
//!
//! # Example
//!
//! ```
//! use cyrus_core::types::i64::I64;
//! use cyrus_core::types::tags::Pos;
//!
//! let flux = I64::<Pos>::new(100).unwrap();
//! assert!(flux.get() > 0);
//! ```

use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

pub use super::tags::{
    Finite, IsFinite, IsMinusOne, IsNegative, IsNonNeg, IsNonPos, IsNonZero, IsOne, IsPositive,
    IsZero, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Zero,
};

// ============================================================================
// I64<Tag> Wrapper
// ============================================================================

/// A type-safe i64 wrapper with compile-time invariant tracking.
///
/// The default tag is `Finite`, so `I64` means `I64<Finite>`.
/// For integers, `Finite` is always satisfied (no NaN/infinity).
///
/// # Representation
///
/// `#[repr(transparent)]` ensures this has the same memory layout as `i64`.
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct I64<Tag = Finite>(pub(crate) i64, pub(crate) PhantomData<Tag>);

// ============================================================================
// Core Impls
// ============================================================================

impl<Tag> I64<Tag> {
    /// Get the underlying i64 value.
    #[inline(always)]
    #[must_use]
    pub fn get(self) -> i64 {
        self.0
    }

    /// Create from a raw i64 without validation.
    ///
    /// Used internally and by macros after compile-time verification.
    #[doc(hidden)]
    #[inline(always)]
    #[must_use]
    pub const fn from_raw(x: i64) -> Self {
        Self(x, PhantomData)
    }
}

impl<Tag> std::fmt::Debug for I64<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let full_name = std::any::type_name::<Tag>();
        let tag_name = full_name.rsplit("::").next().unwrap_or(full_name);
        write!(f, "I64<{tag_name}>({})", self.0)
    }
}

impl<Tag> std::fmt::Display for I64<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<Tag> PartialEq for I64<Tag> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<Tag> Eq for I64<Tag> {}

impl<Tag> Hash for I64<Tag> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl<Tag> PartialOrd for I64<Tag> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Tag> Ord for I64<Tag> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

// ============================================================================
// Constructors
// ============================================================================

impl I64<Finite> {
    /// Create a finite i64. Always succeeds for integers.
    #[must_use]
    pub fn new(x: i64) -> Self {
        Self(x, PhantomData)
    }
}

impl I64<Pos> {
    /// Create a positive i64. Returns `None` if ≤ 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x > 0).then(|| Self(x, PhantomData))
    }
}

impl I64<Neg> {
    /// Create a negative i64. Returns `None` if ≥ 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x < 0).then(|| Self(x, PhantomData))
    }
}

impl I64<NonZero> {
    /// Create a non-zero i64. Returns `None` if = 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x != 0).then(|| Self(x, PhantomData))
    }
}

impl I64<NonNeg> {
    /// Create a non-negative i64. Returns `None` if < 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x >= 0).then(|| Self(x, PhantomData))
    }
}

impl I64<NonPos> {
    /// Create a non-positive i64. Returns `None` if > 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x <= 0).then(|| Self(x, PhantomData))
    }
}

impl I64<Zero> {
    /// The zero constant.
    pub const ZERO: Self = Self(0, PhantomData);

    /// Create a zero. Returns `None` if not exactly 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x == 0).then_some(Self::ZERO)
    }
}

impl I64<One> {
    /// The one constant.
    pub const ONE: Self = Self(1, PhantomData);

    /// Create a one. Returns `None` if not exactly 1.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x == 1).then_some(Self::ONE)
    }
}

impl I64<MinusOne> {
    /// The minus-one constant.
    pub const MINUS_ONE: Self = Self(-1, PhantomData);

    /// Create a minus-one. Returns `None` if not exactly -1.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x == -1).then_some(Self::MINUS_ONE)
    }
}

impl Default for I64<Zero> {
    fn default() -> Self {
        Self::ZERO
    }
}

// ============================================================================
// Trait Implementations
// ============================================================================

impl IsFinite for I64<Finite> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        self
    }
}

impl IsFinite for I64<Pos> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<Neg> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<Zero> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<One> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<MinusOne> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<NonNeg> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<NonZero> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

impl IsFinite for I64<NonPos> {
    type Finite = I64<Finite>;
    fn to_finite(self) -> Self::Finite {
        I64(self.0, PhantomData)
    }
}

// --- IsNonZero ---

impl IsNonZero for I64<NonZero> {
    type NonZero = I64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        self
    }
}

impl IsNonZero for I64<Pos> {
    type NonZero = I64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I64(self.0, PhantomData)
    }
}

impl IsNonZero for I64<Neg> {
    type NonZero = I64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I64(self.0, PhantomData)
    }
}

impl IsNonZero for I64<One> {
    type NonZero = I64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I64(self.0, PhantomData)
    }
}

impl IsNonZero for I64<MinusOne> {
    type NonZero = I64<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I64(self.0, PhantomData)
    }
}

// --- IsPositive ---

impl IsPositive for I64<Pos> {
    type Positive = I64<Pos>;
    fn to_positive(self) -> Self::Positive {
        self
    }
}

impl IsPositive for I64<One> {
    type Positive = I64<Pos>;
    fn to_positive(self) -> Self::Positive {
        I64(self.0, PhantomData)
    }
}

// --- IsNegative ---

impl IsNegative for I64<Neg> {
    type Negative = I64<Neg>;
    fn to_negative(self) -> Self::Negative {
        self
    }
}

impl IsNegative for I64<MinusOne> {
    type Negative = I64<Neg>;
    fn to_negative(self) -> Self::Negative {
        I64(self.0, PhantomData)
    }
}

// --- IsZero ---

impl IsZero for I64<Zero> {}

// --- IsNonNeg ---

impl IsNonNeg for I64<NonNeg> {
    type NonNeg = I64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        self
    }
}

impl IsNonNeg for I64<Pos> {
    type NonNeg = I64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I64(self.0, PhantomData)
    }
}

impl IsNonNeg for I64<One> {
    type NonNeg = I64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I64(self.0, PhantomData)
    }
}

impl IsNonNeg for I64<Zero> {
    type NonNeg = I64<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I64(self.0, PhantomData)
    }
}

// --- IsNonPos ---

impl IsNonPos for I64<NonPos> {
    type NonPos = I64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        self
    }
}

impl IsNonPos for I64<Neg> {
    type NonPos = I64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        I64(self.0, PhantomData)
    }
}

impl IsNonPos for I64<MinusOne> {
    type NonPos = I64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        I64(self.0, PhantomData)
    }
}

impl IsNonPos for I64<Zero> {
    type NonPos = I64<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        I64(self.0, PhantomData)
    }
}

// --- IsOne / IsMinusOne ---

impl IsOne for I64<One> {}
impl IsMinusOne for I64<MinusOne> {}

// ============================================================================
// Compile-time Macros
// ============================================================================

/// Create an `I64<Pos>` from a positive literal. Compile-time verified.
#[macro_export]
macro_rules! i64_pos {
    ($val:expr) => {{
        const V: i64 = $val;
        const _: () = assert!(V > 0, "i64_pos! requires positive value");
        $crate::types::i64::I64::<$crate::types::tags::Pos>::from_raw(V)
    }};
}

/// Create an `I64<Neg>` from a negative literal. Compile-time verified.
#[macro_export]
macro_rules! i64_neg {
    ($val:expr) => {{
        const V: i64 = $val;
        const _: () = assert!(V < 0, "i64_neg! requires negative value");
        $crate::types::i64::I64::<$crate::types::tags::Neg>::from_raw(V)
    }};
}

/// Create an `I64<NonNeg>` from a non-negative literal. Compile-time verified.
#[macro_export]
macro_rules! i64_nonneg {
    ($val:expr) => {{
        const V: i64 = $val;
        const _: () = assert!(V >= 0, "i64_nonneg! requires non-negative value");
        $crate::types::i64::I64::<$crate::types::tags::NonNeg>::from_raw(V)
    }};
}

/// Create an `I64<NonPos>` from a non-positive literal. Compile-time verified.
#[macro_export]
macro_rules! i64_nonpos {
    ($val:expr) => {{
        const V: i64 = $val;
        const _: () = assert!(V <= 0, "i64_nonpos! requires non-positive value");
        $crate::types::i64::I64::<$crate::types::tags::NonPos>::from_raw(V)
    }};
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_i64_size() {
        assert_eq!(std::mem::size_of::<I64>(), std::mem::size_of::<i64>());
        assert_eq!(std::mem::size_of::<I64<Pos>>(), std::mem::size_of::<i64>());
    }

    #[test]
    fn test_pos_new() {
        assert!(I64::<Pos>::new(1).is_some());
        assert!(I64::<Pos>::new(0).is_none());
        assert!(I64::<Pos>::new(-1).is_none());
    }

    #[test]
    fn test_neg_new() {
        assert!(I64::<Neg>::new(-1).is_some());
        assert!(I64::<Neg>::new(0).is_none());
        assert!(I64::<Neg>::new(1).is_none());
    }

    #[test]
    fn test_constants() {
        assert_eq!(I64::<Zero>::ZERO.get(), 0);
        assert_eq!(I64::<One>::ONE.get(), 1);
        assert_eq!(I64::<MinusOne>::MINUS_ONE.get(), -1);
    }

    #[test]
    fn test_macros() {
        let x = i64_pos!(42);
        assert_eq!(x.get(), 42);

        let y = i64_neg!(-5);
        assert_eq!(y.get(), -5);
    }

    #[test]
    fn test_ord() {
        let a = I64::<Pos>::new(1).unwrap();
        let b = I64::<Pos>::new(2).unwrap();
        assert!(a < b);
    }
}
