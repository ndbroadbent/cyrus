//! Type-safe i32 wrapper with compile-time invariant tracking.
//!
//! `I32<Tag>` wraps an `i32` with a phantom type that encodes constraints:
//!
//! - `I32` or `I32<Finite>`: Any value (always finite for integers)
//! - `I32<Pos>`: Strictly positive (> 0)
//! - `I32<Neg>`: Strictly negative (< 0)
//! - `I32<Zero>`: Exactly zero
//! - `I32<One>`: Exactly 1
//! - `I32<MinusOne>`: Exactly -1
//! - `I32<NonNeg>`: Non-negative (≥ 0)
//! - `I32<NonPos>`: Non-positive (≤ 0)
//!
//! # Example
//!
//! ```
//! use cyrus_core::types::i32::I32;
//! use cyrus_core::types::tags::Pos;
//!
//! let h11 = I32::<Pos>::new(4).unwrap();  // Hodge number
//! assert!(h11.get() > 0);
//! ```

use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

pub use super::tags::{
    Finite, GTEOne, IsFinite, IsGTEOne, IsMinusOne, IsNegative, IsNonNeg, IsNonPos, IsNonZero,
    IsOne, IsPositive, IsTwo, IsZero, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
};

// ============================================================================
// I32<Tag> Wrapper
// ============================================================================

/// A type-safe i32 wrapper with compile-time invariant tracking.
///
/// The default tag is `Finite`, so `I32` means `I32<Finite>`.
/// For integers, `Finite` is always satisfied (no NaN/infinity).
///
/// # Representation
///
/// `#[repr(transparent)]` ensures this has the same memory layout as `i32`.
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct I32<Tag = Finite>(pub(crate) i32, pub(crate) PhantomData<Tag>);

// ============================================================================
// Core Impls
// ============================================================================

impl<Tag> I32<Tag> {
    /// Get the underlying i32 value.
    #[inline(always)]
    #[must_use]
    pub fn get(self) -> i32 {
        self.0
    }

    /// Create from a raw i32 without validation.
    ///
    /// Used internally and by macros after compile-time verification.
    #[doc(hidden)]
    #[inline(always)]
    #[must_use]
    pub const fn from_raw(x: i32) -> Self {
        Self(x, PhantomData)
    }
}

impl<Tag> std::fmt::Debug for I32<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let full_name = std::any::type_name::<Tag>();
        let tag_name = full_name.rsplit("::").next().unwrap_or(full_name);
        write!(f, "I32<{tag_name}>({})", self.0)
    }
}

impl<Tag> std::fmt::Display for I32<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<Tag> PartialEq for I32<Tag> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<Tag> Eq for I32<Tag> {}

impl<Tag> Hash for I32<Tag> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl<Tag> PartialOrd for I32<Tag> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Tag> Ord for I32<Tag> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

// ============================================================================
// Constructors
// ============================================================================

impl I32<Finite> {
    /// Create a finite i32. Always succeeds for integers.
    #[must_use]
    pub fn new(x: i32) -> Self {
        Self(x, PhantomData)
    }
}

impl I32<Pos> {
    /// Create a positive i32. Returns `None` if ≤ 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x > 0).then(|| Self(x, PhantomData))
    }
}

impl I32<Neg> {
    /// Create a negative i32. Returns `None` if ≥ 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x < 0).then(|| Self(x, PhantomData))
    }
}

impl I32<NonZero> {
    /// Create a non-zero i32. Returns `None` if = 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x != 0).then(|| Self(x, PhantomData))
    }
}

impl I32<NonNeg> {
    /// Create a non-negative i32. Returns `None` if < 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x >= 0).then(|| Self(x, PhantomData))
    }
}

impl I32<NonPos> {
    /// Create a non-positive i32. Returns `None` if > 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x <= 0).then(|| Self(x, PhantomData))
    }
}

impl I32<Zero> {
    /// The zero constant.
    pub const ZERO: Self = Self(0, PhantomData);

    /// Create a zero. Returns `None` if not exactly 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x == 0).then_some(Self::ZERO)
    }
}

impl I32<One> {
    /// The one constant.
    pub const ONE: Self = Self(1, PhantomData);

    /// Create a one. Returns `None` if not exactly 1.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x == 1).then_some(Self::ONE)
    }
}

impl I32<Two> {
    /// The two constant.
    pub const TWO: Self = Self(2, PhantomData);

    /// Create a two. Returns `None` if not exactly 2.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x == 2).then_some(Self::TWO)
    }
}

impl I32<GTEOne> {
    /// Create a value ≥ 1. Returns `None` if < 1.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x >= 1).then(|| Self(x, PhantomData))
    }
}

impl I32<MinusOne> {
    /// The minus-one constant.
    pub const MINUS_ONE: Self = Self(-1, PhantomData);

    /// Create a minus-one. Returns `None` if not exactly -1.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x == -1).then_some(Self::MINUS_ONE)
    }
}

impl Default for I32<Zero> {
    fn default() -> Self {
        Self::ZERO
    }
}

// ============================================================================
// Trait Implementations
// ============================================================================

impl IsFinite for I32<Finite> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        self
    }
}

impl IsFinite for I32<Pos> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<Neg> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<Zero> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<One> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<Two> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<MinusOne> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<NonZero> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<NonNeg> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<NonPos> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

impl IsFinite for I32<GTEOne> {
    type Finite = I32<Finite>;
    fn to_finite(self) -> Self::Finite {
        I32(self.0, PhantomData)
    }
}

// --- IsNonZero ---

impl IsNonZero for I32<NonZero> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        self
    }
}

impl IsNonZero for I32<Pos> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I32(self.0, PhantomData)
    }
}

impl IsNonZero for I32<Neg> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I32(self.0, PhantomData)
    }
}

impl IsNonZero for I32<One> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I32(self.0, PhantomData)
    }
}

impl IsNonZero for I32<Two> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I32(self.0, PhantomData)
    }
}

impl IsNonZero for I32<MinusOne> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I32(self.0, PhantomData)
    }
}

impl IsNonZero for I32<GTEOne> {
    type NonZero = I32<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        I32(self.0, PhantomData)
    }
}

// --- IsPositive ---

impl IsPositive for I32<Pos> {
    type Positive = I32<Pos>;
    fn to_positive(self) -> Self::Positive {
        self
    }
}

impl IsPositive for I32<One> {
    type Positive = I32<Pos>;
    fn to_positive(self) -> Self::Positive {
        I32(self.0, PhantomData)
    }
}

impl IsPositive for I32<Two> {
    type Positive = I32<Pos>;
    fn to_positive(self) -> Self::Positive {
        I32(self.0, PhantomData)
    }
}

impl IsPositive for I32<GTEOne> {
    type Positive = I32<Pos>;
    fn to_positive(self) -> Self::Positive {
        I32(self.0, PhantomData)
    }
}

// --- IsNegative ---

impl IsNegative for I32<Neg> {
    type Negative = I32<Neg>;
    fn to_negative(self) -> Self::Negative {
        self
    }
}

impl IsNegative for I32<MinusOne> {
    type Negative = I32<Neg>;
    fn to_negative(self) -> Self::Negative {
        I32(self.0, PhantomData)
    }
}

// --- IsZero ---

impl IsZero for I32<Zero> {}

// --- IsNonNeg ---

impl IsNonNeg for I32<NonNeg> {
    type NonNeg = I32<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        self
    }
}

impl IsNonNeg for I32<Pos> {
    type NonNeg = I32<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I32(self.0, PhantomData)
    }
}

impl IsNonNeg for I32<One> {
    type NonNeg = I32<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I32(self.0, PhantomData)
    }
}

impl IsNonNeg for I32<Two> {
    type NonNeg = I32<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I32(self.0, PhantomData)
    }
}

impl IsNonNeg for I32<Zero> {
    type NonNeg = I32<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I32(self.0, PhantomData)
    }
}

impl IsNonNeg for I32<GTEOne> {
    type NonNeg = I32<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        I32(self.0, PhantomData)
    }
}

// --- IsNonPos ---

impl IsNonPos for I32<NonPos> {
    type NonPos = I32<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        self
    }
}

impl IsNonPos for I32<Neg> {
    type NonPos = I32<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        I32(self.0, PhantomData)
    }
}

impl IsNonPos for I32<MinusOne> {
    type NonPos = I32<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        I32(self.0, PhantomData)
    }
}

impl IsNonPos for I32<Zero> {
    type NonPos = I32<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        I32(self.0, PhantomData)
    }
}

// --- IsOne / IsMinusOne ---

impl IsOne for I32<One> {}
impl IsTwo for I32<Two> {}
impl IsMinusOne for I32<MinusOne> {}
impl IsGTEOne for I32<GTEOne> {}
impl IsGTEOne for I32<One> {}
impl IsGTEOne for I32<Two> {}

// ============================================================================
// Arithmetic Operations (generated from algebra rules)
// ============================================================================

crate::impl_ops!(I32, i32);

// ============================================================================
// Compile-time Macros
// ============================================================================

/// Create an `I32<Pos>` from a positive literal. Compile-time verified.
#[macro_export]
macro_rules! i32_pos {
    ($val:expr) => {{
        const V: i32 = $val;
        const _: () = assert!(V > 0, "i32_pos! requires positive value");
        $crate::types::i32::I32::<$crate::types::tags::Pos>::from_raw(V)
    }};
}

/// Create an `I32<Neg>` from a negative literal. Compile-time verified.
#[macro_export]
macro_rules! i32_neg {
    ($val:expr) => {{
        const V: i32 = $val;
        const _: () = assert!(V < 0, "i32_neg! requires negative value");
        $crate::types::i32::I32::<$crate::types::tags::Neg>::from_raw(V)
    }};
}

/// Create an `I32<NonNeg>` from a non-negative literal. Compile-time verified.
#[macro_export]
macro_rules! i32_nonneg {
    ($val:expr) => {{
        const V: i32 = $val;
        const _: () = assert!(V >= 0, "i32_nonneg! requires non-negative value");
        $crate::types::i32::I32::<$crate::types::tags::NonNeg>::from_raw(V)
    }};
}

/// Create an `I32<NonPos>` from a non-positive literal. Compile-time verified.
#[macro_export]
macro_rules! i32_nonpos {
    ($val:expr) => {{
        const V: i32 = $val;
        const _: () = assert!(V <= 0, "i32_nonpos! requires non-positive value");
        $crate::types::i32::I32::<$crate::types::tags::NonPos>::from_raw(V)
    }};
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_i32_size() {
        assert_eq!(std::mem::size_of::<I32>(), std::mem::size_of::<i32>());
        assert_eq!(std::mem::size_of::<I32<Pos>>(), std::mem::size_of::<i32>());
    }

    #[test]
    fn test_pos_new() {
        assert!(I32::<Pos>::new(1).is_some());
        assert!(I32::<Pos>::new(0).is_none());
        assert!(I32::<Pos>::new(-1).is_none());
    }

    #[test]
    fn test_neg_new() {
        assert!(I32::<Neg>::new(-1).is_some());
        assert!(I32::<Neg>::new(0).is_none());
        assert!(I32::<Neg>::new(1).is_none());
    }

    #[test]
    fn test_constants() {
        assert_eq!(I32::<Zero>::ZERO.get(), 0);
        assert_eq!(I32::<One>::ONE.get(), 1);
        assert_eq!(I32::<MinusOne>::MINUS_ONE.get(), -1);
    }

    #[test]
    fn test_macros() {
        let x = i32_pos!(42);
        assert_eq!(x.get(), 42);

        let y = i32_neg!(-5);
        assert_eq!(y.get(), -5);
    }

    #[test]
    fn test_ord() {
        let a = I32::<Pos>::new(1).unwrap();
        let b = I32::<Pos>::new(2).unwrap();
        assert!(a < b);
    }
}
