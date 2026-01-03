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
    #[inline]
    #[must_use]
    pub const fn get(self) -> i32 {
        self.0
    }

    /// Create from a raw i32 without validation.
    ///
    /// Used internally and by macros after compile-time verification.
    #[doc(hidden)]
    #[inline]
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
    pub const fn new(x: i32) -> Self {
        Self(x, PhantomData)
    }
}

impl I32<Pos> {
    /// Create a positive i32. Returns `None` if ≤ 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x > 0).then_some(Self(x, PhantomData))
    }
}

impl I32<Neg> {
    /// Create a negative i32. Returns `None` if ≥ 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x < 0).then_some(Self(x, PhantomData))
    }
}

impl I32<NonZero> {
    /// Create a non-zero i32. Returns `None` if = 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x != 0).then_some(Self(x, PhantomData))
    }
}

impl I32<NonNeg> {
    /// Create a non-negative i32. Returns `None` if < 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x >= 0).then_some(Self(x, PhantomData))
    }
}

impl I32<NonPos> {
    /// Create a non-positive i32. Returns `None` if > 0.
    #[must_use]
    pub fn new(x: i32) -> Option<Self> {
        (x <= 0).then_some(Self(x, PhantomData))
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
        (x >= 1).then_some(Self(x, PhantomData))
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
// Marker Trait Implementations
// ============================================================================

impl IsFinite for I32<Finite> {}
impl IsFinite for I32<Pos> {}
impl IsFinite for I32<Neg> {}
impl IsFinite for I32<Zero> {}
impl IsFinite for I32<One> {}
impl IsFinite for I32<Two> {}
impl IsFinite for I32<MinusOne> {}
impl IsFinite for I32<NonZero> {}
impl IsFinite for I32<NonNeg> {}
impl IsFinite for I32<NonPos> {}
impl IsFinite for I32<GTEOne> {}

impl IsNonZero for I32<NonZero> {}
impl IsNonZero for I32<Pos> {}
impl IsNonZero for I32<Neg> {}
impl IsNonZero for I32<One> {}
impl IsNonZero for I32<Two> {}
impl IsNonZero for I32<MinusOne> {}
impl IsNonZero for I32<GTEOne> {}

impl IsPositive for I32<Pos> {}
impl IsPositive for I32<One> {}
impl IsPositive for I32<Two> {}
impl IsPositive for I32<GTEOne> {}

impl IsNegative for I32<Neg> {}
impl IsNegative for I32<MinusOne> {}

impl IsZero for I32<Zero> {}

impl IsNonNeg for I32<NonNeg> {}
impl IsNonNeg for I32<Pos> {}
impl IsNonNeg for I32<One> {}
impl IsNonNeg for I32<Two> {}
impl IsNonNeg for I32<Zero> {}
impl IsNonNeg for I32<GTEOne> {}

impl IsNonPos for I32<NonPos> {}
impl IsNonPos for I32<Neg> {}
impl IsNonPos for I32<MinusOne> {}
impl IsNonPos for I32<Zero> {}

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
    use std::collections::hash_map::DefaultHasher;

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
    fn test_nonzero_new() {
        assert!(I32::<NonZero>::new(1).is_some());
        assert!(I32::<NonZero>::new(-1).is_some());
        assert!(I32::<NonZero>::new(0).is_none());
    }

    #[test]
    fn test_nonneg_new() {
        assert!(I32::<NonNeg>::new(0).is_some());
        assert!(I32::<NonNeg>::new(1).is_some());
        assert!(I32::<NonNeg>::new(-1).is_none());
    }

    #[test]
    fn test_nonpos_new() {
        assert!(I32::<NonPos>::new(0).is_some());
        assert!(I32::<NonPos>::new(-1).is_some());
        assert!(I32::<NonPos>::new(1).is_none());
    }

    #[test]
    fn test_zero_new() {
        assert!(I32::<Zero>::new(0).is_some());
        assert!(I32::<Zero>::new(1).is_none());
        assert!(I32::<Zero>::new(-1).is_none());
    }

    #[test]
    fn test_one_new() {
        assert!(I32::<One>::new(1).is_some());
        assert!(I32::<One>::new(0).is_none());
        assert!(I32::<One>::new(2).is_none());
    }

    #[test]
    fn test_two_new() {
        assert!(I32::<Two>::new(2).is_some());
        assert!(I32::<Two>::new(1).is_none());
        assert!(I32::<Two>::new(0).is_none());
    }

    #[test]
    fn test_gteone_new() {
        assert!(I32::<GTEOne>::new(1).is_some());
        assert!(I32::<GTEOne>::new(2).is_some());
        assert!(I32::<GTEOne>::new(100).is_some());
        assert!(I32::<GTEOne>::new(0).is_none());
        assert!(I32::<GTEOne>::new(-1).is_none());
    }

    #[test]
    fn test_minusone_new() {
        assert!(I32::<MinusOne>::new(-1).is_some());
        assert!(I32::<MinusOne>::new(0).is_none());
        assert!(I32::<MinusOne>::new(-2).is_none());
    }

    #[test]
    fn test_finite_new() {
        let x = I32::<Finite>::new(42);
        assert_eq!(x.get(), 42);
        let y = I32::<Finite>::new(-100);
        assert_eq!(y.get(), -100);
    }

    #[test]
    fn test_constants() {
        assert_eq!(I32::<Zero>::ZERO.get(), 0);
        assert_eq!(I32::<One>::ONE.get(), 1);
        assert_eq!(I32::<Two>::TWO.get(), 2);
        assert_eq!(I32::<MinusOne>::MINUS_ONE.get(), -1);
    }

    #[test]
    fn test_debug_display() {
        let x = I32::<Pos>::new(42).unwrap();
        assert!(format!("{x:?}").contains("Pos"));
        assert!(format!("{x:?}").contains("42"));
        assert_eq!(format!("{x}"), "42");

        let y = I32::<Finite>::new(-5);
        assert!(format!("{y:?}").contains("Finite"));
    }

    #[test]
    fn test_eq_hash() {
        let a = I32::<Pos>::new(5).unwrap();
        let b = I32::<Pos>::new(5).unwrap();
        let c = I32::<Pos>::new(6).unwrap();
        assert_eq!(a, b);
        assert_ne!(a, c);

        let mut h1 = DefaultHasher::new();
        let mut h2 = DefaultHasher::new();
        a.hash(&mut h1);
        b.hash(&mut h2);
        assert_eq!(h1.finish(), h2.finish());
    }

    #[test]
    fn test_default_zero() {
        let z: I32<Zero> = I32::default();
        assert_eq!(z.get(), 0);
    }

    #[test]
    fn test_macros() {
        let x = i32_pos!(42);
        assert_eq!(x.get(), 42);

        let y = i32_neg!(-5);
        assert_eq!(y.get(), -5);

        let z = i32_nonneg!(0);
        assert_eq!(z.get(), 0);
        let z2 = i32_nonneg!(10);
        assert_eq!(z2.get(), 10);

        let w = i32_nonpos!(0);
        assert_eq!(w.get(), 0);
        let w2 = i32_nonpos!(-10);
        assert_eq!(w2.get(), -10);
    }

    #[test]
    fn test_ord() {
        let a = I32::<Pos>::new(1).unwrap();
        let b = I32::<Pos>::new(2).unwrap();
        assert!(a < b);
        assert!(b > a);
        assert!(a <= b);
        assert!(b >= a);
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
    }

    #[test]
    fn test_from_raw() {
        let x = I32::<Pos>::from_raw(100);
        assert_eq!(x.get(), 100);
    }
}
