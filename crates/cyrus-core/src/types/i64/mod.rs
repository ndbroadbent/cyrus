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
    /// The zero constant.
    pub const ZERO: Self = Self(0, PhantomData);

    /// Create a finite i64. Always succeeds for integers.
    #[must_use]
    pub fn new(x: i64) -> Self {
        Self(x, PhantomData)
    }

    /// Convert to F64<Finite>.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<Finite> {
        crate::types::f64::F64::from_raw(self.0 as f64)
    }
}

impl I64<Pos> {
    /// Create a positive i64. Returns `None` if ≤ 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x > 0).then(|| Self(x, PhantomData))
    }

    /// Convert to F64<Pos>. Preserves positivity.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<Pos> {
        crate::types::f64::F64::from_raw(self.0 as f64)
    }
}

impl I64<Neg> {
    /// Create a negative i64. Returns `None` if ≥ 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x < 0).then(|| Self(x, PhantomData))
    }

    /// Convert to F64<Neg>. Preserves negativity.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<Neg> {
        crate::types::f64::F64::from_raw(self.0 as f64)
    }
}

impl I64<NonZero> {
    /// Create a non-zero i64. Returns `None` if = 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x != 0).then(|| Self(x, PhantomData))
    }

    /// Convert to F64<NonZero>. Preserves non-zero.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<NonZero> {
        crate::types::f64::F64::from_raw(self.0 as f64)
    }
}

impl I64<NonNeg> {
    /// Create a non-negative i64. Returns `None` if < 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x >= 0).then(|| Self(x, PhantomData))
    }

    /// Convert to F64<NonNeg>. Preserves non-negativity.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<NonNeg> {
        crate::types::f64::F64::from_raw(self.0 as f64)
    }
}

impl I64<NonPos> {
    /// Create a non-positive i64. Returns `None` if > 0.
    #[must_use]
    pub fn new(x: i64) -> Option<Self> {
        (x <= 0).then(|| Self(x, PhantomData))
    }

    /// Convert to F64<NonPos>. Preserves non-positivity.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<NonPos> {
        crate::types::f64::F64::from_raw(self.0 as f64)
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

    /// Convert to F64<Zero>.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<Zero> {
        crate::types::f64::F64::from_raw(0.0)
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

    /// Convert to F64<One>.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<One> {
        crate::types::f64::F64::from_raw(1.0)
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

    /// Convert to F64<MinusOne>.
    #[must_use]
    pub fn to_f64(self) -> crate::types::f64::F64<MinusOne> {
        crate::types::f64::F64::from_raw(-1.0)
    }
}

impl Default for I64<Zero> {
    fn default() -> Self {
        Self::ZERO
    }
}

// ============================================================================
// Marker Trait Implementations
// ============================================================================

impl IsFinite for I64<Finite> {}
impl IsFinite for I64<Pos> {}
impl IsFinite for I64<Neg> {}
impl IsFinite for I64<Zero> {}
impl IsFinite for I64<One> {}
impl IsFinite for I64<MinusOne> {}
impl IsFinite for I64<NonNeg> {}
impl IsFinite for I64<NonZero> {}
impl IsFinite for I64<NonPos> {}

impl IsNonZero for I64<NonZero> {}
impl IsNonZero for I64<Pos> {}
impl IsNonZero for I64<Neg> {}
impl IsNonZero for I64<One> {}
impl IsNonZero for I64<MinusOne> {}

impl IsPositive for I64<Pos> {}
impl IsPositive for I64<One> {}

impl IsNegative for I64<Neg> {}
impl IsNegative for I64<MinusOne> {}

impl IsZero for I64<Zero> {}

impl IsNonNeg for I64<NonNeg> {}
impl IsNonNeg for I64<Pos> {}
impl IsNonNeg for I64<One> {}
impl IsNonNeg for I64<Zero> {}

impl IsNonPos for I64<NonPos> {}
impl IsNonPos for I64<Neg> {}
impl IsNonPos for I64<MinusOne> {}
impl IsNonPos for I64<Zero> {}

impl IsOne for I64<One> {}
impl IsMinusOne for I64<MinusOne> {}

// ============================================================================
// Arithmetic Operations (generated from algebra rules)
// ============================================================================

crate::impl_ops!(I64, i64);

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

#[cfg(test)]
mod tests;
