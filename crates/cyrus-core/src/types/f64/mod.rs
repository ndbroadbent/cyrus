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
//! use cyrus_core::{f64_pos, f64_finite};
//! use cyrus_core::types::f64::F64;
//!
//! let x = f64_pos!(4.0);     // Compile-time verified F64<Pos>
//! let y = x * x;             // Pos * Pos = Pos
//! assert!((y.get() - 16.0).abs() < 1e-10);
//!
//! let z = f64_finite!(-3.5); // Compile-time verified F64 (F64<Finite>)
//! ```

use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

pub use super::tags::{
    Finite, GTEOne, IsFinite, IsGTEOne, IsMinusOne, IsNegative, IsNonNeg, IsNonPos, IsNonZero,
    IsOne, IsPositive, IsTwo, IsZero, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
    tag,
};

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
    #[inline]
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }

    /// Create from a raw f64 without validation.
    ///
    /// Used internally and by macros after compile-time verification.
    #[doc(hidden)]
    #[inline]
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
        x.is_finite().then_some(Self(x, PhantomData))
    }

    /// Zero constant.
    pub const ZERO: Self = Self(0.0, PhantomData);

    /// Compute absolute value. Finite.abs() is NonNeg (could be zero).
    #[must_use]
    pub const fn abs(self) -> F64<NonNeg> {
        F64(self.0.abs(), PhantomData)
    }

    /// Compute square. x² is always non-negative.
    #[must_use]
    pub fn square(self) -> F64<NonNeg> {
        F64(self.0 * self.0, PhantomData)
    }

    /// Try to narrow to positive. Returns `None` if ≤ 0.
    #[must_use]
    pub fn try_to_pos(self) -> Option<F64<Pos>> {
        (self.0 > 0.0).then_some(F64(self.0, PhantomData))
    }

    /// Try to narrow to negative. Returns `None` if ≥ 0.
    #[must_use]
    pub fn try_to_neg(self) -> Option<F64<Neg>> {
        (self.0 < 0.0).then_some(F64(self.0, PhantomData))
    }

    /// Try to narrow to non-zero. Returns `None` if = 0.
    #[must_use]
    pub fn try_to_non_zero(self) -> Option<F64<NonZero>> {
        (self.0 != 0.0).then_some(F64(self.0, PhantomData))
    }
}

impl F64<Pos> {
    /// Create a positive f64. Returns `None` if ≤ 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x > 0.0).then_some(Self(x, PhantomData))
    }

    /// Compute absolute value. Pos.abs() is Pos (already positive).
    #[must_use]
    pub const fn abs(self) -> Self {
        self
    }

    /// Compute square root. sqrt(Pos) is Pos.
    #[must_use]
    pub fn sqrt(self) -> Self {
        Self(self.0.sqrt(), PhantomData)
    }

    /// Compute natural logarithm. ln(Pos) is Finite (any real number).
    #[must_use]
    pub fn ln(self) -> F64<Finite> {
        F64(self.0.ln(), PhantomData)
    }

    /// Compute reciprocal. 1/Pos is Pos.
    #[must_use]
    pub fn recip(self) -> Self {
        Self(1.0 / self.0, PhantomData)
    }
}

impl F64<Neg> {
    /// Create a negative f64. Returns `None` if ≥ 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x < 0.0).then_some(Self(x, PhantomData))
    }

    /// Compute absolute value. Neg.abs() is Pos.
    #[must_use]
    pub const fn abs(self) -> F64<Pos> {
        F64(self.0.abs(), PhantomData)
    }
}

impl F64<NonZero> {
    /// Create a non-zero f64. Returns `None` if = 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x != 0.0).then_some(Self(x, PhantomData))
    }

    /// Compute absolute value. NonZero.abs() is Pos (not just NonNeg).
    #[must_use]
    pub const fn abs(self) -> F64<Pos> {
        F64(self.0.abs(), PhantomData)
    }
}

impl F64<NonNeg> {
    /// The zero constant.
    pub const ZERO: Self = Self(0.0, PhantomData);

    /// Create a non-negative f64. Returns `None` if < 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x >= 0.0).then_some(Self(x, PhantomData))
    }

    /// Compute square root. sqrt(NonNeg) is NonNeg.
    #[must_use]
    pub fn sqrt(self) -> Self {
        Self(self.0.sqrt(), PhantomData)
    }
}

impl F64<NonPos> {
    /// Create a non-positive f64. Returns `None` if > 0, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x <= 0.0).then_some(Self(x, PhantomData))
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
    #[allow(clippy::float_cmp)] // Intentional exact comparison for special value
    pub fn new(x: f64) -> Option<Self> {
        (x == 1.0).then_some(Self::ONE)
    }
}

impl F64<MinusOne> {
    /// The minus-one constant.
    pub const MINUS_ONE: Self = Self(-1.0, PhantomData);

    /// Create a minus-one. Returns `None` if not exactly -1.0.
    #[must_use]
    #[allow(clippy::float_cmp)] // Intentional exact comparison for special value
    pub fn new(x: f64) -> Option<Self> {
        (x == -1.0).then_some(Self::MINUS_ONE)
    }
}

impl F64<Two> {
    /// The two constant.
    pub const TWO: Self = Self(2.0, PhantomData);

    /// Create a two. Returns `None` if not exactly 2.0.
    #[must_use]
    #[allow(clippy::float_cmp)] // Intentional exact comparison for special value
    pub fn new(x: f64) -> Option<Self> {
        (x == 2.0).then_some(Self::TWO)
    }
}

impl F64<GTEOne> {
    /// Create a value ≥ 1. Returns `None` if < 1, NaN, or infinite.
    #[must_use]
    pub fn new(x: f64) -> Option<Self> {
        (x.is_finite() && x >= 1.0).then_some(Self(x, PhantomData))
    }
}

impl Default for F64<Zero> {
    fn default() -> Self {
        Self::ZERO
    }
}

// ============================================================================
// Marker Trait Implementations
// ============================================================================

// IsFinite - all tags are finite
impl IsFinite for F64<Finite> {}
impl IsFinite for F64<Pos> {}
impl IsFinite for F64<Neg> {}
impl IsFinite for F64<Zero> {}
impl IsFinite for F64<One> {}
impl IsFinite for F64<MinusOne> {}
impl IsFinite for F64<NonZero> {}
impl IsFinite for F64<NonNeg> {}
impl IsFinite for F64<NonPos> {}
impl IsFinite for F64<Two> {}
impl IsFinite for F64<GTEOne> {}

// IsNonZero
impl IsNonZero for F64<NonZero> {}
impl IsNonZero for F64<Pos> {}
impl IsNonZero for F64<Neg> {}
impl IsNonZero for F64<One> {}
impl IsNonZero for F64<MinusOne> {}
impl IsNonZero for F64<Two> {}
impl IsNonZero for F64<GTEOne> {}

// IsPositive
impl IsPositive for F64<Pos> {}
impl IsPositive for F64<One> {}
impl IsPositive for F64<Two> {}
impl IsPositive for F64<GTEOne> {}

// IsNegative
impl IsNegative for F64<Neg> {}
impl IsNegative for F64<MinusOne> {}

// IsZero
impl IsZero for F64<Zero> {}

// IsNonNeg
impl IsNonNeg for F64<NonNeg> {}
impl IsNonNeg for F64<Pos> {}
impl IsNonNeg for F64<One> {}
impl IsNonNeg for F64<Zero> {}
impl IsNonNeg for F64<Two> {}
impl IsNonNeg for F64<GTEOne> {}

// IsNonPos
impl IsNonPos for F64<NonPos> {}
impl IsNonPos for F64<Neg> {}
impl IsNonPos for F64<MinusOne> {}
impl IsNonPos for F64<Zero> {}

// Exact values
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

#[cfg(test)]
mod tests;
