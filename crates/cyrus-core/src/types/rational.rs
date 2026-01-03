//! Type-safe malachite Rational wrapper with compile-time invariant tracking.
//!
//! `Rational<Tag>` wraps a `malachite::Rational` with a phantom type that encodes constraints:
//!
//! - `Rational` or `Rational<Finite>`: Any finite rational (always true for rationals)
//! - `Rational<Pos>`: Strictly positive (> 0)
//! - `Rational<Neg>`: Strictly negative (< 0)
//! - `Rational<Zero>`: Exactly zero
//! - `Rational<NonNeg>`: Non-negative (≥ 0)
//! - `Rational<NonPos>`: Non-positive (≤ 0)

use malachite::Rational as MalachiteRational;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

pub use super::tags::{
    Finite, IsFinite, IsNegative, IsNonNeg, IsNonPos, IsNonZero, IsPositive, IsZero, Neg, NonNeg,
    NonPos, NonZero, Pos, Zero,
};

// ============================================================================
// Rational<Tag> Wrapper
// ============================================================================

/// A type-safe Rational wrapper with compile-time invariant tracking.
///
/// The default tag is `Finite`, so `Rational` means `Rational<Finite>`.
/// For rationals, `Finite` is always satisfied (no NaN/infinity).
#[derive(Clone)]
pub struct Rational<Tag = Finite>(pub(crate) MalachiteRational, pub(crate) PhantomData<Tag>);

// ============================================================================
// Core Impls
// ============================================================================

impl<Tag> Rational<Tag> {
    /// Get the underlying malachite Rational value.
    #[inline(always)]
    #[must_use]
    pub fn get(&self) -> &MalachiteRational {
        &self.0
    }

    /// Consume and return the inner malachite Rational.
    #[inline(always)]
    #[must_use]
    pub fn into_inner(self) -> MalachiteRational {
        self.0
    }

    /// Create from a raw Rational without validation.
    #[doc(hidden)]
    #[inline(always)]
    #[must_use]
    pub fn from_raw(x: MalachiteRational) -> Self {
        Self(x, PhantomData)
    }
}

impl<Tag> std::fmt::Debug for Rational<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let full_name = std::any::type_name::<Tag>();
        let tag_name = full_name.rsplit("::").next().unwrap_or(full_name);
        write!(f, "Rational<{tag_name}>({})", self.0)
    }
}

impl<Tag> std::fmt::Display for Rational<Tag> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<Tag> PartialEq for Rational<Tag> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<Tag> Eq for Rational<Tag> {}

impl<Tag> Hash for Rational<Tag> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl<Tag> PartialOrd for Rational<Tag> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Tag> Ord for Rational<Tag> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

// ============================================================================
// Constructors
// ============================================================================

impl Rational<Finite> {
    /// Create a finite Rational. Always succeeds for rationals.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Self {
        Self(x, PhantomData)
    }
}

impl Rational<Pos> {
    /// Create a positive Rational. Returns `None` if ≤ 0.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Option<Self> {
        use malachite::num::basic::traits::Zero;
        (x > MalachiteRational::ZERO).then(|| Self(x, PhantomData))
    }
}

impl Rational<Neg> {
    /// Create a negative Rational. Returns `None` if ≥ 0.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Option<Self> {
        use malachite::num::basic::traits::Zero;
        (x < MalachiteRational::ZERO).then(|| Self(x, PhantomData))
    }
}

impl Rational<NonZero> {
    /// Create a non-zero Rational. Returns `None` if = 0.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Option<Self> {
        use malachite::num::basic::traits::Zero;
        (x != MalachiteRational::ZERO).then(|| Self(x, PhantomData))
    }
}

impl Rational<NonNeg> {
    /// Create a non-negative Rational. Returns `None` if < 0.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Option<Self> {
        use malachite::num::basic::traits::Zero;
        (x >= MalachiteRational::ZERO).then(|| Self(x, PhantomData))
    }
}

impl Rational<NonPos> {
    /// Create a non-positive Rational. Returns `None` if > 0.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Option<Self> {
        use malachite::num::basic::traits::Zero;
        (x <= MalachiteRational::ZERO).then(|| Self(x, PhantomData))
    }
}

impl Rational<Zero> {
    /// Create a zero. Returns `None` if not exactly 0.
    #[must_use]
    pub fn new(x: MalachiteRational) -> Option<Self> {
        use malachite::num::basic::traits::Zero;
        (x == MalachiteRational::ZERO).then(|| Self(MalachiteRational::ZERO, PhantomData))
    }

    /// The zero constant.
    #[must_use]
    pub fn zero() -> Self {
        use malachite::num::basic::traits::Zero;
        Self(MalachiteRational::ZERO, PhantomData)
    }
}

impl Default for Rational<Zero> {
    fn default() -> Self {
        Self::zero()
    }
}

// ============================================================================
// Trait Implementations
// ============================================================================

impl IsFinite for Rational<Finite> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        self
    }
}

impl IsFinite for Rational<Pos> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        Rational(self.0, PhantomData)
    }
}

impl IsFinite for Rational<Neg> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        Rational(self.0, PhantomData)
    }
}

impl IsFinite for Rational<Zero> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        Rational(self.0, PhantomData)
    }
}

impl IsFinite for Rational<NonNeg> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        Rational(self.0, PhantomData)
    }
}

impl IsFinite for Rational<NonZero> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        Rational(self.0, PhantomData)
    }
}

impl IsFinite for Rational<NonPos> {
    type Finite = Rational<Finite>;
    fn to_finite(self) -> Self::Finite {
        Rational(self.0, PhantomData)
    }
}

// --- IsNonZero ---

impl IsNonZero for Rational<NonZero> {
    type NonZero = Rational<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        self
    }
}

impl IsNonZero for Rational<Pos> {
    type NonZero = Rational<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        Rational(self.0, PhantomData)
    }
}

impl IsNonZero for Rational<Neg> {
    type NonZero = Rational<NonZero>;
    fn to_non_zero(self) -> Self::NonZero {
        Rational(self.0, PhantomData)
    }
}

// --- IsPositive ---

impl IsPositive for Rational<Pos> {
    type Positive = Rational<Pos>;
    fn to_positive(self) -> Self::Positive {
        self
    }
}

// --- IsNegative ---

impl IsNegative for Rational<Neg> {
    type Negative = Rational<Neg>;
    fn to_negative(self) -> Self::Negative {
        self
    }
}

// --- IsZero ---

impl IsZero for Rational<Zero> {}

// --- IsNonNeg ---

impl IsNonNeg for Rational<NonNeg> {
    type NonNeg = Rational<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        self
    }
}

impl IsNonNeg for Rational<Pos> {
    type NonNeg = Rational<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        Rational(self.0, PhantomData)
    }
}

impl IsNonNeg for Rational<Zero> {
    type NonNeg = Rational<NonNeg>;
    fn to_non_neg(self) -> Self::NonNeg {
        Rational(self.0, PhantomData)
    }
}

// --- IsNonPos ---

impl IsNonPos for Rational<NonPos> {
    type NonPos = Rational<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        self
    }
}

impl IsNonPos for Rational<Neg> {
    type NonPos = Rational<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        Rational(self.0, PhantomData)
    }
}

impl IsNonPos for Rational<Zero> {
    type NonPos = Rational<NonPos>;
    fn to_non_pos(self) -> Self::NonPos {
        Rational(self.0, PhantomData)
    }
}

// ============================================================================
// Conversions
// ============================================================================

impl Rational<Pos> {
    /// Convert to f64. Positive rational → positive f64.
    #[must_use]
    pub fn to_f64(&self) -> crate::types::f64::F64<Pos> {
        use malachite::num::conversion::traits::RoundingFrom;
        use malachite::rounding_modes::RoundingMode;
        let (value, _) = f64::rounding_from(&self.0, RoundingMode::Nearest);
        // SAFETY: positive rational rounds to positive or zero f64
        // In edge cases of very small positives, this could round to 0
        // but for physics values this won't happen
        crate::types::f64::F64::<Pos>::new(value)
            .expect("positive rational should convert to positive f64")
    }
}

impl Rational<Neg> {
    /// Convert to f64. Negative rational → negative f64.
    #[must_use]
    pub fn to_f64(&self) -> crate::types::f64::F64<Neg> {
        use malachite::num::conversion::traits::RoundingFrom;
        use malachite::rounding_modes::RoundingMode;
        let (value, _) = f64::rounding_from(&self.0, RoundingMode::Nearest);
        crate::types::f64::F64::<Neg>::new(value)
            .expect("negative rational should convert to negative f64")
    }
}

impl Rational<Finite> {
    /// Convert to f64.
    #[must_use]
    pub fn to_f64(&self) -> crate::types::f64::F64<Finite> {
        use malachite::num::conversion::traits::RoundingFrom;
        use malachite::rounding_modes::RoundingMode;
        let (value, _) = f64::rounding_from(&self.0, RoundingMode::Nearest);
        crate::types::f64::F64::<Finite>::new(value)
            .expect("rational should convert to finite f64")
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Rational as MR;

    #[test]
    fn test_pos_new() {
        assert!(Rational::<Pos>::new(MR::from(1)).is_some());
        assert!(Rational::<Pos>::new(MR::from(0)).is_none());
        assert!(Rational::<Pos>::new(MR::from(-1)).is_none());
    }

    #[test]
    fn test_neg_new() {
        assert!(Rational::<Neg>::new(MR::from(-1)).is_some());
        assert!(Rational::<Neg>::new(MR::from(0)).is_none());
        assert!(Rational::<Neg>::new(MR::from(1)).is_none());
    }

    #[test]
    fn test_to_f64() {
        let r = Rational::<Pos>::new(MR::from(42)).unwrap();
        let f = r.to_f64();
        assert_eq!(f.get(), 42.0);
    }

    #[test]
    fn test_display() {
        let r = Rational::<Pos>::new(MR::from_signeds(3, 4)).unwrap();
        assert_eq!(format!("{r}"), "3/4");
    }
}
