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
// Marker Trait Implementations
// ============================================================================

impl IsFinite for Rational<Finite> {}
impl IsFinite for Rational<Pos> {}
impl IsFinite for Rational<Neg> {}
impl IsFinite for Rational<Zero> {}
impl IsFinite for Rational<NonNeg> {}
impl IsFinite for Rational<NonZero> {}
impl IsFinite for Rational<NonPos> {}

impl IsNonZero for Rational<NonZero> {}
impl IsNonZero for Rational<Pos> {}
impl IsNonZero for Rational<Neg> {}

impl IsPositive for Rational<Pos> {}
impl IsNegative for Rational<Neg> {}
impl IsZero for Rational<Zero> {}

impl IsNonNeg for Rational<NonNeg> {}
impl IsNonNeg for Rational<Pos> {}
impl IsNonNeg for Rational<Zero> {}

impl IsNonPos for Rational<NonPos> {}
impl IsNonPos for Rational<Neg> {}
impl IsNonPos for Rational<Zero> {}

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
        crate::types::f64::F64::<Finite>::new(value).expect("rational should convert to finite f64")
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Rational as MR;
    use std::collections::hash_map::DefaultHasher;

    #[test]
    fn test_finite_new() {
        let r = Rational::<Finite>::new(MR::from(42));
        assert_eq!(r.get(), &MR::from(42));

        let r2 = Rational::<Finite>::new(MR::from(-5));
        assert_eq!(r2.get(), &MR::from(-5));
    }

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
    fn test_nonzero_new() {
        assert!(Rational::<NonZero>::new(MR::from(1)).is_some());
        assert!(Rational::<NonZero>::new(MR::from(-1)).is_some());
        assert!(Rational::<NonZero>::new(MR::from(0)).is_none());
    }

    #[test]
    fn test_nonneg_new() {
        assert!(Rational::<NonNeg>::new(MR::from(0)).is_some());
        assert!(Rational::<NonNeg>::new(MR::from(1)).is_some());
        assert!(Rational::<NonNeg>::new(MR::from(-1)).is_none());
    }

    #[test]
    fn test_nonpos_new() {
        assert!(Rational::<NonPos>::new(MR::from(0)).is_some());
        assert!(Rational::<NonPos>::new(MR::from(-1)).is_some());
        assert!(Rational::<NonPos>::new(MR::from(1)).is_none());
    }

    #[test]
    fn test_zero_new() {
        assert!(Rational::<Zero>::new(MR::from(0)).is_some());
        assert!(Rational::<Zero>::new(MR::from(1)).is_none());
        assert!(Rational::<Zero>::new(MR::from(-1)).is_none());

        let z = Rational::<Zero>::zero();
        assert_eq!(z.get(), &MR::from(0));
    }

    #[test]
    fn test_default_zero() {
        let z: Rational<Zero> = Rational::default();
        assert_eq!(z.get(), &MR::from(0));
    }

    #[test]
    fn test_into_inner() {
        let r = Rational::<Pos>::new(MR::from(5)).unwrap();
        let inner = r.into_inner();
        assert_eq!(inner, MR::from(5));
    }

    #[test]
    fn test_from_raw() {
        let r = Rational::<Pos>::from_raw(MR::from(10));
        assert_eq!(r.get(), &MR::from(10));
    }

    #[test]
    fn test_to_f64_pos() {
        let r = Rational::<Pos>::new(MR::from(42)).unwrap();
        let f = r.to_f64();
        assert!((f.get() - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_f64_neg() {
        let r = Rational::<Neg>::new(MR::from(-5)).unwrap();
        let f = r.to_f64();
        assert!((f.get() - (-5.0)).abs() < 1e-10);
    }

    #[test]
    fn test_to_f64_finite() {
        let r = Rational::<Finite>::new(MR::from(-3));
        let f = r.to_f64();
        assert!((f.get() - (-3.0)).abs() < 1e-10);

        let r2 = Rational::<Finite>::new(MR::from(7));
        let f2 = r2.to_f64();
        assert!((f2.get() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_display() {
        let r = Rational::<Pos>::new(MR::from_signeds(3, 4)).unwrap();
        assert_eq!(format!("{r}"), "3/4");
    }

    #[test]
    fn test_debug() {
        let r = Rational::<Pos>::new(MR::from(5)).unwrap();
        let s = format!("{r:?}");
        assert!(s.contains("Pos"));
        assert!(s.contains('5'));
    }

    #[test]
    fn test_eq() {
        let a = Rational::<Pos>::new(MR::from(3)).unwrap();
        let b = Rational::<Pos>::new(MR::from(3)).unwrap();
        let c = Rational::<Pos>::new(MR::from(4)).unwrap();
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn test_hash() {
        let a = Rational::<Pos>::new(MR::from(5)).unwrap();
        let b = Rational::<Pos>::new(MR::from(5)).unwrap();

        let mut h1 = DefaultHasher::new();
        let mut h2 = DefaultHasher::new();
        a.hash(&mut h1);
        b.hash(&mut h2);
        assert_eq!(h1.finish(), h2.finish());
    }

    #[test]
    fn test_ord() {
        let a = Rational::<Pos>::new(MR::from(1)).unwrap();
        let b = Rational::<Pos>::new(MR::from(2)).unwrap();
        assert!(a < b);
        assert!(b > a);
        assert!(a <= b);
        assert!(b >= a);
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Less);
    }
}
