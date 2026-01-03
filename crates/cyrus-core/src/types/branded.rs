//! Branded types for compile-time dimension safety.
//!
//! Uses phantom lifetimes to tie moduli vectors to their originating
//! intersection tensor, making dimension mismatches impossible.

use rand::Rng;
use std::marker::PhantomData;

use super::f64::F64;
use super::tags::Pos;

/// Dimension handle branded to a specific intersection tensor.
///
/// Can only be obtained from a `NonEmptyIntersection`.
/// Used to generate or validate moduli vectors of the correct length.
pub struct Dim<'a>(pub(crate) usize, pub(crate) PhantomData<&'a ()>);

impl<'a> Dim<'a> {
    /// The number of moduli (dimension of Kähler moduli space).
    #[must_use]
    pub fn len(&self) -> usize {
        self.0
    }

    /// Generate random moduli using the provided RNG.
    ///
    /// The result is guaranteed to have the correct length and positive values.
    #[must_use]
    pub fn generate<R: Rng>(&self, rng: &mut R, min: F64<Pos>, max: F64<Pos>) -> Moduli<'a> {
        debug_assert!(max.get() > min.get(), "max must be greater than min");

        let values: Vec<F64<Pos>> = (0..self.0)
            .map(|_| {
                let val = rng.gen_range(min.get()..max.get());
                // Safe: range is positive because min and max are F64<Pos>
                F64::<Pos>::new(val).unwrap()
            })
            .collect();

        Moduli {
            values,
            _brand: PhantomData,
        }
    }

    /// Create moduli from existing positive values.
    ///
    /// Returns `None` if length doesn't match.
    #[must_use]
    pub fn from_values(&self, values: Vec<F64<Pos>>) -> Option<Moduli<'a>> {
        if values.len() != self.0 {
            return None;
        }
        Some(Moduli {
            values,
            _brand: PhantomData,
        })
    }

    /// Validate raw f64 values into moduli.
    ///
    /// Returns `None` if length doesn't match or any value is non-positive.
    #[must_use]
    pub fn validate(&self, raw: &[f64]) -> Option<Moduli<'a>> {
        if raw.len() != self.0 {
            return None;
        }
        let values: Vec<F64<Pos>> = raw
            .iter()
            .map(|&x| F64::<Pos>::new(x))
            .collect::<Option<_>>()?;
        Some(Moduli {
            values,
            _brand: PhantomData,
        })
    }
}

/// Moduli vector branded to a specific intersection tensor.
///
/// Can only be created through `Dim`, guaranteeing correct length.
/// Can only be used with the `NonEmptyIntersection` it was created from.
pub struct Moduli<'a> {
    pub(crate) values: Vec<F64<Pos>>,
    pub(crate) _brand: PhantomData<&'a ()>,
}

impl<'a> Moduli<'a> {
    /// Access the underlying values.
    #[must_use]
    pub fn values(&self) -> &[F64<Pos>] {
        &self.values
    }
}

impl<'a> std::ops::Index<usize> for Moduli<'a> {
    type Output = F64<Pos>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.values[index]
    }
}
