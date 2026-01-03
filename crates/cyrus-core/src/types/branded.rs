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
    pub const fn len(&self) -> usize {
        self.0
    }

    /// Returns true if the dimension is zero.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.0 == 0
    }

    /// Generate random moduli using the provided RNG.
    ///
    /// The result is guaranteed to have the correct length and positive values.
    ///
    /// # Panics
    /// Panics in debug mode if `max <= min`.
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

impl Moduli<'_> {
    /// Access the underlying values.
    #[must_use]
    pub fn values(&self) -> &[F64<Pos>] {
        &self.values
    }
}

impl std::ops::Index<usize> for Moduli<'_> {
    type Output = F64<Pos>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.values[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::f64_pos;
    use rand::SeedableRng;

    // Helper to create a Dim for testing (simulate getting it from an intersection)
    fn test_dim<'a>(n: usize) -> Dim<'a> {
        Dim(n, PhantomData)
    }

    #[test]
    fn test_dim_len() {
        let dim = test_dim(5);
        assert_eq!(dim.len(), 5);

        let dim2 = test_dim(10);
        assert_eq!(dim2.len(), 10);
    }

    #[test]
    fn test_generate() {
        let dim = test_dim(3);
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let min = f64_pos!(0.1);
        let max = f64_pos!(10.0);

        let moduli = dim.generate(&mut rng, min, max);
        assert_eq!(moduli.values().len(), 3);

        // Check all values are within range
        for v in moduli.values() {
            assert!(v.get() >= 0.1);
            assert!(v.get() <= 10.0);
        }
    }

    #[test]
    fn test_from_values_success() {
        let dim = test_dim(2);
        let values = vec![f64_pos!(1.0), f64_pos!(2.0)];

        let moduli = dim.from_values(values).unwrap();
        assert_eq!(moduli.values().len(), 2);
        assert!((moduli[0].get() - 1.0).abs() < 1e-10);
        assert!((moduli[1].get() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_from_values_wrong_length() {
        let dim = test_dim(3);
        let values = vec![f64_pos!(1.0), f64_pos!(2.0)]; // Only 2 values

        assert!(dim.from_values(values).is_none());
    }

    #[test]
    fn test_validate_success() {
        let dim = test_dim(2);
        let raw = [1.5, 2.5];

        let moduli = dim.validate(&raw).unwrap();
        assert_eq!(moduli.values().len(), 2);
        assert!((moduli[0].get() - 1.5).abs() < 1e-10);
        assert!((moduli[1].get() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_validate_wrong_length() {
        let dim = test_dim(3);
        let raw = [1.0, 2.0]; // Only 2 values

        assert!(dim.validate(&raw).is_none());
    }

    #[test]
    fn test_validate_non_positive() {
        let dim = test_dim(2);
        let raw = [1.0, -0.5]; // Negative value

        assert!(dim.validate(&raw).is_none());
    }

    #[test]
    fn test_validate_zero() {
        let dim = test_dim(2);
        let raw = [1.0, 0.0]; // Zero value

        assert!(dim.validate(&raw).is_none());
    }

    #[test]
    fn test_moduli_index() {
        let dim = test_dim(3);
        let values = vec![f64_pos!(1.0), f64_pos!(2.0), f64_pos!(3.0)];
        let moduli = dim.from_values(values).unwrap();

        assert!((moduli[0].get() - 1.0).abs() < 1e-10);
        assert!((moduli[1].get() - 2.0).abs() < 1e-10);
        assert!((moduli[2].get() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_moduli_values() {
        let dim = test_dim(2);
        let values = vec![f64_pos!(5.0), f64_pos!(10.0)];
        let moduli = dim.from_values(values).unwrap();

        let vals = moduli.values();
        assert_eq!(vals.len(), 2);
        assert!((vals[0].get() - 5.0).abs() < 1e-10);
        assert!((vals[1].get() - 10.0).abs() < 1e-10);
    }
}
