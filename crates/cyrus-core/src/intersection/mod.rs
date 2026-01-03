//! Intersection number computation.
//!
//! Implements the clean-room algorithm for computing intersection numbers
//! of divisors in a Calabi-Yau threefold.
//!
//! Algorithm:
//! 1. Identify distinct intersections from the triangulation (facets).
//! 2. Use GLSM linear relations to constrain the remaining intersections.
//! 3. Solve the linear system to find all entries of the intersection tensor.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 1.2

mod compute;

use crate::f64_pos;
use crate::types::f64::F64;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::{Finite, IsFinite, Pos};
use malachite::Rational;
use std::collections::HashMap;

pub use compute::compute_intersection_numbers;

/// Intersection tensor `κ_ijk`.
///
/// Stores intersection numbers for divisors `D_i, D_j, D_k`.
/// Symmetric in all indices. Values can be positive or negative.
#[derive(Debug, Clone)]
pub struct Intersection {
    /// Dimension of the divisor space (number of divisors).
    dim: usize,
    /// Storage for unique entries (i ≤ j ≤ k).
    /// Key: (i, j, k), Value: Intersection number (any finite value).
    entries: HashMap<(usize, usize, usize), TypedRational<Finite>>,
}

impl Intersection {
    /// Create a new empty intersection tensor.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            entries: HashMap::new(),
        }
    }

    /// Get the dimension (number of divisors).
    pub const fn dim(&self) -> usize {
        self.dim
    }

    /// Get the number of non-zero entries.
    pub fn num_nonzero(&self) -> usize {
        self.entries.len()
    }

    /// Get value at indices (i, j, k). Returns zero if not set.
    pub fn get(&self, i: usize, j: usize, k: usize) -> TypedRational<Finite> {
        let key = canonical_key(i, j, k);
        self.entries
            .get(&key)
            .cloned()
            .unwrap_or_else(|| TypedRational::<Finite>::from_raw(Rational::from(0)))
    }

    /// Set value at indices (i, j, k).
    ///
    /// Accepts any rational type (Pos, Neg, Finite, etc.).
    pub fn set<T>(&mut self, i: usize, j: usize, k: usize, val: TypedRational<T>)
    where
        TypedRational<T>: IsFinite,
    {
        let key = canonical_key(i, j, k);
        // All IsFinite types have valid finite values - just re-wrap
        let val_finite = TypedRational::<Finite>::from_raw(val.into_inner());
        if *val_finite.get() == 0 {
            self.entries.remove(&key);
        } else {
            self.entries.insert(key, val_finite);
        }
    }

    /// Iterator over all non-zero entries.
    /// Returns ((i, j, k), val) where i ≤ j ≤ k.
    pub fn iter(&self) -> impl Iterator<Item = (&(usize, usize, usize), &TypedRational<Finite>)> {
        self.entries.iter()
    }

    /// Parallel iterator over non-zero entries.
    pub fn par_iter(
        &self,
    ) -> impl rayon::iter::ParallelIterator<Item = (&(usize, usize, usize), &TypedRational<Finite>)>
    where
        HashMap<(usize, usize, usize), TypedRational<Finite>>:
            rayon::iter::IntoParallelRefIterator<'static>,
    {
        use rayon::prelude::*;
        self.entries.par_iter()
    }

    /// Non-parallel iterator over non-zero entries.
    pub fn iter_entries(
        &self,
    ) -> std::collections::hash_map::Iter<'_, (usize, usize, usize), TypedRational<Finite>> {
        self.entries.iter()
    }

    /// Contract with positive moduli `t` to compute `κ_ijk t^i t^j t^k`.
    ///
    /// Input moduli must be positive (F64<Pos>) - for Kähler moduli.
    /// Returns `None` if the contraction is not positive (invalid physics).
    ///
    /// # Panics
    /// Panics if `t.len() != dim`.
    pub fn contract_triple(&self, t: &[F64<Pos>]) -> Option<F64<Pos>> {
        assert_eq!(t.len(), self.dim, "Vector dimension mismatch");

        // Accumulate as Finite (sum of terms could be any sign)
        // Type algebra: Pos * Finite = Finite, Finite + Finite = Finite
        let mut sum = F64::<Finite>::ZERO;
        for (&(i, j, k), val) in &self.entries {
            let mult = symmetry_multiplicity(i, j, k);
            let kappa = val.to_f64();
            // mult: Pos, kappa: Finite, t[i]: Pos → Pos * Finite * Pos * Pos * Pos = Finite
            let term = mult * kappa * t[i] * t[j] * t[k];
            sum = sum + term;
        }

        sum.try_to_pos()
    }

    /// Contract with finite-valued vector `t` to compute `κ_ijk t^i t^j t^k`.
    ///
    /// For cases where the vector can have any sign (e.g., flat directions).
    /// Returns `None` if the contraction is not positive.
    ///
    /// # Panics
    /// Panics if `t.len() != dim`.
    pub fn contract_triple_finite(&self, t: &[F64<Finite>]) -> Option<F64<Pos>> {
        assert_eq!(t.len(), self.dim, "Vector dimension mismatch");

        // Type algebra: Pos * Finite = Finite, Finite * Finite = Finite
        let mut sum = F64::<Finite>::ZERO;
        for (&(i, j, k), val) in &self.entries {
            let mult = symmetry_multiplicity(i, j, k);
            let kappa = val.to_f64();
            // mult: Pos, kappa: Finite, t[i]: Finite → all Finite
            let term = mult * kappa * t[i] * t[j] * t[k];
            sum = sum + term;
        }

        sum.try_to_pos()
    }

    /// Convert to NonEmptyIntersection if this tensor has entries.
    ///
    /// Returns `None` if the tensor is empty (no non-zero entries).
    #[must_use]
    pub fn into_non_empty(self) -> Option<NonEmptyIntersection> {
        if self.entries.is_empty() {
            None
        } else {
            Some(NonEmptyIntersection(self))
        }
    }
}

/// Non-empty intersection tensor.
///
/// Wrapper around `Intersection` that guarantees at least one non-zero entry.
/// Provides branded dimension for type-safe moduli handling.
#[derive(Debug, Clone)]
pub struct NonEmptyIntersection(Intersection);

impl NonEmptyIntersection {
    /// Get the dimension handle branded to this tensor.
    ///
    /// Use this to create `Moduli` vectors of the correct length.
    #[must_use]
    pub const fn dim(&self) -> crate::types::branded::Dim<'_> {
        crate::types::branded::Dim(self.0.dim, std::marker::PhantomData)
    }

    /// Contract with validated moduli to compute `κ_ijk t^i t^j t^k`.
    ///
    /// # Panics
    /// Panics if the contraction is non-positive (moduli outside Kähler cone).
    #[must_use]
    pub fn contract_triple(&self, t: &crate::types::branded::Moduli<'_>) -> F64<Pos> {
        // Type algebra: accumulate as Finite, narrow at boundary
        let mut sum = F64::<Finite>::ZERO;
        for (&(i, j, k), val) in &self.0.entries {
            let mult = symmetry_multiplicity(i, j, k);
            let kappa = val.to_f64();
            let term = mult * kappa * t[i] * t[j] * t[k];
            sum = sum + term;
        }

        sum.try_to_pos()
            .expect("contraction must be positive - moduli outside Kähler cone")
    }

    /// Get the underlying intersection tensor.
    #[must_use]
    pub const fn inner(&self) -> &Intersection {
        &self.0
    }
}

/// Compute canonical key (i ≤ j ≤ k) for symmetric tensor.
pub(crate) fn canonical_key(i: usize, j: usize, k: usize) -> (usize, usize, usize) {
    let mut v = [i, j, k];
    v.sort_unstable();
    v.into()
}

/// Compute symmetry multiplicity for entry (i, j, k) with i ≤ j ≤ k.
/// Returns F64<Pos> (1, 3, or 6 are all positive).
const fn symmetry_multiplicity(i: usize, j: usize, k: usize) -> F64<Pos> {
    if i == j && j == k {
        f64_pos!(1.0)
    } else if i == j || j == k || i == k {
        f64_pos!(3.0)
    } else {
        f64_pos!(6.0)
    }
}

#[cfg(test)]
mod tests;
