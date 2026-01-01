//! Intersection numbers for Calabi-Yau manifolds.
//!
//! The intersection tensor `κ_ijk` encodes the triple intersection numbers
//! of divisors on a Calabi-Yau threefold. It is fully symmetric:
//! `κ_ijk = κ_jik = κ_kij` etc.
//!
//! Reference: arXiv:2107.09064

use std::collections::HashMap;

/// Sparse representation of intersection numbers.
///
/// Only stores unique entries (i ≤ j ≤ k) since the tensor is symmetric.
#[derive(Debug, Clone, Default)]
pub struct Intersection {
    /// Dimension (h11).
    dim: usize,
    /// Map from (i, j, k) with i ≤ j ≤ k to intersection number.
    entries: HashMap<(usize, usize, usize), i64>,
}

impl Intersection {
    /// Create a new empty intersection tensor of given dimension.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            entries: HashMap::new(),
        }
    }

    /// Get the dimension (h11).
    pub const fn dim(&self) -> usize {
        self.dim
    }

    /// Set an intersection number `κ_ijk`.
    ///
    /// Automatically handles symmetry - only stores canonical form (i ≤ j ≤ k).
    pub fn set(&mut self, i: usize, j: usize, k: usize, value: i64) {
        let key = canonical_key(i, j, k);
        if value == 0 {
            self.entries.remove(&key);
        } else {
            self.entries.insert(key, value);
        }
    }

    /// Get an intersection number `κ_ijk`.
    pub fn get(&self, i: usize, j: usize, k: usize) -> i64 {
        let key = canonical_key(i, j, k);
        self.entries.get(&key).copied().unwrap_or(0)
    }

    /// Number of non-zero entries (in canonical form).
    pub fn num_nonzero(&self) -> usize {
        self.entries.len()
    }

    /// Iterate over non-zero entries as ((i, j, k), value) with i ≤ j ≤ k.
    pub fn iter(&self) -> impl Iterator<Item = ((usize, usize, usize), i64)> + '_ {
        self.entries.iter().map(|(&k, &v)| (k, v))
    }

    /// Compute `κ_ijk t^i t^j t^k` (triple contraction).
    ///
    /// This is used in volume computations: V = (1/6) `κ_ijk` t^i t^j t^k
    ///
    /// # Panics
    /// Panics if `t.len() != self.dim()`.
    #[allow(clippy::cast_precision_loss)] // mult is at most 6, val fits in f64 mantissa for physics
    pub fn contract_triple(&self, t: &[f64]) -> f64 {
        assert_eq!(t.len(), self.dim, "dimension mismatch");

        let mut result = 0.0;
        for (&(i, j, k), &val) in &self.entries {
            // Count multiplicity based on symmetry
            let mult = symmetry_multiplicity(i, j, k);
            result += f64::from(mult) * (val as f64) * t[i] * t[j] * t[k];
        }
        result
    }
}

/// Compute canonical key (i ≤ j ≤ k) for symmetric tensor.
fn canonical_key(i: usize, j: usize, k: usize) -> (usize, usize, usize) {
    let mut v = [i, j, k];
    v.sort_unstable();
    v.into()
}

/// Compute symmetry multiplicity for entry (i, j, k) with i ≤ j ≤ k.
///
/// - If all equal: multiplicity 1 (iii)
/// - If two equal: multiplicity 3 (iij, iji, jii or ijj, jij, jji)
/// - If all different: multiplicity 6 (all permutations)
const fn symmetry_multiplicity(i: usize, j: usize, k: usize) -> u8 {
    if i == j && j == k {
        1
    } else if i == j || j == k || i == k {
        3
    } else {
        6
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_get() {
        let mut kappa = Intersection::new(3);
        kappa.set(0, 1, 2, 5);

        // All permutations should return the same value
        assert_eq!(kappa.get(0, 1, 2), 5);
        assert_eq!(kappa.get(0, 2, 1), 5);
        assert_eq!(kappa.get(1, 0, 2), 5);
        assert_eq!(kappa.get(1, 2, 0), 5);
        assert_eq!(kappa.get(2, 0, 1), 5);
        assert_eq!(kappa.get(2, 1, 0), 5);
    }

    #[test]
    fn test_symmetry_multiplicity() {
        assert_eq!(symmetry_multiplicity(0, 0, 0), 1);
        assert_eq!(symmetry_multiplicity(0, 0, 1), 3);
        assert_eq!(symmetry_multiplicity(0, 1, 1), 3);
        assert_eq!(symmetry_multiplicity(0, 1, 2), 6);
    }

    #[test]
    fn test_contract_triple() {
        // Simple case: κ_000 = 6 means V = (1/6) * 6 * t^3 = t^3
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let t = vec![2.0];
        let result = kappa.contract_triple(&t);
        // κ_000 t^0 t^0 t^0 = 6 * 8 = 48
        assert!((result - 48.0).abs() < 1e-10);
    }

    #[test]
    fn test_contract_mixed() {
        // κ_012 = 1, all indices different
        let mut kappa = Intersection::new(3);
        kappa.set(0, 1, 2, 1);

        let t = vec![1.0, 2.0, 3.0];
        let result = kappa.contract_triple(&t);
        // 6 permutations, each contributes 1 * 1 * 2 * 3 = 6
        assert!((result - 36.0).abs() < 1e-10);
    }
}
