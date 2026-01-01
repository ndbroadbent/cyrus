//! Flat direction computation for flux compactifications.
//!
//! Given flux vectors K and M and intersection numbers `κ_ijk`,
//! compute the flat direction p and Kähler potential factor `e^{K₀}`.
//!
//! The flat direction satisfies:
//! ```text
//! N_ab p^b = K_a  where  N_ab = κ_abc M^c
//! ```
//!
//! The Kähler potential factor is:
//! ```text
//! e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}
//! ```
//!
//! Reference: arXiv:2107.09064, Section 6

use crate::intersection::Intersection;

/// Compute the N matrix: `N_ab = κ_abc M^c`.
///
/// Contracts the intersection tensor with flux vector M.
///
/// # Panics
/// Panics if `m.len() != kappa.dim()`.
#[must_use]
#[allow(clippy::cast_precision_loss)] // Physics values fit in f64 mantissa
pub fn compute_n_matrix(kappa: &Intersection, m: &[i64]) -> Vec<Vec<f64>> {
    let dim = kappa.dim();
    assert_eq!(m.len(), dim, "dimension mismatch");

    let mut n = vec![vec![0.0; dim]; dim];

    for ((i, j, k), val) in kappa.iter() {
        let val_f = val as f64;
        // Contract with M^c
        // N_ab = κ_abc M^c
        // Since κ is symmetric, we need to account for all permutations
        // κ_{ijk} contributes to N_{ij} when we contract with M_k

        // For canonical entry (i,j,k) with i <= j <= k:
        // - If all different: contributes to N_{ij}, N_{ik}, N_{jk} etc.
        // - Handle each unique pair

        let m_i = m[i] as f64;
        let m_j = m[j] as f64;
        let m_k = m[k] as f64;

        if i == j && j == k {
            // All same: (i,i,i) - only N_{ii} += κ_{iii} M^i
            n[i][i] += val_f * m_i;
        } else if i == j {
            // Two same: (i,i,k) with i < k
            // N_{ii} += κ_{iik} M^k
            // N_{ik} += κ_{iik} M^i × 2 (from iik and iki)
            n[i][i] += val_f * m_k;
            n[i][k] += 2.0 * val_f * m_i;
            n[k][i] += 2.0 * val_f * m_i;
        } else if j == k {
            // Two same: (i,j,j) with i < j
            // N_{jj} += κ_{ijj} M^i
            // N_{ij} += κ_{ijj} M^j × 2
            n[j][j] += val_f * m_i;
            n[i][j] += 2.0 * val_f * m_j;
            n[j][i] += 2.0 * val_f * m_j;
        } else if i == k {
            // This shouldn't happen with canonical ordering (i <= j <= k)
            unreachable!("canonical ordering violated: i == k but i != j");
        } else {
            // All different: (i,j,k) with i < j < k
            // Each pair gets contributions from contracting with the third index
            // N_{ij} += κ_{ijk} M^k × 2 (from ijk and jik)
            // N_{ik} += κ_{ijk} M^j × 2 (from ijk and kij)
            // N_{jk} += κ_{ijk} M^i × 2 (from jki and kjj)
            n[i][j] += 2.0 * val_f * m_k;
            n[j][i] += 2.0 * val_f * m_k;
            n[i][k] += 2.0 * val_f * m_j;
            n[k][i] += 2.0 * val_f * m_j;
            n[j][k] += 2.0 * val_f * m_i;
            n[k][j] += 2.0 * val_f * m_i;
        }
    }

    n
}

/// Solve `N p = K` for p using Gaussian elimination with partial pivoting.
///
/// Returns `None` if the matrix is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
#[allow(clippy::cast_precision_loss)] // Physics values fit in f64 mantissa
pub fn solve_linear_system(n: &[Vec<f64>], k: &[i64]) -> Option<Vec<f64>> {
    let dim = n.len();
    assert!(n.iter().all(|row| row.len() == dim), "N must be square");
    assert_eq!(k.len(), dim, "K dimension mismatch");

    // Create augmented matrix
    let mut aug: Vec<Vec<f64>> = n
        .iter()
        .zip(k.iter())
        .map(|(row, &ki)| {
            let mut r = row.clone();
            r.push(ki as f64);
            r
        })
        .collect();

    // Gaussian elimination with partial pivoting
    for col in 0..dim {
        // Find pivot
        let (max_row, max_val) = aug
            .iter()
            .enumerate()
            .skip(col)
            .map(|(row, r)| (row, r[col].abs()))
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();

        if max_val < 1e-14 {
            return None; // Singular matrix
        }

        // Swap rows
        aug.swap(col, max_row);

        // Eliminate column
        let pivot = aug[col][col];
        for row in (col + 1)..dim {
            let factor = aug[row][col] / pivot;
            for c in col..=dim {
                aug[row][c] -= factor * aug[col][c];
            }
        }
    }

    // Back substitution
    let mut p = vec![0.0; dim];
    for i in (0..dim).rev() {
        let mut sum = aug[i][dim];
        for (j, &pj) in p.iter().enumerate().skip(i + 1) {
            sum -= aug[i][j] * pj;
        }
        p[i] = sum / aug[i][i];
    }

    Some(p)
}

/// Compute flat direction p = N⁻¹ K where `N_ab = κ_abc M^c`.
///
/// Returns `None` if N is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
pub fn compute_flat_direction(kappa: &Intersection, k: &[i64], m: &[i64]) -> Option<Vec<f64>> {
    let n = compute_n_matrix(kappa, m);
    solve_linear_system(&n, k)
}

/// Compute `e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}`.
///
/// # Panics
/// Panics if `p.len() != kappa.dim()`.
#[must_use]
pub fn compute_ek0(kappa: &Intersection, p: &[f64]) -> f64 {
    let kappa_ppp = kappa.contract_triple(p);
    1.0 / ((4.0 / 3.0) * kappa_ppp)
}

/// Result of flat direction computation.
#[derive(Debug, Clone)]
pub struct FlatDirectionResult {
    /// N matrix: `N_ab = κ_abc M^c`.
    pub n_matrix: Vec<Vec<f64>>,
    /// Flat direction: p = N⁻¹ K.
    pub p: Vec<f64>,
    /// Kähler potential factor: `e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}`.
    pub ek0: f64,
}

/// Compute flat direction with full breakdown.
///
/// Returns `None` if N is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
pub fn compute_flat_direction_full(
    kappa: &Intersection,
    k: &[i64],
    m: &[i64],
) -> Option<FlatDirectionResult> {
    let n_matrix = compute_n_matrix(kappa, m);
    let p = solve_linear_system(&n_matrix, k)?;
    let ek0 = compute_ek0(kappa, &p);

    Some(FlatDirectionResult { n_matrix, p, ek0 })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_n_matrix_simple() {
        // κ_000 = 6 means N_00 = 6 * M^0
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let m = vec![2];
        let n = compute_n_matrix(&kappa, &m);

        assert_eq!(n.len(), 1);
        assert!((n[0][0] - 12.0).abs() < 1e-10); // 6 * 2 = 12
    }

    #[test]
    fn test_solve_linear_system() {
        // Simple 2x2: [2, 1; 1, 3] x = [5, 5]
        // Solution: x = [2, 1]
        let n = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let k = vec![5, 5];

        let p = solve_linear_system(&n, &k).unwrap();

        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!((p[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_flat_direction() {
        // With κ_000 = 6, M = [1], K = [3]:
        // N_00 = 6, so p = 6⁻¹ × 3 = 0.5
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let k = vec![3];
        let m = vec![1];

        let p = compute_flat_direction(&kappa, &k, &m).unwrap();
        assert!((p[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_compute_ek0() {
        // κ_000 = 6, p = [1]
        // κ_ppp = 6 × 1 = 6 (note: multiplicity 1 for iii)
        // e^K₀ = 1 / (4/3 × 6) = 1/8 = 0.125
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let p = vec![1.0];
        let ek0 = compute_ek0(&kappa, &p);

        assert!((ek0 - 0.125).abs() < 1e-10);
    }
}
