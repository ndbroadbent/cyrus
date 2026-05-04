//! Matrix utility functions.
//!
//! Provides common matrix operations like rank computation, inversion,
//! and GCD/LCM for integers.

use malachite::Integer;
use malachite::Rational;
use malachite::num::arithmetic::traits::Abs;

/// Compute the rank of a matrix.
pub fn matrix_rank(mat: &[Vec<Rational>]) -> usize {
    if mat.is_empty() {
        return 0;
    }

    let n_rows = mat.len();
    let n_cols = mat[0].len();
    let mut work = mat.to_vec();

    let mut rank = 0;
    let mut col = 0;

    while rank < n_rows && col < n_cols {
        // Find pivot
        let mut pivot_row = None;
        for r in rank..n_rows {
            if work[r][col] != 0 {
                pivot_row = Some(r);
                break;
            }
        }

        if let Some(pr) = pivot_row {
            work.swap(rank, pr);
            let pivot_val = work[rank][col].clone();
            for c in col..n_cols {
                work[rank][c] /= &pivot_val;
            }
            for r in 0..n_rows {
                if r != rank && work[r][col] != 0 {
                    let factor = work[r][col].clone();
                    for c in col..n_cols {
                        let sub = &factor * &work[rank][c];
                        work[r][c] -= sub;
                    }
                }
            }
            rank += 1;
        }
        col += 1;
    }

    rank
}

/// Invert a square matrix of rationals.
pub fn invert_matrix(mat: &[Vec<Rational>]) -> Option<Vec<Vec<Rational>>> {
    let n = mat.len();
    if n == 0 || mat[0].len() != n {
        return None;
    }

    // Augmented matrix [A | I]
    let mut aug: Vec<Vec<Rational>> = vec![vec![Rational::from(0); 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = mat[i][j].clone();
        }
        aug[i][n + i] = Rational::from(1);
    }

    // Gaussian elimination
    for i in 0..n {
        // Find pivot
        let mut pivot = i;
        while pivot < n && aug[pivot][i] == 0 {
            pivot += 1;
        }
        if pivot == n {
            return None; // Singular
        }

        aug.swap(i, pivot);

        let pivot_val = aug[i][i].clone();
        for j in 0..(2 * n) {
            aug[i][j] /= &pivot_val;
        }

        for r in 0..n {
            if r != i && aug[r][i] != 0 {
                let factor = aug[r][i].clone();
                for j in 0..(2 * n) {
                    let sub = &factor * &aug[i][j];
                    aug[r][j] -= sub;
                }
            }
        }
    }

    // Extract inverse
    let inv: Vec<Vec<Rational>> = aug.into_iter().map(|row| row[n..].to_vec()).collect();

    Some(inv)
}

/// Compute LCM of two integers.
pub fn lcm_integer(a: &Integer, b: &Integer) -> Integer {
    if *a == 0 || *b == 0 {
        return Integer::from(0);
    }
    let gcd = gcd_integer(a, b);
    (a / &gcd) * b
}

/// Compute GCD of two integers.
pub fn gcd_integer(a: &Integer, b: &Integer) -> Integer {
    if *b == 0 {
        a.clone().abs()
    } else {
        gcd_integer(b, &(a % b))
    }
}

/// Convert a rational matrix to integers by scaling each row.
pub fn rational_matrix_to_integer(mat: &[Vec<Rational>]) -> Vec<Vec<Integer>> {
    let mut result = Vec::with_capacity(mat.len());

    for row in mat {
        // Find LCM of denominators
        let mut lcm = Integer::from(1);
        for val in row {
            let (_, denom) = val.clone().into_numerator_and_denominator();
            lcm = lcm_integer(&lcm, &Integer::from(denom));
        }

        // Scale row by LCM
        let int_row: Vec<Integer> = row
            .iter()
            .map(|val| {
                let scaled = val * Rational::from(&lcm);
                Integer::try_from(scaled).expect("Non-integer after scaling")
            })
            .collect();

        result.push(int_row);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_matrix_to_integer_preserves_negative_entries() {
        let mat = vec![vec![
            Rational::from(-1) / Rational::from(2),
            Rational::from(3) / Rational::from(2),
        ]];

        assert_eq!(
            rational_matrix_to_integer(&mat),
            vec![vec![Integer::from(-1), Integer::from(3)]]
        );
    }
}
