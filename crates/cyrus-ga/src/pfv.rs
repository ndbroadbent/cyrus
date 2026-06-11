//! Constructive perturbatively-flat-vacuum sampling.
//!
//! The PFV orthogonality condition `K . p = 0` with `p = N(M)^{-1} K` is
//! the statement that K is an isotropic vector of the rational quadratic
//! form `N^{-1}`: `K^T N^{-1} K = 0`. Random integer mutation almost never
//! lands on this measure-zero set (observed: half-million evaluations
//! across 149 polytopes plateaued at |K.p| ~ 1e-3), and greedy integer
//! repair cannot descend below the lattice step. But the condition is
//! checkable EXACTLY in rational arithmetic, so the search is inverted:
//! sample M, build the exact form once, then test candidate K vectors
//! exactly - accepting only true isotropic pairs into the population.

use malachite::Rational;
use rand::Rng;
use rand_chacha::ChaCha8Rng;

use cyrus_core::Intersection;
use cyrus_core::integer_math::invert_matrix;

use crate::genome::Genome;
use crate::geometry::GaGeometry;

/// Exact `N(M)^{-1}` for a flux vector M, or `None` when N is singular.
#[must_use]
pub fn n_inverse_for_m(kappa: &Intersection, m: &[i64]) -> Option<Vec<Vec<Rational>>> {
    let dim = m.len();
    let mut n = vec![vec![Rational::from(0); dim]; dim];
    for (a, row) in n.iter_mut().enumerate() {
        for (b, entry) in row.iter_mut().enumerate() {
            let mut sum = Rational::from(0);
            for (c, &mc) in m.iter().enumerate() {
                if mc != 0 {
                    sum += kappa.get(a, b, c).into_inner() * Rational::from(mc);
                }
            }
            *entry = sum;
        }
    }
    invert_matrix(&n)
}

/// Exact isotropy test: `K^T N^{-1} K == 0` in rational arithmetic.
#[must_use]
pub fn is_isotropic(n_inverse: &[Vec<Rational>], k: &[i64]) -> bool {
    let mut total = Rational::from(0);
    for (a, &ka) in k.iter().enumerate() {
        if ka == 0 {
            continue;
        }
        for (b, &kb) in k.iter().enumerate() {
            if kb == 0 {
                continue;
            }
            total += &n_inverse[a][b] * Rational::from(ka) * Rational::from(kb);
        }
    }
    total == 0
}

/// Sample a genome on the PFV constraint manifold.
///
/// Random M (with invertible N), then random K candidates tested against
/// the exact quadratic form. Returns `None` if the try budgets run out.
#[must_use]
pub fn sample_isotropic_genome(
    geom: &GaGeometry,
    rng: &mut ChaCha8Rng,
    flux_range: i64,
    max_m_tries: usize,
    max_k_tries: usize,
) -> Option<Genome> {
    let dim = geom.h21_primal;
    for _ in 0..max_m_tries {
        let m: Vec<i64> = (0..dim)
            .map(|_| rng.gen_range(-flux_range..=flux_range))
            .collect();
        if m.iter().all(|&x| x == 0) {
            continue;
        }
        let Some(n_inverse) = n_inverse_for_m(&geom.kappa_basis, &m) else {
            continue;
        };
        for _ in 0..max_k_tries {
            let k: Vec<i64> = (0..dim)
                .map(|_| rng.gen_range(-flux_range..=flux_range))
                .collect();
            if k.iter().all(|&x| x == 0) {
                continue;
            }
            if is_isotropic(&n_inverse, &k) {
                return Some(Genome { k, m });
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use cyrus_core::types::rational::Rational as TypedRational;
    use cyrus_core::types::tags::Finite;

    fn diag_kappa() -> Intersection {
        // kappa_abc with kappa_112 = kappa_121 = ... structure giving
        // N(M) = [[m2, m1], [m1, 0]] for M = (m1, m2) under
        // kappa_112 = 1, kappa_111 = ... keep it simple:
        // kappa_111 = 0, kappa_112 = 1, kappa_122 = 0, kappa_222 = 0
        // N_ab = kappa_abc m^c:
        //   N_11 = kappa_111 m1 + kappa_112 m2 = m2
        //   N_12 = kappa_121 m1 + kappa_122 m2 = m1
        //   N_22 = kappa_221 m1 + kappa_222 m2 = 0
        let mut kappa = Intersection::new(2);
        let one = TypedRational::<Finite>::new(Rational::from(1));
        kappa.set(0, 0, 1, one);
        kappa
    }

    #[test]
    fn isotropy_test_is_exact() {
        let kappa = diag_kappa();
        // M = (1, 0): N = [[0, 1], [1, 0]], N^{-1} = same.
        let n_inv = n_inverse_for_m(&kappa, &[1, 0]).expect("invertible");
        // Q(K) = 2 k1 k2: isotropic iff k1 == 0 or k2 == 0.
        assert!(is_isotropic(&n_inv, &[3, 0]));
        assert!(is_isotropic(&n_inv, &[0, -7]));
        assert!(!is_isotropic(&n_inv, &[1, 1]));
        assert!(!is_isotropic(&n_inv, &[5, -3]));
    }

    #[test]
    fn singular_n_returns_none() {
        let kappa = diag_kappa();
        // M = (0, 0) is skipped by samplers, but the inverse must be None
        // for singular N anyway (M = (0, 1) gives N = [[1, 0], [0, 0]]).
        assert!(n_inverse_for_m(&kappa, &[0, 1]).is_none());
    }

    #[test]
    fn sampler_finds_exact_isotropic_pairs() {
        use rand::SeedableRng;
        let kappa = diag_kappa();
        let mut rng = ChaCha8Rng::seed_from_u64(123);
        // Build a minimal geometry around the toy kappa.
        // (Only kappa_basis and h21_primal are used by the sampler.)
        let geom = GaGeometry {
            kappa_basis: kappa,
            mori: cyrus_core::MoriCone::new(vec![]),
            gv: vec![],
            h21_primal: 2,
            mirror_h11: cyrus_core::types::i32::I32::<cyrus_core::types::tags::GTEOne>::new(2)
                .expect("h11"),
            mirror_h21: cyrus_core::types::i32::I32::<cyrus_core::types::tags::NonNeg>::new(2)
                .expect("h21"),
            dual_basis: vec![0, 1],
            dual_glsm: vec![],
            q_d3: 100.0,
        };
        let genome = sample_isotropic_genome(&geom, &mut rng, 10, 64, 256).expect("finds a sample");
        let n_inv = n_inverse_for_m(&geom.kappa_basis, &genome.m).expect("invertible");
        assert!(is_isotropic(&n_inv, &genome.k));
    }
}
