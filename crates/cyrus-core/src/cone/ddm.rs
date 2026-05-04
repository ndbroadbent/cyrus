//! Double Description Method (DDM) for cone dualization.
//!
//! Converts between H-representation (hyperplanes) and V-representation (rays).
//!
//! Reference: PPL (Parma Polyhedra Library), cddlib

use std::collections::HashMap;

mod init;
mod order;
mod rank;
mod step;
mod types;
mod util;

use init::{full_space_initialization, initialize_ddm};
use order::choose_next_hyperplane;
use step::{deduplicate_rays, log_step, process_hyperplane, should_log_step};

#[cfg(test)]
mod tests;

#[cfg(test)]
use crate::integer_math::invert_matrix;
#[cfg(test)]
use malachite::{Integer, Rational};
#[cfg(test)]
use order::{choose_next_hyperplane_with_window, order_remaining_hyperplanes};
#[cfg(test)]
use rank::{
    RANK_MODULUS, RankContext, RankTracker, active_rank_at_least, active_rank_mod_prime_at_least,
    are_adjacent, integer_mod_prime, quotient_rank_mod_prime_at_least,
};
#[cfg(test)]
use step::{compute_intersection_ray, parse_pair_log_interval};
#[cfg(test)]
use types::{DdmHyperplane, DdmRay, RankStats};
#[cfg(test)]
use util::{
    active_set_for_intersection, dot_product, hyperplanes_from_dense,
    normalize_ray_preserving_direction,
};

/// Convert hyperplanes to rays (or vice versa).
///
/// Given H = hyperplanes defining cone {x : H @ x >= 0},
/// compute rays R such that cone = {R^T @ λ : λ >= 0}.
///
/// The algorithm is symmetric: dualize(rays) gives hyperplanes,
/// dualize(hyperplanes) gives rays.
pub fn dualize(matrix: &[Vec<i128>], ambient_dim: usize) -> Vec<Vec<i128>> {
    if matrix.is_empty() {
        // No constraints = full space
        // Return ±e_i for each dimension
        let mut rays = Vec::new();
        for i in 0..ambient_dim {
            let mut pos = vec![0i128; ambient_dim];
            let mut neg = vec![0i128; ambient_dim];
            pos[i] = 1;
            neg[i] = -1;
            rays.push(pos);
            rays.push(neg);
        }
        return rays;
    }

    if ambient_dim == 0 {
        return Vec::new();
    }

    // Use the incremental DDM algorithm
    ddm_incremental(matrix, ambient_dim)
}

/// Incremental Double Description Method.
///
/// Start with the full space (identity rays), then add constraints one by one.
fn ddm_incremental(hyperplanes: &[Vec<i128>], dim: usize) -> Vec<Vec<i128>> {
    let init = initialize_ddm(hyperplanes, dim).unwrap_or_else(|| {
        eprintln!("[WARN] ddm: full-rank initialization failed; starting from full space");
        full_space_initialization(hyperplanes, dim)
    });
    let (mut rays, mut ordered_hyperplanes, start_idx, mut processed_rank, mut rank_context) = init;
    let mut adjacency_cache = HashMap::<Vec<usize>, bool>::new();

    // Process each hyperplane
    for idx in start_idx..ordered_hyperplanes.len() {
        let best_idx = choose_next_hyperplane(&ordered_hyperplanes, idx, &rays);
        if best_idx != idx {
            ordered_hyperplanes.swap(idx, best_idx);
        }
        let h = &ordered_hyperplanes[idx];
        let pointed_before = processed_rank.rank() == dim;
        let (next_rays, stats) = process_hyperplane(
            &rays,
            h,
            idx,
            &ordered_hyperplanes[..idx],
            dim,
            pointed_before,
            &mut rank_context,
            &mut adjacency_cache,
        );
        if stats.redundant {
            if should_log_step(idx, &stats) {
                log_step(idx, ordered_hyperplanes.len(), rays.len(), &stats);
            }
            continue;
        }
        rays = next_rays;
        if rays.len() > 2000 || idx % 20 == 0 {
            rays = deduplicate_rays(rays);
        }
        processed_rank.add_row(h.dense());
        if should_log_step(idx, &stats) {
            log_step(idx, ordered_hyperplanes.len(), rays.len(), &stats);
        }

        // Early termination if cone becomes empty
        if rays.is_empty() {
            break;
        }
    }

    // Remove duplicates and normalize
    deduplicate_rays(rays)
        .into_iter()
        .map(|ray| ray.coeffs)
        .collect()
}
