//! Chamber search: walk the secondary fan toward GV-covered phases.
//!
//! The full KKLT verification selects the small toric curves at its phase-1
//! Kahler point and needs a computed GV invariant for every one of them.
//! Coverage is chamber-dependent: the toric two-face formulas cover the
//! curves of the chamber's own Mori cap, and a hostile chamber can leave
//! hundreds of exotic circuits uncovered (observed: 166 on a GA candidate's
//! default chamber, degrees to 130). McAllister hand-picked heights per
//! example; this module automates that choice by best-first search over
//! bistellar flips, scoring each neighbor chamber by how many of its
//! selected small curves the two-face formulas miss.
//!
//! Scoring approximation (documented, conservative in the right direction):
//! the selection Kahler point is held FIXED across neighbors (it is
//! expressed in the chamber-independent divisor basis), while the Mori cap
//! and two-face coverage are recomputed per chamber. The runner's extra
//! rescue paths (minimal-degree, general fallback, extremal-wall rays) are
//! not modelled, so a chamber scored zero may still need them - but the
//! definitive check is always the full runner re-run in the winning
//! chamber, which fails loudly if the score misled.

use std::collections::{BTreeSet, BinaryHeap};

use crate::Point;
use crate::error::{Error, Result};
use crate::gv::{
    CurvePruningStrategy, compute_mori_cone_cap_rays, compute_toric_two_face_curve_gv_invariants,
    prune_decomposable_curve_candidates, subcutoff_toric_curve_candidates,
};
use crate::polytope::Polytope;
use crate::triangulation::{
    Triangulation, flip_circuit_in_triangulation, secondary_cone_hyperplanes_native,
    triangulation_heights_from_secondary_cone,
};
use crate::types::f64::F64;
use crate::types::tags::{Finite, Pos};

/// Two-face GV coverage of one chamber's selected small curves.
#[derive(Debug, Clone)]
pub struct ChamberScore {
    /// Selected (subcutoff, pruned) small-curve classes.
    pub selected: usize,
    /// Selected classes with no toric two-face GV value.
    pub uncovered: usize,
    /// The uncovered classes themselves (ambient coordinates).
    pub uncovered_classes: Vec<Vec<i64>>,
}

/// Result of a chamber search.
#[derive(Debug, Clone)]
pub struct ChamberSearchOutcome {
    /// Simplices of the best chamber found.
    pub simplices: Vec<Vec<usize>>,
    /// Interior heights certifying the best chamber (secondary-cone point).
    pub heights: Vec<f64>,
    /// Uncovered-curve count in the best chamber.
    pub uncovered: usize,
    /// Uncovered-curve count in the starting chamber.
    pub start_uncovered: usize,
    /// Flip distance of the best chamber from the start.
    pub flips_from_start: usize,
    /// Chambers expanded (scored neighbor sets) before stopping.
    pub expanded: usize,
}

/// Score one chamber: select subcutoff curves of ITS Mori cap at the
/// supplied Kahler point, prune decomposables, and count how many of the
/// survivors the chamber's toric two-face formulas do not cover.
///
/// # Errors
/// Returns an error if the Mori cap, selection, pruning, or two-face GV
/// computation fails for this chamber.
pub fn score_chamber_two_face_coverage(
    polytope: &Polytope,
    points: &[Point],
    tri: &Triangulation,
    basis: &[usize],
    selection_t: &[F64<Finite>],
    cutoff: F64<Pos>,
    pruning: CurvePruningStrategy,
) -> Result<ChamberScore> {
    let ambient_rays = compute_mori_cone_cap_rays(tri, points, polytope, false, false, None)?;
    let candidates = subcutoff_toric_curve_candidates(&ambient_rays, basis, selection_t, cutoff)?;
    let pruned = prune_decomposable_curve_candidates(&candidates, pruning)?;
    let two_face = compute_toric_two_face_curve_gv_invariants(tri, points, polytope)?;
    let covered: BTreeSet<&[i64]> = two_face.iter().map(|inv| inv.class.as_slice()).collect();
    let uncovered_classes: Vec<Vec<i64>> = pruned
        .iter()
        .filter(|candidate| !covered.contains(candidate.class.as_slice()))
        .map(|candidate| candidate.class.clone())
        .collect();
    Ok(ChamberScore {
        selected: pruned.len(),
        uncovered: uncovered_classes.len(),
        uncovered_classes,
    })
}

/// A triangulation is an FRST candidate when it is fine (every point used)
/// and star (the origin in every simplex). Regularity is certified
/// separately via the secondary cone.
fn is_fine_star(simplices: &[Vec<usize>], n_points: usize, origin_idx: usize) -> bool {
    let mut used = vec![false; n_points];
    for simplex in simplices {
        if !simplex.contains(&origin_idx) {
            return false;
        }
        for &idx in simplex {
            if idx >= n_points {
                return false;
            }
            used[idx] = true;
        }
    }
    used.iter().all(|&u| u)
}

fn canonical_key(simplices: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut key: Vec<Vec<usize>> = simplices
        .iter()
        .map(|s| {
            let mut s = s.clone();
            s.sort_unstable();
            s
        })
        .collect();
    key.sort();
    key
}

/// Best-first search over bistellar flips for a chamber whose selected
/// small curves are fully covered by the toric two-face formulas.
///
/// Expands up to `max_expansions` chambers (scoring every valid FRST
/// neighbor of each) and stops early on a fully covered chamber. Returns
/// the best chamber found either way; the caller decides whether a nonzero
/// residual is acceptable. Fails loudly if even the START chamber cannot
/// be scored.
///
/// `reachable(neighbor_tri, neighbor_heights)` must decide whether the
/// KKLT solve still converges with the neighbor chamber's own
/// intersection numbers - coverage without reachability is useless
/// (observed: an 11-flip chamber covered 27 more curves but evicted the
/// phase-1 solution entirely). Unreachable neighbors are pruned from the
/// search frontier.
///
/// # Errors
/// Returns an error if the starting chamber cannot be scored or certified.
#[allow(clippy::too_many_arguments)] // one orchestration entry point
pub fn search_covered_chamber(
    polytope: &Polytope,
    points: &[Point],
    start_simplices: &[Vec<usize>],
    basis: &[usize],
    selection_t: &[F64<Finite>],
    cutoff: F64<Pos>,
    pruning: CurvePruningStrategy,
    max_expansions: usize,
    reachable: &(dyn Fn(&Triangulation, &[f64]) -> bool + Sync),
) -> Result<ChamberSearchOutcome> {
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .ok_or_else(|| Error::InvalidInput("origin not among points".into()))?;

    let start_tri = Triangulation::new(start_simplices.to_vec());
    let start_score = score_chamber_two_face_coverage(
        polytope,
        points,
        &start_tri,
        basis,
        selection_t,
        cutoff,
        pruning,
    )?;

    let mut visited: BTreeSet<Vec<Vec<usize>>> = BTreeSet::new();
    visited.insert(canonical_key(start_simplices));

    // Max-heap on Reverse(uncovered, depth): fewest uncovered first, then
    // shallowest. Entries carry their simplices and whether the expensive
    // reachability check has run. Reachability is LAZY: scoring thousands
    // of neighbors with a full KKLT solve each made deep searches cost
    // hours, but only chambers we actually adopt or expand need the check
    // - so it runs at pop time (O(expansions), not O(scored)).
    type QueueEntry = (std::cmp::Reverse<(usize, usize)>, Vec<Vec<usize>>, bool);
    let mut queue: BinaryHeap<QueueEntry> = BinaryHeap::new();
    queue.push((
        std::cmp::Reverse((start_score.uncovered, 0)),
        start_simplices.to_vec(),
        true, // the start chamber is where the failure happened: reachable
    ));

    let mut best = (
        start_score.uncovered,
        start_simplices.to_vec(),
        0usize, // depth
    );
    let mut expanded = 0usize;

    while let Some((std::cmp::Reverse((uncovered, depth)), simplices, checked)) = queue.pop() {
        if !checked {
            let tri = Triangulation::new(simplices.clone());
            let Ok(Some(heights)) = triangulation_heights_from_secondary_cone(points, &tri) else {
                continue;
            };
            if !reachable(&tri, &heights) {
                continue; // unreachable chamber: drop, never adopt or expand
            }
        }
        if uncovered < best.0 {
            best = (uncovered, simplices.clone(), depth);
        }
        if uncovered == 0 {
            break;
        }
        if expanded >= max_expansions {
            break;
        }
        expanded += 1;

        let tri = Triangulation::new(simplices.clone());
        let circuits = secondary_cone_hyperplanes_native(points, &tri)?;

        // Serial phase: flips are cheap. Enumerate, filter to fresh valid
        // FRST candidates, and mark them visited before the parallel phase.
        let mut seen_circuits: BTreeSet<Vec<i64>> = BTreeSet::new();
        let mut fresh_neighbors: Vec<Vec<Vec<usize>>> = Vec::new();
        for circuit_dense in circuits {
            if !seen_circuits.insert(circuit_dense.clone()) {
                continue;
            }
            let circuit_sparse: Vec<(usize, i64)> = circuit_dense
                .iter()
                .enumerate()
                .filter(|&(_, &c)| c != 0)
                .map(|(idx, &c)| (idx, c))
                .collect();
            let Ok(flip) = flip_circuit_in_triangulation(&simplices, &circuit_sparse) else {
                continue; // not a flippable wall of this chamber
            };
            let neighbor = flip.simplices;
            if !is_fine_star(&neighbor, points.len(), origin_idx) {
                continue;
            }
            let key = canonical_key(&neighbor);
            if visited.contains(&key) {
                continue;
            }
            visited.insert(key);
            fresh_neighbors.push(neighbor);
        }

        // Parallel phase: regularity certificates and coverage scoring are
        // the expensive parts; each neighbor is independent.
        use rayon::prelude::*;
        let scored: Vec<(usize, Vec<Vec<usize>>)> = fresh_neighbors
            .into_par_iter()
            .filter_map(|neighbor| {
                let neighbor_tri = Triangulation::new(neighbor.clone());
                // Only certified-regular triangulations are chambers.
                match triangulation_heights_from_secondary_cone(points, &neighbor_tri) {
                    Ok(Some(_)) => {}
                    _ => return None,
                }
                score_chamber_two_face_coverage(
                    polytope,
                    points,
                    &neighbor_tri,
                    basis,
                    selection_t,
                    cutoff,
                    pruning,
                )
                .ok() // unscorable neighbor: skip, never guess
                .map(|score| (score.uncovered, neighbor))
            })
            .collect();
        for (uncovered, neighbor) in scored {
            queue.push((std::cmp::Reverse((uncovered, depth + 1)), neighbor, false));
        }
    }

    let best_tri = Triangulation::new(best.1.clone());
    let heights = triangulation_heights_from_secondary_cone(points, &best_tri)?
        .ok_or_else(|| Error::InvalidInput("best chamber lost regularity certificate".into()))?;
    Ok(ChamberSearchOutcome {
        simplices: best.1,
        heights,
        uncovered: best.0,
        start_uncovered: start_score.uncovered,
        flips_from_start: best.2,
        expanded,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fine_star_detects_missing_point_and_non_star() {
        // 3 points, origin = 0.
        assert!(is_fine_star(&[vec![0, 1], vec![0, 2]], 3, 0));
        // point 2 unused -> not fine
        assert!(!is_fine_star(&[vec![0, 1]], 3, 0));
        // simplex without origin -> not star
        assert!(!is_fine_star(&[vec![0, 1], vec![1, 2]], 3, 0));
    }

    #[test]
    fn canonical_key_is_order_invariant() {
        let a = canonical_key(&[vec![2, 0, 1], vec![0, 3, 2]]);
        let b = canonical_key(&[vec![0, 2, 3], vec![0, 1, 2]]);
        assert_eq!(a, b);
    }
}
