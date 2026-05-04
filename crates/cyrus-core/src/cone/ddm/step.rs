use std::collections::HashMap;
use std::env;
use std::sync::OnceLock;

use super::rank::{RankContext, are_adjacent_cached};
use super::types::{DdmHyperplane, DdmRay, RankStats, StepStats};
use super::util::{
    active_set_for_intersection, gcd_vec, insert_active, merge_active,
    normalize_ray_preserving_direction,
};

pub(super) fn process_hyperplane(
    rays: &[DdmRay],
    h: &DdmHyperplane,
    hyperplane_idx: usize,
    processed_hyperplanes: &[DdmHyperplane],
    dim: usize,
    pointed_before: bool,
    rank_context: &mut RankContext,
    adjacency_cache: &mut HashMap<Vec<usize>, bool>,
) -> (Vec<DdmRay>, StepStats) {
    let (positive, zero, negative) = partition_rays(rays, h, hyperplane_idx);
    let mut new_rays = surviving_rays(&positive, zero);
    let mut stats = StepStats::from_partitions(
        positive.len(),
        new_rays.len() - positive.len(),
        negative.len(),
    );
    if negative.is_empty() {
        stats.redundant = true;
        return (rays.to_vec(), stats);
    }
    log_step_start(hyperplane_idx, rays.len(), &stats);
    add_intersection_rays(
        &mut new_rays,
        &mut stats,
        &positive,
        &negative,
        h,
        hyperplane_idx,
        processed_hyperplanes,
        dim,
        pointed_before,
        rank_context,
        adjacency_cache,
    );
    stats.rank_cache_entries = adjacency_cache.len();

    (new_rays, stats)
}

fn partition_rays(
    rays: &[DdmRay],
    h: &DdmHyperplane,
    hyperplane_idx: usize,
) -> (Vec<DdmRay>, Vec<DdmRay>, Vec<DdmRay>) {
    let mut positive = Vec::new();
    let mut zero = Vec::new();
    let mut negative = Vec::new();

    for ray in rays {
        match h.dot(&ray.coeffs).cmp(&0) {
            std::cmp::Ordering::Greater => positive.push(ray.clone()),
            std::cmp::Ordering::Equal => {
                let mut ray = ray.clone();
                insert_active(&mut ray.active, hyperplane_idx);
                zero.push(ray);
            }
            std::cmp::Ordering::Less => negative.push(ray.clone()),
        }
    }

    (positive, zero, negative)
}

fn surviving_rays(positive: &[DdmRay], zero: Vec<DdmRay>) -> Vec<DdmRay> {
    let mut out = positive.to_vec();
    out.extend(zero);
    out
}

fn log_step_start(hyperplane_idx: usize, rays_in: usize, stats: &StepStats) {
    if stats.pairs >= step_start_pair_threshold() {
        eprintln!(
            "[DEBUG] ddm step start: hyperplane={}, rays_in={}, pos={}, zero={}, neg={}, pairs={}",
            hyperplane_idx + 1,
            rays_in,
            stats.positive,
            stats.zero,
            stats.negative,
            stats.pairs
        );
    }
}

fn step_start_pair_threshold() -> usize {
    static THRESHOLD: OnceLock<usize> = OnceLock::new();
    *THRESHOLD.get_or_init(|| {
        env::var("CYRUS_DDM_STEP_START_PAIR_THRESHOLD")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(100_000)
    })
}

#[allow(clippy::too_many_arguments)]
fn add_intersection_rays(
    new_rays: &mut Vec<DdmRay>,
    stats: &mut StepStats,
    positive: &[DdmRay],
    negative: &[DdmRay],
    h: &DdmHyperplane,
    hyperplane_idx: usize,
    processed_hyperplanes: &[DdmHyperplane],
    dim: usize,
    pointed_before: bool,
    rank_context: &mut RankContext,
    adjacency_cache: &mut HashMap<Vec<usize>, bool>,
) {
    let mut processed_pairs = 0usize;
    for p in positive {
        for n in negative {
            processed_pairs += 1;
            if pair_is_redundant(
                p,
                n,
                processed_hyperplanes,
                dim,
                pointed_before,
                rank_context,
                adjacency_cache,
                &mut stats.rank,
            ) {
                log_pair_progress(stats, processed_pairs, adjacency_cache.len());
                continue;
            }
            stats.adjacent += 1;
            if let Some(intersection) = compute_intersection_ray(&p.coeffs, &n.coeffs, h) {
                let active = active_set_for_intersection(p, n, hyperplane_idx);
                new_rays.push(DdmRay {
                    coeffs: intersection,
                    active,
                });
                stats.created += 1;
            }
            log_pair_progress(stats, processed_pairs, adjacency_cache.len());
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn pair_is_redundant(
    p: &DdmRay,
    n: &DdmRay,
    processed_hyperplanes: &[DdmHyperplane],
    dim: usize,
    pointed_before: bool,
    rank_context: &mut RankContext,
    adjacency_cache: &mut HashMap<Vec<usize>, bool>,
    rank_stats: &mut RankStats,
) -> bool {
    pointed_before
        && !are_adjacent_cached(
            p,
            n,
            processed_hyperplanes,
            dim,
            rank_context,
            adjacency_cache,
            rank_stats,
        )
}

/// Compute the intersection ray for a (positive, negative) pair.
///
/// Given rays r+ and r- with h·r+ > 0 and h·r- < 0,
/// find the ray r on the hyperplane h·r = 0 that lies between them.
///
/// r = (h·r+) * r- - (h·r-) * r+
///
/// This makes h·r = (h·r+)(h·r-) - (h·r-)(h·r+) = 0.
pub(super) fn compute_intersection_ray(
    pos: &[i128],
    neg: &[i128],
    h: &DdmHyperplane,
) -> Option<Vec<i128>> {
    let dp = h.dot(pos); // > 0
    let dn = h.dot(neg); // < 0

    // r = dp * neg - dn * pos
    let mut result: Vec<i128> = pos
        .iter()
        .zip(neg.iter())
        .map(|(&p, &n)| dp * n - dn * p)
        .collect();

    // Normalize by GCD
    let g = gcd_vec(&result);
    if g == 0 {
        return None;
    }

    for x in &mut result {
        *x /= g;
    }

    Some(result)
}

/// Remove duplicate rays (after normalization).
pub(super) fn deduplicate_rays(rays: Vec<DdmRay>) -> Vec<DdmRay> {
    let mut by_coeffs: HashMap<Vec<i128>, Vec<usize>> = HashMap::new();
    for ray in rays {
        let normalized = normalize_ray_preserving_direction(ray.coeffs);
        if normalized.iter().all(|&x| x == 0) {
            continue;
        }
        by_coeffs
            .entry(normalized)
            .and_modify(|active| merge_active(active, &ray.active))
            .or_insert(ray.active);
    }

    by_coeffs
        .into_iter()
        .map(|(coeffs, active)| DdmRay { coeffs, active })
        .collect()
}

fn log_pair_progress(stats: &StepStats, processed_pairs: usize, rank_cache_entries: usize) {
    let Some(interval) = pair_log_interval() else {
        return;
    };
    if stats.pairs > interval && processed_pairs % interval == 0 {
        eprintln!(
            "[DEBUG] ddm pairs: {}/{} processed, adjacent={}, created={}, rank_cache={}, rank={}",
            processed_pairs,
            stats.pairs,
            stats.adjacent,
            stats.created,
            rank_cache_entries,
            rank_summary(&stats.rank)
        );
    }
}

fn pair_log_interval() -> Option<usize> {
    static INTERVAL: OnceLock<Option<usize>> = OnceLock::new();
    *INTERVAL.get_or_init(|| {
        env::var("CYRUS_DDM_PAIR_LOG_INTERVAL")
            .ok()
            .map_or(Some(100_000), |value| parse_pair_log_interval(&value))
    })
}

pub(super) fn parse_pair_log_interval(value: &str) -> Option<usize> {
    match value.parse::<usize>() {
        Ok(0) | Err(_) => None,
        Ok(interval) => Some(interval),
    }
}

pub(super) const fn should_log_step(idx: usize, stats: &StepStats) -> bool {
    if stats.redundant {
        return idx % 1000 == 0;
    }
    idx % 25 == 0 || stats.created > 10_000 || stats.pairs > 1_000
}

pub(super) fn log_step(idx: usize, total: usize, rays: usize, stats: &StepStats) {
    eprintln!(
        "[DEBUG] ddm: hyperplane {}/{}, rays={}, pos={}, zero={}, neg={}, pairs={}, adjacent={}, created={}, redundant={}, rank_cache={}, rank={}",
        idx + 1,
        total,
        rays,
        stats.positive,
        stats.zero,
        stats.negative,
        stats.pairs,
        stats.adjacent,
        stats.created,
        stats.redundant,
        stats.rank_cache_entries,
        rank_summary(&stats.rank)
    );
}

fn rank_summary(stats: &RankStats) -> String {
    format!(
        "checks={}, short={}, cache_hits={}, basis={}, mod={}, quotient_checks={}, quotient={}/{}, int_checks={}, int={}/{}",
        stats.checks,
        stats.too_short,
        stats.cache_hits,
        stats.basis_hits,
        stats.modular_hits,
        stats.quotient_checks,
        stats.quotient_true,
        stats.quotient_false,
        stats.integer_checks,
        stats.integer_true,
        stats.integer_false
    )
}
