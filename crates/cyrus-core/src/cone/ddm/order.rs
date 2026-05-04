use std::env;
use std::sync::OnceLock;

use super::types::{DdmHyperplane, DdmRay};

pub(super) fn order_remaining_hyperplanes(
    hyperplanes: &[Vec<i128>],
    selected: &[bool],
    initial_rays: &[DdmRay],
) -> (Vec<DdmHyperplane>, usize) {
    let mut pruned = 0usize;
    let mut remaining = Vec::new();
    for (idx, row) in hyperplanes.iter().enumerate() {
        if selected[idx] {
            continue;
        }
        let hyperplane = DdmHyperplane::new(row.clone());
        let Some(order) = remaining_hyperplane_order(&hyperplane, initial_rays, idx) else {
            pruned += 1;
            continue;
        };
        remaining.push((order, hyperplane));
    }
    remaining.sort_by(|(left, _), (right, _)| left.cmp(right));
    (
        remaining
            .into_iter()
            .map(|(_, hyperplane)| hyperplane)
            .collect(),
        pruned,
    )
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct RemainingHyperplaneOrder {
    category: u8,
    pair_count: usize,
    sparse_len: usize,
    original_idx: usize,
}

fn remaining_hyperplane_order(
    hyperplane: &DdmHyperplane,
    initial_rays: &[DdmRay],
    original_idx: usize,
) -> Option<RemainingHyperplaneOrder> {
    let SignCounts { positive, negative } = count_ray_signs(hyperplane, initial_rays);

    if negative == 0 {
        return None;
    }

    Some(RemainingHyperplaneOrder {
        // Facet cuts with no initial positive ray collapse the current cone to
        // a face and should happen before cuts that create intersections.
        category: u8::from(positive > 0),
        pair_count: positive * negative,
        sparse_len: hyperplane.sparse.len(),
        original_idx,
    })
}

pub(super) fn choose_next_hyperplane(
    hyperplanes: &[DdmHyperplane],
    start_idx: usize,
    rays: &[DdmRay],
) -> usize {
    choose_next_hyperplane_with_window(
        hyperplanes,
        start_idx,
        rays,
        ddm_order_window(),
        ddm_order_max_window(),
        ddm_order_expand_pair_threshold(),
    )
}

pub(super) fn choose_next_hyperplane_with_window(
    hyperplanes: &[DdmHyperplane],
    start_idx: usize,
    rays: &[DdmRay],
    base_window: usize,
    max_window: usize,
    expand_pair_threshold: usize,
) -> usize {
    if base_window <= 1 {
        return start_idx;
    }

    let max_window = max_window.max(base_window);
    let max_end = (start_idx + max_window).min(hyperplanes.len());
    let mut end = (start_idx + base_window).min(hyperplanes.len());
    let (mut best_idx, mut best_order) =
        best_hyperplane_in_window(hyperplanes, start_idx, end, rays);

    while best_order.category != 0 && best_order.pair_count > expand_pair_threshold && end < max_end
    {
        let current_width = end - start_idx;
        end = (start_idx + current_width * 2).min(max_end);
        (best_idx, best_order) = best_hyperplane_in_window(hyperplanes, start_idx, end, rays);
    }

    best_idx
}

fn best_hyperplane_in_window(
    hyperplanes: &[DdmHyperplane],
    start_idx: usize,
    end: usize,
    rays: &[DdmRay],
) -> (usize, CurrentHyperplaneOrder) {
    (start_idx..end)
        .map(|idx| {
            (
                idx,
                current_hyperplane_order(&hyperplanes[idx], rays, idx - start_idx),
            )
        })
        .min_by(|(_, left), (_, right)| left.cmp(right))
        .expect("window always includes start_idx")
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SignCounts {
    positive: usize,
    negative: usize,
}

fn count_ray_signs(hyperplane: &DdmHyperplane, rays: &[DdmRay]) -> SignCounts {
    let mut positive = 0usize;
    let mut negative = 0usize;
    for ray in rays {
        match hyperplane.dot_sign(&ray.coeffs) {
            std::cmp::Ordering::Greater => positive += 1,
            std::cmp::Ordering::Equal => {}
            std::cmp::Ordering::Less => negative += 1,
        }
    }
    SignCounts { positive, negative }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct CurrentHyperplaneOrder {
    category: u8,
    pair_count: usize,
    sparse_len: usize,
    offset: usize,
}

fn current_hyperplane_order(
    hyperplane: &DdmHyperplane,
    rays: &[DdmRay],
    offset: usize,
) -> CurrentHyperplaneOrder {
    let SignCounts { positive, negative } = count_ray_signs(hyperplane, rays);
    let category = match (positive, negative) {
        (_, 0) => 2,
        (0, _) => 0,
        _ => 1,
    };
    CurrentHyperplaneOrder {
        category,
        pair_count: positive * negative,
        sparse_len: hyperplane.sparse.len(),
        offset,
    }
}

fn ddm_order_window() -> usize {
    static WINDOW: OnceLock<usize> = OnceLock::new();
    *WINDOW.get_or_init(|| {
        env::var("CYRUS_DDM_ORDER_WINDOW")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(256)
    })
}

fn ddm_order_max_window() -> usize {
    static WINDOW: OnceLock<usize> = OnceLock::new();
    *WINDOW.get_or_init(|| {
        env::var("CYRUS_DDM_ORDER_MAX_WINDOW")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(8192)
    })
}

fn ddm_order_expand_pair_threshold() -> usize {
    static THRESHOLD: OnceLock<usize> = OnceLock::new();
    *THRESHOLD.get_or_init(|| {
        env::var("CYRUS_DDM_ORDER_EXPAND_PAIR_THRESHOLD")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(100_000)
    })
}
