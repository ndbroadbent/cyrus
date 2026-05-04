use malachite::{Integer, Rational};

use crate::integer_math::invert_matrix;

use super::order::order_remaining_hyperplanes;
use super::rank::{RankContext, RankTracker};
use super::types::{DdmHyperplane, DdmRay};
use super::util::{active_set_for_existing_ray, rational_vector_to_i128};

pub(super) type DdmInitialization = (
    Vec<DdmRay>,
    Vec<DdmHyperplane>,
    usize,
    RankTracker,
    RankContext,
);

pub(super) fn initialize_ddm(hyperplanes: &[Vec<i128>], dim: usize) -> Option<DdmInitialization> {
    let mut rank = RankTracker::new(dim);
    let mut selected = vec![false; hyperplanes.len()];
    let mut basis_rows = Vec::with_capacity(dim);

    for (idx, row) in hyperplanes.iter().enumerate() {
        let before = rank.rank();
        rank.add_row(row);
        if rank.rank() > before {
            selected[idx] = true;
            basis_rows.push(row.clone());
            if rank.rank() == dim {
                break;
            }
        }
    }

    if basis_rows.len() != dim {
        return None;
    }

    let basis_rational: Vec<Vec<Rational>> = basis_rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|&x| Rational::from(Integer::from(x)))
                .collect()
        })
        .collect();
    let inv = invert_matrix(&basis_rational)?;
    let rank_context = RankContext::with_basis_inverse(&inv)?;

    let mut rays = Vec::with_capacity(dim);
    for col in 0..dim {
        let column = inv.iter().map(|row| row[col].clone()).collect::<Vec<_>>();
        let coeffs = rational_vector_to_i128(&column)?;
        let active = active_set_for_existing_ray(&coeffs, &basis_rows);
        rays.push(DdmRay { coeffs, active });
    }

    let mut ordered_hyperplanes: Vec<DdmHyperplane> =
        basis_rows.into_iter().map(DdmHyperplane::new).collect();
    let (mut remaining, pruned) = order_remaining_hyperplanes(hyperplanes, &selected, &rays);
    if pruned > 0 {
        eprintln!(
            "[DEBUG] ddm init: pruned {pruned} constraints redundant against initial basis cone"
        );
    }
    ordered_hyperplanes.append(&mut remaining);

    let mut processed_rank = RankTracker::new(dim);
    for row in &ordered_hyperplanes[..dim] {
        processed_rank.add_row(row.dense());
    }

    Some((rays, ordered_hyperplanes, dim, processed_rank, rank_context))
}

pub(super) fn full_space_initialization(
    hyperplanes: &[Vec<i128>],
    dim: usize,
) -> DdmInitialization {
    let rays: Vec<DdmRay> = (0..dim)
        .flat_map(|i| {
            let mut pos = vec![0i128; dim];
            let mut neg = vec![0i128; dim];
            pos[i] = 1;
            neg[i] = -1;
            vec![
                DdmRay {
                    coeffs: pos,
                    active: Vec::new(),
                },
                DdmRay {
                    coeffs: neg,
                    active: Vec::new(),
                },
            ]
        })
        .collect();
    (
        rays,
        hyperplanes
            .iter()
            .cloned()
            .map(DdmHyperplane::new)
            .collect(),
        0,
        RankTracker::new(dim),
        RankContext::default(),
    )
}
