//! General (HKTY-series) GV fallback for selected primal small curves the
//! toric two-face formulas and the minimal-degree pre-pass leave uncovered.
//!
//! This is the opt-in `--primal-gv-min-points` / `--primal-gv-max-deg` layer:
//! it computes GV invariants in the divisor basis up to the requested degree
//! and maps them back to ambient curve classes, covering exactly the missing
//! classes the cheaper layers could not reach.

use std::collections::HashMap;

use malachite::Integer;

use crate::{
    compute_grading_vector, compute_gv_invariants, compute_mori_cone_cap_rays,
    map_basis_gv_invariants_to_ambient,
};

use super::gv_basis::vector_gv_basis_data;
use super::{PrimalGeom, PrimalIntersection};

/// Degree statistics for the missing primal GV classes, in the GV grading.
#[derive(Clone, Debug, PartialEq, Eq)]
struct RequiredGvDegreeSummary {
    count: usize,
    min_degree: i128,
    max_degree: i128,
    first_over_max: Option<Vec<i64>>,
}

/// Compute exact general GV invariants for the missing ambient classes via the
/// divisor-basis HKTY series, then map them back to ambient curve coordinates.
///
/// Requires at least one of `min_points` / `max_deg`; `max_deg` must cover the
/// maximal grading degree among the missing classes (otherwise the series
/// would silently drop curves, so this fails loudly instead).
pub fn compute_primal_general_gv_by_ambient_class(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    required_ambient_classes: &[Vec<i64>],
    ambient_mori_rays: Option<&[Vec<i64>]>,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<HashMap<Vec<i64>, Integer>, String> {
    if min_points.is_none() && max_deg.is_none() {
        return Err(
            "primal general GV fallback requires --primal-gv-min-points or --primal-gv-max-deg"
                .into(),
        );
    }

    let computed_ambient_rays;
    let ambient_rays = if let Some(ambient_mori_rays) = ambient_mori_rays {
        ambient_mori_rays
    } else {
        computed_ambient_rays = compute_mori_cone_cap_rays(
            &geom.triangulation,
            &geom.triangulation_points,
            &geom.polytope,
            false,
            false,
            None,
        )
        .map_err(|e| format!("failed to compute primal ambient Mori-cap rays: {e}"))?;
        &computed_ambient_rays
    };
    let gv_basis_data = vector_gv_basis_data(
        ambient_rays,
        &intersection.linrels,
        &intersection.basis,
        "primal",
    )?;
    let rays = &gv_basis_data.mori_rays;
    let grading = compute_grading_vector(rays)
        .ok_or_else(|| "failed to compute primal GV grading vector".to_string())?;
    let degree_summary = summarize_required_gv_degrees(
        required_ambient_classes,
        &intersection.basis,
        &grading,
        max_deg,
    )?;
    if degree_summary.count > 0 {
        eprintln!(
            "[INFO] missing primal GV curve degree range: min={} max={} count={}",
            degree_summary.min_degree, degree_summary.max_degree, degree_summary.count
        );
    }
    if let Some(max_deg) = max_deg
        && degree_summary.max_degree > i128::from(max_deg)
    {
        return Err(format!(
            "primal general GV max_deg={max_deg} cannot cover all selected missing curves: required degree range {}..{} ({} curves), first_over_max={:?}",
            degree_summary.min_degree,
            degree_summary.max_degree,
            degree_summary.count,
            degree_summary.first_over_max
        ));
    }
    let general_gvs = compute_gv_invariants(
        rays,
        &grading,
        &gv_basis_data.q_matrix,
        &intersection.kappa_basis,
        min_points,
        max_deg,
    )
    .map_err(|e| format!("failed to compute primal general GV invariants: {e}"))?;
    let ambient_gvs =
        map_basis_gv_invariants_to_ambient(&general_gvs, &gv_basis_data.curve_basis_matrix)
            .map_err(|e| format!("failed to map primal general GV invariants to ambient: {e}"))?;

    let mut out = HashMap::with_capacity(ambient_gvs.len());
    for (class, gv) in ambient_gvs {
        match out.entry(class) {
            std::collections::hash_map::Entry::Occupied(existing) => {
                if existing.get() != &gv {
                    return Err(format!(
                        "conflicting general GV values for duplicate ambient curve: {} vs {gv}",
                        existing.get()
                    ));
                }
            }
            std::collections::hash_map::Entry::Vacant(slot) => {
                slot.insert(gv);
            }
        }
    }
    Ok(out)
}

/// Summarize the GV grading degrees of the missing ambient classes, flagging
/// the first class whose degree exceeds `max_deg`.
fn summarize_required_gv_degrees(
    ambient_classes: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    max_deg: Option<u32>,
) -> Result<RequiredGvDegreeSummary, String> {
    if basis.len() != grading.len() {
        return Err(format!(
            "basis length {} does not match grading vector length {}",
            basis.len(),
            grading.len()
        ));
    }
    if ambient_classes.is_empty() {
        return Ok(RequiredGvDegreeSummary {
            count: 0,
            min_degree: 0,
            max_degree: 0,
            first_over_max: None,
        });
    }

    let mut min_degree = i128::MAX;
    let mut max_degree = i128::MIN;
    let mut first_over_max = None;
    for class in ambient_classes {
        let degree = ambient_curve_grading_degree(class, basis, grading)?;
        min_degree = min_degree.min(degree);
        max_degree = max_degree.max(degree);
        if degree <= 0 {
            return Err(format!(
                "required primal GV curve has non-positive grading degree {degree}: {class:?}"
            ));
        }
        if let Some(max_deg) = max_deg
            && degree > i128::from(max_deg)
            && first_over_max.is_none()
        {
            first_over_max = Some(class.clone());
        }
    }

    Ok(RequiredGvDegreeSummary {
        count: ambient_classes.len(),
        min_degree,
        max_degree,
        first_over_max,
    })
}

/// Grading degree of an ambient curve class: the basis-restricted dot product
/// of the class with the GV grading vector.
fn ambient_curve_grading_degree(
    ambient_class: &[i64],
    basis: &[usize],
    grading: &[i64],
) -> Result<i128, String> {
    let mut degree = 0i128;
    for (&idx, &weight) in basis.iter().zip(grading.iter()) {
        let Some(&coefficient) = ambient_class.get(idx) else {
            return Err(format!(
                "basis index {idx} is out of bounds for ambient curve dimension {}",
                ambient_class.len()
            ));
        };
        degree += i128::from(coefficient) * i128::from(weight);
    }
    Ok(degree)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn required_gv_degree_summary_projects_ambient_classes_to_basis() {
        let classes = vec![vec![99, 2, 0], vec![99, 0, 3]];
        let summary = summarize_required_gv_degrees(&classes, &[1, 2], &[5, 7], Some(10)).unwrap();

        assert_eq!(
            summary,
            RequiredGvDegreeSummary {
                count: 2,
                min_degree: 10,
                max_degree: 21,
                first_over_max: Some(vec![99, 0, 3]),
            }
        );
    }

    #[test]
    fn required_gv_degree_summary_empty_is_zero() {
        let summary = summarize_required_gv_degrees(&[], &[1], &[1], None).unwrap();
        assert_eq!(summary.count, 0);
        assert_eq!(summary.first_over_max, None);
    }

    #[test]
    fn required_gv_degree_summary_rejects_invalid_projection_inputs() {
        assert!(
            summarize_required_gv_degrees(&[vec![1, 2]], &[0, 1], &[1], None)
                .unwrap_err()
                .contains("basis length")
        );
        assert!(
            summarize_required_gv_degrees(&[vec![1]], &[1], &[1], None)
                .unwrap_err()
                .contains("out of bounds")
        );
        assert!(
            summarize_required_gv_degrees(&[vec![0, -1]], &[1], &[1], None)
                .unwrap_err()
                .contains("non-positive")
        );
    }
}
