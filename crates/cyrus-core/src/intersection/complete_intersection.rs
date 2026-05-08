//! Complete-intersection intersection-number reduction.
//!
//! This ports the complete-intersection branch of
//! `cytools.calabiyau.CalabiYau.intersection_numbers`: ambient top-form
//! intersections are successively intersected with each nef-partition divisor
//! class until a Calabi-Yau threefold triple-intersection tensor remains.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use malachite::Integer;
use malachite::Rational as MalachiteRational;
use malachite::num::basic::traits::Zero;

use crate::Point;
use crate::error::{Error, Result};
use crate::geometry::ConvexHull;
use crate::integer_math::matrix_rank;
use crate::triangulation::Triangulation;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;
use crate::utils::{gcd_integers, lll_reduce};

use super::Intersection;
use super::cytools_algorithm::compute_ambient_intersections_cytools;

/// CYTools-style nef-partition data certified from polytope checks.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FullDimensionalNefPartitionCertificate {
    /// Nef parts converted to triangulation-point indices, excluding origin.
    pub nef_parts: Vec<Vec<usize>>,
    /// Number of lattice points in each part polytope, including origin.
    pub part_lattice_point_counts: Vec<usize>,
    /// Number of vertices of the Minkowski sum of the part polytopes.
    pub minkowski_sum_vertex_count: usize,
}

/// Alias for the general CYTools-style nef certificate.
pub type CytoolsNefPartitionCertificate = FullDimensionalNefPartitionCertificate;

/// Certify a CYTools-style nef partition.
///
/// This ports the explicit checks in `CalabiYau.__init__`:
///
/// 1. the convex hull of the union of partition points equals the ambient hull;
/// 2. the Minkowski sum of the part polytopes is reflexive;
/// 3. each part is converted from part-polytope lattice points to
///    triangulation-point indices, excluding the origin.
///
/// Lower-dimensional part polytopes use CYTools' translate + LLL-reduce + drop
/// zero directions convention before lattice-point enumeration, then map the
/// resulting lattice points back to the original ambient coordinates.
///
/// # Errors
/// Returns an error if indices are invalid, the partition does not cover the
/// ambient hull, the affine reduction is inconsistent, or the Minkowski sum is
/// not reflexive.
pub fn certify_nef_partition_cytools_style(
    ambient_points: &[Point],
    triangulation_points: &[Point],
    nef_partition: &[Vec<usize>],
) -> Result<CytoolsNefPartitionCertificate> {
    validate_nef_partition_point_indices(ambient_points, nef_partition)?;
    validate_same_point_dimension(ambient_points, "ambient points")?;
    validate_same_point_dimension(triangulation_points, "triangulation points")?;
    if ambient_points.first().map(Point::dim) != triangulation_points.first().map(Point::dim) {
        return Err(Error::InvalidInput(
            "ambient and triangulation points have different dimensions".into(),
        ));
    }
    validate_nef_partition_union_hull(ambient_points, nef_partition)?;

    let certified =
        certify_nef_partition_parts_cytools(ambient_points, triangulation_points, nef_partition)?;

    let minkowski_vertices = minkowski_sum_cytools_vertices(&certified.part_vertex_sets)?;
    if !is_reflexive_cytools_polytope(&minkowski_vertices, "nef partition Minkowski sum")? {
        return Err(Error::InvalidInput(
            "nef partition Minkowski sum is not reflexive".into(),
        ));
    }

    Ok(FullDimensionalNefPartitionCertificate {
        nef_parts: certified.nef_parts,
        part_lattice_point_counts: certified.part_lattice_point_counts,
        minkowski_sum_vertex_count: minkowski_vertices.len(),
    })
}

/// Certify a CYTools-style nef partition, requiring every part to be full-dimensional.
///
/// This is a stricter wrapper around [`certify_nef_partition_cytools_style`]
/// used by tests and callers that need to avoid affine-span reduction.
///
/// # Errors
/// Returns an error if any nef part is lower-dimensional, or if the general
/// CYTools-style nef checks fail.
pub fn certify_nef_partition_cytools_style_full_dim(
    ambient_points: &[Point],
    triangulation_points: &[Point],
    nef_partition: &[Vec<usize>],
) -> Result<FullDimensionalNefPartitionCertificate> {
    validate_nef_partition_point_indices(ambient_points, nef_partition)?;
    let ambient_dim = ambient_points
        .first()
        .ok_or_else(|| Error::InvalidInput("ambient point set is empty".into()))?
        .dim();
    for (part_idx, part) in nef_partition.iter().enumerate() {
        let part_input = nef_part_input_points(ambient_points, part);
        if affine_rank_i64_points(&part_input) != ambient_dim {
            return Err(Error::InvalidInput(format!(
                "nef part {part_idx} polytope must be full-dimensional"
            )));
        }
    }
    certify_nef_partition_cytools_style(ambient_points, triangulation_points, nef_partition)
}

struct CertifiedCytoolsNefParts {
    nef_parts: Vec<Vec<usize>>,
    part_lattice_point_counts: Vec<usize>,
    part_vertex_sets: Vec<Vec<Point>>,
}

struct CytoolsPolytopeModel {
    dimension: usize,
    input_lattice_points: Vec<Point>,
    input_vertices: Vec<Point>,
    optimal_vertices: Vec<Point>,
}

/// Compute a CICY threefold triple-intersection tensor from ambient toric data.
///
/// This composes the CYTools generic ambient top-form solver with the CYTools
/// complete-intersection reduction. It still requires a source-derived
/// `nef_parts` input; it does not infer or choose a nef partition.
///
/// # Errors
/// Returns an error if ambient top-form computation or CICY reduction fails.
pub fn compute_complete_intersection_cy3_from_ambient_cytools(
    tri: &Triangulation,
    points: &[Point],
    linear_relations_no_origin: &[Vec<i64>],
    nef_parts: &[Vec<usize>],
) -> Result<Intersection> {
    let ambient = compute_ambient_intersections_cytools(tri, points, linear_relations_no_origin)?;
    let ambient_entries = ambient
        .iter()
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect::<Vec<_>>();
    compute_complete_intersection_cy3_intersection_numbers(
        &ambient_entries,
        nef_parts,
        ambient.divisor_count(),
    )
}

fn validate_nef_partition_point_indices(
    ambient_points: &[Point],
    nef_partition: &[Vec<usize>],
) -> Result<()> {
    if ambient_points.is_empty() {
        return Err(Error::InvalidInput("ambient point set is empty".into()));
    }
    if nef_partition.is_empty() {
        return Err(Error::InvalidInput("nef partition is empty".into()));
    }
    let mut seen = HashSet::new();
    for part in nef_partition {
        if part.is_empty() {
            return Err(Error::InvalidInput(
                "nef partition contains an empty part".into(),
            ));
        }
        for &idx in part {
            if idx == 0 {
                return Err(Error::InvalidInput(
                    "nef partition parts must exclude the origin index".into(),
                ));
            }
            if idx >= ambient_points.len() {
                return Err(Error::InvalidInput(format!(
                    "nef partition point index {idx} is out of range"
                )));
            }
            if !seen.insert(idx) {
                return Err(Error::InvalidInput(format!(
                    "nef partition point index {idx} appears in more than one part"
                )));
            }
        }
    }
    Ok(())
}

fn validate_nef_partition_union_hull(
    ambient_points: &[Point],
    nef_partition: &[Vec<usize>],
) -> Result<()> {
    let union_points = nef_partition
        .iter()
        .flat_map(|part| part.iter().copied())
        .map(|idx| ambient_points[idx].clone())
        .collect::<Vec<_>>();
    if cytools_polytope_vertex_set(&union_points, "nef partition union")?
        != cytools_polytope_vertex_set(ambient_points, "ambient polytope")?
    {
        return Err(Error::InvalidInput(
            "nef partition union hull does not equal ambient polytope hull".into(),
        ));
    }
    Ok(())
}

fn certify_nef_partition_parts_cytools(
    ambient_points: &[Point],
    triangulation_points: &[Point],
    nef_partition: &[Vec<usize>],
) -> Result<CertifiedCytoolsNefParts> {
    let triangulation_index_by_point = triangulation_point_index_by_coords(triangulation_points)?;
    let mut part_vertex_sets = Vec::with_capacity(nef_partition.len());
    let mut nef_parts = Vec::with_capacity(nef_partition.len());
    let mut part_lattice_point_counts = Vec::with_capacity(nef_partition.len());
    for (part_idx, part) in nef_partition.iter().enumerate() {
        let part_input = nef_part_input_points(ambient_points, part);
        let part_model = cytools_polytope_model(&part_input, &format!("nef part {part_idx}"))?;
        let part_lattice_points = part_model.input_lattice_points;
        part_lattice_point_counts.push(part_lattice_points.len());
        nef_parts.push(triangulation_indices_for_part_lattice_points(
            &part_lattice_points,
            &triangulation_index_by_point,
        ));
        part_vertex_sets.push(part_model.input_vertices);
    }
    Ok(CertifiedCytoolsNefParts {
        nef_parts,
        part_lattice_point_counts,
        part_vertex_sets,
    })
}

fn nef_part_input_points(ambient_points: &[Point], part: &[usize]) -> Vec<Point> {
    std::iter::once(0)
        .chain(part.iter().copied())
        .map(|idx| ambient_points[idx].clone())
        .collect()
}

fn validate_same_point_dimension(points: &[Point], name: &str) -> Result<()> {
    let Some(first) = points.first() else {
        return Err(Error::InvalidInput(format!("{name} are empty")));
    };
    let dim = first.dim();
    if points.iter().any(|point| point.dim() != dim) {
        return Err(Error::InvalidInput(format!(
            "{name} have inconsistent dimensions"
        )));
    }
    Ok(())
}

fn triangulation_point_index_by_coords(points: &[Point]) -> Result<HashMap<Vec<i64>, usize>> {
    let mut out = HashMap::new();
    for (idx, point) in points.iter().enumerate() {
        if out.insert(point.coords().to_vec(), idx).is_some() {
            return Err(Error::InvalidInput(
                "triangulation points contain duplicate coordinates".into(),
            ));
        }
    }
    Ok(out)
}

fn triangulation_indices_for_part_lattice_points(
    part_lattice_points: &[Point],
    triangulation_index_by_point: &HashMap<Vec<i64>, usize>,
) -> Vec<usize> {
    let mut out = part_lattice_points
        .iter()
        .filter(|point| point.coords().iter().any(|&coord| coord != 0))
        .filter_map(|point| triangulation_index_by_point.get(point.coords()).copied())
        .collect::<Vec<_>>();
    out.sort_unstable();
    out.dedup();
    out
}

fn cytools_polytope_vertex_set(points: &[Point], context: &str) -> Result<BTreeSet<Vec<i64>>> {
    Ok(cytools_polytope_model(points, context)?
        .input_vertices
        .into_iter()
        .map(|point| point.coords().to_vec())
        .collect())
}

fn cytools_polytope_model(points: &[Point], context: &str) -> Result<CytoolsPolytopeModel> {
    validate_same_point_dimension(points, context)?;
    let ambient_dim = points
        .first()
        .expect("validated point set is nonempty")
        .dim();
    let dimension = affine_rank_i64_points(points);
    if dimension == 0 {
        return Ok(cytools_zero_dim_polytope_model(points));
    }
    cytools_positive_dim_polytope_model(points, context, dimension, ambient_dim)
}

fn cytools_zero_dim_polytope_model(points: &[Point]) -> CytoolsPolytopeModel {
    let point = points
        .first()
        .expect("validated point set is nonempty")
        .clone();
    CytoolsPolytopeModel {
        dimension: 0,
        input_lattice_points: vec![point.clone()],
        input_vertices: vec![point],
        optimal_vertices: vec![Point::new(Vec::new())],
    }
}

fn cytools_positive_dim_polytope_model(
    points: &[Point],
    context: &str,
    dimension: usize,
    ambient_dim: usize,
) -> Result<CytoolsPolytopeModel> {
    let translation = cytools_affine_translation(points, dimension, ambient_dim);
    let translated = translate_points(points, &translation)?;
    let lll = lll_reduce(&translated, true);
    let dim_diff = ambient_dim
        .checked_sub(dimension)
        .ok_or_else(|| Error::InvalidInput(format!("{context} has impossible affine rank")))?;
    let optimal_points = cytools_optimal_points_from_lll(&lll.reduced, dim_diff, context)?;
    let optimal_lattice_points =
        lattice_points_in_full_dim_polytope(&optimal_points, &format!("{context} optimal"))?;
    let optimal_vertices =
        convex_hull_vertex_points(&optimal_points, &format!("{context} optimal"))?;
    let transform_inv = lll.transform_inv.ok_or_else(|| {
        Error::InvalidInput(format!("{context} LLL transform inverse was not computed"))
    })?;
    let input_lattice_points = cytools_optimal_points_to_input(
        &optimal_lattice_points,
        ambient_dim,
        dim_diff,
        &transform_inv,
        &translation,
    )?;
    let input_vertices = cytools_optimal_points_to_input(
        &optimal_vertices,
        ambient_dim,
        dim_diff,
        &transform_inv,
        &translation,
    )?;
    Ok(CytoolsPolytopeModel {
        dimension,
        input_lattice_points,
        input_vertices,
        optimal_vertices,
    })
}

fn affine_rank_i64_points(points: &[Point]) -> usize {
    let Some(base) = points.first() else {
        return 0;
    };
    let rows = points
        .iter()
        .skip(1)
        .map(|point| {
            point
                .coords()
                .iter()
                .zip(base.coords().iter())
                .map(|(&coord, &base_coord)| MalachiteRational::from(coord - base_coord))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    matrix_rank(&rows)
}

fn cytools_affine_translation(points: &[Point], dimension: usize, ambient_dim: usize) -> Vec<i64> {
    if dimension == ambient_dim {
        vec![0; ambient_dim]
    } else {
        points
            .first()
            .expect("validated point set is nonempty")
            .coords()
            .to_vec()
    }
}

fn translate_points(points: &[Point], translation: &[i64]) -> Result<Vec<Vec<i64>>> {
    points
        .iter()
        .map(|point| {
            if point.dim() != translation.len() {
                return Err(Error::InvalidInput(
                    "point dimension does not match affine translation".into(),
                ));
            }
            Ok(point
                .coords()
                .iter()
                .zip(translation.iter())
                .map(|(&coord, &shift)| coord - shift)
                .collect::<Vec<_>>())
        })
        .collect()
}

fn cytools_optimal_points_from_lll(
    reduced: &[Vec<i64>],
    dim_diff: usize,
    context: &str,
) -> Result<Vec<Point>> {
    let mut optimal = Vec::with_capacity(reduced.len());
    for row in reduced {
        if row.len() < dim_diff {
            return Err(Error::InvalidInput(format!(
                "{context} LLL row is shorter than dropped dimension count"
            )));
        }
        if row[..dim_diff].iter().any(|&coord| coord != 0) {
            return Err(Error::InvalidInput(format!(
                "{context} LLL reduction did not isolate affine zero directions"
            )));
        }
        optimal.push(Point::new(row[dim_diff..].to_vec()));
    }
    Ok(dedup_points(optimal))
}

fn dedup_points(points: Vec<Point>) -> Vec<Point> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for point in points {
        if seen.insert(point.coords().to_vec()) {
            out.push(point);
        }
    }
    out
}

fn cytools_optimal_points_to_input(
    optimal_points: &[Point],
    ambient_dim: usize,
    dim_diff: usize,
    transform_inv: &[Vec<i64>],
    translation: &[i64],
) -> Result<Vec<Point>> {
    if transform_inv.len() != ambient_dim
        || transform_inv.iter().any(|row| row.len() != ambient_dim)
    {
        return Err(Error::InvalidInput(
            "LLL transform inverse has wrong dimensions".into(),
        ));
    }
    optimal_points
        .iter()
        .map(|point| {
            if point.dim() + dim_diff != ambient_dim {
                return Err(Error::InvalidInput(
                    "optimal point dimension does not match affine embedding".into(),
                ));
            }
            let mut padded = vec![0i64; ambient_dim];
            for (idx, &coord) in point.coords().iter().enumerate() {
                padded[dim_diff + idx] = coord;
            }
            let coords = transform_inv
                .iter()
                .zip(translation.iter())
                .map(|(row, &shift)| {
                    row.iter()
                        .zip(padded.iter())
                        .map(|(&matrix_entry, &coord)| matrix_entry * coord)
                        .sum::<i64>()
                        + shift
                })
                .collect::<Vec<_>>();
            Ok(Point::new(coords))
        })
        .collect()
}

fn convex_hull_vertex_points(points: &[Point], context: &str) -> Result<Vec<Point>> {
    let coords = point_coords(points);
    let hull = ConvexHull::compute(&coords).ok_or_else(|| {
        Error::InvalidInput(format!(
            "{context} must be full-dimensional for Cyrus' current CYTools-style nef validation"
        ))
    })?;
    Ok(hull
        .vertex_indices
        .into_iter()
        .map(|idx| Point::new(coords[idx].clone()))
        .collect())
}

fn point_coords(points: &[Point]) -> Vec<Vec<i64>> {
    points.iter().map(|point| point.coords().to_vec()).collect()
}

fn lattice_points_in_full_dim_polytope(points: &[Point], context: &str) -> Result<Vec<Point>> {
    let coords = point_coords(points);
    let hull = ConvexHull::compute(&coords).ok_or_else(|| {
        Error::InvalidInput(format!(
            "{context} must be full-dimensional for lattice-point enumeration"
        ))
    })?;
    let (min_coords, max_coords) = bounding_box(&coords)?;
    let mut lattice_points = Vec::new();
    let mut candidate = min_coords.clone();
    loop {
        if point_satisfies_hull_facets(&candidate, &hull) {
            lattice_points.push(Point::new(candidate.clone()));
        }
        if !advance_lattice_box_point(&mut candidate, &min_coords, &max_coords) {
            break;
        }
    }
    Ok(lattice_points)
}

fn bounding_box(points: &[Vec<i64>]) -> Result<(Vec<i64>, Vec<i64>)> {
    let Some(first) = points.first() else {
        return Err(Error::InvalidInput("cannot bound empty point set".into()));
    };
    let mut min_coords = first.clone();
    let mut max_coords = first.clone();
    for point in &points[1..] {
        for (idx, &coord) in point.iter().enumerate() {
            min_coords[idx] = min_coords[idx].min(coord);
            max_coords[idx] = max_coords[idx].max(coord);
        }
    }
    Ok((min_coords, max_coords))
}

fn point_satisfies_hull_facets(point: &[i64], hull: &ConvexHull) -> bool {
    hull.facets.iter().all(|facet| {
        let dot = facet
            .normal
            .iter()
            .zip(point.iter())
            .map(|(normal, &coord)| normal * Integer::from(coord))
            .sum::<Integer>();
        dot <= facet.constant
    })
}

fn advance_lattice_box_point(
    candidate: &mut [i64],
    min_coords: &[i64],
    max_coords: &[i64],
) -> bool {
    for idx in 0..candidate.len() {
        candidate[idx] += 1;
        if candidate[idx] <= max_coords[idx] {
            return true;
        }
        candidate[idx] = min_coords[idx];
    }
    false
}

fn minkowski_sum_cytools_vertices(part_vertex_sets: &[Vec<Point>]) -> Result<Vec<Point>> {
    let Some(first_part) = part_vertex_sets.first() else {
        return Err(Error::InvalidInput(
            "Minkowski sum requires at least one part".into(),
        ));
    };
    let dim = first_part
        .first()
        .ok_or_else(|| Error::InvalidInput("Minkowski sum part has no vertices".into()))?
        .dim();
    let mut sums = vec![vec![0i64; dim]];
    for vertices in part_vertex_sets {
        if vertices.is_empty() {
            return Err(Error::InvalidInput(
                "Minkowski sum part has no vertices".into(),
            ));
        }
        let mut next = Vec::new();
        for sum in &sums {
            for vertex in vertices {
                if vertex.dim() != dim {
                    return Err(Error::InvalidInput(
                        "Minkowski sum vertices have inconsistent dimensions".into(),
                    ));
                }
                next.push(
                    sum.iter()
                        .zip(vertex.coords().iter())
                        .map(|(&lhs, &rhs)| lhs + rhs)
                        .collect::<Vec<_>>(),
                );
            }
        }
        next.sort();
        next.dedup();
        sums = next;
    }
    let sum_points = sums.into_iter().map(Point::new).collect::<Vec<_>>();
    Ok(cytools_polytope_model(&sum_points, "nef partition Minkowski sum")?.input_vertices)
}

fn is_reflexive_cytools_polytope(points: &[Point], context: &str) -> Result<bool> {
    let model = cytools_polytope_model(points, context)?;
    if model.dimension == 0 {
        return Ok(true);
    }
    is_reflexive_full_dim_by_facets(&model.optimal_vertices, &format!("{context} optimal"))
}

fn is_reflexive_full_dim_by_facets(points: &[Point], context: &str) -> Result<bool> {
    let coords = point_coords(points);
    let hull = ConvexHull::compute(&coords).ok_or_else(|| {
        Error::InvalidInput(format!(
            "{context} must be full-dimensional for reflexivity check"
        ))
    })?;
    for facet in &hull.facets {
        if facet.constant <= Integer::from(0) {
            return Ok(false);
        }
        let normal_gcd = gcd_integers(&facet.normal);
        let primitive_constant = &facet.constant / normal_gcd;
        if primitive_constant != Integer::from(1) {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Reduce ambient top-form intersections to a CICY threefold intersection tensor.
///
/// `ambient_intersections` must contain sorted ambient intersection tuples of
/// length `3 + nef_parts.len()`. Each nef part is a set of ambient divisor
/// indices whose sum defines one complete-intersection divisor class.
///
/// The reduction is the CYTools algorithm: for each nef part and each current
/// tuple, remove one occurrence of any divisor index belonging to that part,
/// canonicalize the remaining tuple, and add the ambient coefficient once per
/// distinct resulting tuple.
///
/// # Errors
/// Returns an error if tuple lengths, divisor indices, or nef partition data are
/// inconsistent.
pub fn compute_complete_intersection_cy3_intersection_numbers(
    ambient_intersections: &[(Vec<usize>, TypedRational<Finite>)],
    nef_parts: &[Vec<usize>],
    divisor_count: usize,
) -> Result<Intersection> {
    if divisor_count == 0 {
        return Err(Error::InvalidInput("divisor count is zero".into()));
    }
    if nef_parts.is_empty() {
        return Err(Error::InvalidInput(
            "complete-intersection reduction requires at least one nef part".into(),
        ));
    }
    validate_nef_parts(nef_parts, divisor_count)?;

    let ambient_tuple_len = nef_parts
        .len()
        .checked_add(3)
        .ok_or_else(|| Error::InvalidInput("ambient tuple length overflowed".into()))?;
    let mut current =
        ambient_intersection_map(ambient_intersections, ambient_tuple_len, divisor_count)?;
    let nef_part_sets = nef_parts
        .iter()
        .map(|part| part.iter().copied().collect::<HashSet<_>>())
        .collect::<Vec<_>>();
    for (step, part) in nef_part_sets.iter().enumerate() {
        let tuple_len = ambient_tuple_len
            .checked_sub(step)
            .expect("step never exceeds ambient tuple length");
        current = apply_nef_part_reduction(&current, part, tuple_len)?;
    }
    intersection_from_reduced_triples(current, divisor_count)
}

fn ambient_intersection_map(
    ambient_intersections: &[(Vec<usize>, TypedRational<Finite>)],
    ambient_tuple_len: usize,
    divisor_count: usize,
) -> Result<BTreeMap<Vec<usize>, MalachiteRational>> {
    let mut current = BTreeMap::<Vec<usize>, MalachiteRational>::new();
    for (indices, value) in ambient_intersections {
        validate_ambient_tuple(indices, ambient_tuple_len, divisor_count)?;
        if *value.get() == MalachiteRational::ZERO {
            continue;
        }
        let mut key = indices.clone();
        key.sort_unstable();
        let entry = current
            .entry(key)
            .or_insert_with(|| MalachiteRational::ZERO);
        *entry += value.get().clone();
    }
    Ok(current)
}

fn apply_nef_part_reduction(
    current: &BTreeMap<Vec<usize>, MalachiteRational>,
    part: &HashSet<usize>,
    tuple_len: usize,
) -> Result<BTreeMap<Vec<usize>, MalachiteRational>> {
    let mut next = BTreeMap::<Vec<usize>, MalachiteRational>::new();
    for (indices, value) in current {
        if indices.len() != tuple_len {
            return Err(Error::InvalidInput(format!(
                "complete-intersection reduction expected tuple length {tuple_len}, got {}",
                indices.len()
            )));
        }
        for choice in remove_one_nef_part_index_choices(indices, part) {
            let entry = next
                .entry(choice)
                .or_insert_with(|| MalachiteRational::ZERO);
            *entry += value.clone();
        }
    }
    next.retain(|_, value| *value != MalachiteRational::ZERO);
    Ok(next)
}

fn intersection_from_reduced_triples(
    triples: BTreeMap<Vec<usize>, MalachiteRational>,
    divisor_count: usize,
) -> Result<Intersection> {
    let mut out = Intersection::new(divisor_count);
    for (indices, value) in triples {
        if indices.len() != 3 {
            return Err(Error::InvalidInput(format!(
                "complete-intersection reduction should produce triples, got length {}",
                indices.len()
            )));
        }
        out.set(
            indices[0],
            indices[1],
            indices[2],
            TypedRational::<Finite>::from_raw(value),
        );
    }
    Ok(out)
}

fn validate_nef_parts(nef_parts: &[Vec<usize>], divisor_count: usize) -> Result<()> {
    let mut seen = HashSet::new();
    for part in nef_parts {
        if part.is_empty() {
            return Err(Error::InvalidInput(
                "nef partition contains an empty part".into(),
            ));
        }
        for &idx in part {
            if idx >= divisor_count {
                return Err(Error::InvalidInput(format!(
                    "nef partition index {idx} is out of range for divisor count {divisor_count}"
                )));
            }
            if !seen.insert(idx) {
                return Err(Error::InvalidInput(format!(
                    "nef partition index {idx} appears in more than one part"
                )));
            }
        }
    }
    Ok(())
}

fn validate_ambient_tuple(
    indices: &[usize],
    ambient_tuple_len: usize,
    divisor_count: usize,
) -> Result<()> {
    if indices.len() != ambient_tuple_len {
        return Err(Error::InvalidInput(format!(
            "ambient intersection tuple length must be {ambient_tuple_len}, got {}",
            indices.len()
        )));
    }
    for &idx in indices {
        if idx >= divisor_count {
            return Err(Error::InvalidInput(format!(
                "ambient intersection index {idx} is out of range for divisor count {divisor_count}"
            )));
        }
    }
    Ok(())
}

fn remove_one_nef_part_index_choices(
    indices: &[usize],
    part: &HashSet<usize>,
) -> BTreeSet<Vec<usize>> {
    let mut choices = BTreeSet::new();
    for remove_pos in 0..indices.len() {
        if !part.contains(&indices[remove_pos]) {
            continue;
        }
        let mut reduced = indices
            .iter()
            .enumerate()
            .filter_map(|(idx, &value)| (idx != remove_pos).then_some(value))
            .collect::<Vec<_>>();
        reduced.sort_unstable();
        choices.insert(reduced);
    }
    choices
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(value: i64) -> TypedRational<Finite> {
        TypedRational::from_raw(MalachiteRational::from(value))
    }

    #[test]
    fn cicy_reduction_applies_single_nef_part_once_per_distinct_choice() {
        let ambient = vec![
            (vec![0, 0, 0, 3], q(-1)),
            (vec![0, 1, 1, 3], q(2)),
            (vec![0, 1, 3, 3], q(4)),
            (vec![1, 1, 2, 2], q(7)),
        ];

        let reduced =
            compute_complete_intersection_cy3_intersection_numbers(&ambient, &[vec![3]], 4)
                .unwrap();

        assert_eq!(reduced.get(0, 0, 0).get(), &MalachiteRational::from(-1));
        assert_eq!(reduced.get(0, 1, 1).get(), &MalachiteRational::from(2));
        assert_eq!(reduced.get(0, 1, 3).get(), &MalachiteRational::from(4));
        assert_eq!(reduced.get(1, 1, 2).get(), &MalachiteRational::ZERO);
    }

    #[test]
    fn cicy_reduction_threads_multiple_nef_parts_in_order() {
        let ambient = vec![
            (vec![0, 1, 2, 3, 4], q(6)),
            (vec![0, 1, 3, 4, 4], q(2)),
            (vec![0, 2, 2, 3, 3], q(9)),
        ];

        let reduced = compute_complete_intersection_cy3_intersection_numbers(
            &ambient,
            &[vec![3], vec![4]],
            5,
        )
        .unwrap();

        assert_eq!(reduced.get(0, 1, 2).get(), &MalachiteRational::from(6));
        assert_eq!(reduced.get(0, 1, 4).get(), &MalachiteRational::from(2));
        assert_eq!(reduced.get(0, 2, 2).get(), &MalachiteRational::ZERO);
    }

    #[test]
    fn cicy_reduction_rejects_overlapping_nef_parts() {
        let ambient = vec![(vec![0, 1, 2, 3, 4], q(1))];

        let err = compute_complete_intersection_cy3_intersection_numbers(
            &ambient,
            &[vec![3], vec![3]],
            5,
        )
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("nef partition index 3 appears in more than one part")
        );
    }

    #[test]
    fn full_dim_nef_certificate_extracts_quintic_anticanonical_part() {
        let points = vec![
            Point::new(vec![0, 0, 0, 0]),
            Point::new(vec![-1, -1, -1, -1]),
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![0, 0, 0, 1]),
        ];

        let certificate =
            certify_nef_partition_cytools_style_full_dim(&points, &points, &[vec![1, 2, 3, 4, 5]])
                .unwrap();

        assert_eq!(certificate.nef_parts, vec![vec![1, 2, 3, 4, 5]]);
        assert_eq!(certificate.part_lattice_point_counts, vec![6]);
        assert_eq!(certificate.minkowski_sum_vertex_count, 5);
    }

    #[test]
    fn full_dim_nef_certificate_rejects_union_hull_mismatch() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![-1, 0]),
            Point::new(vec![0, -1]),
        ];

        let err = certify_nef_partition_cytools_style_full_dim(&points, &points, &[vec![1, 2, 3]])
            .unwrap_err();

        assert!(
            err.to_string()
                .contains("nef partition union hull does not equal ambient polytope hull")
        );
    }

    #[test]
    fn full_dim_nef_certificate_rejects_lower_dimensional_part_gap() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![-1, 0]),
            Point::new(vec![0, -1]),
        ];

        let err = certify_nef_partition_cytools_style_full_dim(
            &points,
            &points,
            &[vec![1, 3], vec![2, 4]],
        )
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("nef part 0 polytope must be full-dimensional")
        );
    }

    #[test]
    fn cytools_nef_certificate_accepts_lower_dimensional_segment_parts() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![-1, 0]),
            Point::new(vec![0, -1]),
        ];

        let certificate =
            certify_nef_partition_cytools_style(&points, &points, &[vec![1, 3], vec![2, 4]])
                .unwrap();

        assert_eq!(certificate.nef_parts, vec![vec![1, 3], vec![2, 4]]);
        assert_eq!(certificate.part_lattice_point_counts, vec![3, 3]);
        assert_eq!(certificate.minkowski_sum_vertex_count, 4);
    }

    #[test]
    fn cicy_from_ambient_cytools_matches_quintic_triples() {
        let points = vec![
            Point::new(vec![0, 0, 0, 0]),
            Point::new(vec![-1, -1, -1, -1]),
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![0, 0, 0, 1]),
        ];
        let tri = Triangulation::new(vec![
            vec![0, 1, 2, 3, 4],
            vec![0, 1, 2, 3, 5],
            vec![0, 1, 2, 4, 5],
            vec![0, 1, 3, 4, 5],
            vec![0, 2, 3, 4, 5],
        ]);
        let linrels = vec![
            vec![-1, 1, 0, 0, 0],
            vec![-1, 0, 1, 0, 0],
            vec![-1, 0, 0, 1, 0],
            vec![-1, 0, 0, 0, 1],
        ];

        let quintic = compute_complete_intersection_cy3_from_ambient_cytools(
            &tri,
            &points,
            &linrels,
            &[vec![1, 2, 3, 4, 5]],
        )
        .unwrap();

        assert_eq!(quintic.get(1, 1, 1).get(), &MalachiteRational::from(5));
        assert_eq!(quintic.get(1, 2, 3).get(), &MalachiteRational::from(5));
        assert_eq!(quintic.get(2, 4, 5).get(), &MalachiteRational::from(5));
    }
}
