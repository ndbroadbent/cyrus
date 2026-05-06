//! Consume corrected-chamber GV context exports.
//!
//! This is an opt-in diagnostic binary for the JSON written by
//! `mcallister_first_principles --dump-corrected-chamber-gv-context`. It does
//! not read McAllister downstream GV rows. Its job is to validate the
//! CYTools/cygv-shaped context and, when requested, run small explicit
//! semigroup HKTY checks for missing targets whose exact active-generator
//! decomposition is an integer semigroup.

use malachite::Rational as MalachiteRational;
use nalgebra::{DMatrix, RowDVector};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

use cyrus_core::types::rational::Rational;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    Intersection, compute_gv_invariants_with_explicit_semigroup,
    compute_gv_invariants_with_provided_generators, cygv_pair_reduced_seed_generators,
};

#[derive(Debug, Deserialize)]
struct CorrectedChamberGvContext {
    schema_version: u32,
    ambient_rays: usize,
    toric_gv_missing_count: usize,
    remaining_gv_missing_count: usize,
    basis_mori_ray_count: Option<usize>,
    basis_mori_rays_for_missing_degree_bound: Option<i128>,
    basis_mori_rays_for_missing_degree_bounded: Option<Vec<Vec<i64>>>,
    gv_q_matrix_for_missing: Option<Vec<Vec<i64>>>,
    grading_for_missing: Option<Vec<i64>>,
    corrected_kappa_basis_for_missing: Option<Vec<SparseIntersectionEntry>>,
    missing_target_stats: Option<MissingGvTargetStats>,
}

#[derive(Debug, Deserialize)]
struct SparseIntersectionEntry {
    indices: [usize; 3],
    value: String,
}

#[derive(Debug, Deserialize)]
struct MissingGvTargetStats {
    target_count: usize,
    real_cone_decomposition_exact_kind_counts: HashMap<String, usize>,
    sample: Vec<MissingGvTargetSample>,
}

#[derive(Debug, Deserialize)]
struct MissingGvTargetSample {
    degree: i128,
    generators_le_degree: usize,
    is_mori_generator: bool,
    origin_circuit_pattern: Option<String>,
    origin_circuit_witness_count: Option<usize>,
    origin_circuit_first_witness: Option<OriginCircuitWitnessSample>,
    origin_circuit_affine_support: Option<OriginCircuitAffineSupportSample>,
    cms_general_divisor_shape_candidates: Option<Vec<CmsGeneralDivisorShapeCandidate>>,
    cms_general_divisor_intersection_checks: Option<Vec<CmsGeneralDivisorIntersectionCheck>>,
    branch_diagnostic: Option<MissingGvBranchDiagnostic>,
    real_cone_decomposable_by_other_generators: bool,
    real_cone_decomposition_active_generators: Option<usize>,
    real_cone_decomposition_active_generator_basis_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    real_cone_decomposition_exact_coefficients: Option<Vec<String>>,
    real_cone_decomposition_exact_kind: Option<String>,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct MissingGvBranchDiagnostic {
    q_dot_t: String,
    parity: i128,
    parity_mod2: i128,
    q_dot_bucket: String,
    dilog_status: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct OriginCircuitAffineSupportSample {
    affine_rank: usize,
    coefficient_counts: BTreeMap<i64, usize>,
    local_charge_basis: Vec<Vec<i64>>,
    local_coordinates_2d: Option<Vec<OriginCircuitLocalCoordinate2DSample>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct OriginCircuitLocalCoordinate2DSample {
    point_index: usize,
    coordinates: [i64; 2],
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct CmsGeneralDivisorShapeCandidate {
    shrinking_divisor_index: usize,
    shrinking_divisor_coefficient: i64,
    shrinking_divisor_coordinates: Vec<i64>,
    inferred_other_normal_degree: i64,
    toric_gv1_formula_value: Option<i64>,
    all_non_origin_relation_points_are_two_face: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct CmsGeneralDivisorIntersectionCheck {
    shrinking_divisor_index: usize,
    has_rational_divisor_solution: bool,
    solution_basis_support_len: Option<usize>,
    solution_is_integral: Option<bool>,
    computed_other_normal_degree: Option<String>,
    matches_inferred_other_normal_degree: Option<bool>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct OriginCircuitWitnessSample {
    first_facet_exclusive_point: usize,
    second_facet_exclusive_point: usize,
    shared_two_simplex: Vec<usize>,
    first_facet_size: usize,
    second_facet_size: usize,
    sparse_relation: Vec<(usize, i64)>,
    relation_points: Vec<OriginCircuitRelationPointSample>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct OriginCircuitRelationPointSample {
    point_index: usize,
    coefficient: i64,
    coordinates: Vec<i64>,
    face_dimension: Option<usize>,
}

#[derive(Debug, Serialize)]
struct ContextReport {
    schema_version: u32,
    dimension: usize,
    ambient_rays: usize,
    basis_mori_ray_count: Option<usize>,
    degree_bound: i128,
    degree_bounded_ray_count: usize,
    q_rows: usize,
    q_cols: usize,
    kappa_nonzero_entries: usize,
    toric_gv_missing_count: usize,
    remaining_gv_missing_count: usize,
    missing_target_count: usize,
    exact_kind_counts: HashMap<String, usize>,
    targets: Vec<TargetReport>,
}

#[derive(Debug, Serialize)]
struct TargetReport {
    index: usize,
    degree: i128,
    generators_le_degree: usize,
    is_mori_generator: bool,
    origin_circuit_pattern: Option<String>,
    origin_circuit_witness_count: Option<usize>,
    origin_circuit_first_witness: Option<OriginCircuitWitnessSample>,
    origin_circuit_affine_support: Option<OriginCircuitAffineSupportSample>,
    local_cygv_hypersurface_shape: Option<LocalCygvHypersurfaceShape>,
    cms_general_divisor_shape_candidates: Option<Vec<CmsGeneralDivisorShapeCandidate>>,
    cms_general_divisor_intersection_checks: Option<Vec<CmsGeneralDivisorIntersectionCheck>>,
    branch_diagnostic: Option<MissingGvBranchDiagnostic>,
    real_cone_decomposable_by_other_generators: bool,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
    target_cygv_negative_intersections: Option<usize>,
    target_cygv_omega_bucket: Option<String>,
    target_cygv_series_coordinate: Option<usize>,
    target_cygv_series_coordinate_kappa_pair_count: Option<usize>,
    target_cygv_nonzero_coordinate_kappa_pair_counts: Vec<CygvCoordinateKappaSupport>,
    exact_kind: Option<String>,
    active_generator_count: Option<usize>,
    integer_term_count: Option<usize>,
    diamond_element_count: Option<usize>,
    status: String,
    gv: Option<String>,
    error: Option<String>,
    active_support_generator_count: Option<usize>,
    active_support_status: Option<String>,
    active_support_gv: Option<String>,
    active_support_error: Option<String>,
    degree_bounded_candidate_count: usize,
    support_overlap_generator_counts: Vec<SupportOverlapCount>,
    support_closure_layer_counts: Vec<SupportClosureLayerCount>,
    support_overlap_min_for_run: Option<usize>,
    support_overlap_pair_reduce_for_run: bool,
    support_overlap_run_generator_count: Option<usize>,
    support_overlap_run_status: Option<String>,
    support_overlap_run_gv: Option<String>,
    support_overlap_run_error: Option<String>,
    cygv_semigroup_measure_status: Option<String>,
    cygv_semigroup_seed_count: Option<usize>,
    cygv_semigroup_reduced_seed_count: Option<usize>,
    cygv_semigroup_target_is_seed: Option<bool>,
    cygv_semigroup_target_is_reduced_seed: Option<bool>,
    cygv_semigroup_seed_negative_histogram: Option<CygvNegativeIntersectionHistogram>,
    cygv_semigroup_reduced_seed_negative_histogram: Option<CygvNegativeIntersectionHistogram>,
    cygv_semigroup_element_count: Option<usize>,
    cygv_semigroup_max_degree: Option<u32>,
    cygv_semigroup_error: Option<String>,
    cygv_path_history_probe: Option<CygvPathHistoryProbe>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvHypersurfaceShape {
    q_rows: usize,
    q_cols: usize,
    cy_codim: usize,
    ambient_dim: i64,
    cy_dim: i64,
    charge_sums: Vec<i64>,
    is_calabi_yau_charge: bool,
    is_compact_threefold_hypersurface_shape: bool,
}

#[derive(Debug, Serialize)]
struct SupportOverlapCount {
    min_overlap: usize,
    generator_count: usize,
}

#[derive(Debug, Serialize)]
struct SupportClosureLayerCount {
    layer: usize,
    generator_count: usize,
    support_size: usize,
}

#[derive(Clone, Debug, Serialize)]
struct CygvNegativeIntersectionHistogram {
    neg0: usize,
    neg1: usize,
    neg2: usize,
    gt2: usize,
}

#[derive(Clone, Debug, Serialize)]
struct CygvCoordinateKappaSupport {
    coordinate: usize,
    kappa_pair_count: usize,
}

#[derive(Clone, Debug, Serialize)]
struct CygvPathHistoryProbe {
    status: String,
    closure_element_count: Option<usize>,
    closure_degree_counts: BTreeMap<i128, usize>,
    target_in_closure: Option<bool>,
    previous_level_count: usize,
    previous_window_degrees: Vec<i128>,
    previous_window_degree_count: Option<usize>,
    previous_window_element_count: Option<usize>,
    predecessor_difference_count: Option<usize>,
    improving_predecessor_difference_count: Option<usize>,
    closest_series_distance: Option<String>,
    closest_series_predecessor_nonzero: Option<Vec<(usize, i64)>>,
    closest_series_difference_nonzero: Option<Vec<(usize, i64)>>,
}

struct CygvSemigroupMeasurement {
    status: String,
    seed_count: usize,
    reduced_seed_count: usize,
    target_is_seed: bool,
    target_is_reduced_seed: bool,
    seed_negative_histogram: CygvNegativeIntersectionHistogram,
    reduced_seed_negative_histogram: CygvNegativeIntersectionHistogram,
    max_degree: u32,
    element_count: Option<usize>,
}

struct CygvSemigroupDegreeMeasurement {
    status: String,
    seed_count: usize,
    reduced_seed_count: usize,
    seed_set: HashSet<Vec<i64>>,
    reduced_seed_set: HashSet<Vec<i64>>,
    seed_negative_histogram: CygvNegativeIntersectionHistogram,
    reduced_seed_negative_histogram: CygvNegativeIntersectionHistogram,
    max_degree: u32,
    element_count: Option<usize>,
}

fn parse_arg_value<T: std::str::FromStr>(flag: &str) -> Option<T> {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == flag {
            return args.next().and_then(|value| value.parse::<T>().ok());
        }
    }
    None
}

fn parse_flag(flag: &str) -> bool {
    std::env::args().any(|arg| arg == flag)
}

fn panic_payload_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else {
        "non-string panic payload".to_string()
    }
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &PathBuf) -> Result<T, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    serde_json::from_str(&content).map_err(|e| format!("failed to parse {}: {e}", path.display()))
}

fn dense_from_sparse(entries: &[(usize, i64)], dimension: usize) -> Result<Vec<i64>, String> {
    let mut out = vec![0i64; dimension];
    for &(idx, value) in entries {
        let Some(slot) = out.get_mut(idx) else {
            return Err(format!(
                "sparse coordinate {idx} is out of bounds for dimension {dimension}"
            ));
        };
        if *slot != 0 {
            return Err(format!("duplicate sparse coordinate {idx}"));
        }
        *slot = value;
    }
    Ok(out)
}

fn sparse_from_dense(values: &[i64]) -> Vec<(usize, i64)> {
    values
        .iter()
        .enumerate()
        .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
        .collect()
}

fn q_intersections(curve: &[i64], q_matrix: &[Vec<i64>]) -> Result<Vec<i128>, String> {
    if q_matrix.is_empty() {
        return Err("q-matrix is empty".to_string());
    }
    if q_matrix.len() != curve.len() {
        return Err(format!(
            "q-matrix row count {} does not match curve dimension {}",
            q_matrix.len(),
            curve.len()
        ));
    }
    let q_cols = q_matrix[0].len();
    if q_matrix.iter().any(|row| row.len() != q_cols) {
        return Err("q-matrix rows have inconsistent lengths".to_string());
    }
    let mut intersections = vec![0i128; q_cols];
    for (coefficient, row) in curve.iter().zip(q_matrix.iter()) {
        for (slot, &charge) in intersections.iter_mut().zip(row.iter()) {
            *slot = slot
                .checked_add(i128::from(*coefficient) * i128::from(charge))
                .ok_or_else(|| "q-intersection overflowed i128".to_string())?;
        }
    }
    Ok(intersections)
}

fn cygv_negative_intersection_count(curve: &[i64], q_matrix: &[Vec<i64>]) -> Result<usize, String> {
    Ok(q_intersections(curve, q_matrix)?
        .into_iter()
        .filter(|value| *value < 0)
        .count())
}

fn cygv_omega_bucket(negative_intersections: usize) -> String {
    match negative_intersections {
        0 => "neg0".to_string(),
        1 => "neg1".to_string(),
        2 => "neg2".to_string(),
        _ => "ignored_gt2".to_string(),
    }
}

fn cygv_negative_intersection_histogram(
    curves: impl IntoIterator<Item = Vec<i64>>,
    q_matrix: &[Vec<i64>],
) -> Result<CygvNegativeIntersectionHistogram, String> {
    let mut histogram = CygvNegativeIntersectionHistogram {
        neg0: 0,
        neg1: 0,
        neg2: 0,
        gt2: 0,
    };
    for curve in curves {
        match cygv_negative_intersection_count(&curve, q_matrix)? {
            0 => histogram.neg0 += 1,
            1 => histogram.neg1 += 1,
            2 => histogram.neg2 += 1,
            _ => histogram.gt2 += 1,
        }
    }
    Ok(histogram)
}

fn kappa_pair_count_for_series_coordinate(
    intersection: &Intersection,
    coordinate: usize,
) -> Result<usize, String> {
    if coordinate >= intersection.dim() {
        return Err(format!(
            "series coordinate {coordinate} is out of bounds for intersection dimension {}",
            intersection.dim()
        ));
    }
    let mut count = 0usize;
    for lhs in 0..intersection.dim() {
        for rhs in lhs..intersection.dim() {
            if *intersection.get(coordinate, lhs, rhs).get() != 0 {
                count += 1;
            }
        }
    }
    Ok(count)
}

fn cygv_series_coordinate_support(
    curve: &[i64],
    intersection: &Intersection,
) -> Result<
    (
        Option<usize>,
        Option<usize>,
        Vec<CygvCoordinateKappaSupport>,
    ),
    String,
> {
    if curve.len() != intersection.dim() {
        return Err(format!(
            "curve dimension {} does not match intersection dimension {}",
            curve.len(),
            intersection.dim()
        ));
    }
    let mut nonzero_coordinate_counts = Vec::new();
    for (coordinate, value) in curve.iter().enumerate() {
        if *value == 0 {
            continue;
        }
        nonzero_coordinate_counts.push(CygvCoordinateKappaSupport {
            coordinate,
            kappa_pair_count: kappa_pair_count_for_series_coordinate(intersection, coordinate)?,
        });
    }
    let series_coordinate = nonzero_coordinate_counts
        .first()
        .map(|entry| entry.coordinate);
    let series_coordinate_kappa_pair_count = nonzero_coordinate_counts
        .first()
        .map(|entry| entry.kappa_pair_count);
    Ok((
        series_coordinate,
        series_coordinate_kappa_pair_count,
        nonzero_coordinate_counts,
    ))
}

fn parse_rational(value: &str) -> Result<MalachiteRational, String> {
    value
        .parse::<MalachiteRational>()
        .map_err(|_| format!("failed to parse rational coefficient {value}"))
}

fn rational_to_nonnegative_integer(value: &MalachiteRational) -> Result<Option<usize>, String> {
    let text = value.to_string();
    if text.contains('/') {
        return Ok(None);
    }
    let parsed = text
        .parse::<i128>()
        .map_err(|e| format!("failed to parse integer rational {text}: {e}"))?;
    if parsed < 0 {
        return Err(format!(
            "integer semigroup coefficient is negative: {parsed}"
        ));
    }
    usize::try_from(parsed)
        .map(Some)
        .map_err(|_| format!("integer coefficient {parsed} does not fit in usize"))
}

fn integer_decomposition_terms(
    generators: &[Vec<i64>],
    coefficients: &[String],
) -> Result<Option<Vec<Vec<i64>>>, String> {
    if generators.len() != coefficients.len() {
        return Err(format!(
            "active generator count {} does not match coefficient count {}",
            generators.len(),
            coefficients.len()
        ));
    }

    let mut terms = Vec::new();
    for (generator, coefficient) in generators.iter().zip(coefficients.iter()) {
        let rational = parse_rational(coefficient)?;
        let Some(multiplicity) = rational_to_nonnegative_integer(&rational)? else {
            return Ok(None);
        };
        for _ in 0..multiplicity {
            terms.push(generator.clone());
        }
    }
    Ok(Some(terms))
}

fn decomposition_diamond_elements(
    terms: &[Vec<i64>],
    target: &[i64],
) -> Result<Vec<Vec<i64>>, String> {
    let dimension = target.len();
    let zero = vec![0i64; dimension];
    let mut elements = vec![zero.clone()];
    let mut seen = HashSet::from([zero]);
    for term in terms {
        if term.len() != dimension {
            return Err(format!(
                "decomposition term dimension {} does not match target dimension {dimension}",
                term.len()
            ));
        }
        let existing = elements.clone();
        for element in existing {
            let mut sum = Vec::with_capacity(dimension);
            for (&lhs, &rhs) in element.iter().zip(term.iter()) {
                sum.push(lhs.checked_add(rhs).ok_or_else(|| {
                    "decomposition diamond coordinate overflowed i64".to_string()
                })?);
            }
            if seen.insert(sum.clone()) {
                elements.push(sum);
            }
        }
    }
    if !seen.contains(target) {
        return Err("decomposition diamond does not contain target".to_string());
    }
    elements.sort();
    Ok(elements)
}

fn reconstruct_intersection(
    dimension: usize,
    entries: &[SparseIntersectionEntry],
) -> Result<Intersection, String> {
    let mut intersection = Intersection::new(dimension);
    for entry in entries {
        if entry.indices.iter().any(|&idx| idx >= dimension) {
            return Err(format!(
                "intersection index {:?} is out of bounds for dimension {dimension}",
                entry.indices
            ));
        }
        let value = parse_rational(&entry.value)?;
        intersection.set(
            entry.indices[0],
            entry.indices[1],
            entry.indices[2],
            Rational::<Finite>::new(value),
        );
    }
    Ok(intersection)
}

fn target_active_support(
    sample: &MissingGvTargetSample,
    dimension: usize,
) -> Result<HashSet<usize>, String> {
    let target = dense_from_sparse(&sample.basis_nonzero, dimension)?;
    let mut support = target
        .iter()
        .enumerate()
        .filter_map(|(idx, &value)| (value != 0).then_some(idx))
        .collect::<HashSet<_>>();
    if let Some(active_generators) = sample
        .real_cone_decomposition_active_generator_basis_nonzero
        .as_ref()
    {
        for generator in active_generators {
            for &(idx, value) in generator {
                if idx >= dimension {
                    return Err(format!(
                        "active generator coordinate {idx} is out of bounds for dimension {dimension}"
                    ));
                }
                if value != 0 {
                    support.insert(idx);
                }
            }
        }
    }
    Ok(support)
}

fn support_window_stats(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
) -> Result<
    (
        usize,
        Vec<SupportOverlapCount>,
        Vec<SupportClosureLayerCount>,
    ),
    String,
> {
    let seed_support = target_active_support(sample, context.dimension)?;
    let mut eligible_supports = Vec::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > sample.degree {
            continue;
        }
        let support = ray
            .iter()
            .enumerate()
            .filter_map(|(idx, &value)| (value != 0).then_some(idx))
            .collect::<HashSet<_>>();
        eligible_supports.push(support);
    }
    let overlap_counts = (1..=4)
        .map(|min_overlap| SupportOverlapCount {
            min_overlap,
            generator_count: eligible_supports
                .iter()
                .filter(|support| support.intersection(&seed_support).count() >= min_overlap)
                .count(),
        })
        .collect::<Vec<_>>();

    let mut closure_support = seed_support;
    let mut closure_counts = Vec::new();
    for layer in 1..=4 {
        let selected = eligible_supports
            .iter()
            .filter(|support| !support.is_disjoint(&closure_support))
            .collect::<Vec<_>>();
        for support in &selected {
            closure_support.extend(support.iter().copied());
        }
        closure_counts.push(SupportClosureLayerCount {
            layer,
            generator_count: selected.len(),
            support_size: closure_support.len(),
        });
    }
    Ok((eligible_supports.len(), overlap_counts, closure_counts))
}

fn local_cygv_hypersurface_shape(
    sample: &MissingGvTargetSample,
) -> Result<Option<LocalCygvHypersurfaceShape>, String> {
    let Some(support) = sample.origin_circuit_affine_support.as_ref() else {
        return Ok(None);
    };
    let Some(first_row) = support.local_charge_basis.first() else {
        return Err("origin-circuit affine support has no local charge rows".to_string());
    };
    let q_rows = first_row.len();
    if q_rows == 0 {
        return Err("origin-circuit local charge row is empty".to_string());
    }
    if support
        .local_charge_basis
        .iter()
        .any(|row| row.len() != q_rows)
    {
        return Err("origin-circuit local charge rows have inconsistent lengths".to_string());
    }
    let q_cols = support.local_charge_basis.len();
    let q_rows_i64 =
        i64::try_from(q_rows).map_err(|_| "local q row count does not fit in i64".to_string())?;
    let q_cols_i64 = i64::try_from(q_cols)
        .map_err(|_| "local q column count does not fit in i64".to_string())?;
    let cy_codim = 1usize;
    let ambient_dim = q_rows_i64 - q_cols_i64;
    let cy_dim = ambient_dim - i64::try_from(cy_codim).expect("cy_codim fits in i64");
    let charge_sums = support
        .local_charge_basis
        .iter()
        .map(|row| row.iter().sum())
        .collect::<Vec<i64>>();
    let is_calabi_yau_charge = charge_sums.iter().all(|&sum| sum == 0);
    Ok(Some(LocalCygvHypersurfaceShape {
        q_rows,
        q_cols,
        cy_codim,
        ambient_dim,
        cy_dim,
        charge_sums,
        is_calabi_yau_charge,
        is_compact_threefold_hypersurface_shape: is_calabi_yau_charge && cy_dim == 3,
    }))
}

fn validate_context<'a>(
    context: &'a CorrectedChamberGvContext,
) -> Result<ValidatedContext<'a>, String> {
    if context.schema_version != 1 {
        return Err(format!(
            "unsupported corrected-chamber GV context schema {}",
            context.schema_version
        ));
    }
    let grading = context
        .grading_for_missing
        .as_ref()
        .ok_or_else(|| "context is missing grading_for_missing".to_string())?;
    let dimension = grading.len();
    if dimension == 0 {
        return Err("grading_for_missing is empty".to_string());
    }
    let q_matrix = context
        .gv_q_matrix_for_missing
        .as_ref()
        .ok_or_else(|| "context is missing gv_q_matrix_for_missing".to_string())?;
    if q_matrix.len() != dimension {
        return Err(format!(
            "q-matrix row count {} does not match dimension {dimension}",
            q_matrix.len()
        ));
    }
    let Some(q_cols) = q_matrix.first().map(Vec::len) else {
        return Err("q-matrix is empty".to_string());
    };
    if q_matrix.iter().any(|row| row.len() != q_cols) {
        return Err("q-matrix rows have inconsistent lengths".to_string());
    }
    let degree_bounded_rays = context
        .basis_mori_rays_for_missing_degree_bounded
        .as_ref()
        .ok_or_else(|| "context is missing degree-bounded Mori rays".to_string())?;
    for ray in degree_bounded_rays {
        if ray.len() != dimension {
            return Err(format!(
                "degree-bounded Mori ray dimension {} does not match {dimension}",
                ray.len()
            ));
        }
    }
    let degree_bound = context
        .basis_mori_rays_for_missing_degree_bound
        .ok_or_else(|| "context is missing Mori-ray degree bound".to_string())?;
    let kappa_entries = context
        .corrected_kappa_basis_for_missing
        .as_ref()
        .ok_or_else(|| "context is missing corrected kappa entries".to_string())?;
    let intersection = reconstruct_intersection(dimension, kappa_entries)?;
    let stats = context
        .missing_target_stats
        .as_ref()
        .ok_or_else(|| "context is missing missing_target_stats".to_string())?;
    if stats.sample.len() != stats.target_count {
        return Err(format!(
            "missing target sample is incomplete: sample={} target_count={}",
            stats.sample.len(),
            stats.target_count
        ));
    }
    for (idx, sample) in stats.sample.iter().enumerate() {
        if let (Some(declared), Some(active_generators)) = (
            sample.real_cone_decomposition_active_generators,
            sample
                .real_cone_decomposition_active_generator_basis_nonzero
                .as_ref(),
        ) {
            if declared != active_generators.len() {
                return Err(format!(
                    "missing target sample {idx} declares {declared} active generators but contains {} vectors",
                    active_generators.len()
                ));
            }
        }
    }
    Ok(ValidatedContext {
        dimension,
        degree_bound,
        q_cols,
        grading,
        q_matrix,
        degree_bounded_rays,
        intersection,
        stats,
    })
}

struct ValidatedContext<'a> {
    dimension: usize,
    degree_bound: i128,
    q_cols: usize,
    grading: &'a [i64],
    q_matrix: &'a [Vec<i64>],
    degree_bounded_rays: &'a [Vec<i64>],
    intersection: Intersection,
    stats: &'a MissingGvTargetStats,
}

fn report_target(
    index: usize,
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    run_integer_diamonds: bool,
    run_active_support_generators: bool,
    support_overlap_min_for_run: Option<usize>,
    support_overlap_max_target_degree: Option<i128>,
    support_overlap_pair_reduce_for_run: bool,
    measure_cygv_semigroups: bool,
    semigroup_measure_max_target_degree: Option<i128>,
    semigroup_measure_max_seed_count: Option<usize>,
    semigroup_measurement_cache: &mut HashMap<i128, Result<CygvSemigroupDegreeMeasurement, String>>,
    element_limit: usize,
) -> TargetReport {
    let exact_kind = sample.real_cone_decomposition_exact_kind.clone();
    let active_generator_count = sample
        .real_cone_decomposition_active_generator_basis_nonzero
        .as_ref()
        .map(Vec::len)
        .or(sample.real_cone_decomposition_active_generators);
    let local_cygv_hypersurface_shape = match local_cygv_hypersurface_shape(sample) {
        Ok(shape) => shape,
        Err(error) => {
            return TargetReport {
                index,
                degree: sample.degree,
                generators_le_degree: sample.generators_le_degree,
                is_mori_generator: sample.is_mori_generator,
                origin_circuit_pattern: sample.origin_circuit_pattern.clone(),
                origin_circuit_witness_count: sample.origin_circuit_witness_count,
                origin_circuit_first_witness: sample.origin_circuit_first_witness.clone(),
                origin_circuit_affine_support: sample.origin_circuit_affine_support.clone(),
                local_cygv_hypersurface_shape: None,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                branch_diagnostic: sample.branch_diagnostic.clone(),
                real_cone_decomposable_by_other_generators: sample
                    .real_cone_decomposable_by_other_generators,
                ambient_nonzero: sample.ambient_nonzero.clone(),
                basis_nonzero: sample.basis_nonzero.clone(),
                target_cygv_negative_intersections: None,
                target_cygv_omega_bucket: None,
                target_cygv_series_coordinate: None,
                target_cygv_series_coordinate_kappa_pair_count: None,
                target_cygv_nonzero_coordinate_kappa_pair_counts: Vec::new(),
                exact_kind,
                active_generator_count,
                integer_term_count: None,
                diamond_element_count: None,
                status: "error".to_string(),
                gv: None,
                error: Some(error),
                active_support_generator_count: None,
                active_support_status: None,
                active_support_gv: None,
                active_support_error: None,
                degree_bounded_candidate_count: 0,
                support_overlap_generator_counts: Vec::new(),
                support_closure_layer_counts: Vec::new(),
                support_overlap_min_for_run,
                support_overlap_pair_reduce_for_run,
                support_overlap_run_generator_count: None,
                support_overlap_run_status: None,
                support_overlap_run_gv: None,
                support_overlap_run_error: None,
                cygv_semigroup_measure_status: None,
                cygv_semigroup_seed_count: None,
                cygv_semigroup_reduced_seed_count: None,
                cygv_semigroup_target_is_seed: None,
                cygv_semigroup_target_is_reduced_seed: None,
                cygv_semigroup_seed_negative_histogram: None,
                cygv_semigroup_reduced_seed_negative_histogram: None,
                cygv_semigroup_element_count: None,
                cygv_semigroup_max_degree: None,
                cygv_semigroup_error: None,
                cygv_path_history_probe: None,
            };
        }
    };
    let target = match dense_from_sparse(&sample.basis_nonzero, context.dimension) {
        Ok(target) => target,
        Err(error) => {
            return TargetReport {
                index,
                degree: sample.degree,
                generators_le_degree: sample.generators_le_degree,
                is_mori_generator: sample.is_mori_generator,
                origin_circuit_pattern: sample.origin_circuit_pattern.clone(),
                origin_circuit_witness_count: sample.origin_circuit_witness_count,
                origin_circuit_first_witness: sample.origin_circuit_first_witness.clone(),
                origin_circuit_affine_support: sample.origin_circuit_affine_support.clone(),
                local_cygv_hypersurface_shape,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                branch_diagnostic: sample.branch_diagnostic.clone(),
                real_cone_decomposable_by_other_generators: sample
                    .real_cone_decomposable_by_other_generators,
                ambient_nonzero: sample.ambient_nonzero.clone(),
                basis_nonzero: sample.basis_nonzero.clone(),
                target_cygv_negative_intersections: None,
                target_cygv_omega_bucket: None,
                target_cygv_series_coordinate: None,
                target_cygv_series_coordinate_kappa_pair_count: None,
                target_cygv_nonzero_coordinate_kappa_pair_counts: Vec::new(),
                exact_kind,
                active_generator_count,
                integer_term_count: None,
                diamond_element_count: None,
                status: "error".to_string(),
                gv: None,
                error: Some(error),
                active_support_generator_count: None,
                active_support_status: None,
                active_support_gv: None,
                active_support_error: None,
                degree_bounded_candidate_count: 0,
                support_overlap_generator_counts: Vec::new(),
                support_closure_layer_counts: Vec::new(),
                support_overlap_min_for_run,
                support_overlap_pair_reduce_for_run,
                support_overlap_run_generator_count: None,
                support_overlap_run_status: None,
                support_overlap_run_gv: None,
                support_overlap_run_error: None,
                cygv_semigroup_measure_status: None,
                cygv_semigroup_seed_count: None,
                cygv_semigroup_reduced_seed_count: None,
                cygv_semigroup_target_is_seed: None,
                cygv_semigroup_target_is_reduced_seed: None,
                cygv_semigroup_seed_negative_histogram: None,
                cygv_semigroup_reduced_seed_negative_histogram: None,
                cygv_semigroup_element_count: None,
                cygv_semigroup_max_degree: None,
                cygv_semigroup_error: None,
                cygv_path_history_probe: None,
            };
        }
    };
    let negative_intersections = match cygv_negative_intersection_count(&target, context.q_matrix) {
        Ok(count) => count,
        Err(error) => {
            return TargetReport {
                index,
                degree: sample.degree,
                generators_le_degree: sample.generators_le_degree,
                is_mori_generator: sample.is_mori_generator,
                origin_circuit_pattern: sample.origin_circuit_pattern.clone(),
                origin_circuit_witness_count: sample.origin_circuit_witness_count,
                origin_circuit_first_witness: sample.origin_circuit_first_witness.clone(),
                origin_circuit_affine_support: sample.origin_circuit_affine_support.clone(),
                local_cygv_hypersurface_shape,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                branch_diagnostic: sample.branch_diagnostic.clone(),
                real_cone_decomposable_by_other_generators: sample
                    .real_cone_decomposable_by_other_generators,
                ambient_nonzero: sample.ambient_nonzero.clone(),
                basis_nonzero: sample.basis_nonzero.clone(),
                target_cygv_negative_intersections: None,
                target_cygv_omega_bucket: None,
                target_cygv_series_coordinate: None,
                target_cygv_series_coordinate_kappa_pair_count: None,
                target_cygv_nonzero_coordinate_kappa_pair_counts: Vec::new(),
                exact_kind,
                active_generator_count,
                integer_term_count: None,
                diamond_element_count: None,
                status: "error".to_string(),
                gv: None,
                error: Some(error),
                active_support_generator_count: None,
                active_support_status: None,
                active_support_gv: None,
                active_support_error: None,
                degree_bounded_candidate_count: 0,
                support_overlap_generator_counts: Vec::new(),
                support_closure_layer_counts: Vec::new(),
                support_overlap_min_for_run,
                support_overlap_pair_reduce_for_run,
                support_overlap_run_generator_count: None,
                support_overlap_run_status: None,
                support_overlap_run_gv: None,
                support_overlap_run_error: None,
                cygv_semigroup_measure_status: None,
                cygv_semigroup_seed_count: None,
                cygv_semigroup_reduced_seed_count: None,
                cygv_semigroup_target_is_seed: None,
                cygv_semigroup_target_is_reduced_seed: None,
                cygv_semigroup_seed_negative_histogram: None,
                cygv_semigroup_reduced_seed_negative_histogram: None,
                cygv_semigroup_element_count: None,
                cygv_semigroup_max_degree: None,
                cygv_semigroup_error: None,
                cygv_path_history_probe: None,
            };
        }
    };
    let (
        target_cygv_series_coordinate,
        target_cygv_series_coordinate_kappa_pair_count,
        target_cygv_nonzero_coordinate_kappa_pair_counts,
    ) = match cygv_series_coordinate_support(&target, &context.intersection) {
        Ok(support) => support,
        Err(error) => {
            return TargetReport {
                index,
                degree: sample.degree,
                generators_le_degree: sample.generators_le_degree,
                is_mori_generator: sample.is_mori_generator,
                origin_circuit_pattern: sample.origin_circuit_pattern.clone(),
                origin_circuit_witness_count: sample.origin_circuit_witness_count,
                origin_circuit_first_witness: sample.origin_circuit_first_witness.clone(),
                origin_circuit_affine_support: sample.origin_circuit_affine_support.clone(),
                local_cygv_hypersurface_shape,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                branch_diagnostic: sample.branch_diagnostic.clone(),
                real_cone_decomposable_by_other_generators: sample
                    .real_cone_decomposable_by_other_generators,
                ambient_nonzero: sample.ambient_nonzero.clone(),
                basis_nonzero: sample.basis_nonzero.clone(),
                target_cygv_negative_intersections: Some(negative_intersections),
                target_cygv_omega_bucket: Some(cygv_omega_bucket(negative_intersections)),
                target_cygv_series_coordinate: None,
                target_cygv_series_coordinate_kappa_pair_count: None,
                target_cygv_nonzero_coordinate_kappa_pair_counts: Vec::new(),
                exact_kind,
                active_generator_count,
                integer_term_count: None,
                diamond_element_count: None,
                status: "error".to_string(),
                gv: None,
                error: Some(error),
                active_support_generator_count: None,
                active_support_status: None,
                active_support_gv: None,
                active_support_error: None,
                degree_bounded_candidate_count: 0,
                support_overlap_generator_counts: Vec::new(),
                support_closure_layer_counts: Vec::new(),
                support_overlap_min_for_run,
                support_overlap_pair_reduce_for_run,
                support_overlap_run_generator_count: None,
                support_overlap_run_status: None,
                support_overlap_run_gv: None,
                support_overlap_run_error: None,
                cygv_semigroup_measure_status: None,
                cygv_semigroup_seed_count: None,
                cygv_semigroup_reduced_seed_count: None,
                cygv_semigroup_target_is_seed: None,
                cygv_semigroup_target_is_reduced_seed: None,
                cygv_semigroup_seed_negative_histogram: None,
                cygv_semigroup_reduced_seed_negative_histogram: None,
                cygv_semigroup_element_count: None,
                cygv_semigroup_max_degree: None,
                cygv_semigroup_error: None,
                cygv_path_history_probe: None,
            };
        }
    };
    let base = TargetReport {
        index,
        degree: sample.degree,
        generators_le_degree: sample.generators_le_degree,
        is_mori_generator: sample.is_mori_generator,
        origin_circuit_pattern: sample.origin_circuit_pattern.clone(),
        origin_circuit_witness_count: sample.origin_circuit_witness_count,
        origin_circuit_first_witness: sample.origin_circuit_first_witness.clone(),
        origin_circuit_affine_support: sample.origin_circuit_affine_support.clone(),
        local_cygv_hypersurface_shape,
        cms_general_divisor_shape_candidates: sample.cms_general_divisor_shape_candidates.clone(),
        cms_general_divisor_intersection_checks: sample
            .cms_general_divisor_intersection_checks
            .clone(),
        branch_diagnostic: sample.branch_diagnostic.clone(),
        real_cone_decomposable_by_other_generators: sample
            .real_cone_decomposable_by_other_generators,
        ambient_nonzero: sample.ambient_nonzero.clone(),
        basis_nonzero: sample.basis_nonzero.clone(),
        target_cygv_negative_intersections: Some(negative_intersections),
        target_cygv_omega_bucket: Some(cygv_omega_bucket(negative_intersections)),
        target_cygv_series_coordinate,
        target_cygv_series_coordinate_kappa_pair_count,
        target_cygv_nonzero_coordinate_kappa_pair_counts,
        exact_kind,
        active_generator_count,
        integer_term_count: None,
        diamond_element_count: None,
        status: String::new(),
        gv: None,
        error: None,
        active_support_generator_count: None,
        active_support_status: None,
        active_support_gv: None,
        active_support_error: None,
        degree_bounded_candidate_count: 0,
        support_overlap_generator_counts: Vec::new(),
        support_closure_layer_counts: Vec::new(),
        support_overlap_min_for_run,
        support_overlap_pair_reduce_for_run,
        support_overlap_run_generator_count: None,
        support_overlap_run_status: None,
        support_overlap_run_gv: None,
        support_overlap_run_error: None,
        cygv_semigroup_measure_status: None,
        cygv_semigroup_seed_count: None,
        cygv_semigroup_reduced_seed_count: None,
        cygv_semigroup_target_is_seed: None,
        cygv_semigroup_target_is_reduced_seed: None,
        cygv_semigroup_seed_negative_histogram: None,
        cygv_semigroup_reduced_seed_negative_histogram: None,
        cygv_semigroup_element_count: None,
        cygv_semigroup_max_degree: None,
        cygv_semigroup_error: None,
        cygv_path_history_probe: None,
    };
    let (
        degree_bounded_candidate_count,
        support_overlap_generator_counts,
        support_closure_layer_counts,
    ) = match support_window_stats(sample, context) {
        Ok(stats) => stats,
        Err(error) => {
            return TargetReport {
                status: "error".to_string(),
                error: Some(error),
                ..base
            };
        }
    };

    let (
        active_support_generator_count,
        active_support_status,
        active_support_gv,
        active_support_error,
    ) = if run_active_support_generators {
        match active_support_generator_gv(sample, context) {
            Ok((count, status, gv, error)) => (Some(count), Some(status), gv, error),
            Err(error) => (None, Some("error".to_string()), None, Some(error)),
        }
    } else {
        (None, None, None, None)
    };
    let (
        support_overlap_run_generator_count,
        support_overlap_run_status,
        support_overlap_run_gv,
        support_overlap_run_error,
    ) = if let Some(min_overlap) = support_overlap_min_for_run {
        if support_overlap_max_target_degree.is_some_and(|max_degree| sample.degree > max_degree) {
            (
                None,
                Some("skipped_target_degree_limit".to_string()),
                None,
                None,
            )
        } else {
            match support_overlap_generator_gv(
                sample,
                context,
                min_overlap,
                support_overlap_pair_reduce_for_run,
            ) {
                Ok((count, status, gv, error)) => (Some(count), Some(status), gv, error),
                Err(error) => (None, Some("error".to_string()), None, Some(error)),
            }
        }
    } else {
        (None, None, None, None)
    };
    let (
        cygv_semigroup_measure_status,
        cygv_semigroup_seed_count,
        cygv_semigroup_reduced_seed_count,
        cygv_semigroup_target_is_seed,
        cygv_semigroup_target_is_reduced_seed,
        cygv_semigroup_seed_negative_histogram,
        cygv_semigroup_reduced_seed_negative_histogram,
        cygv_semigroup_element_count,
        cygv_semigroup_max_degree,
        cygv_semigroup_error,
    ) = if measure_cygv_semigroups {
        if semigroup_measure_max_target_degree.is_some_and(|max_degree| sample.degree > max_degree)
        {
            (
                Some("skipped_target_degree_limit".to_string()),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
        } else {
            match measure_cygv_semigroup(
                sample,
                context,
                semigroup_measure_max_seed_count,
                semigroup_measurement_cache,
            ) {
                Ok(measurement) => (
                    Some(measurement.status),
                    Some(measurement.seed_count),
                    Some(measurement.reduced_seed_count),
                    Some(measurement.target_is_seed),
                    Some(measurement.target_is_reduced_seed),
                    Some(measurement.seed_negative_histogram),
                    Some(measurement.reduced_seed_negative_histogram),
                    measurement.element_count,
                    Some(measurement.max_degree),
                    None,
                ),
                Err(error) => (
                    Some("error".to_string()),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    Some(error),
                ),
            }
        }
    } else {
        (None, None, None, None, None, None, None, None, None, None)
    };
    let cygv_path_history_probe = if measure_cygv_semigroups {
        if semigroup_measure_max_target_degree.is_some_and(|max_degree| sample.degree > max_degree)
        {
            None
        } else {
            Some(cygv_path_history_probe(
                sample,
                context,
                &target,
                element_limit,
            ))
        }
    } else {
        None
    };

    match report_target_inner(sample, context, run_integer_diamonds, element_limit) {
        Ok((status, term_count, element_count, gv, error)) => TargetReport {
            status,
            integer_term_count: term_count,
            diamond_element_count: element_count,
            gv,
            error,
            active_support_generator_count,
            active_support_status,
            active_support_gv,
            active_support_error,
            support_overlap_run_generator_count,
            support_overlap_run_status,
            support_overlap_run_gv,
            support_overlap_run_error,
            cygv_semigroup_measure_status,
            cygv_semigroup_seed_count,
            cygv_semigroup_reduced_seed_count,
            cygv_semigroup_target_is_seed,
            cygv_semigroup_target_is_reduced_seed,
            cygv_semigroup_seed_negative_histogram,
            cygv_semigroup_reduced_seed_negative_histogram,
            cygv_semigroup_element_count,
            cygv_semigroup_max_degree,
            cygv_semigroup_error,
            cygv_path_history_probe,
            degree_bounded_candidate_count,
            support_overlap_generator_counts,
            support_closure_layer_counts,
            ..base
        },
        Err(error) => TargetReport {
            status: "error".to_string(),
            error: Some(error),
            active_support_generator_count,
            active_support_status,
            active_support_gv,
            active_support_error,
            support_overlap_run_generator_count,
            support_overlap_run_status,
            support_overlap_run_gv,
            support_overlap_run_error,
            cygv_semigroup_measure_status,
            cygv_semigroup_seed_count,
            cygv_semigroup_reduced_seed_count,
            cygv_semigroup_target_is_seed,
            cygv_semigroup_target_is_reduced_seed,
            cygv_semigroup_seed_negative_histogram,
            cygv_semigroup_reduced_seed_negative_histogram,
            cygv_semigroup_element_count,
            cygv_semigroup_max_degree,
            cygv_semigroup_error,
            cygv_path_history_probe,
            degree_bounded_candidate_count,
            support_overlap_generator_counts,
            support_closure_layer_counts,
            ..base
        },
    }
}

fn report_target_inner(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    run_integer_diamonds: bool,
    element_limit: usize,
) -> Result<
    (
        String,
        Option<usize>,
        Option<usize>,
        Option<String>,
        Option<String>,
    ),
    String,
> {
    if sample.real_cone_decomposition_exact_kind.as_deref() != Some("integer_semigroup") {
        return Ok((
            "skipped_non_integer_decomposition".to_string(),
            None,
            None,
            None,
            None,
        ));
    }
    let active_generators = sample
        .real_cone_decomposition_active_generator_basis_nonzero
        .as_ref()
        .ok_or_else(|| "integer target is missing active generator vectors".to_string())?
        .iter()
        .map(|sparse| dense_from_sparse(sparse, context.dimension))
        .collect::<Result<Vec<_>, _>>()?;
    let coefficients = sample
        .real_cone_decomposition_exact_coefficients
        .as_ref()
        .ok_or_else(|| "integer target is missing exact coefficients".to_string())?;
    let Some(terms) = integer_decomposition_terms(&active_generators, coefficients)? else {
        return Ok((
            "skipped_rational_coefficients".to_string(),
            None,
            None,
            None,
            None,
        ));
    };
    let target = dense_from_sparse(&sample.basis_nonzero, context.dimension)?;
    let elements = decomposition_diamond_elements(&terms, &target)?;
    if elements.len() > element_limit {
        return Ok((
            format!("skipped_element_limit_{element_limit}"),
            Some(terms.len()),
            Some(elements.len()),
            None,
            None,
        ));
    }
    if !run_integer_diamonds {
        return Ok((
            "ready_integer_decomposition_diamond".to_string(),
            Some(terms.len()),
            Some(elements.len()),
            None,
            None,
        ));
    }
    if cfg!(panic = "abort") {
        return Err(
            "running cygv explicit-semigroup HKTY requires a panic=unwind build for diagnostics"
                .to_string(),
        );
    }
    let target_i32 = target
        .iter()
        .map(|&value| {
            i32::try_from(value).map_err(|_| "target coordinate does not fit in i32".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_gv_invariants_with_explicit_semigroup(
            &elements,
            context.grading,
            context.q_matrix,
            &context.intersection,
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    let gvs = match gvs_result {
        Ok(Ok(gvs)) => gvs,
        Ok(Err(error)) => {
            return Ok((
                "hkty_error".to_string(),
                Some(terms.len()),
                Some(elements.len()),
                None,
                Some(format!("explicit-semigroup HKTY failed: {error}")),
            ));
        }
        Err(payload) => {
            return Ok((
                "hkty_panic".to_string(),
                Some(terms.len()),
                Some(elements.len()),
                None,
                Some(format!(
                    "explicit-semigroup HKTY panicked: {}",
                    panic_payload_message(payload.as_ref())
                )),
            ));
        }
    };
    let gv = gvs
        .into_iter()
        .find_map(|(curve, value)| (curve == target_i32).then(|| value.to_string()))
        .unwrap_or_else(|| "0".to_string());
    Ok((
        "computed_integer_decomposition_diamond".to_string(),
        Some(terms.len()),
        Some(elements.len()),
        Some(gv),
        None,
    ))
}

fn measure_cygv_semigroup(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    max_seed_count: Option<usize>,
    cache: &mut HashMap<i128, Result<CygvSemigroupDegreeMeasurement, String>>,
) -> Result<CygvSemigroupMeasurement, String> {
    let target = dense_from_sparse(&sample.basis_nonzero, context.dimension)?;
    if !cache.contains_key(&sample.degree) {
        let measurement = measure_cygv_semigroup_degree(sample.degree, context, max_seed_count);
        cache.insert(sample.degree, measurement);
    }
    let degree_measurement = cache
        .get(&sample.degree)
        .expect("measurement was inserted above")
        .as_ref()
        .map_err(Clone::clone)?;
    Ok(CygvSemigroupMeasurement {
        status: degree_measurement.status.clone(),
        seed_count: degree_measurement.seed_count,
        reduced_seed_count: degree_measurement.reduced_seed_count,
        target_is_seed: degree_measurement.seed_set.contains(&target),
        target_is_reduced_seed: degree_measurement.reduced_seed_set.contains(&target),
        seed_negative_histogram: degree_measurement.seed_negative_histogram.clone(),
        reduced_seed_negative_histogram: degree_measurement.reduced_seed_negative_histogram.clone(),
        max_degree: degree_measurement.max_degree,
        element_count: degree_measurement.element_count,
    })
}

fn measure_cygv_semigroup_degree(
    target_degree: i128,
    context: &ValidatedContext<'_>,
    max_seed_count: Option<usize>,
) -> Result<CygvSemigroupDegreeMeasurement, String> {
    let max_deg = u32::try_from(target_degree)
        .map_err(|_| format!("target degree {target_degree} does not fit in u32"))?;
    let mut seeds = Vec::new();
    let mut seen = HashSet::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > target_degree {
            continue;
        }
        if seen.insert(ray.clone()) {
            seeds.push(ray.clone());
        }
    }
    let reduced_seeds = cygv_pair_reduced_seed_generators(&seeds)
        .map_err(|error| format!("cygv seed reduction failed: {error}"))?;
    let reduced_seed_set = reduced_seeds.into_iter().collect::<HashSet<_>>();
    let seed_negative_histogram =
        cygv_negative_intersection_histogram(seeds.iter().cloned(), context.q_matrix)?;
    let reduced_seed_negative_histogram =
        cygv_negative_intersection_histogram(reduced_seed_set.iter().cloned(), context.q_matrix)?;
    if max_seed_count.is_some_and(|limit| seeds.len() > limit) {
        return Ok(CygvSemigroupDegreeMeasurement {
            status: "skipped_seed_limit".to_string(),
            seed_count: seeds.len(),
            reduced_seed_count: reduced_seed_set.len(),
            seed_set: seen,
            reduced_seed_set,
            seed_negative_histogram,
            reduced_seed_negative_histogram,
            max_degree: max_deg,
            element_count: None,
        });
    }
    measure_cygv_semigroup_size(&seeds, context.grading, max_deg).map(|element_count| {
        CygvSemigroupDegreeMeasurement {
            status: "measured_cygv_semigroup".to_string(),
            seed_count: seeds.len(),
            reduced_seed_count: reduced_seed_set.len(),
            seed_set: seen,
            reduced_seed_set,
            seed_negative_histogram,
            reduced_seed_negative_histogram,
            max_degree: max_deg,
            element_count: Some(element_count),
        }
    })
}

fn measure_cygv_semigroup_size(
    seeds: &[Vec<i64>],
    grading_vector: &[i64],
    max_deg: u32,
) -> Result<usize, String> {
    if seeds.is_empty() {
        return Err("cygv semigroup seed set is empty".to_string());
    }
    let dimension = grading_vector.len();
    if dimension == 0 {
        return Err("grading vector is empty".to_string());
    }
    if seeds.iter().any(|row| row.len() != dimension) {
        return Err("cygv semigroup seed dimensions do not match grading".to_string());
    }
    let grading = RowDVector::from_row_slice(
        &grading_vector
            .iter()
            .map(|&value| {
                i32::try_from(value).map_err(|_| "grading entry does not fit in i32".to_string())
            })
            .collect::<Result<Vec<_>, _>>()?,
    );
    let mut data = Vec::with_capacity(dimension * seeds.len());
    for col in 0..seeds.len() {
        for row in 0..dimension {
            data.push(
                i32::try_from(seeds[col][row])
                    .map_err(|_| "semigroup seed entry does not fit in i32".to_string())?,
            );
        }
    }
    let generators = DMatrix::from_column_slice(dimension, seeds.len(), &data);
    let semigroup = cygv::Semigroup::with_max_degree(generators, grading, max_deg)
        .map_err(|error| format!("cygv semigroup construction failed: {error}"))?;
    Ok(semigroup.elements.ncols())
}

struct BoundedCygvClosure {
    status: String,
    elements: HashSet<Vec<i64>>,
    degree_counts: BTreeMap<i128, usize>,
    completed: bool,
}

fn cygv_path_history_probe(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    target: &[i64],
    element_limit: usize,
) -> CygvPathHistoryProbe {
    match cygv_path_history_probe_inner(sample, context, target, element_limit) {
        Ok(probe) => probe,
        Err(error) => CygvPathHistoryProbe {
            status: format!("error: {error}"),
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            target_in_closure: None,
            previous_level_count: cygv_previous_level_count(context.dimension),
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
        },
    }
}

fn cygv_path_history_probe_inner(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    target: &[i64],
    element_limit: usize,
) -> Result<CygvPathHistoryProbe, String> {
    let previous_level_count = cygv_previous_level_count(context.dimension);
    let mut seeds = Vec::new();
    let mut seen = HashSet::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > sample.degree {
            continue;
        }
        if seen.insert(ray.clone()) {
            seeds.push(ray.clone());
        }
    }
    if seeds.is_empty() {
        return Ok(CygvPathHistoryProbe {
            status: "skipped_empty_seed_set".to_string(),
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            target_in_closure: None,
            previous_level_count,
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
        });
    }

    let closure =
        bounded_cygv_semigroup_closure(&seeds, context.grading, sample.degree, element_limit)?;
    let target_in_closure = closure.elements.contains(target);
    if !closure.completed {
        return Ok(CygvPathHistoryProbe {
            status: closure.status,
            closure_element_count: Some(closure.elements.len()),
            closure_degree_counts: closure.degree_counts,
            target_in_closure: Some(target_in_closure),
            previous_level_count,
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
        });
    }

    let mut lower_degrees = closure
        .degree_counts
        .keys()
        .copied()
        .filter(|degree| *degree > 0 && *degree < sample.degree)
        .collect::<Vec<_>>();
    lower_degrees.sort_unstable();
    let mut selected_degree_vec = lower_degrees
        .into_iter()
        .rev()
        .take(previous_level_count)
        .collect::<Vec<_>>();
    selected_degree_vec.sort_unstable();
    let selected_degrees = selected_degree_vec.iter().copied().collect::<HashSet<_>>();

    let mut previous_window_element_count = 0usize;
    let mut predecessor_difference_count = 0usize;
    let mut improving_predecessor_difference_count = 0usize;
    let mut closest_distance = cygv_series_distance(target);
    let mut closest_predecessor = None;
    let mut closest_difference = None;
    let mut sorted_elements = closure.elements.iter().collect::<Vec<_>>();
    sorted_elements.sort();
    for element in sorted_elements {
        let degree = curve_degree(element, context.grading)?;
        if !selected_degrees.contains(&degree) {
            continue;
        }
        previous_window_element_count += 1;
        let difference = checked_vector_difference(target, element)?;
        if !closure.elements.contains(&difference) {
            continue;
        }
        predecessor_difference_count += 1;
        let distance = cygv_series_distance(&difference);
        if distance < closest_distance {
            closest_distance = distance;
            improving_predecessor_difference_count += 1;
            closest_predecessor = Some((*element).clone());
            closest_difference = Some(difference);
        }
    }

    Ok(CygvPathHistoryProbe {
        status: "completed_bounded_closure".to_string(),
        closure_element_count: Some(closure.elements.len()),
        closure_degree_counts: closure.degree_counts,
        target_in_closure: Some(target_in_closure),
        previous_level_count,
        previous_window_degrees: selected_degree_vec,
        previous_window_degree_count: Some(selected_degrees.len()),
        previous_window_element_count: Some(previous_window_element_count),
        predecessor_difference_count: Some(predecessor_difference_count),
        improving_predecessor_difference_count: Some(improving_predecessor_difference_count),
        closest_series_distance: Some(format!("{closest_distance:.6}")),
        closest_series_predecessor_nonzero: closest_predecessor.as_deref().map(sparse_from_dense),
        closest_series_difference_nonzero: closest_difference.as_deref().map(sparse_from_dense),
    })
}

fn cygv_previous_level_count(dimension: usize) -> usize {
    if dimension < 4 {
        2
    } else if dimension < 10 {
        5
    } else {
        10
    }
}

fn bounded_cygv_semigroup_closure(
    seeds: &[Vec<i64>],
    grading_vector: &[i64],
    target_degree: i128,
    element_limit: usize,
) -> Result<BoundedCygvClosure, String> {
    if element_limit == 0 {
        return Err("bounded cygv closure element limit must be positive".to_string());
    }
    let dimension = grading_vector.len();
    if dimension == 0 {
        return Err("grading vector is empty".to_string());
    }
    if seeds.iter().any(|row| row.len() != dimension) {
        return Err("cygv closure seed dimensions do not match grading".to_string());
    }
    let generators = cygv_pair_reduced_seed_generators(seeds)
        .map_err(|error| format!("cygv seed reduction failed: {error}"))?;
    let zero = vec![0i64; dimension];
    let mut elements = HashSet::new();
    let mut starting_elements = HashSet::new();
    elements.insert(zero);
    for seed in seeds {
        elements.insert(seed.clone());
        starting_elements.insert(seed.clone());
    }
    let degree_counts = cygv_closure_degree_counts(&elements, grading_vector)?;
    if elements.len() > element_limit {
        return Ok(BoundedCygvClosure {
            status: format!("exceeded_element_limit_initial_{element_limit}"),
            elements,
            degree_counts,
            completed: false,
        });
    }

    loop {
        let mut new_elements = HashSet::new();
        for generator in &generators {
            for element in &starting_elements {
                let sum = checked_vector_sum(generator, element)?;
                let degree = curve_degree(&sum, grading_vector)?;
                if degree <= target_degree && !elements.contains(&sum) {
                    new_elements.insert(sum);
                }
            }
        }
        if new_elements.is_empty() {
            let degree_counts = cygv_closure_degree_counts(&elements, grading_vector)?;
            return Ok(BoundedCygvClosure {
                status: "completed_bounded_closure".to_string(),
                elements,
                degree_counts,
                completed: true,
            });
        }
        if elements.len() + new_elements.len() > element_limit {
            for element in new_elements {
                if elements.len() >= element_limit {
                    break;
                }
                elements.insert(element);
            }
            let degree_counts = cygv_closure_degree_counts(&elements, grading_vector)?;
            return Ok(BoundedCygvClosure {
                status: format!("exceeded_element_limit_{element_limit}"),
                elements,
                degree_counts,
                completed: false,
            });
        }
        for element in &new_elements {
            elements.insert(element.clone());
        }
        starting_elements = new_elements;
    }
}

fn cygv_closure_degree_counts(
    elements: &HashSet<Vec<i64>>,
    grading_vector: &[i64],
) -> Result<BTreeMap<i128, usize>, String> {
    let mut counts = BTreeMap::new();
    for element in elements {
        *counts
            .entry(curve_degree(element, grading_vector)?)
            .or_insert(0) += 1;
    }
    Ok(counts)
}

fn checked_vector_sum(lhs: &[i64], rhs: &[i64]) -> Result<Vec<i64>, String> {
    if lhs.len() != rhs.len() {
        return Err("vector dimensions do not match for sum".to_string());
    }
    lhs.iter()
        .zip(rhs.iter())
        .map(|(&left, &right)| {
            left.checked_add(right)
                .ok_or_else(|| "vector sum overflowed i64".to_string())
        })
        .collect()
}

fn checked_vector_difference(lhs: &[i64], rhs: &[i64]) -> Result<Vec<i64>, String> {
    if lhs.len() != rhs.len() {
        return Err("vector dimensions do not match for difference".to_string());
    }
    lhs.iter()
        .zip(rhs.iter())
        .map(|(&left, &right)| {
            left.checked_sub(right)
                .ok_or_else(|| "vector difference overflowed i64".to_string())
        })
        .collect()
}

fn cygv_series_distance(curve: &[i64]) -> f64 {
    curve
        .iter()
        .map(|value| {
            if *value == 0 {
                0.0
            } else {
                (*value as f64).abs().log2() + 1.0
            }
        })
        .sum()
}

fn active_support_generator_gv(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
) -> Result<(usize, String, Option<String>, Option<String>), String> {
    if cfg!(panic = "abort") {
        return Ok((
            0,
            "skipped_panic_abort".to_string(),
            None,
            Some("cygv diagnostic requires a panic=unwind build".to_string()),
        ));
    }
    let target = dense_from_sparse(&sample.basis_nonzero, context.dimension)?;
    let support = target_active_support(sample, context.dimension)?;
    if support.is_empty() {
        return Err("active-support generator window is empty".to_string());
    }

    let mut generators = Vec::new();
    let mut seen = HashSet::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > sample.degree {
            continue;
        }
        if ray
            .iter()
            .enumerate()
            .all(|(idx, &value)| value == 0 || support.contains(&idx))
            && seen.insert(ray.clone())
        {
            generators.push(ray.clone());
        }
    }
    if !generators
        .iter()
        .any(|ray| ray.as_slice() == target.as_slice())
    {
        generators.push(target.clone());
    }
    generators.sort();
    generators.dedup();

    let max_deg = u32::try_from(sample.degree)
        .map_err(|_| format!("target degree {} does not fit in u32", sample.degree))?;
    let target_i32 = target
        .iter()
        .map(|&value| {
            i32::try_from(value).map_err(|_| "target coordinate does not fit in i32".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_gv_invariants_with_provided_generators(
            &generators,
            context.grading,
            context.q_matrix,
            &context.intersection,
            None,
            Some(max_deg),
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    let gvs = match gvs_result {
        Ok(Ok(gvs)) => gvs,
        Ok(Err(error)) => {
            return Ok((
                generators.len(),
                "hkty_error".to_string(),
                None,
                Some(format!(
                    "active-support provided-generator HKTY failed: {error}"
                )),
            ));
        }
        Err(payload) => {
            return Ok((
                generators.len(),
                "hkty_panic".to_string(),
                None,
                Some(format!(
                    "active-support provided-generator HKTY panicked: {}",
                    panic_payload_message(payload.as_ref())
                )),
            ));
        }
    };
    let gv = gvs
        .into_iter()
        .find_map(|(curve, value)| (curve == target_i32).then(|| value.to_string()))
        .unwrap_or_else(|| "0".to_string());
    Ok((
        generators.len(),
        "computed_active_support_generators".to_string(),
        Some(gv),
        None,
    ))
}

fn support_overlap_generator_gv(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    min_overlap: usize,
    pair_reduce: bool,
) -> Result<(usize, String, Option<String>, Option<String>), String> {
    if cfg!(panic = "abort") {
        return Ok((
            0,
            "skipped_panic_abort".to_string(),
            None,
            Some("cygv diagnostic requires a panic=unwind build".to_string()),
        ));
    }
    let target = dense_from_sparse(&sample.basis_nonzero, context.dimension)?;
    let support = if min_overlap == 0 {
        None
    } else {
        Some(target_active_support(sample, context.dimension)?)
    };
    let mut generators = Vec::new();
    let mut seen = HashSet::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > sample.degree {
            continue;
        }
        let include = support.as_ref().is_none_or(|support| {
            ray.iter()
                .enumerate()
                .filter(|(idx, value)| **value != 0 && support.contains(idx))
                .count()
                >= min_overlap
        });
        if include && seen.insert(ray.clone()) {
            generators.push(ray.clone());
        }
    }
    if !generators
        .iter()
        .any(|ray| ray.as_slice() == target.as_slice())
    {
        generators.push(target.clone());
    }
    generators.sort();
    generators.dedup();
    if pair_reduce {
        generators = cygv_pair_reduced_seed_generators(&generators).map_err(|error| {
            format!("cygv pair-reducing support-overlap generators failed: {error}")
        })?;
    }
    let label = match (min_overlap, pair_reduce) {
        (0, false) => "degree_bounded_generators",
        (0, true) => "degree_bounded_pair_reduced_generators",
        (_, false) => "support_overlap_generators",
        (_, true) => "support_overlap_pair_reduced_generators",
    };
    run_provided_generator_target_gv(&generators, &target, sample.degree, context, label)
}

fn run_provided_generator_target_gv(
    generators: &[Vec<i64>],
    target: &[i64],
    target_degree: i128,
    context: &ValidatedContext<'_>,
    label: &str,
) -> Result<(usize, String, Option<String>, Option<String>), String> {
    let max_deg = u32::try_from(target_degree)
        .map_err(|_| format!("target degree {target_degree} does not fit in u32"))?;
    let target_i32 = target
        .iter()
        .map(|&value| {
            i32::try_from(value).map_err(|_| "target coordinate does not fit in i32".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_gv_invariants_with_provided_generators(
            generators,
            context.grading,
            context.q_matrix,
            &context.intersection,
            None,
            Some(max_deg),
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    let gvs = match gvs_result {
        Ok(Ok(gvs)) => gvs,
        Ok(Err(error)) => {
            return Ok((
                generators.len(),
                "hkty_error".to_string(),
                None,
                Some(format!("{label} HKTY failed: {error}")),
            ));
        }
        Err(payload) => {
            return Ok((
                generators.len(),
                "hkty_panic".to_string(),
                None,
                Some(format!(
                    "{label} HKTY panicked: {}",
                    panic_payload_message(payload.as_ref())
                )),
            ));
        }
    };
    let gv = gvs
        .into_iter()
        .find_map(|(curve, value)| (curve == target_i32).then(|| value.to_string()))
        .unwrap_or_else(|| "0".to_string());
    Ok((
        generators.len(),
        format!("computed_{label}"),
        Some(gv),
        None,
    ))
}

fn curve_degree(curve: &[i64], grading: &[i64]) -> Result<i128, String> {
    if curve.len() != grading.len() {
        return Err(format!(
            "curve dimension {} does not match grading dimension {}",
            curve.len(),
            grading.len()
        ));
    }
    Ok(curve
        .iter()
        .zip(grading.iter())
        .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
        .sum())
}

fn build_report(
    context: &CorrectedChamberGvContext,
    validated: &ValidatedContext<'_>,
    run_integer_diamonds: bool,
    run_active_support_generators: bool,
    support_overlap_min_for_run: Option<usize>,
    support_overlap_max_target_degree: Option<i128>,
    support_overlap_pair_reduce_for_run: bool,
    measure_cygv_semigroups: bool,
    semigroup_measure_max_target_degree: Option<i128>,
    semigroup_measure_max_seed_count: Option<usize>,
    element_limit: usize,
) -> ContextReport {
    let mut semigroup_measurement_cache = HashMap::new();
    let mut targets = Vec::with_capacity(validated.stats.sample.len());
    for (idx, sample) in validated.stats.sample.iter().enumerate() {
        targets.push(report_target(
            idx,
            sample,
            validated,
            run_integer_diamonds,
            run_active_support_generators,
            support_overlap_min_for_run,
            support_overlap_max_target_degree,
            support_overlap_pair_reduce_for_run,
            measure_cygv_semigroups,
            semigroup_measure_max_target_degree,
            semigroup_measure_max_seed_count,
            &mut semigroup_measurement_cache,
            element_limit,
        ));
    }
    ContextReport {
        schema_version: context.schema_version,
        dimension: validated.dimension,
        ambient_rays: context.ambient_rays,
        basis_mori_ray_count: context.basis_mori_ray_count,
        degree_bound: validated.degree_bound,
        degree_bounded_ray_count: validated.degree_bounded_rays.len(),
        q_rows: validated.q_matrix.len(),
        q_cols: validated.q_cols,
        kappa_nonzero_entries: validated.intersection.num_nonzero(),
        toric_gv_missing_count: context.toric_gv_missing_count,
        remaining_gv_missing_count: context.remaining_gv_missing_count,
        missing_target_count: validated.stats.target_count,
        exact_kind_counts: validated
            .stats
            .real_cone_decomposition_exact_kind_counts
            .clone(),
        targets,
    }
}

fn main() {
    let Some(context_path) = parse_arg_value::<PathBuf>("--context") else {
        eprintln!(
            "[ERROR] usage: mcallister_gv_context --context path [--run-integer-diamonds] [--run-active-support-generators] [--run-support-overlap-generators N] [--pair-reduce-support-overlap-generators] [--support-overlap-max-target-degree N] [--measure-cygv-semigroups] [--semigroup-measure-max-target-degree N] [--semigroup-measure-max-seeds N] [--element-limit N] [--out path]\n       use --run-support-overlap-generators 0 to try all degree-bounded generators up to each target degree"
        );
        std::process::exit(2);
    };
    let run_integer_diamonds = parse_flag("--run-integer-diamonds");
    let run_active_support_generators = parse_flag("--run-active-support-generators");
    let support_overlap_min_for_run = parse_arg_value::<usize>("--run-support-overlap-generators");
    let support_overlap_pair_reduce_for_run =
        parse_flag("--pair-reduce-support-overlap-generators");
    let support_overlap_max_target_degree =
        parse_arg_value::<i128>("--support-overlap-max-target-degree");
    let measure_cygv_semigroups = parse_flag("--measure-cygv-semigroups");
    let semigroup_measure_max_target_degree =
        parse_arg_value::<i128>("--semigroup-measure-max-target-degree");
    let semigroup_measure_max_seed_count =
        parse_arg_value::<usize>("--semigroup-measure-max-seeds");
    let element_limit = parse_arg_value::<usize>("--element-limit").unwrap_or(256);
    let out_path = parse_arg_value::<PathBuf>("--out");

    let context = load_json::<CorrectedChamberGvContext>(&context_path).unwrap_or_else(|e| {
        eprintln!("[ERROR] {e}");
        std::process::exit(2);
    });
    let validated = validate_context(&context).unwrap_or_else(|e| {
        eprintln!("[ERROR] invalid corrected-chamber GV context: {e}");
        std::process::exit(2);
    });
    let report = build_report(
        &context,
        &validated,
        run_integer_diamonds,
        run_active_support_generators,
        support_overlap_min_for_run,
        support_overlap_max_target_degree,
        support_overlap_pair_reduce_for_run,
        measure_cygv_semigroups,
        semigroup_measure_max_target_degree,
        semigroup_measure_max_seed_count,
        element_limit,
    );
    let content = serde_json::to_string_pretty(&report).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to serialize context report: {e}");
        std::process::exit(2);
    });
    if let Some(path) = out_path {
        std::fs::write(&path, content).unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to write {}: {e}", path.display());
            std::process::exit(2);
        });
    } else {
        println!("{content}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_from_sparse_rejects_duplicate_or_out_of_bounds_entries() {
        assert_eq!(
            dense_from_sparse(&[(0, 2), (3, -1)], 4).unwrap(),
            vec![2, 0, 0, -1]
        );
        assert_eq!(sparse_from_dense(&[2, 0, 0, -1]), vec![(0, 2), (3, -1)]);
        assert!(dense_from_sparse(&[(4, 1)], 4).is_err());
        assert!(dense_from_sparse(&[(1, 1), (1, 2)], 4).is_err());
    }

    #[test]
    fn integer_decomposition_terms_expands_integral_coefficients_only() {
        let generators = vec![vec![1, 0], vec![0, 1]];
        let terms = integer_decomposition_terms(&generators, &["2".to_string(), "1".to_string()])
            .unwrap()
            .unwrap();
        assert_eq!(terms, vec![vec![1, 0], vec![1, 0], vec![0, 1]]);
        assert!(
            integer_decomposition_terms(&generators, &["1/2".to_string(), "1".to_string()])
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn decomposition_diamond_contains_all_partial_sums() {
        let terms = vec![vec![1, 0], vec![0, 1], vec![1, 0]];
        let target = vec![2, 1];
        let elements = decomposition_diamond_elements(&terms, &target).unwrap();
        assert_eq!(
            elements,
            vec![
                vec![0, 0],
                vec![0, 1],
                vec![1, 0],
                vec![1, 1],
                vec![2, 0],
                vec![2, 1],
            ]
        );
    }

    #[test]
    fn curve_degree_rejects_dimension_mismatch() {
        assert_eq!(curve_degree(&[2, -1], &[3, 5]).unwrap(), 1);
        assert!(curve_degree(&[1], &[1, 2]).is_err());
    }

    #[test]
    fn q_intersection_bucket_counts_cygv_negative_intersections() {
        let q_matrix = vec![vec![1, -2, 0], vec![-1, 1, 3]];
        assert_eq!(q_intersections(&[2, 1], &q_matrix).unwrap(), vec![1, -3, 3]);
        let negative_count = cygv_negative_intersection_count(&[2, 1], &q_matrix).unwrap();
        assert_eq!(negative_count, 1);
        assert_eq!(cygv_omega_bucket(negative_count), "neg1");
        assert!(q_intersections(&[1], &q_matrix).is_err());
    }

    #[test]
    fn q_intersection_histogram_tracks_ignored_gt2_bucket() {
        let q_matrix = vec![vec![-1, -1, -1], vec![2, 0, 0]];
        let histogram = cygv_negative_intersection_histogram(
            vec![vec![1, 0], vec![0, 1], vec![1, 1]],
            &q_matrix,
        )
        .unwrap();
        assert_eq!(histogram.neg0, 1);
        assert_eq!(histogram.neg1, 0);
        assert_eq!(histogram.neg2, 1);
        assert_eq!(histogram.gt2, 1);
    }

    #[test]
    fn series_coordinate_support_counts_kappa_pairs_for_cygv_coordinate() {
        let mut intersection = Intersection::new(3);
        intersection.set(1, 0, 0, Rational::<Finite>::new(MalachiteRational::from(1)));
        intersection.set(1, 1, 2, Rational::<Finite>::new(MalachiteRational::from(1)));
        intersection.set(2, 2, 2, Rational::<Finite>::new(MalachiteRational::from(1)));

        let (series_coordinate, pair_count, all_counts) =
            cygv_series_coordinate_support(&[0, 1, 2], &intersection).unwrap();
        assert_eq!(series_coordinate, Some(1));
        assert_eq!(pair_count, Some(2));
        assert_eq!(all_counts.len(), 2);
        assert_eq!(all_counts[0].coordinate, 1);
        assert_eq!(all_counts[0].kappa_pair_count, 2);
        assert_eq!(all_counts[1].coordinate, 2);
        assert_eq!(all_counts[1].kappa_pair_count, 2);
    }

    #[test]
    fn cygv_semigroup_size_closes_seed_generators_to_degree() {
        let seeds = vec![vec![1, 0], vec![0, 1]];
        let size = measure_cygv_semigroup_size(&seeds, &[1, 1], 2).unwrap();
        assert_eq!(size, 6);
    }

    #[test]
    fn bounded_cygv_closure_mirrors_degree_limited_seed_closure() {
        let seeds = vec![vec![1, 0], vec![0, 1]];
        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 2, 16).unwrap();
        assert!(closure.completed);
        assert_eq!(closure.elements.len(), 6);
        assert_eq!(closure.degree_counts.get(&0), Some(&1));
        assert_eq!(closure.degree_counts.get(&1), Some(&2));
        assert_eq!(closure.degree_counts.get(&2), Some(&3));
        assert!(closure.elements.contains(&vec![1, 1]));
    }

    #[test]
    fn bounded_cygv_closure_matches_actual_cygv_after_pair_seed_reduction() {
        let seeds = vec![vec![1, 0], vec![0, 1], vec![1, 1], vec![2, 0]];
        let reduced = cygv_pair_reduced_seed_generators(&seeds).unwrap();
        assert_eq!(reduced, vec![vec![0, 1], vec![1, 0]]);

        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 2, 16).unwrap();
        assert!(closure.completed);

        let dim = 2;
        let mut generator_data = Vec::new();
        for seed in &seeds {
            for &entry in seed {
                generator_data.push(i32::try_from(entry).unwrap());
            }
        }
        let generators = DMatrix::from_column_slice(dim, seeds.len(), &generator_data);
        let grading = RowDVector::from_row_slice(&[1, 1]);
        let semigroup = cygv::Semigroup::with_max_degree(generators, grading, 2).unwrap();
        let actual = (0..semigroup.elements.ncols())
            .map(|col| {
                semigroup
                    .elements
                    .column(col)
                    .iter()
                    .map(|&entry| i64::from(entry))
                    .collect::<Vec<_>>()
            })
            .collect::<HashSet<_>>();

        assert_eq!(closure.elements, actual);
    }

    #[test]
    fn path_history_probe_counts_cygv_monomial_map_predecessors() {
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 2,
                generators_le_degree: 2,
                is_mori_generator: false,
                origin_circuit_pattern: None,
                origin_circuit_witness_count: None,
                origin_circuit_first_witness: None,
                origin_circuit_affine_support: None,
                cms_general_divisor_shape_candidates: None,
                cms_general_divisor_intersection_checks: None,
                branch_diagnostic: None,
                real_cone_decomposable_by_other_generators: false,
                real_cone_decomposition_active_generators: None,
                real_cone_decomposition_active_generator_basis_nonzero: None,
                real_cone_decomposition_exact_coefficients: None,
                real_cone_decomposition_exact_kind: None,
                ambient_nonzero: vec![(0, 1), (1, 1)],
                basis_nonzero: vec![(0, 1), (1, 1)],
            }],
        };
        let grading = vec![1, 1];
        let q_matrix = vec![vec![1, 0], vec![0, 1]];
        let degree_bounded_rays = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 2,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            intersection: Intersection::new(2),
            stats: &stats,
        };

        let target = vec![1, 1];
        let probe = cygv_path_history_probe_inner(&stats.sample[0], &context, &target, 16).unwrap();
        assert_eq!(probe.status, "completed_bounded_closure");
        assert_eq!(probe.target_in_closure, Some(true));
        assert_eq!(probe.previous_level_count, 2);
        assert_eq!(probe.previous_window_degrees, vec![1]);
        assert_eq!(probe.previous_window_degree_count, Some(1));
        assert_eq!(probe.previous_window_element_count, Some(2));
        assert_eq!(probe.predecessor_difference_count, Some(2));
        assert_eq!(probe.improving_predecessor_difference_count, Some(1));
        assert_eq!(probe.closest_series_distance.as_deref(), Some("1.000000"));
        assert_eq!(probe.closest_series_predecessor_nonzero, Some(vec![(1, 1)]));
        assert_eq!(probe.closest_series_difference_nonzero, Some(vec![(0, 1)]));
    }

    #[test]
    fn target_active_support_merges_target_and_active_generator_supports() {
        let sample = MissingGvTargetSample {
            degree: 3,
            generators_le_degree: 2,
            is_mori_generator: false,
            origin_circuit_pattern: None,
            origin_circuit_witness_count: None,
            origin_circuit_first_witness: None,
            origin_circuit_affine_support: None,
            cms_general_divisor_shape_candidates: None,
            cms_general_divisor_intersection_checks: None,
            branch_diagnostic: None,
            real_cone_decomposable_by_other_generators: true,
            real_cone_decomposition_active_generators: Some(2),
            real_cone_decomposition_active_generator_basis_nonzero: Some(vec![
                vec![(0, 1), (2, -1)],
                vec![(3, 0), (4, 2)],
            ]),
            real_cone_decomposition_exact_coefficients: None,
            real_cone_decomposition_exact_kind: None,
            ambient_nonzero: Vec::new(),
            basis_nonzero: vec![(1, -2), (2, 1)],
        };
        let support = target_active_support(&sample, 5).unwrap();
        assert_eq!(support, HashSet::from([0, 1, 2, 4]));
        assert!(target_active_support(&sample, 4).is_err());
    }

    #[test]
    fn target_report_preserves_origin_circuit_context() {
        let origin_witness = OriginCircuitWitnessSample {
            first_facet_exclusive_point: 7,
            second_facet_exclusive_point: 11,
            shared_two_simplex: vec![2, 3, 5],
            first_facet_size: 4,
            second_facet_size: 4,
            sparse_relation: vec![(7, 1), (11, 1), (0, -2)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 7,
                    coefficient: 1,
                    coordinates: vec![1, 0, 0, 0],
                    face_dimension: Some(3),
                },
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -2,
                    coordinates: vec![0, 0, 0, 0],
                    face_dimension: None,
                },
            ],
        };
        let affine_support = OriginCircuitAffineSupportSample {
            affine_rank: 3,
            coefficient_counts: BTreeMap::from([(-2, 1), (1, 2)]),
            local_charge_basis: vec![vec![1, -2, 1]],
            local_coordinates_2d: None,
        };
        let branch_diagnostic = MissingGvBranchDiagnostic {
            q_dot_t: "1/10".to_string(),
            parity: 1,
            parity_mod2: 1,
            q_dot_bucket: "positive".to_string(),
            dilog_status: "real_ok".to_string(),
        };
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 1,
                generators_le_degree: 1,
                is_mori_generator: false,
                origin_circuit_pattern: Some("-2:1,1:2".to_string()),
                origin_circuit_witness_count: Some(1),
                origin_circuit_first_witness: Some(origin_witness),
                origin_circuit_affine_support: Some(affine_support),
                cms_general_divisor_shape_candidates: Some(vec![CmsGeneralDivisorShapeCandidate {
                    shrinking_divisor_index: 7,
                    shrinking_divisor_coefficient: 1,
                    shrinking_divisor_coordinates: vec![1, 0],
                    inferred_other_normal_degree: 2,
                    toric_gv1_formula_value: Some(4),
                    all_non_origin_relation_points_are_two_face: true,
                }]),
                cms_general_divisor_intersection_checks: Some(vec![
                    CmsGeneralDivisorIntersectionCheck {
                        shrinking_divisor_index: 7,
                        has_rational_divisor_solution: true,
                        solution_basis_support_len: Some(2),
                        solution_is_integral: Some(false),
                        computed_other_normal_degree: Some("3/2".to_string()),
                        matches_inferred_other_normal_degree: Some(false),
                    },
                ]),
                branch_diagnostic: Some(branch_diagnostic),
                real_cone_decomposable_by_other_generators: false,
                real_cone_decomposition_active_generators: None,
                real_cone_decomposition_active_generator_basis_nonzero: None,
                real_cone_decomposition_exact_coefficients: None,
                real_cone_decomposition_exact_kind: None,
                ambient_nonzero: vec![(4, 1)],
                basis_nonzero: vec![(0, 1)],
            }],
        };
        let grading = vec![1, 1];
        let q_matrix = vec![vec![1, 0], vec![0, 1]];
        let degree_bounded_rays = vec![vec![1, 0]];
        let intersection = Intersection::new(2);
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 1,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            intersection,
            stats: &stats,
        };

        let mut semigroup_measurement_cache = HashMap::new();
        let report = report_target(
            0,
            &stats.sample[0],
            &context,
            false,
            false,
            None,
            None,
            false,
            false,
            None,
            None,
            &mut semigroup_measurement_cache,
            256,
        );

        assert_eq!(report.origin_circuit_pattern.as_deref(), Some("-2:1,1:2"));
        assert_eq!(report.origin_circuit_witness_count, Some(1));
        assert_eq!(
            report
                .origin_circuit_first_witness
                .as_ref()
                .unwrap()
                .sparse_relation,
            vec![(7, 1), (11, 1), (0, -2)]
        );
        assert_eq!(
            report
                .origin_circuit_affine_support
                .as_ref()
                .unwrap()
                .affine_rank,
            3
        );
        let shape = report.local_cygv_hypersurface_shape.as_ref().unwrap();
        assert_eq!(shape.q_rows, 3);
        assert_eq!(shape.q_cols, 1);
        assert_eq!(shape.cy_dim, 1);
        assert_eq!(shape.charge_sums, vec![0]);
        assert!(!shape.is_compact_threefold_hypersurface_shape);
        assert_eq!(
            report
                .cms_general_divisor_intersection_checks
                .as_ref()
                .unwrap()[0]
                .computed_other_normal_degree
                .as_deref(),
            Some("3/2")
        );
        assert_eq!(
            report.branch_diagnostic.as_ref().unwrap().dilog_status,
            "real_ok"
        );
    }

    #[test]
    fn support_overlap_run_respects_target_degree_limit() {
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 5,
                generators_le_degree: 1,
                is_mori_generator: true,
                origin_circuit_pattern: None,
                origin_circuit_witness_count: None,
                origin_circuit_first_witness: None,
                origin_circuit_affine_support: None,
                cms_general_divisor_shape_candidates: None,
                cms_general_divisor_intersection_checks: None,
                branch_diagnostic: None,
                real_cone_decomposable_by_other_generators: false,
                real_cone_decomposition_active_generators: None,
                real_cone_decomposition_active_generator_basis_nonzero: None,
                real_cone_decomposition_exact_coefficients: None,
                real_cone_decomposition_exact_kind: None,
                ambient_nonzero: vec![(0, 1)],
                basis_nonzero: vec![(0, 1)],
            }],
        };
        let grading = vec![1];
        let q_matrix = vec![vec![1]];
        let degree_bounded_rays = vec![vec![1]];
        let intersection = Intersection::new(1);
        let context = ValidatedContext {
            dimension: 1,
            degree_bound: 5,
            q_cols: 1,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            intersection,
            stats: &stats,
        };
        let mut semigroup_measurement_cache = HashMap::new();

        let report = report_target(
            0,
            &stats.sample[0],
            &context,
            false,
            false,
            Some(0),
            Some(4),
            false,
            false,
            None,
            None,
            &mut semigroup_measurement_cache,
            256,
        );

        assert_eq!(
            report.support_overlap_run_status.as_deref(),
            Some("skipped_target_degree_limit")
        );
        assert_eq!(report.support_overlap_run_generator_count, None);
    }
}
