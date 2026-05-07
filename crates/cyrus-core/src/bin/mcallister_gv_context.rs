//! Consume corrected-chamber GV context exports.
//!
//! This is an opt-in diagnostic binary for the JSON written by
//! `mcallister_first_principles --dump-corrected-chamber-gv-context`. It does
//! not read McAllister downstream GV rows. Its job is to validate the
//! CYTools/cygv-shaped context and, when requested, run small explicit
//! semigroup HKTY checks for missing targets whose exact active-generator
//! decomposition is an integer semigroup.

use malachite::{Integer, Rational as MalachiteRational};
use nalgebra::{DMatrix, RowDVector};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

use cyrus_core::gv::{
    SupportingMoriFaceLpSearchOptions, certify_supporting_mori_face_by_exact_kernel,
    certify_supporting_mori_face_by_lp_search, cygv_pair_reduced_seed_generators,
    find_extremal_mori_ray_separator,
};
use cyrus_core::types::rational::Rational;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    Intersection, Point, compute_gv_invariants_with_explicit_semigroup,
    compute_gv_invariants_with_provided_generators, curve_row_span_rank,
    diagnose_affine_toric_circuit, integer_math::solve_linear_system_rational, utils::gcd_list_int,
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
    degree_bounded_mori_ray_context_for_missing: Option<Vec<DegreeBoundedMoriRayContextSample>>,
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
struct DegreeBoundedMoriRayContextSample {
    degree: i128,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
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
    #[serde(default)]
    local_coordinates: Vec<OriginCircuitLocalCoordinateSample>,
    local_coordinates_2d: Option<Vec<OriginCircuitLocalCoordinate2DSample>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct OriginCircuitLocalCoordinateSample {
    point_index: usize,
    coordinates: Vec<i64>,
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
    #[serde(default)]
    first_facet: Vec<usize>,
    #[serde(default)]
    second_facet: Vec<usize>,
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
    degree_bounded_mori_ray_context_count: Option<usize>,
    degree_bounded_mori_ray_context_status: String,
    q_rows: usize,
    q_cols: usize,
    kappa_nonzero_entries: usize,
    toric_gv_missing_count: usize,
    remaining_gv_missing_count: usize,
    missing_target_count: usize,
    exact_kind_counts: HashMap<String, usize>,
    local_cygv_charge_signature_counts: BTreeMap<String, usize>,
    local_cygv_target_candidate_status_counts: BTreeMap<String, usize>,
    local_cygv_actual_call_readiness_counts: BTreeMap<String, usize>,
    local_cygv_missing_source_input_counts: BTreeMap<String, usize>,
    origin_circuit_facet_context_status_counts: BTreeMap<String, usize>,
    active_support_status_counts: BTreeMap<String, usize>,
    active_support_face_certificate_status_counts: BTreeMap<String, usize>,
    target_extremal_ray_certificate_status_counts: BTreeMap<String, usize>,
    origin_relation_support_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_shared_facet_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_facet_union_face_certificate_status_counts: BTreeMap<String, usize>,
    local_cygv_q_matrix_orientation_status_counts: BTreeMap<String, usize>,
    local_cygv_q_matrix_layout_status_counts: BTreeMap<String, usize>,
    local_cygv_origin_point_status_counts: BTreeMap<String, usize>,
    local_cytools_origin_circuit_status_counts: BTreeMap<String, usize>,
    local_cygv_grading_vector_status_counts: BTreeMap<String, usize>,
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
    local_cygv_input_skeleton: Option<LocalCygvInputSkeleton>,
    cms_general_divisor_shape_candidates: Option<Vec<CmsGeneralDivisorShapeCandidate>>,
    cms_general_divisor_intersection_checks: Option<Vec<CmsGeneralDivisorIntersectionCheck>>,
    branch_diagnostic: Option<MissingGvBranchDiagnostic>,
    real_cone_decomposable_by_other_generators: bool,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
    target_extremal_ray_certificate: Option<TargetExtremalRayCertificateProbe>,
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
    origin_relation_support_generator_count: Option<usize>,
    origin_shared_facet_generator_count: Option<usize>,
    origin_facet_union_generator_count: Option<usize>,
    origin_relation_support_face_certificate_status: Option<String>,
    origin_shared_facet_face_certificate_status: Option<String>,
    origin_facet_union_face_certificate_status: Option<String>,
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
    cygv_semigroup_degree_ladder: Option<Vec<CygvSemigroupDegreeLadderStep>>,
    cygv_path_history_probe: Option<CygvPathHistoryProbe>,
}

#[derive(Clone, Debug, Serialize)]
struct TargetExtremalRayCertificateProbe {
    status: String,
    same_ray_generator_count: Option<usize>,
    zero_other_generator_count: Option<usize>,
    positive_other_generator_count: Option<usize>,
    separator_normal_nonzero: Option<Vec<(usize, i64)>>,
    decomposition_kind: Option<String>,
    decomposition_active_generator_count: Option<usize>,
    decomposition_exact_coefficients: Option<Vec<String>>,
    decomposition_active_generators_nonzero: Option<Vec<Vec<(usize, i64)>>>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvHypersurfaceShape {
    q_rows: usize,
    q_cols: usize,
    cy_codim: usize,
    ambient_dim: i64,
    cy_dim: i64,
    charge_sums: Vec<i64>,
    charge_row_permutation_signatures: Vec<Vec<i64>>,
    charge_row_multiplicities: Vec<Vec<LocalChargeMultiplicity>>,
    is_calabi_yau_charge: bool,
    is_compact_threefold_hypersurface_shape: bool,
    cygv_compact_input_status: String,
    cygv_compact_input_missing: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalChargeMultiplicity {
    charge: i64,
    count: usize,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvInputSkeleton {
    support_point_indices: Vec<usize>,
    support_contains_origin_point: bool,
    local_cygv_origin_point_status: String,
    origin_point_relation_coefficient: Option<i64>,
    local_cytools_origin_circuit_status: String,
    local_q_matrix_rows: Vec<Vec<i64>>,
    target_relation_coefficients: Option<Vec<i64>>,
    target_relation_in_charge_basis: Option<Vec<i64>>,
    target_relation_status: String,
    orientation_candidates: Vec<LocalCygvOrientationCandidate>,
    local_q_matrix_orientation_candidate: Option<i64>,
    local_q_matrix_orientation_status: String,
    local_cygv_q_matrix_rows_candidate: Option<Vec<Vec<i64>>>,
    local_cygv_wrapper_q_matrix_candidate: Option<Vec<Vec<i64>>>,
    local_cygv_q_matrix_layout_status: String,
    local_grading_vector_candidate: Option<Vec<i64>>,
    local_grading_vector_status: String,
    remaining_uncertified_inputs: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvOrientationCandidate {
    overall_charge_basis_sign: i64,
    local_q_matrix_rows: Vec<Vec<i64>>,
    target_coordinate: Option<Vec<i64>>,
    target_candidate_status: String,
    target_coordinate_is_nonnegative: Option<bool>,
    target_coordinate_gcd: Option<i64>,
    target_coordinate_is_primitive: Option<bool>,
    target_primitive_direction: Option<Vec<i64>>,
    positive_unit_generator_negative_intersections: Option<usize>,
    positive_unit_generator_omega_bucket: Option<String>,
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
    seed_count: Option<usize>,
    reduced_seed_count: Option<usize>,
    closure_element_count: Option<usize>,
    closure_degree_counts: BTreeMap<i128, usize>,
    target_in_closure: Option<bool>,
    previous_level_count: usize,
    previous_window_degrees: Vec<i128>,
    previous_window_degree_count: Option<usize>,
    previous_window_element_count: Option<usize>,
    predecessor_counts_complete: bool,
    predecessor_difference_count: Option<usize>,
    improving_predecessor_difference_count: Option<usize>,
    closest_series_distance: Option<String>,
    closest_series_predecessor_nonzero: Option<Vec<(usize, i64)>>,
    closest_series_difference_nonzero: Option<Vec<(usize, i64)>>,
    lower_seed_decomposition_max_terms: usize,
    lower_seed_decomposition_status: String,
    lower_seed_decomposition_term_count: Option<usize>,
    lower_seed_decomposition_terms_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    lower_seed_decomposition_error: Option<String>,
    lower_seed_diamond_status: Option<String>,
    lower_seed_diamond_element_count: Option<usize>,
    lower_seed_diamond_gv: Option<String>,
    lower_seed_diamond_error: Option<String>,
}

struct CygvPathPredecessorStats {
    previous_window_element_count: usize,
    predecessor_difference_count: usize,
    improving_predecessor_difference_count: usize,
    closest_distance: f64,
    closest_predecessor: Option<Vec<i64>>,
    closest_difference: Option<Vec<i64>>,
}

struct LowerSeedDecompositionProbe {
    status: String,
    term_count: Option<usize>,
    terms_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    terms: Option<Vec<Vec<i64>>>,
    error: Option<String>,
}

struct LowerSeedDiamondProbe {
    status: Option<String>,
    element_count: Option<usize>,
    gv: Option<String>,
    error: Option<String>,
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

#[derive(Clone, Debug, Serialize)]
struct CygvSemigroupDegreeLadderStep {
    degree: i128,
    effective_seed_count: usize,
    reduced_seed_count: Option<usize>,
    status: String,
    element_count: Option<usize>,
    elapsed_ms: Option<u128>,
    error: Option<String>,
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

fn lower_seed_decomposition_probe(
    target: &[i64],
    seeds: &[Vec<i64>],
    max_terms: usize,
) -> LowerSeedDecompositionProbe {
    match bounded_seed_decomposition(target, seeds, max_terms) {
        Ok(Some(terms)) => LowerSeedDecompositionProbe {
            status: "found_lower_seed_decomposition".to_string(),
            term_count: Some(terms.len()),
            terms_nonzero: Some(
                terms
                    .iter()
                    .map(|term| sparse_from_dense(term))
                    .collect::<Vec<_>>(),
            ),
            terms: Some(terms),
            error: None,
        },
        Ok(None) => LowerSeedDecompositionProbe {
            status: format!("not_found_up_to_{max_terms}"),
            term_count: None,
            terms_nonzero: None,
            terms: None,
            error: None,
        },
        Err(error) => LowerSeedDecompositionProbe {
            status: "error".to_string(),
            term_count: None,
            terms_nonzero: None,
            terms: None,
            error: Some(error),
        },
    }
}

fn lower_seed_diamond_probe(
    target: &[i64],
    lower_seed_decomposition: &LowerSeedDecompositionProbe,
    context: &ValidatedContext<'_>,
    run_lower_seed_diamonds: bool,
    element_limit: usize,
) -> LowerSeedDiamondProbe {
    if !run_lower_seed_diamonds {
        return LowerSeedDiamondProbe {
            status: None,
            element_count: None,
            gv: None,
            error: None,
        };
    }
    let Some(terms) = lower_seed_decomposition.terms.as_deref() else {
        return LowerSeedDiamondProbe {
            status: Some("skipped_no_lower_seed_decomposition".to_string()),
            element_count: None,
            gv: None,
            error: lower_seed_decomposition.error.clone(),
        };
    };
    match lower_seed_diamond_probe_inner(target, terms, context, element_limit) {
        Ok(probe) => probe,
        Err(error) => LowerSeedDiamondProbe {
            status: Some("error".to_string()),
            element_count: None,
            gv: None,
            error: Some(error),
        },
    }
}

fn lower_seed_diamond_probe_inner(
    target: &[i64],
    terms: &[Vec<i64>],
    context: &ValidatedContext<'_>,
    element_limit: usize,
) -> Result<LowerSeedDiamondProbe, String> {
    let elements = decomposition_diamond_elements(terms, target)?;
    if elements.len() > element_limit {
        return Ok(LowerSeedDiamondProbe {
            status: Some(format!("skipped_element_limit_{element_limit}")),
            element_count: Some(elements.len()),
            gv: None,
            error: None,
        });
    }
    if cfg!(panic = "abort") {
        return Ok(LowerSeedDiamondProbe {
            status: Some("hkty_unavailable_panic_abort".to_string()),
            element_count: Some(elements.len()),
            gv: None,
            error: Some(
                "running cygv explicit-semigroup HKTY requires a panic=unwind build for diagnostics"
                    .to_string(),
            ),
        });
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
            return Ok(LowerSeedDiamondProbe {
                status: Some("hkty_error".to_string()),
                element_count: Some(elements.len()),
                gv: None,
                error: Some(format!("explicit-semigroup cygv HKTY failed: {error}")),
            });
        }
        Err(payload) => {
            return Ok(LowerSeedDiamondProbe {
                status: Some("hkty_panic".to_string()),
                element_count: Some(elements.len()),
                gv: None,
                error: Some(format!(
                    "explicit-semigroup cygv HKTY panicked: {}",
                    panic_payload_message(payload.as_ref())
                )),
            });
        }
    };
    let gv = gvs
        .into_iter()
        .find_map(|(curve, value)| (curve == target_i32).then(|| value.to_string()))
        .unwrap_or_else(|| "0".to_string());
    Ok(LowerSeedDiamondProbe {
        status: Some("computed_lower_seed_decomposition_diamond".to_string()),
        element_count: Some(elements.len()),
        gv: Some(gv),
        error: None,
    })
}

fn bounded_seed_decomposition(
    target: &[i64],
    seeds: &[Vec<i64>],
    max_terms: usize,
) -> Result<Option<Vec<Vec<i64>>>, String> {
    if max_terms < 2 || seeds.is_empty() {
        return Ok(None);
    }
    if max_terms > 4 {
        return Err("bounded seed decomposition currently supports max_terms <= 4".to_string());
    }
    if seeds.iter().any(|seed| seed.len() != target.len()) {
        return Err("bounded seed decomposition seed dimension mismatch".to_string());
    }

    let pair_sums = seed_pair_sums(seeds)?;
    if let Some(pairs) = pair_sums.get(target)
        && let Some(&(i, j)) = pairs.first()
    {
        return Ok(Some(sorted_decomposition(vec![
            seeds[i].clone(),
            seeds[j].clone(),
        ])));
    }
    if max_terms < 3 {
        return Ok(None);
    }

    for seed in seeds {
        let remainder = checked_vector_difference(target, seed)?;
        let Some(pairs) = pair_sums.get(&remainder) else {
            continue;
        };
        if let Some(&(i, j)) = pairs.first() {
            return Ok(Some(sorted_decomposition(vec![
                seed.clone(),
                seeds[i].clone(),
                seeds[j].clone(),
            ])));
        }
    }
    if max_terms < 4 {
        return Ok(None);
    }

    let mut pair_items = pair_sums.iter().collect::<Vec<_>>();
    pair_items.sort_by(|(left, _), (right, _)| left.cmp(right));
    for (first_sum, first_pairs) in pair_items {
        let remainder = checked_vector_difference(target, first_sum)?;
        let Some(second_pairs) = pair_sums.get(&remainder) else {
            continue;
        };
        let Some(&(i, j)) = first_pairs.first() else {
            continue;
        };
        let Some(&(k, l)) = second_pairs.first() else {
            continue;
        };
        return Ok(Some(sorted_decomposition(vec![
            seeds[i].clone(),
            seeds[j].clone(),
            seeds[k].clone(),
            seeds[l].clone(),
        ])));
    }

    Ok(None)
}

fn seed_pair_sums(seeds: &[Vec<i64>]) -> Result<HashMap<Vec<i64>, Vec<(usize, usize)>>, String> {
    let mut pair_sums: HashMap<Vec<i64>, Vec<(usize, usize)>> = HashMap::new();
    for i in 0..seeds.len() {
        for j in i..seeds.len() {
            let sum = checked_vector_sum(&seeds[i], &seeds[j])?;
            pair_sums.entry(sum).or_default().push((i, j));
        }
    }
    Ok(pair_sums)
}

fn sorted_decomposition(mut terms: Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    terms.sort();
    terms
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

struct OriginCircuitAmbientGeneratorCounts {
    relation_support: Option<usize>,
    shared_facet: Option<usize>,
    facet_union: Option<usize>,
}

struct OriginCircuitAmbientSupportCertificateStatuses {
    relation_support: Option<String>,
    shared_facet: Option<String>,
    facet_union: Option<String>,
}

struct OriginCircuitAmbientSupportSets {
    relation_support: HashSet<usize>,
    shared_facet: HashSet<usize>,
    facet_union: HashSet<usize>,
}

fn origin_circuit_ambient_generator_counts(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
) -> OriginCircuitAmbientGeneratorCounts {
    let Some(ray_context) = context.degree_bounded_ray_context else {
        return OriginCircuitAmbientGeneratorCounts {
            relation_support: None,
            shared_facet: None,
            facet_union: None,
        };
    };
    let Some(supports) = origin_circuit_ambient_support_sets(sample) else {
        return OriginCircuitAmbientGeneratorCounts {
            relation_support: None,
            shared_facet: None,
            facet_union: None,
        };
    };

    OriginCircuitAmbientGeneratorCounts {
        relation_support: Some(degree_bounded_ray_context_support_count(
            ray_context,
            sample.degree,
            &supports.relation_support,
        )),
        shared_facet: Some(degree_bounded_ray_context_support_count(
            ray_context,
            sample.degree,
            &supports.shared_facet,
        )),
        facet_union: Some(degree_bounded_ray_context_support_count(
            ray_context,
            sample.degree,
            &supports.facet_union,
        )),
    }
}

fn origin_circuit_ambient_support_sets(
    sample: &MissingGvTargetSample,
) -> Option<OriginCircuitAmbientSupportSets> {
    let witness = sample.origin_circuit_first_witness.as_ref()?;

    let relation_support = witness
        .relation_points
        .iter()
        .map(|point| point.point_index)
        .collect::<HashSet<_>>();
    let first_facet = witness.first_facet.iter().copied().collect::<HashSet<_>>();
    let second_facet = witness.second_facet.iter().copied().collect::<HashSet<_>>();
    let mut shared_facet = first_facet
        .intersection(&second_facet)
        .copied()
        .collect::<HashSet<_>>();
    shared_facet.insert(0);
    shared_facet.insert(witness.first_facet_exclusive_point);
    shared_facet.insert(witness.second_facet_exclusive_point);
    let mut facet_union = first_facet
        .union(&second_facet)
        .copied()
        .collect::<HashSet<_>>();
    facet_union.insert(0);

    Some(OriginCircuitAmbientSupportSets {
        relation_support,
        shared_facet,
        facet_union,
    })
}

fn degree_bounded_ray_context_support_count(
    ray_context: &[DegreeBoundedMoriRayContextSample],
    max_degree: i128,
    allowed_ambient_support: &HashSet<usize>,
) -> usize {
    ray_context
        .iter()
        .filter(|ray| ray.degree > 0 && ray.degree <= max_degree)
        .filter(|ray| {
            ray.ambient_nonzero
                .iter()
                .all(|(idx, _)| allowed_ambient_support.contains(idx))
        })
        .count()
}

fn origin_circuit_ambient_support_face_certificate_statuses(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    generator_limit: usize,
) -> OriginCircuitAmbientSupportCertificateStatuses {
    let Some(ray_context) = context.degree_bounded_ray_context else {
        return OriginCircuitAmbientSupportCertificateStatuses {
            relation_support: Some("missing_degree_bounded_mori_ray_context".to_string()),
            shared_facet: Some("missing_degree_bounded_mori_ray_context".to_string()),
            facet_union: Some("missing_degree_bounded_mori_ray_context".to_string()),
        };
    };
    let Some(supports) = origin_circuit_ambient_support_sets(sample) else {
        return OriginCircuitAmbientSupportCertificateStatuses {
            relation_support: Some("no_origin_circuit_witness".to_string()),
            shared_facet: Some("no_origin_circuit_witness".to_string()),
            facet_union: Some("no_origin_circuit_witness".to_string()),
        };
    };

    OriginCircuitAmbientSupportCertificateStatuses {
        relation_support: Some(origin_circuit_ambient_support_face_certificate_status(
            ray_context,
            sample.degree,
            &supports.relation_support,
            context,
            generator_limit,
        )),
        shared_facet: Some(origin_circuit_ambient_support_face_certificate_status(
            ray_context,
            sample.degree,
            &supports.shared_facet,
            context,
            generator_limit,
        )),
        facet_union: Some(origin_circuit_ambient_support_face_certificate_status(
            ray_context,
            sample.degree,
            &supports.facet_union,
            context,
            generator_limit,
        )),
    }
}

fn origin_circuit_ambient_support_face_certificate_status(
    ray_context: &[DegreeBoundedMoriRayContextSample],
    max_degree: i128,
    allowed_ambient_support: &HashSet<usize>,
    context: &ValidatedContext<'_>,
    generator_limit: usize,
) -> String {
    let generators = degree_bounded_ray_context_support_generators(
        ray_context,
        max_degree,
        allowed_ambient_support,
        context.dimension,
    );
    let generators = match generators {
        Ok(generators) => generators,
        Err(error) => {
            return format!(
                "origin_support_certificate_error_{}",
                status_error_fragment(&error)
            );
        }
    };
    if generators.is_empty() {
        return "origin_support_no_generators".to_string();
    }
    if generators.len() > generator_limit {
        return format!(
            "origin_support_skipped_generator_limit_{generator_limit}_actual_{}",
            generators.len()
        );
    }
    let rank = match curve_row_span_rank(&generators) {
        Ok(rank) => rank,
        Err(error) => {
            return format!(
                "origin_support_certificate_error_{}",
                status_error_fragment(&error.to_string())
            );
        }
    };
    if rank != context.dimension.saturating_sub(1) {
        let options = SupportingMoriFaceLpSearchOptions::default();
        return match certify_supporting_mori_face_by_lp_search(
            &generators,
            context.degree_bounded_rays,
            &options,
        ) {
            Ok(Some(certificate)) => format!(
                "origin_support_certified_lp_containing_face_rank_{rank}_dim_{}_zero_{}_positive_{}",
                context.dimension,
                certificate.zero_generator_count,
                certificate.positive_generator_count
            ),
            Ok(None) => format!(
                "origin_support_lp_no_certificate_rank_{rank}_dim_{}",
                context.dimension
            ),
            Err(error) => format!(
                "origin_support_certificate_error_{}",
                status_error_fragment(&error.to_string())
            ),
        };
    }
    match certify_supporting_mori_face_by_exact_kernel(&generators, context.degree_bounded_rays) {
        Ok(Some(_)) => "origin_support_certified_codimension_one_face".to_string(),
        Ok(None) => "origin_support_codimension_one_but_not_supporting".to_string(),
        Err(error) => format!(
            "origin_support_certificate_error_{}",
            status_error_fragment(&error.to_string())
        ),
    }
}

fn status_error_fragment(error: &str) -> String {
    let mut out = String::new();
    for ch in error.chars() {
        let mapped = if ch.is_ascii_alphanumeric() {
            ch.to_ascii_lowercase()
        } else {
            '_'
        };
        if out.ends_with('_') && mapped == '_' {
            continue;
        }
        out.push(mapped);
        if out.len() >= 160 {
            break;
        }
    }
    out.trim_matches('_').to_string()
}

fn target_extremal_ray_certificate_probe(
    sample: &MissingGvTargetSample,
    target: &[i64],
    context: &ValidatedContext<'_>,
    run: bool,
    generator_limit: usize,
    max_target_degree: Option<i128>,
) -> Option<TargetExtremalRayCertificateProbe> {
    if !run {
        return None;
    }
    if let Some(probe) = target_exact_non_extremal_decomposition_probe(sample, target, context) {
        return Some(probe);
    }
    if !sample.is_mori_generator {
        return Some(TargetExtremalRayCertificateProbe {
            status: "skipped_non_mori_generator".to_string(),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: None,
            decomposition_active_generator_count: None,
            decomposition_exact_coefficients: None,
            decomposition_active_generators_nonzero: None,
        });
    }
    if max_target_degree.is_some_and(|max_degree| sample.degree > max_degree) {
        return Some(TargetExtremalRayCertificateProbe {
            status: "skipped_target_degree_limit".to_string(),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: None,
            decomposition_active_generator_count: None,
            decomposition_exact_coefficients: None,
            decomposition_active_generators_nonzero: None,
        });
    }
    if context.degree_bounded_rays.len() > generator_limit {
        return Some(TargetExtremalRayCertificateProbe {
            status: format!(
                "skipped_generator_limit_{generator_limit}_actual_{}",
                context.degree_bounded_rays.len()
            ),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: None,
            decomposition_active_generator_count: None,
            decomposition_exact_coefficients: None,
            decomposition_active_generators_nonzero: None,
        });
    }
    match find_extremal_mori_ray_separator(target, context.degree_bounded_rays) {
        Ok(Some(certificate)) => Some(TargetExtremalRayCertificateProbe {
            status: "certified_exact_extremal_ray".to_string(),
            same_ray_generator_count: Some(certificate.same_ray_generator_count),
            zero_other_generator_count: Some(certificate.zero_other_generator_count),
            positive_other_generator_count: Some(certificate.positive_other_generator_count),
            separator_normal_nonzero: Some(sparse_from_dense(&certificate.separator_normal)),
            decomposition_kind: None,
            decomposition_active_generator_count: None,
            decomposition_exact_coefficients: None,
            decomposition_active_generators_nonzero: None,
        }),
        Ok(None) => Some(TargetExtremalRayCertificateProbe {
            status: "not_certified_as_extremal_ray".to_string(),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: None,
            decomposition_active_generator_count: None,
            decomposition_exact_coefficients: None,
            decomposition_active_generators_nonzero: None,
        }),
        Err(error) => Some(TargetExtremalRayCertificateProbe {
            status: format!("error_{}", status_error_fragment(&error.to_string())),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: None,
            decomposition_active_generator_count: None,
            decomposition_exact_coefficients: None,
            decomposition_active_generators_nonzero: None,
        }),
    }
}

fn target_exact_non_extremal_decomposition_probe(
    sample: &MissingGvTargetSample,
    target: &[i64],
    context: &ValidatedContext<'_>,
) -> Option<TargetExtremalRayCertificateProbe> {
    if !sample.real_cone_decomposable_by_other_generators {
        return None;
    }
    let Some(kind) = sample.real_cone_decomposition_exact_kind.as_deref() else {
        return None;
    };
    let Some(active_generators_sparse) = sample
        .real_cone_decomposition_active_generator_basis_nonzero
        .as_ref()
    else {
        return None;
    };
    let Some(coefficients) = sample.real_cone_decomposition_exact_coefficients.as_ref() else {
        return None;
    };

    match verify_exact_positive_decomposition(
        target,
        active_generators_sparse,
        coefficients,
        context.dimension,
    ) {
        Ok(()) => Some(TargetExtremalRayCertificateProbe {
            status: format!("not_extremal_by_exact_{kind}_decomposition"),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: Some(kind.to_string()),
            decomposition_active_generator_count: Some(active_generators_sparse.len()),
            decomposition_exact_coefficients: Some(coefficients.clone()),
            decomposition_active_generators_nonzero: Some(active_generators_sparse.clone()),
        }),
        Err(error) => Some(TargetExtremalRayCertificateProbe {
            status: format!(
                "error_invalid_exact_decomposition_{}",
                status_error_fragment(&error)
            ),
            same_ray_generator_count: None,
            zero_other_generator_count: None,
            positive_other_generator_count: None,
            separator_normal_nonzero: None,
            decomposition_kind: Some(kind.to_string()),
            decomposition_active_generator_count: Some(active_generators_sparse.len()),
            decomposition_exact_coefficients: Some(coefficients.clone()),
            decomposition_active_generators_nonzero: Some(active_generators_sparse.clone()),
        }),
    }
}

fn verify_exact_positive_decomposition(
    target: &[i64],
    active_generators_sparse: &[Vec<(usize, i64)>],
    coefficients: &[String],
    dimension: usize,
) -> Result<(), String> {
    if target.len() != dimension {
        return Err(format!(
            "target dimension {} does not match context dimension {dimension}",
            target.len()
        ));
    }
    if active_generators_sparse.len() != coefficients.len() {
        return Err(format!(
            "active generator count {} does not match coefficient count {}",
            active_generators_sparse.len(),
            coefficients.len()
        ));
    }

    let zero = MalachiteRational::from(0);
    let mut reconstructed = vec![zero.clone(); dimension];
    for (generator_sparse, coefficient) in active_generators_sparse.iter().zip(coefficients) {
        let coefficient = parse_rational(coefficient)?;
        if coefficient <= zero {
            return Err("decomposition coefficient is not positive".to_string());
        }
        let generator = dense_from_sparse(generator_sparse, dimension)?;
        if generator == target {
            return Err("decomposition uses the target ray itself".to_string());
        }
        for (slot, &entry) in reconstructed.iter_mut().zip(generator.iter()) {
            *slot += coefficient.clone() * MalachiteRational::from(Integer::from(entry));
        }
    }

    for (idx, (actual, &expected)) in reconstructed.iter().zip(target.iter()).enumerate() {
        let expected = MalachiteRational::from(Integer::from(expected));
        if actual != &expected {
            return Err(format!(
                "decomposition coordinate {idx} reconstructs {actual} instead of {expected}"
            ));
        }
    }
    Ok(())
}

fn degree_bounded_ray_context_support_generators(
    ray_context: &[DegreeBoundedMoriRayContextSample],
    max_degree: i128,
    allowed_ambient_support: &HashSet<usize>,
    dimension: usize,
) -> Result<Vec<Vec<i64>>, String> {
    let mut generators = Vec::new();
    let mut seen = HashSet::new();
    for ray in ray_context {
        if ray.degree <= 0 || ray.degree > max_degree {
            continue;
        }
        if !ray
            .ambient_nonzero
            .iter()
            .all(|(idx, _)| allowed_ambient_support.contains(idx))
        {
            continue;
        }
        let basis_ray = dense_from_sparse(&ray.basis_nonzero, dimension)?;
        if seen.insert(basis_ray.clone()) {
            generators.push(basis_ray);
        }
    }
    generators.sort();
    Ok(generators)
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
    let charge_row_permutation_signatures =
        local_charge_row_permutation_signatures(&support.local_charge_basis);
    let charge_row_multiplicities = support
        .local_charge_basis
        .iter()
        .map(|row| local_charge_multiplicities(row))
        .collect::<Vec<_>>();
    let is_calabi_yau_charge = charge_sums.iter().all(|&sum| sum == 0);
    let is_compact_threefold_hypersurface_shape = is_calabi_yau_charge && cy_dim == 3;
    let (cygv_compact_input_status, cygv_compact_input_missing) =
        cygv_compact_input_readiness(is_compact_threefold_hypersurface_shape);
    Ok(Some(LocalCygvHypersurfaceShape {
        q_rows,
        q_cols,
        cy_codim,
        ambient_dim,
        cy_dim,
        charge_sums,
        charge_row_permutation_signatures,
        charge_row_multiplicities,
        is_calabi_yau_charge,
        is_compact_threefold_hypersurface_shape,
        cygv_compact_input_status,
        cygv_compact_input_missing,
    }))
}

fn local_charge_row_permutation_signatures(local_charge_basis: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let mut rows = local_charge_basis
        .iter()
        .map(|row| {
            let mut sorted = row.clone();
            sorted.sort_unstable();
            sorted
        })
        .collect::<Vec<_>>();
    rows.sort();
    rows
}

fn local_charge_signature_key(local_charge_basis: &[Vec<i64>]) -> String {
    local_charge_row_permutation_signatures(local_charge_basis)
        .iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join("|")
}

fn local_charge_multiplicities(row: &[i64]) -> Vec<LocalChargeMultiplicity> {
    let mut counts = BTreeMap::new();
    for &charge in row {
        *counts.entry(charge).or_insert(0usize) += 1;
    }
    counts
        .into_iter()
        .map(|(charge, count)| LocalChargeMultiplicity { charge, count })
        .collect()
}

fn cygv_compact_input_readiness(
    is_compact_threefold_hypersurface_shape: bool,
) -> (String, Vec<String>) {
    if !is_compact_threefold_hypersurface_shape {
        return (
            "not_compact_threefold_hypersurface_shape".to_string(),
            Vec::new(),
        );
    }
    (
        "shape_only_missing_source_derived_cygv_inputs".to_string(),
        vec![
            "local_semigroup_generators".to_string(),
            "local_grading_vector".to_string(),
            "local_q_matrix_orientation_and_phase".to_string(),
            "local_intersection_tensor".to_string(),
            "local_chamber_certificate".to_string(),
            "target_class_to_local_semigroup_coordinate".to_string(),
        ],
    )
}

fn local_cygv_input_skeleton(
    sample: &MissingGvTargetSample,
    support: Option<&OriginCircuitAffineSupportSample>,
) -> Result<Option<LocalCygvInputSkeleton>, String> {
    let Some(support) = support else {
        return Ok(None);
    };
    let Some(first_row) = support.local_charge_basis.first() else {
        return Ok(None);
    };
    let support_len = first_row.len();
    if support_len == 0 {
        return Ok(None);
    }
    if support
        .local_charge_basis
        .iter()
        .any(|row| row.len() != support_len)
    {
        return Err("local charge basis rows have inconsistent lengths".to_string());
    }
    let local_q_matrix_rows = transpose_local_charge_basis(&support.local_charge_basis);
    let support_point_indices = sample
        .origin_circuit_first_witness
        .as_ref()
        .map(|witness| {
            witness
                .relation_points
                .iter()
                .map(|point| point.point_index)
                .collect::<Vec<_>>()
        })
        .filter(|indices| indices.len() == support_len)
        .unwrap_or_else(|| {
            support
                .local_coordinates
                .iter()
                .map(|point| point.point_index)
                .collect()
        });
    let support_contains_origin_point = support_point_indices.contains(&0);
    let local_cygv_origin_point_status = if support_contains_origin_point {
        "local_gkz_relation_includes_origin_point_requires_phase_mapping"
    } else {
        "local_support_has_no_origin_point"
    }
    .to_string();
    let target_relation_coefficients =
        sample
            .origin_circuit_first_witness
            .as_ref()
            .and_then(|witness| {
                (witness.relation_points.len() == support_len).then(|| {
                    witness
                        .relation_points
                        .iter()
                        .map(|point| point.coefficient)
                        .collect::<Vec<_>>()
                })
            });
    let (origin_point_relation_coefficient, local_cytools_origin_circuit_status) =
        local_cytools_origin_circuit_status(
            &support_point_indices,
            target_relation_coefficients.as_deref(),
        );
    let target_relation_in_charge_basis = target_relation_coefficients
        .as_ref()
        .map(|target| {
            relation_coordinates_in_local_charge_basis(target, &support.local_charge_basis)
        })
        .transpose()?
        .flatten();
    let target_relation_status = match (
        target_relation_coefficients.as_ref(),
        target_relation_in_charge_basis.as_ref(),
    ) {
        (Some(_), Some(_)) => "target_relation_integral_in_local_charge_basis",
        (Some(_), None) => "target_relation_not_in_local_charge_basis",
        (None, _) => "target_relation_unavailable",
    }
    .to_string();
    let orientation_candidates = local_cygv_orientation_candidates(
        &local_q_matrix_rows,
        target_relation_in_charge_basis.as_deref(),
    );
    let (local_q_matrix_orientation_candidate, local_q_matrix_orientation_status) =
        local_cygv_q_matrix_orientation_candidate(&orientation_candidates);
    let (
        local_cygv_q_matrix_rows_candidate,
        local_cygv_wrapper_q_matrix_candidate,
        local_cygv_q_matrix_layout_status,
    ) = local_cygv_q_matrix_layout_candidate(
        &orientation_candidates,
        local_q_matrix_orientation_candidate,
    );
    let (local_grading_vector_candidate, local_grading_vector_status) =
        local_cygv_grading_vector_candidate(&orientation_candidates);
    let mut remaining_uncertified_inputs = vec!["local_semigroup_generators".to_string()];
    if local_grading_vector_candidate.is_none() {
        remaining_uncertified_inputs.push("local_grading_vector".to_string());
    }
    if local_q_matrix_orientation_candidate.is_some() {
        remaining_uncertified_inputs.push("local_q_matrix_phase".to_string());
    } else {
        remaining_uncertified_inputs.push("local_q_matrix_orientation_and_phase".to_string());
    }
    remaining_uncertified_inputs.extend([
        "local_intersection_tensor".to_string(),
        "local_chamber_certificate".to_string(),
    ]);
    Ok(Some(LocalCygvInputSkeleton {
        support_point_indices,
        support_contains_origin_point,
        local_cygv_origin_point_status,
        origin_point_relation_coefficient,
        local_cytools_origin_circuit_status,
        local_q_matrix_rows,
        target_relation_coefficients,
        target_relation_in_charge_basis,
        target_relation_status,
        orientation_candidates,
        local_q_matrix_orientation_candidate,
        local_q_matrix_orientation_status,
        local_cygv_q_matrix_rows_candidate,
        local_cygv_wrapper_q_matrix_candidate,
        local_cygv_q_matrix_layout_status,
        local_grading_vector_candidate,
        local_grading_vector_status,
        remaining_uncertified_inputs,
    }))
}

fn local_cygv_q_matrix_layout_candidate(
    orientation_candidates: &[LocalCygvOrientationCandidate],
    orientation_candidate: Option<i64>,
) -> (Option<Vec<Vec<i64>>>, Option<Vec<Vec<i64>>>, String) {
    let Some(sign) = orientation_candidate else {
        return (
            None,
            None,
            "local_q_matrix_layout_blocked_no_orientation".to_string(),
        );
    };
    let Some(candidate) = orientation_candidates
        .iter()
        .find(|candidate| candidate.overall_charge_basis_sign == sign)
    else {
        return (
            None,
            None,
            "local_q_matrix_layout_missing_selected_orientation".to_string(),
        );
    };
    if candidate.local_q_matrix_rows.is_empty() {
        return (None, None, "local_q_matrix_layout_empty".to_string());
    }
    let Some(width) = candidate.local_q_matrix_rows.first().map(Vec::len) else {
        return (None, None, "local_q_matrix_layout_empty".to_string());
    };
    if width == 0 {
        return (None, None, "local_q_matrix_layout_zero_rank".to_string());
    }
    if candidate
        .local_q_matrix_rows
        .iter()
        .any(|row| row.len() != width)
    {
        return (
            None,
            None,
            "local_q_matrix_layout_inconsistent_row_widths".to_string(),
        );
    }
    let wrapper_q_matrix = transpose_local_charge_basis(&candidate.local_q_matrix_rows);
    (
        Some(candidate.local_q_matrix_rows.clone()),
        Some(wrapper_q_matrix),
        "source_derived_oriented_q_matrix_layout".to_string(),
    )
}

fn local_cytools_origin_circuit_status(
    support_point_indices: &[usize],
    target_relation_coefficients: Option<&[i64]>,
) -> (Option<i64>, String) {
    let Some(origin_position) = support_point_indices.iter().position(|&point| point == 0) else {
        return (None, "not_cytools_origin_circuit_support".to_string());
    };
    let Some(coefficients) = target_relation_coefficients else {
        return (None, "origin_relation_coefficients_unavailable".to_string());
    };
    if coefficients.len() != support_point_indices.len() {
        return (
            None,
            "origin_relation_coefficients_support_mismatch".to_string(),
        );
    }
    let origin_coefficient = coefficients[origin_position];
    let status = if origin_coefficient < 0 {
        "source_cytools_retains_negative_origin_coefficient"
    } else if origin_coefficient == 0 {
        "source_cytools_skips_zero_origin_coefficient"
    } else {
        "source_cytools_rejects_positive_origin_coefficient"
    };
    (Some(origin_coefficient), status.to_string())
}

fn local_cygv_q_matrix_orientation_candidate(
    orientation_candidates: &[LocalCygvOrientationCandidate],
) -> (Option<i64>, String) {
    let eligible = orientation_candidates
        .iter()
        .filter(|candidate| {
            candidate.target_coordinate_is_nonnegative == Some(true)
                && candidate.target_coordinate_is_primitive == Some(true)
        })
        .collect::<Vec<_>>();
    match eligible.as_slice() {
        [] => {
            if orientation_candidates
                .iter()
                .any(|candidate| candidate.target_coordinate.is_none())
            {
                (
                    None,
                    "local_q_orientation_target_coordinate_unavailable".to_string(),
                )
            } else {
                (
                    None,
                    "local_q_orientation_blocked_no_positive_primitive_target".to_string(),
                )
            }
        }
        [candidate] => (
            Some(candidate.overall_charge_basis_sign),
            "source_derived_target_positive_orientation".to_string(),
        ),
        _ => (
            None,
            "local_q_orientation_ambiguous_multiple_positive_primitive_targets".to_string(),
        ),
    }
}

fn local_cygv_grading_vector_candidate(
    orientation_candidates: &[LocalCygvOrientationCandidate],
) -> (Option<Vec<i64>>, String) {
    for candidate in orientation_candidates {
        if candidate.target_coordinate_is_nonnegative != Some(true)
            || candidate.target_coordinate_is_primitive != Some(true)
        {
            continue;
        }
        let Some(direction) = candidate.target_primitive_direction.as_ref() else {
            continue;
        };
        if direction.len() == 1 && direction[0] > 0 {
            return (
                Some(vec![1]),
                "source_derived_primitive_one_parameter_grading".to_string(),
            );
        }
        return (
            None,
            "local_grading_requires_higher_rank_dual_cone_certificate".to_string(),
        );
    }
    if orientation_candidates
        .iter()
        .any(|candidate| candidate.target_coordinate.is_none())
    {
        (
            None,
            "local_grading_target_coordinate_unavailable".to_string(),
        )
    } else {
        (
            None,
            "local_grading_blocked_no_positive_primitive_target".to_string(),
        )
    }
}

fn local_cygv_orientation_candidates(
    local_q_matrix_rows: &[Vec<i64>],
    target_relation_in_charge_basis: Option<&[i64]>,
) -> Vec<LocalCygvOrientationCandidate> {
    [-1, 1]
        .into_iter()
        .map(|sign| {
            let oriented_q = local_q_matrix_rows
                .iter()
                .map(|row| row.iter().map(|&entry| sign * entry).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let target_coordinate = target_relation_in_charge_basis
                .map(|coordinate| coordinate.iter().map(|&entry| sign * entry).collect());
            let target_coordinate_is_nonnegative = target_coordinate
                .as_ref()
                .map(|coordinate: &Vec<i64>| coordinate.iter().all(|&entry| entry >= 0));
            let target_coordinate_gcd = target_coordinate
                .as_ref()
                .map(|coordinate: &Vec<i64>| gcd_list_int(coordinate).abs());
            let target_coordinate_is_primitive = target_coordinate_gcd.map(|gcd| gcd == 1);
            let target_primitive_direction = target_coordinate
                .as_ref()
                .zip(target_coordinate_gcd)
                .and_then(|(coordinate, gcd)| {
                    (gcd > 0).then(|| {
                        coordinate
                            .iter()
                            .map(|&entry| entry / gcd)
                            .collect::<Vec<_>>()
                    })
                });
            let positive_unit_generator_negative_intersections =
                (oriented_q.first().map_or(0, Vec::len) == 1)
                    .then(|| oriented_q.iter().filter(|row| row[0] < 0).count());
            let positive_unit_generator_omega_bucket =
                positive_unit_generator_negative_intersections.map(cygv_omega_bucket);
            let target_candidate_status = local_cygv_target_candidate_status(
                target_coordinate_is_nonnegative,
                target_coordinate_is_primitive,
                positive_unit_generator_omega_bucket.as_deref(),
            );
            LocalCygvOrientationCandidate {
                overall_charge_basis_sign: sign,
                local_q_matrix_rows: oriented_q,
                target_coordinate,
                target_candidate_status,
                target_coordinate_is_nonnegative,
                target_coordinate_gcd,
                target_coordinate_is_primitive,
                target_primitive_direction,
                positive_unit_generator_negative_intersections,
                positive_unit_generator_omega_bucket,
            }
        })
        .collect()
}

fn local_cygv_target_candidate_status(
    target_coordinate_is_nonnegative: Option<bool>,
    target_coordinate_is_primitive: Option<bool>,
    positive_unit_generator_omega_bucket: Option<&str>,
) -> String {
    let Some(nonnegative) = target_coordinate_is_nonnegative else {
        return "target_coordinate_unavailable".to_string();
    };
    if !nonnegative {
        return "target_not_in_nonnegative_local_semigroup".to_string();
    }
    if target_coordinate_is_primitive == Some(false) {
        return "target_nonprimitive_local_multiple".to_string();
    }
    match positive_unit_generator_omega_bucket {
        Some(bucket) if bucket.starts_with("ignored") => {
            "target_positive_but_ignored_by_cygv_omega_bucket".to_string()
        }
        Some(_) => "target_primitive_positive_supported_by_cygv_omega_bucket".to_string(),
        None => "target_positive_but_omega_bucket_unavailable".to_string(),
    }
}

fn transpose_local_charge_basis(local_charge_basis: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let Some(width) = local_charge_basis.first().map(Vec::len) else {
        return Vec::new();
    };
    (0..width)
        .map(|row| local_charge_basis.iter().map(|basis| basis[row]).collect())
        .collect()
}

fn relation_coordinates_in_local_charge_basis(
    target: &[i64],
    rows: &[Vec<i64>],
) -> Result<Option<Vec<i64>>, String> {
    if rows.is_empty() {
        return Ok(target.iter().all(|&value| value == 0).then(Vec::new));
    }
    let dim = rows[0].len();
    if target.len() != dim {
        return Err("target relation dimension does not match local charge basis".to_string());
    }
    if rows.iter().any(|row| row.len() != dim) {
        return Err("local charge basis rows have inconsistent dimensions".to_string());
    }
    for columns in column_combinations(dim, rows.len()) {
        let matrix = columns
            .iter()
            .map(|&col| {
                rows.iter()
                    .map(|row| MalachiteRational::from(Integer::from(row[col])))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let rhs = columns
            .iter()
            .map(|&col| MalachiteRational::from(Integer::from(target[col])))
            .collect::<Vec<_>>();
        let Some(coordinates) = solve_linear_system_rational(&matrix, &rhs) else {
            continue;
        };
        if !local_charge_coordinates_match_target(&coordinates, target, rows) {
            continue;
        }
        let coordinates = coordinates
            .into_iter()
            .map(|coordinate| {
                if coordinate.denominator_ref() != &1u32 {
                    return Err(
                        "target relation has non-integral local charge coordinates".to_string()
                    );
                }
                let integer = Integer::try_from(coordinate).map_err(|_| {
                    "target relation local charge coordinate is not integral".to_string()
                })?;
                i64::try_from(&integer).map_err(|_| {
                    "target relation local charge coordinate does not fit in i64".to_string()
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        return Ok(Some(coordinates));
    }
    Ok(None)
}

fn local_charge_coordinates_match_target(
    coordinates: &[MalachiteRational],
    target: &[i64],
    rows: &[Vec<i64>],
) -> bool {
    (0..target.len()).all(|col| {
        let reconstructed = coordinates.iter().zip(rows.iter()).fold(
            MalachiteRational::from(0),
            |acc, (coordinate, row)| {
                acc + coordinate * MalachiteRational::from(Integer::from(row[col]))
            },
        );
        reconstructed == MalachiteRational::from(Integer::from(target[col]))
    })
}

fn column_combinations(n_columns: usize, size: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current = Vec::with_capacity(size);
    push_column_combinations(0, n_columns, size, &mut current, &mut out);
    out
}

fn push_column_combinations(
    start: usize,
    n_columns: usize,
    size: usize,
    current: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if current.len() == size {
        out.push(current.clone());
        return;
    }
    let remaining = size - current.len();
    if remaining > n_columns.saturating_sub(start) {
        return;
    }
    for col in start..=n_columns - remaining {
        current.push(col);
        push_column_combinations(col + 1, n_columns, size, current, out);
        current.pop();
    }
}

fn origin_circuit_affine_support_with_coordinates(
    sample: &MissingGvTargetSample,
) -> Result<Option<OriginCircuitAffineSupportSample>, String> {
    let Some(mut support) = sample.origin_circuit_affine_support.clone() else {
        return Ok(None);
    };
    if !support.local_coordinates.is_empty() {
        return Ok(Some(support));
    }
    let Some(witness) = sample.origin_circuit_first_witness.as_ref() else {
        return Ok(Some(support));
    };
    if witness.relation_points.is_empty() {
        return Ok(Some(support));
    }

    let points = witness
        .relation_points
        .iter()
        .map(|point| Point::new(point.coordinates.clone()))
        .collect::<Vec<_>>();
    let relation = witness
        .relation_points
        .iter()
        .map(|point| point.coefficient)
        .collect::<Vec<_>>();
    let diagnostic = diagnose_affine_toric_circuit(&relation, &points)
        .map_err(|error| format!("failed to reconstruct local affine coordinates: {error}"))?
        .ok_or_else(|| {
            "failed to reconstruct local affine coordinates: relation is not an affine circuit"
                .to_string()
        })?;

    support.local_coordinates = diagnostic
        .local_coordinates
        .iter()
        .map(|point| {
            let ambient_point = witness
                .relation_points
                .get(point.point_index)
                .ok_or_else(|| "local coordinate point index is out of witness bounds".to_string())?
                .point_index;
            Ok(OriginCircuitLocalCoordinateSample {
                point_index: ambient_point,
                coordinates: point.coordinates.clone(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok(Some(support))
}

fn validate_context<'a>(
    context: &'a CorrectedChamberGvContext,
) -> Result<ValidatedContext<'a>, String> {
    if !matches!(context.schema_version, 1 | 2 | 3) {
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
    let degree_bounded_ray_context = context
        .degree_bounded_mori_ray_context_for_missing
        .as_deref();
    if context.schema_version >= 3 && degree_bounded_ray_context.is_none() {
        return Err("schema-3 context is missing degree-bounded Mori ray context".to_string());
    }
    if let Some(ray_context) = degree_bounded_ray_context {
        for (idx, sample) in ray_context.iter().enumerate() {
            let mut ambient_seen = HashSet::new();
            for &(ambient_idx, _) in &sample.ambient_nonzero {
                if !ambient_seen.insert(ambient_idx) {
                    return Err(format!(
                        "degree-bounded Mori ray context sample {idx} has duplicate ambient coordinate {ambient_idx}"
                    ));
                }
            }
            let basis_ray = dense_from_sparse(&sample.basis_nonzero, dimension).map_err(|e| {
                format!(
                    "degree-bounded Mori ray context sample {idx} has invalid basis support: {e}"
                )
            })?;
            let computed_degree = curve_degree(&basis_ray, grading).map_err(|e| {
                format!(
                    "degree-bounded Mori ray context sample {idx} has invalid basis degree: {e}"
                )
            })?;
            if computed_degree != sample.degree {
                return Err(format!(
                    "degree-bounded Mori ray context sample {idx} declares degree {} but computes {computed_degree}",
                    sample.degree
                ));
            }
            if sample.degree > degree_bound {
                return Err(format!(
                    "degree-bounded Mori ray context sample {idx} has degree {} above bound {degree_bound}",
                    sample.degree
                ));
            }
        }
    }
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
        degree_bounded_ray_context,
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
    degree_bounded_ray_context: Option<&'a [DegreeBoundedMoriRayContextSample]>,
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
    certify_origin_support_domains: bool,
    origin_support_certificate_limit: usize,
    certify_target_extremal_rays: bool,
    target_extremal_generator_limit: usize,
    target_extremal_max_degree: Option<i128>,
    measure_cygv_semigroups: bool,
    probe_cygv_path_history: bool,
    run_lower_seed_diamonds: bool,
    measure_cygv_degree_ladder: bool,
    cygv_degree_ladder_max_degree: Option<i128>,
    semigroup_measure_max_target_degree: Option<i128>,
    semigroup_measure_max_seed_count: Option<usize>,
    semigroup_measurement_cache: &mut HashMap<i128, Result<CygvSemigroupDegreeMeasurement, String>>,
    semigroup_ladder_cache: &mut HashMap<
        (i128, i128, Option<usize>),
        Result<Vec<CygvSemigroupDegreeLadderStep>, String>,
    >,
    element_limit: usize,
) -> TargetReport {
    let exact_kind = sample.real_cone_decomposition_exact_kind.clone();
    let active_generator_count = sample
        .real_cone_decomposition_active_generator_basis_nonzero
        .as_ref()
        .map(Vec::len)
        .or(sample.real_cone_decomposition_active_generators);
    let origin_circuit_affine_support = match origin_circuit_affine_support_with_coordinates(sample)
    {
        Ok(support) => support,
        Err(error) => {
            eprintln!(
                "[WARN] target {index}: {error}; preserving serialized origin-circuit support"
            );
            sample.origin_circuit_affine_support.clone()
        }
    };
    let local_cygv_input_skeleton =
        match local_cygv_input_skeleton(sample, origin_circuit_affine_support.as_ref()) {
            Ok(skeleton) => skeleton,
            Err(error) => {
                eprintln!("[WARN] target {index}: {error}; omitting local cygv input skeleton");
                None
            }
        };
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
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape: None,
                local_cygv_input_skeleton,
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
                target_extremal_ray_certificate: None,
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
                origin_relation_support_generator_count: None,
                origin_shared_facet_generator_count: None,
                origin_facet_union_generator_count: None,
                origin_relation_support_face_certificate_status: None,
                origin_shared_facet_face_certificate_status: None,
                origin_facet_union_face_certificate_status: None,
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
                cygv_semigroup_degree_ladder: None,
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
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape,
                local_cygv_input_skeleton,
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
                target_extremal_ray_certificate: None,
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
                origin_relation_support_generator_count: None,
                origin_shared_facet_generator_count: None,
                origin_facet_union_generator_count: None,
                origin_relation_support_face_certificate_status: None,
                origin_shared_facet_face_certificate_status: None,
                origin_facet_union_face_certificate_status: None,
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
                cygv_semigroup_degree_ladder: None,
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
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape,
                local_cygv_input_skeleton,
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
                target_extremal_ray_certificate: None,
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
                origin_relation_support_generator_count: None,
                origin_shared_facet_generator_count: None,
                origin_facet_union_generator_count: None,
                origin_relation_support_face_certificate_status: None,
                origin_shared_facet_face_certificate_status: None,
                origin_facet_union_face_certificate_status: None,
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
                cygv_semigroup_degree_ladder: None,
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
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape,
                local_cygv_input_skeleton,
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
                target_extremal_ray_certificate: None,
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
                origin_relation_support_generator_count: None,
                origin_shared_facet_generator_count: None,
                origin_facet_union_generator_count: None,
                origin_relation_support_face_certificate_status: None,
                origin_shared_facet_face_certificate_status: None,
                origin_facet_union_face_certificate_status: None,
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
                cygv_semigroup_degree_ladder: None,
                cygv_path_history_probe: None,
            };
        }
    };
    let origin_counts = origin_circuit_ambient_generator_counts(sample, context);
    let origin_certificates = if certify_origin_support_domains {
        origin_circuit_ambient_support_face_certificate_statuses(
            sample,
            context,
            origin_support_certificate_limit,
        )
    } else {
        OriginCircuitAmbientSupportCertificateStatuses {
            relation_support: None,
            shared_facet: None,
            facet_union: None,
        }
    };
    let target_extremal_ray_certificate = target_extremal_ray_certificate_probe(
        sample,
        &target,
        context,
        certify_target_extremal_rays,
        target_extremal_generator_limit,
        target_extremal_max_degree,
    );
    let base = TargetReport {
        index,
        degree: sample.degree,
        generators_le_degree: sample.generators_le_degree,
        is_mori_generator: sample.is_mori_generator,
        origin_circuit_pattern: sample.origin_circuit_pattern.clone(),
        origin_circuit_witness_count: sample.origin_circuit_witness_count,
        origin_circuit_first_witness: sample.origin_circuit_first_witness.clone(),
        origin_circuit_affine_support,
        local_cygv_hypersurface_shape,
        local_cygv_input_skeleton,
        cms_general_divisor_shape_candidates: sample.cms_general_divisor_shape_candidates.clone(),
        cms_general_divisor_intersection_checks: sample
            .cms_general_divisor_intersection_checks
            .clone(),
        branch_diagnostic: sample.branch_diagnostic.clone(),
        real_cone_decomposable_by_other_generators: sample
            .real_cone_decomposable_by_other_generators,
        ambient_nonzero: sample.ambient_nonzero.clone(),
        basis_nonzero: sample.basis_nonzero.clone(),
        target_extremal_ray_certificate,
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
        origin_relation_support_generator_count: origin_counts.relation_support,
        origin_shared_facet_generator_count: origin_counts.shared_facet,
        origin_facet_union_generator_count: origin_counts.facet_union,
        origin_relation_support_face_certificate_status: origin_certificates.relation_support,
        origin_shared_facet_face_certificate_status: origin_certificates.shared_facet,
        origin_facet_union_face_certificate_status: origin_certificates.facet_union,
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
        cygv_semigroup_degree_ladder: None,
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
    let cygv_semigroup_degree_ladder = if measure_cygv_degree_ladder {
        if semigroup_measure_max_target_degree.is_some_and(|max_degree| sample.degree > max_degree)
        {
            Some(vec![CygvSemigroupDegreeLadderStep {
                degree: sample.degree,
                effective_seed_count: 0,
                reduced_seed_count: None,
                status: "skipped_target_degree_limit".to_string(),
                element_count: None,
                elapsed_ms: None,
                error: None,
            }])
        } else {
            match cygv_degree_ladder_max_degree {
                Some(max_ladder_degree) => Some(measure_cygv_semigroup_degree_ladder(
                    sample,
                    context,
                    max_ladder_degree,
                    semigroup_measure_max_seed_count,
                    semigroup_ladder_cache,
                )),
                None => Some(vec![CygvSemigroupDegreeLadderStep {
                    degree: sample.degree,
                    effective_seed_count: 0,
                    reduced_seed_count: None,
                    status: "skipped_missing_ladder_max_degree".to_string(),
                    element_count: None,
                    elapsed_ms: None,
                    error: Some(
                        "--measure-cygv-degree-ladder requires --cygv-degree-ladder-max-degree"
                            .to_string(),
                    ),
                }]),
            }
        }
    } else {
        None
    };
    let cygv_path_history_probe = if measure_cygv_semigroups || probe_cygv_path_history {
        if semigroup_measure_max_target_degree.is_some_and(|max_degree| sample.degree > max_degree)
        {
            None
        } else {
            Some(cygv_path_history_probe(
                sample,
                context,
                &target,
                run_lower_seed_diamonds,
                element_limit,
                semigroup_measure_max_seed_count,
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
            cygv_semigroup_degree_ladder,
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
            cygv_semigroup_degree_ladder,
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

fn measure_cygv_semigroup_degree_ladder(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    max_ladder_degree: i128,
    max_seed_count: Option<usize>,
    cache: &mut HashMap<
        (i128, i128, Option<usize>),
        Result<Vec<CygvSemigroupDegreeLadderStep>, String>,
    >,
) -> Vec<CygvSemigroupDegreeLadderStep> {
    if max_ladder_degree <= 0 {
        return vec![CygvSemigroupDegreeLadderStep {
            degree: max_ladder_degree,
            effective_seed_count: 0,
            reduced_seed_count: None,
            status: "error".to_string(),
            element_count: None,
            elapsed_ms: None,
            error: Some("cygv degree ladder max degree must be positive".to_string()),
        }];
    }
    let effective_max = sample.degree.min(max_ladder_degree);
    let key = (sample.degree, effective_max, max_seed_count);
    if !cache.contains_key(&key) {
        let measurement = measure_cygv_semigroup_degree_ladder_uncached(
            sample.degree,
            effective_max,
            context,
            max_seed_count,
        );
        cache.insert(key, measurement);
    }
    cache
        .get(&key)
        .expect("ladder measurement was inserted above")
        .clone()
        .unwrap_or_else(|error| {
            vec![CygvSemigroupDegreeLadderStep {
                degree: effective_max,
                effective_seed_count: 0,
                reduced_seed_count: None,
                status: "error".to_string(),
                element_count: None,
                elapsed_ms: None,
                error: Some(error),
            }]
        })
}

fn measure_cygv_semigroup_degree_ladder_uncached(
    target_degree: i128,
    max_ladder_degree: i128,
    context: &ValidatedContext<'_>,
    max_seed_count: Option<usize>,
) -> Result<Vec<CygvSemigroupDegreeLadderStep>, String> {
    let mut all_seeds = Vec::new();
    let mut seen = HashSet::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > target_degree {
            continue;
        }
        if seen.insert(ray.clone()) {
            all_seeds.push((degree, ray.clone()));
        }
    }

    let mut steps = Vec::new();
    for degree in 1..=max_ladder_degree {
        let seeds = all_seeds
            .iter()
            .filter_map(|(seed_degree, seed)| (*seed_degree <= degree).then(|| seed.clone()))
            .collect::<Vec<_>>();
        let effective_seed_count = seeds.len();
        if effective_seed_count == 0 {
            steps.push(CygvSemigroupDegreeLadderStep {
                degree,
                effective_seed_count,
                reduced_seed_count: None,
                status: "skipped_empty_seed_set".to_string(),
                element_count: None,
                elapsed_ms: None,
                error: None,
            });
            continue;
        }
        let reduced_seed_count = match cygv_pair_reduced_seed_generators(&seeds) {
            Ok(reduced) => Some(reduced.len()),
            Err(error) => {
                steps.push(CygvSemigroupDegreeLadderStep {
                    degree,
                    effective_seed_count,
                    reduced_seed_count: None,
                    status: "error".to_string(),
                    element_count: None,
                    elapsed_ms: None,
                    error: Some(format!("cygv seed reduction failed: {error}")),
                });
                continue;
            }
        };
        if max_seed_count.is_some_and(|limit| effective_seed_count > limit) {
            steps.push(CygvSemigroupDegreeLadderStep {
                degree,
                effective_seed_count,
                reduced_seed_count,
                status: "skipped_seed_limit".to_string(),
                element_count: None,
                elapsed_ms: None,
                error: None,
            });
            continue;
        }
        let max_deg = match u32::try_from(degree) {
            Ok(max_deg) => max_deg,
            Err(_) => {
                steps.push(CygvSemigroupDegreeLadderStep {
                    degree,
                    effective_seed_count,
                    reduced_seed_count,
                    status: "error".to_string(),
                    element_count: None,
                    elapsed_ms: None,
                    error: Some(format!("degree {degree} does not fit in u32")),
                });
                continue;
            }
        };
        let started = std::time::Instant::now();
        match measure_cygv_semigroup_size(&seeds, context.grading, max_deg) {
            Ok(element_count) => steps.push(CygvSemigroupDegreeLadderStep {
                degree,
                effective_seed_count,
                reduced_seed_count,
                status: "measured_cygv_semigroup".to_string(),
                element_count: Some(element_count),
                elapsed_ms: Some(started.elapsed().as_millis()),
                error: None,
            }),
            Err(error) => steps.push(CygvSemigroupDegreeLadderStep {
                degree,
                effective_seed_count,
                reduced_seed_count,
                status: "error".to_string(),
                element_count: None,
                elapsed_ms: Some(started.elapsed().as_millis()),
                error: Some(error),
            }),
        }
    }
    Ok(steps)
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
    run_lower_seed_diamonds: bool,
    element_limit: usize,
    max_seed_count: Option<usize>,
) -> CygvPathHistoryProbe {
    match cygv_path_history_probe_inner(
        sample,
        context,
        target,
        run_lower_seed_diamonds,
        element_limit,
        max_seed_count,
    ) {
        Ok(probe) => probe,
        Err(error) => CygvPathHistoryProbe {
            status: format!("error: {error}"),
            seed_count: None,
            reduced_seed_count: None,
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            target_in_closure: None,
            previous_level_count: cygv_previous_level_count(context.dimension),
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_counts_complete: false,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: "not_run".to_string(),
            lower_seed_decomposition_term_count: None,
            lower_seed_decomposition_terms_nonzero: None,
            lower_seed_decomposition_error: Some(
                "path-history probe failed before seed decomposition".to_string(),
            ),
            lower_seed_diamond_status: None,
            lower_seed_diamond_element_count: None,
            lower_seed_diamond_gv: None,
            lower_seed_diamond_error: None,
        },
    }
}

fn cygv_path_history_probe_inner(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    target: &[i64],
    run_lower_seed_diamonds: bool,
    element_limit: usize,
    max_seed_count: Option<usize>,
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
            seed_count: Some(0),
            reduced_seed_count: Some(0),
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            target_in_closure: None,
            previous_level_count,
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_counts_complete: false,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: "skipped_empty_seed_set".to_string(),
            lower_seed_decomposition_term_count: None,
            lower_seed_decomposition_terms_nonzero: None,
            lower_seed_decomposition_error: None,
            lower_seed_diamond_status: None,
            lower_seed_diamond_element_count: None,
            lower_seed_diamond_gv: None,
            lower_seed_diamond_error: None,
        });
    }
    if max_seed_count.is_some_and(|limit| seeds.len() > limit) {
        return Ok(CygvPathHistoryProbe {
            status: "skipped_seed_limit".to_string(),
            seed_count: Some(seeds.len()),
            reduced_seed_count: None,
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            target_in_closure: None,
            previous_level_count,
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_counts_complete: false,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: "skipped_seed_limit".to_string(),
            lower_seed_decomposition_term_count: None,
            lower_seed_decomposition_terms_nonzero: None,
            lower_seed_decomposition_error: None,
            lower_seed_diamond_status: Some("skipped_seed_limit".to_string()),
            lower_seed_diamond_element_count: None,
            lower_seed_diamond_gv: None,
            lower_seed_diamond_error: None,
        });
    }
    let reduced_seed_count = cygv_pair_reduced_seed_generators(&seeds)
        .map_err(|error| format!("cygv seed reduction failed: {error}"))?
        .len();
    let lower_seed_decomposition = lower_seed_decomposition_probe(target, &seeds, 4);
    let lower_seed_diamond = lower_seed_diamond_probe(
        target,
        &lower_seed_decomposition,
        context,
        run_lower_seed_diamonds,
        element_limit,
    );

    let closure =
        bounded_cygv_semigroup_closure(&seeds, context.grading, sample.degree, element_limit)?;
    let target_in_closure = closure.elements.contains(target);
    let selected_degree_vec =
        cygv_previous_window_degrees(&closure.degree_counts, sample.degree, previous_level_count);
    let selected_degrees = selected_degree_vec.iter().copied().collect::<HashSet<_>>();
    let predecessor_stats = cygv_path_predecessor_stats(
        &closure.elements,
        context.grading,
        target,
        &selected_degrees,
    )?;
    if !closure.completed {
        return Ok(CygvPathHistoryProbe {
            status: closure.status,
            seed_count: Some(seeds.len()),
            reduced_seed_count: Some(reduced_seed_count),
            closure_element_count: Some(closure.elements.len()),
            previous_window_degrees: selected_degree_vec,
            closure_degree_counts: closure.degree_counts,
            target_in_closure: Some(target_in_closure),
            previous_level_count,
            previous_window_degree_count: Some(selected_degrees.len()),
            previous_window_element_count: Some(predecessor_stats.previous_window_element_count),
            predecessor_counts_complete: false,
            predecessor_difference_count: Some(predecessor_stats.predecessor_difference_count),
            improving_predecessor_difference_count: Some(
                predecessor_stats.improving_predecessor_difference_count,
            ),
            closest_series_distance: Some(format!("{:.6}", predecessor_stats.closest_distance)),
            closest_series_predecessor_nonzero: predecessor_stats
                .closest_predecessor
                .as_deref()
                .map(sparse_from_dense),
            closest_series_difference_nonzero: predecessor_stats
                .closest_difference
                .as_deref()
                .map(sparse_from_dense),
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: lower_seed_decomposition.status,
            lower_seed_decomposition_term_count: lower_seed_decomposition.term_count,
            lower_seed_decomposition_terms_nonzero: lower_seed_decomposition.terms_nonzero,
            lower_seed_decomposition_error: lower_seed_decomposition.error,
            lower_seed_diamond_status: lower_seed_diamond.status,
            lower_seed_diamond_element_count: lower_seed_diamond.element_count,
            lower_seed_diamond_gv: lower_seed_diamond.gv,
            lower_seed_diamond_error: lower_seed_diamond.error,
        });
    }

    Ok(CygvPathHistoryProbe {
        status: "completed_bounded_closure".to_string(),
        seed_count: Some(seeds.len()),
        reduced_seed_count: Some(reduced_seed_count),
        closure_element_count: Some(closure.elements.len()),
        closure_degree_counts: closure.degree_counts,
        target_in_closure: Some(target_in_closure),
        previous_level_count,
        previous_window_degrees: selected_degree_vec,
        previous_window_degree_count: Some(selected_degrees.len()),
        previous_window_element_count: Some(predecessor_stats.previous_window_element_count),
        predecessor_counts_complete: true,
        predecessor_difference_count: Some(predecessor_stats.predecessor_difference_count),
        improving_predecessor_difference_count: Some(
            predecessor_stats.improving_predecessor_difference_count,
        ),
        closest_series_distance: Some(format!("{:.6}", predecessor_stats.closest_distance)),
        closest_series_predecessor_nonzero: predecessor_stats
            .closest_predecessor
            .as_deref()
            .map(sparse_from_dense),
        closest_series_difference_nonzero: predecessor_stats
            .closest_difference
            .as_deref()
            .map(sparse_from_dense),
        lower_seed_decomposition_max_terms: 4,
        lower_seed_decomposition_status: lower_seed_decomposition.status,
        lower_seed_decomposition_term_count: lower_seed_decomposition.term_count,
        lower_seed_decomposition_terms_nonzero: lower_seed_decomposition.terms_nonzero,
        lower_seed_decomposition_error: lower_seed_decomposition.error,
        lower_seed_diamond_status: lower_seed_diamond.status,
        lower_seed_diamond_element_count: lower_seed_diamond.element_count,
        lower_seed_diamond_gv: lower_seed_diamond.gv,
        lower_seed_diamond_error: lower_seed_diamond.error,
    })
}

fn cygv_path_predecessor_stats(
    elements: &HashSet<Vec<i64>>,
    grading: &[i64],
    target: &[i64],
    selected_degrees: &HashSet<i128>,
) -> Result<CygvPathPredecessorStats, String> {
    let mut previous_window_element_count = 0usize;
    let mut predecessor_difference_count = 0usize;
    let mut improving_predecessor_difference_count = 0usize;
    let mut closest_distance = cygv_series_distance(target);
    let mut closest_predecessor = None;
    let mut closest_difference = None;
    let mut sorted_elements = elements.iter().collect::<Vec<_>>();
    sorted_elements.sort();
    for element in sorted_elements {
        let degree = curve_degree(element, grading)?;
        if !selected_degrees.contains(&degree) {
            continue;
        }
        previous_window_element_count += 1;
        let difference = checked_vector_difference(target, element)?;
        if !elements.contains(&difference) {
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
    Ok(CygvPathPredecessorStats {
        previous_window_element_count,
        predecessor_difference_count,
        improving_predecessor_difference_count,
        closest_distance,
        closest_predecessor,
        closest_difference,
    })
}

fn cygv_previous_window_degrees(
    degree_counts: &BTreeMap<i128, usize>,
    target_degree: i128,
    previous_level_count: usize,
) -> Vec<i128> {
    let mut selected = degree_counts
        .keys()
        .copied()
        .filter(|degree| *degree > 0 && *degree < target_degree)
        .rev()
        .take(previous_level_count)
        .collect::<Vec<_>>();
    selected.sort_unstable();
    selected
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
        let mut sorted_new_elements = new_elements.into_iter().collect::<Vec<_>>();
        sorted_new_elements.sort();
        if elements.len() + sorted_new_elements.len() > element_limit {
            for element in sorted_new_elements {
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
        for element in &sorted_new_elements {
            elements.insert(element.clone());
        }
        starting_elements = sorted_new_elements.into_iter().collect();
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
    let (target, generators) = active_support_generators(sample, context)?;

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

fn active_support_generators(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
) -> Result<(Vec<i64>, Vec<Vec<i64>>), String> {
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
    Ok((target, generators))
}

fn active_support_face_certificate_status(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
) -> String {
    let Ok((_, generators)) = active_support_generators(sample, context) else {
        return "active_support_face_certificate_error".to_string();
    };
    match certify_supporting_mori_face_by_exact_kernel(&generators, context.degree_bounded_rays) {
        Ok(Some(_)) => "active_support_certified_codimension_one_face".to_string(),
        Ok(None) => "active_support_not_certified_as_codimension_one_face".to_string(),
        Err(_) => "active_support_face_certificate_error".to_string(),
    }
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
    target_index_filter: Option<usize>,
    run_integer_diamonds: bool,
    run_active_support_generators: bool,
    support_overlap_min_for_run: Option<usize>,
    support_overlap_max_target_degree: Option<i128>,
    support_overlap_pair_reduce_for_run: bool,
    certify_origin_support_domains: bool,
    origin_support_certificate_limit: usize,
    certify_target_extremal_rays: bool,
    target_extremal_generator_limit: usize,
    target_extremal_max_degree: Option<i128>,
    measure_cygv_semigroups: bool,
    probe_cygv_path_history: bool,
    run_lower_seed_diamonds: bool,
    measure_cygv_degree_ladder: bool,
    cygv_degree_ladder_max_degree: Option<i128>,
    semigroup_measure_max_target_degree: Option<i128>,
    semigroup_measure_max_seed_count: Option<usize>,
    element_limit: usize,
) -> ContextReport {
    let mut semigroup_measurement_cache = HashMap::new();
    let mut semigroup_ladder_cache = HashMap::new();
    let mut targets = Vec::with_capacity(validated.stats.sample.len());
    for (idx, sample) in validated.stats.sample.iter().enumerate() {
        if !target_index_selected(idx, target_index_filter) {
            continue;
        }
        targets.push(report_target(
            idx,
            sample,
            validated,
            run_integer_diamonds,
            run_active_support_generators,
            support_overlap_min_for_run,
            support_overlap_max_target_degree,
            support_overlap_pair_reduce_for_run,
            certify_origin_support_domains,
            origin_support_certificate_limit,
            certify_target_extremal_rays,
            target_extremal_generator_limit,
            target_extremal_max_degree,
            measure_cygv_semigroups,
            probe_cygv_path_history,
            run_lower_seed_diamonds,
            measure_cygv_degree_ladder,
            cygv_degree_ladder_max_degree,
            semigroup_measure_max_target_degree,
            semigroup_measure_max_seed_count,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            element_limit,
        ));
    }
    let local_cygv_charge_signature_counts =
        local_cygv_charge_signature_counts(&validated.stats.sample, target_index_filter);
    let local_cygv_target_candidate_status_counts = local_cygv_target_candidate_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cygv_actual_call_readiness_counts = local_cygv_actual_call_readiness_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cygv_missing_source_input_counts = local_cygv_missing_source_input_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let origin_circuit_facet_context_status_counts =
        origin_circuit_facet_context_status_counts(&validated.stats.sample, target_index_filter);
    let active_support_status_counts = optional_status_counts(
        targets
            .iter()
            .map(|target| target.active_support_status.as_deref()),
        "not_run",
    );
    let active_support_face_certificate_status_counts =
        active_support_face_certificate_status_counts(
            &validated.stats.sample,
            &validated,
            target_index_filter,
        );
    let target_extremal_ray_certificate_status_counts = optional_status_counts(
        targets.iter().map(|target| {
            target
                .target_extremal_ray_certificate
                .as_ref()
                .map(|probe| probe.status.as_str())
        }),
        "not_run",
    );
    let origin_relation_support_face_certificate_status_counts = optional_status_counts(
        targets.iter().map(|target| {
            target
                .origin_relation_support_face_certificate_status
                .as_deref()
        }),
        "not_run",
    );
    let origin_shared_facet_face_certificate_status_counts = optional_status_counts(
        targets.iter().map(|target| {
            target
                .origin_shared_facet_face_certificate_status
                .as_deref()
        }),
        "not_run",
    );
    let origin_facet_union_face_certificate_status_counts = optional_status_counts(
        targets
            .iter()
            .map(|target| target.origin_facet_union_face_certificate_status.as_deref()),
        "not_run",
    );
    let local_cygv_q_matrix_orientation_status_counts =
        local_cygv_q_matrix_orientation_status_counts(
            targets
                .iter()
                .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
        );
    let local_cygv_q_matrix_layout_status_counts = local_cygv_q_matrix_layout_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cygv_origin_point_status_counts = local_cygv_origin_point_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cytools_origin_circuit_status_counts = local_cytools_origin_circuit_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cygv_grading_vector_status_counts = local_cygv_grading_vector_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    ContextReport {
        schema_version: context.schema_version,
        dimension: validated.dimension,
        ambient_rays: context.ambient_rays,
        basis_mori_ray_count: context.basis_mori_ray_count,
        degree_bound: validated.degree_bound,
        degree_bounded_ray_count: validated.degree_bounded_rays.len(),
        degree_bounded_mori_ray_context_count: validated
            .degree_bounded_ray_context
            .map(<[DegreeBoundedMoriRayContextSample]>::len),
        degree_bounded_mori_ray_context_status: degree_bounded_mori_ray_context_status(validated),
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
        local_cygv_charge_signature_counts,
        local_cygv_target_candidate_status_counts,
        local_cygv_actual_call_readiness_counts,
        local_cygv_missing_source_input_counts,
        origin_circuit_facet_context_status_counts,
        active_support_status_counts,
        active_support_face_certificate_status_counts,
        target_extremal_ray_certificate_status_counts,
        origin_relation_support_face_certificate_status_counts,
        origin_shared_facet_face_certificate_status_counts,
        origin_facet_union_face_certificate_status_counts,
        local_cygv_q_matrix_orientation_status_counts,
        local_cygv_q_matrix_layout_status_counts,
        local_cygv_origin_point_status_counts,
        local_cytools_origin_circuit_status_counts,
        local_cygv_grading_vector_status_counts,
        targets,
    }
}

fn degree_bounded_mori_ray_context_status(context: &ValidatedContext<'_>) -> String {
    match context.degree_bounded_ray_context {
        Some(samples) if samples.is_empty() => {
            "source_derived_degree_bounded_mori_ray_context_empty".to_string()
        }
        Some(_) => "source_derived_ambient_and_basis_degree_bounded_mori_ray_context".to_string(),
        None => "missing_degree_bounded_mori_ray_context".to_string(),
    }
}

fn local_cygv_charge_signature_counts(
    samples: &[MissingGvTargetSample],
    target_index_filter: Option<usize>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for (idx, sample) in samples.iter().enumerate() {
        if !target_index_selected(idx, target_index_filter) {
            continue;
        }
        let Some(support) = sample.origin_circuit_affine_support.as_ref() else {
            continue;
        };
        let key = local_charge_signature_key(&support.local_charge_basis);
        *counts.entry(key).or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_target_candidate_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        for candidate in &skeleton.orientation_candidates {
            *counts
                .entry(candidate.target_candidate_status.clone())
                .or_insert(0usize) += 1;
        }
    }
    counts
}

fn local_cygv_actual_call_readiness_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(local_cygv_actual_call_readiness(skeleton))
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_missing_source_input_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        for input in &skeleton.remaining_uncertified_inputs {
            *counts.entry(input.clone()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn optional_status_counts<'a>(
    statuses: impl IntoIterator<Item = Option<&'a str>>,
    missing_status: &str,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for status in statuses {
        *counts
            .entry(status.unwrap_or(missing_status).to_string())
            .or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_facet_context_status_counts(
    samples: &[MissingGvTargetSample],
    target_index_filter: Option<usize>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for (idx, sample) in samples.iter().enumerate() {
        if target_index_filter.is_some_and(|filter| filter != idx) {
            continue;
        }
        let status = origin_circuit_facet_context_status(sample);
        *counts.entry(status).or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_facet_context_status(sample: &MissingGvTargetSample) -> String {
    let Some(witness) = sample.origin_circuit_first_witness.as_ref() else {
        return "no_origin_circuit_witness".to_string();
    };
    if witness.first_facet.is_empty() || witness.second_facet.is_empty() {
        return "origin_circuit_missing_full_facet_context".to_string();
    }
    if witness.first_facet.len() != witness.first_facet_size
        || witness.second_facet.len() != witness.second_facet_size
    {
        return "origin_circuit_facet_context_size_mismatch".to_string();
    }

    let first_facet = witness.first_facet.iter().copied().collect::<HashSet<_>>();
    let second_facet = witness.second_facet.iter().copied().collect::<HashSet<_>>();
    let shared_two_simplex_is_common = witness
        .shared_two_simplex
        .iter()
        .all(|point| first_facet.contains(point) && second_facet.contains(point));
    let exclusive_points_are_exclusive = first_facet.contains(&witness.first_facet_exclusive_point)
        && !second_facet.contains(&witness.first_facet_exclusive_point)
        && second_facet.contains(&witness.second_facet_exclusive_point)
        && !first_facet.contains(&witness.second_facet_exclusive_point);

    if shared_two_simplex_is_common && exclusive_points_are_exclusive {
        "source_derived_full_facet_context".to_string()
    } else {
        "origin_circuit_facet_context_inconsistent".to_string()
    }
}

fn active_support_face_certificate_status_counts(
    samples: &[MissingGvTargetSample],
    context: &ValidatedContext<'_>,
    target_index_filter: Option<usize>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for (idx, sample) in samples.iter().enumerate() {
        if target_index_filter.is_some_and(|filter| filter != idx) {
            continue;
        }
        let status = active_support_face_certificate_status(sample, context);
        *counts.entry(status).or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_grading_vector_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(skeleton.local_grading_vector_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_q_matrix_orientation_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(skeleton.local_q_matrix_orientation_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_q_matrix_layout_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(skeleton.local_cygv_q_matrix_layout_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_origin_point_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(skeleton.local_cygv_origin_point_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cytools_origin_circuit_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(skeleton.local_cytools_origin_circuit_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_actual_call_readiness(skeleton: &LocalCygvInputSkeleton) -> String {
    if !skeleton.remaining_uncertified_inputs.is_empty() {
        return "blocked_missing_source_derived_inputs".to_string();
    }
    if skeleton.orientation_candidates.iter().any(|candidate| {
        candidate.target_candidate_status
            == "target_primitive_positive_supported_by_cygv_omega_bucket"
    }) {
        "ready_for_actual_cygv_call".to_string()
    } else {
        "blocked_no_supported_target_orientation".to_string()
    }
}

fn target_index_selected(index: usize, target_index_filter: Option<usize>) -> bool {
    target_index_filter.is_none_or(|selected| selected == index)
}

fn main() {
    let Some(context_path) = parse_arg_value::<PathBuf>("--context") else {
        eprintln!(
            "[ERROR] usage: mcallister_gv_context --context path [--target-index N] [--run-integer-diamonds] [--run-active-support-generators] [--run-support-overlap-generators N] [--pair-reduce-support-overlap-generators] [--support-overlap-max-target-degree N] [--certify-origin-support-domains] [--origin-support-certificate-limit N] [--certify-target-extremal-rays] [--target-extremal-generator-limit N] [--target-extremal-max-degree N] [--measure-cygv-semigroups] [--probe-cygv-path-history] [--run-lower-seed-diamonds] [--measure-cygv-degree-ladder --cygv-degree-ladder-max-degree N] [--semigroup-measure-max-target-degree N] [--semigroup-measure-max-seeds N] [--element-limit N] [--out path]\n       use --run-support-overlap-generators 0 to try all degree-bounded generators up to each target degree"
        );
        std::process::exit(2);
    };
    let target_index_filter = parse_arg_value::<usize>("--target-index");
    let run_integer_diamonds = parse_flag("--run-integer-diamonds");
    let run_active_support_generators = parse_flag("--run-active-support-generators");
    let support_overlap_min_for_run = parse_arg_value::<usize>("--run-support-overlap-generators");
    let support_overlap_pair_reduce_for_run =
        parse_flag("--pair-reduce-support-overlap-generators");
    let support_overlap_max_target_degree =
        parse_arg_value::<i128>("--support-overlap-max-target-degree");
    let certify_origin_support_domains = parse_flag("--certify-origin-support-domains");
    let origin_support_certificate_limit =
        parse_arg_value::<usize>("--origin-support-certificate-limit").unwrap_or(256);
    let certify_target_extremal_rays = parse_flag("--certify-target-extremal-rays");
    let target_extremal_generator_limit =
        parse_arg_value::<usize>("--target-extremal-generator-limit").unwrap_or(256);
    let target_extremal_max_degree = parse_arg_value::<i128>("--target-extremal-max-degree");
    let measure_cygv_semigroups = parse_flag("--measure-cygv-semigroups");
    let probe_cygv_path_history = parse_flag("--probe-cygv-path-history");
    let run_lower_seed_diamonds = parse_flag("--run-lower-seed-diamonds");
    let measure_cygv_degree_ladder = parse_flag("--measure-cygv-degree-ladder");
    let cygv_degree_ladder_max_degree = parse_arg_value::<i128>("--cygv-degree-ladder-max-degree");
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
        target_index_filter,
        run_integer_diamonds,
        run_active_support_generators,
        support_overlap_min_for_run,
        support_overlap_max_target_degree,
        support_overlap_pair_reduce_for_run,
        certify_origin_support_domains,
        origin_support_certificate_limit,
        certify_target_extremal_rays,
        target_extremal_generator_limit,
        target_extremal_max_degree,
        measure_cygv_semigroups,
        probe_cygv_path_history,
        run_lower_seed_diamonds,
        measure_cygv_degree_ladder,
        cygv_degree_ladder_max_degree,
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

    fn minimal_corrected_context(
        schema_version: u32,
        ray_context: Option<Vec<DegreeBoundedMoriRayContextSample>>,
    ) -> CorrectedChamberGvContext {
        CorrectedChamberGvContext {
            schema_version,
            ambient_rays: 2,
            toric_gv_missing_count: 0,
            remaining_gv_missing_count: 0,
            basis_mori_ray_count: Some(2),
            basis_mori_rays_for_missing_degree_bound: Some(2),
            basis_mori_rays_for_missing_degree_bounded: Some(vec![vec![1, 0], vec![0, 1]]),
            degree_bounded_mori_ray_context_for_missing: ray_context,
            gv_q_matrix_for_missing: Some(vec![vec![1, 0], vec![0, 1]]),
            grading_for_missing: Some(vec![1, 1]),
            corrected_kappa_basis_for_missing: Some(Vec::new()),
            missing_target_stats: Some(MissingGvTargetStats {
                target_count: 0,
                real_cone_decomposition_exact_kind_counts: HashMap::new(),
                sample: Vec::new(),
            }),
        }
    }

    #[test]
    fn validate_context_accepts_schema3_degree_bounded_ray_context() {
        let context = minimal_corrected_context(
            3,
            Some(vec![
                DegreeBoundedMoriRayContextSample {
                    degree: 1,
                    ambient_nonzero: vec![(5, 1)],
                    basis_nonzero: vec![(0, 1)],
                },
                DegreeBoundedMoriRayContextSample {
                    degree: 1,
                    ambient_nonzero: vec![(8, 1)],
                    basis_nonzero: vec![(1, 1)],
                },
            ]),
        );
        let validated = validate_context(&context).unwrap();

        assert_eq!(validated.degree_bounded_ray_context.unwrap().len(), 2);
        assert_eq!(
            degree_bounded_mori_ray_context_status(&validated),
            "source_derived_ambient_and_basis_degree_bounded_mori_ray_context"
        );
    }

    #[test]
    fn validate_context_requires_schema3_degree_bounded_ray_context() {
        let context = minimal_corrected_context(3, None);
        let err = match validate_context(&context) {
            Ok(_) => panic!("schema-3 context without ray context should be rejected"),
            Err(err) => err,
        };
        assert!(
            err.contains("schema-3 context is missing degree-bounded Mori ray context"),
            "{err}"
        );
    }

    #[test]
    fn validate_context_checks_degree_bounded_ray_context_degrees() {
        let context = minimal_corrected_context(
            3,
            Some(vec![DegreeBoundedMoriRayContextSample {
                degree: 2,
                ambient_nonzero: vec![(5, 1)],
                basis_nonzero: vec![(0, 1)],
            }]),
        );
        let err = match validate_context(&context) {
            Ok(_) => panic!("ray context with wrong degree should be rejected"),
            Err(err) => err,
        };
        assert!(err.contains("declares degree 2 but computes 1"), "{err}");
    }

    #[test]
    fn origin_circuit_ambient_generator_counts_filter_source_supports() {
        let stats = MissingGvTargetStats {
            target_count: 0,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: Vec::new(),
        };
        let grading = vec![1, 1];
        let q_matrix = vec![vec![1, 0], vec![0, 1]];
        let degree_bounded_rays = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        let ray_context = vec![
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, -1), (1, 1), (4, 1), (5, -1)],
                basis_nonzero: vec![(0, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 2,
                ambient_nonzero: vec![(2, 1), (3, -1)],
                basis_nonzero: vec![(1, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 2,
                ambient_nonzero: vec![(6, 1)],
                basis_nonzero: vec![(0, 1), (1, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 2,
                ambient_nonzero: vec![(8, 1)],
                basis_nonzero: vec![(0, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 3,
                ambient_nonzero: vec![(0, -1), (1, 1), (4, 1), (5, -1)],
                basis_nonzero: vec![(0, 2), (1, 1)],
            },
        ];
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 3,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: Some(&ray_context),
            intersection: Intersection::new(2),
            stats: &stats,
        };
        let sample = MissingGvTargetSample {
            degree: 2,
            generators_le_degree: 0,
            is_mori_generator: true,
            origin_circuit_pattern: None,
            origin_circuit_witness_count: Some(1),
            origin_circuit_first_witness: Some(OriginCircuitWitnessSample {
                first_facet_exclusive_point: 4,
                second_facet_exclusive_point: 5,
                shared_two_simplex: vec![1, 2, 3],
                first_facet: vec![1, 2, 3, 4, 6],
                second_facet: vec![1, 2, 3, 5, 7],
                first_facet_size: 5,
                second_facet_size: 5,
                sparse_relation: Vec::new(),
                relation_points: vec![
                    OriginCircuitRelationPointSample {
                        point_index: 0,
                        coefficient: -1,
                        coordinates: vec![0, 0],
                        face_dimension: None,
                    },
                    OriginCircuitRelationPointSample {
                        point_index: 1,
                        coefficient: -1,
                        coordinates: vec![1, 0],
                        face_dimension: None,
                    },
                    OriginCircuitRelationPointSample {
                        point_index: 4,
                        coefficient: 1,
                        coordinates: vec![0, 1],
                        face_dimension: None,
                    },
                    OriginCircuitRelationPointSample {
                        point_index: 5,
                        coefficient: 1,
                        coordinates: vec![1, 1],
                        face_dimension: None,
                    },
                ],
            }),
            origin_circuit_affine_support: None,
            cms_general_divisor_shape_candidates: None,
            cms_general_divisor_intersection_checks: None,
            branch_diagnostic: None,
            real_cone_decomposable_by_other_generators: false,
            real_cone_decomposition_active_generators: None,
            real_cone_decomposition_active_generator_basis_nonzero: None,
            real_cone_decomposition_exact_coefficients: None,
            real_cone_decomposition_exact_kind: None,
            ambient_nonzero: vec![(0, -1), (1, 1), (4, 1), (5, -1)],
            basis_nonzero: vec![(0, 1)],
        };

        let counts = origin_circuit_ambient_generator_counts(&sample, &context);
        assert_eq!(counts.relation_support, Some(1));
        assert_eq!(counts.shared_facet, Some(2));
        assert_eq!(counts.facet_union, Some(3));

        let certificates =
            origin_circuit_ambient_support_face_certificate_statuses(&sample, &context, 2);
        assert_eq!(
            certificates.relation_support.as_deref(),
            Some("origin_support_certified_codimension_one_face")
        );
        assert_eq!(
            certificates.shared_facet.as_deref(),
            Some("origin_support_lp_no_certificate_rank_2_dim_2")
        );
        assert_eq!(
            certificates.facet_union.as_deref(),
            Some("origin_support_skipped_generator_limit_2_actual_3")
        );
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
    fn bounded_seed_decomposition_finds_short_lower_seed_sums() {
        let seeds = vec![vec![0, 1], vec![1, 0], vec![2, 0]];

        assert_eq!(
            bounded_seed_decomposition(&[2, 1], &seeds, 2).unwrap(),
            Some(vec![vec![0, 1], vec![2, 0]])
        );
        assert_eq!(
            bounded_seed_decomposition(&[3, 1], &seeds, 3).unwrap(),
            Some(vec![vec![0, 1], vec![1, 0], vec![2, 0]])
        );
        assert_eq!(
            bounded_seed_decomposition(&[5, 1], &seeds, 3).unwrap(),
            None
        );
        assert_eq!(
            bounded_seed_decomposition(&[5, 1], &seeds, 4).unwrap(),
            Some(vec![vec![0, 1], vec![1, 0], vec![2, 0], vec![2, 0]])
        );
    }

    #[test]
    fn local_charge_signature_groups_permuted_rows() {
        let first = vec![vec![1, -2, -1, -1, 3]];
        let second = vec![vec![1, -1, -1, -2, 3]];

        assert_eq!(
            local_charge_row_permutation_signatures(&first),
            vec![vec![-2, -1, -1, 1, 3]]
        );
        assert_eq!(
            local_charge_signature_key(&first),
            local_charge_signature_key(&second)
        );
        assert_eq!(
            local_charge_multiplicities(&first[0])
                .into_iter()
                .map(|entry| (entry.charge, entry.count))
                .collect::<Vec<_>>(),
            vec![(-2, 1), (-1, 2), (1, 1), (3, 1)]
        );
    }

    #[test]
    fn compact_threefold_shape_is_not_enough_for_cygv_input() {
        let (status, missing) = cygv_compact_input_readiness(true);

        assert_eq!(status, "shape_only_missing_source_derived_cygv_inputs");
        assert_eq!(
            missing,
            vec![
                "local_semigroup_generators",
                "local_grading_vector",
                "local_q_matrix_orientation_and_phase",
                "local_intersection_tensor",
                "local_chamber_certificate",
                "target_class_to_local_semigroup_coordinate",
            ]
        );
    }

    #[test]
    fn origin_circuit_affine_support_reconstructs_full_coordinates_for_old_dumps() {
        let relation_points = vec![
            OriginCircuitRelationPointSample {
                point_index: 0,
                coefficient: -1,
                coordinates: vec![0, 0, 0, 0],
                face_dimension: Some(4),
            },
            OriginCircuitRelationPointSample {
                point_index: 46,
                coefficient: 2,
                coordinates: vec![1, 2, 0, 2],
                face_dimension: Some(1),
            },
            OriginCircuitRelationPointSample {
                point_index: 208,
                coefficient: 1,
                coordinates: vec![2, 2, 1, 2],
                face_dimension: Some(2),
            },
            OriginCircuitRelationPointSample {
                point_index: 211,
                coefficient: -3,
                coordinates: vec![2, 3, 1, 3],
                face_dimension: Some(2),
            },
            OriginCircuitRelationPointSample {
                point_index: 214,
                coefficient: 1,
                coordinates: vec![2, 3, 2, 3],
                face_dimension: Some(2),
            },
        ];
        let sample = MissingGvTargetSample {
            degree: 10,
            generators_le_degree: 720,
            is_mori_generator: true,
            origin_circuit_pattern: Some("-3:1,-1:1,1:2,2:1".to_string()),
            origin_circuit_witness_count: Some(1),
            origin_circuit_first_witness: Some(OriginCircuitWitnessSample {
                first_facet_exclusive_point: 208,
                second_facet_exclusive_point: 214,
                shared_two_simplex: vec![46, 55, 211],
                first_facet: Vec::new(),
                second_facet: Vec::new(),
                first_facet_size: 23,
                second_facet_size: 191,
                sparse_relation: relation_points
                    .iter()
                    .map(|point| (point.point_index, point.coefficient))
                    .collect(),
                relation_points: relation_points.clone(),
            }),
            origin_circuit_affine_support: Some(OriginCircuitAffineSupportSample {
                affine_rank: 3,
                coefficient_counts: BTreeMap::from([(-3, 1), (-1, 1), (1, 2), (2, 1)]),
                local_charge_basis: vec![vec![1, -2, -1, 3, -1]],
                local_coordinates: Vec::new(),
                local_coordinates_2d: None,
            }),
            cms_general_divisor_shape_candidates: None,
            cms_general_divisor_intersection_checks: None,
            branch_diagnostic: None,
            real_cone_decomposable_by_other_generators: false,
            real_cone_decomposition_active_generators: None,
            real_cone_decomposition_active_generator_basis_nonzero: None,
            real_cone_decomposition_exact_coefficients: None,
            real_cone_decomposition_exact_kind: None,
            ambient_nonzero: Vec::new(),
            basis_nonzero: Vec::new(),
        };

        let support = origin_circuit_affine_support_with_coordinates(&sample)
            .unwrap()
            .expect("support should be present");

        assert_eq!(support.local_coordinates.len(), 5);
        for point in &support.local_coordinates {
            assert_eq!(point.coordinates.len(), 3);
        }
        let coefficients = relation_points
            .iter()
            .map(|point| (point.point_index, point.coefficient))
            .collect::<HashMap<_, _>>();
        for coordinate_idx in 0..3 {
            let weighted_sum: i64 = support
                .local_coordinates
                .iter()
                .map(|point| coefficients[&point.point_index] * point.coordinates[coordinate_idx])
                .sum();
            assert_eq!(weighted_sum, 0);
        }

        let skeleton = local_cygv_input_skeleton(&sample, Some(&support))
            .unwrap()
            .expect("local skeleton should be present");
        assert_eq!(skeleton.support_point_indices, vec![0, 46, 208, 211, 214]);
        assert!(skeleton.support_contains_origin_point);
        assert_eq!(
            skeleton.local_cygv_origin_point_status,
            "local_gkz_relation_includes_origin_point_requires_phase_mapping"
        );
        assert_eq!(skeleton.origin_point_relation_coefficient, Some(-1));
        assert_eq!(
            skeleton.local_cytools_origin_circuit_status,
            "source_cytools_retains_negative_origin_coefficient"
        );
        assert_eq!(
            skeleton.local_q_matrix_rows,
            vec![vec![1], vec![-2], vec![-1], vec![3], vec![-1]]
        );
        assert_eq!(
            skeleton.target_relation_coefficients,
            Some(vec![-1, 2, 1, -3, 1])
        );
        assert_eq!(skeleton.target_relation_in_charge_basis, Some(vec![-1]));
        assert_eq!(
            skeleton.target_relation_status,
            "target_relation_integral_in_local_charge_basis"
        );
        assert_eq!(skeleton.local_q_matrix_orientation_candidate, Some(-1));
        assert_eq!(
            skeleton.local_q_matrix_orientation_status,
            "source_derived_target_positive_orientation"
        );
        assert_eq!(
            skeleton.local_cygv_q_matrix_rows_candidate,
            Some(vec![vec![-1], vec![2], vec![1], vec![-3], vec![1]])
        );
        assert_eq!(
            skeleton.local_cygv_wrapper_q_matrix_candidate,
            Some(vec![vec![-1, 2, 1, -3, 1]])
        );
        assert_eq!(
            skeleton.local_cygv_q_matrix_layout_status,
            "source_derived_oriented_q_matrix_layout"
        );
        assert_eq!(skeleton.local_grading_vector_candidate, Some(vec![1]));
        assert_eq!(
            skeleton.local_grading_vector_status,
            "source_derived_primitive_one_parameter_grading"
        );
        assert!(
            !skeleton
                .remaining_uncertified_inputs
                .iter()
                .any(|input| input == "local_grading_vector")
        );
        assert!(
            !skeleton
                .remaining_uncertified_inputs
                .iter()
                .any(|input| input == "local_q_matrix_orientation_and_phase")
        );
        assert!(
            skeleton
                .remaining_uncertified_inputs
                .iter()
                .any(|input| input == "local_q_matrix_phase")
        );
        assert_eq!(skeleton.orientation_candidates.len(), 2);
        let positive_target = &skeleton.orientation_candidates[0];
        assert_eq!(positive_target.overall_charge_basis_sign, -1);
        assert_eq!(positive_target.target_coordinate, Some(vec![1]));
        assert_eq!(
            positive_target.target_candidate_status,
            "target_primitive_positive_supported_by_cygv_omega_bucket"
        );
        assert_eq!(positive_target.target_coordinate_is_nonnegative, Some(true));
        assert_eq!(positive_target.target_coordinate_gcd, Some(1));
        assert_eq!(positive_target.target_coordinate_is_primitive, Some(true));
        assert_eq!(positive_target.target_primitive_direction, Some(vec![1]));
        assert_eq!(
            positive_target.positive_unit_generator_negative_intersections,
            Some(2)
        );
        assert_eq!(
            positive_target
                .positive_unit_generator_omega_bucket
                .as_deref(),
            Some("neg2")
        );
        let original_orientation = &skeleton.orientation_candidates[1];
        assert_eq!(original_orientation.overall_charge_basis_sign, 1);
        assert_eq!(original_orientation.target_coordinate, Some(vec![-1]));
        assert_eq!(
            original_orientation.target_candidate_status,
            "target_not_in_nonnegative_local_semigroup"
        );
        assert_eq!(original_orientation.target_coordinate_gcd, Some(1));
        assert_eq!(
            original_orientation.target_primitive_direction,
            Some(vec![-1])
        );
        assert_eq!(
            original_orientation.positive_unit_generator_negative_intersections,
            Some(3)
        );
        assert_eq!(
            original_orientation
                .positive_unit_generator_omega_bucket
                .as_deref(),
            Some("ignored_gt2")
        );
    }

    #[test]
    fn local_cygv_orientation_status_marks_positive_omega_ignored_candidates() {
        let candidates = local_cygv_orientation_candidates(
            &[vec![2], vec![1], vec![2], vec![-1], vec![-2], vec![-2]],
            Some(&[-1]),
        );

        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].overall_charge_basis_sign, -1);
        assert_eq!(candidates[0].target_coordinate, Some(vec![1]));
        assert_eq!(
            candidates[0]
                .positive_unit_generator_omega_bucket
                .as_deref(),
            Some("ignored_gt2")
        );
        assert_eq!(
            candidates[0].target_candidate_status,
            "target_positive_but_ignored_by_cygv_omega_bucket"
        );
        assert_eq!(candidates[1].overall_charge_basis_sign, 1);
        assert_eq!(
            candidates[1].target_candidate_status,
            "target_not_in_nonnegative_local_semigroup"
        );
    }

    #[test]
    fn local_cygv_target_status_counts_aggregate_orientation_candidates() {
        let compact_threefold_like = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: local_cygv_orientation_candidates(
                &[vec![1], vec![-2], vec![-1], vec![3], vec![-1]],
                Some(&[-1]),
            ),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };
        let fourfold_like = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: local_cygv_orientation_candidates(
                &[vec![2], vec![1], vec![2], vec![-1], vec![-2], vec![-2]],
                Some(&[-1]),
            ),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };

        let counts =
            local_cygv_target_candidate_status_counts([&compact_threefold_like, &fourfold_like]);

        assert_eq!(
            counts
                .get("target_primitive_positive_supported_by_cygv_omega_bucket")
                .copied(),
            Some(1)
        );
        assert_eq!(
            counts
                .get("target_positive_but_ignored_by_cygv_omega_bucket")
                .copied(),
            Some(1)
        );
        assert_eq!(
            counts
                .get("target_not_in_nonnegative_local_semigroup")
                .copied(),
            Some(2)
        );
    }

    #[test]
    fn local_cygv_actual_call_readiness_requires_all_source_inputs() {
        let supported_candidate = local_cygv_orientation_candidates(
            &[vec![1], vec![-2], vec![-1], vec![3], vec![-1]],
            Some(&[-1]),
        );
        let blocked_missing_inputs = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: supported_candidate.clone(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: vec!["local_intersection_tensor".to_string()],
        };
        let ready = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: supported_candidate,
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };
        let blocked_orientation = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: local_cygv_orientation_candidates(
                &[vec![2], vec![1], vec![2], vec![-1], vec![-2], vec![-2]],
                Some(&[1]),
            ),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };

        assert_eq!(
            local_cygv_actual_call_readiness(&blocked_missing_inputs),
            "blocked_missing_source_derived_inputs"
        );
        assert_eq!(
            local_cygv_actual_call_readiness(&blocked_orientation),
            "blocked_no_supported_target_orientation"
        );
        assert_eq!(
            local_cygv_actual_call_readiness(&ready),
            "ready_for_actual_cygv_call"
        );

        let counts = local_cygv_actual_call_readiness_counts([
            &blocked_missing_inputs,
            &ready,
            &blocked_orientation,
        ]);
        assert_eq!(
            counts.get("blocked_missing_source_derived_inputs").copied(),
            Some(1)
        );
        assert_eq!(
            counts
                .get("blocked_no_supported_target_orientation")
                .copied(),
            Some(1)
        );
        assert_eq!(counts.get("ready_for_actual_cygv_call").copied(), Some(1));
    }

    #[test]
    fn local_cygv_missing_source_input_counts_aggregate_uncertified_inputs() {
        let first = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: vec![
                "local_semigroup_generators".to_string(),
                "local_intersection_tensor".to_string(),
            ],
        };
        let second = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: vec![
                "local_semigroup_generators".to_string(),
                "local_chamber_certificate".to_string(),
            ],
        };

        let counts = local_cygv_missing_source_input_counts([&first, &second]);

        assert_eq!(counts.get("local_semigroup_generators").copied(), Some(2));
        assert_eq!(counts.get("local_intersection_tensor").copied(), Some(1));
        assert_eq!(counts.get("local_chamber_certificate").copied(), Some(1));
    }

    #[test]
    fn optional_status_counts_aggregates_missing_and_present_statuses() {
        let computed = "computed_active_support_generators".to_string();
        let hkty_error = "hkty_error".to_string();
        let counts = optional_status_counts(
            [
                Some(computed.as_str()),
                None,
                Some(hkty_error.as_str()),
                Some(computed.as_str()),
            ],
            "not_run",
        );

        assert_eq!(
            counts.get("computed_active_support_generators").copied(),
            Some(2)
        );
        assert_eq!(counts.get("hkty_error").copied(), Some(1));
        assert_eq!(counts.get("not_run").copied(), Some(1));
    }

    #[test]
    fn local_cygv_grading_vector_status_counts_aggregate_skeletons() {
        let first = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: Some(vec![1]),
            local_grading_vector_status: "source_derived_primitive_one_parameter_grading"
                .to_string(),
            remaining_uncertified_inputs: Vec::new(),
        };
        let second = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: "local_grading_blocked_no_positive_primitive_target"
                .to_string(),
            remaining_uncertified_inputs: Vec::new(),
        };
        let third = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_grading_vector_candidate: Some(vec![1]),
            local_grading_vector_status: "source_derived_primitive_one_parameter_grading"
                .to_string(),
            remaining_uncertified_inputs: Vec::new(),
        };

        let counts = local_cygv_grading_vector_status_counts([&first, &second, &third]);

        assert_eq!(
            counts
                .get("source_derived_primitive_one_parameter_grading")
                .copied(),
            Some(2)
        );
        assert_eq!(
            counts
                .get("local_grading_blocked_no_positive_primitive_target")
                .copied(),
            Some(1)
        );
    }

    #[test]
    fn local_cygv_q_matrix_orientation_status_counts_aggregate_skeletons() {
        let first = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: true,
            local_cygv_origin_point_status:
                "local_gkz_relation_includes_origin_point_requires_phase_mapping".to_string(),
            origin_point_relation_coefficient: Some(-1),
            local_cytools_origin_circuit_status:
                "source_cytools_retains_negative_origin_coefficient".to_string(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: Some(-1),
            local_q_matrix_orientation_status: "source_derived_target_positive_orientation"
                .to_string(),
            local_cygv_q_matrix_rows_candidate: Some(vec![vec![-1], vec![2]]),
            local_cygv_wrapper_q_matrix_candidate: Some(vec![vec![-1, 2]]),
            local_cygv_q_matrix_layout_status: "source_derived_oriented_q_matrix_layout"
                .to_string(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };
        let second = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: false,
            local_cygv_origin_point_status: "local_support_has_no_origin_point".to_string(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: "not_cytools_origin_circuit_support".to_string(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status:
                "local_q_orientation_blocked_no_positive_primitive_target".to_string(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: "local_q_matrix_layout_blocked_no_orientation"
                .to_string(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };
        let third = LocalCygvInputSkeleton {
            support_point_indices: Vec::new(),
            support_contains_origin_point: true,
            local_cygv_origin_point_status:
                "local_gkz_relation_includes_origin_point_requires_phase_mapping".to_string(),
            origin_point_relation_coefficient: Some(-1),
            local_cytools_origin_circuit_status:
                "source_cytools_retains_negative_origin_coefficient".to_string(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: Some(-1),
            local_q_matrix_orientation_status: "source_derived_target_positive_orientation"
                .to_string(),
            local_cygv_q_matrix_rows_candidate: Some(vec![vec![-1], vec![1]]),
            local_cygv_wrapper_q_matrix_candidate: Some(vec![vec![-1, 1]]),
            local_cygv_q_matrix_layout_status: "source_derived_oriented_q_matrix_layout"
                .to_string(),
            local_grading_vector_candidate: None,
            local_grading_vector_status: String::new(),
            remaining_uncertified_inputs: Vec::new(),
        };

        let counts = local_cygv_q_matrix_orientation_status_counts([&first, &second, &third]);

        assert_eq!(
            counts
                .get("source_derived_target_positive_orientation")
                .copied(),
            Some(2)
        );
        assert_eq!(
            counts
                .get("local_q_orientation_blocked_no_positive_primitive_target")
                .copied(),
            Some(1)
        );

        let layout_counts = local_cygv_q_matrix_layout_status_counts([&first, &second, &third]);
        assert_eq!(
            layout_counts
                .get("source_derived_oriented_q_matrix_layout")
                .copied(),
            Some(2)
        );
        assert_eq!(
            layout_counts
                .get("local_q_matrix_layout_blocked_no_orientation")
                .copied(),
            Some(1)
        );

        let origin_counts = local_cygv_origin_point_status_counts([&first, &second, &third]);
        assert_eq!(
            origin_counts
                .get("local_gkz_relation_includes_origin_point_requires_phase_mapping")
                .copied(),
            Some(2)
        );
        assert_eq!(
            origin_counts
                .get("local_support_has_no_origin_point")
                .copied(),
            Some(1)
        );

        let cytools_origin_counts =
            local_cytools_origin_circuit_status_counts([&first, &second, &third]);
        assert_eq!(
            cytools_origin_counts
                .get("source_cytools_retains_negative_origin_coefficient")
                .copied(),
            Some(2)
        );
        assert_eq!(
            cytools_origin_counts
                .get("not_cytools_origin_circuit_support")
                .copied(),
            Some(1)
        );
    }

    #[test]
    fn lower_seed_diamond_probe_is_opt_in_and_respects_element_limit() {
        let stats = MissingGvTargetStats {
            target_count: 0,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: Vec::new(),
        };
        let grading = vec![1, 1];
        let q_matrix = vec![vec![1, 0], vec![0, 1]];
        let degree_bounded_rays = Vec::new();
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 2,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };
        let decomposition = LowerSeedDecompositionProbe {
            status: "found_lower_seed_decomposition".to_string(),
            term_count: Some(2),
            terms_nonzero: Some(vec![vec![(0, 1)], vec![(1, 1)]]),
            terms: Some(vec![vec![1, 0], vec![0, 1]]),
            error: None,
        };

        let disabled = lower_seed_diamond_probe(&[1, 1], &decomposition, &context, false, 2);
        assert_eq!(disabled.status, None);
        assert_eq!(disabled.element_count, None);

        let limited = lower_seed_diamond_probe(&[1, 1], &decomposition, &context, true, 2);
        assert_eq!(limited.status.as_deref(), Some("skipped_element_limit_2"));
        assert_eq!(limited.element_count, Some(4));
    }

    #[test]
    fn cygv_degree_ladder_uses_actual_semigroup_counts() {
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 3,
                generators_le_degree: 3,
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
                ambient_nonzero: vec![(0, 2), (1, 1)],
                basis_nonzero: vec![(0, 2), (1, 1)],
            }],
        };
        let grading = vec![1, 1];
        let q_matrix = vec![vec![1, 0], vec![0, 1]];
        let degree_bounded_rays = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 3,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };
        let mut cache = HashMap::new();

        let steps =
            measure_cygv_semigroup_degree_ladder(&stats.sample[0], &context, 2, None, &mut cache);

        assert_eq!(steps.len(), 2);
        assert_eq!(steps[0].degree, 1);
        assert_eq!(steps[0].effective_seed_count, 2);
        assert_eq!(steps[0].reduced_seed_count, Some(2));
        assert_eq!(steps[0].status, "measured_cygv_semigroup");
        assert_eq!(steps[0].element_count, Some(3));
        assert_eq!(steps[1].degree, 2);
        assert_eq!(steps[1].effective_seed_count, 3);
        assert_eq!(steps[1].reduced_seed_count, Some(2));
        assert_eq!(steps[1].status, "measured_cygv_semigroup");
        assert_eq!(steps[1].element_count, Some(6));
    }

    #[test]
    fn curve_degree_rejects_dimension_mismatch() {
        assert_eq!(curve_degree(&[2, -1], &[3, 5]).unwrap(), 1);
        assert!(curve_degree(&[1], &[1, 2]).is_err());
    }

    #[test]
    fn target_index_filter_selects_only_requested_target() {
        assert!(target_index_selected(3, None));
        assert!(target_index_selected(3, Some(3)));
        assert!(!target_index_selected(3, Some(4)));
    }

    #[test]
    fn status_error_fragment_is_stable_for_report_keys() {
        assert_eq!(
            status_error_fragment("supporting Mori face normal LP failed: infeasible?"),
            "supporting_mori_face_normal_lp_failed_infeasible"
        );
    }

    #[test]
    fn target_extremal_ray_probe_uses_exact_separator() {
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 1,
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
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };

        let probe = target_extremal_ray_certificate_probe(
            &stats.sample[0],
            &[1, 0],
            &context,
            true,
            3,
            None,
        )
        .expect("enabled probe should report a status");
        assert_eq!(probe.status, "certified_exact_extremal_ray");
        assert_eq!(probe.same_ray_generator_count, Some(1));
        assert_eq!(probe.zero_other_generator_count, Some(1));
        assert_eq!(probe.positive_other_generator_count, Some(1));
        assert_eq!(probe.separator_normal_nonzero, Some(vec![(0, -1), (1, 1)]));

        let limited = target_extremal_ray_certificate_probe(
            &stats.sample[0],
            &[1, 0],
            &context,
            true,
            2,
            None,
        )
        .expect("enabled probe should report a limit status");
        assert_eq!(limited.status, "skipped_generator_limit_2_actual_3");
    }

    #[test]
    fn target_extremal_ray_probe_short_circuits_exact_decomposition() {
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 2,
                generators_le_degree: 3,
                is_mori_generator: true,
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
                    vec![(0, 1)],
                    vec![(1, 1)],
                ]),
                real_cone_decomposition_exact_coefficients: Some(vec![
                    "1".to_string(),
                    "1".to_string(),
                ]),
                real_cone_decomposition_exact_kind: Some("integer_semigroup".to_string()),
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
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };

        let probe = target_extremal_ray_certificate_probe(
            &stats.sample[0],
            &[1, 1],
            &context,
            true,
            0,
            None,
        )
        .expect("enabled probe should report the exact decomposition");
        assert_eq!(
            probe.status,
            "not_extremal_by_exact_integer_semigroup_decomposition"
        );
        assert_eq!(probe.decomposition_active_generator_count, Some(2));
        assert_eq!(
            probe.decomposition_exact_coefficients,
            Some(vec!["1".to_string(), "1".to_string()])
        );
        assert_eq!(probe.separator_normal_nonzero, None);
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
    fn bounded_cygv_closure_truncates_new_elements_deterministically() {
        let seeds = vec![vec![1, 0], vec![0, 1]];
        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 2, 4).unwrap();

        assert_eq!(closure.status, "exceeded_element_limit_4");
        assert!(!closure.completed);
        assert_eq!(
            closure.elements,
            HashSet::from([vec![0, 0], vec![0, 1], vec![1, 0], vec![0, 2]])
        );
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
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };

        let target = vec![1, 1];
        let probe =
            cygv_path_history_probe_inner(&stats.sample[0], &context, &target, false, 16, None)
                .unwrap();
        assert_eq!(probe.status, "completed_bounded_closure");
        assert_eq!(probe.seed_count, Some(3));
        assert_eq!(probe.reduced_seed_count, Some(2));
        assert_eq!(probe.target_in_closure, Some(true));
        assert_eq!(probe.previous_level_count, 2);
        assert_eq!(probe.previous_window_degrees, vec![1]);
        assert_eq!(probe.previous_window_degree_count, Some(1));
        assert_eq!(probe.previous_window_element_count, Some(2));
        assert!(probe.predecessor_counts_complete);
        assert_eq!(probe.predecessor_difference_count, Some(2));
        assert_eq!(probe.improving_predecessor_difference_count, Some(1));
        assert_eq!(probe.closest_series_distance.as_deref(), Some("1.000000"));
        assert_eq!(probe.closest_series_predecessor_nonzero, Some(vec![(1, 1)]));
        assert_eq!(probe.closest_series_difference_nonzero, Some(vec![(0, 1)]));
    }

    #[test]
    fn path_history_probe_reports_previous_degrees_for_incomplete_closure() {
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
        let degree_bounded_rays = vec![vec![1, 0], vec![0, 1]];
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 2,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };

        let target = vec![1, 1];
        let probe =
            cygv_path_history_probe_inner(&stats.sample[0], &context, &target, false, 2, None)
                .unwrap();

        assert_eq!(probe.status, "exceeded_element_limit_initial_2");
        assert_eq!(probe.seed_count, Some(2));
        assert_eq!(probe.reduced_seed_count, Some(2));
        assert_eq!(probe.previous_window_degrees, vec![1]);
        assert_eq!(probe.previous_window_degree_count, Some(1));
        assert!(!probe.predecessor_counts_complete);
        assert_eq!(probe.predecessor_difference_count, Some(2));
        assert_eq!(probe.closest_series_distance.as_deref(), Some("1.000000"));
    }

    #[test]
    fn path_history_probe_respects_seed_limit_before_closure() {
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 2,
                generators_le_degree: 3,
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
            degree_bounded_ray_context: None,
            intersection: Intersection::new(2),
            stats: &stats,
        };

        let target = vec![1, 1];
        let probe =
            cygv_path_history_probe_inner(&stats.sample[0], &context, &target, false, 16, Some(2))
                .unwrap();

        assert_eq!(probe.status, "skipped_seed_limit");
        assert_eq!(probe.seed_count, Some(3));
        assert_eq!(probe.reduced_seed_count, None);
        assert_eq!(probe.lower_seed_decomposition_status, "skipped_seed_limit");
        assert_eq!(
            probe.lower_seed_diamond_status.as_deref(),
            Some("skipped_seed_limit")
        );
        assert_eq!(probe.closure_element_count, None);
        assert_eq!(probe.previous_window_element_count, None);
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
            first_facet: vec![2, 3, 5, 7],
            second_facet: vec![2, 3, 5, 11],
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
                    point_index: 11,
                    coefficient: 1,
                    coordinates: vec![-1, 0, 0, 0],
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
            local_coordinates: Vec::new(),
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
            degree_bounded_ray_context: None,
            intersection,
            stats: &stats,
        };

        let mut semigroup_measurement_cache = HashMap::new();
        let mut semigroup_ladder_cache = HashMap::new();
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
            256,
            false,
            256,
            None,
            false,
            false,
            false,
            false,
            None,
            None,
            None,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
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
            origin_circuit_facet_context_status(&stats.sample[0]),
            "source_derived_full_facet_context"
        );
        assert_eq!(
            origin_circuit_facet_context_status_counts(&stats.sample, None)
                .get("source_derived_full_facet_context")
                .copied(),
            Some(1)
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
        assert_eq!(
            shape.charge_row_permutation_signatures,
            vec![vec![-2, 1, 1]]
        );
        assert!(!shape.is_compact_threefold_hypersurface_shape);
        assert_eq!(
            shape.cygv_compact_input_status,
            "not_compact_threefold_hypersurface_shape"
        );
        assert!(shape.cygv_compact_input_missing.is_empty());
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
            degree_bounded_ray_context: None,
            intersection,
            stats: &stats,
        };
        let mut semigroup_measurement_cache = HashMap::new();
        let mut semigroup_ladder_cache = HashMap::new();

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
            256,
            false,
            256,
            None,
            false,
            false,
            false,
            false,
            None,
            None,
            None,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            256,
        );

        assert_eq!(
            report.support_overlap_run_status.as_deref(),
            Some("skipped_target_degree_limit")
        );
        assert_eq!(report.support_overlap_run_generator_count, None);

        let path_report = report_target(
            0,
            &stats.sample[0],
            &context,
            false,
            false,
            None,
            None,
            false,
            false,
            256,
            false,
            256,
            None,
            false,
            true,
            false,
            false,
            None,
            None,
            None,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            256,
        );

        assert_eq!(path_report.cygv_semigroup_measure_status, None);
        let path = path_report
            .cygv_path_history_probe
            .as_ref()
            .expect("path-history flag should populate the bounded probe");
        assert_eq!(path.status, "completed_bounded_closure");
        assert_eq!(path.target_in_closure, Some(true));
    }
}
