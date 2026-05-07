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
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
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

const CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT: usize = 32;
const ORIGIN_CIRCUIT_WITNESS_DOMAIN_UNRESOLVED_SAMPLE_LIMIT: usize = 64;
const ORIGIN_CIRCUIT_WITNESS_DOMAIN_OCCURRENCE_SAMPLE_LIMIT: usize = 8;

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
    covered_toric_gv_context_for_missing: Option<Vec<CoveredToricGvContextSample>>,
    #[serde(default)]
    uncovered_source_ray_toric_diagnostic_sample: Option<Vec<ToricGvDiagnosticContextSample>>,
    #[serde(default)]
    degree_bounded_toric_gv_diagnostic_context_for_missing:
        Option<Vec<ToricGvDiagnosticContextSample>>,
    uncovered_source_ray_stats_for_missing: Option<MissingGvTargetStats>,
    #[serde(default)]
    shared_facet_unresolved_source_ray_stats_for_missing: Option<MissingGvTargetStats>,
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

#[derive(Clone, Debug, Deserialize, Serialize)]
struct CoveredToricGvContextSample {
    degree: i128,
    gv: String,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ToricGvDiagnosticContextSample {
    degree: i128,
    gv: String,
    source_bucket: String,
    source_summary: String,
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
    #[serde(default)]
    origin_circuit_witnesses: Option<Vec<OriginCircuitWitnessSample>>,
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
    #[serde(default)]
    solution_basis_nonzero: Option<Vec<(usize, String)>>,
    #[serde(default)]
    solution_ambient_basis_nonzero: Option<Vec<(usize, String)>>,
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
    covered_toric_gv_context_count: Option<usize>,
    degree_bounded_toric_gv_diagnostic_context_count: Option<usize>,
    q_rows: usize,
    q_cols: usize,
    kappa_nonzero_entries: usize,
    toric_gv_missing_count: usize,
    remaining_gv_missing_count: usize,
    missing_target_count: usize,
    exact_kind_counts: HashMap<String, usize>,
    target_status_counts: BTreeMap<String, usize>,
    active_decomposition_generator_source_status_counts: BTreeMap<String, usize>,
    active_decomposition_unresolved_source_leaf_sample: Vec<ActiveDecompositionSourceLeafSummary>,
    active_decomposition_source_leaf_unit_phase_probe_status_counts: BTreeMap<String, usize>,
    active_decomposition_source_leaf_origin_omitted_unit_phase_probe_status_counts:
        BTreeMap<String, usize>,
    active_decomposition_source_leaf_cms_solution_status_counts: BTreeMap<String, usize>,
    active_decomposition_source_leaf_cms_candidate_status_counts: BTreeMap<String, usize>,
    local_cygv_charge_signature_counts: BTreeMap<String, usize>,
    local_cygv_target_candidate_status_counts: BTreeMap<String, usize>,
    local_cygv_actual_call_readiness_counts: BTreeMap<String, usize>,
    local_cygv_missing_source_input_counts: BTreeMap<String, usize>,
    cms_general_divisor_candidate_status_counts: BTreeMap<String, usize>,
    cms_general_divisor_intersection_check_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_relation_status_counts: BTreeMap<String, usize>,
    uncovered_source_ray_origin_circuit_witness_relation_status_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_target_count: Option<usize>,
    shared_facet_unresolved_source_ray_origin_circuit_pattern_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_origin_circuit_witness_relation_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_local_cygv_charge_signature_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_local_cygv_actual_call_readiness_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_local_cygv_missing_source_input_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_general_divisor_candidate_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_general_divisor_intersection_check_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_solution_status_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_intersection_tensor_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_primitive_probe_status_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_primitive_probe_gv_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_cms_solution_summary_error_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_unit_phase_probe_status_counts: BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_origin_omitted_unit_phase_probe_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_unit_effective_tensor_requirement_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_origin_omitted_effective_tensor_requirement_status_counts:
        BTreeMap<String, usize>,
    shared_facet_unresolved_source_ray_unit_phase_probe_sample:
        Vec<LocalCygvTargetUnitPhaseProbeSummary>,
    shared_facet_unresolved_source_ray_sample: Vec<TargetReport>,
    origin_circuit_facet_context_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_relation_support_face_profile_counts: BTreeMap<String, usize>,
    origin_circuit_witness_shared_facet_face_profile_counts: BTreeMap<String, usize>,
    origin_circuit_witness_facet_union_face_profile_counts: BTreeMap<String, usize>,
    origin_circuit_witness_relation_support_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_shared_facet_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_facet_union_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_domain_sample: Vec<OriginCircuitWitnessDomainSummary>,
    origin_circuit_witness_domain_unresolved_generator_unique_count: usize,
    origin_circuit_witness_domain_unresolved_generator_occurrence_count: usize,
    origin_circuit_witness_domain_unresolved_generator_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_domain_unresolved_generator_degree_counts: BTreeMap<i128, usize>,
    origin_circuit_witness_domain_unresolved_generator_sample:
        Vec<OriginCircuitWitnessDomainUnresolvedGeneratorSummary>,
    origin_circuit_witness_shared_facet_unresolved_generator_unique_count: usize,
    origin_circuit_witness_shared_facet_unresolved_generator_occurrence_count: usize,
    origin_circuit_witness_shared_facet_unresolved_generator_status_counts: BTreeMap<String, usize>,
    origin_circuit_witness_shared_facet_unresolved_generator_degree_counts: BTreeMap<i128, usize>,
    origin_circuit_witness_shared_facet_unresolved_generator_sample:
        Vec<OriginCircuitWitnessDomainUnresolvedGeneratorSummary>,
    active_support_status_counts: BTreeMap<String, usize>,
    active_support_face_certificate_status_counts: BTreeMap<String, usize>,
    target_extremal_ray_certificate_status_counts: BTreeMap<String, usize>,
    origin_relation_support_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_shared_facet_face_certificate_status_counts: BTreeMap<String, usize>,
    origin_facet_union_face_certificate_status_counts: BTreeMap<String, usize>,
    cygv_path_history_status_counts: BTreeMap<String, usize>,
    cygv_lower_seed_decomposition_status_counts: BTreeMap<String, usize>,
    cygv_lower_seed_diamond_status_counts: BTreeMap<String, usize>,
    path_support_uncovered_source_ray_unique_count: usize,
    path_support_uncovered_source_ray_occurrence_count: usize,
    path_support_uncovered_source_ray_degree_counts: BTreeMap<i128, usize>,
    path_support_uncovered_source_ray_lookup_status_counts: BTreeMap<String, usize>,
    path_support_uncovered_source_ray_path_support_gv_counts: BTreeMap<String, usize>,
    path_support_uncovered_source_ray_local_cygv_primitive_probe_status_counts:
        BTreeMap<String, usize>,
    path_support_uncovered_source_ray_sample: Vec<CygvPathSupportSourceRaySummary>,
    local_cygv_q_matrix_orientation_status_counts: BTreeMap<String, usize>,
    local_cygv_q_matrix_layout_status_counts: BTreeMap<String, usize>,
    local_cygv_q_matrix_phase_status_counts: BTreeMap<String, usize>,
    local_cygv_origin_point_status_counts: BTreeMap<String, usize>,
    local_cygv_origin_omitted_compact_shape_status_counts: BTreeMap<String, usize>,
    local_cygv_target_unit_phase_probe_status_counts: BTreeMap<String, usize>,
    local_cygv_target_origin_omitted_unit_phase_probe_status_counts: BTreeMap<String, usize>,
    local_cygv_target_unit_effective_tensor_requirement_status_counts: BTreeMap<String, usize>,
    local_cygv_target_origin_omitted_effective_tensor_requirement_status_counts:
        BTreeMap<String, usize>,
    local_cygv_target_unit_formula_sum_probe_status_counts: BTreeMap<String, usize>,
    local_cygv_target_origin_omitted_unit_formula_sum_probe_status_counts: BTreeMap<String, usize>,
    local_cygv_target_unit_formula_sum_effective_tensor_requirement_status_counts:
        BTreeMap<String, usize>,
    local_cygv_target_origin_omitted_unit_formula_sum_effective_tensor_requirement_status_counts:
        BTreeMap<String, usize>,
    local_cygv_target_unit_phase_probe_sample: Vec<LocalCygvTargetUnitPhaseProbeSummary>,
    local_cygv_target_integer_tensor_scan_status_counts: BTreeMap<String, usize>,
    local_cygv_target_integer_tensor_scan_sample: Vec<LocalCygvIntegerTensorScanSummary>,
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
    origin_circuit_witnesses: Option<Vec<OriginCircuitWitnessSample>>,
    origin_circuit_affine_support: Option<OriginCircuitAffineSupportSample>,
    local_cygv_hypersurface_shape: Option<LocalCygvHypersurfaceShape>,
    local_cygv_input_skeleton: Option<LocalCygvInputSkeleton>,
    cms_general_divisor_shape_candidates: Option<Vec<CmsGeneralDivisorShapeCandidate>>,
    cms_general_divisor_intersection_checks: Option<Vec<CmsGeneralDivisorIntersectionCheck>>,
    cms_general_divisor_solution_summaries: Vec<CmsGeneralDivisorSolutionSummary>,
    cms_general_divisor_solution_summary_error: Option<String>,
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
    active_support_face_certificate_status: Option<String>,
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
    local_semigroup_generators_candidate: Option<Vec<Vec<i64>>>,
    local_semigroup_generator_status: String,
    orientation_candidates: Vec<LocalCygvOrientationCandidate>,
    local_q_matrix_orientation_candidate: Option<i64>,
    local_q_matrix_orientation_status: String,
    local_cygv_q_matrix_rows_candidate: Option<Vec<Vec<i64>>>,
    local_cygv_wrapper_q_matrix_candidate: Option<Vec<Vec<i64>>>,
    local_cygv_q_matrix_layout_status: String,
    local_cygv_phase_q_matrix_candidate: Option<Vec<Vec<i64>>>,
    local_cygv_q_matrix_phase_status: String,
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
    closure_generation_counts: Vec<CygvClosureGenerationCount>,
    target_in_closure: Option<bool>,
    previous_level_count: usize,
    previous_window_degrees: Vec<i128>,
    previous_window_degree_count: Option<usize>,
    previous_window_element_count: Option<usize>,
    predecessor_counts_complete: bool,
    predecessor_difference_count: Option<usize>,
    improving_predecessor_difference_count: Option<usize>,
    predecessor_toric_coverage_counts: BTreeMap<String, usize>,
    predecessor_known_qn_history_counts: BTreeMap<String, usize>,
    closest_series_distance: Option<String>,
    closest_series_predecessor_nonzero: Option<Vec<(usize, i64)>>,
    closest_series_difference_nonzero: Option<Vec<(usize, i64)>>,
    closest_known_qn_predecessor: Option<CygvClosestKnownQnPredecessor>,
    closest_known_qn_residual_predecessor: Option<CygvClosestKnownQnPredecessor>,
    predecessor_candidate_sample_limit: usize,
    predecessor_candidate_sample: Vec<CygvPathPredecessorCandidate>,
    lower_seed_decomposition_max_terms: usize,
    lower_seed_decomposition_status: String,
    lower_seed_decomposition_term_count: Option<usize>,
    lower_seed_decomposition_terms_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    lower_seed_decomposition_error: Option<String>,
    lower_seed_diamond_status: Option<String>,
    lower_seed_diamond_element_count: Option<usize>,
    lower_seed_diamond_gv: Option<String>,
    lower_seed_diamond_error: Option<String>,
    pair_expanded_lower_seed_decomposition_status: Option<String>,
    pair_expanded_lower_seed_decomposition_term_count: Option<usize>,
    pair_expanded_lower_seed_decomposition_terms_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    pair_expanded_lower_seed_diamond_status: Option<String>,
    pair_expanded_lower_seed_diamond_element_count: Option<usize>,
    pair_expanded_lower_seed_diamond_gv: Option<String>,
    pair_expanded_lower_seed_diamond_error: Option<String>,
    path_support_size: Option<usize>,
    path_support_generator_count: Option<usize>,
    path_support_status: Option<String>,
    path_support_gv: Option<String>,
    path_support_error: Option<String>,
    path_support_lookup_status_counts: BTreeMap<String, usize>,
    path_support_source_class_status_counts: BTreeMap<String, usize>,
    path_support_lookup_sample: Vec<CygvPathSupportLookup>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvClosureGenerationCount {
    generation: usize,
    starting_element_count: usize,
    new_element_count: usize,
    total_element_count_after_full_generation: usize,
    truncated_at_limit: bool,
}

#[derive(Clone, Debug, Serialize)]
struct CygvPathPredecessorCandidate {
    predecessor_degree: i128,
    difference_degree: i128,
    series_distance: String,
    predecessor_toric_gv: Option<String>,
    difference_toric_gv: Option<String>,
    predecessor_source_derived_gv: Option<String>,
    difference_source_derived_gv: Option<String>,
    predecessor_known_qn_history_status: String,
    difference_known_qn_history_status: String,
    known_qn_history_pair_status: String,
    predecessor_is_seed: bool,
    difference_is_seed: bool,
    predecessor_is_reduced_seed: bool,
    difference_is_reduced_seed: bool,
    predecessor_first_generation_seed_sum: Option<CygvSeedSumDecomposition>,
    difference_first_generation_seed_sum: Option<CygvSeedSumDecomposition>,
    predecessor_nonzero: Vec<(usize, i64)>,
    difference_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvClosestKnownQnPredecessor {
    predecessor_degree: i128,
    difference_degree: i128,
    series_distance: String,
    predecessor_toric_gv: Option<String>,
    predecessor_source_derived_gv: Option<String>,
    predecessor_known_qn_history_status: String,
    difference_toric_gv: Option<String>,
    difference_source_derived_gv: Option<String>,
    difference_known_qn_history_status: String,
    predecessor_nonzero: Vec<(usize, i64)>,
    difference_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvSeedSumDecomposition {
    reduced_seed_degree: i128,
    seed_degree: i128,
    reduced_seed_toric_gv: Option<String>,
    seed_toric_gv: Option<String>,
    reduced_seed_source_derived_gv: Option<String>,
    seed_source_derived_gv: Option<String>,
    reduced_seed_known_qn_history_status: String,
    seed_known_qn_history_status: String,
    seed_is_reduced_seed: bool,
    seed_pair_reduction_sum: Option<CygvSeedPairDecomposition>,
    reduced_seed_nonzero: Vec<(usize, i64)>,
    seed_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvSeedPairDecomposition {
    lhs_degree: i128,
    rhs_degree: i128,
    lhs_toric_gv: Option<String>,
    rhs_toric_gv: Option<String>,
    lhs_source_derived_gv: Option<String>,
    rhs_source_derived_gv: Option<String>,
    lhs_known_qn_history_status: String,
    rhs_known_qn_history_status: String,
    lhs_is_reduced_seed: bool,
    rhs_is_reduced_seed: bool,
    lhs_nonzero: Vec<(usize, i64)>,
    rhs_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvPathSupportLookup {
    candidate_index: usize,
    side: String,
    degree: i128,
    known_qn_history_status: String,
    toric_gv: Option<String>,
    source_derived_gv: Option<String>,
    path_support_gv: Option<String>,
    path_support_lookup_status: String,
    source_class_status: String,
    source_ray_ambient_nonzero: Option<Vec<(usize, i64)>>,
    matching_missing_target_index: Option<usize>,
    matching_missing_target_degree: Option<i128>,
    matching_missing_target_origin_circuit_pattern: Option<String>,
    matching_missing_target_exact_kind: Option<String>,
    matching_uncovered_source_ray_index: Option<usize>,
    matching_uncovered_source_ray_degree: Option<i128>,
    matching_uncovered_source_ray_origin_circuit_pattern: Option<String>,
    matching_uncovered_source_ray_exact_kind: Option<String>,
    matching_uncovered_source_ray_cms_check_status_counts: BTreeMap<String, usize>,
    matching_uncovered_source_ray_cms_solution_summaries: Vec<CmsGeneralDivisorSolutionSummary>,
    matching_uncovered_source_ray_local_charge_signature: Option<String>,
    matching_uncovered_source_ray_local_cygv_readiness: Option<String>,
    matching_uncovered_source_ray_local_missing_inputs: Vec<String>,
    curve_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct CmsGeneralDivisorSolutionSummary {
    shrinking_divisor_index: usize,
    cms_check_status: String,
    solution_basis_nonzero: Vec<(usize, String)>,
    solution_ambient_basis_nonzero: Vec<(usize, String)>,
    computed_other_normal_degree: String,
    solution_basis_cubic_self_intersection: String,
    local_intersection_tensor_candidate_status: String,
    local_intersection_tensor_candidate: Option<Vec<LocalCygvIntersectionTensorEntry>>,
    local_cygv_primitive_probe: Option<LocalCygvPrimitiveProbe>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct LocalCygvIntersectionTensorEntry {
    indices: [usize; 3],
    value: String,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct LocalCygvPrimitiveProbe {
    status: String,
    q_matrix: Vec<Vec<i64>>,
    grading_vector: Vec<i64>,
    semigroup_elements: Vec<Vec<i64>>,
    candidate_gv: Option<String>,
    unit_tensor_candidate_gv: Option<String>,
    unit_tensor_probe_status: String,
    origin_omitted_unit_tensor_candidate_gv: Option<String>,
    origin_omitted_unit_tensor_probe_status: String,
    expected_toric_gv1_formula_value: Option<String>,
    error: Option<String>,
    unit_tensor_error: Option<String>,
    origin_omitted_unit_tensor_error: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvTargetUnitPhaseProbeSummary {
    target_index: usize,
    degree: i128,
    probe: LocalCygvUnitPhaseProbe,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvUnitPhaseProbe {
    q_matrix: Option<Vec<Vec<i64>>>,
    unit_tensor_candidate_gv: Option<String>,
    unit_tensor_probe_status: String,
    unit_tensor_error: Option<String>,
    unit_tensor_effective_tensor_requirements: Vec<LocalCygvEffectiveTensorRequirement>,
    origin_omitted_q_matrix: Option<Vec<Vec<i64>>>,
    origin_omitted_unit_tensor_candidate_gv: Option<String>,
    origin_omitted_unit_tensor_probe_status: String,
    origin_omitted_unit_tensor_error: Option<String>,
    origin_omitted_unit_tensor_effective_tensor_requirements:
        Vec<LocalCygvEffectiveTensorRequirement>,
    expected_toric_gv1_formula_values: Vec<String>,
    expected_toric_gv1_formula_value_sum: Option<String>,
    unit_tensor_formula_sum_probe_status: String,
    unit_tensor_formula_sum_effective_tensor_requirement:
        Option<LocalCygvEffectiveTensorRequirement>,
    origin_omitted_unit_tensor_formula_sum_probe_status: String,
    origin_omitted_unit_tensor_formula_sum_effective_tensor_requirement:
        Option<LocalCygvEffectiveTensorRequirement>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvEffectiveTensorRequirement {
    expected_gv: String,
    unit_tensor_gv: Option<String>,
    required_tensor_value: Option<String>,
    status: String,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvIntegerTensorScanSummary {
    target_index: usize,
    degree: i128,
    status: String,
    scan_bound: i64,
    expected_toric_gv1_formula_values: Vec<String>,
    entries: Vec<LocalCygvIntegerTensorScanEntry>,
}

#[derive(Clone, Debug, Serialize)]
struct LocalCygvIntegerTensorScanEntry {
    tensor_value: i64,
    candidate_gv: Option<String>,
    status: String,
    error: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvPathSupportSourceRaySummary {
    degree: i128,
    occurrence_count: usize,
    known_qn_history_status_counts: BTreeMap<String, usize>,
    path_support_lookup_status_counts: BTreeMap<String, usize>,
    path_support_gv_counts: BTreeMap<String, usize>,
    uncovered_source_ray_origin_circuit_pattern: Option<String>,
    uncovered_source_ray_exact_kind: Option<String>,
    uncovered_source_ray_cms_check_status_counts: BTreeMap<String, usize>,
    uncovered_source_ray_cms_solution_summaries: Vec<CmsGeneralDivisorSolutionSummary>,
    uncovered_source_ray_local_charge_signature: Option<String>,
    uncovered_source_ray_local_cygv_readiness_counts: BTreeMap<String, usize>,
    uncovered_source_ray_local_missing_input_counts: BTreeMap<String, usize>,
    source_ray_ambient_nonzero: Vec<(usize, i64)>,
    curve_nonzero: Vec<(usize, i64)>,
    occurrences: Vec<CygvPathSupportSourceRayOccurrence>,
}

#[derive(Clone, Debug, Serialize)]
struct CygvPathSupportSourceRayOccurrence {
    target_index: usize,
    candidate_index: usize,
    side: String,
}

#[derive(Clone, Debug, Serialize)]
struct ActiveDecompositionSourceLeafSummary {
    source_status: String,
    occurrence_count: usize,
    degree: Option<i128>,
    curve_nonzero: Vec<(usize, i64)>,
    source_ray_ambient_nonzero: Option<Vec<(usize, i64)>>,
    source_ray_ambient_origin_relation_pattern: Option<String>,
    matching_missing_target_index: Option<usize>,
    matching_missing_target_degree: Option<i128>,
    matching_missing_target_exact_kind: Option<String>,
    matching_uncovered_source_ray_index: Option<usize>,
    matching_uncovered_source_ray_degree: Option<i128>,
    matching_uncovered_source_ray_origin_circuit_pattern: Option<String>,
    matching_uncovered_source_ray_exact_kind: Option<String>,
    matching_uncovered_source_ray_cms_check_status_counts: BTreeMap<String, usize>,
    matching_uncovered_source_ray_cms_solution_summaries: Vec<CmsGeneralDivisorSolutionSummary>,
    matching_uncovered_source_ray_local_charge_signature: Option<String>,
    matching_uncovered_source_ray_local_cygv_readiness: Option<String>,
    matching_uncovered_source_ray_local_missing_inputs: Vec<String>,
    matching_uncovered_source_ray_local_cygv_input_skeleton: Option<LocalCygvInputSkeleton>,
    matching_uncovered_source_ray_local_unit_phase_probe: Option<LocalCygvUnitPhaseProbe>,
    occurrences: Vec<ActiveDecompositionSourceLeafOccurrence>,
}

#[derive(Clone, Debug, Serialize)]
struct ActiveDecompositionSourceLeafOccurrence {
    target_index: usize,
    target_degree: i128,
    target_exact_kind: Option<String>,
    active_generator_index: usize,
}

struct CygvPathSupportSourceRaySummaryBuilder {
    degree: i128,
    known_qn_history_status_counts: BTreeMap<String, usize>,
    path_support_lookup_status_counts: BTreeMap<String, usize>,
    path_support_gv_counts: BTreeMap<String, usize>,
    uncovered_source_ray_origin_circuit_pattern: Option<String>,
    uncovered_source_ray_exact_kind: Option<String>,
    uncovered_source_ray_cms_check_status_counts: BTreeMap<String, usize>,
    uncovered_source_ray_cms_solution_summaries: Vec<CmsGeneralDivisorSolutionSummary>,
    uncovered_source_ray_local_charge_signature: Option<String>,
    uncovered_source_ray_local_cygv_readiness_counts: BTreeMap<String, usize>,
    uncovered_source_ray_local_missing_input_counts: BTreeMap<String, usize>,
    source_ray_ambient_nonzero: Vec<(usize, i64)>,
    curve_nonzero: Vec<(usize, i64)>,
    occurrences: Vec<CygvPathSupportSourceRayOccurrence>,
}

struct PathSupportSourceClassContext {
    status: String,
    source_ray_ambient_nonzero: Option<Vec<(usize, i64)>>,
    matching_missing_target_index: Option<usize>,
    matching_missing_target_degree: Option<i128>,
    matching_missing_target_origin_circuit_pattern: Option<String>,
    matching_missing_target_exact_kind: Option<String>,
    matching_uncovered_source_ray_index: Option<usize>,
    matching_uncovered_source_ray_degree: Option<i128>,
    matching_uncovered_source_ray_origin_circuit_pattern: Option<String>,
    matching_uncovered_source_ray_exact_kind: Option<String>,
    matching_uncovered_source_ray_cms_check_status_counts: BTreeMap<String, usize>,
    matching_uncovered_source_ray_cms_solution_summaries: Vec<CmsGeneralDivisorSolutionSummary>,
    matching_uncovered_source_ray_local_charge_signature: Option<String>,
    matching_uncovered_source_ray_local_cygv_readiness: Option<String>,
    matching_uncovered_source_ray_local_missing_inputs: Vec<String>,
    matching_uncovered_source_ray_local_cygv_input_skeleton: Option<LocalCygvInputSkeleton>,
    matching_uncovered_source_ray_local_unit_phase_probe: Option<LocalCygvUnitPhaseProbe>,
}

struct CygvPathPredecessorStats {
    previous_window_element_count: usize,
    predecessor_difference_count: usize,
    improving_predecessor_difference_count: usize,
    toric_coverage_counts: BTreeMap<String, usize>,
    known_qn_history_counts: BTreeMap<String, usize>,
    closest_distance: f64,
    closest_predecessor: Option<Vec<i64>>,
    closest_difference: Option<Vec<i64>>,
    closest_known_qn_predecessor: Option<CygvClosestKnownQnPredecessor>,
    candidate_sample: Vec<CygvPathPredecessorCandidate>,
}

struct LowerSeedDecompositionProbe {
    status: String,
    term_count: Option<usize>,
    terms_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    terms: Option<Vec<Vec<i64>>>,
    error: Option<String>,
}

struct PairExpandedLowerSeedProbe {
    status: String,
    term_count: Option<usize>,
    terms_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    terms: Option<Vec<Vec<i64>>>,
    error: Option<String>,
}

struct PathSupportGeneratorProbe {
    support_size: Option<usize>,
    generator_count: Option<usize>,
    status: Option<String>,
    gv: Option<String>,
    error: Option<String>,
    lookup_status_counts: BTreeMap<String, usize>,
    source_class_status_counts: BTreeMap<String, usize>,
    lookup_sample: Vec<CygvPathSupportLookup>,
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

fn curve_i64_to_i32(curve: &[i64], label: &str) -> Result<Vec<i32>, String> {
    curve
        .iter()
        .map(|&value| {
            i32::try_from(value).map_err(|_| format!("{label} coordinate does not fit in i32"))
        })
        .collect()
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

fn pair_expanded_lower_seed_diamond_probe(
    target: &[i64],
    pair_expanded: &PairExpandedLowerSeedProbe,
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
    let Some(terms) = pair_expanded.terms.as_deref() else {
        return LowerSeedDiamondProbe {
            status: Some("skipped_no_pair_expanded_lower_seed_decomposition".to_string()),
            element_count: None,
            gv: None,
            error: pair_expanded.error.clone(),
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

fn path_support_generator_probe(
    target: &[i64],
    target_degree: i128,
    candidates: &[CygvPathPredecessorCandidate],
    context: &ValidatedContext<'_>,
    run_path_support_generators: bool,
) -> PathSupportGeneratorProbe {
    if !run_path_support_generators {
        return PathSupportGeneratorProbe {
            support_size: None,
            generator_count: None,
            status: None,
            gv: None,
            error: None,
            lookup_status_counts: BTreeMap::new(),
            source_class_status_counts: BTreeMap::new(),
            lookup_sample: Vec::new(),
        };
    }
    match path_support_generator_probe_inner(target, target_degree, candidates, context) {
        Ok(probe) => probe,
        Err(error) => PathSupportGeneratorProbe {
            support_size: None,
            generator_count: None,
            status: Some("error".to_string()),
            gv: None,
            error: Some(error),
            lookup_status_counts: BTreeMap::new(),
            source_class_status_counts: BTreeMap::new(),
            lookup_sample: Vec::new(),
        },
    }
}

fn path_support_generator_probe_inner(
    target: &[i64],
    target_degree: i128,
    candidates: &[CygvPathPredecessorCandidate],
    context: &ValidatedContext<'_>,
) -> Result<PathSupportGeneratorProbe, String> {
    let support = path_candidate_support(target, candidates);
    if support.is_empty() {
        return Ok(PathSupportGeneratorProbe {
            support_size: Some(0),
            generator_count: Some(0),
            status: Some("skipped_empty_path_support".to_string()),
            gv: None,
            error: None,
            lookup_status_counts: BTreeMap::new(),
            source_class_status_counts: BTreeMap::new(),
            lookup_sample: Vec::new(),
        });
    }
    let mut generators = Vec::new();
    let mut seen = HashSet::new();
    for ray in context.degree_bounded_rays {
        let degree = curve_degree(ray, context.grading)?;
        if degree <= 0 || degree > target_degree {
            continue;
        }
        let supported = ray
            .iter()
            .enumerate()
            .all(|(idx, &value)| value == 0 || support.contains(&idx));
        if supported && seen.insert(ray.clone()) {
            generators.push(ray.clone());
        }
    }
    if !generators.iter().any(|ray| ray.as_slice() == target) {
        generators.push(target.to_vec());
    }
    generators.sort();
    generators.dedup();
    let max_deg = u32::try_from(target_degree)
        .map_err(|_| format!("target degree {target_degree} does not fit in u32"))?;
    let target_i32 = curve_i64_to_i32(target, "path-support target")?;
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
            return Ok(PathSupportGeneratorProbe {
                support_size: Some(support.len()),
                generator_count: Some(generators.len()),
                status: Some("hkty_error".to_string()),
                gv: None,
                error: Some(format!("path_support_generators HKTY failed: {error}")),
                lookup_status_counts: BTreeMap::new(),
                source_class_status_counts: BTreeMap::new(),
                lookup_sample: Vec::new(),
            });
        }
        Err(payload) => {
            return Ok(PathSupportGeneratorProbe {
                support_size: Some(support.len()),
                generator_count: Some(generators.len()),
                status: Some("hkty_panic".to_string()),
                gv: None,
                error: Some(format!(
                    "path_support_generators HKTY panicked: {}",
                    panic_payload_message(payload.as_ref())
                )),
                lookup_status_counts: BTreeMap::new(),
                source_class_status_counts: BTreeMap::new(),
                lookup_sample: Vec::new(),
            });
        }
    };
    let gvs_by_curve = gvs
        .into_iter()
        .map(|(curve, value)| (curve, value.to_string()))
        .collect::<HashMap<_, _>>();
    let gv = gvs_by_curve
        .get(&target_i32)
        .cloned()
        .unwrap_or_else(|| "0".to_string());
    let (lookup_status_counts, source_class_status_counts, lookup_sample) =
        path_support_gv_lookup_sample(candidates, &gvs_by_curve, context)?;
    Ok(PathSupportGeneratorProbe {
        support_size: Some(support.len()),
        generator_count: Some(generators.len()),
        status: Some("computed_path_support_generators".to_string()),
        gv: Some(gv),
        error: None,
        lookup_status_counts,
        source_class_status_counts,
        lookup_sample,
    })
}

fn path_support_gv_lookup_sample(
    candidates: &[CygvPathPredecessorCandidate],
    gvs_by_curve: &HashMap<Vec<i32>, String>,
    context: &ValidatedContext<'_>,
) -> Result<
    (
        BTreeMap<String, usize>,
        BTreeMap<String, usize>,
        Vec<CygvPathSupportLookup>,
    ),
    String,
> {
    let mut lookup_status_counts = BTreeMap::new();
    let mut source_class_status_counts = BTreeMap::new();
    let mut sample = Vec::new();
    for (candidate_index, candidate) in candidates.iter().enumerate() {
        push_path_support_lookup(
            candidate_index,
            "predecessor",
            candidate.predecessor_degree,
            &candidate.predecessor_known_qn_history_status,
            candidate.predecessor_toric_gv.as_ref(),
            candidate.predecessor_source_derived_gv.as_ref(),
            &candidate.predecessor_nonzero,
            gvs_by_curve,
            context,
            &mut lookup_status_counts,
            &mut source_class_status_counts,
            &mut sample,
        )?;
        push_path_support_lookup(
            candidate_index,
            "difference",
            candidate.difference_degree,
            &candidate.difference_known_qn_history_status,
            candidate.difference_toric_gv.as_ref(),
            candidate.difference_source_derived_gv.as_ref(),
            &candidate.difference_nonzero,
            gvs_by_curve,
            context,
            &mut lookup_status_counts,
            &mut source_class_status_counts,
            &mut sample,
        )?;
    }
    Ok((lookup_status_counts, source_class_status_counts, sample))
}

#[allow(clippy::too_many_arguments)]
fn push_path_support_lookup(
    candidate_index: usize,
    side: &str,
    degree: i128,
    known_qn_history_status: &str,
    toric_gv: Option<&String>,
    source_derived_gv: Option<&String>,
    curve_nonzero: &[(usize, i64)],
    gvs_by_curve: &HashMap<Vec<i32>, String>,
    context: &ValidatedContext<'_>,
    lookup_status_counts: &mut BTreeMap<String, usize>,
    source_class_status_counts: &mut BTreeMap<String, usize>,
    sample: &mut Vec<CygvPathSupportLookup>,
) -> Result<(), String> {
    let curve = dense_from_sparse(curve_nonzero, context.dimension)?;
    let curve_i32 = curve_i64_to_i32(&curve, "path-support lookup curve")?;
    let path_support_gv = gvs_by_curve
        .get(&curve_i32)
        .cloned()
        .unwrap_or_else(|| "0".to_string());
    let path_support_lookup_status =
        path_support_lookup_status(known_qn_history_status, &path_support_gv)?;
    *lookup_status_counts
        .entry(path_support_lookup_status.clone())
        .or_insert(0) += 1;
    let source_context = path_support_source_class_context(&curve, context)?;
    *source_class_status_counts
        .entry(source_context.status.clone())
        .or_insert(0) += 1;
    sample.push(CygvPathSupportLookup {
        candidate_index,
        side: side.to_string(),
        degree,
        known_qn_history_status: known_qn_history_status.to_string(),
        toric_gv: toric_gv.cloned(),
        source_derived_gv: source_derived_gv.cloned(),
        path_support_gv: Some(path_support_gv),
        path_support_lookup_status,
        source_class_status: source_context.status,
        source_ray_ambient_nonzero: source_context.source_ray_ambient_nonzero,
        matching_missing_target_index: source_context.matching_missing_target_index,
        matching_missing_target_degree: source_context.matching_missing_target_degree,
        matching_missing_target_origin_circuit_pattern: source_context
            .matching_missing_target_origin_circuit_pattern,
        matching_missing_target_exact_kind: source_context.matching_missing_target_exact_kind,
        matching_uncovered_source_ray_index: source_context.matching_uncovered_source_ray_index,
        matching_uncovered_source_ray_degree: source_context.matching_uncovered_source_ray_degree,
        matching_uncovered_source_ray_origin_circuit_pattern: source_context
            .matching_uncovered_source_ray_origin_circuit_pattern,
        matching_uncovered_source_ray_exact_kind: source_context
            .matching_uncovered_source_ray_exact_kind,
        matching_uncovered_source_ray_cms_check_status_counts: source_context
            .matching_uncovered_source_ray_cms_check_status_counts,
        matching_uncovered_source_ray_cms_solution_summaries: source_context
            .matching_uncovered_source_ray_cms_solution_summaries,
        matching_uncovered_source_ray_local_charge_signature: source_context
            .matching_uncovered_source_ray_local_charge_signature,
        matching_uncovered_source_ray_local_cygv_readiness: source_context
            .matching_uncovered_source_ray_local_cygv_readiness,
        matching_uncovered_source_ray_local_missing_inputs: source_context
            .matching_uncovered_source_ray_local_missing_inputs,
        curve_nonzero: curve_nonzero.to_vec(),
    });
    Ok(())
}

fn path_support_source_class_context(
    curve: &[i64],
    context: &ValidatedContext<'_>,
) -> Result<PathSupportSourceClassContext, String> {
    let matching_missing_target = matching_missing_target_for_curve(curve, context)?;
    let matching_uncovered_source_ray = matching_uncovered_source_ray_for_curve(curve, context)?;
    let matching_uncovered_source_ray_cms_check_status_counts =
        missing_sample_cms_check_status_counts(
            matching_uncovered_source_ray
                .as_ref()
                .and_then(|entry| entry.1.cms_general_divisor_intersection_checks.as_deref()),
        );
    let matching_uncovered_source_ray_cms_solution_summaries =
        cms_general_divisor_solution_summaries(
            matching_uncovered_source_ray.as_ref().map(|entry| entry.1),
            &context.intersection,
        )?;
    let (
        matching_uncovered_source_ray_local_charge_signature,
        matching_uncovered_source_ray_local_cygv_readiness,
        matching_uncovered_source_ray_local_missing_inputs,
        matching_uncovered_source_ray_local_cygv_input_skeleton,
        matching_uncovered_source_ray_local_unit_phase_probe,
    ) = uncovered_source_ray_local_cygv_context(matching_uncovered_source_ray.as_ref())?;
    let Some(ray_context) = context.degree_bounded_ray_context else {
        return Ok(PathSupportSourceClassContext {
            status: if context.source_derived_gv_by_basis.contains_key(curve) {
                "source_ray_known_source_derived_gv".to_string()
            } else if matching_missing_target.is_some() {
                "source_ray_context_missing_but_matches_missing_target".to_string()
            } else {
                "source_ray_context_missing".to_string()
            },
            source_ray_ambient_nonzero: None,
            matching_missing_target_index: matching_missing_target.as_ref().map(|entry| entry.0),
            matching_missing_target_degree: matching_missing_target
                .as_ref()
                .map(|entry| entry.1.degree),
            matching_missing_target_origin_circuit_pattern: matching_missing_target
                .as_ref()
                .and_then(|entry| entry.1.origin_circuit_pattern.clone()),
            matching_missing_target_exact_kind: matching_missing_target
                .as_ref()
                .and_then(|entry| entry.1.real_cone_decomposition_exact_kind.clone()),
            matching_uncovered_source_ray_index: matching_uncovered_source_ray
                .as_ref()
                .map(|entry| entry.0),
            matching_uncovered_source_ray_degree: matching_uncovered_source_ray
                .as_ref()
                .map(|entry| entry.1.degree),
            matching_uncovered_source_ray_origin_circuit_pattern: matching_uncovered_source_ray
                .as_ref()
                .and_then(|entry| entry.1.origin_circuit_pattern.clone()),
            matching_uncovered_source_ray_exact_kind: matching_uncovered_source_ray
                .as_ref()
                .and_then(|entry| entry.1.real_cone_decomposition_exact_kind.clone()),
            matching_uncovered_source_ray_cms_check_status_counts,
            matching_uncovered_source_ray_cms_solution_summaries,
            matching_uncovered_source_ray_local_charge_signature,
            matching_uncovered_source_ray_local_cygv_readiness,
            matching_uncovered_source_ray_local_missing_inputs,
            matching_uncovered_source_ray_local_cygv_input_skeleton,
            matching_uncovered_source_ray_local_unit_phase_probe,
        });
    };

    let source_ray = ray_context
        .iter()
        .find(|sample| sparse_matches_dense(&sample.basis_nonzero, curve));
    let status = if context.covered_toric_gv_by_basis.contains_key(curve) {
        "source_ray_known_toric_covered"
    } else if context.source_derived_gv_by_basis.contains_key(curve) {
        "source_ray_known_source_derived_gv"
    } else if matching_missing_target.is_some() {
        "source_ray_matches_missing_target"
    } else if source_ray.is_some() {
        "source_ray_not_toric_covered"
    } else {
        "not_source_degree_bounded_ray"
    };
    Ok(PathSupportSourceClassContext {
        status: status.to_string(),
        source_ray_ambient_nonzero: source_ray.map(|sample| sample.ambient_nonzero.clone()),
        matching_missing_target_index: matching_missing_target.as_ref().map(|entry| entry.0),
        matching_missing_target_degree: matching_missing_target
            .as_ref()
            .map(|entry| entry.1.degree),
        matching_missing_target_origin_circuit_pattern: matching_missing_target
            .as_ref()
            .and_then(|entry| entry.1.origin_circuit_pattern.clone()),
        matching_missing_target_exact_kind: matching_missing_target
            .as_ref()
            .and_then(|entry| entry.1.real_cone_decomposition_exact_kind.clone()),
        matching_uncovered_source_ray_index: matching_uncovered_source_ray
            .as_ref()
            .map(|entry| entry.0),
        matching_uncovered_source_ray_degree: matching_uncovered_source_ray
            .as_ref()
            .map(|entry| entry.1.degree),
        matching_uncovered_source_ray_origin_circuit_pattern: matching_uncovered_source_ray
            .as_ref()
            .and_then(|entry| entry.1.origin_circuit_pattern.clone()),
        matching_uncovered_source_ray_exact_kind: matching_uncovered_source_ray
            .as_ref()
            .and_then(|entry| entry.1.real_cone_decomposition_exact_kind.clone()),
        matching_uncovered_source_ray_cms_check_status_counts,
        matching_uncovered_source_ray_cms_solution_summaries,
        matching_uncovered_source_ray_local_charge_signature,
        matching_uncovered_source_ray_local_cygv_readiness,
        matching_uncovered_source_ray_local_missing_inputs,
        matching_uncovered_source_ray_local_cygv_input_skeleton,
        matching_uncovered_source_ray_local_unit_phase_probe,
    })
}

fn uncovered_source_ray_local_cygv_context(
    matching_uncovered_source_ray: Option<&(usize, &MissingGvTargetSample)>,
) -> Result<
    (
        Option<String>,
        Option<String>,
        Vec<String>,
        Option<LocalCygvInputSkeleton>,
        Option<LocalCygvUnitPhaseProbe>,
    ),
    String,
> {
    let Some((_, sample)) = matching_uncovered_source_ray else {
        return Ok((None, None, Vec::new(), None, None));
    };
    let signature = sample
        .origin_circuit_affine_support
        .as_ref()
        .map(|support| local_charge_signature_key(&support.local_charge_basis));
    let skeleton =
        local_cygv_input_skeleton(sample, sample.origin_circuit_affine_support.as_ref())?;
    let readiness = skeleton.as_ref().map(local_cygv_actual_call_readiness);
    let missing_inputs = skeleton
        .as_ref()
        .map(|skeleton| skeleton.remaining_uncertified_inputs.clone())
        .unwrap_or_default();
    let unit_phase_probe = skeleton.as_ref().map(|skeleton| {
        local_cygv_unit_phase_probe_from_skeleton(
            skeleton,
            sample_expected_toric_gv1_formula_values(sample),
            sample_expected_toric_gv1_formula_value_sum(sample),
        )
    });
    Ok((
        signature,
        readiness,
        missing_inputs,
        skeleton,
        unit_phase_probe,
    ))
}

fn missing_sample_cms_check_status_counts(
    checks: Option<&[CmsGeneralDivisorIntersectionCheck]>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    let Some(checks) = checks else {
        return counts;
    };
    for check in checks {
        *counts
            .entry(cms_general_divisor_intersection_check_status(check).to_string())
            .or_insert(0) += 1;
    }
    counts
}

fn cms_general_divisor_solution_summaries(
    sample: Option<&MissingGvTargetSample>,
    intersection: &Intersection,
) -> Result<Vec<CmsGeneralDivisorSolutionSummary>, String> {
    let Some(sample) = sample else {
        return Ok(Vec::new());
    };
    let checks = sample.cms_general_divisor_intersection_checks.as_deref();
    let Some(checks) = checks else {
        return Ok(Vec::new());
    };
    let mut summaries = Vec::new();
    for check in checks {
        let cms_check_status = cms_general_divisor_intersection_check_status(check).to_string();
        let Some(solution_basis_nonzero) = check.solution_basis_nonzero.clone() else {
            continue;
        };
        let Some(solution_ambient_basis_nonzero) = check.solution_ambient_basis_nonzero.clone()
        else {
            continue;
        };
        let Some(computed_other_normal_degree) = check.computed_other_normal_degree.clone() else {
            continue;
        };
        let solution_basis_cubic_self_intersection =
            divisor_cubic_self_intersection(intersection, &solution_basis_nonzero)?;
        let expected_toric_gv1_formula_value =
            cms_matching_toric_gv1_formula_value(sample, check.shrinking_divisor_index)?;
        let (
            local_intersection_tensor_candidate,
            local_intersection_tensor_candidate_status,
            local_cygv_primitive_probe,
        ) = local_cygv_intersection_tensor_candidate_from_cms_divisor(
            sample,
            &solution_basis_cubic_self_intersection,
            expected_toric_gv1_formula_value,
            &cms_check_status,
        )?;
        summaries.push(CmsGeneralDivisorSolutionSummary {
            shrinking_divisor_index: check.shrinking_divisor_index,
            cms_check_status,
            solution_basis_cubic_self_intersection,
            local_intersection_tensor_candidate_status,
            local_intersection_tensor_candidate,
            local_cygv_primitive_probe,
            solution_basis_nonzero,
            solution_ambient_basis_nonzero,
            computed_other_normal_degree,
        });
    }
    summaries.sort();
    summaries.dedup();
    Ok(summaries)
}

fn local_cygv_intersection_tensor_candidate_from_cms_divisor(
    sample: &MissingGvTargetSample,
    divisor_cubic_self_intersection: &str,
    expected_toric_gv1_formula_value: Option<i64>,
    cms_check_status: &str,
) -> Result<
    (
        Option<Vec<LocalCygvIntersectionTensorEntry>>,
        String,
        Option<LocalCygvPrimitiveProbe>,
    ),
    String,
> {
    let Some(support) = sample.origin_circuit_affine_support.as_ref() else {
        return Ok((
            None,
            "local_intersection_candidate_missing_affine_support".to_string(),
            None,
        ));
    };
    let Some(skeleton) = local_cygv_input_skeleton(sample, Some(support))? else {
        return Ok((
            None,
            "local_intersection_candidate_missing_skeleton".to_string(),
            None,
        ));
    };
    let has_unit_semigroup = skeleton
        .local_semigroup_generators_candidate
        .as_ref()
        .is_some_and(|generators| generators.as_slice() == [vec![1]]);
    let has_unit_grading = skeleton
        .local_grading_vector_candidate
        .as_ref()
        .is_some_and(|grading| grading.as_slice() == [1]);
    let has_one_parameter_q_layout = skeleton
        .local_cygv_wrapper_q_matrix_candidate
        .as_ref()
        .is_some_and(|q_matrix| q_matrix.len() == 1);
    if !(has_unit_semigroup && has_unit_grading && has_one_parameter_q_layout) {
        return Ok((
            None,
            "local_intersection_candidate_not_one_parameter_source_derived_skeleton".to_string(),
            None,
        ));
    }
    let local_cygv_primitive_probe = Some(local_cygv_primitive_probe_from_cms_divisor(
        &skeleton,
        divisor_cubic_self_intersection,
        expected_toric_gv1_formula_value,
    )?);
    let candidate_status =
        if cms_check_status == "cms_general_divisor_integral_solution_matches_inferred_degree" {
            "candidate_from_cms_divisor_cubic_needs_phase_and_chamber_certificate".to_string()
        } else {
            format!("diagnostic_from_{cms_check_status}_not_promoted")
        };
    Ok((
        Some(vec![LocalCygvIntersectionTensorEntry {
            indices: [0, 0, 0],
            value: divisor_cubic_self_intersection.to_string(),
        }]),
        candidate_status,
        local_cygv_primitive_probe,
    ))
}

fn cms_matching_toric_gv1_formula_value(
    sample: &MissingGvTargetSample,
    shrinking_divisor_index: usize,
) -> Result<Option<i64>, String> {
    let Some(candidates) = sample.cms_general_divisor_shape_candidates.as_deref() else {
        return Ok(None);
    };
    let mut values = candidates
        .iter()
        .filter(|candidate| candidate.shrinking_divisor_index == shrinking_divisor_index)
        .filter_map(|candidate| candidate.toric_gv1_formula_value)
        .collect::<Vec<_>>();
    values.sort_unstable();
    values.dedup();
    match values.as_slice() {
        [] => Ok(None),
        [value] => Ok(Some(*value)),
        _ => Err(format!(
            "CMS divisor {shrinking_divisor_index} has conflicting toric GV1 formula values {values:?}"
        )),
    }
}

fn local_cygv_primitive_probe_from_cms_divisor(
    skeleton: &LocalCygvInputSkeleton,
    divisor_cubic_self_intersection: &str,
    expected_toric_gv1_formula_value: Option<i64>,
) -> Result<LocalCygvPrimitiveProbe, String> {
    let q_matrix = skeleton
        .local_cygv_phase_q_matrix_candidate
        .clone()
        .ok_or_else(|| {
            "local primitive cygv probe is missing phase-selected q-matrix candidate".to_string()
        })?;
    let grading_vector = skeleton
        .local_grading_vector_candidate
        .clone()
        .ok_or_else(|| "local primitive cygv probe is missing grading vector".to_string())?;
    let semigroup_elements = vec![vec![0], vec![1]];
    let expected_toric_gv1_formula_value =
        expected_toric_gv1_formula_value.map(|value| value.to_string());

    let cubic = parse_rational(divisor_cubic_self_intersection)?;
    if cubic.denominator_ref() != &1u32 {
        return Ok(LocalCygvPrimitiveProbe {
            status: "primitive_cygv_probe_blocked_nonintegral_raw_cubic".to_string(),
            q_matrix,
            grading_vector,
            semigroup_elements,
            candidate_gv: None,
            unit_tensor_candidate_gv: None,
            unit_tensor_probe_status: "unit_tensor_probe_not_run_raw_cubic_nonintegral".to_string(),
            origin_omitted_unit_tensor_candidate_gv: None,
            origin_omitted_unit_tensor_probe_status:
                "origin_omitted_unit_tensor_probe_not_run_raw_cubic_nonintegral".to_string(),
            expected_toric_gv1_formula_value,
            error: Some(format!(
                "raw cubic candidate {divisor_cubic_self_intersection} is not integral"
            )),
            unit_tensor_error: None,
            origin_omitted_unit_tensor_error: None,
        });
    }

    let result =
        one_parameter_primitive_cygv_value(&q_matrix, &grading_vector, &semigroup_elements, cubic);
    let gvs = match result {
        Ok(gv) => gv,
        Err(error) => {
            return Ok(LocalCygvPrimitiveProbe {
                status: "primitive_cygv_probe_hkty_error".to_string(),
                q_matrix,
                grading_vector,
                semigroup_elements,
                candidate_gv: None,
                unit_tensor_candidate_gv: None,
                unit_tensor_probe_status: "unit_tensor_probe_not_run_raw_cubic_failed".to_string(),
                origin_omitted_unit_tensor_candidate_gv: None,
                origin_omitted_unit_tensor_probe_status:
                    "origin_omitted_unit_tensor_probe_not_run_raw_cubic_failed".to_string(),
                expected_toric_gv1_formula_value,
                error: Some(error),
                unit_tensor_error: None,
                origin_omitted_unit_tensor_error: None,
            });
        }
    };
    let status = match expected_toric_gv1_formula_value.as_deref() {
        Some(expected) if expected == gvs => {
            "primitive_cygv_probe_matches_expected_formula_but_chamber_uncertified"
        }
        Some(_) => "primitive_cygv_probe_mismatch_raw_cubic_is_not_certified_tensor",
        None => "primitive_cygv_probe_computed_without_expected_formula",
    };
    let unit_result = one_parameter_primitive_cygv_value(
        &q_matrix,
        &grading_vector,
        &semigroup_elements,
        MalachiteRational::from(1),
    );
    let (unit_tensor_candidate_gv, unit_tensor_probe_status, unit_tensor_error) = match unit_result
    {
        Ok(unit_gv) => {
            let unit_status = match expected_toric_gv1_formula_value.as_deref() {
                Some(expected) if expected == unit_gv => {
                    "unit_tensor_probe_matches_expected_formula_but_chamber_uncertified"
                }
                Some(_) => "unit_tensor_probe_mismatch",
                None => "unit_tensor_probe_computed_without_expected_formula",
            };
            (Some(unit_gv), unit_status.to_string(), None)
        }
        Err(error) => (
            None,
            "unit_tensor_probe_hkty_error".to_string(),
            Some(error),
        ),
    };
    let origin_omitted_q_matrix = local_origin_omitted_wrapper_q_matrix(skeleton)?;
    let (
        origin_omitted_unit_tensor_candidate_gv,
        origin_omitted_unit_tensor_probe_status,
        origin_omitted_unit_tensor_error,
    ) = match origin_omitted_q_matrix {
        Some(origin_omitted_q_matrix) => match one_parameter_primitive_cygv_value(
            &origin_omitted_q_matrix,
            &grading_vector,
            &semigroup_elements,
            MalachiteRational::from(1),
        ) {
            Ok(origin_omitted_gv) => {
                let status = match expected_toric_gv1_formula_value.as_deref() {
                    Some(expected) if expected == origin_omitted_gv => {
                        "origin_omitted_unit_tensor_probe_matches_expected_formula_but_phase_uncertified"
                    }
                    Some(_) => "origin_omitted_unit_tensor_probe_mismatch",
                    None => "origin_omitted_unit_tensor_probe_computed_without_expected_formula",
                };
                (Some(origin_omitted_gv), status.to_string(), None)
            }
            Err(error) => (
                None,
                "origin_omitted_unit_tensor_probe_hkty_error".to_string(),
                Some(error),
            ),
        },
        None => (
            None,
            "origin_omitted_unit_tensor_probe_not_run_no_origin_column".to_string(),
            None,
        ),
    };
    Ok(LocalCygvPrimitiveProbe {
        status: status.to_string(),
        q_matrix,
        grading_vector,
        semigroup_elements,
        candidate_gv: Some(gvs),
        unit_tensor_candidate_gv,
        unit_tensor_probe_status,
        origin_omitted_unit_tensor_candidate_gv,
        origin_omitted_unit_tensor_probe_status,
        expected_toric_gv1_formula_value,
        error: None,
        unit_tensor_error,
        origin_omitted_unit_tensor_error,
    })
}

fn local_origin_omitted_wrapper_q_matrix(
    skeleton: &LocalCygvInputSkeleton,
) -> Result<Option<Vec<Vec<i64>>>, String> {
    let Some(q_matrix) = skeleton.local_cygv_wrapper_q_matrix_candidate.as_ref() else {
        return Ok(None);
    };
    local_origin_omitted_wrapper_q_matrix_from_parts(&skeleton.support_point_indices, q_matrix)
}

fn local_origin_omitted_wrapper_q_matrix_from_parts(
    support_point_indices: &[usize],
    q_matrix: &[Vec<i64>],
) -> Result<Option<Vec<Vec<i64>>>, String> {
    let Some(origin_position) = support_point_indices
        .iter()
        .position(|&point_index| point_index == 0)
    else {
        return Ok(None);
    };
    if q_matrix
        .iter()
        .any(|row| row.len() != support_point_indices.len())
    {
        return Err("local origin-omitted q matrix support width mismatch".to_string());
    }
    let mut omitted = Vec::with_capacity(q_matrix.len());
    for row in q_matrix {
        let mut new_row = Vec::with_capacity(row.len().saturating_sub(1));
        for (idx, &entry) in row.iter().enumerate() {
            if idx != origin_position {
                new_row.push(entry);
            }
        }
        if new_row.is_empty() {
            return Ok(None);
        }
        omitted.push(new_row);
    }
    Ok(Some(omitted))
}

fn local_cygv_unit_phase_probe_from_skeleton(
    skeleton: &LocalCygvInputSkeleton,
    expected_toric_gv1_formula_values: Vec<String>,
    expected_toric_gv1_formula_value_sum: Option<String>,
) -> LocalCygvUnitPhaseProbe {
    let q_matrix = skeleton.local_cygv_phase_q_matrix_candidate.clone();
    let grading_vector = skeleton.local_grading_vector_candidate.clone();
    let semigroup_elements = skeleton
        .local_semigroup_generators_candidate
        .as_ref()
        .filter(|generators| generators.as_slice() == [vec![1]])
        .map(|_| vec![vec![0], vec![1]]);

    let (unit_tensor_candidate_gv, unit_tensor_probe_status, unit_tensor_error) =
        match (&q_matrix, &grading_vector, &semigroup_elements) {
            (Some(q_matrix), Some(grading_vector), Some(semigroup_elements))
                if q_matrix.len() == 1 =>
            {
                match one_parameter_primitive_cygv_value(
                    q_matrix,
                    grading_vector,
                    semigroup_elements,
                    MalachiteRational::from(1),
                ) {
                    Ok(value) => {
                        let status =
                            unit_phase_probe_status(&value, &expected_toric_gv1_formula_values);
                        (Some(value), status, None)
                    }
                    Err(error) => (
                        None,
                        "unit_tensor_probe_hkty_error".to_string(),
                        Some(error),
                    ),
                }
            }
            (Some(_), Some(_), Some(_)) => (
                None,
                "unit_tensor_probe_not_run_not_one_parameter_q_matrix".to_string(),
                None,
            ),
            _ => (
                None,
                "unit_tensor_probe_not_run_missing_source_inputs".to_string(),
                None,
            ),
        };

    let origin_omitted_q_matrix = local_origin_omitted_wrapper_q_matrix(skeleton)
        .ok()
        .flatten();
    let (
        origin_omitted_unit_tensor_candidate_gv,
        origin_omitted_unit_tensor_probe_status,
        origin_omitted_unit_tensor_error,
    ) = match (
        &origin_omitted_q_matrix,
        &grading_vector,
        &semigroup_elements,
    ) {
        (Some(q_matrix), Some(grading_vector), Some(semigroup_elements)) if q_matrix.len() == 1 => {
            match one_parameter_primitive_cygv_value(
                q_matrix,
                grading_vector,
                semigroup_elements,
                MalachiteRational::from(1),
            ) {
                Ok(value) => {
                    let status = origin_omitted_unit_phase_probe_status(
                        &value,
                        &expected_toric_gv1_formula_values,
                    );
                    (Some(value), status, None)
                }
                Err(error) => (
                    None,
                    "origin_omitted_unit_tensor_probe_hkty_error".to_string(),
                    Some(error),
                ),
            }
        }
        (Some(_), Some(_), Some(_)) => (
            None,
            "origin_omitted_unit_tensor_probe_not_run_not_one_parameter_q_matrix".to_string(),
            None,
        ),
        (None, _, _) => (
            None,
            "origin_omitted_unit_tensor_probe_not_run_no_origin_column".to_string(),
            None,
        ),
        _ => (
            None,
            "origin_omitted_unit_tensor_probe_not_run_missing_source_inputs".to_string(),
            None,
        ),
    };

    let unit_tensor_effective_tensor_requirements = effective_tensor_requirements_for_unit_probe(
        unit_tensor_candidate_gv.as_deref(),
        &expected_toric_gv1_formula_values,
    );
    let origin_omitted_unit_tensor_effective_tensor_requirements =
        effective_tensor_requirements_for_unit_probe(
            origin_omitted_unit_tensor_candidate_gv.as_deref(),
            &expected_toric_gv1_formula_values,
        );
    let unit_tensor_formula_sum_probe_status = formula_sum_probe_status(
        "unit_tensor",
        unit_tensor_candidate_gv.as_deref(),
        expected_toric_gv1_formula_value_sum.as_deref(),
    );
    let origin_omitted_unit_tensor_formula_sum_probe_status = formula_sum_probe_status(
        "origin_omitted_unit_tensor",
        origin_omitted_unit_tensor_candidate_gv.as_deref(),
        expected_toric_gv1_formula_value_sum.as_deref(),
    );
    let unit_tensor_formula_sum_effective_tensor_requirement = expected_toric_gv1_formula_value_sum
        .as_deref()
        .map(|expected| {
            effective_tensor_requirement_for_unit_probe(
                unit_tensor_candidate_gv.as_deref(),
                expected,
            )
        });
    let origin_omitted_unit_tensor_formula_sum_effective_tensor_requirement =
        expected_toric_gv1_formula_value_sum
            .as_deref()
            .map(|expected| {
                effective_tensor_requirement_for_unit_probe(
                    origin_omitted_unit_tensor_candidate_gv.as_deref(),
                    expected,
                )
            });

    LocalCygvUnitPhaseProbe {
        q_matrix,
        unit_tensor_candidate_gv,
        unit_tensor_probe_status,
        unit_tensor_error,
        unit_tensor_effective_tensor_requirements,
        origin_omitted_q_matrix,
        origin_omitted_unit_tensor_candidate_gv,
        origin_omitted_unit_tensor_probe_status,
        origin_omitted_unit_tensor_error,
        origin_omitted_unit_tensor_effective_tensor_requirements,
        expected_toric_gv1_formula_values,
        expected_toric_gv1_formula_value_sum,
        unit_tensor_formula_sum_probe_status,
        unit_tensor_formula_sum_effective_tensor_requirement,
        origin_omitted_unit_tensor_formula_sum_probe_status,
        origin_omitted_unit_tensor_formula_sum_effective_tensor_requirement,
    }
}

fn effective_tensor_requirements_for_unit_probe(
    unit_tensor_gv: Option<&str>,
    expected_values: &[String],
) -> Vec<LocalCygvEffectiveTensorRequirement> {
    expected_values
        .iter()
        .map(|expected| effective_tensor_requirement_for_unit_probe(unit_tensor_gv, expected))
        .collect()
}

fn effective_tensor_requirement_for_unit_probe(
    unit_tensor_gv: Option<&str>,
    expected: &str,
) -> LocalCygvEffectiveTensorRequirement {
    let blocked = |status: &str| LocalCygvEffectiveTensorRequirement {
        expected_gv: expected.to_string(),
        unit_tensor_gv: unit_tensor_gv.map(ToString::to_string),
        required_tensor_value: None,
        status: status.to_string(),
    };
    let Some(unit_tensor_gv) = unit_tensor_gv else {
        return blocked("effective_tensor_not_computed_missing_unit_probe");
    };
    let Ok(unit) = parse_rational(unit_tensor_gv) else {
        return blocked("effective_tensor_invalid_unit_gv");
    };
    let Ok(expected) = parse_rational(expected) else {
        return blocked("effective_tensor_invalid_expected_gv");
    };
    let zero = MalachiteRational::from(0);
    if unit == zero {
        let status = if expected == zero {
            "effective_tensor_underdetermined_unit_gv_zero_expected_zero"
        } else {
            "effective_tensor_no_scalar_match_unit_gv_zero"
        };
        return blocked(status);
    }
    let expected_gv = expected.to_string();
    let required = expected / unit;
    let status = if required.denominator_ref() == &1u32 {
        "effective_tensor_integral_candidate_but_uncertified"
    } else {
        "effective_tensor_nonintegral_candidate_rejected_by_cygv_threefold_intnums"
    };
    LocalCygvEffectiveTensorRequirement {
        expected_gv,
        unit_tensor_gv: Some(unit_tensor_gv.to_string()),
        required_tensor_value: Some(required.to_string()),
        status: status.to_string(),
    }
}

fn unit_phase_probe_status(value: &str, expected_values: &[String]) -> String {
    if expected_values.is_empty() {
        return "unit_tensor_probe_computed_without_expected_formula".to_string();
    }
    if expected_values.iter().any(|expected| expected == value) {
        "unit_tensor_probe_matches_expected_formula_set_but_uncertified".to_string()
    } else {
        "unit_tensor_probe_mismatch_expected_formula_set".to_string()
    }
}

fn origin_omitted_unit_phase_probe_status(value: &str, expected_values: &[String]) -> String {
    if expected_values.is_empty() {
        return "origin_omitted_unit_tensor_probe_computed_without_expected_formula".to_string();
    }
    if expected_values.iter().any(|expected| expected == value) {
        "origin_omitted_unit_tensor_probe_matches_expected_formula_set_but_phase_uncertified"
            .to_string()
    } else {
        "origin_omitted_unit_tensor_probe_mismatch_expected_formula_set".to_string()
    }
}

fn formula_sum_probe_status(
    label: &str,
    unit_tensor_gv: Option<&str>,
    expected_sum: Option<&str>,
) -> String {
    let Some(expected_sum) = expected_sum else {
        return format!("{label}_formula_sum_not_available");
    };
    let Some(unit_tensor_gv) = unit_tensor_gv else {
        return format!("{label}_formula_sum_probe_not_run_missing_unit_probe");
    };
    if unit_tensor_gv == expected_sum {
        format!("{label}_matches_formula_sum_but_uncertified")
    } else {
        format!("{label}_mismatch_formula_sum")
    }
}

fn one_parameter_primitive_cygv_value(
    q_matrix: &[Vec<i64>],
    grading_vector: &[i64],
    semigroup_elements: &[Vec<i64>],
    tensor_value: MalachiteRational,
) -> Result<String, String> {
    if tensor_value.denominator_ref() != &1u32 {
        return Err(format!(
            "one-parameter primitive cygv tensor value {tensor_value} is not integral"
        ));
    }
    let mut intnums = Intersection::new(1);
    intnums.set(0, 0, 0, Rational::<Finite>::new(tensor_value));
    let gvs = compute_gv_invariants_with_explicit_semigroup(
        semigroup_elements,
        grading_vector,
        q_matrix,
        &intnums,
    )
    .map_err(|error| error.to_string())?;
    Ok(gvs
        .iter()
        .find(|(curve, _)| curve.as_slice() == [1])
        .map_or_else(|| "0".to_string(), |(_, value)| value.to_string()))
}

fn divisor_cubic_self_intersection(
    intersection: &Intersection,
    solution_basis_nonzero: &[(usize, String)],
) -> Result<String, String> {
    let zero = MalachiteRational::from(0);
    let mut coefficients = Vec::new();
    let mut seen = HashSet::new();
    for (idx, value) in solution_basis_nonzero {
        if *idx >= intersection.dim() {
            return Err(format!(
                "CMS solution basis index {idx} is out of bounds for intersection dimension {}",
                intersection.dim()
            ));
        }
        if !seen.insert(*idx) {
            return Err(format!(
                "CMS solution basis index {idx} appears more than once"
            ));
        }
        let coefficient = parse_rational(value)?;
        if coefficient != zero {
            coefficients.push((*idx, coefficient));
        }
    }
    let mut value = zero.clone();
    for (i, coeff_i) in &coefficients {
        for (j, coeff_j) in &coefficients {
            for (k, coeff_k) in &coefficients {
                value += coeff_i * coeff_j * coeff_k * intersection.get(*i, *j, *k).get();
            }
        }
    }
    Ok(value.to_string())
}

fn matching_missing_target_for_curve<'a>(
    curve: &[i64],
    context: &'a ValidatedContext<'_>,
) -> Result<Option<(usize, &'a MissingGvTargetSample)>, String> {
    for (index, sample) in context.stats.sample.iter().enumerate() {
        let sample_curve = dense_from_sparse(&sample.basis_nonzero, context.dimension)?;
        if sample_curve == curve {
            return Ok(Some((index, sample)));
        }
    }
    Ok(None)
}

fn matching_uncovered_source_ray_for_curve<'a>(
    curve: &[i64],
    context: &'a ValidatedContext<'_>,
) -> Result<Option<(usize, &'a MissingGvTargetSample)>, String> {
    let Some(stats) = context.uncovered_source_ray_stats else {
        return Ok(None);
    };
    for (index, sample) in stats.sample.iter().enumerate() {
        let sample_curve = dense_from_sparse(&sample.basis_nonzero, context.dimension)?;
        if sample_curve == curve {
            return Ok(Some((index, sample)));
        }
    }
    Ok(None)
}

fn sparse_matches_dense(sparse: &[(usize, i64)], dense: &[i64]) -> bool {
    let mut seen = HashSet::new();
    for &(idx, value) in sparse {
        if idx >= dense.len() || value == 0 || dense[idx] != value || !seen.insert(idx) {
            return false;
        }
    }
    dense
        .iter()
        .enumerate()
        .all(|(idx, &value)| value == 0 || seen.contains(&idx))
}

fn ambient_origin_relation_pattern(ambient_nonzero: &[(usize, i64)]) -> Option<String> {
    let origin_coefficient = ambient_nonzero
        .iter()
        .find_map(|&(idx, value)| (idx == 0).then_some(value))?;
    let mut negative_counts = BTreeMap::new();
    let mut positive_counts = BTreeMap::new();
    for &(idx, coefficient) in ambient_nonzero {
        if idx == 0 || coefficient == 0 {
            continue;
        }
        if coefficient < 0 {
            *negative_counts.entry(coefficient).or_insert(0usize) += 1;
        } else {
            *positive_counts.entry(coefficient).or_insert(0usize) += 1;
        }
    }
    Some(format!(
        "origin={origin_coefficient};neg={};pos={}",
        coefficient_count_map_label(&negative_counts),
        coefficient_count_map_label(&positive_counts)
    ))
}

fn coefficient_count_map_label(counts: &BTreeMap<i64, usize>) -> String {
    if counts.is_empty() {
        return "{}".to_string();
    }
    let body = counts
        .iter()
        .map(|(coefficient, count)| format!("{coefficient}: {count}"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("{{{body}}}")
}

fn active_decomposition_generator_source_status_counts(
    samples: &[MissingGvTargetSample],
    context: &ValidatedContext<'_>,
    target_index_filter: Option<usize>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for (idx, sample) in samples.iter().enumerate() {
        if !target_index_selected(idx, target_index_filter) {
            continue;
        }
        let Some(active_generators) = sample
            .real_cone_decomposition_active_generator_basis_nonzero
            .as_ref()
        else {
            continue;
        };
        for generator in active_generators {
            let status = match dense_from_sparse(generator, context.dimension)
                .and_then(|curve| active_decomposition_generator_source_status(&curve, context))
            {
                Ok(status) => status,
                Err(error) => format!(
                    "active_generator_source_status_error_{}",
                    status_error_fragment(&error)
                ),
            };
            *counts.entry(status).or_insert(0) += 1;
        }
    }
    counts
}

fn active_decomposition_generator_source_status(
    curve: &[i64],
    context: &ValidatedContext<'_>,
) -> Result<String, String> {
    if context.covered_toric_gv_by_basis.contains_key(curve) {
        return Ok("active_generator_known_toric_covered".to_string());
    }
    if context.source_derived_gv_by_basis.contains_key(curve) {
        return Ok("active_generator_known_source_derived_gv".to_string());
    }
    if matching_missing_target_for_curve(curve, context)?.is_some() {
        return Ok("active_generator_matches_missing_target".to_string());
    }
    if matching_uncovered_source_ray_for_curve(curve, context)?.is_some() {
        return Ok("active_generator_matches_uncovered_source_ray".to_string());
    }
    let Some(ray_context) = context.degree_bounded_ray_context else {
        return Ok("active_generator_source_context_missing".to_string());
    };
    if ray_context
        .iter()
        .any(|sample| sparse_matches_dense(&sample.basis_nonzero, curve))
    {
        return Ok("active_generator_source_ray_not_toric_covered".to_string());
    }
    Ok("active_generator_not_source_degree_bounded_ray".to_string())
}

fn active_decomposition_unresolved_source_leaf_summaries(
    samples: &[MissingGvTargetSample],
    context: &ValidatedContext<'_>,
    target_index_filter: Option<usize>,
) -> Vec<ActiveDecompositionSourceLeafSummary> {
    let mut summaries: BTreeMap<Vec<i64>, ActiveDecompositionSourceLeafSummary> = BTreeMap::new();
    for (target_index, sample) in samples.iter().enumerate() {
        if !target_index_selected(target_index, target_index_filter) {
            continue;
        }
        let Some(active_generators) = sample
            .real_cone_decomposition_active_generator_basis_nonzero
            .as_ref()
        else {
            continue;
        };
        for (active_generator_index, generator) in active_generators.iter().enumerate() {
            let Ok(curve) = dense_from_sparse(generator, context.dimension) else {
                continue;
            };
            let Ok(source_status) = active_decomposition_generator_source_status(&curve, context)
            else {
                continue;
            };
            if source_status == "active_generator_known_toric_covered"
                || source_status == "active_generator_known_source_derived_gv"
            {
                continue;
            }
            let occurrence = ActiveDecompositionSourceLeafOccurrence {
                target_index,
                target_degree: sample.degree,
                target_exact_kind: sample.real_cone_decomposition_exact_kind.clone(),
                active_generator_index,
            };
            match summaries.entry(curve.clone()) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    let summary = entry.get_mut();
                    summary.occurrence_count += 1;
                    summary.occurrences.push(occurrence);
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    let source_context = path_support_source_class_context(&curve, context).ok();
                    let source_ray_ambient_nonzero = source_context
                        .as_ref()
                        .and_then(|context| context.source_ray_ambient_nonzero.clone());
                    let source_ray_ambient_origin_relation_pattern = source_ray_ambient_nonzero
                        .as_deref()
                        .and_then(ambient_origin_relation_pattern);
                    entry.insert(ActiveDecompositionSourceLeafSummary {
                        source_status,
                        occurrence_count: 1,
                        degree: curve_degree(&curve, context.grading).ok(),
                        curve_nonzero: sparse_from_dense(&curve),
                        source_ray_ambient_nonzero,
                        source_ray_ambient_origin_relation_pattern,
                        matching_missing_target_index: source_context
                            .as_ref()
                            .and_then(|context| context.matching_missing_target_index),
                        matching_missing_target_degree: source_context
                            .as_ref()
                            .and_then(|context| context.matching_missing_target_degree),
                        matching_missing_target_exact_kind: source_context
                            .as_ref()
                            .and_then(|context| context.matching_missing_target_exact_kind.clone()),
                        matching_uncovered_source_ray_index: source_context
                            .as_ref()
                            .and_then(|context| context.matching_uncovered_source_ray_index),
                        matching_uncovered_source_ray_degree: source_context
                            .as_ref()
                            .and_then(|context| context.matching_uncovered_source_ray_degree),
                        matching_uncovered_source_ray_origin_circuit_pattern: source_context
                            .as_ref()
                            .and_then(|context| {
                                context
                                    .matching_uncovered_source_ray_origin_circuit_pattern
                                    .clone()
                            }),
                        matching_uncovered_source_ray_exact_kind: source_context.as_ref().and_then(
                            |context| context.matching_uncovered_source_ray_exact_kind.clone(),
                        ),
                        matching_uncovered_source_ray_cms_check_status_counts: source_context
                            .as_ref()
                            .map(|context| {
                                context
                                    .matching_uncovered_source_ray_cms_check_status_counts
                                    .clone()
                            })
                            .unwrap_or_default(),
                        matching_uncovered_source_ray_cms_solution_summaries: source_context
                            .as_ref()
                            .map(|context| {
                                context
                                    .matching_uncovered_source_ray_cms_solution_summaries
                                    .clone()
                            })
                            .unwrap_or_default(),
                        matching_uncovered_source_ray_local_charge_signature: source_context
                            .as_ref()
                            .and_then(|context| {
                                context
                                    .matching_uncovered_source_ray_local_charge_signature
                                    .clone()
                            }),
                        matching_uncovered_source_ray_local_cygv_readiness: source_context
                            .as_ref()
                            .and_then(|context| {
                                context
                                    .matching_uncovered_source_ray_local_cygv_readiness
                                    .clone()
                            }),
                        matching_uncovered_source_ray_local_missing_inputs: source_context
                            .as_ref()
                            .map(|context| {
                                context
                                    .matching_uncovered_source_ray_local_missing_inputs
                                    .clone()
                            })
                            .unwrap_or_default(),
                        matching_uncovered_source_ray_local_cygv_input_skeleton: source_context
                            .as_ref()
                            .and_then(|context| {
                                context
                                    .matching_uncovered_source_ray_local_cygv_input_skeleton
                                    .clone()
                            }),
                        matching_uncovered_source_ray_local_unit_phase_probe: source_context
                            .as_ref()
                            .and_then(|context| {
                                context
                                    .matching_uncovered_source_ray_local_unit_phase_probe
                                    .clone()
                            }),
                        occurrences: vec![occurrence],
                    });
                }
            }
        }
    }
    summaries.into_values().collect()
}

fn active_decomposition_source_leaf_unit_phase_probe_status_counts(
    summaries: &[ActiveDecompositionSourceLeafSummary],
) -> BTreeMap<String, usize> {
    optional_status_counts(
        summaries.iter().map(|summary| {
            summary
                .matching_uncovered_source_ray_local_unit_phase_probe
                .as_ref()
                .map(|probe| probe.unit_tensor_probe_status.as_str())
        }),
        "not_available",
    )
}

fn active_decomposition_source_leaf_origin_omitted_unit_phase_probe_status_counts(
    summaries: &[ActiveDecompositionSourceLeafSummary],
) -> BTreeMap<String, usize> {
    optional_status_counts(
        summaries.iter().map(|summary| {
            summary
                .matching_uncovered_source_ray_local_unit_phase_probe
                .as_ref()
                .map(|probe| probe.origin_omitted_unit_tensor_probe_status.as_str())
        }),
        "not_available",
    )
}

fn active_decomposition_source_leaf_cms_solution_status_counts(
    summaries: &[ActiveDecompositionSourceLeafSummary],
) -> BTreeMap<String, usize> {
    active_decomposition_source_leaf_cms_summary_counts(summaries, |summary| {
        summary.cms_check_status.as_str()
    })
}

fn active_decomposition_source_leaf_cms_candidate_status_counts(
    summaries: &[ActiveDecompositionSourceLeafSummary],
) -> BTreeMap<String, usize> {
    active_decomposition_source_leaf_cms_summary_counts(summaries, |summary| {
        summary.local_intersection_tensor_candidate_status.as_str()
    })
}

fn active_decomposition_source_leaf_cms_summary_counts(
    summaries: &[ActiveDecompositionSourceLeafSummary],
    key: impl Fn(&CmsGeneralDivisorSolutionSummary) -> &str,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        if summary
            .matching_uncovered_source_ray_cms_solution_summaries
            .is_empty()
        {
            *counts
                .entry("no_cms_solution_summary".to_string())
                .or_insert(0) += 1;
            continue;
        }
        for cms_summary in &summary.matching_uncovered_source_ray_cms_solution_summaries {
            *counts.entry(key(cms_summary).to_string()).or_insert(0) += 1;
        }
    }
    counts
}

fn path_support_lookup_status(
    known_qn_history_status: &str,
    path_support_gv: &str,
) -> Result<String, String> {
    let path_support_nonzero = parse_rational(path_support_gv)? != MalachiteRational::from(0);
    let status = match (known_qn_history_status, path_support_nonzero) {
        ("known_nonzero_toric_gv", true) => "path_support_matches_known_nonzero_toric_gv",
        ("known_nonzero_toric_gv", false) => "path_support_misses_known_nonzero_toric_gv",
        ("known_zero_toric_gv", true) => "path_support_contradicts_known_zero_toric_gv",
        ("known_zero_toric_gv", false) => "path_support_matches_known_zero_toric_gv",
        ("known_nonzero_source_gv", true) => "path_support_matches_known_nonzero_source_gv",
        ("known_nonzero_source_gv", false) => "path_support_misses_known_nonzero_source_gv",
        ("known_zero_source_gv", true) => "path_support_contradicts_known_zero_source_gv",
        ("known_zero_source_gv", false) => "path_support_matches_known_zero_source_gv",
        ("unknown_not_toric_covered", true) => "path_support_nonzero_unknown_not_toric_covered",
        ("unknown_not_toric_covered", false) => {
            "path_support_zero_or_absent_unknown_not_toric_covered"
        }
        (_, true) => "path_support_nonzero_unrecognized_history_status",
        (_, false) => "path_support_zero_or_absent_unrecognized_history_status",
    };
    Ok(status.to_string())
}

fn path_support_uncovered_source_ray_summaries(
    targets: &[TargetReport],
) -> Vec<CygvPathSupportSourceRaySummary> {
    let mut summaries: BTreeMap<Vec<(usize, i64)>, CygvPathSupportSourceRaySummaryBuilder> =
        BTreeMap::new();
    for target in targets {
        let Some(probe) = target.cygv_path_history_probe.as_ref() else {
            continue;
        };
        for lookup in &probe.path_support_lookup_sample {
            add_path_support_uncovered_source_ray_lookup(&mut summaries, target.index, lookup);
        }
    }
    let mut out = summaries
        .into_values()
        .map(|builder| CygvPathSupportSourceRaySummary {
            degree: builder.degree,
            occurrence_count: builder.occurrences.len(),
            known_qn_history_status_counts: builder.known_qn_history_status_counts,
            path_support_lookup_status_counts: builder.path_support_lookup_status_counts,
            path_support_gv_counts: builder.path_support_gv_counts,
            uncovered_source_ray_origin_circuit_pattern: builder
                .uncovered_source_ray_origin_circuit_pattern,
            uncovered_source_ray_exact_kind: builder.uncovered_source_ray_exact_kind,
            uncovered_source_ray_cms_check_status_counts: builder
                .uncovered_source_ray_cms_check_status_counts,
            uncovered_source_ray_cms_solution_summaries: builder
                .uncovered_source_ray_cms_solution_summaries,
            uncovered_source_ray_local_charge_signature: builder
                .uncovered_source_ray_local_charge_signature,
            uncovered_source_ray_local_cygv_readiness_counts: builder
                .uncovered_source_ray_local_cygv_readiness_counts,
            uncovered_source_ray_local_missing_input_counts: builder
                .uncovered_source_ray_local_missing_input_counts,
            source_ray_ambient_nonzero: builder.source_ray_ambient_nonzero,
            curve_nonzero: builder.curve_nonzero,
            occurrences: builder.occurrences,
        })
        .collect::<Vec<_>>();
    out.sort_by(|lhs, rhs| {
        lhs.degree
            .cmp(&rhs.degree)
            .then_with(|| lhs.curve_nonzero.cmp(&rhs.curve_nonzero))
    });
    out
}

fn add_path_support_uncovered_source_ray_lookup(
    summaries: &mut BTreeMap<Vec<(usize, i64)>, CygvPathSupportSourceRaySummaryBuilder>,
    target_index: usize,
    lookup: &CygvPathSupportLookup,
) {
    if lookup.source_class_status != "source_ray_not_toric_covered" {
        return;
    }
    let entry = summaries
        .entry(lookup.curve_nonzero.clone())
        .or_insert_with(|| CygvPathSupportSourceRaySummaryBuilder {
            degree: lookup.degree,
            known_qn_history_status_counts: BTreeMap::new(),
            path_support_lookup_status_counts: BTreeMap::new(),
            path_support_gv_counts: BTreeMap::new(),
            uncovered_source_ray_origin_circuit_pattern: lookup
                .matching_uncovered_source_ray_origin_circuit_pattern
                .clone(),
            uncovered_source_ray_exact_kind: lookup
                .matching_uncovered_source_ray_exact_kind
                .clone(),
            uncovered_source_ray_cms_check_status_counts: BTreeMap::new(),
            uncovered_source_ray_cms_solution_summaries: Vec::new(),
            uncovered_source_ray_local_charge_signature: lookup
                .matching_uncovered_source_ray_local_charge_signature
                .clone(),
            uncovered_source_ray_local_cygv_readiness_counts: BTreeMap::new(),
            uncovered_source_ray_local_missing_input_counts: BTreeMap::new(),
            source_ray_ambient_nonzero: lookup
                .source_ray_ambient_nonzero
                .clone()
                .unwrap_or_default(),
            curve_nonzero: lookup.curve_nonzero.clone(),
            occurrences: Vec::new(),
        });
    debug_assert_eq!(entry.degree, lookup.degree);
    if entry.uncovered_source_ray_origin_circuit_pattern.is_none() {
        entry.uncovered_source_ray_origin_circuit_pattern = lookup
            .matching_uncovered_source_ray_origin_circuit_pattern
            .clone();
    }
    if entry.uncovered_source_ray_exact_kind.is_none() {
        entry.uncovered_source_ray_exact_kind =
            lookup.matching_uncovered_source_ray_exact_kind.clone();
    }
    if entry.uncovered_source_ray_local_charge_signature.is_none() {
        entry.uncovered_source_ray_local_charge_signature = lookup
            .matching_uncovered_source_ray_local_charge_signature
            .clone();
    }
    for (status, count) in &lookup.matching_uncovered_source_ray_cms_check_status_counts {
        *entry
            .uncovered_source_ray_cms_check_status_counts
            .entry(status.clone())
            .or_insert(0) += count;
    }
    for solution in &lookup.matching_uncovered_source_ray_cms_solution_summaries {
        if !entry
            .uncovered_source_ray_cms_solution_summaries
            .contains(solution)
        {
            entry
                .uncovered_source_ray_cms_solution_summaries
                .push(solution.clone());
        }
    }
    entry.uncovered_source_ray_cms_solution_summaries.sort();
    entry.uncovered_source_ray_cms_solution_summaries.dedup();
    if let Some(readiness) = &lookup.matching_uncovered_source_ray_local_cygv_readiness {
        *entry
            .uncovered_source_ray_local_cygv_readiness_counts
            .entry(readiness.clone())
            .or_insert(0) += 1;
    }
    for input in &lookup.matching_uncovered_source_ray_local_missing_inputs {
        *entry
            .uncovered_source_ray_local_missing_input_counts
            .entry(input.clone())
            .or_insert(0) += 1;
    }
    *entry
        .known_qn_history_status_counts
        .entry(lookup.known_qn_history_status.clone())
        .or_insert(0) += 1;
    *entry
        .path_support_lookup_status_counts
        .entry(lookup.path_support_lookup_status.clone())
        .or_insert(0) += 1;
    *entry
        .path_support_gv_counts
        .entry(
            lookup
                .path_support_gv
                .clone()
                .unwrap_or_else(|| "missing_path_support_gv".to_string()),
        )
        .or_insert(0) += 1;
    entry.occurrences.push(CygvPathSupportSourceRayOccurrence {
        target_index,
        candidate_index: lookup.candidate_index,
        side: lookup.side.clone(),
    });
}

fn path_support_uncovered_source_ray_degree_counts(
    summaries: &[CygvPathSupportSourceRaySummary],
) -> BTreeMap<i128, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        *counts.entry(summary.degree).or_insert(0) += 1;
    }
    counts
}

fn path_support_uncovered_source_ray_lookup_status_counts(
    summaries: &[CygvPathSupportSourceRaySummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        for (status, count) in &summary.path_support_lookup_status_counts {
            *counts.entry(status.clone()).or_insert(0) += count;
        }
    }
    counts
}

fn path_support_uncovered_source_ray_path_support_gv_counts(
    summaries: &[CygvPathSupportSourceRaySummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        for (gv, count) in &summary.path_support_gv_counts {
            *counts.entry(gv.clone()).or_insert(0) += count;
        }
    }
    counts
}

fn path_support_uncovered_source_ray_local_cygv_primitive_probe_status_counts(
    summaries: &[CygvPathSupportSourceRaySummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        for solution in &summary.uncovered_source_ray_cms_solution_summaries {
            let status = solution
                .local_cygv_primitive_probe
                .as_ref()
                .map_or("local_cygv_primitive_probe_not_run", |probe| {
                    probe.status.as_str()
                });
            *counts.entry(status.to_string()).or_insert(0) += 1;
        }
    }
    counts
}

fn target_cms_solution_status_counts(
    targets: &[TargetReport],
    key: impl Fn(&CmsGeneralDivisorSolutionSummary) -> &str,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for target in targets {
        if target.cms_general_divisor_solution_summaries.is_empty() {
            *counts
                .entry("no_cms_solution_summary".to_string())
                .or_insert(0usize) += 1;
            continue;
        }
        for summary in &target.cms_general_divisor_solution_summaries {
            *counts.entry(key(summary).to_string()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn target_cms_primitive_probe_gv_counts(targets: &[TargetReport]) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for target in targets {
        for summary in &target.cms_general_divisor_solution_summaries {
            let key = summary
                .local_cygv_primitive_probe
                .as_ref()
                .and_then(|probe| probe.candidate_gv.clone())
                .unwrap_or_else(|| "missing_cms_primitive_probe_gv".to_string());
            *counts.entry(key).or_insert(0usize) += 1;
        }
    }
    counts
}

fn target_cms_solution_summary_error_counts(targets: &[TargetReport]) -> BTreeMap<String, usize> {
    optional_status_counts(
        targets
            .iter()
            .map(|target| target.cms_general_divisor_solution_summary_error.as_deref()),
        "no_cms_solution_summary_error",
    )
}

fn path_candidate_support(
    target: &[i64],
    candidates: &[CygvPathPredecessorCandidate],
) -> HashSet<usize> {
    let mut support = sparse_from_dense(target)
        .into_iter()
        .map(|(idx, _)| idx)
        .collect::<HashSet<_>>();
    for candidate in candidates {
        extend_support_from_sparse(&mut support, &candidate.predecessor_nonzero);
        extend_support_from_sparse(&mut support, &candidate.difference_nonzero);
        if let Some(decomposition) = candidate.predecessor_first_generation_seed_sum.as_ref() {
            extend_support_from_seed_sum(&mut support, decomposition);
        }
        if let Some(decomposition) = candidate.difference_first_generation_seed_sum.as_ref() {
            extend_support_from_seed_sum(&mut support, decomposition);
        }
    }
    support
}

fn extend_support_from_seed_sum(
    support: &mut HashSet<usize>,
    decomposition: &CygvSeedSumDecomposition,
) {
    extend_support_from_sparse(support, &decomposition.reduced_seed_nonzero);
    extend_support_from_sparse(support, &decomposition.seed_nonzero);
    if let Some(pair) = decomposition.seed_pair_reduction_sum.as_ref() {
        extend_support_from_sparse(support, &pair.lhs_nonzero);
        extend_support_from_sparse(support, &pair.rhs_nonzero);
    }
}

fn extend_support_from_sparse(support: &mut HashSet<usize>, sparse: &[(usize, i64)]) {
    for &(idx, value) in sparse {
        if value != 0 {
            support.insert(idx);
        }
    }
}

fn pair_expanded_lower_seed_probe(
    lower_seed_decomposition: &LowerSeedDecompositionProbe,
    seed_set: &HashSet<Vec<i64>>,
    reduced_seed_set: &HashSet<Vec<i64>>,
    max_terms: usize,
    max_depth: usize,
) -> PairExpandedLowerSeedProbe {
    let Some(terms) = lower_seed_decomposition.terms.as_deref() else {
        return PairExpandedLowerSeedProbe {
            status: "skipped_no_lower_seed_decomposition".to_string(),
            term_count: None,
            terms_nonzero: None,
            terms: None,
            error: lower_seed_decomposition.error.clone(),
        };
    };
    match pair_expand_terms_to_reduced_seeds(
        terms,
        seed_set,
        reduced_seed_set,
        max_terms,
        max_depth,
    ) {
        Ok(Some(expanded_terms)) => PairExpandedLowerSeedProbe {
            status: "expanded_to_pair_reduced_seed_terms".to_string(),
            term_count: Some(expanded_terms.len()),
            terms_nonzero: Some(
                expanded_terms
                    .iter()
                    .map(|term| sparse_from_dense(term))
                    .collect(),
            ),
            terms: Some(expanded_terms),
            error: None,
        },
        Ok(None) => PairExpandedLowerSeedProbe {
            status: "already_pair_reduced_seed_terms".to_string(),
            term_count: Some(terms.len()),
            terms_nonzero: Some(terms.iter().map(|term| sparse_from_dense(term)).collect()),
            terms: Some(terms.to_vec()),
            error: None,
        },
        Err(error) => PairExpandedLowerSeedProbe {
            status: "error".to_string(),
            term_count: None,
            terms_nonzero: None,
            terms: None,
            error: Some(error),
        },
    }
}

fn pair_expand_terms_to_reduced_seeds(
    terms: &[Vec<i64>],
    seed_set: &HashSet<Vec<i64>>,
    reduced_seed_set: &HashSet<Vec<i64>>,
    max_terms: usize,
    max_depth: usize,
) -> Result<Option<Vec<Vec<i64>>>, String> {
    let mut expanded = Vec::new();
    let mut changed = false;
    for term in terms {
        let mut term_expansion = Vec::new();
        pair_expand_seed_to_reduced_terms(
            term,
            seed_set,
            reduced_seed_set,
            max_depth,
            &mut term_expansion,
        )?;
        if term_expansion.len() != 1 || term_expansion.first() != Some(term) {
            changed = true;
        }
        expanded.extend(term_expansion);
        if expanded.len() > max_terms {
            return Err(format!(
                "pair-expanded lower-seed decomposition exceeds term limit {max_terms}"
            ));
        }
    }
    if !changed {
        return Ok(None);
    }
    Ok(Some(sorted_decomposition(expanded)))
}

fn pair_expand_seed_to_reduced_terms(
    seed: &[i64],
    seed_set: &HashSet<Vec<i64>>,
    reduced_seed_set: &HashSet<Vec<i64>>,
    remaining_depth: usize,
    out: &mut Vec<Vec<i64>>,
) -> Result<(), String> {
    if reduced_seed_set.contains(seed) {
        out.push(seed.to_vec());
        return Ok(());
    }
    if remaining_depth == 0 {
        out.push(seed.to_vec());
        return Ok(());
    }
    let Some((lhs, rhs)) = seed_pair_reduction_terms(seed, seed_set)? else {
        out.push(seed.to_vec());
        return Ok(());
    };
    pair_expand_seed_to_reduced_terms(&lhs, seed_set, reduced_seed_set, remaining_depth - 1, out)?;
    pair_expand_seed_to_reduced_terms(&rhs, seed_set, reduced_seed_set, remaining_depth - 1, out)?;
    Ok(())
}

fn seed_pair_reduction_terms(
    seed: &[i64],
    seed_set: &HashSet<Vec<i64>>,
) -> Result<Option<(Vec<i64>, Vec<i64>)>, String> {
    let mut sorted_seeds = seed_set.iter().collect::<Vec<_>>();
    sorted_seeds.sort();
    for lhs in sorted_seeds {
        let rhs = checked_vector_difference(seed, lhs)?;
        if seed_set.contains(&rhs) {
            return Ok(Some((lhs.clone(), rhs)));
        }
    }
    Ok(None)
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

#[derive(Clone, Debug, Serialize)]
struct OriginCircuitWitnessDomainSummary {
    target_index: usize,
    degree: i128,
    witness_index: usize,
    first_facet_exclusive_point: usize,
    second_facet_exclusive_point: usize,
    facet_context_status: String,
    relation_support_size: usize,
    shared_facet_support_size: usize,
    facet_union_support_size: usize,
    relation_generator_count: Option<usize>,
    shared_facet_generator_count: Option<usize>,
    facet_union_generator_count: Option<usize>,
    relation_rank: Option<usize>,
    shared_facet_rank: Option<usize>,
    facet_union_rank: Option<usize>,
    relation_support_face_profile: String,
    shared_facet_face_profile: String,
    facet_union_face_profile: String,
    relation_source_status_counts: BTreeMap<String, usize>,
    shared_facet_source_status_counts: BTreeMap<String, usize>,
    facet_union_source_status_counts: BTreeMap<String, usize>,
    relation_support_face_certificate_status: String,
    shared_facet_face_certificate_status: String,
    facet_union_face_certificate_status: String,
}

#[derive(Clone, Debug, Serialize)]
struct OriginCircuitWitnessDomainUnresolvedGeneratorSummary {
    degree: i128,
    source_status: String,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
    occurrence_count: usize,
    occurrences: Vec<OriginCircuitWitnessDomainGeneratorOccurrence>,
}

#[derive(Clone, Debug, Serialize)]
struct OriginCircuitWitnessDomainGeneratorOccurrence {
    target_index: usize,
    witness_index: usize,
    domain_kind: String,
    target_degree: i128,
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
    let witnesses = origin_circuit_witnesses(sample);
    if witnesses.is_empty() {
        return None;
    }

    let mut relation_support = HashSet::new();
    let mut shared_facet = HashSet::new();
    let mut facet_union = HashSet::new();
    for witness in witnesses {
        relation_support.extend(
            witness
                .relation_points
                .iter()
                .map(|point| point.point_index),
        );
        let first_facet = witness.first_facet.iter().copied().collect::<HashSet<_>>();
        let second_facet = witness.second_facet.iter().copied().collect::<HashSet<_>>();
        shared_facet.extend(first_facet.intersection(&second_facet).copied());
        shared_facet.insert(0);
        shared_facet.insert(witness.first_facet_exclusive_point);
        shared_facet.insert(witness.second_facet_exclusive_point);
        facet_union.extend(first_facet.union(&second_facet).copied());
        facet_union.insert(0);
    }

    Some(OriginCircuitAmbientSupportSets {
        relation_support,
        shared_facet,
        facet_union,
    })
}

fn origin_circuit_witness_ambient_support_sets(
    witness: &OriginCircuitWitnessSample,
) -> OriginCircuitAmbientSupportSets {
    let relation_support = origin_circuit_witness_relation_signature(witness)
        .into_iter()
        .map(|(idx, _)| idx)
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

    OriginCircuitAmbientSupportSets {
        relation_support,
        shared_facet,
        facet_union,
    }
}

fn origin_circuit_witness_domain_summaries(
    samples: &[MissingGvTargetSample],
    context: &ValidatedContext<'_>,
    target_index_filter: Option<usize>,
    certify_domains: bool,
    generator_limit: usize,
) -> Vec<OriginCircuitWitnessDomainSummary> {
    let mut summaries = Vec::new();
    for (target_index, sample) in samples.iter().enumerate() {
        if target_index_filter.is_some_and(|filter| filter != target_index) {
            continue;
        }
        for (witness_index, witness) in origin_circuit_witnesses(sample).into_iter().enumerate() {
            summaries.push(origin_circuit_witness_domain_summary(
                target_index,
                sample.degree,
                witness_index,
                witness,
                context,
                certify_domains,
                generator_limit,
            ));
        }
    }
    summaries
}

fn origin_circuit_witness_domain_summary(
    target_index: usize,
    degree: i128,
    witness_index: usize,
    witness: &OriginCircuitWitnessSample,
    context: &ValidatedContext<'_>,
    certify_domains: bool,
    generator_limit: usize,
) -> OriginCircuitWitnessDomainSummary {
    let supports = origin_circuit_witness_ambient_support_sets(witness);
    let relation = origin_circuit_witness_domain_stats(
        context,
        degree,
        &supports.relation_support,
        certify_domains,
        generator_limit,
    );
    let shared = origin_circuit_witness_domain_stats(
        context,
        degree,
        &supports.shared_facet,
        certify_domains,
        generator_limit,
    );
    let union = origin_circuit_witness_domain_stats(
        context,
        degree,
        &supports.facet_union,
        certify_domains,
        generator_limit,
    );
    OriginCircuitWitnessDomainSummary {
        target_index,
        degree,
        witness_index,
        first_facet_exclusive_point: witness.first_facet_exclusive_point,
        second_facet_exclusive_point: witness.second_facet_exclusive_point,
        facet_context_status: origin_circuit_witness_facet_context_status(witness).to_string(),
        relation_support_size: supports.relation_support.len(),
        shared_facet_support_size: supports.shared_facet.len(),
        facet_union_support_size: supports.facet_union.len(),
        relation_generator_count: relation.generator_count,
        shared_facet_generator_count: shared.generator_count,
        facet_union_generator_count: union.generator_count,
        relation_rank: relation.rank,
        shared_facet_rank: shared.rank,
        facet_union_rank: union.rank,
        relation_support_face_profile: relation.support_face_profile,
        shared_facet_face_profile: shared.support_face_profile,
        facet_union_face_profile: union.support_face_profile,
        relation_source_status_counts: relation.source_status_counts,
        shared_facet_source_status_counts: shared.source_status_counts,
        facet_union_source_status_counts: union.source_status_counts,
        relation_support_face_certificate_status: relation.certificate_status,
        shared_facet_face_certificate_status: shared.certificate_status,
        facet_union_face_certificate_status: union.certificate_status,
    }
}

struct OriginCircuitWitnessDomainStats {
    generator_count: Option<usize>,
    rank: Option<usize>,
    support_face_profile: String,
    source_status_counts: BTreeMap<String, usize>,
    certificate_status: String,
}

fn origin_circuit_witness_domain_stats(
    context: &ValidatedContext<'_>,
    degree: i128,
    allowed_ambient_support: &HashSet<usize>,
    certify_domain: bool,
    generator_limit: usize,
) -> OriginCircuitWitnessDomainStats {
    let Some(ray_context) = context.degree_bounded_ray_context else {
        return OriginCircuitWitnessDomainStats {
            generator_count: None,
            rank: None,
            support_face_profile: "missing_degree_bounded_mori_ray_context".to_string(),
            source_status_counts: BTreeMap::new(),
            certificate_status: "missing_degree_bounded_mori_ray_context".to_string(),
        };
    };
    let generators = degree_bounded_ray_context_support_generators(
        ray_context,
        degree,
        allowed_ambient_support,
        context.dimension,
    );
    let generators = match generators {
        Ok(generators) => generators,
        Err(error) => {
            return OriginCircuitWitnessDomainStats {
                generator_count: None,
                rank: None,
                support_face_profile: format!(
                    "origin_witness_domain_error_{}",
                    status_error_fragment(&error)
                ),
                source_status_counts: BTreeMap::new(),
                certificate_status: format!(
                    "origin_witness_domain_error_{}",
                    status_error_fragment(&error)
                ),
            };
        }
    };
    let rank = if generators.is_empty() {
        None
    } else {
        match curve_row_span_rank(&generators) {
            Ok(rank) => Some(rank),
            Err(error) => {
                return OriginCircuitWitnessDomainStats {
                    generator_count: Some(generators.len()),
                    rank: None,
                    support_face_profile: format!(
                        "origin_witness_domain_error_{}",
                        status_error_fragment(&error.to_string())
                    ),
                    source_status_counts: BTreeMap::new(),
                    certificate_status: format!(
                        "origin_witness_domain_error_{}",
                        status_error_fragment(&error.to_string())
                    ),
                };
            }
        }
    };
    let support_face_profile =
        support_face_pre_lp_profile(generators.len(), rank, context.dimension, generator_limit);
    let source_status_counts = source_status_counts_for_generators(&generators, context);
    let certificate_status = if certify_domain {
        origin_circuit_ambient_support_face_certificate_status(
            ray_context,
            degree,
            allowed_ambient_support,
            context,
            generator_limit,
        )
    } else {
        "not_run".to_string()
    };
    OriginCircuitWitnessDomainStats {
        generator_count: Some(generators.len()),
        rank,
        support_face_profile,
        source_status_counts,
        certificate_status,
    }
}

fn support_face_pre_lp_profile(
    generator_count: usize,
    rank: Option<usize>,
    dimension: usize,
    generator_limit: usize,
) -> String {
    if generator_count == 0 {
        return "support_face_empty_domain".to_string();
    }
    let Some(rank) = rank else {
        return "support_face_rank_error".to_string();
    };
    if rank > dimension {
        return format!("support_face_rank_{rank}_exceeds_dim_{dimension}");
    }
    let generator_bucket = generator_count_bucket(generator_count);
    let codimension = dimension - rank;
    if codimension == 0 {
        return format!("support_face_full_dimensional_generators_{generator_bucket}");
    }
    if codimension == 1 {
        return format!(
            "support_face_codim_1_exact_kernel_candidate_generators_{generator_bucket}"
        );
    }
    if generator_count > generator_limit {
        return format!(
            "support_face_codim_{codimension}_lp_skipped_generator_limit_{generator_limit}_generators_{generator_bucket}"
        );
    }
    format!("support_face_codim_{codimension}_lp_candidate_generators_{generator_bucket}")
}

fn generator_count_bucket(generator_count: usize) -> &'static str {
    match generator_count {
        0 => "0",
        1 => "1",
        2..=4 => "2_4",
        5..=16 => "5_16",
        17..=64 => "17_64",
        65..=256 => "65_256",
        257..=1024 => "257_1024",
        _ => "gt_1024",
    }
}

fn source_status_counts_for_generators(
    generators: &[Vec<i64>],
    context: &ValidatedContext<'_>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for generator in generators {
        let status = active_decomposition_generator_source_status(generator, context)
            .unwrap_or_else(|error| {
                format!(
                    "generator_source_status_error_{}",
                    status_error_fragment(&error)
                )
            });
        *counts.entry(status).or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_witness_domain_status_counts(
    summaries: &[OriginCircuitWitnessDomainSummary],
    status: impl Fn(&OriginCircuitWitnessDomainSummary) -> &str,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        *counts.entry(status(summary).to_string()).or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_witness_domain_unresolved_generator_summaries(
    samples: &[MissingGvTargetSample],
    context: &ValidatedContext<'_>,
    target_index_filter: Option<usize>,
    domain_kind_filter: Option<&str>,
) -> Vec<OriginCircuitWitnessDomainUnresolvedGeneratorSummary> {
    let Some(ray_context) = context.degree_bounded_ray_context else {
        return Vec::new();
    };
    let mut summaries: BTreeMap<Vec<i64>, OriginCircuitWitnessDomainUnresolvedGeneratorSummary> =
        BTreeMap::new();
    for (target_index, sample) in samples.iter().enumerate() {
        if target_index_filter.is_some_and(|filter| filter != target_index) {
            continue;
        }
        for (witness_index, witness) in origin_circuit_witnesses(sample).into_iter().enumerate() {
            let supports = origin_circuit_witness_ambient_support_sets(witness);
            let domains = [
                ("relation_support", &supports.relation_support),
                ("shared_facet", &supports.shared_facet),
                ("facet_union", &supports.facet_union),
            ];
            for (domain_kind, support) in domains {
                if domain_kind_filter.is_some_and(|filter| filter != domain_kind) {
                    continue;
                }
                let mut seen_in_domain = HashSet::new();
                for ray in ray_context {
                    if ray.degree <= 0 || ray.degree > sample.degree {
                        continue;
                    }
                    if !ray
                        .ambient_nonzero
                        .iter()
                        .all(|(idx, _)| support.contains(idx))
                    {
                        continue;
                    }
                    let Ok(basis_ray) = dense_from_sparse(&ray.basis_nonzero, context.dimension)
                    else {
                        continue;
                    };
                    if !seen_in_domain.insert(basis_ray.clone()) {
                        continue;
                    }
                    let source_status =
                        match active_decomposition_generator_source_status(&basis_ray, context) {
                            Ok(status) => status,
                            Err(error) => format!(
                                "generator_source_status_error_{}",
                                status_error_fragment(&error)
                            ),
                        };
                    if source_status == "active_generator_known_toric_covered"
                        || source_status == "active_generator_known_source_derived_gv"
                    {
                        continue;
                    }
                    let occurrence = OriginCircuitWitnessDomainGeneratorOccurrence {
                        target_index,
                        witness_index,
                        domain_kind: domain_kind.to_string(),
                        target_degree: sample.degree,
                    };
                    match summaries.entry(basis_ray) {
                        std::collections::btree_map::Entry::Occupied(mut entry) => {
                            let summary = entry.get_mut();
                            summary.occurrence_count += 1;
                            if summary.occurrences.len()
                                < ORIGIN_CIRCUIT_WITNESS_DOMAIN_OCCURRENCE_SAMPLE_LIMIT
                            {
                                summary.occurrences.push(occurrence);
                            }
                        }
                        std::collections::btree_map::Entry::Vacant(entry) => {
                            entry.insert(OriginCircuitWitnessDomainUnresolvedGeneratorSummary {
                                degree: ray.degree,
                                source_status,
                                ambient_nonzero: ray.ambient_nonzero.clone(),
                                basis_nonzero: ray.basis_nonzero.clone(),
                                occurrence_count: 1,
                                occurrences: vec![occurrence],
                            });
                        }
                    }
                }
            }
        }
    }
    let mut summaries = summaries.into_values().collect::<Vec<_>>();
    summaries.sort_by(|left, right| {
        (
            left.degree,
            left.source_status.as_str(),
            left.basis_nonzero.as_slice(),
        )
            .cmp(&(
                right.degree,
                right.source_status.as_str(),
                right.basis_nonzero.as_slice(),
            ))
    });
    summaries
}

fn origin_circuit_witness_domain_unresolved_generator_status_counts_for(
    summaries: &[OriginCircuitWitnessDomainUnresolvedGeneratorSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        *counts
            .entry(summary.source_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_witness_domain_unresolved_generator_degree_counts_for(
    summaries: &[OriginCircuitWitnessDomainUnresolvedGeneratorSummary],
) -> BTreeMap<i128, usize> {
    let mut counts = BTreeMap::new();
    for summary in summaries {
        *counts.entry(summary.degree).or_insert(0usize) += 1;
    }
    counts
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
    let codimension = context.dimension.saturating_sub(rank);
    if generators.len() == 1 && codimension > 1 {
        return format!(
            "origin_support_skipped_single_generator_higher_codimension_rank_{rank}_dim_{}",
            context.dimension
        );
    }
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
    if let Some(probe) = target_exact_non_extremal_decomposition_probe(sample, target, context) {
        return Some(probe);
    }
    if !run {
        return None;
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
    let (local_semigroup_generators_candidate, local_semigroup_generator_status) =
        local_cygv_semigroup_generators_candidate(
            &orientation_candidates,
            local_q_matrix_orientation_candidate,
        );
    let (
        local_cygv_q_matrix_rows_candidate,
        local_cygv_wrapper_q_matrix_candidate,
        local_cygv_q_matrix_layout_status,
    ) = local_cygv_q_matrix_layout_candidate(
        &orientation_candidates,
        local_q_matrix_orientation_candidate,
    );
    let (local_cygv_phase_q_matrix_candidate, local_cygv_q_matrix_phase_status) =
        local_cygv_q_matrix_phase_candidate(
            &support_point_indices,
            local_cygv_wrapper_q_matrix_candidate.as_deref(),
        );
    let (local_grading_vector_candidate, local_grading_vector_status) =
        local_cygv_grading_vector_candidate(&orientation_candidates);
    let mut remaining_uncertified_inputs = Vec::new();
    if local_semigroup_generators_candidate.is_none() {
        remaining_uncertified_inputs.push("local_semigroup_generators".to_string());
    }
    if local_grading_vector_candidate.is_none() {
        remaining_uncertified_inputs.push("local_grading_vector".to_string());
    }
    if local_cygv_phase_q_matrix_candidate.is_none() {
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
        local_semigroup_generators_candidate,
        local_semigroup_generator_status,
        orientation_candidates,
        local_q_matrix_orientation_candidate,
        local_q_matrix_orientation_status,
        local_cygv_q_matrix_rows_candidate,
        local_cygv_wrapper_q_matrix_candidate,
        local_cygv_q_matrix_layout_status,
        local_cygv_phase_q_matrix_candidate,
        local_cygv_q_matrix_phase_status,
        local_grading_vector_candidate,
        local_grading_vector_status,
        remaining_uncertified_inputs,
    }))
}

fn local_cygv_semigroup_generators_candidate(
    orientation_candidates: &[LocalCygvOrientationCandidate],
    orientation_candidate: Option<i64>,
) -> (Option<Vec<Vec<i64>>>, String) {
    let Some(sign) = orientation_candidate else {
        return (
            None,
            "local_semigroup_generators_blocked_no_orientation".to_string(),
        );
    };
    let Some(candidate) = orientation_candidates
        .iter()
        .find(|candidate| candidate.overall_charge_basis_sign == sign)
    else {
        return (
            None,
            "local_semigroup_generators_missing_selected_orientation".to_string(),
        );
    };
    let Some(direction) = candidate.target_primitive_direction.as_ref() else {
        return (
            None,
            "local_semigroup_generators_blocked_no_target_direction".to_string(),
        );
    };
    if direction.len() != 1 || direction[0] != 1 {
        return (
            None,
            "local_semigroup_generators_blocked_not_one_parameter_unit".to_string(),
        );
    }
    (
        Some(vec![vec![1]]),
        "source_derived_one_parameter_unit_semigroup".to_string(),
    )
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

fn local_cygv_q_matrix_phase_candidate(
    support_point_indices: &[usize],
    wrapper_q_matrix: Option<&[Vec<i64>]>,
) -> (Option<Vec<Vec<i64>>>, String) {
    let Some(wrapper_q_matrix) = wrapper_q_matrix else {
        return (
            None,
            "local_q_matrix_phase_blocked_no_wrapper_q_matrix".to_string(),
        );
    };
    let mut candidates = Vec::new();
    if local_cygv_q_matrix_cy_dimension(wrapper_q_matrix) == Some(3) {
        candidates.push((
            "including_origin",
            wrapper_q_matrix
                .iter()
                .map(std::clone::Clone::clone)
                .collect::<Vec<_>>(),
        ));
    }
    if let Ok(Some(origin_omitted)) =
        local_origin_omitted_wrapper_q_matrix_from_parts(support_point_indices, wrapper_q_matrix)
    {
        if local_cygv_q_matrix_cy_dimension(&origin_omitted) == Some(3) {
            candidates.push(("omitting_origin", origin_omitted));
        }
    }

    match candidates.as_slice() {
        [] => (
            None,
            "local_q_matrix_phase_blocked_no_compact_threefold_phase".to_string(),
        ),
        [(label, matrix)] => (
            Some(matrix.clone()),
            format!("source_derived_unique_compact_threefold_phase_{label}"),
        ),
        _ => (
            None,
            "local_q_matrix_phase_ambiguous_multiple_compact_threefold_phases".to_string(),
        ),
    }
}

fn local_cygv_q_matrix_cy_dimension(q_matrix: &[Vec<i64>]) -> Option<i64> {
    let first_row = q_matrix.first()?;
    if first_row.is_empty() || q_matrix.iter().any(|row| row.len() != first_row.len()) {
        return None;
    }
    let h11 = i64::try_from(q_matrix.len()).ok()?;
    let h11pd = i64::try_from(first_row.len()).ok()?;
    Some(h11pd - h11 - 1)
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
    let Some(witness) = origin_circuit_witnesses(sample)
        .into_iter()
        .find(|witness| !witness.relation_points.is_empty())
    else {
        return Ok(Some(support));
    };

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
    let mut covered_toric_gv_by_basis = HashMap::new();
    if let Some(covered_context) = context.covered_toric_gv_context_for_missing.as_ref() {
        for (idx, sample) in covered_context.iter().enumerate() {
            let basis_class = dense_from_sparse(&sample.basis_nonzero, dimension).map_err(|e| {
                format!("covered toric GV context sample {idx} has invalid basis support: {e}")
            })?;
            let computed_degree = curve_degree(&basis_class, grading).map_err(|e| {
                format!("covered toric GV context sample {idx} has invalid basis degree: {e}")
            })?;
            if computed_degree != sample.degree {
                return Err(format!(
                    "covered toric GV context sample {idx} declares degree {} but computes {computed_degree}",
                    sample.degree
                ));
            }
            if sample.degree > degree_bound {
                return Err(format!(
                    "covered toric GV context sample {idx} has degree {} above bound {degree_bound}",
                    sample.degree
                ));
            }
            if covered_toric_gv_by_basis
                .insert(basis_class, sample.gv.clone())
                .is_some()
            {
                return Err(format!(
                    "covered toric GV context sample {idx} duplicates a basis class"
                ));
            }
        }
    }
    let kappa_entries = context
        .corrected_kappa_basis_for_missing
        .as_ref()
        .ok_or_else(|| "context is missing corrected kappa entries".to_string())?;
    let intersection = reconstruct_intersection(dimension, kappa_entries)?;
    let mut source_derived_gv_by_basis = source_derived_gv_by_basis(
        context.uncovered_source_ray_stats_for_missing.as_ref(),
        dimension,
    )?;
    add_toric_diagnostic_context_gvs(
        &mut source_derived_gv_by_basis,
        context
            .uncovered_source_ray_toric_diagnostic_sample
            .as_deref(),
        dimension,
        grading,
        degree_bound,
    )?;
    add_toric_diagnostic_context_gvs(
        &mut source_derived_gv_by_basis,
        context
            .degree_bounded_toric_gv_diagnostic_context_for_missing
            .as_deref(),
        dimension,
        grading,
        degree_bound,
    )?;
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
        covered_toric_gv_by_basis,
        source_derived_gv_by_basis,
        intersection,
        stats,
        uncovered_source_ray_stats: context.uncovered_source_ray_stats_for_missing.as_ref(),
        shared_facet_unresolved_source_ray_stats: context
            .shared_facet_unresolved_source_ray_stats_for_missing
            .as_ref(),
    })
}

fn source_derived_gv_by_basis(
    stats: Option<&MissingGvTargetStats>,
    dimension: usize,
) -> Result<HashMap<Vec<i64>, String>, String> {
    let mut out = HashMap::new();
    let Some(stats) = stats else {
        return Ok(out);
    };
    for (index, sample) in stats.sample.iter().enumerate() {
        let basis_class = dense_from_sparse(&sample.basis_nonzero, dimension).map_err(|e| {
            format!("uncovered source-ray stats sample {index} has invalid basis support: {e}")
        })?;
        let Some(gv) = source_derived_gv_for_sample(sample)? else {
            continue;
        };
        if let Some(previous) = out.insert(basis_class, gv.clone())
            && previous != gv
        {
            return Err(format!(
                "uncovered source-ray stats sample {index} gives conflicting source-derived GV values {previous} and {gv}"
            ));
        }
    }
    Ok(out)
}

fn add_toric_diagnostic_context_gvs(
    out: &mut HashMap<Vec<i64>, String>,
    samples: Option<&[ToricGvDiagnosticContextSample]>,
    dimension: usize,
    grading: &[i64],
    degree_bound: i128,
) -> Result<(), String> {
    let Some(samples) = samples else {
        return Ok(());
    };
    for (index, sample) in samples.iter().enumerate() {
        let basis_class = dense_from_sparse(&sample.basis_nonzero, dimension).map_err(|e| {
            format!("uncovered source-ray toric diagnostic sample {index} has invalid basis support: {e}")
        })?;
        let computed_degree = curve_degree(&basis_class, grading).map_err(|e| {
            format!(
                "uncovered source-ray toric diagnostic sample {index} has invalid basis degree: {e}"
            )
        })?;
        if computed_degree != sample.degree {
            return Err(format!(
                "uncovered source-ray toric diagnostic sample {index} declares degree {} but computes {computed_degree}",
                sample.degree
            ));
        }
        if sample.degree > degree_bound {
            return Err(format!(
                "uncovered source-ray toric diagnostic sample {index} has degree {} above context bound {degree_bound}",
                sample.degree
            ));
        }
        if let Some(previous) = out.insert(basis_class, sample.gv.clone())
            && previous != sample.gv
        {
            return Err(format!(
                "uncovered source-ray toric diagnostic sample {index} gives conflicting source-derived GV values {previous} and {}",
                sample.gv
            ));
        }
    }
    Ok(())
}

fn source_derived_gv_for_sample(sample: &MissingGvTargetSample) -> Result<Option<String>, String> {
    let Some(candidates) = sample.cms_general_divisor_shape_candidates.as_deref() else {
        return Ok(None);
    };
    let Some(checks) = sample.cms_general_divisor_intersection_checks.as_deref() else {
        return Ok(None);
    };
    let mut values = Vec::new();
    for check in checks {
        if cms_general_divisor_intersection_check_status(check)
            != "cms_general_divisor_integral_solution_matches_inferred_degree"
        {
            continue;
        }
        let Some(candidate) = candidates
            .iter()
            .find(|candidate| candidate.shrinking_divisor_index == check.shrinking_divisor_index)
        else {
            return Err(format!(
                "CMS divisor check for shrinking divisor {} has no matching shape candidate",
                check.shrinking_divisor_index
            ));
        };
        if let Some(value) = candidate.toric_gv1_formula_value {
            values.push(value);
        }
    }
    values.sort_unstable();
    values.dedup();
    match values.as_slice() {
        [] => Ok(None),
        [value] => Ok(Some(value.to_string())),
        _ => Err(format!(
            "source-derived CMS checks produce conflicting GV formula values {values:?}"
        )),
    }
}

struct ValidatedContext<'a> {
    dimension: usize,
    degree_bound: i128,
    q_cols: usize,
    grading: &'a [i64],
    q_matrix: &'a [Vec<i64>],
    degree_bounded_rays: &'a [Vec<i64>],
    degree_bounded_ray_context: Option<&'a [DegreeBoundedMoriRayContextSample]>,
    covered_toric_gv_by_basis: HashMap<Vec<i64>, String>,
    source_derived_gv_by_basis: HashMap<Vec<i64>, String>,
    intersection: Intersection,
    stats: &'a MissingGvTargetStats,
    uncovered_source_ray_stats: Option<&'a MissingGvTargetStats>,
    shared_facet_unresolved_source_ray_stats: Option<&'a MissingGvTargetStats>,
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
    run_path_support_generators: bool,
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
    closure_generation_limit: Option<usize>,
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
    let (cms_general_divisor_solution_summaries, cms_general_divisor_solution_summary_error) =
        match cms_general_divisor_solution_summaries(Some(sample), &context.intersection) {
            Ok(summaries) => (summaries, None),
            Err(error) => (Vec::new(), Some(error)),
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
                origin_circuit_witnesses: sample.origin_circuit_witnesses.clone(),
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape: None,
                local_cygv_input_skeleton,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                cms_general_divisor_solution_summaries,
                cms_general_divisor_solution_summary_error,
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
                active_support_face_certificate_status: None,
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
                origin_circuit_witnesses: sample.origin_circuit_witnesses.clone(),
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape,
                local_cygv_input_skeleton,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                cms_general_divisor_solution_summaries,
                cms_general_divisor_solution_summary_error,
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
                active_support_face_certificate_status: None,
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
                origin_circuit_witnesses: sample.origin_circuit_witnesses.clone(),
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape,
                local_cygv_input_skeleton,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                cms_general_divisor_solution_summaries,
                cms_general_divisor_solution_summary_error,
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
                active_support_face_certificate_status: None,
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
                origin_circuit_witnesses: sample.origin_circuit_witnesses.clone(),
                origin_circuit_affine_support,
                local_cygv_hypersurface_shape,
                local_cygv_input_skeleton,
                cms_general_divisor_shape_candidates: sample
                    .cms_general_divisor_shape_candidates
                    .clone(),
                cms_general_divisor_intersection_checks: sample
                    .cms_general_divisor_intersection_checks
                    .clone(),
                cms_general_divisor_solution_summaries,
                cms_general_divisor_solution_summary_error,
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
                active_support_face_certificate_status: None,
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
        origin_circuit_witnesses: sample.origin_circuit_witnesses.clone(),
        origin_circuit_affine_support,
        local_cygv_hypersurface_shape,
        local_cygv_input_skeleton,
        cms_general_divisor_shape_candidates: sample.cms_general_divisor_shape_candidates.clone(),
        cms_general_divisor_intersection_checks: sample
            .cms_general_divisor_intersection_checks
            .clone(),
        cms_general_divisor_solution_summaries,
        cms_general_divisor_solution_summary_error,
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
        active_support_face_certificate_status: Some(active_support_face_certificate_status(
            sample, context,
        )),
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
                run_path_support_generators,
                element_limit,
                semigroup_measure_max_seed_count,
                closure_generation_limit,
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
    generation_counts: Vec<CygvClosureGenerationCount>,
    completed: bool,
}

fn cygv_path_history_probe(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    target: &[i64],
    run_lower_seed_diamonds: bool,
    run_path_support_generators: bool,
    element_limit: usize,
    max_seed_count: Option<usize>,
    closure_generation_limit: Option<usize>,
) -> CygvPathHistoryProbe {
    match cygv_path_history_probe_inner(
        sample,
        context,
        target,
        run_lower_seed_diamonds,
        run_path_support_generators,
        element_limit,
        max_seed_count,
        closure_generation_limit,
    ) {
        Ok(probe) => probe,
        Err(error) => CygvPathHistoryProbe {
            status: format!("error: {error}"),
            seed_count: None,
            reduced_seed_count: None,
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            closure_generation_counts: Vec::new(),
            target_in_closure: None,
            previous_level_count: cygv_previous_level_count(context.dimension),
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_counts_complete: false,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            predecessor_toric_coverage_counts: BTreeMap::new(),
            predecessor_known_qn_history_counts: BTreeMap::new(),
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
            closest_known_qn_predecessor: None,
            closest_known_qn_residual_predecessor: None,
            predecessor_candidate_sample_limit: CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT,
            predecessor_candidate_sample: Vec::new(),
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
            pair_expanded_lower_seed_decomposition_status: None,
            pair_expanded_lower_seed_decomposition_term_count: None,
            pair_expanded_lower_seed_decomposition_terms_nonzero: None,
            pair_expanded_lower_seed_diamond_status: None,
            pair_expanded_lower_seed_diamond_element_count: None,
            pair_expanded_lower_seed_diamond_gv: None,
            pair_expanded_lower_seed_diamond_error: None,
            path_support_size: None,
            path_support_generator_count: None,
            path_support_status: None,
            path_support_gv: None,
            path_support_error: None,
            path_support_lookup_status_counts: BTreeMap::new(),
            path_support_source_class_status_counts: BTreeMap::new(),
            path_support_lookup_sample: Vec::new(),
        },
    }
}

fn cygv_path_history_probe_inner(
    sample: &MissingGvTargetSample,
    context: &ValidatedContext<'_>,
    target: &[i64],
    run_lower_seed_diamonds: bool,
    run_path_support_generators: bool,
    element_limit: usize,
    max_seed_count: Option<usize>,
    closure_generation_limit: Option<usize>,
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
            closure_generation_counts: Vec::new(),
            target_in_closure: None,
            previous_level_count,
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_counts_complete: false,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            predecessor_toric_coverage_counts: BTreeMap::new(),
            predecessor_known_qn_history_counts: BTreeMap::new(),
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
            closest_known_qn_predecessor: None,
            closest_known_qn_residual_predecessor: None,
            predecessor_candidate_sample_limit: CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT,
            predecessor_candidate_sample: Vec::new(),
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: "skipped_empty_seed_set".to_string(),
            lower_seed_decomposition_term_count: None,
            lower_seed_decomposition_terms_nonzero: None,
            lower_seed_decomposition_error: None,
            lower_seed_diamond_status: None,
            lower_seed_diamond_element_count: None,
            lower_seed_diamond_gv: None,
            lower_seed_diamond_error: None,
            pair_expanded_lower_seed_decomposition_status: None,
            pair_expanded_lower_seed_decomposition_term_count: None,
            pair_expanded_lower_seed_decomposition_terms_nonzero: None,
            pair_expanded_lower_seed_diamond_status: None,
            pair_expanded_lower_seed_diamond_element_count: None,
            pair_expanded_lower_seed_diamond_gv: None,
            pair_expanded_lower_seed_diamond_error: None,
            path_support_size: None,
            path_support_generator_count: None,
            path_support_status: None,
            path_support_gv: None,
            path_support_error: None,
            path_support_lookup_status_counts: BTreeMap::new(),
            path_support_source_class_status_counts: BTreeMap::new(),
            path_support_lookup_sample: Vec::new(),
        });
    }
    if max_seed_count.is_some_and(|limit| seeds.len() > limit) {
        return Ok(CygvPathHistoryProbe {
            status: "skipped_seed_limit".to_string(),
            seed_count: Some(seeds.len()),
            reduced_seed_count: None,
            closure_element_count: None,
            closure_degree_counts: BTreeMap::new(),
            closure_generation_counts: Vec::new(),
            target_in_closure: None,
            previous_level_count,
            previous_window_degrees: Vec::new(),
            previous_window_degree_count: None,
            previous_window_element_count: None,
            predecessor_counts_complete: false,
            predecessor_difference_count: None,
            improving_predecessor_difference_count: None,
            predecessor_toric_coverage_counts: BTreeMap::new(),
            predecessor_known_qn_history_counts: BTreeMap::new(),
            closest_series_distance: None,
            closest_series_predecessor_nonzero: None,
            closest_series_difference_nonzero: None,
            closest_known_qn_predecessor: None,
            closest_known_qn_residual_predecessor: None,
            predecessor_candidate_sample_limit: CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT,
            predecessor_candidate_sample: Vec::new(),
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: "skipped_seed_limit".to_string(),
            lower_seed_decomposition_term_count: None,
            lower_seed_decomposition_terms_nonzero: None,
            lower_seed_decomposition_error: None,
            lower_seed_diamond_status: Some("skipped_seed_limit".to_string()),
            lower_seed_diamond_element_count: None,
            lower_seed_diamond_gv: None,
            lower_seed_diamond_error: None,
            pair_expanded_lower_seed_decomposition_status: None,
            pair_expanded_lower_seed_decomposition_term_count: None,
            pair_expanded_lower_seed_decomposition_terms_nonzero: None,
            pair_expanded_lower_seed_diamond_status: None,
            pair_expanded_lower_seed_diamond_element_count: None,
            pair_expanded_lower_seed_diamond_gv: None,
            pair_expanded_lower_seed_diamond_error: None,
            path_support_size: None,
            path_support_generator_count: None,
            path_support_status: None,
            path_support_gv: None,
            path_support_error: None,
            path_support_lookup_status_counts: BTreeMap::new(),
            path_support_source_class_status_counts: BTreeMap::new(),
            path_support_lookup_sample: Vec::new(),
        });
    }
    let reduced_seeds = cygv_pair_reduced_seed_generators(&seeds)
        .map_err(|error| format!("cygv seed reduction failed: {error}"))?
        .into_iter()
        .collect::<HashSet<_>>();
    let reduced_seed_count = reduced_seeds.len();
    let lower_seed_decomposition = lower_seed_decomposition_probe(target, &seeds, 4);
    let lower_seed_diamond = lower_seed_diamond_probe(
        target,
        &lower_seed_decomposition,
        context,
        run_lower_seed_diamonds,
        element_limit,
    );
    let pair_expanded_lower_seed =
        pair_expanded_lower_seed_probe(&lower_seed_decomposition, &seen, &reduced_seeds, 16, 8);
    let pair_expanded_lower_seed_diamond = pair_expanded_lower_seed_diamond_probe(
        target,
        &pair_expanded_lower_seed,
        context,
        run_lower_seed_diamonds,
        element_limit,
    );

    let closure = bounded_cygv_semigroup_closure(
        &seeds,
        context.grading,
        sample.degree,
        element_limit,
        closure_generation_limit,
    )?;
    let target_in_closure = closure.elements.contains(target);
    let selected_degree_vec =
        cygv_previous_window_degrees(&closure.degree_counts, sample.degree, previous_level_count);
    let selected_degrees = selected_degree_vec.iter().copied().collect::<HashSet<_>>();
    let predecessor_stats = cygv_path_predecessor_stats(
        &closure.elements,
        context.grading,
        target,
        &selected_degrees,
        &context.covered_toric_gv_by_basis,
        &context.source_derived_gv_by_basis,
        &seen,
        &reduced_seeds,
    )?;
    let closest_known_qn_residual_predecessor = cygv_closest_known_qn_residual_predecessor(
        &closure.elements,
        &closure.degree_counts,
        context.grading,
        context.dimension,
        previous_level_count,
        predecessor_stats.closest_known_qn_predecessor.as_ref(),
        &context.covered_toric_gv_by_basis,
        &context.source_derived_gv_by_basis,
    )?;
    let path_support_generators = path_support_generator_probe(
        target,
        sample.degree,
        &predecessor_stats.candidate_sample,
        context,
        run_path_support_generators,
    );
    if !closure.completed {
        return Ok(CygvPathHistoryProbe {
            status: closure.status,
            seed_count: Some(seeds.len()),
            reduced_seed_count: Some(reduced_seed_count),
            closure_element_count: Some(closure.elements.len()),
            previous_window_degrees: selected_degree_vec,
            closure_degree_counts: closure.degree_counts,
            closure_generation_counts: closure.generation_counts,
            target_in_closure: Some(target_in_closure),
            previous_level_count,
            previous_window_degree_count: Some(selected_degrees.len()),
            previous_window_element_count: Some(predecessor_stats.previous_window_element_count),
            predecessor_counts_complete: false,
            predecessor_difference_count: Some(predecessor_stats.predecessor_difference_count),
            improving_predecessor_difference_count: Some(
                predecessor_stats.improving_predecessor_difference_count,
            ),
            predecessor_toric_coverage_counts: predecessor_stats.toric_coverage_counts,
            predecessor_known_qn_history_counts: predecessor_stats.known_qn_history_counts,
            closest_series_distance: Some(format!("{:.6}", predecessor_stats.closest_distance)),
            closest_series_predecessor_nonzero: predecessor_stats
                .closest_predecessor
                .as_deref()
                .map(sparse_from_dense),
            closest_series_difference_nonzero: predecessor_stats
                .closest_difference
                .as_deref()
                .map(sparse_from_dense),
            closest_known_qn_predecessor: predecessor_stats.closest_known_qn_predecessor,
            closest_known_qn_residual_predecessor,
            predecessor_candidate_sample_limit: CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT,
            predecessor_candidate_sample: predecessor_stats.candidate_sample,
            lower_seed_decomposition_max_terms: 4,
            lower_seed_decomposition_status: lower_seed_decomposition.status,
            lower_seed_decomposition_term_count: lower_seed_decomposition.term_count,
            lower_seed_decomposition_terms_nonzero: lower_seed_decomposition.terms_nonzero,
            lower_seed_decomposition_error: lower_seed_decomposition.error,
            lower_seed_diamond_status: lower_seed_diamond.status,
            lower_seed_diamond_element_count: lower_seed_diamond.element_count,
            lower_seed_diamond_gv: lower_seed_diamond.gv,
            lower_seed_diamond_error: lower_seed_diamond.error,
            pair_expanded_lower_seed_decomposition_status: Some(pair_expanded_lower_seed.status),
            pair_expanded_lower_seed_decomposition_term_count: pair_expanded_lower_seed.term_count,
            pair_expanded_lower_seed_decomposition_terms_nonzero: pair_expanded_lower_seed
                .terms_nonzero,
            pair_expanded_lower_seed_diamond_status: pair_expanded_lower_seed_diamond.status,
            pair_expanded_lower_seed_diamond_element_count: pair_expanded_lower_seed_diamond
                .element_count,
            pair_expanded_lower_seed_diamond_gv: pair_expanded_lower_seed_diamond.gv,
            pair_expanded_lower_seed_diamond_error: pair_expanded_lower_seed_diamond.error,
            path_support_size: path_support_generators.support_size,
            path_support_generator_count: path_support_generators.generator_count,
            path_support_status: path_support_generators.status,
            path_support_gv: path_support_generators.gv,
            path_support_error: path_support_generators.error,
            path_support_lookup_status_counts: path_support_generators.lookup_status_counts,
            path_support_source_class_status_counts: path_support_generators
                .source_class_status_counts,
            path_support_lookup_sample: path_support_generators.lookup_sample,
        });
    }

    Ok(CygvPathHistoryProbe {
        status: "completed_bounded_closure".to_string(),
        seed_count: Some(seeds.len()),
        reduced_seed_count: Some(reduced_seed_count),
        closure_element_count: Some(closure.elements.len()),
        closure_degree_counts: closure.degree_counts,
        closure_generation_counts: closure.generation_counts,
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
        predecessor_toric_coverage_counts: predecessor_stats.toric_coverage_counts,
        predecessor_known_qn_history_counts: predecessor_stats.known_qn_history_counts,
        closest_series_distance: Some(format!("{:.6}", predecessor_stats.closest_distance)),
        closest_series_predecessor_nonzero: predecessor_stats
            .closest_predecessor
            .as_deref()
            .map(sparse_from_dense),
        closest_series_difference_nonzero: predecessor_stats
            .closest_difference
            .as_deref()
            .map(sparse_from_dense),
        closest_known_qn_predecessor: predecessor_stats.closest_known_qn_predecessor,
        closest_known_qn_residual_predecessor,
        predecessor_candidate_sample_limit: CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT,
        predecessor_candidate_sample: predecessor_stats.candidate_sample,
        lower_seed_decomposition_max_terms: 4,
        lower_seed_decomposition_status: lower_seed_decomposition.status,
        lower_seed_decomposition_term_count: lower_seed_decomposition.term_count,
        lower_seed_decomposition_terms_nonzero: lower_seed_decomposition.terms_nonzero,
        lower_seed_decomposition_error: lower_seed_decomposition.error,
        lower_seed_diamond_status: lower_seed_diamond.status,
        lower_seed_diamond_element_count: lower_seed_diamond.element_count,
        lower_seed_diamond_gv: lower_seed_diamond.gv,
        lower_seed_diamond_error: lower_seed_diamond.error,
        pair_expanded_lower_seed_decomposition_status: Some(pair_expanded_lower_seed.status),
        pair_expanded_lower_seed_decomposition_term_count: pair_expanded_lower_seed.term_count,
        pair_expanded_lower_seed_decomposition_terms_nonzero: pair_expanded_lower_seed
            .terms_nonzero,
        pair_expanded_lower_seed_diamond_status: pair_expanded_lower_seed_diamond.status,
        pair_expanded_lower_seed_diamond_element_count: pair_expanded_lower_seed_diamond
            .element_count,
        pair_expanded_lower_seed_diamond_gv: pair_expanded_lower_seed_diamond.gv,
        pair_expanded_lower_seed_diamond_error: pair_expanded_lower_seed_diamond.error,
        path_support_size: path_support_generators.support_size,
        path_support_generator_count: path_support_generators.generator_count,
        path_support_status: path_support_generators.status,
        path_support_gv: path_support_generators.gv,
        path_support_error: path_support_generators.error,
        path_support_lookup_status_counts: path_support_generators.lookup_status_counts,
        path_support_source_class_status_counts: path_support_generators.source_class_status_counts,
        path_support_lookup_sample: path_support_generators.lookup_sample,
    })
}

fn known_qn_history_status(
    toric_gv: Option<&str>,
    source_derived_gv: Option<&str>,
) -> Result<&'static str, String> {
    let toric_parsed = toric_gv.map(parse_rational).transpose()?;
    let source_parsed = source_derived_gv.map(parse_rational).transpose()?;
    if let (Some(toric), Some(source)) = (&toric_parsed, &source_parsed)
        && toric != source
    {
        return Err(format!(
            "toric GV value {toric} conflicts with source-derived GV value {source}"
        ));
    }
    let zero = MalachiteRational::from(0);
    if let Some(toric) = toric_parsed {
        if toric == zero {
            Ok("known_zero_toric_gv")
        } else {
            Ok("known_nonzero_toric_gv")
        }
    } else if let Some(source) = source_parsed {
        if source == zero {
            Ok("known_zero_source_gv")
        } else {
            Ok("known_nonzero_source_gv")
        }
    } else {
        Ok("unknown_not_toric_covered")
    }
}

fn known_qn_history_pair_status(predecessor_status: &str, difference_status: &str) -> String {
    format!("predecessor_{predecessor_status}__difference_{difference_status}")
}

fn is_known_nonzero_qn_history_status(status: &str) -> bool {
    matches!(status, "known_nonzero_toric_gv" | "known_nonzero_source_gv")
}

fn cygv_closest_known_qn_residual_predecessor(
    elements: &HashSet<Vec<i64>>,
    degree_counts: &BTreeMap<i128, usize>,
    grading: &[i64],
    dimension: usize,
    previous_level_count: usize,
    closest_known_qn_predecessor: Option<&CygvClosestKnownQnPredecessor>,
    covered_toric_gv_by_basis: &HashMap<Vec<i64>, String>,
    source_derived_gv_by_basis: &HashMap<Vec<i64>, String>,
) -> Result<Option<CygvClosestKnownQnPredecessor>, String> {
    let Some(closest_known_qn_predecessor) = closest_known_qn_predecessor else {
        return Ok(None);
    };
    if closest_known_qn_predecessor.difference_known_qn_history_status
        != "unknown_not_toric_covered"
    {
        return Ok(None);
    }
    let residual = dense_from_sparse(&closest_known_qn_predecessor.difference_nonzero, dimension)?;
    let residual_degree = curve_degree(&residual, grading)?;
    let selected_degree_vec =
        cygv_previous_window_degrees(degree_counts, residual_degree, previous_level_count);
    let selected_degrees = selected_degree_vec.into_iter().collect::<HashSet<_>>();
    cygv_closest_known_qn_predecessor(
        elements,
        grading,
        &residual,
        &selected_degrees,
        covered_toric_gv_by_basis,
        source_derived_gv_by_basis,
    )
}

fn cygv_closest_known_qn_predecessor(
    elements: &HashSet<Vec<i64>>,
    grading: &[i64],
    target: &[i64],
    selected_degrees: &HashSet<i128>,
    covered_toric_gv_by_basis: &HashMap<Vec<i64>, String>,
    source_derived_gv_by_basis: &HashMap<Vec<i64>, String>,
) -> Result<Option<CygvClosestKnownQnPredecessor>, String> {
    let target_series_distance = cygv_series_distance(target);
    let mut closest = None;
    let mut sorted_elements = elements.iter().collect::<Vec<_>>();
    sorted_elements.sort();
    for element in sorted_elements {
        let degree = curve_degree(element, grading)?;
        if !selected_degrees.contains(&degree) {
            continue;
        }
        let difference = checked_vector_difference(target, element)?;
        if !elements.contains(&difference) {
            continue;
        }
        let predecessor_toric_gv = covered_toric_gv_by_basis.get(element);
        let predecessor_source_derived_gv = source_derived_gv_by_basis.get(element);
        let predecessor_known_qn_history_status = known_qn_history_status(
            predecessor_toric_gv.map(String::as_str),
            predecessor_source_derived_gv.map(String::as_str),
        )?;
        if !is_known_nonzero_qn_history_status(predecessor_known_qn_history_status) {
            continue;
        }
        let distance = cygv_series_distance(&difference);
        if distance >= target_series_distance
            || closest
                .as_ref()
                .is_some_and(|(best_distance, _)| distance >= *best_distance)
        {
            continue;
        }
        let difference_toric_gv = covered_toric_gv_by_basis.get(&difference);
        let difference_source_derived_gv = source_derived_gv_by_basis.get(&difference);
        let difference_known_qn_history_status = known_qn_history_status(
            difference_toric_gv.map(String::as_str),
            difference_source_derived_gv.map(String::as_str),
        )?;
        closest = Some((
            distance,
            CygvClosestKnownQnPredecessor {
                predecessor_degree: degree,
                difference_degree: curve_degree(&difference, grading)?,
                series_distance: format!("{distance:.6}"),
                predecessor_toric_gv: predecessor_toric_gv.cloned(),
                predecessor_source_derived_gv: predecessor_source_derived_gv.cloned(),
                predecessor_known_qn_history_status: predecessor_known_qn_history_status
                    .to_string(),
                difference_toric_gv: difference_toric_gv.cloned(),
                difference_source_derived_gv: difference_source_derived_gv.cloned(),
                difference_known_qn_history_status: difference_known_qn_history_status.to_string(),
                predecessor_nonzero: sparse_from_dense(element),
                difference_nonzero: sparse_from_dense(&difference),
            },
        ));
    }
    Ok(closest.map(|(_, predecessor)| predecessor))
}

fn cygv_path_predecessor_stats(
    elements: &HashSet<Vec<i64>>,
    grading: &[i64],
    target: &[i64],
    selected_degrees: &HashSet<i128>,
    covered_toric_gv_by_basis: &HashMap<Vec<i64>, String>,
    source_derived_gv_by_basis: &HashMap<Vec<i64>, String>,
    seed_set: &HashSet<Vec<i64>>,
    reduced_seed_set: &HashSet<Vec<i64>>,
) -> Result<CygvPathPredecessorStats, String> {
    let mut previous_window_element_count = 0usize;
    let mut predecessor_difference_count = 0usize;
    let mut improving_predecessor_difference_count = 0usize;
    let target_series_distance = cygv_series_distance(target);
    let mut closest_distance = target_series_distance;
    let mut closest_predecessor = None;
    let mut closest_difference = None;
    let mut candidate_sample = Vec::new();
    let mut toric_coverage_counts = BTreeMap::new();
    let mut known_qn_history_counts = BTreeMap::new();
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
        let predecessor_is_toric_covered = covered_toric_gv_by_basis.contains_key(element);
        let difference_is_toric_covered = covered_toric_gv_by_basis.contains_key(&difference);
        let predecessor_toric_gv = covered_toric_gv_by_basis.get(element);
        let difference_toric_gv = covered_toric_gv_by_basis.get(&difference);
        let predecessor_source_derived_gv = source_derived_gv_by_basis.get(element);
        let difference_source_derived_gv = source_derived_gv_by_basis.get(&difference);
        let predecessor_known_qn_history_status = known_qn_history_status(
            predecessor_toric_gv.map(String::as_str),
            predecessor_source_derived_gv.map(String::as_str),
        )?;
        let difference_known_qn_history_status = known_qn_history_status(
            difference_toric_gv.map(String::as_str),
            difference_source_derived_gv.map(String::as_str),
        )?;
        let known_qn_history_pair_status = known_qn_history_pair_status(
            predecessor_known_qn_history_status,
            difference_known_qn_history_status,
        );
        let coverage_key = match (predecessor_is_toric_covered, difference_is_toric_covered) {
            (true, true) => "both_toric_covered",
            (true, false) => "predecessor_only_toric_covered",
            (false, true) => "difference_only_toric_covered",
            (false, false) => "neither_toric_covered",
        };
        *toric_coverage_counts
            .entry(coverage_key.to_string())
            .or_insert(0) += 1;
        *known_qn_history_counts
            .entry(known_qn_history_pair_status.clone())
            .or_insert(0) += 1;
        let difference_degree = curve_degree(&difference, grading)?;
        let distance = cygv_series_distance(&difference);
        if candidate_sample.len() < CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT {
            candidate_sample.push((
                distance,
                CygvPathPredecessorCandidate {
                    predecessor_degree: degree,
                    difference_degree,
                    series_distance: format!("{distance:.6}"),
                    predecessor_toric_gv: predecessor_toric_gv.cloned(),
                    difference_toric_gv: difference_toric_gv.cloned(),
                    predecessor_source_derived_gv: predecessor_source_derived_gv.cloned(),
                    difference_source_derived_gv: difference_source_derived_gv.cloned(),
                    predecessor_known_qn_history_status: predecessor_known_qn_history_status
                        .to_string(),
                    difference_known_qn_history_status: difference_known_qn_history_status
                        .to_string(),
                    known_qn_history_pair_status: known_qn_history_pair_status.clone(),
                    predecessor_is_seed: seed_set.contains(element),
                    difference_is_seed: seed_set.contains(&difference),
                    predecessor_is_reduced_seed: reduced_seed_set.contains(element),
                    difference_is_reduced_seed: reduced_seed_set.contains(&difference),
                    predecessor_first_generation_seed_sum: first_generation_seed_sum_decomposition(
                        element,
                        grading,
                        seed_set,
                        reduced_seed_set,
                        covered_toric_gv_by_basis,
                        source_derived_gv_by_basis,
                    )?,
                    difference_first_generation_seed_sum: first_generation_seed_sum_decomposition(
                        &difference,
                        grading,
                        seed_set,
                        reduced_seed_set,
                        covered_toric_gv_by_basis,
                        source_derived_gv_by_basis,
                    )?,
                    predecessor_nonzero: sparse_from_dense(element),
                    difference_nonzero: sparse_from_dense(&difference),
                },
            ));
            candidate_sample.sort_by(|(lhs_distance, lhs), (rhs_distance, rhs)| {
                lhs_distance
                    .total_cmp(rhs_distance)
                    .then_with(|| lhs.predecessor_nonzero.cmp(&rhs.predecessor_nonzero))
                    .then_with(|| lhs.difference_nonzero.cmp(&rhs.difference_nonzero))
            });
        } else if candidate_sample
            .last()
            .is_some_and(|(worst_distance, _)| distance < *worst_distance)
        {
            candidate_sample.pop();
            candidate_sample.push((
                distance,
                CygvPathPredecessorCandidate {
                    predecessor_degree: degree,
                    difference_degree,
                    series_distance: format!("{distance:.6}"),
                    predecessor_toric_gv: predecessor_toric_gv.cloned(),
                    difference_toric_gv: difference_toric_gv.cloned(),
                    predecessor_source_derived_gv: predecessor_source_derived_gv.cloned(),
                    difference_source_derived_gv: difference_source_derived_gv.cloned(),
                    predecessor_known_qn_history_status: predecessor_known_qn_history_status
                        .to_string(),
                    difference_known_qn_history_status: difference_known_qn_history_status
                        .to_string(),
                    known_qn_history_pair_status: known_qn_history_pair_status.clone(),
                    predecessor_is_seed: seed_set.contains(element),
                    difference_is_seed: seed_set.contains(&difference),
                    predecessor_is_reduced_seed: reduced_seed_set.contains(element),
                    difference_is_reduced_seed: reduced_seed_set.contains(&difference),
                    predecessor_first_generation_seed_sum: first_generation_seed_sum_decomposition(
                        element,
                        grading,
                        seed_set,
                        reduced_seed_set,
                        covered_toric_gv_by_basis,
                        source_derived_gv_by_basis,
                    )?,
                    difference_first_generation_seed_sum: first_generation_seed_sum_decomposition(
                        &difference,
                        grading,
                        seed_set,
                        reduced_seed_set,
                        covered_toric_gv_by_basis,
                        source_derived_gv_by_basis,
                    )?,
                    predecessor_nonzero: sparse_from_dense(element),
                    difference_nonzero: sparse_from_dense(&difference),
                },
            ));
            candidate_sample.sort_by(|(lhs_distance, lhs), (rhs_distance, rhs)| {
                lhs_distance
                    .total_cmp(rhs_distance)
                    .then_with(|| lhs.predecessor_nonzero.cmp(&rhs.predecessor_nonzero))
                    .then_with(|| lhs.difference_nonzero.cmp(&rhs.difference_nonzero))
            });
        }
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
        toric_coverage_counts,
        known_qn_history_counts,
        closest_distance,
        closest_predecessor,
        closest_difference,
        closest_known_qn_predecessor: cygv_closest_known_qn_predecessor(
            elements,
            grading,
            target,
            selected_degrees,
            covered_toric_gv_by_basis,
            source_derived_gv_by_basis,
        )?,
        candidate_sample: candidate_sample
            .into_iter()
            .map(|(_, candidate)| candidate)
            .collect(),
    })
}

fn first_generation_seed_sum_decomposition(
    target: &[i64],
    grading: &[i64],
    seed_set: &HashSet<Vec<i64>>,
    reduced_seed_set: &HashSet<Vec<i64>>,
    covered_toric_gv_by_basis: &HashMap<Vec<i64>, String>,
    source_derived_gv_by_basis: &HashMap<Vec<i64>, String>,
) -> Result<Option<CygvSeedSumDecomposition>, String> {
    let mut reduced_seeds = reduced_seed_set.iter().collect::<Vec<_>>();
    reduced_seeds.sort();
    for reduced_seed in reduced_seeds {
        let seed = checked_vector_difference(target, reduced_seed)?;
        if !seed_set.contains(&seed) {
            continue;
        }
        let reduced_seed_toric_gv = covered_toric_gv_by_basis.get(reduced_seed);
        let seed_toric_gv = covered_toric_gv_by_basis.get(&seed);
        let reduced_seed_source_derived_gv = source_derived_gv_by_basis.get(reduced_seed);
        let seed_source_derived_gv = source_derived_gv_by_basis.get(&seed);
        return Ok(Some(CygvSeedSumDecomposition {
            reduced_seed_degree: curve_degree(reduced_seed, grading)?,
            seed_degree: curve_degree(&seed, grading)?,
            reduced_seed_toric_gv: reduced_seed_toric_gv.cloned(),
            seed_toric_gv: seed_toric_gv.cloned(),
            reduced_seed_source_derived_gv: reduced_seed_source_derived_gv.cloned(),
            seed_source_derived_gv: seed_source_derived_gv.cloned(),
            reduced_seed_known_qn_history_status: known_qn_history_status(
                reduced_seed_toric_gv.map(String::as_str),
                reduced_seed_source_derived_gv.map(String::as_str),
            )?
            .to_string(),
            seed_known_qn_history_status: known_qn_history_status(
                seed_toric_gv.map(String::as_str),
                seed_source_derived_gv.map(String::as_str),
            )?
            .to_string(),
            seed_is_reduced_seed: reduced_seed_set.contains(&seed),
            seed_pair_reduction_sum: seed_pair_reduction_sum_decomposition(
                &seed,
                grading,
                seed_set,
                reduced_seed_set,
                covered_toric_gv_by_basis,
                source_derived_gv_by_basis,
            )?,
            reduced_seed_nonzero: sparse_from_dense(reduced_seed),
            seed_nonzero: sparse_from_dense(&seed),
        }));
    }
    Ok(None)
}

fn seed_pair_reduction_sum_decomposition(
    seed: &[i64],
    grading: &[i64],
    seed_set: &HashSet<Vec<i64>>,
    reduced_seed_set: &HashSet<Vec<i64>>,
    covered_toric_gv_by_basis: &HashMap<Vec<i64>, String>,
    source_derived_gv_by_basis: &HashMap<Vec<i64>, String>,
) -> Result<Option<CygvSeedPairDecomposition>, String> {
    if reduced_seed_set.contains(seed) {
        return Ok(None);
    }
    let mut sorted_seeds = seed_set.iter().collect::<Vec<_>>();
    sorted_seeds.sort();
    for lhs in sorted_seeds {
        let rhs = checked_vector_difference(seed, lhs)?;
        if !seed_set.contains(&rhs) {
            continue;
        }
        let lhs_toric_gv = covered_toric_gv_by_basis.get(lhs);
        let rhs_toric_gv = covered_toric_gv_by_basis.get(&rhs);
        let lhs_source_derived_gv = source_derived_gv_by_basis.get(lhs);
        let rhs_source_derived_gv = source_derived_gv_by_basis.get(&rhs);
        return Ok(Some(CygvSeedPairDecomposition {
            lhs_degree: curve_degree(lhs, grading)?,
            rhs_degree: curve_degree(&rhs, grading)?,
            lhs_toric_gv: lhs_toric_gv.cloned(),
            rhs_toric_gv: rhs_toric_gv.cloned(),
            lhs_source_derived_gv: lhs_source_derived_gv.cloned(),
            rhs_source_derived_gv: rhs_source_derived_gv.cloned(),
            lhs_known_qn_history_status: known_qn_history_status(
                lhs_toric_gv.map(String::as_str),
                lhs_source_derived_gv.map(String::as_str),
            )?
            .to_string(),
            rhs_known_qn_history_status: known_qn_history_status(
                rhs_toric_gv.map(String::as_str),
                rhs_source_derived_gv.map(String::as_str),
            )?
            .to_string(),
            lhs_is_reduced_seed: reduced_seed_set.contains(lhs),
            rhs_is_reduced_seed: reduced_seed_set.contains(&rhs),
            lhs_nonzero: sparse_from_dense(lhs),
            rhs_nonzero: sparse_from_dense(&rhs),
        }));
    }
    Ok(None)
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
    generation_limit: Option<usize>,
) -> Result<BoundedCygvClosure, String> {
    if element_limit == 0 {
        return Err("bounded cygv closure element limit must be positive".to_string());
    }
    if generation_limit == Some(0) {
        return Err("bounded cygv closure generation limit must be positive".to_string());
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
    let mut generation_counts = Vec::new();
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
            generation_counts,
            completed: false,
        });
    }

    let mut generation = 0usize;
    loop {
        generation += 1;
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
                generation_counts,
                completed: true,
            });
        }
        let mut sorted_new_elements = new_elements.into_iter().collect::<Vec<_>>();
        sorted_new_elements.sort();
        let total_after_full_generation = elements
            .len()
            .checked_add(sorted_new_elements.len())
            .ok_or_else(|| "cygv closure generation count overflowed usize".to_string())?;
        if total_after_full_generation > element_limit {
            generation_counts.push(CygvClosureGenerationCount {
                generation,
                starting_element_count: starting_elements.len(),
                new_element_count: sorted_new_elements.len(),
                total_element_count_after_full_generation: total_after_full_generation,
                truncated_at_limit: true,
            });
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
                generation_counts,
                completed: false,
            });
        }
        generation_counts.push(CygvClosureGenerationCount {
            generation,
            starting_element_count: starting_elements.len(),
            new_element_count: sorted_new_elements.len(),
            total_element_count_after_full_generation: total_after_full_generation,
            truncated_at_limit: false,
        });
        for element in &sorted_new_elements {
            elements.insert(element.clone());
        }
        starting_elements = sorted_new_elements.into_iter().collect();
        if generation_limit.is_some_and(|limit| generation >= limit) {
            let degree_counts = cygv_closure_degree_counts(&elements, grading_vector)?;
            return Ok(BoundedCygvClosure {
                status: format!("stopped_generation_limit_{generation}"),
                elements,
                degree_counts,
                generation_counts,
                completed: false,
            });
        }
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
    certify_origin_witness_domains: bool,
    origin_support_certificate_limit: usize,
    certify_target_extremal_rays: bool,
    target_extremal_generator_limit: usize,
    target_extremal_max_degree: Option<i128>,
    measure_cygv_semigroups: bool,
    probe_cygv_path_history: bool,
    run_lower_seed_diamonds: bool,
    run_path_support_generators: bool,
    measure_cygv_degree_ladder: bool,
    cygv_degree_ladder_max_degree: Option<i128>,
    semigroup_measure_max_target_degree: Option<i128>,
    semigroup_measure_max_seed_count: Option<usize>,
    scan_local_integer_tensors: bool,
    local_tensor_scan_bound: i64,
    element_limit: usize,
    closure_generation_limit: Option<usize>,
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
            run_path_support_generators,
            measure_cygv_degree_ladder,
            cygv_degree_ladder_max_degree,
            semigroup_measure_max_target_degree,
            semigroup_measure_max_seed_count,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            element_limit,
            closure_generation_limit,
        ));
    }
    let missing_local_cygv_charge_signature_counts =
        local_cygv_charge_signature_counts(&validated.stats.sample, target_index_filter);
    let local_cygv_target_candidate_status_counts = local_cygv_target_candidate_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let missing_local_cygv_actual_call_readiness_counts = local_cygv_actual_call_readiness_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let missing_local_cygv_missing_source_input_counts = local_cygv_missing_source_input_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let missing_cms_general_divisor_candidate_status_counts =
        cms_general_divisor_candidate_status_counts(&targets);
    let missing_cms_general_divisor_intersection_check_status_counts =
        cms_general_divisor_intersection_check_status_counts(&targets);
    let missing_origin_circuit_witness_relation_status_counts =
        origin_circuit_witness_relation_status_counts(&validated.stats.sample, target_index_filter);
    let uncovered_source_ray_origin_circuit_witness_relation_status_counts = validated
        .uncovered_source_ray_stats
        .map(|stats| origin_circuit_witness_relation_status_counts(&stats.sample, None))
        .unwrap_or_default();
    let shared_facet_unresolved_source_ray_target_count = validated
        .shared_facet_unresolved_source_ray_stats
        .map(|stats| stats.target_count);
    let shared_facet_unresolved_source_ray_origin_circuit_pattern_counts = validated
        .shared_facet_unresolved_source_ray_stats
        .map(|stats| origin_circuit_pattern_counts(&stats.sample))
        .unwrap_or_default();
    let shared_facet_unresolved_source_ray_origin_circuit_witness_relation_status_counts =
        validated
            .shared_facet_unresolved_source_ray_stats
            .map(|stats| origin_circuit_witness_relation_status_counts(&stats.sample, None))
            .unwrap_or_default();
    let shared_facet_unresolved_source_ray_local_cygv_charge_signature_counts = validated
        .shared_facet_unresolved_source_ray_stats
        .map(|stats| local_cygv_charge_signature_counts(&stats.sample, None))
        .unwrap_or_default();
    let mut shared_facet_unresolved_source_ray_sample = Vec::new();
    if let Some(stats) = validated.shared_facet_unresolved_source_ray_stats {
        for (idx, sample) in stats.sample.iter().enumerate() {
            shared_facet_unresolved_source_ray_sample.push(report_target(
                idx,
                sample,
                validated,
                false,
                false,
                None,
                None,
                false,
                false,
                origin_support_certificate_limit,
                false,
                target_extremal_generator_limit,
                None,
                false,
                false,
                false,
                false,
                false,
                None,
                None,
                None,
                &mut semigroup_measurement_cache,
                &mut semigroup_ladder_cache,
                element_limit,
                closure_generation_limit,
            ));
        }
    }
    let shared_facet_unresolved_source_ray_local_cygv_actual_call_readiness_counts =
        local_cygv_actual_call_readiness_counts(
            shared_facet_unresolved_source_ray_sample
                .iter()
                .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
        );
    let shared_facet_unresolved_source_ray_local_cygv_missing_source_input_counts =
        local_cygv_missing_source_input_counts(
            shared_facet_unresolved_source_ray_sample
                .iter()
                .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
        );
    let shared_facet_unresolved_source_ray_cms_general_divisor_candidate_status_counts =
        cms_general_divisor_candidate_status_counts(&shared_facet_unresolved_source_ray_sample);
    let shared_facet_unresolved_source_ray_cms_general_divisor_intersection_check_status_counts =
        cms_general_divisor_intersection_check_status_counts(
            &shared_facet_unresolved_source_ray_sample,
        );
    let shared_facet_unresolved_source_ray_cms_solution_status_counts =
        target_cms_solution_status_counts(&shared_facet_unresolved_source_ray_sample, |summary| {
            summary.cms_check_status.as_str()
        });
    let shared_facet_unresolved_source_ray_cms_intersection_tensor_status_counts =
        target_cms_solution_status_counts(&shared_facet_unresolved_source_ray_sample, |summary| {
            summary.local_intersection_tensor_candidate_status.as_str()
        });
    let shared_facet_unresolved_source_ray_cms_primitive_probe_status_counts =
        target_cms_solution_status_counts(&shared_facet_unresolved_source_ray_sample, |summary| {
            summary
                .local_cygv_primitive_probe
                .as_ref()
                .map_or("local_cygv_primitive_probe_not_run", |probe| {
                    probe.status.as_str()
                })
        });
    let shared_facet_unresolved_source_ray_cms_primitive_probe_gv_counts =
        target_cms_primitive_probe_gv_counts(&shared_facet_unresolved_source_ray_sample);
    let shared_facet_unresolved_source_ray_cms_solution_summary_error_counts =
        target_cms_solution_summary_error_counts(&shared_facet_unresolved_source_ray_sample);
    let shared_facet_unresolved_source_ray_unit_phase_probe_sample =
        local_cygv_target_unit_phase_probe_summaries(&shared_facet_unresolved_source_ray_sample);
    let shared_facet_unresolved_source_ray_unit_phase_probe_status_counts =
        local_cygv_target_unit_phase_probe_status_counts(
            &shared_facet_unresolved_source_ray_unit_phase_probe_sample,
        );
    let shared_facet_unresolved_source_ray_origin_omitted_unit_phase_probe_status_counts =
        local_cygv_target_origin_omitted_unit_phase_probe_status_counts(
            &shared_facet_unresolved_source_ray_unit_phase_probe_sample,
        );
    let shared_facet_unresolved_source_ray_unit_effective_tensor_requirement_status_counts =
        local_cygv_target_unit_effective_tensor_requirement_status_counts(
            &shared_facet_unresolved_source_ray_unit_phase_probe_sample,
        );
    let shared_facet_unresolved_source_ray_origin_omitted_effective_tensor_requirement_status_counts =
        local_cygv_target_origin_omitted_effective_tensor_requirement_status_counts(
            &shared_facet_unresolved_source_ray_unit_phase_probe_sample,
        );
    let origin_circuit_facet_context_status_counts =
        origin_circuit_facet_context_status_counts(&validated.stats.sample, target_index_filter);
    let origin_circuit_witness_domain_sample = origin_circuit_witness_domain_summaries(
        &validated.stats.sample,
        validated,
        target_index_filter,
        certify_origin_witness_domains,
        origin_support_certificate_limit,
    );
    let origin_circuit_witness_relation_support_face_profile_counts =
        origin_circuit_witness_domain_status_counts(
            &origin_circuit_witness_domain_sample,
            |summary| &summary.relation_support_face_profile,
        );
    let origin_circuit_witness_shared_facet_face_profile_counts =
        origin_circuit_witness_domain_status_counts(
            &origin_circuit_witness_domain_sample,
            |summary| &summary.shared_facet_face_profile,
        );
    let origin_circuit_witness_facet_union_face_profile_counts =
        origin_circuit_witness_domain_status_counts(
            &origin_circuit_witness_domain_sample,
            |summary| &summary.facet_union_face_profile,
        );
    let origin_circuit_witness_relation_support_face_certificate_status_counts =
        origin_circuit_witness_domain_status_counts(
            &origin_circuit_witness_domain_sample,
            |summary| &summary.relation_support_face_certificate_status,
        );
    let origin_circuit_witness_shared_facet_face_certificate_status_counts =
        origin_circuit_witness_domain_status_counts(
            &origin_circuit_witness_domain_sample,
            |summary| &summary.shared_facet_face_certificate_status,
        );
    let origin_circuit_witness_facet_union_face_certificate_status_counts =
        origin_circuit_witness_domain_status_counts(
            &origin_circuit_witness_domain_sample,
            |summary| &summary.facet_union_face_certificate_status,
        );
    let origin_circuit_witness_domain_unresolved_generators =
        origin_circuit_witness_domain_unresolved_generator_summaries(
            &validated.stats.sample,
            validated,
            target_index_filter,
            None,
        );
    let origin_circuit_witness_domain_unresolved_generator_unique_count =
        origin_circuit_witness_domain_unresolved_generators.len();
    let origin_circuit_witness_domain_unresolved_generator_occurrence_count =
        origin_circuit_witness_domain_unresolved_generators
            .iter()
            .map(|summary| summary.occurrence_count)
            .sum();
    let origin_circuit_witness_domain_unresolved_generator_status_counts =
        origin_circuit_witness_domain_unresolved_generator_status_counts_for(
            &origin_circuit_witness_domain_unresolved_generators,
        );
    let origin_circuit_witness_domain_unresolved_generator_degree_counts =
        origin_circuit_witness_domain_unresolved_generator_degree_counts_for(
            &origin_circuit_witness_domain_unresolved_generators,
        );
    let origin_circuit_witness_domain_unresolved_generator_sample =
        origin_circuit_witness_domain_unresolved_generators
            .into_iter()
            .take(ORIGIN_CIRCUIT_WITNESS_DOMAIN_UNRESOLVED_SAMPLE_LIMIT)
            .collect::<Vec<_>>();
    let origin_circuit_witness_shared_facet_unresolved_generators =
        origin_circuit_witness_domain_unresolved_generator_summaries(
            &validated.stats.sample,
            validated,
            target_index_filter,
            Some("shared_facet"),
        );
    let origin_circuit_witness_shared_facet_unresolved_generator_unique_count =
        origin_circuit_witness_shared_facet_unresolved_generators.len();
    let origin_circuit_witness_shared_facet_unresolved_generator_occurrence_count =
        origin_circuit_witness_shared_facet_unresolved_generators
            .iter()
            .map(|summary| summary.occurrence_count)
            .sum();
    let origin_circuit_witness_shared_facet_unresolved_generator_status_counts =
        origin_circuit_witness_domain_unresolved_generator_status_counts_for(
            &origin_circuit_witness_shared_facet_unresolved_generators,
        );
    let origin_circuit_witness_shared_facet_unresolved_generator_degree_counts =
        origin_circuit_witness_domain_unresolved_generator_degree_counts_for(
            &origin_circuit_witness_shared_facet_unresolved_generators,
        );
    let origin_circuit_witness_shared_facet_unresolved_generator_sample =
        origin_circuit_witness_shared_facet_unresolved_generators
            .into_iter()
            .take(ORIGIN_CIRCUIT_WITNESS_DOMAIN_UNRESOLVED_SAMPLE_LIMIT)
            .collect::<Vec<_>>();
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
    let cygv_path_history_status_counts = optional_status_counts(
        targets.iter().map(|target| {
            target
                .cygv_path_history_probe
                .as_ref()
                .map(|probe| probe.status.as_str())
        }),
        "not_run",
    );
    let cygv_lower_seed_decomposition_status_counts = optional_status_counts(
        targets.iter().map(|target| {
            target
                .cygv_path_history_probe
                .as_ref()
                .map(|probe| probe.lower_seed_decomposition_status.as_str())
        }),
        "not_run",
    );
    let cygv_lower_seed_diamond_status_counts = optional_status_counts(
        targets.iter().map(|target| {
            target
                .cygv_path_history_probe
                .as_ref()
                .and_then(|probe| probe.lower_seed_diamond_status.as_deref())
        }),
        "not_run",
    );
    let path_support_uncovered_source_ray_sample =
        path_support_uncovered_source_ray_summaries(&targets);
    let path_support_uncovered_source_ray_unique_count =
        path_support_uncovered_source_ray_sample.len();
    let path_support_uncovered_source_ray_occurrence_count =
        path_support_uncovered_source_ray_sample
            .iter()
            .map(|summary| summary.occurrence_count)
            .sum();
    let path_support_uncovered_source_ray_degree_counts =
        path_support_uncovered_source_ray_degree_counts(&path_support_uncovered_source_ray_sample);
    let path_support_uncovered_source_ray_lookup_status_counts =
        path_support_uncovered_source_ray_lookup_status_counts(
            &path_support_uncovered_source_ray_sample,
        );
    let path_support_uncovered_source_ray_path_support_gv_counts =
        path_support_uncovered_source_ray_path_support_gv_counts(
            &path_support_uncovered_source_ray_sample,
        );
    let path_support_uncovered_source_ray_local_cygv_primitive_probe_status_counts =
        path_support_uncovered_source_ray_local_cygv_primitive_probe_status_counts(
            &path_support_uncovered_source_ray_sample,
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
    let local_cygv_q_matrix_phase_status_counts = local_cygv_q_matrix_phase_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cygv_origin_point_status_counts = local_cygv_origin_point_status_counts(
        targets
            .iter()
            .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
    );
    let local_cygv_origin_omitted_compact_shape_status_counts =
        local_cygv_origin_omitted_compact_shape_status_counts(
            targets
                .iter()
                .filter_map(|target| target.local_cygv_input_skeleton.as_ref()),
        );
    let local_cygv_target_unit_phase_probe_sample =
        local_cygv_target_unit_phase_probe_summaries(&targets);
    let local_cygv_target_unit_phase_probe_status_counts =
        local_cygv_target_unit_phase_probe_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_origin_omitted_unit_phase_probe_status_counts =
        local_cygv_target_origin_omitted_unit_phase_probe_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_unit_effective_tensor_requirement_status_counts =
        local_cygv_target_unit_effective_tensor_requirement_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_origin_omitted_effective_tensor_requirement_status_counts =
        local_cygv_target_origin_omitted_effective_tensor_requirement_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_unit_formula_sum_probe_status_counts =
        local_cygv_target_unit_formula_sum_probe_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_origin_omitted_unit_formula_sum_probe_status_counts =
        local_cygv_target_origin_omitted_unit_formula_sum_probe_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_unit_formula_sum_effective_tensor_requirement_status_counts =
        local_cygv_target_unit_formula_sum_effective_tensor_requirement_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_origin_omitted_unit_formula_sum_effective_tensor_requirement_status_counts =
        local_cygv_target_origin_omitted_unit_formula_sum_effective_tensor_requirement_status_counts(
            &local_cygv_target_unit_phase_probe_sample,
        );
    let local_cygv_target_integer_tensor_scan_sample = local_cygv_integer_tensor_scan_summaries(
        &targets,
        scan_local_integer_tensors,
        local_tensor_scan_bound,
    );
    let local_cygv_target_integer_tensor_scan_status_counts =
        local_cygv_integer_tensor_scan_status_counts(
            &local_cygv_target_integer_tensor_scan_sample,
            targets.len(),
            scan_local_integer_tensors,
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
    let target_status_counts = optional_status_counts(
        targets.iter().map(|target| Some(target.status.as_str())),
        "missing",
    );
    let active_decomposition_generator_source_status_counts =
        active_decomposition_generator_source_status_counts(
            &validated.stats.sample,
            validated,
            target_index_filter,
        );
    let active_decomposition_unresolved_source_leaf_sample =
        active_decomposition_unresolved_source_leaf_summaries(
            &validated.stats.sample,
            validated,
            target_index_filter,
        );
    let active_decomposition_source_leaf_unit_phase_probe_status_counts =
        active_decomposition_source_leaf_unit_phase_probe_status_counts(
            &active_decomposition_unresolved_source_leaf_sample,
        );
    let active_decomposition_source_leaf_origin_omitted_unit_phase_probe_status_counts =
        active_decomposition_source_leaf_origin_omitted_unit_phase_probe_status_counts(
            &active_decomposition_unresolved_source_leaf_sample,
        );
    let active_decomposition_source_leaf_cms_solution_status_counts =
        active_decomposition_source_leaf_cms_solution_status_counts(
            &active_decomposition_unresolved_source_leaf_sample,
        );
    let active_decomposition_source_leaf_cms_candidate_status_counts =
        active_decomposition_source_leaf_cms_candidate_status_counts(
            &active_decomposition_unresolved_source_leaf_sample,
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
        covered_toric_gv_context_count: context
            .covered_toric_gv_context_for_missing
            .as_ref()
            .map(Vec::len),
        degree_bounded_toric_gv_diagnostic_context_count: context
            .degree_bounded_toric_gv_diagnostic_context_for_missing
            .as_ref()
            .map(Vec::len),
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
        target_status_counts,
        active_decomposition_generator_source_status_counts,
        active_decomposition_unresolved_source_leaf_sample,
        active_decomposition_source_leaf_unit_phase_probe_status_counts,
        active_decomposition_source_leaf_origin_omitted_unit_phase_probe_status_counts,
        active_decomposition_source_leaf_cms_solution_status_counts,
        active_decomposition_source_leaf_cms_candidate_status_counts,
        local_cygv_charge_signature_counts: missing_local_cygv_charge_signature_counts,
        local_cygv_target_candidate_status_counts,
        local_cygv_actual_call_readiness_counts: missing_local_cygv_actual_call_readiness_counts,
        local_cygv_missing_source_input_counts: missing_local_cygv_missing_source_input_counts,
        cms_general_divisor_candidate_status_counts:
            missing_cms_general_divisor_candidate_status_counts,
        cms_general_divisor_intersection_check_status_counts:
            missing_cms_general_divisor_intersection_check_status_counts,
        origin_circuit_witness_relation_status_counts:
            missing_origin_circuit_witness_relation_status_counts,
        uncovered_source_ray_origin_circuit_witness_relation_status_counts,
        shared_facet_unresolved_source_ray_target_count,
        shared_facet_unresolved_source_ray_origin_circuit_pattern_counts,
        shared_facet_unresolved_source_ray_origin_circuit_witness_relation_status_counts,
        shared_facet_unresolved_source_ray_local_cygv_charge_signature_counts,
        shared_facet_unresolved_source_ray_local_cygv_actual_call_readiness_counts,
        shared_facet_unresolved_source_ray_local_cygv_missing_source_input_counts,
        shared_facet_unresolved_source_ray_cms_general_divisor_candidate_status_counts,
        shared_facet_unresolved_source_ray_cms_general_divisor_intersection_check_status_counts,
        shared_facet_unresolved_source_ray_cms_solution_status_counts,
        shared_facet_unresolved_source_ray_cms_intersection_tensor_status_counts,
        shared_facet_unresolved_source_ray_cms_primitive_probe_status_counts,
        shared_facet_unresolved_source_ray_cms_primitive_probe_gv_counts,
        shared_facet_unresolved_source_ray_cms_solution_summary_error_counts,
        shared_facet_unresolved_source_ray_unit_phase_probe_status_counts,
        shared_facet_unresolved_source_ray_origin_omitted_unit_phase_probe_status_counts,
        shared_facet_unresolved_source_ray_unit_effective_tensor_requirement_status_counts,
        shared_facet_unresolved_source_ray_origin_omitted_effective_tensor_requirement_status_counts,
        shared_facet_unresolved_source_ray_unit_phase_probe_sample,
        shared_facet_unresolved_source_ray_sample,
        origin_circuit_facet_context_status_counts,
        origin_circuit_witness_relation_support_face_profile_counts,
        origin_circuit_witness_shared_facet_face_profile_counts,
        origin_circuit_witness_facet_union_face_profile_counts,
        origin_circuit_witness_relation_support_face_certificate_status_counts,
        origin_circuit_witness_shared_facet_face_certificate_status_counts,
        origin_circuit_witness_facet_union_face_certificate_status_counts,
        origin_circuit_witness_domain_sample,
        origin_circuit_witness_domain_unresolved_generator_unique_count,
        origin_circuit_witness_domain_unresolved_generator_occurrence_count,
        origin_circuit_witness_domain_unresolved_generator_status_counts,
        origin_circuit_witness_domain_unresolved_generator_degree_counts,
        origin_circuit_witness_domain_unresolved_generator_sample,
        origin_circuit_witness_shared_facet_unresolved_generator_unique_count,
        origin_circuit_witness_shared_facet_unresolved_generator_occurrence_count,
        origin_circuit_witness_shared_facet_unresolved_generator_status_counts,
        origin_circuit_witness_shared_facet_unresolved_generator_degree_counts,
        origin_circuit_witness_shared_facet_unresolved_generator_sample,
        active_support_status_counts,
        active_support_face_certificate_status_counts,
        target_extremal_ray_certificate_status_counts,
        origin_relation_support_face_certificate_status_counts,
        origin_shared_facet_face_certificate_status_counts,
        origin_facet_union_face_certificate_status_counts,
        cygv_path_history_status_counts,
        cygv_lower_seed_decomposition_status_counts,
        cygv_lower_seed_diamond_status_counts,
        path_support_uncovered_source_ray_unique_count,
        path_support_uncovered_source_ray_occurrence_count,
        path_support_uncovered_source_ray_degree_counts,
        path_support_uncovered_source_ray_lookup_status_counts,
        path_support_uncovered_source_ray_path_support_gv_counts,
        path_support_uncovered_source_ray_local_cygv_primitive_probe_status_counts,
        path_support_uncovered_source_ray_sample,
        local_cygv_q_matrix_orientation_status_counts,
        local_cygv_q_matrix_layout_status_counts,
        local_cygv_q_matrix_phase_status_counts,
        local_cygv_origin_point_status_counts,
        local_cygv_origin_omitted_compact_shape_status_counts,
        local_cygv_target_unit_phase_probe_status_counts,
        local_cygv_target_origin_omitted_unit_phase_probe_status_counts,
        local_cygv_target_unit_effective_tensor_requirement_status_counts,
        local_cygv_target_origin_omitted_effective_tensor_requirement_status_counts,
        local_cygv_target_unit_formula_sum_probe_status_counts,
        local_cygv_target_origin_omitted_unit_formula_sum_probe_status_counts,
        local_cygv_target_unit_formula_sum_effective_tensor_requirement_status_counts,
        local_cygv_target_origin_omitted_unit_formula_sum_effective_tensor_requirement_status_counts,
        local_cygv_target_unit_phase_probe_sample,
        local_cygv_target_integer_tensor_scan_status_counts,
        local_cygv_target_integer_tensor_scan_sample,
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

fn cms_general_divisor_candidate_status_counts(
    targets: &[TargetReport],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for target in targets {
        let status = match target.cms_general_divisor_shape_candidates.as_deref() {
            None => "cms_general_divisor_candidates_not_available",
            Some([]) => "cms_general_divisor_no_shape_candidates",
            Some(candidates)
                if candidates
                    .iter()
                    .any(|candidate| candidate.toric_gv1_formula_value.is_some()) =>
            {
                "cms_general_divisor_has_formula_candidate"
            }
            Some(_) => "cms_general_divisor_no_formula_candidate",
        };
        *counts.entry(status.to_string()).or_insert(0usize) += 1;
    }
    counts
}

fn cms_general_divisor_intersection_check_status_counts(
    targets: &[TargetReport],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for target in targets {
        let Some(checks) = target.cms_general_divisor_intersection_checks.as_deref() else {
            *counts
                .entry("cms_general_divisor_checks_not_available".to_string())
                .or_insert(0usize) += 1;
            continue;
        };
        if checks.is_empty() {
            *counts
                .entry("cms_general_divisor_no_intersection_checks".to_string())
                .or_insert(0usize) += 1;
            continue;
        }
        for check in checks {
            let status = cms_general_divisor_intersection_check_status(check);
            *counts.entry(status.to_string()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn cms_general_divisor_intersection_check_status(
    check: &CmsGeneralDivisorIntersectionCheck,
) -> &'static str {
    if !check.has_rational_divisor_solution {
        return "cms_general_divisor_no_rational_divisor_solution";
    }
    if check.solution_is_integral != Some(true) {
        return "cms_general_divisor_nonintegral_divisor_solution";
    }
    if check.matches_inferred_other_normal_degree != Some(true) {
        return "cms_general_divisor_integral_solution_mismatches_inferred_degree";
    }
    "cms_general_divisor_integral_solution_matches_inferred_degree"
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

fn origin_circuit_witness_relation_status_counts(
    samples: &[MissingGvTargetSample],
    target_index_filter: Option<usize>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for (idx, sample) in samples.iter().enumerate() {
        if target_index_filter.is_some_and(|filter| filter != idx) {
            continue;
        }
        *counts
            .entry(origin_circuit_witness_relation_status(sample))
            .or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_pattern_counts(samples: &[MissingGvTargetSample]) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for sample in samples {
        let key = sample
            .origin_circuit_pattern
            .clone()
            .unwrap_or_else(|| "not_origin_circuit".to_string());
        *counts.entry(key).or_insert(0usize) += 1;
    }
    counts
}

fn origin_circuit_witness_relation_status(sample: &MissingGvTargetSample) -> String {
    let Some(witnesses) = sample
        .origin_circuit_witnesses
        .as_ref()
        .filter(|witnesses| !witnesses.is_empty())
    else {
        return if sample.origin_circuit_first_witness.is_some() {
            if sample.origin_circuit_witness_count.unwrap_or(1) > 1 {
                "multiple_origin_circuit_witnesses_not_serialized".to_string()
            } else {
                "single_origin_circuit_witness".to_string()
            }
        } else {
            "no_origin_circuit_witness".to_string()
        };
    };

    if witnesses.len() == 1 {
        return "single_origin_circuit_witness".to_string();
    }

    let mut relation_signatures = BTreeSet::new();
    for witness in witnesses {
        relation_signatures.insert(origin_circuit_witness_relation_signature(witness));
    }
    if relation_signatures.len() == 1 {
        "all_origin_circuit_witnesses_share_relation".to_string()
    } else {
        "mixed_origin_circuit_witness_relations".to_string()
    }
}

fn origin_circuit_witness_relation_signature(
    witness: &OriginCircuitWitnessSample,
) -> Vec<(usize, i64)> {
    let mut signature = if witness.relation_points.is_empty() {
        witness.sparse_relation.clone()
    } else {
        witness
            .relation_points
            .iter()
            .map(|point| (point.point_index, point.coefficient))
            .collect()
    };
    signature.sort_unstable();
    signature
}

fn origin_circuit_facet_context_status(sample: &MissingGvTargetSample) -> String {
    let witnesses = origin_circuit_witnesses(sample);
    if witnesses.is_empty() {
        return "no_origin_circuit_witness".to_string();
    }

    let mut statuses = BTreeSet::new();
    for witness in witnesses {
        statuses.insert(origin_circuit_witness_facet_context_status(witness));
    }
    if statuses.len() == 1 {
        return statuses
            .into_iter()
            .next()
            .expect("single status exists")
            .to_string();
    }
    format!(
        "mixed_origin_circuit_facet_context:{}",
        statuses.into_iter().collect::<Vec<_>>().join("|")
    )
}

fn origin_circuit_witnesses(sample: &MissingGvTargetSample) -> Vec<&OriginCircuitWitnessSample> {
    if let Some(witnesses) = sample
        .origin_circuit_witnesses
        .as_ref()
        .filter(|witnesses| !witnesses.is_empty())
    {
        return witnesses.iter().collect();
    }
    sample.origin_circuit_first_witness.iter().collect()
}

fn origin_circuit_witness_facet_context_status(
    witness: &OriginCircuitWitnessSample,
) -> &'static str {
    if witness.first_facet.is_empty() || witness.second_facet.is_empty() {
        return "origin_circuit_missing_full_facet_context";
    }
    if witness.first_facet.len() != witness.first_facet_size
        || witness.second_facet.len() != witness.second_facet_size
    {
        return "origin_circuit_facet_context_size_mismatch";
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
        "source_derived_full_facet_context"
    } else {
        "origin_circuit_facet_context_inconsistent"
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

fn local_cygv_q_matrix_phase_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(skeleton.local_cygv_q_matrix_phase_status.clone())
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

fn local_cygv_origin_omitted_compact_shape_status_counts<'a>(
    skeletons: impl IntoIterator<Item = &'a LocalCygvInputSkeleton>,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for skeleton in skeletons {
        *counts
            .entry(local_cygv_origin_omitted_compact_shape_status(skeleton))
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_origin_omitted_compact_shape_status(skeleton: &LocalCygvInputSkeleton) -> String {
    let q_matrix = match local_origin_omitted_wrapper_q_matrix(skeleton) {
        Ok(Some(q_matrix)) => q_matrix,
        Ok(None) => return "origin_omitted_q_matrix_not_available".to_string(),
        Err(error) => return format!("origin_omitted_q_matrix_error:{error}"),
    };
    let Some(first_row) = q_matrix.first() else {
        return "origin_omitted_q_matrix_not_available".to_string();
    };
    if first_row.is_empty() {
        return "origin_omitted_q_matrix_empty".to_string();
    }
    if q_matrix.iter().any(|row| row.len() != first_row.len()) {
        return "origin_omitted_q_matrix_inconsistent_row_widths".to_string();
    }
    let h11 = match i64::try_from(q_matrix.len()) {
        Ok(value) => value,
        Err(_) => return "origin_omitted_q_matrix_rank_overflow".to_string(),
    };
    let h11pd = match i64::try_from(first_row.len()) {
        Ok(value) => value,
        Err(_) => return "origin_omitted_q_matrix_divisor_count_overflow".to_string(),
    };
    let cy_dim = h11pd - h11 - 1;
    if cy_dim == 3 {
        "origin_omitted_compact_threefold_hypersurface_shape".to_string()
    } else {
        format!("origin_omitted_cy_dim_{cy_dim}_not_compact_threefold")
    }
}

fn local_cygv_target_unit_phase_probe_summaries(
    targets: &[TargetReport],
) -> Vec<LocalCygvTargetUnitPhaseProbeSummary> {
    let mut summaries = Vec::new();
    for target in targets {
        let Some(skeleton) = target.local_cygv_input_skeleton.as_ref() else {
            continue;
        };
        let expected_values = target_expected_toric_gv1_formula_values(target);
        let expected_sum = target_expected_toric_gv1_formula_value_sum(target);
        let probe =
            local_cygv_unit_phase_probe_from_skeleton(skeleton, expected_values, expected_sum);
        summaries.push(LocalCygvTargetUnitPhaseProbeSummary {
            target_index: target.index,
            degree: target.degree,
            probe,
        });
    }
    summaries
}

fn target_expected_toric_gv1_formula_values(target: &TargetReport) -> Vec<String> {
    expected_toric_gv1_formula_values(target.cms_general_divisor_shape_candidates.as_deref())
}

fn target_expected_toric_gv1_formula_value_sum(target: &TargetReport) -> Option<String> {
    expected_toric_gv1_formula_value_sum(target.cms_general_divisor_shape_candidates.as_deref())
}

fn sample_expected_toric_gv1_formula_values(sample: &MissingGvTargetSample) -> Vec<String> {
    expected_toric_gv1_formula_values(sample.cms_general_divisor_shape_candidates.as_deref())
}

fn sample_expected_toric_gv1_formula_value_sum(sample: &MissingGvTargetSample) -> Option<String> {
    expected_toric_gv1_formula_value_sum(sample.cms_general_divisor_shape_candidates.as_deref())
}

fn expected_toric_gv1_formula_values(
    candidates: Option<&[CmsGeneralDivisorShapeCandidate]>,
) -> Vec<String> {
    let mut values = candidates
        .unwrap_or_default()
        .iter()
        .filter_map(|candidate| candidate.toric_gv1_formula_value)
        .collect::<Vec<_>>();
    values.sort_unstable();
    values.dedup();
    values.into_iter().map(|value| value.to_string()).collect()
}

fn expected_toric_gv1_formula_value_sum(
    candidates: Option<&[CmsGeneralDivisorShapeCandidate]>,
) -> Option<String> {
    let mut found = false;
    let mut sum = 0i64;
    for value in candidates
        .unwrap_or_default()
        .iter()
        .filter_map(|candidate| candidate.toric_gv1_formula_value)
    {
        found = true;
        sum = sum.checked_add(value)?;
    }
    found.then(|| sum.to_string())
}

fn local_cygv_target_unit_phase_probe_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        *counts
            .entry(probe.probe.unit_tensor_probe_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_target_origin_omitted_unit_phase_probe_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        *counts
            .entry(probe.probe.origin_omitted_unit_tensor_probe_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_target_unit_effective_tensor_requirement_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        for requirement in &probe.probe.unit_tensor_effective_tensor_requirements {
            *counts.entry(requirement.status.clone()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn local_cygv_target_origin_omitted_effective_tensor_requirement_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        for requirement in &probe
            .probe
            .origin_omitted_unit_tensor_effective_tensor_requirements
        {
            *counts.entry(requirement.status.clone()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn local_cygv_target_unit_formula_sum_probe_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        *counts
            .entry(probe.probe.unit_tensor_formula_sum_probe_status.clone())
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_target_origin_omitted_unit_formula_sum_probe_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        *counts
            .entry(
                probe
                    .probe
                    .origin_omitted_unit_tensor_formula_sum_probe_status
                    .clone(),
            )
            .or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_target_unit_formula_sum_effective_tensor_requirement_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        if let Some(requirement) = &probe
            .probe
            .unit_tensor_formula_sum_effective_tensor_requirement
        {
            *counts.entry(requirement.status.clone()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn local_cygv_target_origin_omitted_unit_formula_sum_effective_tensor_requirement_status_counts(
    probes: &[LocalCygvTargetUnitPhaseProbeSummary],
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for probe in probes {
        if let Some(requirement) = &probe
            .probe
            .origin_omitted_unit_tensor_formula_sum_effective_tensor_requirement
        {
            *counts.entry(requirement.status.clone()).or_insert(0usize) += 1;
        }
    }
    counts
}

fn local_cygv_integer_tensor_scan_status_counts(
    scans: &[LocalCygvIntegerTensorScanSummary],
    target_count: usize,
    scan_enabled: bool,
) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    if !scan_enabled {
        counts.insert("not_run".to_string(), target_count);
        return counts;
    }
    for scan in scans {
        *counts.entry(scan.status.clone()).or_insert(0usize) += 1;
    }
    counts
}

fn local_cygv_integer_tensor_scan_summaries(
    targets: &[TargetReport],
    scan_enabled: bool,
    scan_bound: i64,
) -> Vec<LocalCygvIntegerTensorScanSummary> {
    if !scan_enabled {
        return Vec::new();
    }
    targets
        .iter()
        .map(|target| local_cygv_integer_tensor_scan_summary(target, scan_bound))
        .collect()
}

fn local_cygv_integer_tensor_scan_summary(
    target: &TargetReport,
    scan_bound: i64,
) -> LocalCygvIntegerTensorScanSummary {
    let expected_values = target_expected_toric_gv1_formula_values(target);
    let blocked = |status: &str| LocalCygvIntegerTensorScanSummary {
        target_index: target.index,
        degree: target.degree,
        status: status.to_string(),
        scan_bound,
        expected_toric_gv1_formula_values: expected_values.clone(),
        entries: Vec::new(),
    };
    if scan_bound <= 0 {
        return blocked("integer_tensor_scan_invalid_bound");
    }
    if expected_values.is_empty() {
        return blocked("integer_tensor_scan_no_expected_formula_values");
    }
    let Some(skeleton) = target.local_cygv_input_skeleton.as_ref() else {
        return blocked("integer_tensor_scan_missing_skeleton");
    };
    let Some(q_matrix) = skeleton.local_cygv_phase_q_matrix_candidate.as_ref() else {
        return blocked("integer_tensor_scan_missing_q_matrix");
    };
    if q_matrix.len() != 1 {
        return blocked("integer_tensor_scan_not_one_parameter_q_matrix");
    }
    let Some(grading_vector) = skeleton.local_grading_vector_candidate.as_ref() else {
        return blocked("integer_tensor_scan_missing_grading_vector");
    };
    let Some(generators) = skeleton.local_semigroup_generators_candidate.as_ref() else {
        return blocked("integer_tensor_scan_missing_semigroup_generators");
    };
    if generators.as_slice() != [vec![1]] {
        return blocked("integer_tensor_scan_not_unit_primitive_semigroup");
    }

    let semigroup_elements = vec![vec![0], vec![1]];
    let mut entries = Vec::new();
    for tensor_value in -scan_bound..=scan_bound {
        if tensor_value == 0 {
            continue;
        }
        match one_parameter_primitive_cygv_value(
            q_matrix,
            grading_vector,
            &semigroup_elements,
            MalachiteRational::from(tensor_value),
        ) {
            Ok(candidate_gv) => {
                let status = if expected_values
                    .iter()
                    .any(|expected| expected == &candidate_gv)
                {
                    "integer_tensor_scan_matches_expected_formula_but_uncertified"
                } else {
                    "integer_tensor_scan_mismatch"
                };
                entries.push(LocalCygvIntegerTensorScanEntry {
                    tensor_value,
                    candidate_gv: Some(candidate_gv),
                    status: status.to_string(),
                    error: None,
                });
            }
            Err(error) => entries.push(LocalCygvIntegerTensorScanEntry {
                tensor_value,
                candidate_gv: None,
                status: "integer_tensor_scan_hkty_error".to_string(),
                error: Some(error),
            }),
        }
    }
    let status = if entries
        .iter()
        .any(|entry| entry.status == "integer_tensor_scan_matches_expected_formula_but_uncertified")
    {
        "integer_tensor_scan_has_expected_match_but_uncertified"
    } else if entries
        .iter()
        .all(|entry| entry.status == "integer_tensor_scan_hkty_error")
    {
        "integer_tensor_scan_all_hkty_error"
    } else {
        "integer_tensor_scan_no_expected_match"
    };
    LocalCygvIntegerTensorScanSummary {
        target_index: target.index,
        degree: target.degree,
        status: status.to_string(),
        scan_bound,
        expected_toric_gv1_formula_values: expected_values,
        entries,
    }
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
            "[ERROR] usage: mcallister_gv_context --context path [--target-index N] [--run-integer-diamonds] [--run-active-support-generators] [--run-support-overlap-generators N] [--pair-reduce-support-overlap-generators] [--support-overlap-max-target-degree N] [--certify-origin-support-domains] [--certify-origin-witness-domains] [--origin-support-certificate-limit N] [--certify-target-extremal-rays] [--target-extremal-generator-limit N] [--target-extremal-max-degree N] [--measure-cygv-semigroups] [--probe-cygv-path-history] [--run-lower-seed-diamonds] [--run-path-support-generators] [--measure-cygv-degree-ladder --cygv-degree-ladder-max-degree N] [--semigroup-measure-max-target-degree N] [--semigroup-measure-max-seeds N] [--scan-local-integer-tensors] [--local-tensor-scan-bound N] [--element-limit N] [--closure-generation-limit N] [--out path]\n       use --run-support-overlap-generators 0 to try all degree-bounded generators up to each target degree"
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
    let certify_origin_witness_domains = parse_flag("--certify-origin-witness-domains");
    let origin_support_certificate_limit =
        parse_arg_value::<usize>("--origin-support-certificate-limit").unwrap_or(256);
    let certify_target_extremal_rays = parse_flag("--certify-target-extremal-rays");
    let target_extremal_generator_limit =
        parse_arg_value::<usize>("--target-extremal-generator-limit").unwrap_or(256);
    let target_extremal_max_degree = parse_arg_value::<i128>("--target-extremal-max-degree");
    let measure_cygv_semigroups = parse_flag("--measure-cygv-semigroups");
    let probe_cygv_path_history = parse_flag("--probe-cygv-path-history");
    let run_lower_seed_diamonds = parse_flag("--run-lower-seed-diamonds");
    let run_path_support_generators = parse_flag("--run-path-support-generators");
    let measure_cygv_degree_ladder = parse_flag("--measure-cygv-degree-ladder");
    let cygv_degree_ladder_max_degree = parse_arg_value::<i128>("--cygv-degree-ladder-max-degree");
    let semigroup_measure_max_target_degree =
        parse_arg_value::<i128>("--semigroup-measure-max-target-degree");
    let semigroup_measure_max_seed_count =
        parse_arg_value::<usize>("--semigroup-measure-max-seeds");
    let scan_local_integer_tensors = parse_flag("--scan-local-integer-tensors");
    let local_tensor_scan_bound = parse_arg_value::<i64>("--local-tensor-scan-bound").unwrap_or(8);
    let element_limit = parse_arg_value::<usize>("--element-limit").unwrap_or(256);
    let closure_generation_limit = parse_arg_value::<usize>("--closure-generation-limit");
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
        certify_origin_witness_domains,
        origin_support_certificate_limit,
        certify_target_extremal_rays,
        target_extremal_generator_limit,
        target_extremal_max_degree,
        measure_cygv_semigroups,
        probe_cygv_path_history,
        run_lower_seed_diamonds,
        run_path_support_generators,
        measure_cygv_degree_ladder,
        cygv_degree_ladder_max_degree,
        semigroup_measure_max_target_degree,
        semigroup_measure_max_seed_count,
        scan_local_integer_tensors,
        local_tensor_scan_bound,
        element_limit,
        closure_generation_limit,
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
            covered_toric_gv_context_for_missing: None,
            uncovered_source_ray_toric_diagnostic_sample: None,
            degree_bounded_toric_gv_diagnostic_context_for_missing: None,
            uncovered_source_ray_stats_for_missing: None,
            shared_facet_unresolved_source_ray_stats_for_missing: None,
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

    fn minimal_missing_sample(basis_nonzero: Vec<(usize, i64)>) -> MissingGvTargetSample {
        MissingGvTargetSample {
            degree: 1,
            generators_le_degree: 0,
            is_mori_generator: true,
            origin_circuit_pattern: None,
            origin_circuit_witness_count: None,
            origin_circuit_first_witness: None,
            origin_circuit_witnesses: None,
            origin_circuit_affine_support: None,
            cms_general_divisor_shape_candidates: None,
            cms_general_divisor_intersection_checks: None,
            branch_diagnostic: None,
            real_cone_decomposable_by_other_generators: false,
            real_cone_decomposition_active_generators: None,
            real_cone_decomposition_active_generator_basis_nonzero: None,
            real_cone_decomposition_exact_coefficients: None,
            real_cone_decomposition_exact_kind: None,
            ambient_nonzero: basis_nonzero.clone(),
            basis_nonzero,
        }
    }

    fn skeleton_with_origin_phase_probe(
        support_point_indices: Vec<usize>,
        wrapper_q_matrix: Option<Vec<Vec<i64>>>,
    ) -> LocalCygvInputSkeleton {
        let (phase_q_matrix, phase_status) = local_cygv_q_matrix_phase_candidate(
            &support_point_indices,
            wrapper_q_matrix.as_deref(),
        );
        LocalCygvInputSkeleton {
            support_contains_origin_point: support_point_indices.contains(&0),
            support_point_indices,
            local_cygv_origin_point_status: String::new(),
            origin_point_relation_coefficient: None,
            local_cytools_origin_circuit_status: String::new(),
            local_q_matrix_rows: Vec::new(),
            target_relation_coefficients: None,
            target_relation_in_charge_basis: None,
            target_relation_status: String::new(),
            local_semigroup_generators_candidate: Some(vec![vec![1]]),
            local_semigroup_generator_status: "source_derived_one_parameter_unit_semigroup"
                .to_string(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: wrapper_q_matrix,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: phase_q_matrix,
            local_cygv_q_matrix_phase_status: phase_status,
            local_grading_vector_candidate: Some(vec![1]),
            local_grading_vector_status: "source_derived_primitive_one_parameter_grading"
                .to_string(),
            remaining_uncertified_inputs: Vec::new(),
        }
    }

    #[test]
    fn validate_context_accepts_schema3_degree_bounded_ray_context() {
        let mut context = minimal_corrected_context(
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
        context.covered_toric_gv_context_for_missing = Some(vec![CoveredToricGvContextSample {
            degree: 1,
            gv: "42".to_string(),
            ambient_nonzero: vec![(5, 1)],
            basis_nonzero: vec![(0, 1)],
        }]);
        context.uncovered_source_ray_toric_diagnostic_sample =
            Some(vec![ToricGvDiagnosticContextSample {
                degree: 1,
                gv: "-2".to_string(),
                source_bucket: "two_face".to_string(),
                source_summary: "gv=-2;source_count=1".to_string(),
                ambient_nonzero: vec![(8, 1)],
                basis_nonzero: vec![(1, 1)],
            }]);
        context.degree_bounded_toric_gv_diagnostic_context_for_missing =
            Some(vec![ToricGvDiagnosticContextSample {
                degree: 2,
                gv: "1".to_string(),
                source_bucket: "degree_bounded_two_face".to_string(),
                source_summary: "gv=1;source_count=1".to_string(),
                ambient_nonzero: vec![(5, 1), (8, 1)],
                basis_nonzero: vec![(0, 1), (1, 1)],
            }]);
        let validated = validate_context(&context).unwrap();

        assert_eq!(validated.degree_bounded_ray_context.unwrap().len(), 2);
        assert_eq!(
            validated
                .covered_toric_gv_by_basis
                .get(&vec![1, 0])
                .map(String::as_str),
            Some("42")
        );
        assert_eq!(
            validated
                .source_derived_gv_by_basis
                .get(&vec![0, 1])
                .map(String::as_str),
            Some("-2")
        );
        assert_eq!(
            validated
                .source_derived_gv_by_basis
                .get(&vec![1, 1])
                .map(String::as_str),
            Some("1")
        );
        assert_eq!(
            path_support_source_class_context(&[0, 1], &validated)
                .unwrap()
                .status,
            "source_ray_known_source_derived_gv"
        );
        assert_eq!(
            degree_bounded_mori_ray_context_status(&validated),
            "source_derived_ambient_and_basis_degree_bounded_mori_ray_context"
        );
    }

    #[test]
    fn pair_expansion_replaces_nonreduced_seed_terms() {
        let seed_set = [vec![1, 0], vec![0, 1], vec![1, 1]]
            .into_iter()
            .collect::<HashSet<_>>();
        let reduced_seed_set = [vec![1, 0], vec![0, 1]].into_iter().collect::<HashSet<_>>();

        let expanded = pair_expand_terms_to_reduced_seeds(
            &[vec![1, 1], vec![1, 0]],
            &seed_set,
            &reduced_seed_set,
            4,
            4,
        )
        .unwrap()
        .unwrap();

        assert_eq!(expanded, vec![vec![0, 1], vec![1, 0], vec![1, 0]]);
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
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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
            origin_circuit_witnesses: None,
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

        let summaries = origin_circuit_witness_domain_summaries(
            std::slice::from_ref(&sample),
            &context,
            None,
            false,
            2,
        );
        assert_eq!(summaries.len(), 1);
        assert_eq!(
            summaries[0].relation_support_face_profile,
            "support_face_codim_1_exact_kernel_candidate_generators_1"
        );
        assert_eq!(
            summaries[0].shared_facet_face_profile,
            "support_face_full_dimensional_generators_2_4"
        );
        assert_eq!(
            summaries[0].facet_union_face_profile,
            "support_face_full_dimensional_generators_2_4"
        );

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
    fn origin_support_certificate_skips_single_generator_higher_codimension_lp() {
        let stats = MissingGvTargetStats {
            target_count: 0,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: Vec::new(),
        };
        let grading = vec![1, 1, 1];
        let q_matrix = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let degree_bounded_rays = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let ray_context = vec![DegreeBoundedMoriRayContextSample {
            degree: 1,
            ambient_nonzero: vec![(4, 1)],
            basis_nonzero: vec![(0, 1)],
        }];
        let context = ValidatedContext {
            dimension: 3,
            degree_bound: 1,
            q_cols: 3,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: Some(&ray_context),
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(3),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
        };
        let allowed_support = [4usize].into_iter().collect::<HashSet<_>>();

        let status = origin_circuit_ambient_support_face_certificate_status(
            &ray_context,
            1,
            &allowed_support,
            &context,
            16,
        );

        assert_eq!(
            status,
            "origin_support_skipped_single_generator_higher_codimension_rank_1_dim_3"
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
            origin_circuit_witnesses: None,
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
        assert_eq!(
            skeleton.local_semigroup_generators_candidate,
            Some(vec![vec![1]])
        );
        assert_eq!(
            skeleton.local_semigroup_generator_status,
            "source_derived_one_parameter_unit_semigroup"
        );
        assert!(
            !skeleton
                .remaining_uncertified_inputs
                .iter()
                .any(|input| input == "local_semigroup_generators")
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
            !skeleton
                .remaining_uncertified_inputs
                .iter()
                .any(|input| input == "local_q_matrix_phase")
        );
        assert_eq!(
            skeleton.local_cygv_phase_q_matrix_candidate,
            Some(vec![vec![-1, 2, 1, -3, 1]])
        );
        assert_eq!(
            skeleton.local_cygv_q_matrix_phase_status,
            "source_derived_unique_compact_threefold_phase_including_origin"
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: local_cygv_orientation_candidates(
                &[vec![1], vec![-2], vec![-1], vec![3], vec![-1]],
                Some(&[-1]),
            ),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: local_cygv_orientation_candidates(
                &[vec![2], vec![1], vec![2], vec![-1], vec![-2], vec![-2]],
                Some(&[-1]),
            ),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: supported_candidate.clone(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: supported_candidate,
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: local_cygv_orientation_candidates(
                &[vec![2], vec![1], vec![2], vec![-1], vec![-2], vec![-2]],
                Some(&[1]),
            ),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status: String::new(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: String::new(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: Some(-1),
            local_q_matrix_orientation_status: "source_derived_target_positive_orientation"
                .to_string(),
            local_cygv_q_matrix_rows_candidate: Some(vec![vec![-1], vec![2]]),
            local_cygv_wrapper_q_matrix_candidate: Some(vec![vec![-1, 2]]),
            local_cygv_q_matrix_layout_status: "source_derived_oriented_q_matrix_layout"
                .to_string(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: None,
            local_q_matrix_orientation_status:
                "local_q_orientation_blocked_no_positive_primitive_target".to_string(),
            local_cygv_q_matrix_rows_candidate: None,
            local_cygv_wrapper_q_matrix_candidate: None,
            local_cygv_q_matrix_layout_status: "local_q_matrix_layout_blocked_no_orientation"
                .to_string(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
            local_semigroup_generators_candidate: None,
            local_semigroup_generator_status: String::new(),
            orientation_candidates: Vec::new(),
            local_q_matrix_orientation_candidate: Some(-1),
            local_q_matrix_orientation_status: "source_derived_target_positive_orientation"
                .to_string(),
            local_cygv_q_matrix_rows_candidate: Some(vec![vec![-1], vec![1]]),
            local_cygv_wrapper_q_matrix_candidate: Some(vec![vec![-1, 1]]),
            local_cygv_q_matrix_layout_status: "source_derived_oriented_q_matrix_layout"
                .to_string(),
            local_cygv_phase_q_matrix_candidate: None,
            local_cygv_q_matrix_phase_status: String::new(),
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
    fn local_cygv_origin_omitted_compact_shape_counts_phase_candidates() {
        let compact_threefold = skeleton_with_origin_phase_probe(
            vec![0, 1, 2, 3, 4, 5],
            Some(vec![vec![-1, 1, -1, 1, -1, 1]]),
        );
        let too_low_dimensional = skeleton_with_origin_phase_probe(
            vec![0, 1, 2, 3, 4],
            Some(vec![vec![-1, -2, 1, 1, 1]]),
        );
        let no_origin = skeleton_with_origin_phase_probe(vec![1, 2, 3], Some(vec![vec![1, 1, 1]]));

        let counts = local_cygv_origin_omitted_compact_shape_status_counts([
            &compact_threefold,
            &too_low_dimensional,
            &no_origin,
        ]);

        assert_eq!(
            counts
                .get("origin_omitted_compact_threefold_hypersurface_shape")
                .copied(),
            Some(1)
        );
        assert_eq!(
            counts
                .get("origin_omitted_cy_dim_2_not_compact_threefold")
                .copied(),
            Some(1)
        );
        assert_eq!(
            counts.get("origin_omitted_q_matrix_not_available").copied(),
            Some(1)
        );
    }

    #[test]
    fn local_cygv_unit_phase_probe_compares_expected_formula_set() {
        let resolved_like = skeleton_with_origin_phase_probe(
            vec![0, 1, 2, 3, 4, 5],
            Some(vec![vec![-1, 1, -1, 1, -1, 1]]),
        );

        let probe = local_cygv_unit_phase_probe_from_skeleton(
            &resolved_like,
            vec!["1".into()],
            Some("1".into()),
        );

        assert_eq!(probe.unit_tensor_candidate_gv.as_deref(), Some("1"));
        assert_eq!(
            probe.unit_tensor_probe_status,
            "unit_tensor_probe_matches_expected_formula_set_but_uncertified"
        );
        assert_eq!(
            probe.origin_omitted_q_matrix,
            Some(vec![vec![1, -1, 1, -1, 1]])
        );
        assert_eq!(
            probe.origin_omitted_unit_tensor_candidate_gv.as_deref(),
            Some("1")
        );
        assert_eq!(
            probe.origin_omitted_unit_tensor_probe_status,
            "origin_omitted_unit_tensor_probe_matches_expected_formula_set_but_phase_uncertified"
        );
        assert_eq!(
            probe.unit_tensor_effective_tensor_requirements[0].status,
            "effective_tensor_integral_candidate_but_uncertified"
        );
        assert_eq!(
            probe.origin_omitted_unit_tensor_effective_tensor_requirements[0]
                .required_tensor_value
                .as_deref(),
            Some("1")
        );
        assert_eq!(
            probe.origin_omitted_unit_tensor_effective_tensor_requirements[0].status,
            "effective_tensor_integral_candidate_but_uncertified"
        );
        assert_eq!(
            probe.unit_tensor_formula_sum_probe_status,
            "unit_tensor_matches_formula_sum_but_uncertified"
        );
        assert_eq!(
            probe
                .unit_tensor_formula_sum_effective_tensor_requirement
                .as_ref()
                .and_then(|requirement| requirement.required_tensor_value.as_deref()),
            Some("1")
        );
    }

    #[test]
    fn effective_tensor_requirement_flags_nonintegral_local_reductions() {
        let requirement = effective_tensor_requirement_for_unit_probe(Some("6"), "3");

        assert_eq!(requirement.required_tensor_value.as_deref(), Some("1/2"));
        assert_eq!(
            requirement.status,
            "effective_tensor_nonintegral_candidate_rejected_by_cygv_threefold_intnums"
        );
    }

    #[test]
    fn expected_formula_sum_preserves_candidate_multiplicity() {
        let candidates = vec![
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 1,
                shrinking_divisor_coefficient: -1,
                shrinking_divisor_coordinates: vec![0, 1],
                inferred_other_normal_degree: -1,
                toric_gv1_formula_value: Some(1),
                all_non_origin_relation_points_are_two_face: true,
            },
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 2,
                shrinking_divisor_coefficient: -2,
                shrinking_divisor_coordinates: vec![1, 0],
                inferred_other_normal_degree: 0,
                toric_gv1_formula_value: Some(-2),
                all_non_origin_relation_points_are_two_face: true,
            },
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 3,
                shrinking_divisor_coefficient: -1,
                shrinking_divisor_coordinates: vec![1, 1],
                inferred_other_normal_degree: -1,
                toric_gv1_formula_value: Some(1),
                all_non_origin_relation_points_are_two_face: true,
            },
        ];

        assert_eq!(
            expected_toric_gv1_formula_values(Some(&candidates)),
            vec!["-2".to_string(), "1".to_string()]
        );
        assert_eq!(
            expected_toric_gv1_formula_value_sum(Some(&candidates)).as_deref(),
            Some("0")
        );
    }

    #[test]
    fn local_cygv_integer_tensor_scan_status_counts_are_opt_in() {
        let scans = vec![
            LocalCygvIntegerTensorScanSummary {
                target_index: 0,
                degree: 10,
                status: "integer_tensor_scan_has_expected_match_but_uncertified".to_string(),
                scan_bound: 8,
                expected_toric_gv1_formula_values: vec!["3".to_string()],
                entries: Vec::new(),
            },
            LocalCygvIntegerTensorScanSummary {
                target_index: 1,
                degree: 18,
                status: "integer_tensor_scan_no_expected_match".to_string(),
                scan_bound: 8,
                expected_toric_gv1_formula_values: vec!["-2".to_string(), "1".to_string()],
                entries: Vec::new(),
            },
        ];

        let disabled = local_cygv_integer_tensor_scan_status_counts(&scans, 9, false);
        assert_eq!(disabled, BTreeMap::from([("not_run".to_string(), 9)]));

        let enabled = local_cygv_integer_tensor_scan_status_counts(&scans, 9, true);
        assert_eq!(
            enabled,
            BTreeMap::from([
                (
                    "integer_tensor_scan_has_expected_match_but_uncertified".to_string(),
                    1
                ),
                ("integer_tensor_scan_no_expected_match".to_string(), 1),
            ])
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
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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
                origin_circuit_witnesses: None,
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
        let mut covered_toric_gv_by_basis = HashMap::new();
        covered_toric_gv_by_basis.insert(vec![1, 0], "7".to_string());
        covered_toric_gv_by_basis.insert(vec![0, 1], "11".to_string());
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 3,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            covered_toric_gv_by_basis,
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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
                origin_circuit_witnesses: None,
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
        let mut covered_toric_gv_by_basis = HashMap::new();
        covered_toric_gv_by_basis.insert(vec![1, 0], "7".to_string());
        covered_toric_gv_by_basis.insert(vec![0, 1], "11".to_string());
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 2,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            covered_toric_gv_by_basis,
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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
                origin_circuit_witnesses: None,
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
        let mut covered_toric_gv_by_basis = HashMap::new();
        covered_toric_gv_by_basis.insert(vec![1, 0], "7".to_string());
        covered_toric_gv_by_basis.insert(vec![0, 1], "11".to_string());
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 2,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            covered_toric_gv_by_basis,
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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

        let cheap_probe = target_extremal_ray_certificate_probe(
            &stats.sample[0],
            &[1, 1],
            &context,
            false,
            0,
            None,
        )
        .expect("exact non-extremal decomposition should not require separator search");
        assert_eq!(
            cheap_probe.status,
            "not_extremal_by_exact_integer_semigroup_decomposition"
        );
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
        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 2, 16, None).unwrap();
        assert!(closure.completed);
        assert_eq!(closure.elements.len(), 6);
        assert_eq!(closure.degree_counts.get(&0), Some(&1));
        assert_eq!(closure.degree_counts.get(&1), Some(&2));
        assert_eq!(closure.degree_counts.get(&2), Some(&3));
        assert!(closure.elements.contains(&vec![1, 1]));
        assert_eq!(closure.generation_counts.len(), 1);
        assert_eq!(closure.generation_counts[0].generation, 1);
        assert_eq!(closure.generation_counts[0].starting_element_count, 2);
        assert_eq!(closure.generation_counts[0].new_element_count, 3);
        assert_eq!(
            closure.generation_counts[0].total_element_count_after_full_generation,
            6
        );
        assert!(!closure.generation_counts[0].truncated_at_limit);
    }

    #[test]
    fn bounded_cygv_closure_matches_actual_cygv_after_pair_seed_reduction() {
        let seeds = vec![vec![1, 0], vec![0, 1], vec![1, 1], vec![2, 0]];
        let reduced = cygv_pair_reduced_seed_generators(&seeds).unwrap();
        assert_eq!(reduced, vec![vec![0, 1], vec![1, 0]]);

        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 2, 16, None).unwrap();
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
        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 2, 4, None).unwrap();

        assert_eq!(closure.status, "exceeded_element_limit_4");
        assert!(!closure.completed);
        assert_eq!(
            closure.elements,
            HashSet::from([vec![0, 0], vec![0, 1], vec![1, 0], vec![0, 2]])
        );
        assert_eq!(closure.generation_counts.len(), 1);
        assert_eq!(closure.generation_counts[0].new_element_count, 3);
        assert_eq!(
            closure.generation_counts[0].total_element_count_after_full_generation,
            6
        );
        assert!(closure.generation_counts[0].truncated_at_limit);
    }

    #[test]
    fn bounded_cygv_closure_can_stop_after_full_generation() {
        let seeds = vec![vec![1, 0], vec![0, 1]];
        let closure = bounded_cygv_semigroup_closure(&seeds, &[1, 1], 3, 16, Some(1)).unwrap();

        assert_eq!(closure.status, "stopped_generation_limit_1");
        assert!(!closure.completed);
        assert_eq!(closure.elements.len(), 6);
        assert_eq!(closure.degree_counts.get(&2), Some(&3));
        assert_eq!(closure.degree_counts.get(&3), None);
        assert_eq!(closure.generation_counts.len(), 1);
        assert!(!closure.generation_counts[0].truncated_at_limit);
    }

    #[test]
    fn path_candidate_support_merges_target_path_and_seed_sum_terms() {
        let candidate = CygvPathPredecessorCandidate {
            predecessor_degree: 3,
            difference_degree: 4,
            series_distance: "0.000000".to_string(),
            predecessor_toric_gv: None,
            difference_toric_gv: None,
            predecessor_source_derived_gv: None,
            difference_source_derived_gv: None,
            predecessor_known_qn_history_status: "unknown_not_toric_covered".to_string(),
            difference_known_qn_history_status: "unknown_not_toric_covered".to_string(),
            known_qn_history_pair_status:
                "predecessor_unknown_not_toric_covered__difference_unknown_not_toric_covered"
                    .to_string(),
            predecessor_is_seed: false,
            difference_is_seed: false,
            predecessor_is_reduced_seed: false,
            difference_is_reduced_seed: false,
            predecessor_first_generation_seed_sum: Some(CygvSeedSumDecomposition {
                reduced_seed_degree: 1,
                seed_degree: 2,
                reduced_seed_toric_gv: None,
                seed_toric_gv: None,
                reduced_seed_source_derived_gv: None,
                seed_source_derived_gv: None,
                reduced_seed_known_qn_history_status: "unknown_not_toric_covered".to_string(),
                seed_known_qn_history_status: "unknown_not_toric_covered".to_string(),
                seed_is_reduced_seed: false,
                seed_pair_reduction_sum: Some(CygvSeedPairDecomposition {
                    lhs_degree: 1,
                    rhs_degree: 1,
                    lhs_toric_gv: None,
                    rhs_toric_gv: None,
                    lhs_source_derived_gv: None,
                    rhs_source_derived_gv: None,
                    lhs_known_qn_history_status: "unknown_not_toric_covered".to_string(),
                    rhs_known_qn_history_status: "unknown_not_toric_covered".to_string(),
                    lhs_is_reduced_seed: true,
                    rhs_is_reduced_seed: true,
                    lhs_nonzero: vec![(4, 3), (5, 0)],
                    rhs_nonzero: vec![(6, 1)],
                }),
                reduced_seed_nonzero: vec![(2, 1)],
                seed_nonzero: vec![(3, -1)],
            }),
            difference_first_generation_seed_sum: None,
            predecessor_nonzero: vec![(1, 1)],
            difference_nonzero: vec![(3, 0), (4, 1)],
        };

        let support = path_candidate_support(&[0, 5, 0, 0, 0, 0, 0, 1], &[candidate]);

        assert_eq!(support, HashSet::from([1, 2, 3, 4, 6, 7]));
    }

    #[test]
    fn known_qn_history_status_distinguishes_zero_nonzero_and_unknown() {
        assert_eq!(
            known_qn_history_status(Some("0"), None).unwrap(),
            "known_zero_toric_gv"
        );
        assert_eq!(
            known_qn_history_status(Some("-3"), None).unwrap(),
            "known_nonzero_toric_gv"
        );
        assert_eq!(
            known_qn_history_status(None, Some("5")).unwrap(),
            "known_nonzero_source_gv"
        );
        assert_eq!(
            known_qn_history_status(None, Some("0")).unwrap(),
            "known_zero_source_gv"
        );
        assert_eq!(
            known_qn_history_status(Some("7"), Some("7")).unwrap(),
            "known_nonzero_toric_gv"
        );
        assert!(
            known_qn_history_status(Some("7"), Some("5"))
                .unwrap_err()
                .contains("conflicts")
        );
        assert_eq!(
            known_qn_history_status(None, None).unwrap(),
            "unknown_not_toric_covered"
        );
        assert_eq!(
            known_qn_history_pair_status("known_zero_toric_gv", "unknown_not_toric_covered"),
            "predecessor_known_zero_toric_gv__difference_unknown_not_toric_covered"
        );
    }

    #[test]
    fn source_derived_gv_requires_integral_matching_cms_certificate() {
        let mut sample = minimal_missing_sample(vec![(0, 2)]);
        sample.cms_general_divisor_shape_candidates = Some(vec![
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 5,
                shrinking_divisor_coefficient: -2,
                shrinking_divisor_coordinates: vec![1, 0, 0, 0],
                inferred_other_normal_degree: 0,
                toric_gv1_formula_value: Some(-2),
                all_non_origin_relation_points_are_two_face: true,
            },
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 6,
                shrinking_divisor_coefficient: -1,
                shrinking_divisor_coordinates: vec![0, 1, 0, 0],
                inferred_other_normal_degree: -1,
                toric_gv1_formula_value: Some(1),
                all_non_origin_relation_points_are_two_face: true,
            },
        ]);
        sample.cms_general_divisor_intersection_checks = Some(vec![
            CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 5,
                has_rational_divisor_solution: true,
                solution_basis_support_len: Some(1),
                solution_basis_nonzero: Some(vec![(0, "1".to_string())]),
                solution_ambient_basis_nonzero: Some(vec![(5, "1".to_string())]),
                solution_is_integral: Some(true),
                computed_other_normal_degree: Some("0".to_string()),
                matches_inferred_other_normal_degree: Some(true),
            },
            CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 6,
                has_rational_divisor_solution: false,
                solution_basis_support_len: None,
                solution_basis_nonzero: None,
                solution_ambient_basis_nonzero: None,
                solution_is_integral: None,
                computed_other_normal_degree: None,
                matches_inferred_other_normal_degree: None,
            },
        ]);

        assert_eq!(
            source_derived_gv_for_sample(&sample).unwrap().as_deref(),
            Some("-2")
        );

        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![sample],
        };
        let map = source_derived_gv_by_basis(Some(&stats), 2).unwrap();
        assert_eq!(map.get(&vec![2, 0]).map(String::as_str), Some("-2"));
    }

    #[test]
    fn seed_sum_decomposition_reports_source_derived_history_status() {
        let seed_set = [vec![1, 0], vec![0, 1]].into_iter().collect::<HashSet<_>>();
        let reduced_seed_set = seed_set.clone();
        let mut covered_toric_gv_by_basis = HashMap::new();
        covered_toric_gv_by_basis.insert(vec![1, 0], "7".to_string());
        let mut source_derived_gv_by_basis = HashMap::new();
        source_derived_gv_by_basis.insert(vec![0, 1], "5".to_string());

        let decomposition = first_generation_seed_sum_decomposition(
            &[1, 1],
            &[1, 1],
            &seed_set,
            &reduced_seed_set,
            &covered_toric_gv_by_basis,
            &source_derived_gv_by_basis,
        )
        .unwrap()
        .unwrap();

        assert_eq!(
            decomposition.reduced_seed_source_derived_gv.as_deref(),
            Some("5")
        );
        assert_eq!(decomposition.seed_toric_gv.as_deref(), Some("7"));
        assert_eq!(
            decomposition.reduced_seed_known_qn_history_status,
            "known_nonzero_source_gv"
        );
        assert_eq!(
            decomposition.seed_known_qn_history_status,
            "known_nonzero_toric_gv"
        );
    }

    #[test]
    fn closest_known_qn_residual_predecessor_traces_unknown_residual() {
        let elements = [vec![0, 1], vec![2, 0], vec![2, 1]]
            .into_iter()
            .collect::<HashSet<_>>();
        let degree_counts = BTreeMap::from([(1, 1), (2, 1), (3, 1)]);
        let covered_toric_gv_by_basis = HashMap::new();
        let mut source_derived_gv_by_basis = HashMap::new();
        source_derived_gv_by_basis.insert(vec![0, 1], "5".to_string());
        let closest = CygvClosestKnownQnPredecessor {
            predecessor_degree: 1,
            difference_degree: 3,
            series_distance: "3.000000".to_string(),
            predecessor_toric_gv: None,
            predecessor_source_derived_gv: None,
            predecessor_known_qn_history_status: "known_nonzero_toric_gv".to_string(),
            difference_toric_gv: None,
            difference_source_derived_gv: None,
            difference_known_qn_history_status: "unknown_not_toric_covered".to_string(),
            predecessor_nonzero: vec![(0, 1)],
            difference_nonzero: vec![(0, 2), (1, 1)],
        };

        let residual = cygv_closest_known_qn_residual_predecessor(
            &elements,
            &degree_counts,
            &[1, 1],
            2,
            2,
            Some(&closest),
            &covered_toric_gv_by_basis,
            &source_derived_gv_by_basis,
        )
        .unwrap()
        .unwrap();

        assert_eq!(
            residual.predecessor_known_qn_history_status,
            "known_nonzero_source_gv"
        );
        assert_eq!(residual.predecessor_nonzero, vec![(1, 1)]);
        assert_eq!(residual.difference_nonzero, vec![(0, 2)]);
    }

    #[test]
    fn path_support_lookup_status_compares_diagnostic_gv_to_known_history() {
        assert_eq!(
            path_support_lookup_status("known_nonzero_toric_gv", "-2").unwrap(),
            "path_support_matches_known_nonzero_toric_gv"
        );
        assert_eq!(
            path_support_lookup_status("known_nonzero_toric_gv", "0").unwrap(),
            "path_support_misses_known_nonzero_toric_gv"
        );
        assert_eq!(
            path_support_lookup_status("unknown_not_toric_covered", "3").unwrap(),
            "path_support_nonzero_unknown_not_toric_covered"
        );
        assert_eq!(
            path_support_lookup_status("known_nonzero_source_gv", "3").unwrap(),
            "path_support_matches_known_nonzero_source_gv"
        );
        assert_eq!(
            path_support_lookup_status("known_zero_source_gv", "0").unwrap(),
            "path_support_matches_known_zero_source_gv"
        );
        assert_eq!(
            path_support_lookup_status("known_zero_source_gv", "3").unwrap(),
            "path_support_contradicts_known_zero_source_gv"
        );
        assert_eq!(
            path_support_lookup_status("unknown_not_toric_covered", "0").unwrap(),
            "path_support_zero_or_absent_unknown_not_toric_covered"
        );
    }

    #[test]
    fn sparse_matches_dense_rejects_missing_extra_or_duplicate_entries() {
        assert!(sparse_matches_dense(&[(0, 2), (3, -1)], &[2, 0, 0, -1]));
        assert!(!sparse_matches_dense(&[(0, 2)], &[2, 0, 0, -1]));
        assert!(!sparse_matches_dense(&[(0, 2), (4, 1)], &[2, 0, 0, -1]));
        assert!(!sparse_matches_dense(&[(0, 2), (0, 2)], &[2, 0, 0, -1]));
        assert!(!sparse_matches_dense(&[(0, 0)], &[0]));
        assert_eq!(
            ambient_origin_relation_pattern(&[(0, -1), (5, -2), (6, 1), (7, 1)]).as_deref(),
            Some("origin=-1;neg={-2: 1};pos={1: 2}")
        );
        assert_eq!(ambient_origin_relation_pattern(&[(5, -2), (6, 1)]), None);
    }

    #[test]
    fn active_decomposition_source_status_counts_classify_generator_dependencies() {
        let target_sample = MissingGvTargetSample {
            degree: 2,
            generators_le_degree: 5,
            is_mori_generator: false,
            origin_circuit_pattern: None,
            origin_circuit_witness_count: None,
            origin_circuit_first_witness: None,
            origin_circuit_witnesses: None,
            origin_circuit_affine_support: None,
            cms_general_divisor_shape_candidates: None,
            cms_general_divisor_intersection_checks: None,
            branch_diagnostic: None,
            real_cone_decomposable_by_other_generators: true,
            real_cone_decomposition_active_generators: Some(6),
            real_cone_decomposition_active_generator_basis_nonzero: Some(vec![
                vec![(0, 1)],
                vec![(0, 1), (1, 1)],
                vec![(1, 1)],
                vec![(1, 2)],
                vec![(2, 1)],
                vec![(0, 1), (2, 1)],
            ]),
            real_cone_decomposition_exact_coefficients: Some(vec![
                "1".to_string(),
                "1".to_string(),
                "1".to_string(),
                "1".to_string(),
                "1".to_string(),
                "1".to_string(),
            ]),
            real_cone_decomposition_exact_kind: Some("integer_semigroup".to_string()),
            ambient_nonzero: Vec::new(),
            basis_nonzero: vec![(0, 1), (1, 1)],
        };
        let uncovered_sample = MissingGvTargetSample {
            degree: 1,
            generators_le_degree: 1,
            is_mori_generator: false,
            origin_circuit_pattern: None,
            origin_circuit_witness_count: None,
            origin_circuit_first_witness: None,
            origin_circuit_witnesses: None,
            origin_circuit_affine_support: None,
            cms_general_divisor_shape_candidates: None,
            cms_general_divisor_intersection_checks: None,
            branch_diagnostic: None,
            real_cone_decomposable_by_other_generators: false,
            real_cone_decomposition_active_generators: None,
            real_cone_decomposition_active_generator_basis_nonzero: None,
            real_cone_decomposition_exact_coefficients: None,
            real_cone_decomposition_exact_kind: None,
            ambient_nonzero: Vec::new(),
            basis_nonzero: vec![(1, 1)],
        };
        let second_uncovered_sample = MissingGvTargetSample {
            degree: 2,
            generators_le_degree: 1,
            is_mori_generator: false,
            origin_circuit_pattern: None,
            origin_circuit_witness_count: None,
            origin_circuit_first_witness: Some(OriginCircuitWitnessSample {
                first_facet_exclusive_point: 0,
                second_facet_exclusive_point: 1,
                shared_two_simplex: Vec::new(),
                first_facet: Vec::new(),
                second_facet: Vec::new(),
                first_facet_size: 1,
                second_facet_size: 1,
                sparse_relation: vec![(0, 2), (1, -1)],
                relation_points: vec![
                    OriginCircuitRelationPointSample {
                        point_index: 0,
                        coefficient: 2,
                        coordinates: vec![0],
                        face_dimension: None,
                    },
                    OriginCircuitRelationPointSample {
                        point_index: 1,
                        coefficient: -1,
                        coordinates: vec![1],
                        face_dimension: None,
                    },
                ],
            }),
            origin_circuit_witnesses: None,
            origin_circuit_affine_support: Some(OriginCircuitAffineSupportSample {
                affine_rank: 1,
                coefficient_counts: BTreeMap::from([(-1, 1), (2, 1)]),
                local_charge_basis: vec![vec![2, -1]],
                local_coordinates: vec![
                    OriginCircuitLocalCoordinateSample {
                        point_index: 0,
                        coordinates: vec![0],
                    },
                    OriginCircuitLocalCoordinateSample {
                        point_index: 1,
                        coordinates: vec![1],
                    },
                ],
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
            basis_nonzero: vec![(1, 2)],
        };
        let stats = MissingGvTargetStats {
            target_count: 1,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![target_sample],
        };
        let uncovered_stats = MissingGvTargetStats {
            target_count: 2,
            real_cone_decomposition_exact_kind_counts: HashMap::new(),
            sample: vec![uncovered_sample, second_uncovered_sample],
        };
        let ray_context = vec![DegreeBoundedMoriRayContextSample {
            degree: 1,
            ambient_nonzero: vec![(0, -1), (7, 1)],
            basis_nonzero: vec![(2, 1)],
        }];
        let grading = vec![1, 1, 1];
        let q_matrix = Vec::new();
        let degree_bounded_rays = Vec::new();
        let mut covered_toric_gv_by_basis = HashMap::new();
        covered_toric_gv_by_basis.insert(vec![1, 0, 0], "7".to_string());
        let mut source_derived_gv_by_basis = HashMap::new();
        source_derived_gv_by_basis.insert(vec![0, 1, 0], "5".to_string());
        let context = ValidatedContext {
            dimension: 3,
            degree_bound: 2,
            q_cols: 3,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: Some(&ray_context),
            covered_toric_gv_by_basis,
            source_derived_gv_by_basis,
            intersection: Intersection::new(3),
            stats: &stats,
            uncovered_source_ray_stats: Some(&uncovered_stats),
            shared_facet_unresolved_source_ray_stats: None,
        };

        let counts =
            active_decomposition_generator_source_status_counts(&stats.sample, &context, None);
        assert_eq!(
            counts,
            BTreeMap::from([
                ("active_generator_known_toric_covered".to_string(), 1),
                ("active_generator_known_source_derived_gv".to_string(), 1),
                ("active_generator_matches_missing_target".to_string(), 1),
                (
                    "active_generator_matches_uncovered_source_ray".to_string(),
                    1
                ),
                (
                    "active_generator_source_ray_not_toric_covered".to_string(),
                    1
                ),
                (
                    "active_generator_not_source_degree_bounded_ray".to_string(),
                    1
                ),
            ])
        );

        let summaries =
            active_decomposition_unresolved_source_leaf_summaries(&stats.sample, &context, None);
        assert_eq!(summaries.len(), 4);
        let summary_status_counts = summaries.iter().fold(BTreeMap::new(), |mut acc, summary| {
            *acc.entry(summary.source_status.clone()).or_insert(0) += 1;
            acc
        });
        assert_eq!(
            summary_status_counts,
            BTreeMap::from([
                ("active_generator_matches_missing_target".to_string(), 1),
                (
                    "active_generator_matches_uncovered_source_ray".to_string(),
                    1
                ),
                (
                    "active_generator_source_ray_not_toric_covered".to_string(),
                    1
                ),
                (
                    "active_generator_not_source_degree_bounded_ray".to_string(),
                    1
                ),
            ])
        );
        assert!(summaries.iter().any(|summary| {
            summary.source_status == "active_generator_matches_missing_target"
                && summary.matching_missing_target_index == Some(0)
        }));
        assert!(summaries.iter().any(|summary| {
            summary.source_status == "active_generator_matches_uncovered_source_ray"
                && summary.matching_uncovered_source_ray_index == Some(1)
                && summary
                    .matching_uncovered_source_ray_local_charge_signature
                    .as_deref()
                    == Some("-1,2")
                && summary
                    .matching_uncovered_source_ray_local_cygv_input_skeleton
                    .as_ref()
                    .and_then(|skeleton| skeleton.local_cygv_wrapper_q_matrix_candidate.as_ref())
                    == Some(&vec![vec![2, -1]])
                && summary
                    .matching_uncovered_source_ray_local_unit_phase_probe
                    .is_some()
        }));
        assert!(summaries.iter().any(|summary| {
            summary.source_status == "active_generator_source_ray_not_toric_covered"
                && summary
                    .source_ray_ambient_origin_relation_pattern
                    .as_deref()
                    == Some("origin=-1;neg={};pos={1: 1}")
        }));
        let phase_counts =
            active_decomposition_source_leaf_unit_phase_probe_status_counts(&summaries);
        assert_eq!(phase_counts.values().sum::<usize>(), 4);
        assert_eq!(phase_counts.get("not_available").copied(), Some(3));
        let cms_counts = active_decomposition_source_leaf_cms_solution_status_counts(&summaries);
        assert_eq!(cms_counts.get("no_cms_solution_summary").copied(), Some(4));
    }

    #[test]
    fn path_support_uncovered_source_ray_summary_keeps_unique_source_queue() {
        let mut summaries = BTreeMap::new();
        let lookup = CygvPathSupportLookup {
            candidate_index: 2,
            side: "difference".to_string(),
            degree: 6,
            known_qn_history_status: "unknown_not_toric_covered".to_string(),
            toric_gv: None,
            source_derived_gv: None,
            path_support_gv: Some("1".to_string()),
            path_support_lookup_status: "path_support_nonzero_unknown_not_toric_covered"
                .to_string(),
            source_class_status: "source_ray_not_toric_covered".to_string(),
            source_ray_ambient_nonzero: Some(vec![(0, -1), (5, 1)]),
            matching_missing_target_index: None,
            matching_missing_target_degree: None,
            matching_missing_target_origin_circuit_pattern: None,
            matching_missing_target_exact_kind: None,
            matching_uncovered_source_ray_index: None,
            matching_uncovered_source_ray_degree: None,
            matching_uncovered_source_ray_origin_circuit_pattern: None,
            matching_uncovered_source_ray_exact_kind: None,
            matching_uncovered_source_ray_cms_check_status_counts: BTreeMap::new(),
            matching_uncovered_source_ray_cms_solution_summaries: vec![
                CmsGeneralDivisorSolutionSummary {
                    cms_check_status:
                        "cms_general_divisor_integral_solution_matches_inferred_degree".to_string(),
                    shrinking_divisor_index: 5,
                    solution_basis_nonzero: vec![(0, "1".to_string())],
                    solution_ambient_basis_nonzero: vec![(5, "1".to_string())],
                    computed_other_normal_degree: "0".to_string(),
                    solution_basis_cubic_self_intersection: "0".to_string(),
                    local_intersection_tensor_candidate_status:
                        "candidate_from_cms_divisor_cubic_needs_phase_and_chamber_certificate"
                            .to_string(),
                    local_intersection_tensor_candidate: Some(vec![
                        LocalCygvIntersectionTensorEntry {
                            indices: [0, 0, 0],
                            value: "0".to_string(),
                        },
                    ]),
                    local_cygv_primitive_probe: None,
                },
            ],
            matching_uncovered_source_ray_local_charge_signature: Some("-2,-1,1,1,1".to_string()),
            matching_uncovered_source_ray_local_cygv_readiness: Some(
                "blocked_missing_source_derived_inputs".to_string(),
            ),
            matching_uncovered_source_ray_local_missing_inputs: vec![
                "local_semigroup_generators".to_string(),
                "local_chamber_certificate".to_string(),
            ],
            curve_nonzero: vec![(3, -1), (8, 1)],
        };
        add_path_support_uncovered_source_ray_lookup(&mut summaries, 7, &lookup);
        let repeated_lookup = CygvPathSupportLookup {
            candidate_index: 5,
            side: "predecessor".to_string(),
            path_support_gv: Some("-2".to_string()),
            ..lookup.clone()
        };
        add_path_support_uncovered_source_ray_lookup(&mut summaries, 8, &repeated_lookup);
        let ignored_lookup = CygvPathSupportLookup {
            source_class_status: "not_source_degree_bounded_ray".to_string(),
            curve_nonzero: vec![(9, 1)],
            ..lookup
        };
        add_path_support_uncovered_source_ray_lookup(&mut summaries, 9, &ignored_lookup);

        assert_eq!(summaries.len(), 1);
        let summary = summaries.get(&vec![(3, -1), (8, 1)]).unwrap();
        assert_eq!(summary.degree, 6);
        assert_eq!(summary.occurrences.len(), 2);
        assert_eq!(
            summary.known_qn_history_status_counts,
            BTreeMap::from([("unknown_not_toric_covered".to_string(), 2)])
        );
        assert_eq!(
            summary.path_support_gv_counts,
            BTreeMap::from([("1".to_string(), 1), ("-2".to_string(), 1)])
        );
        assert_eq!(
            summary.uncovered_source_ray_cms_solution_summaries,
            vec![CmsGeneralDivisorSolutionSummary {
                cms_check_status: "cms_general_divisor_integral_solution_matches_inferred_degree"
                    .to_string(),
                shrinking_divisor_index: 5,
                solution_basis_nonzero: vec![(0, "1".to_string())],
                solution_ambient_basis_nonzero: vec![(5, "1".to_string())],
                computed_other_normal_degree: "0".to_string(),
                solution_basis_cubic_self_intersection: "0".to_string(),
                local_intersection_tensor_candidate_status:
                    "candidate_from_cms_divisor_cubic_needs_phase_and_chamber_certificate"
                        .to_string(),
                local_intersection_tensor_candidate: Some(vec![LocalCygvIntersectionTensorEntry {
                    indices: [0, 0, 0],
                    value: "0".to_string(),
                }]),
                local_cygv_primitive_probe: None,
            }]
        );
        assert_eq!(
            summary
                .uncovered_source_ray_local_charge_signature
                .as_deref(),
            Some("-2,-1,1,1,1")
        );
        assert_eq!(
            summary.uncovered_source_ray_local_cygv_readiness_counts,
            BTreeMap::from([("blocked_missing_source_derived_inputs".to_string(), 2)])
        );
        assert_eq!(
            summary.uncovered_source_ray_local_missing_input_counts,
            BTreeMap::from([
                ("local_semigroup_generators".to_string(), 2),
                ("local_chamber_certificate".to_string(), 2)
            ])
        );
        assert_eq!(summary.occurrences[0].target_index, 7);
        assert_eq!(summary.occurrences[1].side, "predecessor");
    }

    #[test]
    fn path_support_primitive_probe_status_counts_aggregate_solution_summaries() {
        let summary = CygvPathSupportSourceRaySummary {
            degree: 6,
            occurrence_count: 1,
            known_qn_history_status_counts: BTreeMap::new(),
            path_support_lookup_status_counts: BTreeMap::new(),
            path_support_gv_counts: BTreeMap::new(),
            uncovered_source_ray_origin_circuit_pattern: None,
            uncovered_source_ray_exact_kind: None,
            uncovered_source_ray_cms_check_status_counts: BTreeMap::new(),
            uncovered_source_ray_cms_solution_summaries: vec![
                CmsGeneralDivisorSolutionSummary {
                    cms_check_status:
                        "cms_general_divisor_integral_solution_matches_inferred_degree".to_string(),
                    shrinking_divisor_index: 1,
                    solution_basis_nonzero: vec![(0, "1".to_string())],
                    solution_ambient_basis_nonzero: vec![(1, "1".to_string())],
                    computed_other_normal_degree: "0".to_string(),
                    solution_basis_cubic_self_intersection: "3".to_string(),
                    local_intersection_tensor_candidate_status: String::new(),
                    local_intersection_tensor_candidate: None,
                    local_cygv_primitive_probe: Some(LocalCygvPrimitiveProbe {
                        status: "primitive_cygv_probe_mismatch_raw_cubic_is_not_certified_tensor"
                            .to_string(),
                        q_matrix: vec![vec![-1, -2, 1, 1, 1]],
                        grading_vector: vec![1],
                        semigroup_elements: vec![vec![0], vec![1]],
                        candidate_gv: Some("-6".to_string()),
                        unit_tensor_candidate_gv: Some("-2".to_string()),
                        unit_tensor_probe_status:
                            "unit_tensor_probe_matches_expected_formula_but_chamber_uncertified"
                                .to_string(),
                        origin_omitted_unit_tensor_candidate_gv: None,
                        origin_omitted_unit_tensor_probe_status:
                            "origin_omitted_unit_tensor_probe_not_run_no_origin_column".to_string(),
                        expected_toric_gv1_formula_value: Some("-2".to_string()),
                        error: None,
                        unit_tensor_error: None,
                        origin_omitted_unit_tensor_error: None,
                    }),
                },
                CmsGeneralDivisorSolutionSummary {
                    cms_check_status:
                        "cms_general_divisor_integral_solution_matches_inferred_degree".to_string(),
                    shrinking_divisor_index: 2,
                    solution_basis_nonzero: vec![(0, "1".to_string())],
                    solution_ambient_basis_nonzero: vec![(2, "1".to_string())],
                    computed_other_normal_degree: "0".to_string(),
                    solution_basis_cubic_self_intersection: "0".to_string(),
                    local_intersection_tensor_candidate_status: String::new(),
                    local_intersection_tensor_candidate: None,
                    local_cygv_primitive_probe: None,
                },
            ],
            uncovered_source_ray_local_charge_signature: None,
            uncovered_source_ray_local_cygv_readiness_counts: BTreeMap::new(),
            uncovered_source_ray_local_missing_input_counts: BTreeMap::new(),
            source_ray_ambient_nonzero: Vec::new(),
            curve_nonzero: Vec::new(),
            occurrences: Vec::new(),
        };

        assert_eq!(
            path_support_uncovered_source_ray_local_cygv_primitive_probe_status_counts(&[summary]),
            BTreeMap::from([
                ("local_cygv_primitive_probe_not_run".to_string(), 1),
                (
                    "primitive_cygv_probe_mismatch_raw_cubic_is_not_certified_tensor".to_string(),
                    1
                )
            ])
        );
    }

    #[test]
    fn cms_general_divisor_check_status_tracks_stable_weyl_readiness() {
        let no_solution = CmsGeneralDivisorIntersectionCheck {
            shrinking_divisor_index: 1,
            has_rational_divisor_solution: false,
            solution_basis_support_len: None,
            solution_basis_nonzero: None,
            solution_ambient_basis_nonzero: None,
            solution_is_integral: None,
            computed_other_normal_degree: None,
            matches_inferred_other_normal_degree: None,
        };
        let matching_solution = CmsGeneralDivisorIntersectionCheck {
            shrinking_divisor_index: 2,
            has_rational_divisor_solution: true,
            solution_basis_support_len: Some(3),
            solution_basis_nonzero: Some(vec![(0, "1".to_string())]),
            solution_ambient_basis_nonzero: Some(vec![(4, "1".to_string())]),
            solution_is_integral: Some(true),
            computed_other_normal_degree: Some("1".to_string()),
            matches_inferred_other_normal_degree: Some(true),
        };

        assert_eq!(
            cms_general_divisor_intersection_check_status(&no_solution),
            "cms_general_divisor_no_rational_divisor_solution"
        );
        assert_eq!(
            cms_general_divisor_intersection_check_status(&matching_solution),
            "cms_general_divisor_integral_solution_matches_inferred_degree"
        );
    }

    #[test]
    fn divisor_cubic_self_intersection_contracts_solution_with_kappa() {
        let mut intersection = Intersection::new(2);
        intersection.set(0, 0, 0, Rational::<Finite>::new(MalachiteRational::from(2)));
        intersection.set(0, 0, 1, Rational::<Finite>::new(MalachiteRational::from(3)));
        intersection.set(0, 1, 1, Rational::<Finite>::new(MalachiteRational::from(5)));
        intersection.set(1, 1, 1, Rational::<Finite>::new(MalachiteRational::from(7)));

        let cubic = divisor_cubic_self_intersection(
            &intersection,
            &[(0, "2".to_string()), (1, "-1".to_string())],
        )
        .unwrap();

        assert_eq!(cubic, "3");
    }

    #[test]
    fn cms_solution_summary_adds_one_parameter_intersection_tensor_candidate() {
        let mut sample = minimal_missing_sample(vec![(0, 1)]);
        sample.origin_circuit_first_witness = Some(OriginCircuitWitnessSample {
            first_facet_exclusive_point: 2,
            second_facet_exclusive_point: 4,
            shared_two_simplex: vec![1, 3],
            first_facet: Vec::new(),
            second_facet: Vec::new(),
            first_facet_size: 2,
            second_facet_size: 2,
            sparse_relation: vec![(0, -1), (1, 2), (2, 1), (3, -3), (4, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 1,
                    coefficient: 2,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 2,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 3,
                    coefficient: -3,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
            ],
        });
        sample.origin_circuit_affine_support = Some(OriginCircuitAffineSupportSample {
            affine_rank: 3,
            coefficient_counts: BTreeMap::new(),
            local_charge_basis: vec![vec![1, -2, -1, 3, -1]],
            local_coordinates: Vec::new(),
            local_coordinates_2d: None,
        });
        sample.cms_general_divisor_intersection_checks =
            Some(vec![CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 4,
                has_rational_divisor_solution: true,
                solution_basis_support_len: Some(1),
                solution_basis_nonzero: Some(vec![(0, "1".to_string())]),
                solution_ambient_basis_nonzero: Some(vec![(4, "1".to_string())]),
                solution_is_integral: Some(true),
                computed_other_normal_degree: Some("0".to_string()),
                matches_inferred_other_normal_degree: Some(true),
            }]);
        let mut intersection = Intersection::new(1);
        intersection.set(0, 0, 0, Rational::<Finite>::new(MalachiteRational::from(3)));

        let summaries =
            cms_general_divisor_solution_summaries(Some(&sample), &intersection).unwrap();

        assert_eq!(summaries.len(), 1);
        assert_eq!(summaries[0].solution_basis_cubic_self_intersection, "3");
        assert_eq!(
            summaries[0].local_intersection_tensor_candidate_status,
            "candidate_from_cms_divisor_cubic_needs_phase_and_chamber_certificate"
        );
        assert_eq!(
            summaries[0].local_intersection_tensor_candidate,
            Some(vec![LocalCygvIntersectionTensorEntry {
                indices: [0, 0, 0],
                value: "3".to_string(),
            }])
        );
    }

    #[test]
    fn cms_solution_summary_keeps_non_promotable_rational_solution() {
        let mut sample = minimal_missing_sample(vec![(0, 1)]);
        sample.origin_circuit_first_witness = Some(OriginCircuitWitnessSample {
            first_facet_exclusive_point: 2,
            second_facet_exclusive_point: 4,
            shared_two_simplex: vec![1, 3],
            first_facet: Vec::new(),
            second_facet: Vec::new(),
            first_facet_size: 2,
            second_facet_size: 2,
            sparse_relation: vec![(0, -1), (1, 2), (2, 1), (3, -3), (4, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 1,
                    coefficient: 2,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 2,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 3,
                    coefficient: -3,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
            ],
        });
        sample.origin_circuit_affine_support = Some(OriginCircuitAffineSupportSample {
            affine_rank: 3,
            coefficient_counts: BTreeMap::new(),
            local_charge_basis: vec![vec![1, -2, -1, 3, -1]],
            local_coordinates: Vec::new(),
            local_coordinates_2d: None,
        });
        sample.cms_general_divisor_intersection_checks =
            Some(vec![CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 4,
                has_rational_divisor_solution: true,
                solution_basis_support_len: Some(1),
                solution_basis_nonzero: Some(vec![(0, "1".to_string())]),
                solution_ambient_basis_nonzero: Some(vec![(4, "1".to_string())]),
                solution_is_integral: Some(true),
                computed_other_normal_degree: Some("-2".to_string()),
                matches_inferred_other_normal_degree: Some(false),
            }]);
        let mut intersection = Intersection::new(1);
        intersection.set(0, 0, 0, Rational::<Finite>::new(MalachiteRational::from(3)));

        let summaries =
            cms_general_divisor_solution_summaries(Some(&sample), &intersection).unwrap();

        assert_eq!(summaries.len(), 1);
        assert_eq!(
            summaries[0].cms_check_status,
            "cms_general_divisor_integral_solution_mismatches_inferred_degree"
        );
        assert_eq!(
            summaries[0].local_intersection_tensor_candidate_status,
            "diagnostic_from_cms_general_divisor_integral_solution_mismatches_inferred_degree_not_promoted"
        );
    }

    #[test]
    fn cms_solution_summary_probe_flags_raw_cubic_mismatch() {
        let mut sample = minimal_missing_sample(vec![(0, 1)]);
        sample.origin_circuit_first_witness = Some(OriginCircuitWitnessSample {
            first_facet_exclusive_point: 3,
            second_facet_exclusive_point: 4,
            shared_two_simplex: vec![1, 2],
            first_facet: Vec::new(),
            second_facet: Vec::new(),
            first_facet_size: 2,
            second_facet_size: 2,
            sparse_relation: vec![(0, -1), (1, -2), (2, 1), (3, 1), (4, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 1,
                    coefficient: -2,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 2,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 3,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
            ],
        });
        sample.origin_circuit_affine_support = Some(OriginCircuitAffineSupportSample {
            affine_rank: 3,
            coefficient_counts: BTreeMap::new(),
            local_charge_basis: vec![vec![1, 2, -1, -1, -1]],
            local_coordinates: Vec::new(),
            local_coordinates_2d: None,
        });
        sample.cms_general_divisor_shape_candidates = Some(vec![CmsGeneralDivisorShapeCandidate {
            shrinking_divisor_index: 1,
            shrinking_divisor_coefficient: -2,
            shrinking_divisor_coordinates: Vec::new(),
            inferred_other_normal_degree: 0,
            toric_gv1_formula_value: Some(-2),
            all_non_origin_relation_points_are_two_face: false,
        }]);
        sample.cms_general_divisor_intersection_checks =
            Some(vec![CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 1,
                has_rational_divisor_solution: true,
                solution_basis_support_len: Some(1),
                solution_basis_nonzero: Some(vec![(0, "1".to_string())]),
                solution_ambient_basis_nonzero: Some(vec![(1, "1".to_string())]),
                solution_is_integral: Some(true),
                computed_other_normal_degree: Some("0".to_string()),
                matches_inferred_other_normal_degree: Some(true),
            }]);
        let mut intersection = Intersection::new(1);
        intersection.set(0, 0, 0, Rational::<Finite>::new(MalachiteRational::from(3)));

        let summaries =
            cms_general_divisor_solution_summaries(Some(&sample), &intersection).unwrap();
        let probe = summaries[0]
            .local_cygv_primitive_probe
            .as_ref()
            .expect("one-parameter CMS summary should run primitive cygv probe");

        assert_eq!(
            probe.status,
            "primitive_cygv_probe_mismatch_raw_cubic_is_not_certified_tensor"
        );
        assert_eq!(probe.q_matrix, vec![vec![-1, -2, 1, 1, 1]]);
        assert_eq!(probe.candidate_gv.as_deref(), Some("-6"));
        assert_eq!(probe.unit_tensor_candidate_gv.as_deref(), Some("-2"));
        assert_eq!(
            probe.unit_tensor_probe_status,
            "unit_tensor_probe_matches_expected_formula_but_chamber_uncertified"
        );
        assert_eq!(
            probe.origin_omitted_unit_tensor_probe_status,
            "origin_omitted_unit_tensor_probe_hkty_error"
        );
        assert!(
            probe
                .origin_omitted_unit_tensor_error
                .as_deref()
                .is_some_and(|error| error.contains("dimension of the CY must be at least three"))
        );
        assert_eq!(
            probe.expected_toric_gv1_formula_value.as_deref(),
            Some("-2")
        );
    }

    #[test]
    fn cms_solution_summary_origin_omitted_probe_checks_origin_phase() {
        let mut sample = minimal_missing_sample(vec![(0, 1)]);
        sample.origin_circuit_first_witness = Some(OriginCircuitWitnessSample {
            first_facet_exclusive_point: 3,
            second_facet_exclusive_point: 4,
            shared_two_simplex: vec![1, 2],
            first_facet: Vec::new(),
            second_facet: Vec::new(),
            first_facet_size: 2,
            second_facet_size: 2,
            sparse_relation: vec![(0, -1), (1, 1), (2, -1), (3, 1), (4, -1), (5, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 1,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 2,
                    coefficient: -1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 3,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: -1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 5,
                    coefficient: 1,
                    coordinates: Vec::new(),
                    face_dimension: None,
                },
            ],
        });
        sample.origin_circuit_affine_support = Some(OriginCircuitAffineSupportSample {
            affine_rank: 4,
            coefficient_counts: BTreeMap::new(),
            local_charge_basis: vec![vec![1, -1, 1, -1, 1, -1]],
            local_coordinates: Vec::new(),
            local_coordinates_2d: None,
        });
        sample.cms_general_divisor_shape_candidates = Some(vec![CmsGeneralDivisorShapeCandidate {
            shrinking_divisor_index: 1,
            shrinking_divisor_coefficient: 1,
            shrinking_divisor_coordinates: Vec::new(),
            inferred_other_normal_degree: -1,
            toric_gv1_formula_value: Some(1),
            all_non_origin_relation_points_are_two_face: false,
        }]);
        sample.cms_general_divisor_intersection_checks =
            Some(vec![CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 1,
                has_rational_divisor_solution: true,
                solution_basis_support_len: Some(1),
                solution_basis_nonzero: Some(vec![(0, "1".to_string())]),
                solution_ambient_basis_nonzero: Some(vec![(1, "1".to_string())]),
                solution_is_integral: Some(true),
                computed_other_normal_degree: Some("-1".to_string()),
                matches_inferred_other_normal_degree: Some(true),
            }]);
        let mut intersection = Intersection::new(1);
        intersection.set(
            0,
            0,
            0,
            Rational::<Finite>::new(MalachiteRational::from(26)),
        );

        let summaries =
            cms_general_divisor_solution_summaries(Some(&sample), &intersection).unwrap();
        let probe = summaries[0]
            .local_cygv_primitive_probe
            .as_ref()
            .expect("one-parameter CMS summary should run primitive cygv probe");

        assert_eq!(probe.q_matrix, vec![vec![1, -1, 1, -1, 1]]);
        assert_eq!(probe.candidate_gv.as_deref(), Some("26"));
        assert_eq!(probe.unit_tensor_candidate_gv.as_deref(), Some("1"));
        assert_eq!(
            probe.origin_omitted_unit_tensor_candidate_gv.as_deref(),
            Some("1")
        );
        assert_eq!(
            probe.origin_omitted_unit_tensor_probe_status,
            "origin_omitted_unit_tensor_probe_matches_expected_formula_but_phase_uncertified"
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
                origin_circuit_witnesses: None,
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
        let mut covered_toric_gv_by_basis = HashMap::new();
        covered_toric_gv_by_basis.insert(vec![1, 0], "7".to_string());
        covered_toric_gv_by_basis.insert(vec![0, 1], "11".to_string());
        let context = ValidatedContext {
            dimension: 2,
            degree_bound: 2,
            q_cols: 2,
            grading: &grading,
            q_matrix: &q_matrix,
            degree_bounded_rays: &degree_bounded_rays,
            degree_bounded_ray_context: None,
            covered_toric_gv_by_basis,
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
        };

        let target = vec![1, 1];
        let probe = cygv_path_history_probe_inner(
            &stats.sample[0],
            &context,
            &target,
            false,
            false,
            16,
            None,
            None,
        )
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
        assert_eq!(
            probe
                .predecessor_toric_coverage_counts
                .get("both_toric_covered")
                .copied(),
            Some(2)
        );
        assert_eq!(
            probe
                .predecessor_known_qn_history_counts
                .get("predecessor_known_nonzero_toric_gv__difference_known_nonzero_toric_gv")
                .copied(),
            Some(2)
        );
        assert_eq!(probe.closest_series_distance.as_deref(), Some("1.000000"));
        assert_eq!(probe.closest_series_predecessor_nonzero, Some(vec![(1, 1)]));
        assert_eq!(probe.closest_series_difference_nonzero, Some(vec![(0, 1)]));
        let closest_known = probe.closest_known_qn_predecessor.as_ref().unwrap();
        assert_eq!(closest_known.series_distance, "1.000000");
        assert_eq!(closest_known.predecessor_nonzero, vec![(1, 1)]);
        assert_eq!(closest_known.difference_nonzero, vec![(0, 1)]);
        assert_eq!(
            closest_known.predecessor_known_qn_history_status,
            "known_nonzero_toric_gv"
        );
        assert_eq!(
            probe.predecessor_candidate_sample_limit,
            CYGV_PATH_PREDECESSOR_SAMPLE_LIMIT
        );
        assert_eq!(probe.predecessor_candidate_sample.len(), 2);
        assert_eq!(
            probe.predecessor_candidate_sample[0].predecessor_nonzero,
            vec![(0, 1)]
        );
        assert_eq!(
            probe.predecessor_candidate_sample[0].difference_nonzero,
            vec![(1, 1)]
        );
        assert_eq!(
            probe.predecessor_candidate_sample[0]
                .predecessor_toric_gv
                .as_deref(),
            Some("7")
        );
        assert_eq!(
            probe.predecessor_candidate_sample[0]
                .difference_toric_gv
                .as_deref(),
            Some("11")
        );
        assert_eq!(
            probe.predecessor_candidate_sample[0].predecessor_known_qn_history_status,
            "known_nonzero_toric_gv"
        );
        assert_eq!(
            probe.predecessor_candidate_sample[0].difference_known_qn_history_status,
            "known_nonzero_toric_gv"
        );
        assert_eq!(
            probe.predecessor_candidate_sample[0].known_qn_history_pair_status,
            "predecessor_known_nonzero_toric_gv__difference_known_nonzero_toric_gv"
        );
        assert!(probe.predecessor_candidate_sample[0].predecessor_is_seed);
        assert!(probe.predecessor_candidate_sample[0].difference_is_seed);
        assert!(probe.predecessor_candidate_sample[0].predecessor_is_reduced_seed);
        assert!(probe.predecessor_candidate_sample[0].difference_is_reduced_seed);
        assert!(
            probe.predecessor_candidate_sample[0]
                .predecessor_first_generation_seed_sum
                .is_none()
        );
        assert!(
            probe.predecessor_candidate_sample[0]
                .difference_first_generation_seed_sum
                .is_none()
        );
        assert_eq!(
            probe.predecessor_candidate_sample[1].predecessor_nonzero,
            vec![(1, 1)]
        );
        assert_eq!(
            probe.predecessor_candidate_sample[1].difference_nonzero,
            vec![(0, 1)]
        );
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
                origin_circuit_witnesses: None,
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
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
        };

        let target = vec![1, 1];
        let probe = cygv_path_history_probe_inner(
            &stats.sample[0],
            &context,
            &target,
            false,
            false,
            2,
            None,
            None,
        )
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
                origin_circuit_witnesses: None,
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
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection: Intersection::new(2),
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
        };

        let target = vec![1, 1];
        let probe = cygv_path_history_probe_inner(
            &stats.sample[0],
            &context,
            &target,
            false,
            false,
            16,
            Some(2),
            None,
        )
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
            origin_circuit_witnesses: None,
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
    fn origin_circuit_support_diagnostics_use_all_serialized_witnesses() {
        let missing_witness = OriginCircuitWitnessSample {
            first_facet_exclusive_point: 4,
            second_facet_exclusive_point: 5,
            shared_two_simplex: vec![1, 2],
            first_facet: Vec::new(),
            second_facet: Vec::new(),
            first_facet_size: 3,
            second_facet_size: 3,
            sparse_relation: vec![(0, -1), (4, 1), (5, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: vec![0, 0],
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: vec![1, 0],
                    face_dimension: None,
                },
            ],
        };
        let full_witness = OriginCircuitWitnessSample {
            first_facet_exclusive_point: 6,
            second_facet_exclusive_point: 7,
            shared_two_simplex: vec![1, 2],
            first_facet: vec![1, 2, 6],
            second_facet: vec![1, 2, 7],
            first_facet_size: 3,
            second_facet_size: 3,
            sparse_relation: vec![(0, -1), (6, 1), (7, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: vec![0, 0],
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 7,
                    coefficient: 1,
                    coordinates: vec![0, 1],
                    face_dimension: None,
                },
            ],
        };
        let mut sample = minimal_missing_sample(Vec::new());
        sample.origin_circuit_first_witness = Some(missing_witness.clone());
        sample.origin_circuit_witnesses = Some(vec![missing_witness, full_witness]);

        assert_eq!(
            origin_circuit_facet_context_status(&sample),
            "mixed_origin_circuit_facet_context:origin_circuit_missing_full_facet_context|source_derived_full_facet_context"
        );
        let supports = origin_circuit_ambient_support_sets(&sample)
            .expect("witness list should define ambient supports");
        assert_eq!(supports.relation_support, HashSet::from([0, 4, 7]));
        assert_eq!(supports.shared_facet, HashSet::from([0, 1, 2, 4, 5, 6, 7]));
        assert_eq!(supports.facet_union, HashSet::from([0, 1, 2, 6, 7]));
    }

    #[test]
    fn origin_circuit_witness_domain_supports_are_per_witness() {
        let witness = OriginCircuitWitnessSample {
            first_facet_exclusive_point: 4,
            second_facet_exclusive_point: 5,
            shared_two_simplex: vec![1, 2],
            first_facet: vec![1, 2, 4, 6],
            second_facet: vec![1, 2, 5, 7],
            first_facet_size: 4,
            second_facet_size: 4,
            sparse_relation: vec![(0, -2), (4, 1), (5, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -2,
                    coordinates: vec![0, 0],
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: vec![1, 0],
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 5,
                    coefficient: 1,
                    coordinates: vec![0, 1],
                    face_dimension: None,
                },
            ],
        };

        let supports = origin_circuit_witness_ambient_support_sets(&witness);

        assert_eq!(supports.relation_support, HashSet::from([0, 4, 5]));
        assert_eq!(supports.shared_facet, HashSet::from([0, 1, 2, 4, 5]));
        assert_eq!(supports.facet_union, HashSet::from([0, 1, 2, 4, 5, 6, 7]));
    }

    #[test]
    fn origin_circuit_witness_relation_status_tracks_multi_witness_relation_changes() {
        let first = OriginCircuitWitnessSample {
            first_facet_exclusive_point: 4,
            second_facet_exclusive_point: 5,
            shared_two_simplex: vec![1, 2],
            first_facet: vec![1, 2, 4],
            second_facet: vec![1, 2, 5],
            first_facet_size: 3,
            second_facet_size: 3,
            sparse_relation: vec![(0, -1), (4, 1), (5, 1)],
            relation_points: vec![
                OriginCircuitRelationPointSample {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: vec![0, 0],
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: vec![1, 0],
                    face_dimension: None,
                },
                OriginCircuitRelationPointSample {
                    point_index: 5,
                    coefficient: 1,
                    coordinates: vec![0, 1],
                    face_dimension: None,
                },
            ],
        };
        let mut same_relation = first.clone();
        same_relation.shared_two_simplex = vec![2, 3];
        same_relation.first_facet = vec![2, 3, 4];
        same_relation.second_facet = vec![2, 3, 5];

        let mut shared = minimal_missing_sample(Vec::new());
        shared.origin_circuit_witness_count = Some(2);
        shared.origin_circuit_first_witness = Some(first.clone());
        shared.origin_circuit_witnesses = Some(vec![first.clone(), same_relation]);
        assert_eq!(
            origin_circuit_witness_relation_status(&shared),
            "all_origin_circuit_witnesses_share_relation"
        );

        let mut changed_relation = first.clone();
        changed_relation.relation_points[2].point_index = 6;
        changed_relation.sparse_relation = vec![(0, -1), (4, 1), (6, 1)];
        let mut mixed = minimal_missing_sample(Vec::new());
        mixed.origin_circuit_witness_count = Some(2);
        mixed.origin_circuit_first_witness = Some(first.clone());
        mixed.origin_circuit_witnesses = Some(vec![first.clone(), changed_relation]);
        assert_eq!(
            origin_circuit_witness_relation_status(&mixed),
            "mixed_origin_circuit_witness_relations"
        );

        let mut legacy = minimal_missing_sample(Vec::new());
        legacy.origin_circuit_witness_count = Some(2);
        legacy.origin_circuit_first_witness = Some(first);
        assert_eq!(
            origin_circuit_witness_relation_status(&legacy),
            "multiple_origin_circuit_witnesses_not_serialized"
        );

        let counts = origin_circuit_witness_relation_status_counts(&[shared, mixed, legacy], None);
        assert_eq!(
            counts.get("all_origin_circuit_witnesses_share_relation"),
            Some(&1)
        );
        assert_eq!(
            counts.get("mixed_origin_circuit_witness_relations"),
            Some(&1)
        );
        assert_eq!(
            counts.get("multiple_origin_circuit_witnesses_not_serialized"),
            Some(&1)
        );
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
                origin_circuit_first_witness: Some(origin_witness.clone()),
                origin_circuit_witnesses: Some(vec![origin_witness]),
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
                        solution_basis_nonzero: Some(vec![
                            (0, "1".to_string()),
                            (1, "1/2".to_string()),
                        ]),
                        solution_ambient_basis_nonzero: Some(vec![
                            (7, "1".to_string()),
                            (8, "1/2".to_string()),
                        ]),
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
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection,
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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
            false,
            None,
            None,
            None,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            256,
            None,
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
            report.origin_circuit_witnesses.as_ref().map(Vec::len),
            Some(1)
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
                origin_circuit_witnesses: None,
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
            covered_toric_gv_by_basis: HashMap::new(),
            source_derived_gv_by_basis: HashMap::new(),
            intersection,
            stats: &stats,
            uncovered_source_ray_stats: None,
            shared_facet_unresolved_source_ray_stats: None,
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
            false,
            None,
            None,
            None,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            256,
            None,
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
            false,
            None,
            None,
            None,
            &mut semigroup_measurement_cache,
            &mut semigroup_ladder_cache,
            256,
            None,
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
