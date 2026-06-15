//! Core mathematical primitives for Calabi-Yau manifold computations.
//!
//! This crate provides the foundational algorithms for:
//! - Lattice polytope operations
//! - Intersection number computation (`κ_ijk`)
//! - Kähler cone and volume calculations
//! - Flat direction and `e^{K₀}` computation
//!
//! ## Volume Computation
//!
//! The string frame volume is computed as:
//! ```text
//! V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
//! ```
//!
//! where BBHL = ζ(3) χ(X) / (4(2π)³) is the α' correction.
//!
//! ## Flat Direction
//!
//! Given flux vectors K and M, the flat direction p and Kähler potential factor are:
//! ```text
//! N_ab = κ_abc M^c
//! p = N⁻¹ K
//! e^{K₀} = (4/3 × κ_abc p^a p^b p^c)⁻¹
//! ```

pub mod basis;
pub mod cone;
pub mod config;
pub mod cosmology;
pub mod curve_basis;
pub mod divisor;
pub mod error;
pub mod flat_direction;
pub mod geometry;
pub mod glsm;
pub mod gv;
pub mod height_kahler;
pub mod integer_math;
pub mod intersection;
pub mod kahler;
pub mod lvs;
pub use cone::Cone;
pub use kahler::{MoriCone, compute_mori_generators};
pub mod chamber_search;
pub mod lattice;
pub mod orientifold;
pub mod pipeline;
pub mod policy;
pub mod polytope;
pub mod quintessence;
pub use pipeline::{EvaluationRequest, EvaluationResult, evaluate_vacuum};
pub use policy::{Abort, ForGA, Strict, VacuumPolicy, VolumePolicy};
pub mod racetrack;
pub mod triangulation;
pub mod types;
pub mod utils;
pub mod vacuum;
pub mod volume;

pub use types::{F64, Finite, H11, H21, I32, I64, Neg, Pos};

pub use basis::{basis_change_matrix, compute_divisor_basis, intersection_in_basis, is_unimodular};
pub use curve_basis::{
    DivisorBasis, apply_finite_f64_basis_transform, apply_integer_basis_transform,
    apply_integer_basis_transform_transpose, compute_curve_basis_matrix,
    compute_curve_basis_matrix_for_divisor_basis,
    compute_curve_basis_matrix_from_divisor_basis_matrix, curve_basis_matrix_without_origin_i64,
    curve_basis_q_matrix_for_divisor_basis_i64, divisor_basis_change_matrix,
    divisor_basis_glsm_coordinate_matrix,
};
pub use divisor::{compute_divisor_jacobian, compute_divisor_volumes, compute_kklt_divisor_chi};
pub use error::{Error, Result};
pub use flat_direction::{
    FlatDirectionResult, compute_ek0, compute_flat_direction, compute_flat_direction_full,
    compute_n_matrix, solve_linear_system,
};
pub use glsm::compute_glsm_and_linrels;
pub use glsm::compute_glsm_charge_matrix;
pub use glsm::compute_glsm_linear_relations;
pub use gv::{
    AffineCircuitRelationPoint, AffineToricCircuitDiagnostic, CkyzLocalCausalDomainSpec,
    CkyzLocalDomainProfile, CkyzLocalIntersectionTerm, CkyzLocalSurfaceIdentification,
    CkyzLocalSurfaceKind, CurveDecompositionTerm, CurvePruningStrategy, ExtremalMoriRayCertificate,
    FiniteCutoffGvChargePartition, FiniteGvTableNopClassification, GvDivisorBasisData,
    LocalToricCircuitKind, LocalToricCoordinate, LocalToricCoordinate2D, NilpotentRayCandidate,
    NilpotentRayDegreeSlice, NilpotentRayDivergenceCheck, NilpotentRaySliceComparisonPoint,
    NilpotentRaySliceDistance, NilpotentRayTwoPassNopClassification, OneDimensionalRayGvSeries,
    OriginCircuitCurveDiagnostic, OriginCircuitCurveWitness, OriginCircuitRelationPoint,
    PotentRayConvergence, RankTwoLocalChargeModel, RankTwoLocalChargeModelPoint,
    RankTwoLocalSupportSignature, RankTwoLocalSupportSignatureEntry, SupportingMoriFace,
    SupportingMoriFaceCertificate, SupportingMoriFaceLpSearchOptions, ToricCurveCandidate,
    ToricCurveGvDiagnostic, ToricCurveGvInvariant, ToricCurveGvSource,
    certify_supporting_mori_face_by_lp_search, check_extremal_mori_ray_separator,
    check_supporting_mori_face_normal, ckyz_local_domain_profile_for_degrees,
    ckyz_local_domain_profile_for_degrees_with_causal_domain,
    ckyz_local_surface_causal_domain_spec, ckyz_local_surface_cover_weight_coefficients,
    ckyz_local_surface_domain_profile_for_multiples, ckyz_local_surface_target_degrees,
    classify_nilpotent_rays_from_two_pass_divergence_checks,
    classify_nop_rays_from_finite_gv_table, compute_ambient_one_dimensional_ray_gv_series,
    compute_ckyz_flat_prepotential_period_corrections, compute_ckyz_inverse_mirror_map,
    compute_ckyz_local_gv_invariants, compute_ckyz_local_gv_invariants_for_degrees,
    compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain,
    compute_ckyz_local_gv_invariants_for_degrees_with_predicted_support_domain,
    compute_ckyz_local_instanton_potential_corrections,
    compute_ckyz_local_prepotential_period_corrections,
    compute_ckyz_local_surface_gv_invariants_for_multiples_with_causal_domain,
    compute_ckyz_log_period_corrections, compute_grading_vector, compute_gv_invariants,
    compute_gv_invariants_with_degree_bounded_lattice,
    compute_gv_invariants_with_explicit_semigroup,
    compute_gv_invariants_with_explicit_semigroup_and_nef_partition,
    compute_gv_invariants_with_provided_generators,
    compute_gv_invariants_with_provided_generators_and_nef_partition,
    compute_local_p2_genus_zero_gv_series, compute_local_toric_circuit_gv_series,
    compute_mori_cone_cap_rays, compute_one_dimensional_ray_gv_series,
    compute_origin_circuit_curve_diagnostics, compute_ray_gv_series_with_provided_generators,
    compute_toric_curve_gv_diagnostics, compute_toric_two_face_curve_gv_invariants,
    curve_in_rational_row_span, curve_row_span_rank, curve_volume_in_divisor_basis,
    detect_apparent_nilpotent_ray_from_gv_multiples, detect_apparent_nilpotent_rays_from_gv_table,
    diagnose_affine_toric_circuit, diagnose_extremal_mori_ray_separator_by_lp_search,
    extract_ckyz_local_gv_invariants_from_potential, find_extremal_mori_ray_separator,
    find_extremal_mori_ray_separator_by_lp_search, find_pair_decomposition,
    find_semigroup_decomposition, finite_cutoff_gv_charges_excluding_primitive_rays,
    finite_gv_nonzero_degree_slice_points, gv_divisor_basis_data, identify_ckyz_local_surface,
    intersection_in_divisor_basis, intersection_in_matrix_divisor_basis,
    map_basis_gv_invariants_to_ambient, nilpotent_ray_degree_slice_for_cutoff_fraction,
    nilpotent_ray_divergence_check_from_slice_distances,
    nilpotent_ray_divergence_check_with_explicit_slice_lattices,
    nilpotent_ray_lll_reduced_slice_distance, nilpotent_ray_slice_comparison_points,
    partition_finite_cutoff_gv_charges_by_nilpotence, potent_ray_convergence,
    potent_ray_log_xi_slope, potent_ray_log_xi_terms, project_ambient_curve_to_basis,
    project_ambient_curve_to_basis_matrix, project_mori_cone_cap_rays_for_divisor_basis,
    project_mori_cone_cap_rays_to_basis, project_mori_cone_cap_rays_to_basis_matrix,
    prune_decomposable_curve_candidates, rank_two_local_charge_model,
    rank_two_local_support_signature, remove_pair_decomposable_curve_candidates,
    remove_semigroup_decomposable_curve_candidates, subcutoff_toric_curve_candidates,
    substitute_ckyz_series_in_flat_coordinates, supporting_mori_face_for_curve_from_normal,
    supporting_mori_face_from_normal,
};
pub use height_kahler::{
    effective_prime_divisors_from_curve_basis, heights_to_kahler, kahler_to_heights,
    project_heights_to_kahler,
};
pub use integer_math::{compute_glsm_and_linear_relations, compute_linear_relations_no_origin};
pub use intersection::{
    AmbientIntersectionNumbers, CytoolsNefPartitionCertificate,
    FullDimensionalNefPartitionCertificate, Intersection, certify_nef_partition_cytools_style,
    certify_nef_partition_cytools_style_full_dim,
    compute_complete_intersection_cy3_from_ambient_cytools,
    compute_complete_intersection_cy3_intersection_numbers, compute_intersection_numbers,
    compute_intersection_numbers_with_linear_relations, compute_intersection_numbers_with_offset,
};
pub use intersection::{compute_ambient_intersections_cytools, compute_intersection_cytools};
pub mod kklt;
pub mod kklt_vacuum;

pub use kklt::{
    CertifiedFlopContinuationStep, FlopContinuationState, GvDilogFailure, KkltBranchSearchResult,
    KkltBranchSolution, KkltJacobianDiagnostics, KkltResult, StableWeylCandidateCertificate,
    apply_certified_flop_continuation_sequence, check_stable_weyl_candidate_certificate,
    compute_c_tau as kklt_compute_c_tau, compute_corrected_target_tau,
    compute_divisor_volumes as kklt_compute_tau, compute_gv_target_correction_for_ambient_curves,
    compute_jacobian as kklt_compute_jacobian, compute_jacobian_diagnostics,
    compute_kklt_divisor_volumes, compute_kklt_divisor_volumes_in_divisor_basis,
    compute_kklt_jacobian, compute_kklt_jacobian_in_divisor_basis, compute_target_tau,
    divisor_quadratic_vanishes_on_curve_facet, find_stable_weyl_candidate_certificate,
    flop_reassign_gv_invariants, flop_transform_c2_vector, flop_transform_intersection_numbers,
    generate_scaled_divisor_basis_branch_initializations,
    generate_scaled_kklt_branch_initializations, gv_dilog_from_curve_volume_checked,
    scale_divisor_basis_kklt_branch_initialization_to_target,
    scale_mixed_basis_kklt_branch_initialization_to_target, solve_divisor_basis_path_following,
    solve_divisor_basis_path_following_branch_candidates, solve_mixed_basis_path_following,
    solve_mixed_basis_path_following_branch_candidates, solve_path_following,
    solve_two_phase_mixed_basis_path_following, solve_two_phase_path_following,
    transform_intersection_numbers_by_matrix, weyl_reflection_matches_flop_transform,
    weyl_reflection_matrix,
};
pub use lattice::Point;
pub use polytope::{Polytope, PolytopeFaces4d};
pub use racetrack::{
    GvInvariant, RacetrackResult, RacetrackTerm, build_racetrack_terms, compute_w0,
    compute_w0_from_terms, solve_racetrack,
};
pub use triangulation::{
    Triangulation, circuit_triangulation_choices, compute_delaunay_heights, compute_frst_heights,
    compute_regular_triangulation, expanded_secondary_chamber_choice_count,
    expanded_secondary_chamber_choice_digits,
    expanded_secondary_chamber_count_from_face_inequality_choices,
    expanded_secondary_chamber_hyperplanes_from_choice,
    expanded_secondary_chamber_hyperplanes_from_choice_index,
    expanded_secondary_chamber_hyperplanes_from_face_inequality_choice_index,
    expanded_secondary_chamber_hyperplanes_from_face_triangulation_choice_index,
    expanded_secondary_chamber_hyperplanes_from_face_triangulation_choice_index_with_star_4d,
    expanded_secondary_chamber_hyperplanes_on_polytope_2faces_choice_index_4d,
    expanded_secondary_chamber_hyperplanes_on_polytope_2faces_choice_index_4d_with_sampling,
    expanded_secondary_chamber_hyperplanes_on_polytope_2faces_choice_index_with_star_4d,
    expanded_secondary_chamber_hyperplanes_on_polytope_2faces_choice_index_with_star_4d_with_sampling,
    expanded_secondary_cone_hyperplanes_from_face_triangulations_with_star_4d,
    expanded_secondary_face_choice_indices_containing_height_vector,
    expanded_secondary_face_inequality_choices_from_circuit_faces,
    expanded_secondary_face_inequality_choices_from_triangulations,
    expanded_secondary_face_inequality_choices_from_triangulations_with_star_4d,
    expanded_secondary_face_inequality_choices_on_polytope_2faces_4d,
    expanded_secondary_face_inequality_choices_on_polytope_2faces_4d_with_sampling,
    expanded_secondary_face_inequality_choices_on_polytope_2faces_with_star_4d,
    expanded_secondary_face_inequality_choices_on_polytope_2faces_with_star_4d_with_sampling,
    expanded_secondary_fan_hyperplanes_on_faces,
    expanded_secondary_fan_hyperplanes_on_polytope_2faces_4d,
    expanded_secondary_group_boring_chamber_choices,
    expanded_secondary_regular_triangulation_from_face_triangulations,
    expanded_secondary_star_hyperplanes_for_face_triangulation_4d,
    fine_regular_triangulation_choices_on_polytope_2faces_4d,
    fine_regular_triangulation_choices_on_polytope_2faces_4d_with_sampling,
    fine_regular_triangulations_of_face_2d, grow_fine_triangulation_of_face_2d,
    sample_fine_regular_triangulations_of_face_2d, secondary_cone_height_pairings,
    secondary_cone_hyperplanes_native, secondary_cone_hyperplanes_native_on_faces,
    secondary_cone_hyperplanes_native_on_polytope_2faces_4d,
    secondary_cone_strictly_contains_height_vector, triangulation_heights_from_secondary_cone,
    triangulation_is_regular_from_secondary_cone,
};
pub use vacuum::{VacuumResult, compute_v0, compute_vacuum};
pub use volume::{VolumeResult, bbhl_correction, compute_volume, volume_classical, volume_string};
