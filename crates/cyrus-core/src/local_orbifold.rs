//! Local orbifold mirror-map diagnostics used by the corrected-chamber GV path.
//!
//! These helpers implement the CCIT `K_{\mathbb P(1,1,2)}` one-parameter
//! mirror map and the adjacent split-bundle diagnostics needed by the
//! McAllister source-history audit. They compute exact rational series
//! coefficients and do not assign compact CY3 GV values.

use std::collections::BTreeMap;

use malachite::{Integer, Rational};
use serde::Serialize;

/// First nonzero half-sector descendant readout for one half-integer degree.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HalfSectorFirstNonzeroDescendant {
    /// Inverse power of `z` where the first nonzero readout appears.
    pub inverse_z_power: Option<String>,
    /// Lambda-polynomial coefficients after dual-basis normalization.
    pub coefficients: Vec<String>,
}

/// Dual-basis-normalized half-sector readout at one inverse power of `z`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct HalfSectorDescendantReadout {
    /// Inverse power of `z` for this readout.
    pub inverse_z_power: String,
    /// Raw lambda-polynomial coefficients before dual-basis normalization.
    pub raw_lambda_coefficients: Vec<String>,
    /// Lambda-polynomial coefficients after division by `2*lambda`, when valid.
    pub dual_basis_coefficients: Option<Vec<String>>,
    /// Normalization status for this readout.
    pub status: String,
}

/// Ordinary-sector dual-basis `p^2` readout at one inverse power of `z`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct OrdinarySectorDualBasisP2ZReadout {
    /// Inverse power of `z` for this readout.
    pub inverse_z_power: String,
    /// Raw ordinary-basis `fun_0` lambda-polynomial coefficients before
    /// dual-basis normalization.
    pub raw_fun0_lambda_coefficients: Vec<String>,
    /// Lambda-polynomial coefficients after division by `2*lambda`, when valid.
    pub dual_basis_p2_coefficients: Option<Vec<String>>,
    /// Normalization status for this readout.
    pub status: String,
}

/// Chen-Ruan source-basis facts for the rank-three
/// `O(-1)+O(-1)+O(-2) -> P(1,1,2)` split-bundle diagnostic.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct WeightedP2RankThreeChenRuanSourceBasisReadout {
    /// Whether the source basis has been identified.
    pub source_status: String,
    /// Boundary between adjacent `K_P(1,1,2)` source data and the split bundle.
    pub source_scope_status: String,
    /// Chen-Ruan inertia-component labels used by the CCIT source.
    pub inertia_component_classes: Vec<String>,
    /// Untwisted hyperplane class on `P(1,1,2)`.
    pub untwisted_hyperplane_class: String,
    /// Stack-normalized square of the hyperplane class.
    pub base_hyperplane_square: String,
    /// Untwisted codimension-two insertion candidate.
    pub untwisted_codim2_candidate_class: String,
    /// Twisted half-sector insertion candidate.
    pub twisted_half_sector_candidate_class: String,
    /// Dual-basis class whose one-point correlator reads the `p^2` insertion.
    pub source_dual_basis_p2_class: String,
    /// Normalization divisor for the `p^2` dual-basis readout.
    pub source_dual_basis_p2_normalization_divisor: String,
    /// Status for reading `p^2` from the dual Chen-Ruan basis.
    pub source_dual_basis_p2_readout_status: String,
    /// Split-bundle twisted dual-basis class for the untwisted `p^2`
    /// insertion, using the inverse equivariant Euler twist.
    pub split_bundle_twisted_dual_basis_p2_class: String,
    /// Fun_0 coefficient divisor for the split-bundle twisted `p^2` dual.
    pub split_bundle_twisted_dual_basis_p2_normalization_divisor: String,
    /// Status for the split-bundle twisted `p^2` pairing normalization.
    pub split_bundle_twisted_dual_basis_p2_pairing_status: String,
    /// Dual-basis half-sector class.
    pub source_dual_basis_half_sector_class: String,
    /// Normalization boundary for stack and orbifold-sector pairings.
    pub stack_pairing_normalization_status: String,
    /// Status for using the adjacent canonical crepant-resolution map.
    pub adjacent_canonical_crepant_resolution_status: String,
    /// Image of `p^2` in the adjacent canonical resolution coordinates.
    pub adjacent_canonical_crepant_resolution_p2_image: String,
    /// Image of the half-sector class in adjacent canonical resolution coordinates.
    pub adjacent_canonical_crepant_resolution_half_sector_image: String,
    /// Label for the adjacent canonical `p^2` correlator.
    pub adjacent_canonical_correlator_p2_label: String,
    /// Label for the adjacent canonical half-sector correlator.
    pub adjacent_canonical_correlator_half_sector_label: String,
    /// Status for the adjacent canonical mirror-map normalization.
    pub adjacent_canonical_mirror_map_status: String,
    /// Source of the split-bundle I-function modification.
    pub split_bundle_ifunction_modification_source: String,
    /// Non-equivariant `p^2` readout status before mirror-map/pairing data.
    pub split_bundle_ordinary_non_equivariant_p2_readout_status: String,
    /// Status for inheriting the adjacent `K_P(1,1,2)` mirror map.
    pub split_bundle_mirror_map_status: String,
    /// Status for extracting split-bundle GW data.
    pub split_bundle_gw_extraction_status: String,
    /// Status for divisor-equation use in the twisted sector.
    pub twisted_sector_divisor_equation_status: String,
    /// Remaining source data needed before promotion.
    pub split_bundle_required_source_data_status: String,
    /// Final promotion boundary for this source readout.
    pub split_bundle_promotion_status: String,
}

/// Source-readiness boundary for a weighted-`P2` rank-three split bundle.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct WeightedP2RankThreeSplitBundleSourceReadiness {
    /// Complex dimension of the weighted-projective base.
    pub base_complex_dimension: i64,
    /// Sum of base weights, the anticanonical degree of the base.
    pub base_anticanonical_degree: i64,
    /// Sum of line-bundle degrees.
    pub bundle_degree_sum: i64,
    /// Total first Chern degree of the split-bundle total space.
    pub total_first_chern_degree: i64,
    /// Stack-normalized square of the base hyperplane class.
    pub base_hyperplane_square: String,
    /// Whether the base tensor is integral in the naive compact CY3 handoff.
    pub base_tensor_status: String,
    /// Rank of the direct-sum bundle.
    pub bundle_rank: i64,
    /// Complex dimension of the visible total-space phase.
    pub total_space_complex_dimension: i64,
    /// Codimension needed to reduce the visible phase to a CY3 invariant.
    pub required_cygv_codimension_for_threefold: i64,
    /// Whether the visible phase is already a numerical CY3 target.
    pub numerical_gv_status: String,
    /// Whether rank-two CKYZ local-surface formulas are applicable.
    pub ckyz_local_surface_source_status: String,
    /// Status for the twisted/vector-bundle I-function source route.
    pub twisted_vector_bundle_ifunction_source_status: String,
    /// Required insertion codimension for the twisted/vector-bundle route.
    pub twisted_vector_bundle_ifunction_required_insertion_complex_codimension: Option<i64>,
    /// Natural candidate insertion class, when the base dimension matches.
    pub twisted_vector_bundle_ifunction_candidate_insertion_class: Option<String>,
    /// Normalization status for the candidate insertion.
    pub twisted_vector_bundle_ifunction_insertion_normalization_status: Option<String>,
    /// Readiness status for promoting the twisted/vector-bundle source route.
    pub twisted_vector_bundle_ifunction_readiness_status: String,
    /// Explicit missing inputs before any numerical GV promotion is allowed.
    pub twisted_vector_bundle_ifunction_missing_inputs: Vec<String>,
}

/// Source-certificate requirements for promoting the weighted-`P2`
/// rank-three split-bundle I-function diagnostic to a numerical CY3
/// contribution.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct WeightedP2RankThreeTwistedIfunctionSourceCertificateRequirements {
    /// Overall theorem-handoff status.
    pub theorem_handoff_status: String,
    /// Source status for the compact/orbifold base small-J input.
    pub base_small_j_source_status: String,
    /// Source status for the split-bundle hypergeometric modification.
    pub split_bundle_modification_source_status: String,
    /// Source status for the Chen-Ruan basis used to read insertions.
    pub chen_ruan_source_basis_status: String,
    /// Reconstruction boundary for small-J/big-J and degree-two generation.
    pub chen_ruan_reconstruction_status: String,
    /// Status for the codimension-two observable needed by the visible CY5.
    pub required_observable_status: String,
    /// Status of the checked primary `z^-2` readout.
    pub primary_readout_status: String,
    /// Status of the first nonzero descendant/equivariant readout.
    pub descendant_readout_status: String,
    /// Minimum inverse `z` power in the checked positive-degree hypergeometric
    /// terms before multiplying by the leading outer `z`.
    pub j_scale_min_checked_hypergeometric_inverse_z_power: Option<i64>,
    /// Minimum inverse `z` power in the checked positive-degree terms after
    /// multiplying by the leading outer `z`.
    pub j_scale_min_checked_after_outer_z_inverse_power: Option<i64>,
    /// Whether checked positive-degree terms can modify the CCIT `F z + G`
    /// normalization data.
    pub j_scale_positive_degree_f_or_g_correction_status: String,
    /// Whether CCIT cone membership has been normalized to a twisted J-function.
    pub twisted_j_normalization_status: String,
    /// Whether the split-bundle mirror map needed for the twisted J-function is known.
    pub mirror_map_inversion_status: String,
    /// Whether the `z^-1` potential/invariant extraction step is available.
    pub z_minus_one_extraction_status: String,
    /// Whether the twisted dual pairing needed for coefficient extraction is known.
    pub twisted_dual_pairing_status: String,
    /// Pairing/residue normalization boundary.
    pub pairing_or_residue_status: String,
    /// Chamber-certificate boundary.
    pub chamber_certificate_status: String,
    /// qN-history boundary.
    pub qn_history_status: String,
    /// CY3 projection boundary.
    pub cy3_projection_status: String,
    /// Final promotion boundary.
    pub promotion_status: String,
    /// Maximum profiled degree, measured as twice the curve degree.
    pub checked_max_degree_twice: i64,
    /// Number of checked untwisted integer sectors.
    pub checked_integer_sector_count: usize,
    /// Number of checked half-degree twisted sectors.
    pub checked_half_sector_count: usize,
    /// Source facts used to justify the boundary.
    pub source_references: Vec<String>,
    /// Explicit missing inputs before any numerical promotion is allowed.
    pub required_inputs: Vec<String>,
}

/// Exact degree-by-degree CCIT I-function diagnostic for the weighted-`P2`
/// rank-three split bundle.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct WeightedP2RankThreeTwistedIfunctionDegreeProfile {
    /// Twice the curve degree. Odd values are half-degree twisted sectors.
    pub degree_twice: i64,
    /// Sector classification for this degree.
    pub sector_status: String,
    /// Weighted-base denominator factor counts.
    pub base_denominator_factor_counts: Vec<i64>,
    /// Split-bundle numerator factor counts.
    pub bundle_numerator_factor_counts: Vec<i64>,
    /// Number of non-equivariant zero factors in the split-bundle numerator.
    pub numerator_zero_factor_order: i64,
    /// Numerator factor count in the adjacent canonical `K_P(1,1,2)` model.
    pub adjacent_canonical_numerator_factor_count: i64,
    /// Zero-factor order in the adjacent canonical model.
    pub adjacent_canonical_numerator_zero_factor_order: i64,
    /// Comparison against the adjacent canonical source.
    pub split_vs_adjacent_canonical_factor_status: String,
    /// Naive non-equivariant insertion coefficient, when meaningful.
    pub split_non_equivariant_candidate_insertion_coefficient: Option<String>,
    /// Status for the naive non-equivariant insertion coefficient.
    pub split_non_equivariant_candidate_insertion_coefficient_status: String,
    /// Truncated non-equivariant numerator coefficients.
    pub split_non_equivariant_bundle_numerator_truncated_coefficients: Option<Vec<String>>,
    /// Status for the truncated non-equivariant numerator.
    pub split_non_equivariant_bundle_numerator_truncation_status: String,
    /// Equivariant candidate insertion lambda-polynomial.
    pub split_equivariant_candidate_insertion_lambda_polynomial: Option<Vec<String>>,
    /// First nonzero lambda order of the equivariant candidate insertion.
    pub split_equivariant_candidate_insertion_lambda_order: Option<i64>,
    /// Status for the equivariant candidate insertion.
    pub split_equivariant_candidate_insertion_status: String,
    /// Numerator-only dual-basis `p^2` lambda-polynomial.
    pub split_equivariant_dual_basis_p2_numerator_lambda_polynomial: Option<Vec<String>>,
    /// First nonzero lambda order of the numerator-only dual-basis readout.
    pub split_equivariant_dual_basis_p2_numerator_lambda_order: Option<i64>,
    /// Status for the numerator-only dual-basis readout.
    pub split_equivariant_dual_basis_p2_readout_status: String,
    /// Weighted-base denominator constant for the scalar hypergeometric readout.
    pub split_equivariant_dual_basis_p2_hypergeometric_denominator_constant: Option<String>,
    /// Scalar-denominator hypergeometric lambda-polynomial.
    pub split_equivariant_dual_basis_p2_hypergeometric_lambda_polynomial: Option<Vec<String>>,
    /// First nonzero lambda order after scalar denominator division.
    pub split_equivariant_dual_basis_p2_hypergeometric_lambda_order: Option<i64>,
    /// Status for the scalar-denominator hypergeometric readout.
    pub split_equivariant_dual_basis_p2_hypergeometric_status: String,
    /// Full weighted-base denominator truncated through the requested `p` power.
    pub split_equivariant_full_hypergeometric_denominator_truncated_coefficients:
        Option<Vec<String>>,
    /// Truncated inverse of the full weighted-base denominator.
    pub split_equivariant_full_hypergeometric_inverse_denominator_truncated_coefficients:
        Option<Vec<String>>,
    /// Full truncated quotient grouped by `p` power, then lambda power.
    pub split_equivariant_full_hypergeometric_truncated_coefficients: Option<Vec<Vec<String>>>,
    /// Full-denominator dual-basis `p^2` lambda-polynomial.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_lambda_polynomial: Option<Vec<String>>,
    /// First nonzero lambda order of the full-denominator dual-basis readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_lambda_order: Option<i64>,
    /// Status for the full-denominator dual-basis readout.
    pub split_equivariant_full_hypergeometric_status: String,
    /// Primary `z^-2` full-denominator dual-basis readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_primary_z2_lambda_polynomial:
        Option<Vec<String>>,
    /// First nonzero inverse `z` power in the full-denominator dual-basis readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_inverse_power:
        Option<i64>,
    /// First nonzero lambda-polynomial in the full-denominator `z` readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_polynomial:
        Option<Vec<String>>,
    /// First nonzero lambda order in the full-denominator `z` readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_order:
        Option<i64>,
    /// Non-equivariant-limit status for the first nonzero `z` readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_non_equivariant_limit_status:
        String,
    /// Status for the full-denominator `z` readout.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_z_readout_status: String,
    /// Whether scalar mirror-map reparametrization can promote a primary.
    pub split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status:
        String,
    /// Visibility status of the natural codimension-two insertion.
    pub candidate_insertion_visibility_status: String,
}

/// Compute the reusable source-readiness boundary for a weighted-`P2`
/// rank-three split bundle.
pub fn weighted_p2_rank_three_split_bundle_source_readiness(
    base_weights: &[i64],
    bundle_degrees: &[i64],
) -> Option<WeightedP2RankThreeSplitBundleSourceReadiness> {
    if base_weights.len() != 3 || bundle_degrees.len() != 3 {
        return None;
    }
    if base_weights.iter().any(|&weight| weight <= 0)
        || bundle_degrees.iter().any(|&degree| degree <= 0)
    {
        return None;
    }
    let base_complex_dimension =
        i64::try_from(base_weights.len() - 1).expect("weighted P2 base weight count fits i64");
    let base_anticanonical_degree = base_weights.iter().sum::<i64>();
    let bundle_degree_sum = bundle_degrees.iter().sum::<i64>();
    let total_first_chern_degree = base_anticanonical_degree - bundle_degree_sum;
    let base_weight_product = base_weights
        .iter()
        .try_fold(1i64, |product, &weight| product.checked_mul(weight))
        .expect("weighted P2 base weight product fits i64");
    let base_hyperplane_square = if base_weight_product == 1 {
        "1".to_string()
    } else {
        format!("1/{base_weight_product}")
    };
    let base_tensor_status = if base_weight_product == 1 {
        "weighted_p2_rank_three_base_hyperplane_square_integral".to_string()
    } else {
        "weighted_p2_rank_three_base_hyperplane_square_fractional_requires_stack_or_source_tensor_normalization"
            .to_string()
    };
    let bundle_rank =
        i64::try_from(bundle_degrees.len()).expect("weighted P2 bundle rank fits i64");
    let total_space_complex_dimension = base_complex_dimension + bundle_rank;
    let required_cygv_codimension_for_threefold = total_space_complex_dimension - 3;
    let numerical_gv_status = if total_space_complex_dimension == 3 {
        "weighted_p2_rank_three_visible_phase_is_threefold"
    } else {
        "weighted_p2_rank_three_visible_phase_is_not_numerical_cy3_requires_source_codim2_or_insertion_history"
    };
    let ckyz_local_surface_source_status = if total_space_complex_dimension == 3 && bundle_rank == 1
    {
        "weighted_p2_rank_three_ckyz_local_surface_source_maybe_applicable_requires_reflexive_polygon_identification"
    } else {
        "weighted_p2_rank_three_ckyz_local_surface_source_not_applicable_visible_phase_not_local_surface_cy3"
    };
    let twisted_vector_bundle_ifunction_source_status = if total_first_chern_degree == 0
        && total_space_complex_dimension != 3
    {
        "weighted_p2_rank_three_twisted_vector_bundle_ifunction_candidate_requires_insertions_and_qn_history"
    } else if total_first_chern_degree == 0 {
        "weighted_p2_rank_three_twisted_vector_bundle_ifunction_candidate_visible_threefold"
    } else {
        "weighted_p2_rank_three_twisted_vector_bundle_ifunction_not_calabi_yau_total_space"
    };
    let twisted_vector_bundle_ifunction_required_insertion_complex_codimension =
        if total_first_chern_degree == 0 && total_space_complex_dimension > 3 {
            Some(total_space_complex_dimension - 3)
        } else {
            None
        };
    let twisted_vector_bundle_ifunction_candidate_insertion_class =
        twisted_vector_bundle_ifunction_required_insertion_complex_codimension.and_then(
            |required_codimension| {
                (required_codimension == base_complex_dimension)
                    .then(|| format!("base_hyperplane_power_{required_codimension}"))
            },
        );
    let twisted_vector_bundle_ifunction_insertion_normalization_status =
        twisted_vector_bundle_ifunction_required_insertion_complex_codimension
            .map(|_| base_tensor_status.clone());
    let (
        twisted_vector_bundle_ifunction_readiness_status,
        twisted_vector_bundle_ifunction_missing_inputs,
    ) = if total_first_chern_degree != 0 {
        (
            "weighted_p2_rank_three_twisted_ifunction_not_calabi_yau_total_space".to_string(),
            vec!["calabi_yau_total_first_chern_zero".to_string()],
        )
    } else if total_space_complex_dimension == 3 {
        (
            "weighted_p2_rank_three_twisted_ifunction_visible_threefold_requires_qn_history"
                .to_string(),
            vec!["twisted_vector_bundle_ifunction_qn_history".to_string()],
        )
    } else if twisted_vector_bundle_ifunction_candidate_insertion_class.is_none() {
        (
            "weighted_p2_rank_three_twisted_ifunction_blocked_missing_source_insertion_class_qn_history".to_string(),
            vec![
                "source_derived_numerical_insertion_class".to_string(),
                "twisted_vector_bundle_ifunction_chamber_certificate".to_string(),
                "twisted_vector_bundle_ifunction_qn_history".to_string(),
            ],
        )
    } else if base_weight_product != 1 {
        (
            "weighted_p2_rank_three_twisted_ifunction_blocked_missing_stack_normalized_codim2_insertion_qn_history".to_string(),
            vec![
                "source_derived_codim2_insertion_or_equivalent_observable".to_string(),
                "stack_normalized_hyperplane_square_tensor".to_string(),
                "orbifold_sector_pairing_data".to_string(),
                "equivariant_residue_or_pairing_normalization".to_string(),
                "twisted_vector_bundle_ifunction_chamber_certificate".to_string(),
                "twisted_vector_bundle_ifunction_qn_history".to_string(),
            ],
        )
    } else {
        (
            "weighted_p2_rank_three_twisted_ifunction_blocked_missing_codim_insertion_qn_history"
                .to_string(),
            vec![
                "source_derived_codim_insertion_or_equivalent_observable".to_string(),
                "twisted_vector_bundle_ifunction_chamber_certificate".to_string(),
                "twisted_vector_bundle_ifunction_qn_history".to_string(),
            ],
        )
    };
    Some(WeightedP2RankThreeSplitBundleSourceReadiness {
        base_complex_dimension,
        base_anticanonical_degree,
        bundle_degree_sum,
        total_first_chern_degree,
        base_hyperplane_square,
        base_tensor_status,
        bundle_rank,
        total_space_complex_dimension,
        required_cygv_codimension_for_threefold,
        numerical_gv_status: numerical_gv_status.to_string(),
        ckyz_local_surface_source_status: ckyz_local_surface_source_status.to_string(),
        twisted_vector_bundle_ifunction_source_status:
            twisted_vector_bundle_ifunction_source_status.to_string(),
        twisted_vector_bundle_ifunction_required_insertion_complex_codimension,
        twisted_vector_bundle_ifunction_candidate_insertion_class,
        twisted_vector_bundle_ifunction_insertion_normalization_status,
        twisted_vector_bundle_ifunction_readiness_status,
        twisted_vector_bundle_ifunction_missing_inputs,
    })
}

/// Compute the exact source-certificate boundary for the weighted-`P2`
/// rank-three split-bundle I-function route.
pub fn weighted_p2_rank_three_twisted_ifunction_source_certificate_requirements(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    max_degree_twice: i64,
) -> Option<WeightedP2RankThreeTwistedIfunctionSourceCertificateRequirements> {
    let readiness =
        weighted_p2_rank_three_split_bundle_source_readiness(base_weights, bundle_degrees)?;
    let source_basis = weighted_p2_rank_three_split_bundle_chen_ruan_source_basis_readout(
        base_weights,
        bundle_degrees,
        readiness.total_first_chern_degree,
    )?;
    let profiles = weighted_p2_rank_three_twisted_ifunction_degree_profiles(
        base_weights,
        bundle_degrees,
        max_degree_twice,
        readiness.twisted_vector_bundle_ifunction_required_insertion_complex_codimension,
    );
    let integer_profiles = profiles
        .iter()
        .filter(|profile| profile.degree_twice % 2 == 0)
        .collect::<Vec<_>>();
    let checked_integer_sector_count = integer_profiles.len();
    let checked_half_sector_count = profiles.len() - checked_integer_sector_count;
    let all_integer_primary_readouts_zero = !integer_profiles.is_empty()
        && integer_profiles.iter().all(|profile| {
            profile
                .split_equivariant_full_hypergeometric_dual_basis_p2_primary_z2_lambda_polynomial
                .as_ref()
                .is_some_and(|polynomial| !lambda_polynomial_has_nonzero_coefficient(polynomial))
        });
    let any_descendant_pairing_blocker = integer_profiles.iter().any(|profile| {
        profile
            .split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_non_equivariant_limit_status
            .contains("requires")
    });
    let j_scale_min_checked_hypergeometric_inverse_z_power = profiles
        .iter()
        .map(checked_hypergeometric_inverse_z_power_bound)
        .min();
    let j_scale_min_checked_after_outer_z_inverse_power =
        j_scale_min_checked_hypergeometric_inverse_z_power
            .map(|inverse_z_power| inverse_z_power - 1);
    let j_scale_positive_degree_f_or_g_correction_status =
        match j_scale_min_checked_hypergeometric_inverse_z_power {
            None => {
                "weighted_p2_rank_three_ccit_j_scale_no_positive_degree_terms_checked".to_string()
            }
            Some(inverse_z_power) if inverse_z_power >= 2 => {
                "weighted_p2_rank_three_ccit_j_scale_checked_positive_degrees_cannot_modify_F_or_G_by_zero_order_bound"
                    .to_string()
            }
            Some(_) => {
                "weighted_p2_rank_three_ccit_j_scale_checked_positive_degrees_may_modify_F_or_G_requires_mirror_map"
                    .to_string()
            }
        };
    let primary_readout_status = if checked_integer_sector_count == 0 {
        "weighted_p2_rank_three_source_certificate_primary_readout_not_checked".to_string()
    } else if all_integer_primary_readouts_zero {
        "weighted_p2_rank_three_source_certificate_primary_z2_p2_readout_zero_to_checked_integer_degrees"
            .to_string()
    } else {
        "weighted_p2_rank_three_source_certificate_primary_z2_p2_readout_nonzero_requires_pairing_qn_history"
            .to_string()
    };
    let descendant_readout_status = if any_descendant_pairing_blocker {
        "weighted_p2_rank_three_source_certificate_first_nonzero_terms_are_descendant_or_equivariant_requires_big_j_pairing"
            .to_string()
    } else {
        "weighted_p2_rank_three_source_certificate_no_descendant_pairing_blocker_seen_in_checked_profile"
            .to_string()
    };
    let required_observable_status = if readiness.required_cygv_codimension_for_threefold
        == readiness.base_complex_dimension
    {
        "weighted_p2_rank_three_source_certificate_requires_base_p2_codim2_observable_or_equivalent"
    } else {
        "weighted_p2_rank_three_source_certificate_requires_source_derived_numerical_observable"
    };
    let mut required_inputs = readiness
        .twisted_vector_bundle_ifunction_missing_inputs
        .clone();
    required_inputs.retain(|input| {
        input != "stack_normalized_hyperplane_square_tensor"
            && input != "equivariant_residue_or_pairing_normalization"
    });
    for input in [
        "twisted_big_j_or_pairing_reconstruction_for_descendant_readout",
        "source_derived_chamber_qn_history_for_selected_phase",
        "cy3_projection_or_codimension_two_local_model",
    ] {
        if !required_inputs.iter().any(|existing| existing == input) {
            required_inputs.push(input.to_string());
        }
    }
    Some(WeightedP2RankThreeTwistedIfunctionSourceCertificateRequirements {
        theorem_handoff_status:
            "weighted_p2_rank_three_ccit_ifunction_lies_on_lagrangian_cone_not_yet_numerical_gv"
                .to_string(),
        base_small_j_source_status:
            "weighted_p2_rank_three_base_small_j_source_required_for_p112_orbifold".to_string(),
        split_bundle_modification_source_status: source_basis
            .split_bundle_ifunction_modification_source
            .clone(),
        chen_ruan_source_basis_status: source_basis.source_status.clone(),
        chen_ruan_reconstruction_status:
            "weighted_p2_rank_three_chen_ruan_degree_two_reconstruction_requires_big_j_restriction_or_pairing_input"
                .to_string(),
        required_observable_status: required_observable_status.to_string(),
        primary_readout_status,
        descendant_readout_status,
        j_scale_min_checked_hypergeometric_inverse_z_power,
        j_scale_min_checked_after_outer_z_inverse_power,
        j_scale_positive_degree_f_or_g_correction_status:
            j_scale_positive_degree_f_or_g_correction_status.clone(),
        twisted_j_normalization_status: if j_scale_positive_degree_f_or_g_correction_status
            .contains("cannot_modify_F_or_G")
        {
            "weighted_p2_rank_three_ccit_j_function_normalization_trivial_to_checked_positive_degrees"
        } else {
            "weighted_p2_rank_three_ccit_j_function_normalization_blocked_missing_split_bundle_Fz_plus_G_mirror_map"
        }
        .to_string(),
        mirror_map_inversion_status: if j_scale_positive_degree_f_or_g_correction_status
            .contains("cannot_modify_F_or_G")
        {
            "weighted_p2_rank_three_ccit_mirror_map_has_no_checked_positive_degree_F_or_G_correction"
        } else {
            "weighted_p2_rank_three_ccit_mirror_map_inversion_blocked_primary_signals_zero_and_descendant_pairing_missing"
        }
        .to_string(),
        z_minus_one_extraction_status: match j_scale_min_checked_hypergeometric_inverse_z_power {
            Some(inverse_z_power) if inverse_z_power > 2 => {
                "weighted_p2_rank_three_ccit_checked_positive_degrees_have_no_z_minus_one_primary_terms_first_nonzero_is_descendant_layer"
            }
            _ => {
                "weighted_p2_rank_three_ccit_z_minus_one_potential_extraction_blocked_missing_twisted_j_and_dual_pairing"
            }
        }
        .to_string(),
        twisted_dual_pairing_status:
            "weighted_p2_rank_three_ccit_untwisted_p2_dual_pairing_computed_but_big_j_pairing_reconstruction_missing"
                .to_string(),
        pairing_or_residue_status:
            "weighted_p2_rank_three_source_certificate_pairing_normalization_computed_for_untwisted_p2_but_descendant_reconstruction_missing"
                .to_string(),
        chamber_certificate_status:
            "weighted_p2_rank_three_source_certificate_blocked_missing_source_chamber_certificate"
                .to_string(),
        qn_history_status:
            "weighted_p2_rank_three_source_certificate_blocked_missing_source_qn_history"
                .to_string(),
        cy3_projection_status: readiness.numerical_gv_status,
        promotion_status:
            "weighted_p2_rank_three_source_certificate_not_promotable_without_pairing_chamber_qn_history"
                .to_string(),
        checked_max_degree_twice: max_degree_twice,
        checked_integer_sector_count,
        checked_half_sector_count,
        source_references: vec![
            "ccit_prop_small_j_reconstructs_lagrangian_cone_only_after_degree_two_or_big_j_input"
                .to_string(),
            "ccit_twisted_gw_total_space_uses_inverse_equivariant_euler_pairing".to_string(),
            "ccit_smalllinebundle_hypergeometric_modification_for_line_bundle".to_string(),
            "ccit_smallvb_direct_sum_modification_for_split_bundle".to_string(),
            "ccit_jscale_requires_Fz_plus_G_expansion_for_twisted_j_identification".to_string(),
            "ccit_local_invariants_extract_z_minus_one_term_using_twisted_dual_pairing"
                .to_string(),
        ],
        required_inputs,
    })
}

/// Compute exact degree profiles for the weighted-`P2` rank-three split-bundle
/// I-function diagnostic.
pub fn weighted_p2_rank_three_twisted_ifunction_degree_profiles(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    max_degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> Vec<WeightedP2RankThreeTwistedIfunctionDegreeProfile> {
    let adjacent_canonical_degree = base_weights.iter().sum::<i64>();
    (1..=max_degree_twice)
        .map(|degree_twice| {
            let base_denominator_factor_counts = base_weights
                .iter()
                .map(|&weight| ceil_positive_half_product(weight, degree_twice))
                .collect::<Vec<_>>();
            let bundle_numerator_factor_counts = bundle_degrees
                .iter()
                .map(|&degree| floor_positive_half_product(degree, degree_twice))
                .collect::<Vec<_>>();
            let numerator_zero_factor_order = bundle_degrees
                .iter()
                .filter(|&&degree| degree * degree_twice > 0 && (degree * degree_twice) % 2 == 0)
                .count() as i64;
            let adjacent_canonical_numerator_factor_count =
                floor_positive_half_product(adjacent_canonical_degree, degree_twice);
            let adjacent_canonical_numerator_zero_factor_order =
                i64::from(adjacent_canonical_degree * degree_twice > 0
                    && (adjacent_canonical_degree * degree_twice) % 2 == 0);
            let sector_status = if degree_twice % 2 == 0 {
                "weighted_p2_rank_three_twisted_ifunction_untwisted_integer_degree_sector"
            } else {
                "weighted_p2_rank_three_twisted_ifunction_half_degree_twisted_sector"
            }
            .to_string();
            let split_vs_adjacent_canonical_factor_status = match (
                degree_twice % 2 == 0,
                numerator_zero_factor_order.cmp(&adjacent_canonical_numerator_zero_factor_order),
            ) {
                (false, std::cmp::Ordering::Equal) => {
                    "weighted_p2_rank_three_split_matches_adjacent_canonical_zero_order_but_twisted_sector_needs_pairing"
                }
                (true, std::cmp::Ordering::Greater) => {
                    "weighted_p2_rank_three_split_integer_zero_order_exceeds_adjacent_canonical_kp112_so_p2_table_not_reusable"
                }
                (_, std::cmp::Ordering::Equal) => {
                    "weighted_p2_rank_three_split_zero_order_matches_adjacent_canonical"
                }
                (_, std::cmp::Ordering::Less) => {
                    "weighted_p2_rank_three_split_zero_order_below_adjacent_canonical_requires_source_check"
                }
                (_, std::cmp::Ordering::Greater) => {
                    "weighted_p2_rank_three_split_zero_order_exceeds_adjacent_canonical_requires_equivariant_normalization"
                }
            }
            .to_string();
            let (
                split_non_equivariant_candidate_insertion_coefficient,
                split_non_equivariant_candidate_insertion_coefficient_status,
            ) = split_non_equivariant_candidate_insertion_coefficient_status(
                degree_twice,
                numerator_zero_factor_order,
                required_insertion_complex_codimension,
            );
            let (
                split_non_equivariant_bundle_numerator_truncated_coefficients,
                split_non_equivariant_bundle_numerator_truncation_status,
            ) = split_non_equivariant_bundle_numerator_truncation_status(
                bundle_degrees,
                degree_twice,
                required_insertion_complex_codimension,
            );
            let (
                split_equivariant_candidate_insertion_lambda_polynomial,
                split_equivariant_candidate_insertion_lambda_order,
                split_equivariant_candidate_insertion_status,
            ) = split_equivariant_candidate_insertion_status(
                bundle_degrees,
                degree_twice,
                required_insertion_complex_codimension,
            );
            let (
                split_equivariant_dual_basis_p2_numerator_lambda_polynomial,
                split_equivariant_dual_basis_p2_numerator_lambda_order,
                split_equivariant_dual_basis_p2_readout_status,
            ) = split_equivariant_dual_basis_p2_readout_status(
                bundle_degrees,
                degree_twice,
                required_insertion_complex_codimension,
            );
            let (
                split_equivariant_dual_basis_p2_hypergeometric_denominator_constant,
                split_equivariant_dual_basis_p2_hypergeometric_lambda_polynomial,
                split_equivariant_dual_basis_p2_hypergeometric_lambda_order,
                split_equivariant_dual_basis_p2_hypergeometric_status,
            ) = split_equivariant_dual_basis_p2_hypergeometric_status(
                base_weights,
                bundle_degrees,
                degree_twice,
                required_insertion_complex_codimension,
            );
            let (
                split_equivariant_full_hypergeometric_denominator_truncated_coefficients,
                split_equivariant_full_hypergeometric_inverse_denominator_truncated_coefficients,
                split_equivariant_full_hypergeometric_truncated_coefficients,
                split_equivariant_full_hypergeometric_dual_basis_p2_lambda_polynomial,
                split_equivariant_full_hypergeometric_dual_basis_p2_lambda_order,
                split_equivariant_full_hypergeometric_status,
            ) = split_equivariant_full_hypergeometric_status(
                base_weights,
                bundle_degrees,
                degree_twice,
                required_insertion_complex_codimension,
                required_insertion_complex_codimension
                    .and_then(|codim| usize::try_from(codim).ok())
                    .unwrap_or(2),
            );
            let (
                split_equivariant_full_hypergeometric_dual_basis_p2_primary_z2_lambda_polynomial,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_inverse_power,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_polynomial,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_order,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_non_equivariant_limit_status,
                split_equivariant_full_hypergeometric_dual_basis_p2_z_readout_status,
            ) = split_equivariant_full_hypergeometric_dual_basis_p2_z_readout_status(
                base_weights,
                bundle_degrees,
                degree_twice,
                required_insertion_complex_codimension,
            );
            let split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status =
                split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status(
                    base_weights,
                    bundle_degrees,
                    degree_twice,
                    required_insertion_complex_codimension,
                );
            let candidate_insertion_visibility_status = match (
                degree_twice % 2 == 0,
                required_insertion_complex_codimension,
            ) {
                (false, _) => {
                    "weighted_p2_rank_three_twisted_ifunction_half_degree_requires_orbifold_sector_pairing"
                        .to_string()
                }
                (true, Some(required_codim))
                    if numerator_zero_factor_order > required_codim =>
                {
                    "weighted_p2_rank_three_twisted_ifunction_untwisted_zero_order_exceeds_required_insertion_codim_requires_equivariant_normalization"
                        .to_string()
                }
                (true, Some(required_codim))
                    if numerator_zero_factor_order == required_codim =>
                {
                    "weighted_p2_rank_three_twisted_ifunction_untwisted_zero_order_matches_required_insertion_codim"
                        .to_string()
                }
                (true, Some(_)) => {
                    "weighted_p2_rank_three_twisted_ifunction_untwisted_zero_order_below_required_insertion_codim_requires_higher_series_terms"
                        .to_string()
                }
                (true, None) => {
                    "weighted_p2_rank_three_twisted_ifunction_untwisted_insertion_codim_not_applicable"
                        .to_string()
                }
            };
            WeightedP2RankThreeTwistedIfunctionDegreeProfile {
                degree_twice,
                sector_status,
                base_denominator_factor_counts,
                bundle_numerator_factor_counts,
                numerator_zero_factor_order,
                adjacent_canonical_numerator_factor_count,
                adjacent_canonical_numerator_zero_factor_order,
                split_vs_adjacent_canonical_factor_status,
                split_non_equivariant_candidate_insertion_coefficient,
                split_non_equivariant_candidate_insertion_coefficient_status,
                split_non_equivariant_bundle_numerator_truncated_coefficients,
                split_non_equivariant_bundle_numerator_truncation_status,
                split_equivariant_candidate_insertion_lambda_polynomial,
                split_equivariant_candidate_insertion_lambda_order,
                split_equivariant_candidate_insertion_status,
                split_equivariant_dual_basis_p2_numerator_lambda_polynomial,
                split_equivariant_dual_basis_p2_numerator_lambda_order,
                split_equivariant_dual_basis_p2_readout_status,
                split_equivariant_dual_basis_p2_hypergeometric_denominator_constant,
                split_equivariant_dual_basis_p2_hypergeometric_lambda_polynomial,
                split_equivariant_dual_basis_p2_hypergeometric_lambda_order,
                split_equivariant_dual_basis_p2_hypergeometric_status,
                split_equivariant_full_hypergeometric_denominator_truncated_coefficients,
                split_equivariant_full_hypergeometric_inverse_denominator_truncated_coefficients,
                split_equivariant_full_hypergeometric_truncated_coefficients,
                split_equivariant_full_hypergeometric_dual_basis_p2_lambda_polynomial,
                split_equivariant_full_hypergeometric_dual_basis_p2_lambda_order,
                split_equivariant_full_hypergeometric_status,
                split_equivariant_full_hypergeometric_dual_basis_p2_primary_z2_lambda_polynomial,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_inverse_power,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_polynomial,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_order,
                split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_non_equivariant_limit_status,
                split_equivariant_full_hypergeometric_dual_basis_p2_z_readout_status,
                split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status,
                candidate_insertion_visibility_status,
            }
        })
        .collect()
}

fn checked_hypergeometric_inverse_z_power_bound(
    profile: &WeightedP2RankThreeTwistedIfunctionDegreeProfile,
) -> i64 {
    let denominator_z_power = profile.base_denominator_factor_counts.iter().sum::<i64>();
    let numerator_factor_count = profile.bundle_numerator_factor_counts.iter().sum::<i64>();
    denominator_z_power - (numerator_factor_count - profile.numerator_zero_factor_order)
}

/// Source-basis facts for the weighted-`P2` rank-three split bundle.
pub fn weighted_p2_rank_three_split_bundle_chen_ruan_source_basis_readout(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    total_first_chern_degree: i64,
) -> Option<WeightedP2RankThreeChenRuanSourceBasisReadout> {
    if base_weights != [1, 1, 2] || bundle_degrees != [1, 1, 2] {
        return None;
    }
    if total_first_chern_degree != 0 {
        return None;
    }
    Some(WeightedP2RankThreeChenRuanSourceBasisReadout {
        source_status: "weighted_p2_rank_three_chern_ruan_source_basis_for_p112_identified"
            .to_string(),
        source_scope_status:
            "weighted_p2_rank_three_chern_ruan_source_is_adjacent_kp112_not_split_bundle_qn_history"
                .to_string(),
        inertia_component_classes: vec!["fun_0".to_string(), "fun_{1/2}".to_string()],
        untwisted_hyperplane_class: "p=c1(O(1))".to_string(),
        base_hyperplane_square: "1/2".to_string(),
        untwisted_codim2_candidate_class: "p^2".to_string(),
        twisted_half_sector_candidate_class: "fun_{1/2}".to_string(),
        source_dual_basis_p2_class: "2*lambda*fun_0-8*p".to_string(),
        source_dual_basis_p2_normalization_divisor: "2*lambda".to_string(),
        source_dual_basis_p2_readout_status:
            "p2_one_point_correlator_is_read_from_dual_basis_phi2_not_from_raw_p2_coefficient"
                .to_string(),
        split_bundle_twisted_dual_basis_p2_class:
            "2*(lambda-p)^2*(lambda-2*p)*fun_0=2*lambda^3*fun_0-8*lambda^2*p+10*lambda*p^2"
                .to_string(),
        split_bundle_twisted_dual_basis_p2_normalization_divisor: "2*lambda^3".to_string(),
        split_bundle_twisted_dual_basis_p2_pairing_status:
            "weighted_p2_rank_three_split_bundle_untwisted_p2_dual_pairing_computed_from_ccit_inverse_euler_twist"
                .to_string(),
        source_dual_basis_half_sector_class: "2*lambda*fun_{1/2}".to_string(),
        stack_pairing_normalization_status:
            "weighted_p2_hyperplane_square_fractional_requires_stack_normalized_source_pairing_and_twisted_sector_pairing"
                .to_string(),
        adjacent_canonical_crepant_resolution_status:
            "weighted_p2_rank_three_adjacent_kp112_crepant_map_reference_only".to_string(),
        adjacent_canonical_crepant_resolution_p2_image: "p1*p2/2".to_string(),
        adjacent_canonical_crepant_resolution_half_sector_image: "i*(p1-2*p2)/2".to_string(),
        adjacent_canonical_correlator_p2_label: "<p^2>_{0,1,d}^{K_P(1,1,2)}".to_string(),
        adjacent_canonical_correlator_half_sector_label:
            "<fun_{1/2}>_{0,1,d}^{K_P(1,1,2)}".to_string(),
        adjacent_canonical_mirror_map_status:
            "kp112_canonical_source_uses_q_equals_x_exp_4h_after_o_minus_4_hypergeometric_modification"
                .to_string(),
        split_bundle_ifunction_modification_source:
            "ccit_smalllinebundle_smallvb_direct_sum_modification_product_over_o_minus_1_o_minus_1_o_minus_2"
                .to_string(),
        split_bundle_ordinary_non_equivariant_p2_readout_status:
            "rank_three_split_ordinary_non_equivariant_p2_readout_vanishes_before_mirror_map_because_codim2_truncation_of_bundle_numerator_is_zero"
                .to_string(),
        split_bundle_mirror_map_status:
            "rank_three_split_integer_terms_have_zero_order_three_so_adjacent_kp112_divisor_mirror_map_is_not_inherited"
                .to_string(),
        split_bundle_gw_extraction_status:
            "rank_three_split_requires_own_j_or_big_j_coefficient_extraction_before_qn_history"
                .to_string(),
        twisted_sector_divisor_equation_status:
            "chen_ruan_twisted_sector_has_no_divisor_equation_requires_big_j_or_pairing_input"
                .to_string(),
        split_bundle_required_source_data_status:
            "rank_three_split_nonzero_contribution_if_any_requires_equivariant_residue_or_twisted_big_j_pairing_not_canonical_kp112_table"
                .to_string(),
        split_bundle_promotion_status:
            "weighted_p2_rank_three_split_bundle_still_requires_twisted_vector_bundle_ifunction_normalization_and_qn_history"
                .to_string(),
    })
}

/// CCIT `K_P(1,1,2)` untwisted `b_d` table through `max_degree`.
pub fn kp112_canonical_b_table_values(max_degree: usize) -> Vec<String> {
    kp112_ordinary_p_primary_lambda_coefficients(&[4], max_degree)
        .iter()
        .map(|degree_primary| {
            let ordinary_p_lambda_coefficient = degree_primary
                .get(&1)
                .cloned()
                .unwrap_or_else(|| "0".to_string())
                .parse::<Rational>()
                .expect("ordinary p coefficient is rational");
            (-ordinary_p_lambda_coefficient / rational_from_i64(16)).to_string()
        })
        .collect()
}

/// Ordinary-`p` primary lambda-polynomial coefficients after the adjacent
/// `K_P(1,1,2)` mirror map for a direct-sum line-bundle modification.
pub fn kp112_ordinary_p_primary_lambda_coefficients(
    bundle_degrees: &[i64],
    max_degree: usize,
) -> Vec<BTreeMap<usize, String>> {
    weighted_p2_ordinary_p_primary_after_kp112_mirror_map(bundle_degrees, max_degree)
        .iter()
        .take(max_degree + 1)
        .skip(1)
        .map(rational_lambda_polynomial_to_strings)
        .collect()
}

/// Split-bundle `[1,1,2]` ordinary-`p` primary coefficients after the
/// adjacent `K_P(1,1,2)` mirror map.
pub fn split_bundle_kp112_mirror_map_primary_p_lambda_coefficients(
    max_degree: usize,
) -> Vec<String> {
    kp112_ordinary_p_primary_lambda_coefficients(&[1, 1, 2], max_degree)
        .into_iter()
        .map(|degree_primary| {
            degree_primary
                .get(&1)
                .cloned()
                .unwrap_or_else(|| "0".to_string())
        })
        .collect()
}

/// Ordinary-sector dual-basis `p^2` readout by inverse `z` power for a
/// weighted-`P2` direct-sum line-bundle modification.
pub fn weighted_p2_ordinary_dual_basis_p2_z_readout_profile(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
) -> Option<Vec<OrdinarySectorDualBasisP2ZReadout>> {
    let raw_by_z = weighted_p2_ordinary_fun0_lambda_polynomials_by_z(
        base_weights,
        bundle_degrees,
        degree_twice,
    )?;
    let max_inverse_z_power = raw_by_z.keys().last().copied().unwrap_or(2).max(2);
    Some(
        (2..=max_inverse_z_power)
            .map(|inverse_z_power| {
                let polynomial = raw_by_z.get(&inverse_z_power).cloned().unwrap_or_default();
                ordinary_sector_dual_basis_p2_z_readout(
                    inverse_z_power,
                    &polynomial,
                    bundle_degrees.len(),
                )
            })
            .collect(),
    )
}

/// Ordinary-sector dual-basis `p^2` lambda-polynomial coefficients by inverse
/// `z` power for a weighted-`P2` direct-sum line-bundle modification.
pub fn weighted_p2_ordinary_dual_basis_p2_z_lambda_polynomials(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
) -> Option<BTreeMap<i64, Vec<String>>> {
    let raw_by_z = weighted_p2_ordinary_fun0_lambda_polynomials_by_z(
        base_weights,
        bundle_degrees,
        degree_twice,
    )?;
    let mut normalized_by_z = BTreeMap::new();
    for (inverse_z_power, polynomial) in raw_by_z {
        let coefficients = divide_rational_lambda_polynomial_by_two_lambda_power_if_possible(
            &polynomial,
            bundle_degrees.len(),
        )
        .map(|coefficients| {
            let coefficients = rational_coefficients_to_strings(&coefficients);
            if coefficients.is_empty() {
                vec!["0".to_string()]
            } else {
                coefficients
            }
        })?;
        normalized_by_z.insert(inverse_z_power, coefficients);
    }
    Some(normalized_by_z)
}

/// CCIT `K_P(1,1,2)` twisted half-sector `c_d` table.
pub fn kp112_canonical_c_table_values(half_sector_count: usize) -> Vec<String> {
    kp112_half_sector_primary_lambda_coefficients(&[4], half_sector_count)
        .iter()
        .map(|half_sector_primary| {
            let fun_half_lambda_coefficient = half_sector_primary
                .get(&1)
                .cloned()
                .unwrap_or_else(|| "0".to_string())
                .parse::<Rational>()
                .expect("half-sector coefficient is rational");
            (fun_half_lambda_coefficient / rational_from_i64(2)).to_string()
        })
        .collect()
}

/// Half-sector primary lambda-polynomial coefficients after the adjacent
/// `K_P(1,1,2)` mirror map for a direct-sum line-bundle modification.
pub fn kp112_half_sector_primary_lambda_coefficients(
    bundle_degrees: &[i64],
    half_sector_count: usize,
) -> Vec<BTreeMap<usize, String>> {
    weighted_p2_half_sector_primary_after_kp112_mirror_map(bundle_degrees, half_sector_count)
        .iter()
        .map(rational_lambda_polynomial_to_strings)
        .collect()
}

/// Split-bundle `[1,1,2]` half-sector primary coefficients after the adjacent
/// `K_P(1,1,2)` mirror map.
pub fn split_bundle_kp112_mirror_map_half_sector_primary_coefficients(
    half_sector_count: usize,
) -> Vec<String> {
    kp112_half_sector_primary_lambda_coefficients(&[1, 1, 2], half_sector_count)
        .into_iter()
        .map(|half_sector_primary| {
            let fun_half_lambda_coefficient = half_sector_primary
                .get(&1)
                .cloned()
                .unwrap_or_else(|| "0".to_string())
                .parse::<Rational>()
                .expect("half-sector coefficient is rational");
            (fun_half_lambda_coefficient / rational_from_i64(2)).to_string()
        })
        .collect()
}

/// Split-bundle `[1,1,2]` first nonzero half-sector descendant readouts after
/// the adjacent `K_P(1,1,2)` mirror map.
pub fn split_bundle_kp112_mirror_map_half_sector_first_nonzero_descendants(
    half_sector_count: usize,
    max_inverse_z_power: i64,
) -> Vec<HalfSectorFirstNonzeroDescendant> {
    kp112_half_sector_first_nonzero_descendants(&[1, 1, 2], half_sector_count, max_inverse_z_power)
}

/// Full split-bundle `[1,1,2]` half-sector descendant profile after the
/// adjacent `K_P(1,1,2)` mirror map, grouped by half-integer degree.
pub fn split_bundle_kp112_mirror_map_half_sector_descendant_profiles(
    half_sector_count: usize,
    max_inverse_z_power: i64,
) -> Vec<Vec<HalfSectorDescendantReadout>> {
    kp112_half_sector_dual_basis_descendant_profiles(
        &[1, 1, 2],
        half_sector_count,
        max_inverse_z_power,
    )
}

/// Full half-sector descendant profiles after the adjacent `K_P(1,1,2)`
/// mirror map for a direct-sum line-bundle modification, grouped by
/// half-integer degree.
pub fn kp112_half_sector_dual_basis_descendant_profiles(
    bundle_degrees: &[i64],
    half_sector_count: usize,
    max_inverse_z_power: i64,
) -> Vec<Vec<HalfSectorDescendantReadout>> {
    weighted_p2_half_sector_fun_component_after_kp112_mirror_map_by_z(
        bundle_degrees,
        half_sector_count,
        max_inverse_z_power,
    )
    .into_iter()
    .map(|by_z| {
        (2..=max_inverse_z_power)
            .map(|inverse_z_power| {
                let polynomial = by_z.get(&inverse_z_power).cloned().unwrap_or_default();
                half_sector_descendant_readout(inverse_z_power, &polynomial)
            })
            .collect()
    })
    .collect()
}

/// First nonzero half-sector descendant readouts after the adjacent
/// `K_P(1,1,2)` mirror map for a direct-sum line-bundle modification.
pub fn kp112_half_sector_first_nonzero_descendants(
    bundle_degrees: &[i64],
    half_sector_count: usize,
    max_inverse_z_power: i64,
) -> Vec<HalfSectorFirstNonzeroDescendant> {
    weighted_p2_half_sector_fun_component_after_kp112_mirror_map_by_z(
        bundle_degrees,
        half_sector_count,
        max_inverse_z_power,
    )
    .into_iter()
    .map(|by_z| {
        let first_nonzero = by_z
            .into_iter()
            .filter_map(|(inverse_z_power, polynomial)| {
                divide_rational_lambda_polynomial_by_two_lambda_if_possible(&polynomial).map(
                    |coefficients| {
                        (
                            inverse_z_power,
                            rational_coefficients_to_strings(&coefficients),
                        )
                    },
                )
            })
            .find(|(_, coefficients)| lambda_polynomial_has_nonzero_coefficient(coefficients));
        match first_nonzero {
            Some((inverse_z_power, coefficients)) => HalfSectorFirstNonzeroDescendant {
                inverse_z_power: Some(inverse_z_power.to_string()),
                coefficients,
            },
            None => HalfSectorFirstNonzeroDescendant {
                inverse_z_power: None,
                coefficients: vec!["0".to_string()],
            },
        }
    })
    .collect()
}

fn half_sector_descendant_readout(
    inverse_z_power: i64,
    polynomial: &BTreeMap<usize, Rational>,
) -> HalfSectorDescendantReadout {
    let raw_lambda_coefficients = dense_rational_lambda_polynomial_to_strings(polynomial);
    let dual_basis_coefficients = divide_rational_lambda_polynomial_by_two_lambda_if_possible(
        polynomial,
    )
    .map(|coefficients| {
        let coefficients = rational_coefficients_to_strings(&coefficients);
        if coefficients.is_empty() {
            vec!["0".to_string()]
        } else {
            coefficients
        }
    });
    let status = match &dual_basis_coefficients {
        Some(coefficients) if lambda_polynomial_has_nonzero_coefficient(coefficients) => {
            "dual_basis_normalized_nonzero"
        }
        Some(_) => "dual_basis_normalized_zero",
        None => "blocked_constant_lambda_term_nonzero",
    };
    HalfSectorDescendantReadout {
        inverse_z_power: inverse_z_power.to_string(),
        raw_lambda_coefficients,
        dual_basis_coefficients,
        status: status.to_string(),
    }
}

fn ordinary_sector_dual_basis_p2_z_readout(
    inverse_z_power: i64,
    polynomial: &BTreeMap<usize, Rational>,
    lambda_power: usize,
) -> OrdinarySectorDualBasisP2ZReadout {
    let raw_fun0_lambda_coefficients = dense_rational_lambda_polynomial_to_strings(polynomial);
    let dual_basis_p2_coefficients =
        divide_rational_lambda_polynomial_by_two_lambda_power_if_possible(polynomial, lambda_power)
            .map(|coefficients| {
                let coefficients = rational_coefficients_to_strings(&coefficients);
                if coefficients.is_empty() {
                    vec!["0".to_string()]
                } else {
                    coefficients
                }
            });
    let status = match &dual_basis_p2_coefficients {
        Some(coefficients) if lambda_polynomial_has_nonzero_coefficient(coefficients) => {
            "dual_basis_p2_normalized_nonzero"
        }
        Some(_) => "dual_basis_p2_normalized_zero",
        None => "blocked_constant_lambda_term_nonzero",
    };
    OrdinarySectorDualBasisP2ZReadout {
        inverse_z_power: inverse_z_power.to_string(),
        raw_fun0_lambda_coefficients,
        dual_basis_p2_coefficients,
        status: status.to_string(),
    }
}

fn rational_lambda_polynomial_to_strings(
    polynomial: &BTreeMap<usize, Rational>,
) -> BTreeMap<usize, String> {
    polynomial
        .iter()
        .map(|(&lambda_power, coefficient)| (lambda_power, coefficient.to_string()))
        .collect()
}

fn dense_rational_lambda_polynomial_to_strings(
    polynomial: &BTreeMap<usize, Rational>,
) -> Vec<String> {
    let max_lambda_power = polynomial.keys().last().copied().unwrap_or(0);
    (0..=max_lambda_power)
        .map(|lambda_power| {
            polynomial
                .get(&lambda_power)
                .cloned()
                .unwrap_or_else(|| Rational::from(0))
                .to_string()
        })
        .collect()
}

fn rational_from_i64(value: i64) -> Rational {
    Rational::from(Integer::from(value))
}

fn factorial_integer(value: usize) -> Integer {
    (1..=value).fold(Integer::from(1), |product, factor| {
        product * Integer::from(factor)
    })
}

fn truncated_series_multiply(
    lhs: &[Rational],
    rhs: &[Rational],
    max_degree: usize,
) -> Vec<Rational> {
    let mut product = vec![Rational::from(0); max_degree + 1];
    for (lhs_degree, lhs_coefficient) in lhs.iter().enumerate() {
        if lhs_coefficient == &Rational::from(0) {
            continue;
        }
        for (rhs_degree, rhs_coefficient) in rhs.iter().enumerate() {
            let total_degree = lhs_degree + rhs_degree;
            if total_degree > max_degree || rhs_coefficient == &Rational::from(0) {
                continue;
            }
            product[total_degree] += lhs_coefficient.clone() * rhs_coefficient.clone();
        }
    }
    product
}

fn truncated_series_compose(
    series: &[Rational],
    argument: &[Rational],
    max_degree: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::from(0); max_degree + 1];
    let mut argument_power = vec![Rational::from(0); max_degree + 1];
    argument_power[0] = Rational::from(1);
    for (degree, coefficient) in series.iter().enumerate().take(max_degree + 1) {
        if degree > 0 {
            argument_power = truncated_series_multiply(&argument_power, argument, max_degree);
        }
        if coefficient == &Rational::from(0) {
            continue;
        }
        for output_degree in 0..=max_degree {
            result[output_degree] += coefficient.clone() * argument_power[output_degree].clone();
        }
    }
    result
}

fn truncated_series_exp(argument: &[Rational], max_degree: usize) -> Vec<Rational> {
    let mut result = vec![Rational::from(0); max_degree + 1];
    result[0] = Rational::from(1);
    let mut argument_power = vec![Rational::from(0); max_degree + 1];
    argument_power[0] = Rational::from(1);
    for power in 1..=max_degree {
        argument_power = truncated_series_multiply(&argument_power, argument, max_degree);
        let denominator = Rational::from(factorial_integer(power));
        for degree in 0..=max_degree {
            result[degree] += argument_power[degree].clone() / denominator.clone();
        }
    }
    result
}

fn truncated_series_power(
    series: &[Rational],
    exponent: usize,
    max_degree: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::from(0); max_degree + 1];
    result[0] = Rational::from(1);
    for _ in 0..exponent {
        result = truncated_series_multiply(&result, series, max_degree);
    }
    result
}

fn kp112_h_series(max_degree: usize) -> Vec<Rational> {
    let mut series = vec![Rational::from(0); max_degree + 1];
    for (degree, coefficient) in series.iter_mut().enumerate().take(max_degree + 1).skip(1) {
        let numerator = factorial_integer(4 * degree - 1);
        let denominator =
            factorial_integer(degree) * factorial_integer(degree) * factorial_integer(2 * degree);
        *coefficient = Rational::from(numerator) / Rational::from(denominator);
    }
    series
}

fn kp112_x_of_q_series(max_degree: usize) -> Vec<Rational> {
    let h_series = kp112_h_series(max_degree);
    let mut q_series = vec![Rational::from(0); max_degree + 1];
    if max_degree >= 1 {
        q_series[1] = Rational::from(1);
    }
    let mut x_series = q_series.clone();
    for _ in 0..max_degree {
        let h_of_x = truncated_series_compose(&h_series, &x_series, max_degree);
        let argument = h_of_x
            .into_iter()
            .map(|coefficient| -rational_from_i64(4) * coefficient)
            .collect::<Vec<_>>();
        let exponential = truncated_series_exp(&argument, max_degree);
        x_series = truncated_series_multiply(&q_series, &exponential, max_degree);
    }
    x_series
}

fn weighted_p2_ordinary_p_primary_after_kp112_mirror_map(
    bundle_degrees: &[i64],
    max_degree: usize,
) -> Vec<BTreeMap<usize, Rational>> {
    let x_of_q = kp112_x_of_q_series(max_degree);
    let h_of_q = truncated_series_compose(&kp112_h_series(max_degree), &x_of_q, max_degree);
    let mut raw_by_inverse_z_power = BTreeMap::<i64, BTreeMap<usize, Vec<Rational>>>::new();
    for degree in 1..=max_degree {
        let degree_twice = i64::try_from(2 * degree).expect("degree fits i64");
        let denominator = integer_sector_base_denominator_constant(&[1, 1, 2], degree_twice)
            .expect("integer weighted P2 denominator exists");
        let denominator_z_power = i64::try_from(4 * degree).expect("degree fits i64");
        let mut numerator = BTreeMap::<(usize, usize, i64), Rational>::new();
        numerator.insert((0, 0, 0), Rational::from(1));
        for &line_degree in bundle_degrees {
            let factor_count = line_degree
                .checked_mul(i64::try_from(degree).expect("degree fits i64"))
                .expect("factor count fits i64");
            for factor_offset in (-factor_count + 1)..=0 {
                let mut next = BTreeMap::<(usize, usize, i64), Rational>::new();
                for (&(p_power, lambda_power, z_power), coefficient) in &numerator {
                    *next
                        .entry((p_power, lambda_power + 1, z_power))
                        .or_insert_with(|| Rational::from(0)) += coefficient.clone();
                    if p_power < 1 {
                        *next
                            .entry((p_power + 1, lambda_power, z_power))
                            .or_insert_with(|| Rational::from(0)) +=
                            -rational_from_i64(line_degree) * coefficient.clone();
                    }
                    if factor_offset != 0 {
                        *next
                            .entry((p_power, lambda_power, z_power + 1))
                            .or_insert_with(|| Rational::from(0)) +=
                            rational_from_i64(factor_offset) * coefficient.clone();
                    }
                }
                next.retain(|_, coefficient| *coefficient != Rational::from(0));
                numerator = next;
            }
        }

        let x_power = truncated_series_power(&x_of_q, degree, max_degree);
        for ((p_power, lambda_power, z_power), coefficient) in numerator {
            if p_power != 1 {
                continue;
            }
            let inverse_z_power = denominator_z_power - z_power;
            let q_series = raw_by_inverse_z_power
                .entry(inverse_z_power)
                .or_default()
                .entry(lambda_power)
                .or_insert_with(|| vec![Rational::from(0); max_degree + 1]);
            let scaled = coefficient / denominator.clone();
            for q_degree in 0..=max_degree {
                q_series[q_degree] += scaled.clone() * x_power[q_degree].clone();
            }
        }
    }

    let exp_terms = kp112_lambda_expansion_terms(&h_of_q, max_degree, 2);
    let mut primary_by_q_degree = vec![BTreeMap::<usize, Rational>::new(); max_degree + 1];
    for (inverse_z_power, lambda_series_by_power) in raw_by_inverse_z_power {
        for (extra_inverse_z_power, exp_series) in &exp_terms {
            if inverse_z_power + extra_inverse_z_power != 2 {
                continue;
            }
            for (lambda_power, q_series) in &lambda_series_by_power {
                let product = truncated_series_multiply(q_series, exp_series, max_degree);
                let output_lambda_power = lambda_power
                    + usize::try_from(*extra_inverse_z_power).expect("power is nonnegative");
                for q_degree in 0..=max_degree {
                    *primary_by_q_degree[q_degree]
                        .entry(output_lambda_power)
                        .or_insert_with(|| Rational::from(0)) += product[q_degree].clone();
                }
            }
        }
    }
    primary_by_q_degree
}

fn weighted_p2_ordinary_fun0_lambda_polynomials_by_z(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
) -> Option<BTreeMap<i64, BTreeMap<usize, Rational>>> {
    if degree_twice % 2 != 0 {
        return None;
    }
    let degree = degree_twice / 2;
    if degree < 0 {
        return None;
    }
    let denominator = integer_sector_base_denominator_constant(base_weights, degree_twice)?;
    let denominator_z_power = base_weights
        .iter()
        .map(|weight| {
            weight
                .checked_mul(degree)
                .expect("weighted P2 denominator z power fits i64")
        })
        .try_fold(0i64, |sum, count| sum.checked_add(count))
        .expect("weighted P2 denominator total z power fits i64");
    let mut numerator = BTreeMap::<(usize, i64), Rational>::new();
    numerator.insert((0, 0), Rational::from(1));
    for &line_degree in bundle_degrees {
        let factor_count = line_degree
            .checked_mul(degree)
            .expect("weighted P2 split bundle factor count fits i64");
        if factor_count < 0 {
            return None;
        }
        for factor_offset in (-factor_count + 1)..=0 {
            numerator = multiply_lambda_z_polynomial_by_rational_factor(
                &numerator,
                Rational::from(1),
                rational_from_i64(factor_offset),
            );
        }
    }

    let mut by_inverse_z_power = BTreeMap::<i64, BTreeMap<usize, Rational>>::new();
    for ((lambda_power, z_power), coefficient) in numerator {
        if coefficient == Rational::from(0) {
            continue;
        }
        let inverse_z_power = denominator_z_power - z_power;
        *by_inverse_z_power
            .entry(inverse_z_power)
            .or_default()
            .entry(lambda_power)
            .or_insert_with(|| Rational::from(0)) += coefficient / denominator.clone();
    }
    Some(by_inverse_z_power)
}

fn weighted_p2_half_sector_primary_after_kp112_mirror_map(
    bundle_degrees: &[i64],
    half_sector_count: usize,
) -> Vec<BTreeMap<usize, Rational>> {
    weighted_p2_half_sector_fun_component_after_kp112_mirror_map_by_z(
        bundle_degrees,
        half_sector_count,
        2,
    )
    .into_iter()
    .map(|by_z| by_z.get(&2).cloned().unwrap_or_default())
    .collect()
}

fn weighted_p2_half_sector_fun_component_after_kp112_mirror_map_by_z(
    bundle_degrees: &[i64],
    half_sector_count: usize,
    max_inverse_z_power: i64,
) -> Vec<BTreeMap<i64, BTreeMap<usize, Rational>>> {
    if half_sector_count == 0 {
        return Vec::new();
    }
    let max_shift = half_sector_count - 1;
    let x_of_q = kp112_x_of_q_series(max_shift);
    let h_of_q = truncated_series_compose(&kp112_h_series(max_shift), &x_of_q, max_shift);
    let mut by_z_by_half_index =
        vec![BTreeMap::<i64, BTreeMap<usize, Rational>>::new(); half_sector_count];
    let exp_terms = kp112_lambda_expansion_terms(
        &h_of_q,
        max_shift,
        usize::try_from(max_inverse_z_power).expect("inverse z power is nonnegative"),
    );

    for half_index in 0..half_sector_count {
        let degree_twice = i64::try_from(2 * half_index + 1).expect("half degree fits i64");
        let Some((denominator, denominator_z_power)) =
            half_sector_base_denominator_constant_and_z_power(&[1, 1, 2], degree_twice)
        else {
            continue;
        };
        let raw_by_numerator_z_power =
            half_sector_fun_component_by_numerator_z_power(bundle_degrees, degree_twice);
        let argument = h_of_q
            .iter()
            .map(|coefficient| -rational_from_i64(2 * degree_twice) * coefficient.clone())
            .collect::<Vec<_>>();
        let x_factor = truncated_series_exp(&argument, max_shift);

        for (numerator_z_power, lambda_polynomial) in raw_by_numerator_z_power {
            let input_inverse_z_power = denominator_z_power - numerator_z_power;
            for (extra_inverse_z_power, exp_series) in &exp_terms {
                let output_inverse_z_power = input_inverse_z_power + extra_inverse_z_power;
                if output_inverse_z_power < 2 || output_inverse_z_power > max_inverse_z_power {
                    continue;
                }
                let q_factor = truncated_series_multiply(&x_factor, exp_series, max_shift);
                for (lambda_power, coefficient) in &lambda_polynomial {
                    let scaled = coefficient.clone() / denominator.clone();
                    let output_lambda_power = lambda_power
                        + usize::try_from(*extra_inverse_z_power).expect("power is nonnegative");
                    for shift in 0..=max_shift {
                        let output_half_index = half_index + shift;
                        if output_half_index >= half_sector_count {
                            continue;
                        }
                        *by_z_by_half_index[output_half_index]
                            .entry(output_inverse_z_power)
                            .or_default()
                            .entry(output_lambda_power)
                            .or_insert_with(|| Rational::from(0)) +=
                            scaled.clone() * q_factor[shift].clone();
                    }
                }
            }
        }
    }
    by_z_by_half_index
}

fn kp112_lambda_expansion_terms(
    h_of_q: &[Rational],
    max_degree: usize,
    max_power: usize,
) -> Vec<(i64, Vec<Rational>)> {
    let mut exp_terms = Vec::<(i64, Vec<Rational>)>::new();
    let mut h_power = vec![Rational::from(0); max_degree + 1];
    h_power[0] = Rational::from(1);
    exp_terms.push((0, h_power.clone()));
    for power in 1..=max_power {
        h_power = truncated_series_multiply(&h_power, h_of_q, max_degree);
        let denominator = Rational::from(factorial_integer(power));
        exp_terms.push((
            i64::try_from(power).expect("power fits i64"),
            h_power
                .iter()
                .map(|coefficient| coefficient.clone() / denominator.clone())
                .collect(),
        ));
    }
    exp_terms
}

fn half_sector_fun_component_by_numerator_z_power(
    bundle_degrees: &[i64],
    degree_twice: i64,
) -> BTreeMap<i64, BTreeMap<usize, Rational>> {
    let mut numerator = BTreeMap::<(usize, i64), Rational>::new();
    numerator.insert((0, 0), Rational::from(1));
    for &line_degree in bundle_degrees {
        for offset_twice in half_sector_numerator_offset_twice_values(line_degree, degree_twice) {
            numerator = multiply_lambda_z_polynomial_by_rational_factor(
                &numerator,
                Rational::from(1),
                Rational::from(Integer::from(offset_twice)) / Rational::from(Integer::from(2)),
            );
        }
    }
    let mut by_numerator_z_power = BTreeMap::<i64, BTreeMap<usize, Rational>>::new();
    for ((lambda_power, z_power), coefficient) in numerator {
        *by_numerator_z_power
            .entry(z_power)
            .or_default()
            .entry(lambda_power)
            .or_insert_with(|| Rational::from(0)) += coefficient;
    }
    by_numerator_z_power
}

fn half_sector_numerator_offset_twice_values(line_degree: i64, degree_twice: i64) -> Vec<i64> {
    let lower_twice = -line_degree
        .checked_mul(degree_twice)
        .expect("weighted P2 half-sector numerator bound fits i64");
    let mut offsets = Vec::new();
    let mut offset_twice = lower_twice + 2;
    while offset_twice <= 0 {
        offsets.push(offset_twice);
        offset_twice += 2;
    }
    offsets
}

fn half_sector_base_denominator_constant_and_z_power(
    base_weights: &[i64],
    degree_twice: i64,
) -> Option<(Rational, i64)> {
    if degree_twice % 2 == 0 || degree_twice < 0 {
        return None;
    }
    let mut denominator = Rational::from(1);
    let mut z_power = 0i64;
    for &weight in base_weights {
        let upper_twice = weight
            .checked_mul(degree_twice)
            .expect("weighted P2 half-sector denominator bound fits i64");
        if upper_twice <= 0 {
            return None;
        }
        let mut factor_twice = if upper_twice % 2 == 0 { 2 } else { 1 };
        while factor_twice <= upper_twice {
            denominator *=
                Rational::from(Integer::from(factor_twice)) / Rational::from(Integer::from(2));
            z_power += 1;
            factor_twice += 2;
        }
    }
    Some((denominator, z_power))
}

fn integer_sector_base_denominator_constant(
    base_weights: &[i64],
    degree_twice: i64,
) -> Option<Rational> {
    if degree_twice % 2 != 0 {
        return None;
    }
    let degree = degree_twice / 2;
    if degree < 0 {
        return None;
    }
    let mut product = Integer::from(1);
    for &weight in base_weights {
        let factor_count = weight
            .checked_mul(degree)
            .expect("weighted P2 denominator factor count fits i64");
        if factor_count < 0 {
            return None;
        }
        for offset in 1..=factor_count {
            product *= Integer::from(offset);
        }
    }
    Some(Rational::from(product))
}

fn multiply_lambda_z_polynomial_by_rational_factor(
    coefficients: &BTreeMap<(usize, i64), Rational>,
    lambda_coefficient: Rational,
    z_coefficient: Rational,
) -> BTreeMap<(usize, i64), Rational> {
    let mut next = BTreeMap::<(usize, i64), Rational>::new();
    for (&(lambda_power, z_power), coefficient) in coefficients {
        if lambda_coefficient != Rational::from(0) {
            *next
                .entry((lambda_power + 1, z_power))
                .or_insert_with(|| Rational::from(0)) +=
                coefficient.clone() * lambda_coefficient.clone();
        }
        if z_coefficient != Rational::from(0) {
            *next
                .entry((lambda_power, z_power + 1))
                .or_insert_with(|| Rational::from(0)) +=
                coefficient.clone() * z_coefficient.clone();
        }
    }
    next.retain(|_, coefficient| *coefficient != Rational::from(0));
    next
}

fn divide_rational_lambda_polynomial_by_two_lambda_if_possible(
    polynomial: &BTreeMap<usize, Rational>,
) -> Option<Vec<Rational>> {
    divide_rational_lambda_polynomial_by_two_lambda_power_if_possible(polynomial, 1)
}

fn divide_rational_lambda_polynomial_by_two_lambda_power_if_possible(
    polynomial: &BTreeMap<usize, Rational>,
    lambda_power: usize,
) -> Option<Vec<Rational>> {
    for lower_power in 0..lambda_power {
        let coefficient = polynomial
            .get(&lower_power)
            .cloned()
            .unwrap_or_else(|| Rational::from(0));
        if coefficient != Rational::from(0) {
            return None;
        }
    }
    let max_lambda_power = polynomial.keys().last().copied().unwrap_or(0);
    let two = Rational::from(Integer::from(2));
    Some(
        (lambda_power..=max_lambda_power)
            .map(|lambda_power| {
                polynomial
                    .get(&lambda_power)
                    .cloned()
                    .unwrap_or_else(|| Rational::from(0))
                    / two.clone()
            })
            .collect(),
    )
}

fn rational_coefficients_to_strings(coefficients: &[Rational]) -> Vec<String> {
    coefficients
        .iter()
        .map(std::string::ToString::to_string)
        .collect()
}

fn lambda_polynomial_has_nonzero_coefficient(polynomial: &[String]) -> bool {
    polynomial.iter().any(|coefficient| coefficient != "0")
}

fn split_non_equivariant_candidate_insertion_coefficient_status(
    degree_twice: i64,
    numerator_zero_factor_order: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> (Option<String>, String) {
    if degree_twice % 2 != 0 {
        return (
            None,
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_p2_coefficient_readout"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            "weighted_p2_rank_three_split_non_equivariant_insertion_codim_not_applicable"
                .to_string(),
        );
    };
    match numerator_zero_factor_order.cmp(&required_codimension) {
        std::cmp::Ordering::Greater => (
            Some("0".to_string()),
            "weighted_p2_rank_three_split_untwisted_non_equivariant_p2_coefficient_zero_before_mirror_map_due_to_excess_zero_order"
                .to_string(),
        ),
        std::cmp::Ordering::Equal => (
            None,
            "weighted_p2_rank_three_split_untwisted_non_equivariant_p2_coefficient_requires_exact_factor_expansion"
                .to_string(),
        ),
        std::cmp::Ordering::Less => (
            None,
            "weighted_p2_rank_three_split_untwisted_non_equivariant_p2_coefficient_requires_higher_series_terms"
                .to_string(),
        ),
    }
}

fn split_non_equivariant_bundle_numerator_truncation_status(
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> (Option<Vec<String>>, String) {
    if degree_twice % 2 != 0 {
        return (
            None,
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_bundle_numerator"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            "weighted_p2_rank_three_split_non_equivariant_bundle_numerator_codim_not_applicable"
                .to_string(),
        );
    };
    if required_codimension < 0 {
        return (
            None,
            "weighted_p2_rank_three_split_non_equivariant_bundle_numerator_invalid_negative_codim"
                .to_string(),
        );
    }
    let coefficients = non_equivariant_split_bundle_numerator_truncated_coefficients(
        bundle_degrees,
        degree_twice,
        required_codimension as usize,
    );
    let status = if coefficients[required_codimension as usize] == "0" {
        "weighted_p2_rank_three_split_non_equivariant_bundle_numerator_candidate_codim_coefficient_zero"
    } else {
        "weighted_p2_rank_three_split_non_equivariant_bundle_numerator_candidate_codim_coefficient_nonzero_requires_denominator_and_mirror_map"
    };
    (Some(coefficients), status.to_string())
}

fn non_equivariant_split_bundle_numerator_truncated_coefficients(
    bundle_degrees: &[i64],
    degree_twice: i64,
    max_power: usize,
) -> Vec<String> {
    assert_eq!(
        degree_twice % 2,
        0,
        "non-equivariant untwisted sector has integer degree"
    );
    let degree = degree_twice / 2;
    let mut coefficients = vec![Rational::from(0); max_power + 1];
    coefficients[0] = Rational::from(1);
    for &line_degree in bundle_degrees {
        let factor_count = line_degree
            .checked_mul(degree)
            .expect("weighted P2 split bundle factor count fits i64");
        for factor_offset in (-factor_count + 1)..=0 {
            multiply_truncated_linear_factor(
                &mut coefficients,
                factor_offset,
                -line_degree,
                max_power,
            );
        }
    }
    rational_coefficients_to_strings(&coefficients)
}

fn multiply_truncated_linear_factor(
    coefficients: &mut Vec<Rational>,
    constant: i64,
    linear: i64,
    max_power: usize,
) {
    let mut next = vec![Rational::from(0); max_power + 1];
    for power in 0..=max_power {
        let coefficient = coefficients[power].clone();
        if coefficient == Rational::from(0) {
            continue;
        }
        next[power] += coefficient.clone() * rational_from_i64(constant);
        if power < max_power {
            next[power + 1] += coefficient * rational_from_i64(linear);
        }
    }
    *coefficients = next;
}

fn split_equivariant_candidate_insertion_status(
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> (Option<Vec<String>>, Option<i64>, String) {
    if degree_twice % 2 != 0 {
        return (
            None,
            None,
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_equivariant_p2_readout"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            None,
            "weighted_p2_rank_three_split_equivariant_insertion_codim_not_applicable".to_string(),
        );
    };
    if required_codimension < 0 {
        return (
            None,
            None,
            "weighted_p2_rank_three_split_equivariant_insertion_invalid_negative_codim".to_string(),
        );
    }
    let polynomial = equivariant_split_bundle_numerator_candidate_lambda_polynomial(
        bundle_degrees,
        degree_twice,
        required_codimension as usize,
    );
    let order = lambda_polynomial_order(&polynomial);
    let status = match order {
        None => "weighted_p2_rank_three_split_equivariant_candidate_p2_coefficient_zero",
        Some(0) => {
            "weighted_p2_rank_three_split_equivariant_candidate_p2_has_nonzero_nonequivariant_term"
        }
        Some(_) => {
            "weighted_p2_rank_three_split_equivariant_candidate_p2_positive_lambda_order_requires_pairing_or_residue"
        }
    };
    (Some(polynomial), order, status.to_string())
}

fn split_equivariant_dual_basis_p2_readout_status(
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> (Option<Vec<String>>, Option<i64>, String) {
    if degree_twice % 2 != 0 {
        return (
            None,
            None,
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_dual_basis_p2_readout"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_readout_codim_not_applicable".to_string(),
        );
    };
    if required_codimension != 2 {
        return (
            None,
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_readout_requires_codim2_observable"
                .to_string(),
        );
    }
    let ordinary_fun0_polynomial = equivariant_split_bundle_numerator_candidate_lambda_polynomial(
        bundle_degrees,
        degree_twice,
        0,
    );
    let Some(quotient) = divide_lambda_polynomial_by_two_lambda_power_if_possible(
        &ordinary_fun0_polynomial,
        bundle_degrees.len(),
    ) else {
        return (
            Some(ordinary_fun0_polynomial),
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_readout_blocked_fun0_coefficient_not_divisible_by_two_lambda"
                .to_string(),
        );
    };
    let order = lambda_polynomial_order(&quotient);
    let status = match order {
        None => "weighted_p2_rank_three_split_dual_basis_p2_numerator_readout_zero",
        Some(0) => {
            "weighted_p2_rank_three_split_dual_basis_p2_numerator_readout_has_nonzero_nonequivariant_term"
        }
        Some(_) => {
            "weighted_p2_rank_three_split_dual_basis_p2_numerator_readout_positive_lambda_order_requires_full_ifunction_pairing_or_residue"
        }
    };
    (Some(quotient), order, status.to_string())
}

fn split_equivariant_dual_basis_p2_hypergeometric_status(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> (Option<String>, Option<Vec<String>>, Option<i64>, String) {
    if degree_twice % 2 != 0 {
        return (
            None,
            None,
            None,
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_dual_basis_p2_hypergeometric_readout"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            None,
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_codim_not_applicable"
                .to_string(),
        );
    };
    if required_codimension != 2 {
        return (
            None,
            None,
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_requires_codim2_observable"
                .to_string(),
        );
    }
    let (numerator_readout, _, _) = split_equivariant_dual_basis_p2_readout_status(
        bundle_degrees,
        degree_twice,
        required_insertion_complex_codimension,
    );
    let Some(numerator_readout) = numerator_readout else {
        return (
            None,
            None,
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_blocked_missing_numerator_readout"
                .to_string(),
        );
    };
    let Some(denominator) = integer_sector_base_denominator_constant(base_weights, degree_twice)
    else {
        return (
            None,
            None,
            None,
            "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_blocked_missing_base_denominator"
                .to_string(),
        );
    };
    let hypergeometric = divide_lambda_polynomial_by_rational(&numerator_readout, &denominator);
    let order = lambda_polynomial_order(&hypergeometric);
    let status = match order {
        None => "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_readout_zero",
        Some(0) => {
            "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_readout_has_nonzero_nonequivariant_term"
        }
        Some(_) => {
            "weighted_p2_rank_three_split_dual_basis_p2_hypergeometric_readout_positive_lambda_order_requires_mirror_map_pairing_or_residue"
        }
    };
    (
        Some(denominator.to_string()),
        Some(hypergeometric),
        order,
        status.to_string(),
    )
}

fn split_equivariant_full_hypergeometric_status(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
    max_p_power: usize,
) -> (
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<Vec<Vec<String>>>,
    Option<Vec<String>>,
    Option<i64>,
    String,
) {
    if degree_twice % 2 != 0 {
        return (
            None,
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_full_hypergeometric_readout"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_full_hypergeometric_codim_not_applicable".to_string(),
        );
    };
    if required_codimension != 2 {
        return (
            None,
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_full_hypergeometric_requires_codim2_observable"
                .to_string(),
        );
    }
    let Some(denominator) = integer_sector_base_denominator_truncated_coefficients(
        base_weights,
        degree_twice,
        max_p_power,
    ) else {
        return (
            None,
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_full_hypergeometric_blocked_missing_base_denominator"
                .to_string(),
        );
    };
    let Some(inverse_denominator) = invert_truncated_univariate_polynomial(&denominator) else {
        return (
            Some(rational_coefficients_to_strings(&denominator)),
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_full_hypergeometric_blocked_noninvertible_base_denominator"
                .to_string(),
        );
    };
    let numerator = equivariant_split_bundle_numerator_bivariate_coefficients(
        bundle_degrees,
        degree_twice,
        max_p_power,
    );
    let quotient =
        multiply_bivariate_by_univariate_truncated(&numerator, &inverse_denominator, max_p_power);
    let quotient_polynomials =
        bivariate_coefficients_by_p_power_to_lambda_polynomials(&quotient, max_p_power);
    let Some(fun0_polynomial) = quotient_polynomials.first() else {
        return (
            Some(rational_coefficients_to_strings(&denominator)),
            Some(rational_coefficients_to_strings(&inverse_denominator)),
            Some(quotient_polynomials),
            None,
            None,
            "weighted_p2_rank_three_split_full_hypergeometric_blocked_missing_fun0_coefficient"
                .to_string(),
        );
    };
    let Some(readout) = divide_lambda_polynomial_by_two_lambda_power_if_possible(
        fun0_polynomial,
        bundle_degrees.len(),
    ) else {
        return (
            Some(rational_coefficients_to_strings(&denominator)),
            Some(rational_coefficients_to_strings(&inverse_denominator)),
            Some(quotient_polynomials),
            None,
            None,
            "weighted_p2_rank_three_split_full_hypergeometric_blocked_fun0_coefficient_not_divisible_by_two_lambda"
                .to_string(),
        );
    };
    let order = lambda_polynomial_order(&readout);
    let status = match order {
        None => "weighted_p2_rank_three_split_full_hypergeometric_dual_basis_p2_readout_zero",
        Some(0) => {
            "weighted_p2_rank_three_split_full_hypergeometric_dual_basis_p2_readout_has_nonzero_nonequivariant_term"
        }
        Some(_) => {
            "weighted_p2_rank_three_split_full_hypergeometric_dual_basis_p2_readout_positive_lambda_order_requires_mirror_map_pairing_or_residue"
        }
    };
    (
        Some(rational_coefficients_to_strings(&denominator)),
        Some(rational_coefficients_to_strings(&inverse_denominator)),
        Some(quotient_polynomials),
        Some(readout),
        order,
        status.to_string(),
    )
}

fn split_equivariant_full_hypergeometric_dual_basis_p2_z_readout_status(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> (
    Option<Vec<String>>,
    Option<i64>,
    Option<Vec<String>>,
    Option<i64>,
    String,
    String,
) {
    if degree_twice % 2 != 0 {
        return (
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_first_nonzero_z_readout_not_applicable_twisted_sector"
                .to_string(),
            "weighted_p2_rank_three_split_half_degree_twisted_sector_no_untwisted_full_hypergeometric_z_readout"
                .to_string(),
        );
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return (
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_first_nonzero_z_readout_codim_not_applicable".to_string(),
            "weighted_p2_rank_three_split_full_hypergeometric_z_readout_codim_not_applicable"
                .to_string(),
        );
    };
    if required_codimension != 2 {
        return (
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_first_nonzero_z_readout_requires_codim2_observable"
                .to_string(),
            "weighted_p2_rank_three_split_full_hypergeometric_z_readout_requires_codim2_observable"
                .to_string(),
        );
    }
    let Some(readout_by_z) = weighted_p2_ordinary_dual_basis_p2_z_lambda_polynomials(
        base_weights,
        bundle_degrees,
        degree_twice,
    ) else {
        return (
            None,
            None,
            None,
            None,
            "weighted_p2_rank_three_split_first_nonzero_z_readout_blocked_missing_z_readout"
                .to_string(),
            "weighted_p2_rank_three_split_full_hypergeometric_z_readout_blocked_not_divisible_by_two_lambda"
                .to_string(),
        );
    };
    let primary_z2 = readout_by_z
        .get(&2)
        .cloned()
        .unwrap_or_else(|| vec!["0".to_string()]);
    let first_nonzero = readout_by_z
        .iter()
        .filter(|(inverse_z_power, _)| **inverse_z_power >= 2)
        .find(|(_, polynomial)| lambda_polynomial_has_nonzero_coefficient(polynomial))
        .map(|(&inverse_z_power, polynomial)| (inverse_z_power, polynomial.clone()));
    let first_nonzero_lambda_order = first_nonzero
        .as_ref()
        .and_then(|(_, polynomial)| lambda_polynomial_order(polynomial));
    let first_nonzero_non_equivariant_status = match (
        first_nonzero
            .as_ref()
            .map(|(inverse_z_power, _)| *inverse_z_power),
        first_nonzero_lambda_order,
    ) {
        (None, _) => "weighted_p2_rank_three_split_first_nonzero_z_readout_absent",
        (Some(2), Some(0)) => {
            "weighted_p2_rank_three_split_first_nonzero_primary_z2_has_nonequivariant_limit"
        }
        (Some(2), Some(_)) => {
            "weighted_p2_rank_three_split_first_nonzero_primary_z2_positive_lambda_order_requires_pairing_or_residue"
        }
        (Some(_), Some(0)) => {
            "weighted_p2_rank_three_split_first_nonzero_descendant_has_nonequivariant_limit_requires_reconstruction"
        }
        (Some(_), Some(_)) => {
            "weighted_p2_rank_three_split_first_nonzero_descendant_positive_lambda_order_equivariant_only_requires_residue_or_pairing"
        }
        (Some(_), None) => "weighted_p2_rank_three_split_first_nonzero_z_readout_order_unresolved",
    };
    let status = if lambda_polynomial_has_nonzero_coefficient(&primary_z2) {
        "weighted_p2_rank_three_split_full_hypergeometric_dual_basis_p2_primary_z2_nonzero_requires_pairing_or_residue"
    } else if first_nonzero.is_some() {
        "weighted_p2_rank_three_split_full_hypergeometric_dual_basis_p2_primary_z2_zero_first_nonzero_descendant_requires_big_j_or_pairing"
    } else {
        "weighted_p2_rank_three_split_full_hypergeometric_dual_basis_p2_z_readout_zero"
    };
    (
        Some(primary_z2),
        first_nonzero
            .as_ref()
            .map(|(inverse_z_power, _)| *inverse_z_power),
        first_nonzero.map(|(_, polynomial)| polynomial),
        first_nonzero_lambda_order,
        first_nonzero_non_equivariant_status.to_string(),
        status.to_string(),
    )
}

fn split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status(
    base_weights: &[i64],
    bundle_degrees: &[i64],
    degree_twice: i64,
    required_insertion_complex_codimension: Option<i64>,
) -> String {
    if degree_twice % 2 != 0 {
        return "weighted_p2_rank_three_split_scalar_mirror_map_primary_not_applicable_twisted_sector"
            .to_string();
    }
    let Some(required_codimension) = required_insertion_complex_codimension else {
        return "weighted_p2_rank_three_split_scalar_mirror_map_primary_codim_not_applicable"
            .to_string();
    };
    if required_codimension != 2 {
        return "weighted_p2_rank_three_split_scalar_mirror_map_primary_requires_codim2_observable"
            .to_string();
    }
    let Some(readout_by_z) = weighted_p2_ordinary_dual_basis_p2_z_lambda_polynomials(
        base_weights,
        bundle_degrees,
        degree_twice,
    ) else {
        return "weighted_p2_rank_three_split_scalar_mirror_map_primary_blocked_missing_z_readout"
            .to_string();
    };
    let current_primary = readout_by_z
        .get(&2)
        .cloned()
        .unwrap_or_else(|| vec!["0".to_string()]);
    if lambda_polynomial_has_nonzero_coefficient(&current_primary) {
        return "weighted_p2_rank_three_split_scalar_mirror_map_primary_already_nonzero_before_reparametrization"
            .to_string();
    }
    for prior_degree_twice in (2..degree_twice).step_by(2) {
        let Some(prior_readout_by_z) = weighted_p2_ordinary_dual_basis_p2_z_lambda_polynomials(
            base_weights,
            bundle_degrees,
            prior_degree_twice,
        ) else {
            return "weighted_p2_rank_three_split_scalar_mirror_map_primary_blocked_missing_lower_degree_z_readout"
                .to_string();
        };
        let prior_primary = prior_readout_by_z
            .get(&2)
            .cloned()
            .unwrap_or_else(|| vec!["0".to_string()]);
        if lambda_polynomial_has_nonzero_coefficient(&prior_primary) {
            return "weighted_p2_rank_three_split_scalar_mirror_map_primary_zero_but_lower_primary_terms_can_mix_by_degree_reparametrization"
                .to_string();
        }
    }
    "weighted_p2_rank_three_split_scalar_mirror_map_primary_zero_preserved_by_scalar_reparametrization_to_this_degree".to_string()
}

fn lambda_polynomial_order(polynomial: &[String]) -> Option<i64> {
    polynomial
        .iter()
        .position(|coefficient| coefficient != "0")
        .map(|order| order as i64)
}

fn divide_lambda_polynomial_by_two_lambda_power_if_possible(
    polynomial: &[String],
    lambda_power: usize,
) -> Option<Vec<String>> {
    for coefficient in polynomial.iter().take(lambda_power) {
        if coefficient != "0" {
            return None;
        }
    }
    let two = Rational::from(Integer::from(2));
    Some(
        polynomial
            .iter()
            .skip(lambda_power)
            .map(|coefficient| {
                let coefficient = coefficient
                    .parse::<Rational>()
                    .expect("lambda polynomial coefficient is rational");
                (coefficient / two.clone()).to_string()
            })
            .collect(),
    )
}

fn integer_sector_base_denominator_truncated_coefficients(
    base_weights: &[i64],
    degree_twice: i64,
    max_p_power: usize,
) -> Option<Vec<Rational>> {
    if degree_twice % 2 != 0 {
        return None;
    }
    let degree = degree_twice / 2;
    if degree < 0 {
        return None;
    }
    let mut coefficients = vec![Rational::from(0); max_p_power + 1];
    coefficients[0] = Rational::from(1);
    for &weight in base_weights {
        let factor_count = weight
            .checked_mul(degree)
            .expect("weighted P2 denominator factor count fits i64");
        if factor_count < 0 {
            return None;
        }
        for offset in 1..=factor_count {
            multiply_univariate_truncated_linear_factor(
                &mut coefficients,
                offset,
                weight,
                max_p_power,
            );
        }
    }
    Some(coefficients)
}

fn multiply_univariate_truncated_linear_factor(
    coefficients: &mut Vec<Rational>,
    constant: i64,
    linear: i64,
    max_power: usize,
) {
    let mut next = vec![Rational::from(0); max_power + 1];
    for power in 0..=max_power {
        let coefficient = coefficients[power].clone();
        if coefficient == Rational::from(0) {
            continue;
        }
        next[power] += coefficient.clone() * rational_from_i64(constant);
        if power < max_power {
            next[power + 1] += coefficient * rational_from_i64(linear);
        }
    }
    *coefficients = next;
}

fn invert_truncated_univariate_polynomial(coefficients: &[Rational]) -> Option<Vec<Rational>> {
    let constant = coefficients.first()?;
    if constant == &Rational::from(0) {
        return None;
    }
    let mut inverse = vec![Rational::from(0); coefficients.len()];
    inverse[0] = Rational::from(1) / constant.clone();
    for power in 1..coefficients.len() {
        let mut sum = Rational::from(0);
        for offset in 1..=power {
            sum += coefficients[offset].clone() * inverse[power - offset].clone();
        }
        inverse[power] = -sum / constant.clone();
    }
    Some(inverse)
}

fn divide_lambda_polynomial_by_rational(
    polynomial: &[String],
    denominator: &Rational,
) -> Vec<String> {
    polynomial
        .iter()
        .map(|coefficient| {
            let coefficient = coefficient
                .parse::<Rational>()
                .expect("lambda polynomial coefficient is rational");
            (coefficient / denominator.clone()).to_string()
        })
        .collect()
}

fn equivariant_split_bundle_numerator_candidate_lambda_polynomial(
    bundle_degrees: &[i64],
    degree_twice: i64,
    candidate_p_power: usize,
) -> Vec<String> {
    let coefficients = equivariant_split_bundle_numerator_bivariate_coefficients(
        bundle_degrees,
        degree_twice,
        candidate_p_power,
    );
    bivariate_coefficients_by_p_power_to_lambda_polynomials(&coefficients, candidate_p_power)
        .into_iter()
        .nth(candidate_p_power)
        .expect("candidate p-power polynomial is present")
}

fn equivariant_split_bundle_numerator_bivariate_coefficients(
    bundle_degrees: &[i64],
    degree_twice: i64,
    max_p_power: usize,
) -> BTreeMap<(usize, usize), Rational> {
    assert_eq!(
        degree_twice % 2,
        0,
        "equivariant untwisted sector has integer degree"
    );
    let degree = degree_twice / 2;
    let mut coefficients = BTreeMap::<(usize, usize), Rational>::new();
    coefficients.insert((0, 0), Rational::from(1));
    for &line_degree in bundle_degrees {
        let factor_count = line_degree
            .checked_mul(degree)
            .expect("weighted P2 split bundle factor count fits i64");
        for factor_offset in (-factor_count + 1)..=0 {
            coefficients = multiply_bivariate_truncated_factor(
                &coefficients,
                factor_offset,
                1,
                -line_degree,
                max_p_power,
            );
        }
    }
    coefficients
}

fn multiply_bivariate_by_univariate_truncated(
    coefficients: &BTreeMap<(usize, usize), Rational>,
    univariate: &[Rational],
    max_p_power: usize,
) -> BTreeMap<(usize, usize), Rational> {
    let mut next = BTreeMap::<(usize, usize), Rational>::new();
    for (&(lambda_power, p_power), coefficient) in coefficients {
        for (factor_p_power, factor_coefficient) in univariate.iter().enumerate() {
            let next_p_power = p_power + factor_p_power;
            if next_p_power > max_p_power {
                continue;
            }
            if factor_coefficient == &Rational::from(0) {
                continue;
            }
            *next
                .entry((lambda_power, next_p_power))
                .or_insert_with(|| Rational::from(0)) +=
                coefficient.clone() * factor_coefficient.clone();
        }
    }
    next.retain(|_, coefficient| *coefficient != Rational::from(0));
    next
}

fn bivariate_coefficients_by_p_power_to_lambda_polynomials(
    coefficients: &BTreeMap<(usize, usize), Rational>,
    max_p_power: usize,
) -> Vec<Vec<String>> {
    (0..=max_p_power)
        .map(|target_p_power| {
            let mut lambda_coefficients = BTreeMap::<usize, Rational>::new();
            for (&(lambda_power, p_power), coefficient) in coefficients {
                if p_power == target_p_power {
                    *lambda_coefficients
                        .entry(lambda_power)
                        .or_insert_with(|| Rational::from(0)) += coefficient.clone();
                }
            }
            let max_lambda_power = lambda_coefficients.keys().last().copied().unwrap_or(0);
            (0..=max_lambda_power)
                .map(|lambda_power| {
                    lambda_coefficients
                        .get(&lambda_power)
                        .cloned()
                        .unwrap_or_else(|| Rational::from(0))
                        .to_string()
                })
                .collect()
        })
        .collect()
}

fn multiply_bivariate_truncated_factor(
    coefficients: &BTreeMap<(usize, usize), Rational>,
    constant: i64,
    lambda: i64,
    p: i64,
    max_p_power: usize,
) -> BTreeMap<(usize, usize), Rational> {
    let factor_terms = [
        ((0usize, 0usize), constant),
        ((1usize, 0usize), lambda),
        ((0usize, 1usize), p),
    ];
    let mut next = BTreeMap::<(usize, usize), Rational>::new();
    for (&(lambda_power, p_power), coefficient) in coefficients {
        for &((factor_lambda_power, factor_p_power), factor_coefficient) in &factor_terms {
            let next_p_power = p_power + factor_p_power;
            if next_p_power > max_p_power {
                continue;
            }
            *next
                .entry((lambda_power + factor_lambda_power, next_p_power))
                .or_insert_with(|| Rational::from(0)) +=
                coefficient.clone() * rational_from_i64(factor_coefficient);
        }
    }
    next.retain(|_, coefficient| *coefficient != Rational::from(0));
    next
}

fn ceil_positive_half_product(factor: i64, degree_twice: i64) -> i64 {
    let product = factor
        .checked_mul(degree_twice)
        .expect("weighted P2 I-function factor product fits i64");
    (product + 1) / 2
}

fn floor_positive_half_product(factor: i64, degree_twice: i64) -> i64 {
    let product = factor
        .checked_mul(degree_twice)
        .expect("weighted P2 I-function factor product fits i64");
    product / 2
}

#[cfg(test)]
mod tests {
    use super::{
        kp112_canonical_b_table_values, kp112_canonical_c_table_values,
        kp112_half_sector_first_nonzero_descendants, kp112_half_sector_primary_lambda_coefficients,
        kp112_ordinary_p_primary_lambda_coefficients,
        split_bundle_kp112_mirror_map_half_sector_descendant_profiles,
        split_bundle_kp112_mirror_map_half_sector_first_nonzero_descendants,
        split_bundle_kp112_mirror_map_half_sector_primary_coefficients,
        split_bundle_kp112_mirror_map_primary_p_lambda_coefficients,
        weighted_p2_ordinary_dual_basis_p2_z_readout_profile,
        weighted_p2_rank_three_split_bundle_chen_ruan_source_basis_readout,
        weighted_p2_rank_three_split_bundle_source_readiness,
        weighted_p2_rank_three_twisted_ifunction_degree_profiles,
        weighted_p2_rank_three_twisted_ifunction_source_certificate_requirements,
    };

    #[test]
    fn canonical_weighted_p2_mirror_map_matches_ccit_b_table() {
        assert_eq!(
            kp112_canonical_b_table_values(6),
            vec![
                "11/4".to_string(),
                "525/16".to_string(),
                "6152/9".to_string(),
                "1146765/64".to_string(),
                "53305261/100".to_string(),
                "51550873/3".to_string(),
            ]
        );
    }

    #[test]
    fn canonical_weighted_p2_mirror_map_matches_ccit_c_table() {
        assert_eq!(
            kp112_canonical_c_table_values(7),
            vec![
                "-2".to_string(),
                "-52/9".to_string(),
                "-2002/25".to_string(),
                "-83004/49".to_string(),
                "-3554552/81".to_string(),
                "-154984300/121".to_string(),
                "-6835086702/169".to_string(),
            ]
        );
    }

    #[test]
    fn generic_weighted_p2_line_bundle_api_exposes_primary_coefficients() {
        let ordinary_primary = kp112_ordinary_p_primary_lambda_coefficients(&[4], 1);
        assert_eq!(ordinary_primary[0].get(&1).map(String::as_str), Some("-44"));

        let half_primary = kp112_half_sector_primary_lambda_coefficients(&[4], 1);
        assert_eq!(half_primary[0].get(&1).map(String::as_str), Some("-4"));
    }

    #[test]
    fn generic_weighted_p2_split_bundle_api_matches_mcallister_wrappers() {
        assert_eq!(
            kp112_half_sector_first_nonzero_descendants(&[1, 1, 2], 4, 4),
            split_bundle_kp112_mirror_map_half_sector_first_nonzero_descendants(4, 4)
        );
    }

    #[test]
    fn ordinary_sector_dual_basis_p2_z_profile_matches_split_bundle_diagnostic() {
        let degree_one_profile =
            weighted_p2_ordinary_dual_basis_p2_z_readout_profile(&[1, 1, 2], &[1, 1, 2], 2)
                .expect("integer sector profile");
        assert_eq!(
            degree_one_profile
                .iter()
                .map(|entry| (
                    entry.inverse_z_power.as_str(),
                    entry.dual_basis_p2_coefficients.clone()
                ))
                .collect::<Vec<_>>(),
            vec![
                ("2", Some(vec!["0".to_string()])),
                ("3", Some(vec!["-1/4".to_string()])),
                ("4", Some(vec!["0".to_string(), "1/4".to_string()])),
            ]
        );

        let degree_two_profile =
            weighted_p2_ordinary_dual_basis_p2_z_readout_profile(&[1, 1, 2], &[1, 1, 2], 4)
                .expect("integer sector profile");
        assert_eq!(
            degree_two_profile
                .iter()
                .find(|entry| entry.inverse_z_power == "3")
                .and_then(|entry| entry.dual_basis_p2_coefficients.clone()),
            Some(vec!["-1/32".to_string()])
        );
    }

    #[test]
    fn ordinary_sector_dual_basis_p2_z_profile_matches_canonical_kp112_first_sector() {
        let profile = weighted_p2_ordinary_dual_basis_p2_z_readout_profile(&[1, 1, 2], &[4], 2)
            .expect("canonical sector profile");
        assert_eq!(
            profile
                .iter()
                .find(|entry| entry.inverse_z_power == "2")
                .and_then(|entry| entry.dual_basis_p2_coefficients.clone()),
            Some(vec!["0".to_string(), "11/4".to_string()])
        );
    }

    #[test]
    fn weighted_p2_rank_three_chen_ruan_source_basis_is_pinned() {
        let readout = weighted_p2_rank_three_split_bundle_chen_ruan_source_basis_readout(
            &[1, 1, 2],
            &[1, 1, 2],
            0,
        )
        .expect("rank-three weighted P2 source basis");
        assert_eq!(
            readout.source_status,
            "weighted_p2_rank_three_chern_ruan_source_basis_for_p112_identified"
        );
        assert_eq!(readout.inertia_component_classes, ["fun_0", "fun_{1/2}"]);
        assert_eq!(readout.base_hyperplane_square, "1/2");
        assert_eq!(readout.source_dual_basis_p2_class, "2*lambda*fun_0-8*p");
        assert_eq!(
            readout.source_dual_basis_p2_normalization_divisor,
            "2*lambda"
        );
        assert_eq!(
            readout.split_bundle_twisted_dual_basis_p2_class,
            "2*(lambda-p)^2*(lambda-2*p)*fun_0=2*lambda^3*fun_0-8*lambda^2*p+10*lambda*p^2"
        );
        assert_eq!(
            readout.split_bundle_twisted_dual_basis_p2_normalization_divisor,
            "2*lambda^3"
        );
        assert_eq!(
            readout.split_bundle_twisted_dual_basis_p2_pairing_status,
            "weighted_p2_rank_three_split_bundle_untwisted_p2_dual_pairing_computed_from_ccit_inverse_euler_twist"
        );
        assert_eq!(
            readout.source_dual_basis_half_sector_class,
            "2*lambda*fun_{1/2}"
        );
        assert_eq!(
            readout.stack_pairing_normalization_status,
            "weighted_p2_hyperplane_square_fractional_requires_stack_normalized_source_pairing_and_twisted_sector_pairing"
        );
        assert_eq!(
            readout.split_bundle_promotion_status,
            "weighted_p2_rank_three_split_bundle_still_requires_twisted_vector_bundle_ifunction_normalization_and_qn_history"
        );
        assert!(
            weighted_p2_rank_three_split_bundle_chen_ruan_source_basis_readout(
                &[1, 1, 2],
                &[4],
                0,
            )
            .is_none()
        );
    }

    #[test]
    fn weighted_p2_rank_three_source_readiness_keeps_stack_pairing_blocker_explicit() {
        let readiness =
            weighted_p2_rank_three_split_bundle_source_readiness(&[1, 1, 2], &[1, 1, 2])
                .expect("rank-three weighted P2 source readiness");
        assert_eq!(readiness.base_complex_dimension, 2);
        assert_eq!(readiness.bundle_rank, 3);
        assert_eq!(readiness.total_first_chern_degree, 0);
        assert_eq!(readiness.total_space_complex_dimension, 5);
        assert_eq!(readiness.required_cygv_codimension_for_threefold, 2);
        assert_eq!(readiness.base_hyperplane_square, "1/2");
        assert_eq!(
            readiness.base_tensor_status,
            "weighted_p2_rank_three_base_hyperplane_square_fractional_requires_stack_or_source_tensor_normalization"
        );
        assert_eq!(
            readiness.twisted_vector_bundle_ifunction_candidate_insertion_class,
            Some("base_hyperplane_power_2".to_string())
        );
        assert_eq!(
            readiness.twisted_vector_bundle_ifunction_readiness_status,
            "weighted_p2_rank_three_twisted_ifunction_blocked_missing_stack_normalized_codim2_insertion_qn_history"
        );
        assert_eq!(
            readiness.twisted_vector_bundle_ifunction_missing_inputs,
            vec![
                "source_derived_codim2_insertion_or_equivalent_observable".to_string(),
                "stack_normalized_hyperplane_square_tensor".to_string(),
                "orbifold_sector_pairing_data".to_string(),
                "equivariant_residue_or_pairing_normalization".to_string(),
                "twisted_vector_bundle_ifunction_chamber_certificate".to_string(),
                "twisted_vector_bundle_ifunction_qn_history".to_string(),
            ]
        );
        assert!(weighted_p2_rank_three_split_bundle_source_readiness(&[1, 1, 2], &[4]).is_none());
    }

    #[test]
    fn weighted_p2_rank_three_source_certificate_pins_theorem_handoff_requirements() {
        let certificate = weighted_p2_rank_three_twisted_ifunction_source_certificate_requirements(
            &[1, 1, 2],
            &[1, 1, 2],
            4,
        )
        .expect("rank-three weighted P2 source certificate");
        assert_eq!(
            certificate.theorem_handoff_status,
            "weighted_p2_rank_three_ccit_ifunction_lies_on_lagrangian_cone_not_yet_numerical_gv"
        );
        assert_eq!(
            certificate.split_bundle_modification_source_status,
            "ccit_smalllinebundle_smallvb_direct_sum_modification_product_over_o_minus_1_o_minus_1_o_minus_2"
        );
        assert_eq!(
            certificate.chen_ruan_source_basis_status,
            "weighted_p2_rank_three_chern_ruan_source_basis_for_p112_identified"
        );
        assert_eq!(
            certificate.required_observable_status,
            "weighted_p2_rank_three_source_certificate_requires_base_p2_codim2_observable_or_equivalent"
        );
        assert_eq!(
            certificate.primary_readout_status,
            "weighted_p2_rank_three_source_certificate_primary_z2_p2_readout_zero_to_checked_integer_degrees"
        );
        assert_eq!(
            certificate.descendant_readout_status,
            "weighted_p2_rank_three_source_certificate_first_nonzero_terms_are_descendant_or_equivariant_requires_big_j_pairing"
        );
        assert_eq!(
            certificate.j_scale_min_checked_hypergeometric_inverse_z_power,
            Some(3)
        );
        assert_eq!(
            certificate.j_scale_min_checked_after_outer_z_inverse_power,
            Some(2)
        );
        assert_eq!(
            certificate.j_scale_positive_degree_f_or_g_correction_status,
            "weighted_p2_rank_three_ccit_j_scale_checked_positive_degrees_cannot_modify_F_or_G_by_zero_order_bound"
        );
        assert_eq!(
            certificate.twisted_j_normalization_status,
            "weighted_p2_rank_three_ccit_j_function_normalization_trivial_to_checked_positive_degrees"
        );
        assert_eq!(
            certificate.mirror_map_inversion_status,
            "weighted_p2_rank_three_ccit_mirror_map_has_no_checked_positive_degree_F_or_G_correction"
        );
        assert_eq!(
            certificate.z_minus_one_extraction_status,
            "weighted_p2_rank_three_ccit_checked_positive_degrees_have_no_z_minus_one_primary_terms_first_nonzero_is_descendant_layer"
        );
        assert_eq!(
            certificate.twisted_dual_pairing_status,
            "weighted_p2_rank_three_ccit_untwisted_p2_dual_pairing_computed_but_big_j_pairing_reconstruction_missing"
        );
        assert_eq!(
            certificate.pairing_or_residue_status,
            "weighted_p2_rank_three_source_certificate_pairing_normalization_computed_for_untwisted_p2_but_descendant_reconstruction_missing"
        );
        assert_eq!(
            certificate.qn_history_status,
            "weighted_p2_rank_three_source_certificate_blocked_missing_source_qn_history"
        );
        assert_eq!(
            certificate.cy3_projection_status,
            "weighted_p2_rank_three_visible_phase_is_not_numerical_cy3_requires_source_codim2_or_insertion_history"
        );
        assert_eq!(
            certificate.promotion_status,
            "weighted_p2_rank_three_source_certificate_not_promotable_without_pairing_chamber_qn_history"
        );
        assert_eq!(certificate.checked_max_degree_twice, 4);
        assert_eq!(certificate.checked_integer_sector_count, 2);
        assert_eq!(certificate.checked_half_sector_count, 2);
        assert_eq!(
            certificate.required_inputs,
            vec![
                "source_derived_codim2_insertion_or_equivalent_observable".to_string(),
                "orbifold_sector_pairing_data".to_string(),
                "twisted_vector_bundle_ifunction_chamber_certificate".to_string(),
                "twisted_vector_bundle_ifunction_qn_history".to_string(),
                "twisted_big_j_or_pairing_reconstruction_for_descendant_readout".to_string(),
                "source_derived_chamber_qn_history_for_selected_phase".to_string(),
                "cy3_projection_or_codimension_two_local_model".to_string(),
            ]
        );
        assert!(certificate.source_references.contains(
            &"ccit_jscale_requires_Fz_plus_G_expansion_for_twisted_j_identification".to_string()
        ));
        assert!(
            weighted_p2_rank_three_twisted_ifunction_source_certificate_requirements(
                &[1, 1, 2],
                &[4],
                4,
            )
            .is_none()
        );
    }

    #[test]
    fn weighted_p2_rank_three_ifunction_profiles_pin_full_denominator_readout() {
        let profiles = weighted_p2_rank_three_twisted_ifunction_degree_profiles(
            &[1, 1, 2],
            &[1, 1, 2],
            4,
            Some(2),
        );
        let degree_one = profiles
            .iter()
            .find(|profile| profile.degree_twice == 2)
            .expect("degree-one profile");
        assert_eq!(
            degree_one.split_equivariant_full_hypergeometric_denominator_truncated_coefficients,
            Some(vec!["2".to_string(), "10".to_string(), "18".to_string()])
        );
        assert_eq!(
            degree_one
                .split_equivariant_full_hypergeometric_inverse_denominator_truncated_coefficients,
            Some(vec!["1/2".to_string(), "-5/2".to_string(), "8".to_string()])
        );
        assert_eq!(
            degree_one.split_equivariant_full_hypergeometric_dual_basis_p2_lambda_polynomial,
            Some(vec!["-1/4".to_string(), "1/4".to_string()])
        );
        assert_eq!(
            degree_one
                .split_equivariant_full_hypergeometric_dual_basis_p2_primary_z2_lambda_polynomial,
            Some(vec!["0".to_string()])
        );
        assert_eq!(
            degree_one
                .split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_inverse_power,
            Some(3)
        );
        assert_eq!(
            degree_one
                .split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_polynomial,
            Some(vec!["-1/4".to_string()])
        );
        assert_eq!(
            degree_one.split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status,
            "weighted_p2_rank_three_split_scalar_mirror_map_primary_zero_preserved_by_scalar_reparametrization_to_this_degree"
        );

        let degree_two = profiles
            .iter()
            .find(|profile| profile.degree_twice == 4)
            .expect("degree-two profile");
        assert_eq!(
            degree_two.split_equivariant_full_hypergeometric_denominator_truncated_coefficients,
            Some(vec![
                "96".to_string(),
                "688".to_string(),
                "2072".to_string(),
            ])
        );
        assert_eq!(
            degree_two
                .split_equivariant_full_hypergeometric_inverse_denominator_truncated_coefficients,
            Some(vec![
                "1/96".to_string(),
                "-43/576".to_string(),
                "67/216".to_string(),
            ])
        );
        assert_eq!(
            degree_two.split_equivariant_full_hypergeometric_dual_basis_p2_lambda_polynomial,
            Some(vec![
                "-1/32".to_string(),
                "23/192".to_string(),
                "-17/96".to_string(),
                "1/8".to_string(),
                "-1/24".to_string(),
                "1/192".to_string(),
            ])
        );
        assert_eq!(
            degree_two
                .split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_lambda_polynomial,
            Some(vec!["-1/32".to_string()])
        );
    }

    #[test]
    fn canonical_weighted_p2_ifunction_profile_has_primary_signal() {
        let profile =
            weighted_p2_rank_three_twisted_ifunction_degree_profiles(&[1, 1, 2], &[4], 2, Some(2))
                .into_iter()
                .find(|profile| profile.degree_twice == 2)
                .expect("canonical degree-one profile");
        assert_eq!(
            profile
                .split_equivariant_full_hypergeometric_dual_basis_p2_primary_z2_lambda_polynomial,
            Some(vec!["0".to_string(), "11/4".to_string()])
        );
        assert_eq!(
            profile
                .split_equivariant_full_hypergeometric_dual_basis_p2_first_nonzero_z_inverse_power,
            Some(2)
        );
        assert_eq!(
            profile.split_equivariant_full_hypergeometric_dual_basis_p2_scalar_mirror_map_primary_status,
            "weighted_p2_rank_three_split_scalar_mirror_map_primary_already_nonzero_before_reparametrization"
        );
    }

    #[test]
    fn split_bundle_primary_signals_are_zero_in_adjacent_kp112_diagnostic() {
        assert_eq!(
            split_bundle_kp112_mirror_map_primary_p_lambda_coefficients(4),
            vec![
                "0".to_string(),
                "0".to_string(),
                "0".to_string(),
                "0".to_string(),
            ]
        );
        assert_eq!(
            split_bundle_kp112_mirror_map_half_sector_primary_coefficients(4),
            vec![
                "0".to_string(),
                "0".to_string(),
                "0".to_string(),
                "0".to_string(),
            ]
        );
    }

    #[test]
    fn split_bundle_half_sector_signal_starts_as_descendant() {
        let descendants = split_bundle_kp112_mirror_map_half_sector_first_nonzero_descendants(4, 4);
        assert_eq!(
            descendants
                .iter()
                .map(|entry| entry.inverse_z_power.clone())
                .collect::<Vec<_>>(),
            vec![
                Some("3".to_string()),
                Some("3".to_string()),
                Some("3".to_string()),
                Some("3".to_string()),
            ]
        );
        assert_eq!(
            descendants
                .iter()
                .map(|entry| entry.coefficients.clone())
                .collect::<Vec<_>>(),
            vec![
                vec!["2".to_string()],
                vec!["-322/27".to_string()],
                vec!["-11744/375".to_string()],
                vec!["-22221448/25725".to_string()],
            ]
        );
    }

    #[test]
    fn split_bundle_half_sector_descendant_profile_preserves_z_levels() {
        let profiles = split_bundle_kp112_mirror_map_half_sector_descendant_profiles(2, 4);
        assert_eq!(profiles.len(), 2);
        assert_eq!(
            profiles[0]
                .iter()
                .map(|entry| entry.inverse_z_power.as_str())
                .collect::<Vec<_>>(),
            vec!["2", "3", "4"]
        );
        assert_eq!(profiles[0][0].status, "dual_basis_normalized_zero");
        assert_eq!(
            profiles[0][0].dual_basis_coefficients,
            Some(vec!["0".to_string()])
        );
        assert_eq!(profiles[0][1].status, "dual_basis_normalized_nonzero");
        assert_eq!(
            profiles[0][1].dual_basis_coefficients,
            Some(vec!["2".to_string()])
        );
        assert_eq!(
            profiles[1][1].dual_basis_coefficients,
            Some(vec!["-322/27".to_string()])
        );
    }
}
