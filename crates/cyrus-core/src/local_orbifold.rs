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
                ordinary_sector_dual_basis_p2_z_readout(inverse_z_power, &polynomial)
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
        let coefficients = divide_rational_lambda_polynomial_by_two_lambda_if_possible(&polynomial)
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
) -> OrdinarySectorDualBasisP2ZReadout {
    let raw_fun0_lambda_coefficients = dense_rational_lambda_polynomial_to_strings(polynomial);
    let dual_basis_p2_coefficients = divide_rational_lambda_polynomial_by_two_lambda_if_possible(
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
    let constant = polynomial
        .get(&0)
        .cloned()
        .unwrap_or_else(|| Rational::from(0));
    if constant != Rational::from(0) {
        return None;
    }
    let max_lambda_power = polynomial.keys().last().copied().unwrap_or(0);
    let two = Rational::from(Integer::from(2));
    Some(
        (1..=max_lambda_power)
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
                (
                    "3",
                    Some(vec!["0".to_string(), "0".to_string(), "-1/4".to_string()])
                ),
                (
                    "4",
                    Some(vec![
                        "0".to_string(),
                        "0".to_string(),
                        "0".to_string(),
                        "1/4".to_string()
                    ])
                ),
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
            Some(vec!["0".to_string(), "0".to_string(), "-1/32".to_string()])
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
