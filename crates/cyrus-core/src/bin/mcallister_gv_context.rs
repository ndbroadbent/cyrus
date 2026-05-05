//! Consume corrected-chamber GV context exports.
//!
//! This is an opt-in diagnostic binary for the JSON written by
//! `mcallister_first_principles --dump-corrected-chamber-gv-context`. It does
//! not read McAllister downstream GV rows. Its job is to validate the
//! CYTools/cygv-shaped context and, when requested, run small explicit
//! semigroup HKTY checks for missing targets whose exact active-generator
//! decomposition is an integer semigroup.

use malachite::Rational as MalachiteRational;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use cyrus_core::types::rational::Rational;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    Intersection, compute_gv_invariants_with_explicit_semigroup,
    compute_gv_invariants_with_provided_generators,
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
    real_cone_decomposition_active_generator_basis_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    real_cone_decomposition_exact_coefficients: Option<Vec<String>>,
    real_cone_decomposition_exact_kind: Option<String>,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
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
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
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
    element_limit: usize,
) -> TargetReport {
    let exact_kind = sample.real_cone_decomposition_exact_kind.clone();
    let active_generator_count = sample
        .real_cone_decomposition_active_generator_basis_nonzero
        .as_ref()
        .map(Vec::len);
    let base = TargetReport {
        index,
        degree: sample.degree,
        ambient_nonzero: sample.ambient_nonzero.clone(),
        basis_nonzero: sample.basis_nonzero.clone(),
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
    element_limit: usize,
) -> ContextReport {
    let targets = validated
        .stats
        .sample
        .iter()
        .enumerate()
        .map(|(idx, sample)| {
            report_target(
                idx,
                sample,
                validated,
                run_integer_diamonds,
                run_active_support_generators,
                element_limit,
            )
        })
        .collect();
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
            "[ERROR] usage: mcallister_gv_context --context path [--run-integer-diamonds] [--run-active-support-generators] [--element-limit N] [--out path]"
        );
        std::process::exit(2);
    };
    let run_integer_diamonds = parse_flag("--run-integer-diamonds");
    let run_active_support_generators = parse_flag("--run-active-support-generators");
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
    fn target_active_support_merges_target_and_active_generator_supports() {
        let sample = MissingGvTargetSample {
            degree: 3,
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
}
