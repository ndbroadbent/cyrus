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
use cyrus_core::{Intersection, compute_gv_invariants_with_explicit_semigroup};

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
    };

    match report_target_inner(sample, context, run_integer_diamonds, element_limit) {
        Ok((status, term_count, element_count, gv, error)) => TargetReport {
            status,
            integer_term_count: term_count,
            diamond_element_count: element_count,
            gv,
            error,
            ..base
        },
        Err(error) => TargetReport {
            status: "error".to_string(),
            error: Some(error),
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
        Err(_) => {
            return Ok((
                "hkty_panic".to_string(),
                Some(terms.len()),
                Some(elements.len()),
                None,
                Some("explicit-semigroup HKTY panicked".to_string()),
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

fn build_report(
    context: &CorrectedChamberGvContext,
    validated: &ValidatedContext<'_>,
    run_integer_diamonds: bool,
    element_limit: usize,
) -> ContextReport {
    let targets = validated
        .stats
        .sample
        .iter()
        .enumerate()
        .map(|(idx, sample)| {
            report_target(idx, sample, validated, run_integer_diamonds, element_limit)
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
            "[ERROR] usage: mcallister_gv_context --context path [--run-integer-diamonds] [--element-limit N] [--out path]"
        );
        std::process::exit(2);
    };
    let run_integer_diamonds = parse_flag("--run-integer-diamonds");
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
    let report = build_report(&context, &validated, run_integer_diamonds, element_limit);
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
}
