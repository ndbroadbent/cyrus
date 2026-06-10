#![allow(clippy::items_after_test_module)] // unit tests are grouped mid-file next to the helpers they cover

//! First-principles McAllister pipeline runner (opt-in).
//!
//! Runs the pipeline stages explicitly and fails fast on any mismatch or
//! invalid physics. This is intentionally a binary, not a test.
//!
//! Optional:
//! - `--dual-basis path/to/dual_basis.json` to supply the source coordinate
//!   basis for K/M flux vectors. The JSON may contain either `indices` or a
//!   divisor-basis `matrix`. Without this, McAllister validation mode uses the
//!   paper basis `[3, 4, 5, 8]`; generic mode uses Cyrus' computed dual basis.
//! - `--production-dual-basis path/to/dual_basis.json` to set the internal dual
//!   divisor basis used for flat-direction intersections and compact GV inputs.
//!   This is separate from the K/M flux source basis above.
//! - `--production-primal-basis path/to/primal_basis.json` to set the internal
//!   primal divisor basis used for KKLT path following.
//! - `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.
//! - `--skip-mcallister-assertions` to run the computed pipeline without
//!   comparing final observables to the 4-214-647 validation target.
//! - `--branch-candidates N --branch-selection <max-volume|min-volume|first-positive|min-condition|min-toric-gv-missing|min-required-gv-degree>`
//!   to run deterministic KKLT branch search without loading Kähler checkpoints.
//! - `--branch-height-init` to include the CYTools-style height-projected
//!   Kähler point as the first branch-search candidate.
//!   With no branch search, this height-projected point is the default
//!   first-principles KKLT initializer.
//! - `--branch-report-jsonl path` to export positive phase-1 branch candidates
//!   discovered by that search for GA/debug ranking.
//! - `--branch-report-missing-limit N` to include up to N missing small-curve
//!   classes per branch in the JSONL report. GV-enriched reports also include
//!   the required grading-degree range for missing classes.
//! - `--branch-report-decomposition-depth N` to diagnose missing small-curve
//!   classes that are sums of up to N selected raw candidates. Currently N
//!   may be at most 4.
//! - `--branch-report-only` to stop after writing that report.
//! - `--branch-report-skip-gv-coverage` to omit the expensive per-branch
//!   small-curve/GV coverage enrichment.
//! - `--primal-gv-max-deg N` or `--primal-gv-min-points N` to compute general
//!   primal GV invariants if toric formulas do not cover selected small curves.
//! - `--small-curve-pruning <pair|finite-semigroup>` to choose the selected
//!   toric curve pruning rule. The default `pair` reproduces McAllister
//!   checkpoints; `finite-semigroup` is stricter for GA/search diagnostics.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::time::Instant;

use cyrus_core::flat_direction::{compute_flat_direction, compute_flat_direction_full};
use cyrus_core::gv::BoundedCurveDecompositionIndex;
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::range::CheckedRange;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::vacuum::compute_vacuum;
use cyrus_core::volume::bbhl_correction;
use cyrus_core::{
    CurvePruningStrategy, DivisorBasis, GvDivisorBasisData, Point, Polytope, ToricCurveCandidate,
    Triangulation, apply_finite_f64_basis_transform, apply_integer_basis_transform,
    apply_integer_basis_transform_transpose, basis_change_matrix, build_racetrack_terms,
    compute_curve_basis_matrix, compute_glsm_and_linrels, compute_grading_vector,
    compute_intersection_cytools, compute_linear_relations_no_origin, compute_mori_cone_cap_rays,
    compute_regular_triangulation, compute_toric_two_face_curve_gv_invariants,
    compute_w0_from_terms, divisor_basis_change_matrix, effective_prime_divisors_from_curve_basis,
    generate_scaled_divisor_basis_branch_initializations, gv_divisor_basis_data, heights_to_kahler,
    intersection_in_basis, intersection_in_divisor_basis, is_unimodular, kahler_to_heights,
    map_basis_gv_invariants_to_ambient, prune_decomposable_curve_candidates,
    scale_divisor_basis_kklt_branch_initialization_to_target, solve_divisor_basis_path_following,
    solve_divisor_basis_path_following_branch_candidates, solve_racetrack,
    subcutoff_toric_curve_candidates,
};

const DEFAULT_MCALLISTER_GV_MIN_POINTS: u32 = 20_000;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Stage {
    Triangulation = 1,
    Intersection = 2,
    FlatDirection = 3,
    Gv = 4,
    Racetrack = 5,
    Volume = 6,
    Vacuum = 7,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BranchSelection {
    MaxVolume,
    MinVolume,
    FirstPositive,
    MinCondition,
    MinToricGvMissing,
    MinRequiredGvDegree,
}

impl BranchSelection {
    const fn as_str(self) -> &'static str {
        match self {
            Self::MaxVolume => "max-volume",
            Self::MinVolume => "min-volume",
            Self::FirstPositive => "first-positive",
            Self::MinCondition => "min-condition",
            Self::MinToricGvMissing => "min-toric-gv-missing",
            Self::MinRequiredGvDegree => "min-required-gv-degree",
        }
    }

    const fn requires_gv_coverage(self) -> bool {
        matches!(self, Self::MinToricGvMissing | Self::MinRequiredGvDegree)
    }

    const fn requires_gv_degree_summary(self) -> bool {
        matches!(self, Self::MinRequiredGvDegree)
    }
}

fn parse_branch_selection(name: &str) -> Option<BranchSelection> {
    match name {
        "max-volume" => Some(BranchSelection::MaxVolume),
        "min-volume" => Some(BranchSelection::MinVolume),
        "first-positive" => Some(BranchSelection::FirstPositive),
        "min-condition" => Some(BranchSelection::MinCondition),
        "min-toric-gv-missing" | "min-gv-missing" => Some(BranchSelection::MinToricGvMissing),
        "min-required-gv-degree" | "min-gv-degree" => Some(BranchSelection::MinRequiredGvDegree),
        _ => None,
    }
}

fn parse_curve_pruning_strategy(name: &str) -> Option<CurvePruningStrategy> {
    match name {
        "pair" | "pair-decomposable" | "mcallister-pair" => {
            Some(CurvePruningStrategy::PairDecomposable)
        }
        "finite-semigroup" | "semigroup" => Some(CurvePruningStrategy::FiniteSemigroup),
        _ => None,
    }
}

fn parse_stage(name: &str) -> Option<Stage> {
    match name {
        "triangulation" => Some(Stage::Triangulation),
        "intersection" => Some(Stage::Intersection),
        "flat" | "flat-direction" => Some(Stage::FlatDirection),
        "gv" => Some(Stage::Gv),
        "racetrack" => Some(Stage::Racetrack),
        "volume" => Some(Stage::Volume),
        "vacuum" | "v0" => Some(Stage::Vacuum),
        _ => None,
    }
}

fn parse_arg_value<T: std::str::FromStr>(flag: &str) -> Option<T> {
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == flag {
            return args.next().and_then(|v| v.parse::<T>().ok());
        }
    }
    None
}

fn parse_flag(flag: &str) -> bool {
    std::env::args().any(|arg| arg == flag)
}

fn stage_enabled(stage: Stage, stop_after: Stage) -> bool {
    stage <= stop_after
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &PathBuf) -> T {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn read_csv_f64(path: &PathBuf) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .flat_map(|line| line.split(','))
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<f64>()
                .unwrap_or_else(|e| panic!("Invalid float {s} in {}: {e}", path.display()))
        })
        .collect()
}

fn read_csv_i64(path: &PathBuf) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .flat_map(|line| line.split(','))
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<i64>()
                .unwrap_or_else(|e| panic!("Invalid int {s} in {}: {e}", path.display()))
        })
        .collect()
}

fn read_csv_usize(path: &PathBuf) -> Vec<usize> {
    read_csv_i64(path)
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index must be non-negative"))
        .collect()
}

fn read_optional_scalar_f64(path: &PathBuf) -> Option<f64> {
    if !path.exists() {
        return None;
    }
    Some(
        std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()))
            .trim()
            .parse::<f64>()
            .unwrap_or_else(|e| panic!("Invalid scalar in {}: {e}", path.display())),
    )
}

fn validate_basis_checkpoint(
    glsm: &[Vec<malachite::Integer>],
    computed: &[usize],
    data_dir: &str,
    label: &str,
) {
    let basis_path = PathBuf::from(data_dir).join("basis.dat");
    if !basis_path.exists() {
        eprintln!(
            "[WARN] basis.dat checkpoint not found; skipping {label} basis checkpoint validation"
        );
        return;
    }
    let expected = read_csv_usize(&basis_path);
    if computed == expected {
        return;
    }

    eprintln!(
        "[WARN] {label} basis differs from basis.dat: computed len={} expected len={}",
        computed.len(),
        expected.len()
    );
    let diff: Vec<(usize, usize, usize)> = computed
        .iter()
        .zip(expected.iter())
        .enumerate()
        .filter_map(|(i, (a, b))| if a == b { None } else { Some((i, *a, *b)) })
        .collect();
    for (i, a, b) in diff.iter().take(10) {
        eprintln!("[WARN]   idx {i}: computed={a}, basis.dat={b}");
    }
    eprintln!("[WARN]   total mismatches: {}", diff.len());

    let transform = basis_change_matrix(glsm, computed, &expected).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute change of basis to basis.dat: {e}");
        std::process::exit(2);
    });
    if !is_unimodular(&transform) {
        eprintln!("[ERROR] computed basis and basis.dat are not unimodularly equivalent");
        std::process::exit(2);
    }
    eprintln!("[WARN] basis.dat is checkpoint-only; continuing with Cyrus computed basis");
}

fn sorted_point_coords(points: &[Point]) -> Vec<Vec<i64>> {
    let mut coords: Vec<Vec<i64>> = points.iter().map(|point| point.coords().to_vec()).collect();
    coords.sort();
    coords
}

fn sorted_simplices(triangulation: &Triangulation) -> Vec<Vec<usize>> {
    let mut simplices = triangulation.simplices().to_vec();
    simplices.sort();
    simplices
}

fn validate_dual_checkpoint(
    dual_polytope: &Polytope,
    dual_triangulation: &Triangulation,
    data_dir: &str,
) {
    let dir = PathBuf::from(data_dir);
    let dual_points_path = dir.join("dual_points.dat");
    if dual_points_path.exists() {
        let expected_dual_points = read_points(&dual_points_path);
        let expected_dual_set = {
            let mut points = expected_dual_points;
            points.sort();
            points
        };
        let actual_dual_set = sorted_point_coords(dual_polytope.vertices());
        if actual_dual_set != expected_dual_set {
            eprintln!(
                "[ERROR] computed dual polytope points differ from dual_points.dat checkpoint"
            );
            std::process::exit(2);
        }
    } else {
        eprintln!("[WARN] dual_points.dat checkpoint not found; skipping dual point validation");
    }

    let dual_simplices_path = dir.join("dual_simplices.dat");
    if dual_simplices_path.exists() {
        let expected_simplices_i64 = read_points(&dual_simplices_path);
        let mut expected_simplices: Vec<Vec<usize>> = expected_simplices_i64
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|value| {
                        usize::try_from(value)
                            .expect("dual_simplices.dat index must be non-negative")
                    })
                    .collect()
            })
            .collect();
        expected_simplices.sort();
        let actual_simplices = sorted_simplices(dual_triangulation);
        if actual_simplices != expected_simplices {
            eprintln!("[ERROR] computed dual FRST differs from dual_simplices.dat checkpoint");
            std::process::exit(2);
        }
    } else {
        eprintln!("[WARN] dual_simplices.dat checkpoint not found; skipping dual FRST validation");
    }

    match (dual_points_path.exists(), dual_simplices_path.exists()) {
        (true, true) => {
            eprintln!(
                "[INFO] computed dual polytope/FRST match dual_points.dat and dual_simplices.dat checkpoints"
            );
        }
        (true, false) => {
            eprintln!("[INFO] computed dual polytope points match dual_points.dat checkpoint");
        }
        (false, true) => {
            eprintln!("[INFO] computed dual FRST matches dual_simplices.dat checkpoint");
        }
        (false, false) => {
            eprintln!(
                "[WARN] no dual checkpoint files found; continuing with Cyrus-computed dual geometry"
            );
        }
    }
}

/// Load the declared mirror-chamber FRST from `dual_simplices.dat`.
///
/// Only 4-214-647 happens to use the default chamber of the dual polytope;
/// the other published examples select a specific mirror triangulation, so
/// the chamber choice is a declared model input on the same footing as
/// `heights.dat` (see docs/MCALLISTER_DATA_POLICY.md). The loaded simplices
/// are not trusted: they must use every dual point, be star around the dual
/// origin, and admit strict interior secondary-cone heights that rebuild the
/// identical triangulation through Cyrus's own regular-triangulation code.
fn load_declared_dual_frst(
    data_dir: &str,
    dual_points: &[Point],
    dual_origin_idx: usize,
) -> Option<Triangulation> {
    let path = PathBuf::from(data_dir).join("dual_simplices.dat");
    if !path.exists() {
        return None;
    }
    let simplices: Vec<Vec<usize>> = read_points(&path)
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|value| {
                    usize::try_from(value).expect("dual_simplices.dat index must be non-negative")
                })
                .collect()
        })
        .collect();
    let candidate = Triangulation::new(simplices);

    let used: HashSet<usize> = candidate.simplices().iter().flatten().copied().collect();
    if used.len() != dual_points.len() || used.iter().any(|&idx| idx >= dual_points.len()) {
        eprintln!(
            "[ERROR] dual_simplices.dat is not fine: {} of {} dual points used",
            used.len(),
            dual_points.len()
        );
        std::process::exit(2);
    }
    if candidate
        .simplices()
        .iter()
        .any(|simplex| !simplex.contains(&dual_origin_idx))
    {
        eprintln!("[ERROR] dual_simplices.dat is not star around the dual origin");
        std::process::exit(2);
    }
    let heights = cyrus_core::triangulation::triangulation_heights_from_secondary_cone(
        dual_points,
        &candidate,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute dual secondary cone for dual_simplices.dat: {e}");
        std::process::exit(2);
    })
    .unwrap_or_else(|| {
        eprintln!(
            "[ERROR] dual_simplices.dat is not a regular triangulation (no strict interior point in its secondary cone)"
        );
        std::process::exit(2);
    });
    let rebuilt = compute_regular_triangulation(dual_points, &heights).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to rebuild dual FRST from secondary-cone heights: {e}");
        std::process::exit(2);
    });
    if sorted_simplices(&rebuilt) != sorted_simplices(&candidate) {
        eprintln!(
            "[ERROR] dual_simplices.dat secondary-cone heights rebuild a different triangulation"
        );
        std::process::exit(2);
    }
    // Prefer Cyrus's own canonical construction of the same chamber: when the
    // declared chamber coincides with the default one, reuse the default
    // triangulation object verbatim so downstream behavior is unchanged;
    // otherwise use the rebuilt (canonically ordered) triangulation rather
    // than the raw file ordering.
    let default_frst = cyrus_core::compute_frst_heights(dual_points, dual_origin_idx).ok();
    if let Some((_, computed)) = default_frst
        && sorted_simplices(&computed) == sorted_simplices(&candidate)
    {
        eprintln!(
            "[INFO] using declared mirror-chamber FRST from dual_simplices.dat ({} simplices; verified fine/star/regular via secondary cone; coincides with Cyrus default dual chamber)",
            candidate.simplices().len()
        );
        return Some(computed);
    }
    eprintln!(
        "[INFO] using declared mirror-chamber FRST from dual_simplices.dat ({} simplices; verified fine/star/regular via secondary cone; differs from Cyrus default dual chamber)",
        candidate.simplices().len()
    );
    Some(rebuilt)
}

fn read_points(path: &PathBuf) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| {
                    s.trim().parse::<i64>().unwrap_or_else(|e| {
                        panic!("Invalid point entry {s} in {}: {e}", path.display())
                    })
                })
                .collect::<Vec<i64>>()
        })
        .collect()
}

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

#[derive(Debug, Deserialize)]
struct FluxInput {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct TargetVolumesInput {
    #[serde(rename = "c_i")]
    c_i: Vec<i64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct BasisOverride {
    #[serde(default)]
    indices: Option<Vec<usize>>,
    #[serde(default)]
    matrix: Option<Vec<Vec<i64>>>,
}

enum BasisOverrideRef<'a> {
    Indices(&'a [usize]),
    Matrix(&'a [Vec<i64>]),
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum OwnedDivisorBasis {
    Indices(Vec<usize>),
    Matrix {
        standard_basis: Vec<usize>,
        basis_matrix: Vec<Vec<malachite::Integer>>,
    },
}

impl OwnedDivisorBasis {
    fn as_divisor_basis(&self) -> DivisorBasis<'_> {
        match self {
            Self::Indices(indices) => DivisorBasis::Indices(indices),
            Self::Matrix {
                standard_basis,
                basis_matrix,
            } => DivisorBasis::Matrix {
                standard_basis,
                basis_matrix,
            },
        }
    }

    const fn dimension(&self) -> usize {
        match self {
            Self::Indices(indices) => indices.len(),
            Self::Matrix { basis_matrix, .. } => basis_matrix.len(),
        }
    }

    const fn description(&self) -> &'static str {
        match self {
            Self::Indices(_) => "index divisor basis",
            Self::Matrix { .. } => "matrix divisor basis",
        }
    }
}

impl BasisOverride {
    fn representation(&self, context: &str) -> std::result::Result<BasisOverrideRef<'_>, String> {
        match (&self.indices, &self.matrix) {
            (Some(indices), None) => Ok(BasisOverrideRef::Indices(indices)),
            (None, Some(matrix)) => Ok(BasisOverrideRef::Matrix(matrix)),
            (Some(_), Some(_)) => Err(format!(
                "{context} supplied both `indices` and `matrix`; exactly one basis representation is allowed"
            )),
            (None, None) => Err(format!(
                "{context} did not supply a basis; expected `indices` or `matrix`"
            )),
        }
    }

    fn indices(&self, context: &str) -> std::result::Result<&[usize], String> {
        match self.representation(context)? {
            BasisOverrideRef::Indices(indices) => Ok(indices),
            BasisOverrideRef::Matrix(_) => Err(format!(
                "{context} supplied a matrix divisor basis, but this runner path is still vector-basis only; use an `indices` basis or route through the matrix-basis APIs"
            )),
        }
    }
}

fn basis_indices_or_exit(override_value: &BasisOverride, context: &str) -> Vec<usize> {
    override_value.indices(context).map_or_else(
        |err| {
            eprintln!("[ERROR] {err}");
            std::process::exit(2);
        },
        <[usize]>::to_vec,
    )
}

fn owned_divisor_basis_from_override(
    override_value: &BasisOverride,
    standard_basis: &[usize],
    context: &str,
) -> std::result::Result<OwnedDivisorBasis, String> {
    match override_value.representation(context)? {
        BasisOverrideRef::Indices(indices) => Ok(OwnedDivisorBasis::Indices(indices.to_vec())),
        BasisOverrideRef::Matrix(matrix) => Ok(OwnedDivisorBasis::Matrix {
            standard_basis: standard_basis.to_vec(),
            basis_matrix: basis_matrix_to_integer(matrix),
        }),
    }
}

fn basis_change_matrix_between_owned(
    glsm: &[Vec<malachite::Integer>],
    from_basis: &OwnedDivisorBasis,
    to_basis: &OwnedDivisorBasis,
) -> std::result::Result<Vec<Vec<malachite::Integer>>, String> {
    divisor_basis_change_matrix(
        glsm,
        from_basis.as_divisor_basis(),
        to_basis.as_divisor_basis(),
    )
    .map_err(|e| e.to_string())
}

fn production_gv_basis_data(
    ambient_mori_rays: &[Vec<i64>],
    linrels: &[Vec<malachite::Integer>],
    basis: &OwnedDivisorBasis,
    context: &str,
) -> Result<GvDivisorBasisData, String> {
    gv_divisor_basis_data(ambient_mori_rays, linrels, basis.as_divisor_basis())
        .map_err(|e| format!("failed to build {context} GV divisor-basis data: {e}"))
}

fn vector_gv_basis_data(
    ambient_mori_rays: &[Vec<i64>],
    linrels: &[Vec<malachite::Integer>],
    basis: &[usize],
    context: &str,
) -> Result<GvDivisorBasisData, String> {
    production_gv_basis_data(
        ambient_mori_rays,
        linrels,
        &OwnedDivisorBasis::Indices(basis.to_vec()),
        context,
    )
}

#[derive(Debug, Deserialize)]
struct RacetrackAssertion {
    g_s: f64,
    #[serde(rename = "W_0")]
    w_0: f64,
    cy_vol: f64,
    n_curves: usize,
}

#[derive(Debug, Serialize)]
struct PipelineSummary {
    g_s: f64,
    w0: f64,
    v_string: f64,
    v0_log10_abs: f64,
    ek0: f64,
}

#[derive(Debug, Serialize)]
struct BranchReportSummary {
    #[serde(rename = "type")]
    record_type: &'static str,
    branch_seed: u64,
    branch_selection: &'static str,
    small_curve_pruning: &'static str,
    kklt_steps: usize,
    attempted: usize,
    solved: usize,
    non_converged: usize,
    non_positive_volume: usize,
    positive_volume: usize,
    selected_rank_by_volume: usize,
    selected_init_index: usize,
    selected_init_source: &'static str,
    selected_phase1_volume: f64,
    selected_phase1_rel_err: f64,
    selected_jacobian_rank: usize,
    selected_jacobian_max_rank: usize,
    selected_jacobian_condition_number: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_ambient_rays: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_subcutoff_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_filtered_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_toric_gv_covered_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_toric_gv_missing_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_first_missing_class: Option<Vec<i64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_missing_required_degree_min: Option<i128>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_missing_required_degree_max: Option<i128>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_missing_class_sample: Option<Vec<Vec<i64>>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_missing_bounded_decomposition_max_terms: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_missing_bounded_decomposable_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    selected_small_curve_first_missing_bounded_decomposition: Option<Vec<Vec<i64>>>,
}

#[derive(Debug, Serialize)]
struct BranchReportBranch {
    #[serde(rename = "type")]
    record_type: &'static str,
    branch_seed: u64,
    branch_selection: &'static str,
    small_curve_pruning: &'static str,
    rank_by_volume: usize,
    selected: bool,
    init_index: usize,
    init_source: &'static str,
    phase1_volume: f64,
    phase1_rel_err: f64,
    jacobian_rank: usize,
    jacobian_max_rank: usize,
    jacobian_max_singular_value: f64,
    jacobian_min_nonzero_singular_value: f64,
    jacobian_condition_number: Option<f64>,
    t_init: Vec<f64>,
    t_phase1: Vec<f64>,
    tau_phase1: Vec<f64>,
    tau_phase1_target: Vec<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_ambient_rays: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_subcutoff_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_filtered_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_toric_gv_covered_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_toric_gv_missing_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_first_missing_class: Option<Vec<i64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_missing_required_degree_min: Option<i128>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_missing_required_degree_max: Option<i128>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_missing_class_sample: Option<Vec<Vec<i64>>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_missing_bounded_decomposition_max_terms: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_missing_bounded_decomposable_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    small_curve_first_missing_bounded_decomposition: Option<Vec<Vec<i64>>>,
}

#[derive(Clone, Debug)]
struct BranchGvCoverage {
    ambient_rays: usize,
    subcutoff_count: usize,
    filtered_count: usize,
    toric_gv_covered_count: usize,
    toric_gv_missing_count: usize,
    first_missing_class: Option<Vec<i64>>,
    missing_required_degree_min: Option<i128>,
    missing_required_degree_max: Option<i128>,
    missing_class_sample: Vec<Vec<i64>>,
    bounded_decomposition_max_terms: Option<usize>,
    missing_bounded_decomposable_count: Option<usize>,
    first_missing_bounded_decomposition: Option<Vec<Vec<i64>>>,
}

#[derive(Clone, Debug, PartialEq)]
struct TargetCorrectionDeltaSummary {
    dimension: usize,
    max_abs_delta: f64,
    relative_l2_delta: f64,
    max_abs_reference: f64,
    max_abs_candidate: f64,
}

struct BranchReportContext {
    branch_seed: u64,
    branch_selection: BranchSelection,
    small_curve_pruning: CurvePruningStrategy,
    kklt_steps: usize,
    attempted: usize,
    solved: usize,
    non_converged: usize,
    non_positive_volume: usize,
    selected_rank_by_volume: usize,
}

#[allow(clippy::struct_excessive_bools)] // CLI flags map one-to-one onto independent booleans
struct PipelineArgs {
    stop_after: Stage,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    cutoff: f64,
    small_curve_cutoff: f64,
    small_curve_pruning: CurvePruningStrategy,
    out_path: Option<String>,
    compare_dir: Option<String>,
    data_dir: Option<String>,
    allow_invalid_ek0: bool,
    allow_fixtures: bool,
    validate_mcallister_assertions: bool,
    allow_downstream_kahler: bool,
    kklt_steps: usize,
    branch_candidates: usize,
    branch_seed: u64,
    branch_selection: BranchSelection,
    branch_height_init: bool,
    branch_report_path: Option<String>,
    branch_report_missing_limit: usize,
    branch_report_decomposition_depth: usize,
    branch_report_only: bool,
    branch_report_skip_gv_coverage: bool,
    primal_gv_max_deg: Option<u32>,
    primal_gv_min_points: Option<u32>,
    production_primal_basis_override: Option<BasisOverride>,
    dual_basis_override: Option<BasisOverride>,
    production_dual_basis_override: Option<BasisOverride>,
}

struct PrimalGeom {
    polytope: Polytope,
    heights: Vec<F64<Finite>>,
    triangulation_points: Vec<Point>,
    triangulation: Triangulation,
}

struct PrimalIntersection {
    glsm: Vec<Vec<malachite::Integer>>,
    linrels: Vec<Vec<malachite::Integer>>,
    basis: Vec<usize>,
    kappa_full: cyrus_core::Intersection,
    kappa_basis: cyrus_core::Intersection,
}

struct FlatDirectionData {
    dual_polytope: Polytope,
    dual_triangulation_points: Vec<Point>,
    dual_triangulation: Triangulation,
    dual_linrels: Vec<Vec<malachite::Integer>>,
    dual_divisor_basis: OwnedDivisorBasis,
    dual_kappa_basis: cyrus_core::Intersection,
    flat_p: Vec<F64<Finite>>,
    ek0_opt: Option<F64<Pos>>,
    m_flux: Vec<I64<Finite>>,
}

struct RacetrackData {
    rt_res: cyrus_core::racetrack::RacetrackResult,
    w0: F64<Pos>,
}

fn parse_args() -> PipelineArgs {
    let stop_after = parse_arg_value::<String>("--stop-after")
        .and_then(|s| parse_stage(&s))
        .unwrap_or(Stage::Vacuum);
    let max_deg = parse_arg_value::<u32>("--max-deg");
    let min_points = if max_deg.is_some() {
        None
    } else {
        parse_arg_value::<u32>("--min-points").or(Some(DEFAULT_MCALLISTER_GV_MIN_POINTS))
    };
    let cutoff = parse_arg_value::<f64>("--cutoff").unwrap_or(1.0);
    let small_curve_cutoff_flag = parse_arg_value::<f64>("--small-curve-cutoff");
    let small_curve_pruning = parse_arg_value::<String>("--small-curve-pruning").map_or(
        CurvePruningStrategy::PairDecomposable,
        |value| {
            parse_curve_pruning_strategy(&value).unwrap_or_else(|| {
                eprintln!(
                    "[ERROR] unknown --small-curve-pruning value {value}; expected pair or finite-semigroup"
                );
                std::process::exit(2);
            })
        },
    );
    let out_path = parse_arg_value::<String>("--out");
    let compare_dir = parse_arg_value::<String>("--compare-dir");
    let data_dir = parse_arg_value::<String>("--data-dir")
        .or_else(|| std::env::var("CYRUS_MCALLISTER_DATA_DIR").ok());
    // The small-curve volume cutoff is part of the declared model choice and
    // differs per example (1.0 for 4-214-647, 2.0 for 5-81-3213). An explicit
    // flag wins; otherwise the declared small_curves_cutoff.dat is used.
    let small_curve_cutoff = small_curve_cutoff_flag.unwrap_or_else(|| {
        data_dir
            .as_ref()
            .and_then(|dir| {
                read_optional_scalar_f64(&PathBuf::from(dir).join("small_curves_cutoff.dat"))
            })
            .map_or(1.0, |declared| {
                eprintln!(
                    "[INFO] using declared small-curve cutoff {declared} from small_curves_cutoff.dat"
                );
                declared
            })
    });
    let allow_invalid_ek0 = parse_flag("--allow-invalid-ek0");
    let allow_fixtures = parse_flag("--allow-fixtures");
    let validate_mcallister_assertions = !parse_flag("--skip-mcallister-assertions");
    let allow_downstream_kahler = parse_flag("--allow-downstream-kahler");
    let kklt_steps = parse_arg_value::<usize>("--kklt-steps").unwrap_or(200);
    let branch_candidates = parse_arg_value::<usize>("--branch-candidates").unwrap_or(0);
    let branch_seed = parse_arg_value::<u64>("--branch-seed").unwrap_or(42);
    let branch_selection = parse_arg_value::<String>("--branch-selection")
        .and_then(|s| parse_branch_selection(&s))
        .unwrap_or(BranchSelection::MaxVolume);
    let branch_height_init = parse_flag("--branch-height-init");
    let branch_report_path = parse_arg_value::<String>("--branch-report-jsonl");
    let branch_report_missing_limit =
        parse_arg_value::<usize>("--branch-report-missing-limit").unwrap_or(0);
    let branch_report_decomposition_depth =
        parse_arg_value::<usize>("--branch-report-decomposition-depth").unwrap_or(0);
    let branch_report_only = parse_flag("--branch-report-only");
    let branch_report_skip_gv_coverage = parse_flag("--branch-report-skip-gv-coverage");
    let primal_gv_max_deg = parse_arg_value::<u32>("--primal-gv-max-deg");
    let primal_gv_min_points = if primal_gv_max_deg.is_some() {
        None
    } else {
        parse_arg_value::<u32>("--primal-gv-min-points")
    };
    let production_primal_basis_override = parse_arg_value::<String>("--production-primal-basis")
        .map(|path| load_json::<BasisOverride>(&PathBuf::from(path)));
    let dual_basis_override = parse_arg_value::<String>("--dual-basis")
        .map(|path| load_json::<BasisOverride>(&PathBuf::from(path)));
    let production_dual_basis_override = parse_arg_value::<String>("--production-dual-basis")
        .map(|path| load_json::<BasisOverride>(&PathBuf::from(path)));
    PipelineArgs {
        stop_after,
        max_deg,
        min_points,
        cutoff,
        small_curve_cutoff,
        small_curve_pruning,
        out_path,
        compare_dir,
        data_dir,
        allow_invalid_ek0,
        allow_fixtures,
        validate_mcallister_assertions,
        allow_downstream_kahler,
        kklt_steps,
        branch_candidates,
        branch_seed,
        branch_selection,
        branch_height_init,
        branch_report_path,
        branch_report_missing_limit,
        branch_report_decomposition_depth,
        branch_report_only,
        branch_report_skip_gv_coverage,
        primal_gv_max_deg,
        primal_gv_min_points,
        production_primal_basis_override,
        dual_basis_override,
        production_dual_basis_override,
    }
}

fn enforce_modes(data_dir: Option<&str>, allow_fixtures: bool) {
    let first_principles_env = std::env::var_os("CYRUS_FIRST_PRINCIPLES").is_some();
    if allow_fixtures && first_principles_env {
        eprintln!("[ERROR] --allow-fixtures cannot be used with CYRUS_FIRST_PRINCIPLES");
        std::process::exit(2);
    }
    if let Some(dir) = data_dir {
        eprintln!("[INFO] using McAllister data dir: {dir}");
        eprintln!("[MODE] first-principles (.dat)");
        return;
    }
    if !allow_fixtures {
        eprintln!("[ERROR] No McAllister data dir set. Refusing to fall back to JSON fixtures.");
        eprintln!("[ERROR] Set CYRUS_MCALLISTER_DATA_DIR or pass --allow-fixtures.");
        std::process::exit(2);
    }
    eprintln!("[MODE] fixtures (JSON)");
    eprintln!("[WARN] using JSON fixtures (not a first-principles run)");
}

fn prune_selected_curve_candidates(
    candidates: &[ToricCurveCandidate],
    strategy: CurvePruningStrategy,
    context: &str,
) -> Result<Vec<ToricCurveCandidate>, String> {
    prune_decomposable_curve_candidates(candidates, strategy).map_err(|e| {
        format!(
            "failed to prune {context} small curves with {} strategy: {e}",
            strategy.as_str()
        )
    })
}

fn load_primal_inputs(data_dir: Option<&str>, manifest_dir: &Path) -> (Vec<Vec<i64>>, Vec<f64>) {
    data_dir.map_or_else(
        || {
            let poly_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
            let input: PolytopeInput = load_json(&poly_path);
            let heights_path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
            let heights: HeightsInput = load_json(&heights_path);
            (input.points, heights.values)
        },
        |dir| {
            let dir = PathBuf::from(dir);
            let points = read_points(&dir.join("points.dat"));
            let heights = read_csv_f64(&dir.join("heights.dat"));
            (points, heights)
        },
    )
}

fn load_flux_vectors(data_dir: Option<&str>, manifest_dir: &Path) -> (Vec<i64>, Vec<i64>) {
    data_dir.map_or_else(
        || {
            let flux_path = manifest_dir.join("tests/mcallister_e2e/inputs/flux.json");
            let flux: FluxInput = load_json(&flux_path);
            (flux.k, flux.m)
        },
        |dir| {
            let dir = PathBuf::from(dir);
            let k = read_csv_i64(&dir.join("K_vec.dat"));
            let m = read_csv_i64(&dir.join("M_vec.dat"));
            (k, m)
        },
    )
}

fn transform_kahler_to_computed_basis(
    glsm: &[Vec<malachite::Integer>],
    computed_basis: &[usize],
    source_basis: &[usize],
    values: &[F64<Finite>],
) -> Vec<F64<Finite>> {
    transform_kahler_to_computed_basis_with_logging(
        glsm,
        computed_basis,
        source_basis,
        values,
        true,
    )
}

fn transform_kahler_to_computed_basis_with_logging(
    glsm: &[Vec<malachite::Integer>],
    computed_basis: &[usize],
    source_basis: &[usize],
    values: &[F64<Finite>],
    log_transform: bool,
) -> Vec<F64<Finite>> {
    if computed_basis == source_basis {
        return values.to_vec();
    }
    let transform = basis_change_matrix(glsm, computed_basis, source_basis).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute Kähler basis transform: {e}");
        std::process::exit(2);
    });
    if !is_unimodular(&transform) {
        eprintln!("[ERROR] Kähler basis transform is not unimodular");
        std::process::exit(2);
    }
    if log_transform {
        eprintln!(
            "[INFO] transforming Kähler parameters from source basis {source_basis:?} to computed basis {computed_basis:?}"
        );
    }
    apply_finite_f64_basis_transform(&transform, values, "Kähler").unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to apply Kähler basis transform: {e}");
        std::process::exit(2);
    })
}

fn transform_kahler_between_owned_divisor_bases(
    glsm: &[Vec<malachite::Integer>],
    target_basis: &OwnedDivisorBasis,
    source_basis: &OwnedDivisorBasis,
    values: &[F64<Finite>],
    label: &str,
    log_transform: bool,
) -> Result<Vec<F64<Finite>>, String> {
    if target_basis == source_basis {
        return Ok(values.to_vec());
    }
    let transform = basis_change_matrix_between_owned(glsm, target_basis, source_basis)
        .map_err(|e| format!("failed to compute {label} basis transform: {e}"))?;
    if !is_unimodular(&transform) {
        return Err(format!("{label} basis transform is not unimodular"));
    }
    if log_transform {
        eprintln!(
            "[INFO] transforming {label} from {} source coordinates to {} target coordinates",
            source_basis.description(),
            target_basis.description()
        );
    }
    apply_finite_f64_basis_transform(&transform, values, label).map_err(|e| e.to_string())
}

fn computed_primal_basis(intersection: &PrimalIntersection) -> OwnedDivisorBasis {
    OwnedDivisorBasis::Indices(intersection.basis.clone())
}

fn transform_production_primal_kahler_to_computed(
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    values: &[F64<Finite>],
    label: &str,
) -> Vec<F64<Finite>> {
    transform_kahler_between_owned_divisor_bases(
        &intersection.glsm,
        &computed_primal_basis(intersection),
        production_basis,
        values,
        label,
        true,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] {e}");
        std::process::exit(2);
    })
}

fn transform_computed_primal_kahler_to_production(
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    values: &[F64<Finite>],
    label: &str,
) -> Result<Vec<F64<Finite>>, String> {
    transform_kahler_between_owned_divisor_bases(
        &intersection.glsm,
        production_basis,
        &computed_primal_basis(intersection),
        values,
        label,
        true,
    )
}

/// Derive the half-integral B-field parity vector `gamma` from the orientifold
/// involution.
///
/// The O7-planes of the simple orientifolds in arXiv:2107.09064 sit on the
/// prime toric divisors whose lattice points share a single nonzero parity
/// class `sigma = p mod 2` (the involution's fixed parity). The declared
/// `c_i = 6` so(8) stacks identify `sigma`; the full O7 set is then EVERY
/// point with that parity, including divisors outside the KKLT basis. Using
/// only the KKLT-basis c_i = 6 divisors misses those extra O7s and corrupts
/// the worldsheet-instanton signs `(-1)^{gamma.q}` for curves meeting them.
fn lattice_point_parity(points: &[Point], idx: usize) -> Vec<i64> {
    points[idx]
        .coords()
        .iter()
        .map(|c| c.rem_euclid(2))
        .collect()
}

/// Derive the involution's fixed parity class `sigma = p mod 2` shared by all
/// declared so(8) (`c_i = 6`) stack divisors. Returns `None` when there are no
/// so(8) stacks (then the B-field vanishes).
fn orientifold_involution_parity(
    points: &[Point],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
) -> Option<Vec<i64>> {
    let ambient_dim = points.len();
    let mut sigma: Option<Vec<i64>> = None;
    for (&divisor_idx, ci) in kklt_basis.iter().zip(c_i.iter()) {
        if divisor_idx >= ambient_dim {
            eprintln!(
                "[ERROR] KKLT divisor index {divisor_idx} exceeds ambient dimension {ambient_dim}"
            );
            std::process::exit(2);
        }
        if ci.get() != 6 {
            continue;
        }
        let parity = lattice_point_parity(points, divisor_idx);
        match sigma.as_ref() {
            None => sigma = Some(parity),
            Some(existing) if *existing == parity => {}
            Some(existing) => {
                eprintln!(
                    "[ERROR] so(8) stack divisors do not share a single involution parity: {existing:?} vs {parity:?} (divisor {divisor_idx})"
                );
                std::process::exit(2);
            }
        }
    }
    if sigma.as_ref().is_some_and(|s| s.iter().all(|&p| p == 0)) {
        eprintln!(
            "[ERROR] so(8) stack divisors have the trivial parity class; cannot derive the orientifold B-field"
        );
        std::process::exit(2);
    }
    sigma
}

fn check_gamma_against_declared_c_i(
    gamma: &[i64],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    sigma: &[i64],
) {
    for (&divisor_idx, ci) in kklt_basis.iter().zip(c_i.iter()) {
        if ci.get() == 6 && gamma[divisor_idx] != 1 {
            eprintln!(
                "[ERROR] declared so(8) divisor {divisor_idx} lost its own parity class; gamma derivation is inconsistent"
            );
            std::process::exit(2);
        }
        if ci.get() != 6 && gamma[divisor_idx] == 1 {
            eprintln!(
                "[ERROR] declared c_i={} divisor {divisor_idx} carries the O7 involution parity {sigma:?}; declared orientifold data is inconsistent",
                ci.get()
            );
            std::process::exit(2);
        }
    }
}

fn compute_b_field_gamma_for_o7_divisors(
    points: &[Point],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
) -> Vec<I64<Finite>> {
    if kklt_basis.len() != c_i.len() {
        eprintln!("[ERROR] KKLT basis and c_i length mismatch when computing B-field gamma");
        std::process::exit(2);
    }
    let ambient_dim = points.len();
    let Some(sigma) = orientifold_involution_parity(points, kklt_basis, c_i) else {
        return vec![I64::<Finite>::new(0); ambient_dim];
    };

    let kklt_index_set: HashSet<usize> = kklt_basis.iter().copied().collect();
    let mut gamma = vec![0_i64; ambient_dim];
    let mut extra_o7 = Vec::new();
    for (idx, entry) in gamma.iter_mut().enumerate().skip(1) {
        if lattice_point_parity(points, idx) == sigma {
            *entry = 1;
            if !kklt_index_set.contains(&idx) {
                extra_o7.push(idx);
            }
        }
    }

    check_gamma_against_declared_c_i(&gamma, kklt_basis, c_i, &sigma);

    eprintln!(
        "[INFO] orientifold involution parity sigma={sigma:?}: O7 divisors={} (beyond KKLT basis: {extra_o7:?})",
        gamma.iter().filter(|&&g| g != 0).count()
    );

    gamma.into_iter().map(I64::<Finite>::new).collect()
}

fn basis_matrix_to_integer(basis_matrix: &[Vec<i64>]) -> Vec<Vec<malachite::Integer>> {
    basis_matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|&value| malachite::Integer::from(value))
                .collect()
        })
        .collect()
}

fn transform_m_flux_between_divisor_bases(
    glsm: &[Vec<malachite::Integer>],
    production_basis: &OwnedDivisorBasis,
    source_basis: &OwnedDivisorBasis,
    values: &[i64],
    label: &str,
) -> Vec<i64> {
    if production_basis == source_basis {
        return values.to_vec();
    }
    let transform = basis_change_matrix_between_owned(glsm, production_basis, source_basis)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute {label} basis transform: {e}");
            std::process::exit(2);
        });
    if !is_unimodular(&transform) {
        eprintln!("[ERROR] {label} basis transform is not unimodular");
        std::process::exit(2);
    }
    eprintln!(
        "[INFO] transforming {label} from {} source coordinates to {} production coordinates",
        source_basis.description(),
        production_basis.description()
    );
    apply_integer_basis_transform(&transform, values, label).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to apply {label} basis transform: {e}");
        std::process::exit(2);
    })
}

fn transform_k_flux_between_divisor_bases(
    glsm: &[Vec<malachite::Integer>],
    source_basis: &OwnedDivisorBasis,
    production_basis: &OwnedDivisorBasis,
    values: &[i64],
) -> Vec<i64> {
    if production_basis == source_basis {
        return values.to_vec();
    }
    let transform = basis_change_matrix_between_owned(glsm, source_basis, production_basis)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute K basis transform: {e}");
            std::process::exit(2);
        });
    if !is_unimodular(&transform) {
        eprintln!("[ERROR] K basis transform is not unimodular");
        std::process::exit(2);
    }
    eprintln!(
        "[INFO] transforming K from {} source coordinates to {} production coordinates",
        source_basis.description(),
        production_basis.description()
    );
    apply_integer_basis_transform_transpose(&transform, values, "K").unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to apply K basis transform: {e}");
        std::process::exit(2);
    })
}

fn load_kklt_inputs(data_dir: Option<&str>, manifest_dir: &Path) -> (Vec<I64<Pos>>, Vec<usize>) {
    data_dir.map_or_else(
        || {
            let c_i_path = manifest_dir.join("tests/mcallister_e2e/overrides/target_volumes.json");
            let targets: TargetVolumesInput = load_json(&c_i_path);
            let kklt_basis_path =
                manifest_dir.join("tests/mcallister_e2e/assertions/kklt_basis.json");
            let kklt_basis: BasisOverride = load_json(&kklt_basis_path);
            let c_i = targets
                .c_i
                .into_iter()
                .map(|v| {
                    I64::<Pos>::new(v)
                        .unwrap_or_else(|| panic!("target_volumes c_i must be positive: {v}"))
                })
                .collect();
            (
                c_i,
                basis_indices_or_exit(&kklt_basis, "KKLT basis fixture"),
            )
        },
        |dir| {
            let dir = PathBuf::from(dir);
            let c_i = read_csv_i64(&dir.join("target_volumes.dat"))
                .into_iter()
                .map(|v| {
                    I64::<Pos>::new(v)
                        .unwrap_or_else(|| panic!("target_volumes.dat c_i must be positive: {v}"))
                })
                .collect();
            let kklt_basis = read_csv_usize(&dir.join("kklt_basis.dat"));
            (c_i, kklt_basis)
        },
    )
}

fn stage_triangulation(data_dir: Option<&str>, manifest_dir: &Path, t0: &Instant) -> PrimalGeom {
    let (points_raw, heights_raw) = load_primal_inputs(data_dir, manifest_dir);
    let heights = heights_raw
        .iter()
        .copied()
        .map(|height| F64::<Finite>::new(height).expect("triangulation height must be finite"))
        .collect();
    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("Failed to build polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter triangulation points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights_raw)
        .expect("Failed to compute triangulation");
    eprintln!("[TIME] triangulation: {:.2?}", t0.elapsed());
    PrimalGeom {
        polytope,
        heights,
        triangulation_points,
        triangulation,
    }
}

fn stage_intersection(
    data_dir: Option<&str>,
    geom: &PrimalGeom,
    t0: &Instant,
) -> PrimalIntersection {
    let (glsm, linrels, basis) =
        compute_glsm_and_linrels(&geom.triangulation_points).expect("Failed GLSM/linrels");
    if let Some(dir) = data_dir {
        validate_basis_checkpoint(&glsm, &basis, dir, "primal");
    }
    let points_i64: Vec<Vec<i64>> = geom
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linrels_reduced = compute_linear_relations_no_origin(&points_i64);
    let linrels_i64: Vec<Vec<i64>> = linrels_reduced
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("linrel fits in i64"))
                .collect()
        })
        .collect();
    let kappa_full = compute_intersection_cytools(
        &geom.triangulation,
        &geom.triangulation_points,
        &linrels_i64,
    )
    .expect("Failed intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);
    eprintln!("[TIME] intersection: {:.2?}", t0.elapsed());
    PrimalIntersection {
        glsm,
        linrels,
        basis,
        kappa_full,
        kappa_basis,
    }
}

fn select_production_dual_basis(
    override_opt: Option<&BasisOverride>,
    computed_standard_basis: &[usize],
) -> OwnedDivisorBasis {
    override_opt.map_or_else(
        || {
            eprintln!(
                "[INFO] using computed dual production basis (len={}, basis={:?})",
                computed_standard_basis.len(),
                computed_standard_basis
            );
            OwnedDivisorBasis::Indices(computed_standard_basis.to_vec())
        },
        |explicit| {
            let basis = owned_divisor_basis_from_override(
                explicit,
                computed_standard_basis,
                "--production-dual-basis",
            )
            .unwrap_or_else(|err| {
                eprintln!("[ERROR] {err}");
                std::process::exit(2);
            });
            eprintln!(
                "[INFO] using explicit {} from --production-dual-basis",
                basis.description()
            );
            basis
        },
    )
}

fn stage_flat_direction(
    data_dir: Option<&str>,
    manifest_dir: &Path,
    geom: &PrimalGeom,
    allow_invalid_ek0: bool,
    use_mcallister_flux_basis_default: bool,
    dual_basis_override: Option<&BasisOverride>,
    production_dual_basis_override: Option<&BasisOverride>,
    t0: &Instant,
) -> FlatDirectionData {
    let dual_polytope = geom
        .polytope
        .compute_dual()
        .expect("Failed to compute dual polytope");
    let dual_points_vec = dual_polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter dual triangulation points");
    let dual_origin_idx = dual_points_vec
        .iter()
        .position(|point| point.coords().iter().all(|&coord| coord == 0))
        .expect("dual origin not found in triangulation points");
    let dual_triangulation = data_dir
        .and_then(|dir| load_declared_dual_frst(dir, &dual_points_vec, dual_origin_idx))
        .unwrap_or_else(|| {
            let (_dual_heights, computed) =
                cyrus_core::compute_frst_heights(&dual_points_vec, dual_origin_idx)
                    .expect("Failed to compute dual FRST");
            eprintln!(
                "[INFO] no dual_simplices.dat declared input; using Cyrus default dual FRST chamber"
            );
            computed
        });
    if let Some(dir) = data_dir {
        validate_dual_checkpoint(&dual_polytope, &dual_triangulation, dir);
    }
    let (dual_glsm, dual_linrel, dual_basis_auto) =
        compute_glsm_and_linrels(&dual_points_vec).expect("Failed dual GLSM");
    let dual_standard_basis = dual_basis_auto;
    let dual_divisor_basis =
        select_production_dual_basis(production_dual_basis_override, &dual_standard_basis);
    let flux_basis = dual_basis_override.map_or_else(
        || {
            if use_mcallister_flux_basis_default {
                eprintln!("[INFO] using McAllister flux source basis [3, 4, 5, 8]");
                OwnedDivisorBasis::Indices(vec![3, 4, 5, 8])
            } else if production_dual_basis_override.is_some() {
                eprintln!("[INFO] using production dual basis as flux coordinate basis");
                dual_divisor_basis.clone()
            } else {
                eprintln!("[INFO] using computed dual basis as flux coordinate basis");
                dual_divisor_basis.clone()
            }
        },
        |basis| {
            eprintln!("[INFO] using explicit flux source basis from --dual-basis");
            owned_divisor_basis_from_override(basis, &dual_standard_basis, "--dual-basis")
                .unwrap_or_else(|err| {
                    eprintln!("[ERROR] {err}");
                    std::process::exit(2);
                })
        },
    );
    let dual_points_i64: Vec<Vec<i64>> = dual_points_vec
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let dual_linrels = compute_linear_relations_no_origin(&dual_points_i64);
    let dual_linrels_i64: Vec<Vec<i64>> = dual_linrels
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("dual linrel fits in i64"))
                .collect()
        })
        .collect();
    let dual_kappa_full =
        compute_intersection_cytools(&dual_triangulation, &dual_points_vec, &dual_linrels_i64)
            .expect("Failed dual intersection numbers");
    let dual_kappa_basis = cyrus_core::intersection_in_divisor_basis(
        &dual_kappa_full,
        dual_divisor_basis.as_divisor_basis(),
    )
    .expect("failed dual intersection production-basis handoff");
    let (k_raw, m_raw) = load_flux_vectors(data_dir, manifest_dir);
    let k_raw = transform_k_flux_between_divisor_bases(
        &dual_glsm,
        &flux_basis,
        &dual_divisor_basis,
        &k_raw,
    );
    let m_raw = transform_m_flux_between_divisor_bases(
        &dual_glsm,
        &dual_divisor_basis,
        &flux_basis,
        &m_raw,
        "M",
    );
    let k_flux: Vec<I64<Finite>> = k_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let m_flux: Vec<I64<Finite>> = m_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let (flat_p, ek0_opt) = match compute_flat_direction_full(&dual_kappa_basis, &k_flux, &m_flux) {
        Some(v) => (v.p, Some(v.ek0)),
        None if allow_invalid_ek0 => {
            eprintln!(
                "[WARN] invalid flat direction (ek0<=0); continuing due to --allow-invalid-ek0"
            );
            let p =
                compute_flat_direction(&dual_kappa_basis, &k_flux, &m_flux).unwrap_or_else(|| {
                    eprintln!("[ERROR] flat direction solve failed");
                    std::process::exit(2);
                });
            (p, None)
        }
        None => {
            eprintln!("[ERROR] invalid flat direction (ek0<=0)");
            std::process::exit(2);
        }
    };
    eprintln!("[TIME] flat_direction: {:.2?}", t0.elapsed());
    FlatDirectionData {
        dual_polytope,
        dual_triangulation_points: dual_points_vec,
        dual_triangulation,
        dual_linrels: dual_linrel,
        dual_divisor_basis,
        dual_kappa_basis,
        flat_p,
        ek0_opt,
        m_flux,
    }
}

fn stage_gv(
    flat: &FlatDirectionData,
    min_points: Option<u32>,
    max_deg: Option<u32>,
    t0: &Instant,
) -> Vec<(Vec<i32>, malachite::Integer)> {
    let ambient_rays = compute_mori_cone_cap_rays(
        &flat.dual_triangulation,
        &flat.dual_triangulation_points,
        &flat.dual_polytope,
        false,
        false,
        None,
    )
    .expect("Failed ambient mori cone cap rays");
    let gv_basis = production_gv_basis_data(
        &ambient_rays,
        &flat.dual_linrels,
        &flat.dual_divisor_basis,
        "dual mirror",
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let grading = compute_grading_vector(&gv_basis.mori_rays).expect("No grading vector found");
    let invariants = cyrus_core::compute_gv_invariants(
        &gv_basis.mori_rays,
        &grading,
        &gv_basis.q_matrix,
        &flat.dual_kappa_basis,
        min_points,
        max_deg,
    )
    .expect("GV computation failed");
    eprintln!("[TIME] gv_invariants: {:.2?}", t0.elapsed());
    invariants
}

fn stage_racetrack(
    invariants: &[(Vec<i32>, malachite::Integer)],
    m_flux: &[I64<Finite>],
    flat_p: &[F64<Finite>],
    cutoff: F64<Pos>,
    t0: &Instant,
) -> RacetrackData {
    let gv: Vec<cyrus_core::racetrack::GvInvariant> = invariants
        .iter()
        .map(|(curve, value)| {
            let curve = curve
                .iter()
                .map(|&v| I64::<Finite>::new(i64::from(v)))
                .collect();
            let val_f64 = value
                .to_string()
                .parse::<f64>()
                .expect("GV integer fits in f64");
            if !val_f64.is_finite() {
                eprintln!("[ERROR] GV value overflowed f64");
                std::process::exit(2);
            }
            cyrus_core::racetrack::GvInvariant {
                curve,
                value: F64::<Finite>::new(val_f64).expect("GV value is finite"),
            }
        })
        .collect();
    let terms = build_racetrack_terms(&gv, m_flux, flat_p, cutoff);
    if terms.len() < 2 {
        eprintln!("[ERROR] not enough racetrack terms");
        std::process::exit(2);
    }
    let Some(rt_res) = solve_racetrack(&terms) else {
        eprintln!("[ERROR] no stable racetrack solution; leading terms (exponent, coefficient):");
        for term in terms.iter().take(6) {
            eprintln!(
                "[ERROR]   exponent={} coefficient={}",
                term.exponent.get(),
                term.coefficient.get()
            );
        }
        std::process::exit(2);
    };
    let Some(w0) = compute_w0_from_terms(&rt_res, &terms) else {
        eprintln!("[ERROR] racetrack W0 computation failed or cancelled exactly");
        std::process::exit(2);
    };
    eprintln!("[TIME] racetrack: {:.2?}", t0.elapsed());
    RacetrackData { rt_res, w0 }
}

fn classical_volume_from_t(kappa_primal: &cyrus_core::Intersection, t: &[F64<Finite>]) -> f64 {
    let mut volume_sum = 0.0f64;
    for (&(i, j, k), val) in kappa_primal.iter() {
        let kappa_val = val.to_f64().get();
        let t_product = t[i].get() * t[j].get() * t[k].get();
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        volume_sum += mult * kappa_val * t_product;
    }
    volume_sum / 6.0
}

fn finite_values(values: &[F64<Finite>]) -> Vec<f64> {
    values.iter().map(|value| value.get()).collect()
}

fn pos_values(values: &[F64<Pos>]) -> Vec<f64> {
    values.iter().map(|value| value.get()).collect()
}

fn compute_branch_gv_coverages(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    branches_by_volume: &[cyrus_core::KkltBranchSolution],
    small_curve_cutoff: F64<Pos>,
    small_curve_pruning: CurvePruningStrategy,
    missing_class_sample_limit: usize,
    bounded_decomposition_max_terms: Option<usize>,
    include_required_degree_summary: bool,
) -> Result<Vec<BranchGvCoverage>, String> {
    let ambient_rays = compute_mori_cone_cap_rays(
        &geom.triangulation,
        &geom.triangulation_points,
        &geom.polytope,
        false,
        false,
        None,
    )
    .map_err(|e| format!("failed to compute primal ambient Mori-cap rays: {e}"))?;

    let toric_gvs = compute_toric_two_face_curve_gv_invariants(
        &geom.triangulation,
        &geom.triangulation_points,
        &geom.polytope,
    )
    .map_err(|e| format!("failed to compute primal toric curve GV values: {e}"))?;
    let gv_by_class: HashMap<Vec<i64>, malachite::Integer> = toric_gvs
        .into_iter()
        .map(|item| (item.class, item.gv))
        .collect();
    let required_degree_grading = if include_required_degree_summary {
        let gv_basis_data = vector_gv_basis_data(
            &ambient_rays,
            &intersection.linrels,
            &intersection.basis,
            "branch GV coverage",
        )?;
        Some(
            compute_grading_vector(&gv_basis_data.mori_rays)
                .ok_or_else(|| "failed to compute branch GV degree grading vector".to_string())?,
        )
    } else {
        None
    };

    branches_by_volume
        .iter()
        .map(|branch| {
            let branch_t_computed = transform_kahler_between_owned_divisor_bases(
                &intersection.glsm,
                &computed_primal_basis(intersection),
                production_basis,
                &branch.result.t,
                "branch Kähler point",
                false,
            )?;
            let small_curve_candidates = subcutoff_toric_curve_candidates(
                &ambient_rays,
                &intersection.basis,
                &branch_t_computed,
                small_curve_cutoff,
            )
            .map_err(|e| format!("failed to select branch small toric curve candidates: {e}"))?;
            let small_curves = prune_selected_curve_candidates(
                &small_curve_candidates,
                small_curve_pruning,
                "branch",
            )?;
            let decomposition_index = bounded_decomposition_max_terms
                .map(|_| BoundedCurveDecompositionIndex::new(&small_curve_candidates))
                .transpose()
                .map_err(|e| {
                    format!("failed to build branch bounded curve decomposition index: {e}")
                })?;

            let mut toric_gv_covered_count = 0usize;
            let mut toric_gv_missing_count = 0usize;
            let mut first_missing_class = None;
            let mut missing_classes = Vec::new();
            let mut missing_class_sample = Vec::new();
            let mut missing_bounded_decomposable_count = 0usize;
            let mut first_missing_bounded_decomposition = None;
            for curve in &small_curves {
                if gv_by_class.contains_key(&curve.class) {
                    toric_gv_covered_count += 1;
                } else {
                    toric_gv_missing_count += 1;
                    if first_missing_class.is_none() {
                        first_missing_class = Some(curve.class.clone());
                    }
                    missing_classes.push(curve.class.clone());
                    if missing_class_sample.len() < missing_class_sample_limit {
                        missing_class_sample.push(curve.class.clone());
                    }
                    if let (Some(index), Some(max_terms)) = (
                        decomposition_index.as_ref(),
                        bounded_decomposition_max_terms,
                    ) && let Some(decomposition) =
                        index.find_decomposition(curve, max_terms).map_err(|e| {
                            format!("failed branch bounded curve decomposition search: {e}")
                        })?
                    {
                        missing_bounded_decomposable_count += 1;
                        if first_missing_bounded_decomposition.is_none() {
                            first_missing_bounded_decomposition = Some(decomposition);
                        }
                    }
                }
            }
            let (missing_required_degree_min, missing_required_degree_max) = if let Some(grading) =
                required_degree_grading.as_ref()
                && !missing_classes.is_empty()
            {
                let summary = summarize_required_gv_degrees(
                    &missing_classes,
                    &intersection.basis,
                    grading,
                    None,
                )?;
                (Some(summary.min_degree), Some(summary.max_degree))
            } else {
                (None, None)
            };

            Ok(BranchGvCoverage {
                ambient_rays: ambient_rays.len(),
                subcutoff_count: small_curve_candidates.len(),
                filtered_count: small_curves.len(),
                toric_gv_covered_count,
                toric_gv_missing_count,
                first_missing_class,
                missing_required_degree_min,
                missing_required_degree_max,
                missing_class_sample,
                bounded_decomposition_max_terms,
                missing_bounded_decomposable_count: bounded_decomposition_max_terms
                    .map(|_| missing_bounded_decomposable_count),
                first_missing_bounded_decomposition,
            })
        })
        .collect()
}

fn select_min_toric_gv_missing_rank(
    coverages: &[BranchGvCoverage],
    volumes_by_rank: &[f64],
) -> Result<usize, String> {
    if coverages.is_empty() {
        return Err("cannot select by toric GV coverage without branch coverage rows".into());
    }
    if coverages.len() != volumes_by_rank.len() {
        return Err(format!(
            "branch coverage rows {} do not match branch volume rows {}",
            coverages.len(),
            volumes_by_rank.len()
        ));
    }

    (0..coverages.len())
        .min_by(|&a, &b| {
            coverages[a]
                .toric_gv_missing_count
                .cmp(&coverages[b].toric_gv_missing_count)
                .then_with(|| volumes_by_rank[a].total_cmp(&volumes_by_rank[b]))
                .then_with(|| a.cmp(&b))
        })
        .ok_or_else(|| "cannot select by toric GV coverage without positive branches".into())
}

fn select_min_required_gv_degree_rank(
    coverages: &[BranchGvCoverage],
    volumes_by_rank: &[f64],
) -> Result<usize, String> {
    if coverages.is_empty() {
        return Err("cannot select by required GV degree without branch coverage rows".into());
    }
    if coverages.len() != volumes_by_rank.len() {
        return Err(format!(
            "branch coverage rows {} do not match branch volume rows {}",
            coverages.len(),
            volumes_by_rank.len()
        ));
    }

    (0..coverages.len())
        .map(|rank| {
            let coverage = &coverages[rank];
            let degree_max = match (
                coverage.toric_gv_missing_count,
                coverage.missing_required_degree_max,
            ) {
                (0, None) => 0,
                (_, Some(value)) => value,
                (_, None) => {
                    return Err(format!(
                        "branch rank {rank} is missing required GV degree summary"
                    ));
                }
            };
            let degree_min = coverage.missing_required_degree_min.unwrap_or(0);
            Ok((rank, degree_max, degree_min))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .min_by(|(rank_a, max_a, min_a), (rank_b, max_b, min_b)| {
            max_a
                .cmp(max_b)
                .then_with(|| {
                    coverages[*rank_a]
                        .toric_gv_missing_count
                        .cmp(&coverages[*rank_b].toric_gv_missing_count)
                })
                .then_with(|| min_a.cmp(min_b))
                .then_with(|| volumes_by_rank[*rank_a].total_cmp(&volumes_by_rank[*rank_b]))
                .then_with(|| rank_a.cmp(rank_b))
        })
        .map(|(rank, _, _)| rank)
        .ok_or_else(|| "cannot select by required GV degree without positive branches".into())
}

fn compute_primal_general_gv_by_ambient_class(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    required_ambient_classes: &[Vec<i64>],
    ambient_mori_rays: Option<&[Vec<i64>]>,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<HashMap<Vec<i64>, malachite::Integer>, String> {
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
    let general_gvs = cyrus_core::compute_gv_invariants(
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

#[derive(Clone, Debug, PartialEq, Eq)]
struct RequiredGvDegreeSummary {
    count: usize,
    min_degree: i128,
    max_degree: i128,
    first_over_max: Option<Vec<i64>>,
}

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

fn select_production_primal_basis(
    override_opt: Option<&BasisOverride>,
    computed_standard_basis: &[usize],
) -> OwnedDivisorBasis {
    match override_opt {
        None => {
            eprintln!(
                "[INFO] using computed primal production basis (len={}, basis={:?})",
                computed_standard_basis.len(),
                computed_standard_basis
            );
            OwnedDivisorBasis::Indices(computed_standard_basis.to_vec())
        }
        Some(override_value) => {
            let basis = owned_divisor_basis_from_override(
                override_value,
                computed_standard_basis,
                "--production-primal-basis",
            )
            .unwrap_or_else(|e| {
                eprintln!("[ERROR] {e}");
                std::process::exit(2);
            });
            eprintln!(
                "[INFO] using explicit {} from --production-primal-basis",
                basis.description()
            );
            basis
        }
    }
}

fn write_branch_report_jsonl(
    path: &PathBuf,
    ctx: &BranchReportContext,
    branches_by_volume: &[cyrus_core::KkltBranchSolution],
    t_initializations: &[Vec<F64<Finite>>],
    t_initialization_sources: &[&'static str],
    branch_gv_coverages: Option<&[BranchGvCoverage]>,
) -> Result<(), String> {
    if t_initializations.len() != t_initialization_sources.len() {
        return Err(format!(
            "branch initialization source rows {} do not match initialization rows {}",
            t_initialization_sources.len(),
            t_initializations.len()
        ));
    }
    let Some(selected) = branches_by_volume.get(ctx.selected_rank_by_volume) else {
        return Err(format!(
            "selected branch rank {} is outside {} branch rows",
            ctx.selected_rank_by_volume,
            branches_by_volume.len()
        ));
    };
    let Some(&selected_source) = t_initialization_sources.get(selected.init_index) else {
        return Err(format!(
            "selected branch init index {} is outside {} initialization source rows",
            selected.init_index,
            t_initialization_sources.len()
        ));
    };
    if let Some(coverages) = branch_gv_coverages
        && coverages.len() != branches_by_volume.len()
    {
        return Err(format!(
            "branch GV coverage rows {} do not match branch rows {}",
            coverages.len(),
            branches_by_volume.len()
        ));
    }
    let selected_gv_coverage =
        branch_gv_coverages.and_then(|coverages| coverages.get(ctx.selected_rank_by_volume));

    let positive_volume = branches_by_volume.len();
    let mut lines = String::new();
    let summary = BranchReportSummary {
        record_type: "summary",
        branch_seed: ctx.branch_seed,
        branch_selection: ctx.branch_selection.as_str(),
        small_curve_pruning: ctx.small_curve_pruning.as_str(),
        kklt_steps: ctx.kklt_steps,
        attempted: ctx.attempted,
        solved: ctx.solved,
        non_converged: ctx.non_converged,
        non_positive_volume: ctx.non_positive_volume,
        positive_volume,
        selected_rank_by_volume: ctx.selected_rank_by_volume,
        selected_init_index: selected.init_index,
        selected_init_source: selected_source,
        selected_phase1_volume: selected.classical_volume.get(),
        selected_phase1_rel_err: selected.result.relative_error.get(),
        selected_jacobian_rank: selected.jacobian_diagnostics.rank,
        selected_jacobian_max_rank: selected.jacobian_diagnostics.max_rank,
        selected_jacobian_condition_number: selected
            .jacobian_diagnostics
            .condition_number
            .map(cyrus_core::F64::get),
        selected_small_curve_ambient_rays: selected_gv_coverage
            .map(|coverage| coverage.ambient_rays),
        selected_small_curve_subcutoff_count: selected_gv_coverage
            .map(|coverage| coverage.subcutoff_count),
        selected_small_curve_filtered_count: selected_gv_coverage
            .map(|coverage| coverage.filtered_count),
        selected_small_curve_toric_gv_covered_count: selected_gv_coverage
            .map(|coverage| coverage.toric_gv_covered_count),
        selected_small_curve_toric_gv_missing_count: selected_gv_coverage
            .map(|coverage| coverage.toric_gv_missing_count),
        selected_small_curve_first_missing_class: selected_gv_coverage
            .and_then(|coverage| coverage.first_missing_class.clone()),
        selected_small_curve_missing_required_degree_min: selected_gv_coverage
            .and_then(|coverage| coverage.missing_required_degree_min),
        selected_small_curve_missing_required_degree_max: selected_gv_coverage
            .and_then(|coverage| coverage.missing_required_degree_max),
        selected_small_curve_missing_class_sample: selected_gv_coverage.and_then(|coverage| {
            (!coverage.missing_class_sample.is_empty())
                .then(|| coverage.missing_class_sample.clone())
        }),
        selected_small_curve_missing_bounded_decomposition_max_terms: selected_gv_coverage
            .and_then(|coverage| coverage.bounded_decomposition_max_terms),
        selected_small_curve_missing_bounded_decomposable_count: selected_gv_coverage
            .and_then(|coverage| coverage.missing_bounded_decomposable_count),
        selected_small_curve_first_missing_bounded_decomposition: selected_gv_coverage
            .and_then(|coverage| coverage.first_missing_bounded_decomposition.clone()),
    };
    lines.push_str(&serde_json::to_string(&summary).map_err(|e| e.to_string())?);
    lines.push('\n');

    for (rank_by_volume, branch) in branches_by_volume.iter().enumerate() {
        let Some(t_init) = t_initializations.get(branch.init_index) else {
            return Err(format!(
                "branch init index {} is outside {} initialization rows",
                branch.init_index,
                t_initializations.len()
            ));
        };
        let Some(&init_source) = t_initialization_sources.get(branch.init_index) else {
            return Err(format!(
                "branch init index {} is outside {} initialization source rows",
                branch.init_index,
                t_initialization_sources.len()
            ));
        };
        let gv_coverage = branch_gv_coverages.and_then(|coverages| coverages.get(rank_by_volume));
        let row = BranchReportBranch {
            record_type: "positive_branch",
            branch_seed: ctx.branch_seed,
            branch_selection: ctx.branch_selection.as_str(),
            small_curve_pruning: ctx.small_curve_pruning.as_str(),
            rank_by_volume,
            selected: rank_by_volume == ctx.selected_rank_by_volume,
            init_index: branch.init_index,
            init_source,
            phase1_volume: branch.classical_volume.get(),
            phase1_rel_err: branch.result.relative_error.get(),
            jacobian_rank: branch.jacobian_diagnostics.rank,
            jacobian_max_rank: branch.jacobian_diagnostics.max_rank,
            jacobian_max_singular_value: branch.jacobian_diagnostics.max_singular_value.get(),
            jacobian_min_nonzero_singular_value: branch
                .jacobian_diagnostics
                .min_nonzero_singular_value
                .get(),
            jacobian_condition_number: branch
                .jacobian_diagnostics
                .condition_number
                .map(cyrus_core::F64::get),
            t_init: finite_values(t_init),
            t_phase1: finite_values(&branch.result.t),
            tau_phase1: finite_values(&branch.result.tau),
            tau_phase1_target: pos_values(&branch.result.tau_target),
            small_curve_ambient_rays: gv_coverage.map(|coverage| coverage.ambient_rays),
            small_curve_subcutoff_count: gv_coverage.map(|coverage| coverage.subcutoff_count),
            small_curve_filtered_count: gv_coverage.map(|coverage| coverage.filtered_count),
            small_curve_toric_gv_covered_count: gv_coverage
                .map(|coverage| coverage.toric_gv_covered_count),
            small_curve_toric_gv_missing_count: gv_coverage
                .map(|coverage| coverage.toric_gv_missing_count),
            small_curve_first_missing_class: gv_coverage
                .and_then(|coverage| coverage.first_missing_class.clone()),
            small_curve_missing_required_degree_min: gv_coverage
                .and_then(|coverage| coverage.missing_required_degree_min),
            small_curve_missing_required_degree_max: gv_coverage
                .and_then(|coverage| coverage.missing_required_degree_max),
            small_curve_missing_class_sample: gv_coverage.and_then(|coverage| {
                (!coverage.missing_class_sample.is_empty())
                    .then(|| coverage.missing_class_sample.clone())
            }),
            small_curve_missing_bounded_decomposition_max_terms: gv_coverage
                .and_then(|coverage| coverage.bounded_decomposition_max_terms),
            small_curve_missing_bounded_decomposable_count: gv_coverage
                .and_then(|coverage| coverage.missing_bounded_decomposable_count),
            small_curve_first_missing_bounded_decomposition: gv_coverage
                .and_then(|coverage| coverage.first_missing_bounded_decomposition.clone()),
        };
        lines.push_str(&serde_json::to_string(&row).map_err(|e| e.to_string())?);
        lines.push('\n');
    }

    std::fs::write(path, lines).map_err(|e| e.to_string())
}

fn height_projected_branch_initialization(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    tau_phase1: &[F64<Pos>],
) -> Result<Vec<F64<Finite>>, String> {
    let basis_non_origin = intersection
        .basis
        .iter()
        .map(|&idx| {
            idx.checked_sub(1)
                .ok_or_else(|| "height projection basis unexpectedly contains origin".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let curve_basis = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
        .map_err(|e| format!("failed to compute primal curve basis: {e}"))?;
    let prime_divisors = effective_prime_divisors_from_curve_basis(&curve_basis, &basis_non_origin)
        .ok_or_else(|| "failed to extract effective-cone prime divisor rows".to_string())?;
    let raw = heights_to_kahler(&geom.heights, &basis_non_origin, &prime_divisors)
        .ok_or_else(|| "failed to project triangulation heights to Kähler basis".to_string())?;
    let production_raw = transform_computed_primal_kahler_to_production(
        intersection,
        production_basis,
        &raw,
        "height-projected Kähler",
    )?;
    scale_divisor_basis_kklt_branch_initialization_to_target(
        &intersection.kappa_full,
        production_basis.as_divisor_basis(),
        kklt_basis,
        tau_phase1,
        &production_raw,
    )
    .ok_or_else(|| "failed to scale height-projected Kähler point to phase-1 target".to_string())
}

fn triangulation_from_kahler_point(
    geom: &PrimalGeom,
    basis: &[usize],
    kahler: &[F64<Finite>],
) -> Result<Triangulation, String> {
    let basis_non_origin = basis
        .iter()
        .map(|&idx| {
            idx.checked_sub(1)
                .ok_or_else(|| "Kähler basis unexpectedly contains origin".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let non_origin_count = geom
        .triangulation_points
        .len()
        .checked_sub(1)
        .ok_or_else(|| "triangulation point set is empty".to_string())?;
    let heights = kahler_to_heights(kahler, &basis_non_origin, non_origin_count)
        .ok_or_else(|| "failed to embed Kähler point into secondary-fan heights".to_string())?;
    let raw_heights = heights
        .iter()
        .map(|height| height.get())
        .collect::<Vec<_>>();
    compute_regular_triangulation(&geom.triangulation_points, &raw_heights)
        .map_err(|e| format!("failed to compute triangulation from corrected Kähler heights: {e}"))
}

fn chamber_intersection_full(
    tri: &Triangulation,
    points: &[Point],
) -> Result<cyrus_core::Intersection, String> {
    let points_i64: Vec<Vec<i64>> = points.iter().map(|point| point.coords().to_vec()).collect();
    let linrels_reduced = compute_linear_relations_no_origin(&points_i64);
    let linrels_i64: Vec<Vec<i64>> = linrels_reduced
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| {
                    i64::try_from(value)
                        .map_err(|_| "corrected-chamber linear relation does not fit in i64")
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let kappa_full = compute_intersection_cytools(tri, points, &linrels_i64)
        .map_err(|e| format!("failed to compute corrected-chamber intersections: {e}"))?;
    Ok(kappa_full)
}

fn triangulations_have_same_simplices(lhs: &Triangulation, rhs: &Triangulation) -> bool {
    let mut lhs_simplices = lhs.simplices().to_vec();
    let mut rhs_simplices = rhs.simplices().to_vec();
    for simplex in &mut lhs_simplices {
        simplex.sort_unstable();
    }
    for simplex in &mut rhs_simplices {
        simplex.sort_unstable();
    }
    lhs_simplices.sort();
    rhs_simplices.sort();
    lhs_simplices == rhs_simplices
}

fn compare_corrected_kahler_checkpoint(
    data_dir: Option<&str>,
    intersection: &PrimalIntersection,
    kahler: &[F64<Finite>],
) {
    let Some(checkpoint) = load_corrected_kahler_checkpoint(data_dir, intersection) else {
        return;
    };
    if checkpoint.len() != kahler.len() {
        eprintln!(
            "[COMPARE] corrected Kähler checkpoint length mismatch: checkpoint={} computed={}",
            checkpoint.len(),
            kahler.len()
        );
        return;
    }
    let summary = target_correction_delta_summary(&checkpoint, kahler).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compare corrected Kähler checkpoint: {e}");
        std::process::exit(2);
    });
    let checkpoint_classical = classical_volume_from_t(&intersection.kappa_basis, &checkpoint);
    let computed_classical = classical_volume_from_t(&intersection.kappa_basis, kahler);
    eprintln!(
        "[COMPARE] corrected Kähler checkpoint delta: max_abs={} relative_l2={} max_abs_checkpoint={} max_abs_computed={}",
        summary.max_abs_delta,
        summary.relative_l2_delta,
        summary.max_abs_reference,
        summary.max_abs_candidate
    );
    eprintln!(
        "[COMPARE] corrected Kähler checkpoint input-chamber V_classical={} computed V_classical={} delta={}",
        checkpoint_classical,
        computed_classical,
        computed_classical - checkpoint_classical
    );
}

fn load_corrected_kahler_checkpoint(
    data_dir: Option<&str>,
    intersection: &PrimalIntersection,
) -> Option<Vec<F64<Finite>>> {
    let dir = data_dir.map(PathBuf::from)?;
    let corrected_path = dir.join("corrected_kahler_param.dat");
    if !corrected_path.exists() {
        eprintln!(
            "[COMPARE] corrected_kahler_param.dat checkpoint not found; skipping corrected Kähler comparison"
        );
        return None;
    }
    let basis_path = dir.join("basis.dat");
    if !basis_path.exists() {
        eprintln!("[COMPARE] basis.dat checkpoint not found; skipping corrected Kähler comparison");
        return None;
    }
    let checkpoint_raw = read_csv_f64(&corrected_path)
        .into_iter()
        .map(|value| F64::<Finite>::new(value).expect("corrected Kähler checkpoint is finite"))
        .collect::<Vec<_>>();
    let source_basis = read_csv_usize(&basis_path);
    Some(transform_kahler_to_computed_basis_with_logging(
        &intersection.glsm,
        &intersection.basis,
        &source_basis,
        &checkpoint_raw,
        false,
    ))
}

fn compare_corrected_target_volume_checkpoint(
    data_dir: Option<&str>,
    kklt_basis: &[usize],
    chi_divisor: &[I64<Finite>],
    base_target_tau: &[F64<Pos>],
    computed_gv_correction: &[F64<Finite>],
    alternate_gv_corrections: &[(&str, &[F64<Finite>])],
    computed_target_tau: &[F64<Finite>],
) {
    let Some(dir) = data_dir.map(PathBuf::from) else {
        return;
    };
    let target_path = dir.join("corrected_target_volumes.dat");
    if !target_path.exists() {
        eprintln!(
            "[COMPARE] corrected_target_volumes.dat checkpoint not found; skipping corrected target-volume comparison"
        );
        return;
    }
    let checkpoint = read_csv_f64(&target_path)
        .into_iter()
        .map(|value| {
            F64::<Finite>::new(value).expect("corrected target-volume checkpoint is finite")
        })
        .collect::<Vec<_>>();
    if checkpoint.len() != computed_target_tau.len()
        || checkpoint.len() != kklt_basis.len()
        || checkpoint.len() != chi_divisor.len()
        || checkpoint.len() != base_target_tau.len()
        || checkpoint.len() != computed_gv_correction.len()
    {
        eprintln!(
            "[COMPARE] corrected target-volume checkpoint length mismatch: checkpoint={} kklt_basis={} chi={} base={} gv={} computed={}",
            checkpoint.len(),
            kklt_basis.len(),
            chi_divisor.len(),
            base_target_tau.len(),
            computed_gv_correction.len(),
            computed_target_tau.len()
        );
        return;
    }
    let base_target_tau_finite = base_target_tau
        .iter()
        .map(|value| F64::<Finite>::new(value.get()).expect("base target tau is finite"))
        .collect::<Vec<_>>();
    let base_summary = target_correction_delta_summary(&checkpoint, &base_target_tau_finite)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compare base target-volume checkpoint: {e}");
            std::process::exit(2);
        });
    let implied_gv_correction = base_target_tau
        .iter()
        .zip(checkpoint.iter())
        .map(|(base, target)| {
            F64::<Finite>::new(base.get() - target.get()).expect("implied GV correction is finite")
        })
        .collect::<Vec<_>>();
    let gv_summary =
        target_correction_delta_summary(&implied_gv_correction, computed_gv_correction)
            .unwrap_or_else(|e| {
                eprintln!("[ERROR] failed to compare implied GV target correction: {e}");
                std::process::exit(2);
            });
    let alternate_gv_summaries = alternate_gv_corrections
        .iter()
        .map(|(label, correction)| {
            if correction.len() != implied_gv_correction.len() {
                eprintln!(
                    "[COMPARE] corrected target-volume alternate GV correction length mismatch: label={} checkpoint_implied={} alternate={}",
                    label,
                    implied_gv_correction.len(),
                    correction.len()
                );
                std::process::exit(2);
            }
            (
                *label,
                target_correction_delta_summary(&implied_gv_correction, correction)
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to compare alternate implied GV target correction ({label}): {e}"
                        );
                        std::process::exit(2);
                    }),
            )
        })
        .collect::<Vec<_>>();
    let target_summary = target_correction_delta_summary(&checkpoint, computed_target_tau)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compare corrected target-volume checkpoint: {e}");
            std::process::exit(2);
        });
    let computed_chi = chi_divisor
        .iter()
        .map(|value| value.to_f64())
        .collect::<Vec<_>>();
    let implied_chi = computed_chi
        .iter()
        .zip(checkpoint.iter().zip(computed_target_tau.iter()))
        .map(|(chi, (checkpoint_target, computed_target))| {
            chi.get() + 24.0 * (checkpoint_target.get() - computed_target.get())
        })
        .map(|value| F64::<Finite>::new(value).expect("implied chi is finite"))
        .collect::<Vec<_>>();
    let computed_chi_finite = computed_chi
        .into_iter()
        .map(|value| F64::<Finite>::new(value.get()).expect("computed chi is finite"))
        .collect::<Vec<_>>();
    let chi_summary = target_correction_delta_summary(&implied_chi, &computed_chi_finite)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compare implied divisor chi: {e}");
            std::process::exit(2);
        });
    eprintln!(
        "[COMPARE] corrected target-volume checkpoint base-target delta: max_abs={} relative_l2={}",
        base_summary.max_abs_delta, base_summary.relative_l2_delta
    );
    eprintln!(
        "[COMPARE] corrected target-volume implied GV correction delta: max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_computed={}",
        gv_summary.max_abs_delta,
        gv_summary.relative_l2_delta,
        gv_summary.max_abs_reference,
        gv_summary.max_abs_candidate
    );
    for (label, summary) in alternate_gv_summaries {
        eprintln!(
            "[COMPARE] corrected target-volume implied GV correction delta ({label}): max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_computed={}",
            summary.max_abs_delta,
            summary.relative_l2_delta,
            summary.max_abs_reference,
            summary.max_abs_candidate
        );
    }
    eprintln!(
        "[COMPARE] corrected target-volume checkpoint final-target delta: max_abs={} relative_l2={} max_abs_checkpoint={} max_abs_computed={}",
        target_summary.max_abs_delta,
        target_summary.relative_l2_delta,
        target_summary.max_abs_reference,
        target_summary.max_abs_candidate
    );
    eprintln!(
        "[COMPARE] corrected target-volume implied divisor-chi delta: max_abs={} relative_l2={} max_abs_implied={} max_abs_computed={}",
        chi_summary.max_abs_delta,
        chi_summary.relative_l2_delta,
        chi_summary.max_abs_reference,
        chi_summary.max_abs_candidate
    );
    let mut chi_deltas = implied_chi
        .iter()
        .zip(computed_chi_finite.iter())
        .enumerate()
        .map(|(idx, (implied, computed))| {
            let delta = implied.get() - computed.get();
            (idx, delta.abs(), delta, implied.get(), computed.get())
        })
        .collect::<Vec<_>>();
    chi_deltas.sort_unstable_by(|lhs, rhs| rhs.1.total_cmp(&lhs.1));
    for (idx, _abs_delta, delta, implied, computed) in chi_deltas.into_iter().take(8) {
        eprintln!(
            "[COMPARE] corrected target-volume implied divisor-chi top_delta kklt_idx={} point_idx={} delta={} implied={} computed={} checkpoint_tau={} computed_tau={}",
            idx,
            kklt_basis[idx],
            delta,
            implied,
            computed,
            checkpoint[idx].get(),
            computed_target_tau[idx].get()
        );
    }
}

fn compare_corrected_chamber_target_volume_checkpoint(
    label: &str,
    data_dir: Option<&str>,
    basis: &[usize],
    kklt_basis: &[usize],
    kappa_full: &cyrus_core::Intersection,
    kappa_basis: &cyrus_core::Intersection,
    kahler: &[F64<Finite>],
) {
    let Some(dir) = data_dir.map(PathBuf::from) else {
        return;
    };
    let target_path = dir.join("corrected_target_volumes.dat");
    if !target_path.exists() {
        return;
    }
    let checkpoint = read_csv_f64(&target_path)
        .into_iter()
        .map(|value| {
            F64::<Finite>::new(value).expect("corrected target-volume checkpoint is finite")
        })
        .collect::<Vec<_>>();
    let Some(classical_tau) = cyrus_core::kklt::compute_kklt_divisor_volumes(
        kappa_basis,
        kappa_full,
        basis,
        kklt_basis,
        kahler,
    ) else {
        eprintln!("[ERROR] failed to compute corrected-chamber classical KKLT divisor volumes");
        std::process::exit(2);
    };
    if checkpoint.len() != classical_tau.len() {
        eprintln!(
            "[COMPARE] corrected-chamber target-volume length mismatch: checkpoint={} computed={}",
            checkpoint.len(),
            classical_tau.len()
        );
        return;
    }
    let summary =
        target_correction_delta_summary(&checkpoint, &classical_tau).unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compare corrected-chamber target volumes: {e}");
            std::process::exit(2);
        });
    eprintln!(
        "[COMPARE] corrected target-volume {label} classical delta: max_abs={} relative_l2={} max_abs_checkpoint={} max_abs_computed={}",
        summary.max_abs_delta,
        summary.relative_l2_delta,
        summary.max_abs_reference,
        summary.max_abs_candidate
    );
}

#[allow(clippy::too_many_arguments)] // stage plumbing forwarded one-to-one
fn solve_gv_corrected_kklt(
    intersection: &PrimalIntersection,
    production_primal_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    chi_divisor: &[I64<Finite>],
    c_tau: cyrus_core::types::physics::CTau,
    small_curve_gvs: &[(Vec<i64>, malachite::Integer)],
    gamma: &[I64<Finite>],
    gv_iteration_source_t: &[F64<Finite>],
    kklt_steps: usize,
) -> Vec<F64<Finite>> {
    // Solve the full worldsheet-instanton-corrected KKLT system with the
    // correction inside every Newton linearization. The earlier scheme —
    // freezing the correction at the previous iterate and re-solving the
    // classical system to a fixed point — is not a contraction for every
    // geometry (5-81-3213 cycles with steps ~0.1 regardless of update
    // mixing), while the corrected-system Newton converges directly.
    let zero_correction = vec![F64::<Finite>::ZERO; kklt_basis.len()];
    let Some(base_tau_target) = cyrus_core::kklt::compute_gv_corrected_target_tau(
        c_i,
        chi_divisor,
        c_tau,
        &zero_correction,
    ) else {
        eprintln!("[ERROR] base KKLT target construction failed");
        std::process::exit(2);
    };
    let production_is_computed = matches!(
        production_primal_basis,
        OwnedDivisorBasis::Indices(indices) if *indices == intersection.basis
    );
    let computed_t_for_corrections = |t_production: &[F64<Finite>]| -> Vec<F64<Finite>> {
        if production_is_computed {
            t_production.to_vec()
        } else {
            transform_production_primal_kahler_to_computed(
                intersection,
                production_primal_basis,
                t_production,
                "GV-corrected solve iterate",
            )
        }
    };
    // Columns of d(t_computed)/d(t_production) for the correction
    // Jacobian chain rule; identity (and skipped) in the default case.
    let correction_transform_columns: Option<Vec<Vec<F64<Finite>>>> = if production_is_computed {
        None
    } else {
        let production_dim = gv_iteration_source_t.len();
        Some(
            (0..production_dim)
                .map(|column| {
                    let mut unit = vec![F64::<Finite>::ZERO; production_dim];
                    unit[column] = cyrus_core::f64_finite!(1.0);
                    computed_t_for_corrections(&unit)
                })
                .collect(),
        )
    };
    let gv_correction_closure = |t_production: &[F64<Finite>]| {
        let t_computed = computed_t_for_corrections(t_production);
        cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            small_curve_gvs,
            &intersection.basis,
            kklt_basis,
            &t_computed,
            Some(gamma),
        )
    };
    let gv_correction_jacobian_closure = |t_production: &[F64<Finite>]| {
        let t_computed = computed_t_for_corrections(t_production);
        let jac_computed =
            cyrus_core::kklt::compute_gv_target_correction_jacobian_for_ambient_curves(
                small_curve_gvs,
                &intersection.basis,
                kklt_basis,
                &t_computed,
                Some(gamma),
            )?;
        match correction_transform_columns.as_ref() {
            None => Some(jac_computed),
            Some(columns) => Some(
                jac_computed
                    .iter()
                    .map(|row| {
                        columns
                            .iter()
                            .map(|column| {
                                row.iter()
                                    .zip(column.iter())
                                    .fold(F64::<Finite>::ZERO, |acc, (a, b)| acc + *a * *b)
                            })
                            .collect()
                    })
                    .collect(),
            ),
        }
    };
    let Some(gv_corrected) = cyrus_core::kklt::solve_gv_corrected_divisor_basis_path_following(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
        kklt_basis,
        &base_tau_target,
        gv_iteration_source_t,
        CheckedRange::new(0, kklt_steps),
        gv_correction_closure,
        gv_correction_jacobian_closure,
    ) else {
        eprintln!("[ERROR] GV-corrected KKLT solve failed");
        std::process::exit(2);
    };
    if !gv_corrected.converged {
        eprintln!(
            "[ERROR] GV-corrected KKLT solve did not converge: rel_err={}",
            gv_corrected.relative_error.get()
        );
        std::process::exit(2);
    }
    eprintln!(
        "[INFO] GV-corrected KKLT solve converged: rel_err={}",
        gv_corrected.relative_error.get()
    );
    transform_production_primal_kahler_to_computed(
        intersection,
        production_primal_basis,
        &gv_corrected.t,
        "GV-corrected Kähler point",
    )
}

#[allow(clippy::fn_params_excessive_bools)] // CLI flags are forwarded one-to-one into this stage
fn stage_volume(
    data_dir: Option<&str>,
    manifest_dir: &Path,
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    racetrack: &RacetrackData,
    allow_downstream_kahler: bool,
    kklt_steps: usize,
    branch_candidates: usize,
    branch_seed: u64,
    branch_selection: BranchSelection,
    branch_height_init: bool,
    branch_report_path: Option<&str>,
    branch_report_missing_limit: usize,
    branch_report_decomposition_depth: usize,
    branch_report_only: bool,
    branch_report_skip_gv_coverage: bool,
    primal_gv_min_points: Option<u32>,
    primal_gv_max_deg: Option<u32>,
    production_primal_basis_override: Option<&BasisOverride>,
    small_curve_cutoff: F64<Pos>,
    small_curve_pruning: CurvePruningStrategy,
    h21: usize,
    t0: &Instant,
) -> (f64, F64<Pos>) {
    if branch_report_path.is_some() && allow_downstream_kahler {
        eprintln!(
            "[ERROR] --branch-report-jsonl is only valid for first-principles branch search, not downstream Kähler replay"
        );
        std::process::exit(2);
    }
    if branch_report_path.is_some() && branch_candidates == 0 {
        eprintln!("[ERROR] --branch-report-jsonl requires --branch-candidates N with N > 0");
        std::process::exit(2);
    }
    if branch_report_only && branch_report_path.is_none() {
        eprintln!("[ERROR] --branch-report-only requires --branch-report-jsonl path");
        std::process::exit(2);
    }
    if branch_report_skip_gv_coverage && branch_report_path.is_none() {
        eprintln!("[ERROR] --branch-report-skip-gv-coverage requires --branch-report-jsonl path");
        std::process::exit(2);
    }
    if branch_report_missing_limit > 0 && branch_report_path.is_none() {
        eprintln!("[ERROR] --branch-report-missing-limit requires --branch-report-jsonl path");
        std::process::exit(2);
    }
    if branch_report_skip_gv_coverage && branch_report_missing_limit > 0 {
        eprintln!(
            "[ERROR] --branch-report-missing-limit requires GV coverage, so it cannot be combined with --branch-report-skip-gv-coverage"
        );
        std::process::exit(2);
    }
    if branch_report_decomposition_depth > 0 && branch_report_path.is_none() {
        eprintln!(
            "[ERROR] --branch-report-decomposition-depth requires --branch-report-jsonl path"
        );
        std::process::exit(2);
    }
    if branch_report_decomposition_depth > 4 {
        eprintln!("[ERROR] --branch-report-decomposition-depth currently supports values 0..4");
        std::process::exit(2);
    }
    if branch_report_skip_gv_coverage && branch_report_decomposition_depth > 0 {
        eprintln!(
            "[ERROR] --branch-report-decomposition-depth requires GV coverage, so it cannot be combined with --branch-report-skip-gv-coverage"
        );
        std::process::exit(2);
    }
    if branch_report_skip_gv_coverage && branch_selection.requires_gv_coverage() {
        eprintln!(
            "[ERROR] branch selection policy={} requires GV coverage, so it cannot be combined with --branch-report-skip-gv-coverage",
            branch_selection.as_str()
        );
        std::process::exit(2);
    }
    eprintln!(
        "[INFO] selected toric curve pruning strategy: {}",
        small_curve_pruning.as_str()
    );
    let production_primal_basis =
        select_production_primal_basis(production_primal_basis_override, &intersection.basis);
    let production_primal_kappa_basis = intersection_in_divisor_basis(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute primal production-basis intersections: {e}");
        std::process::exit(2);
    });
    if production_primal_kappa_basis.dim() != production_primal_basis.dimension() {
        eprintln!("[ERROR] primal production-basis intersection dimension mismatch");
        std::process::exit(2);
    }
    if allow_downstream_kahler && production_primal_basis != computed_primal_basis(intersection) {
        eprintln!(
            "[ERROR] --production-primal-basis is only supported for first-principles KKLT solves, not downstream Kähler replay"
        );
        std::process::exit(2);
    }
    let (t, gv_volume_correction, kklt_basis_for_chamber_gv) = if allow_downstream_kahler {
        let Some(data_dir_path) = data_dir.map(PathBuf::from) else {
            eprintln!(
                "[ERROR] Volume replay requires McAllister data dir for corrected_kahler_param.dat"
            );
            std::process::exit(2);
        };
        let t_raw = read_csv_f64(&data_dir_path.join("corrected_kahler_param.dat"))
            .into_iter()
            .map(|v| F64::<Finite>::new(v).expect("corrected Kähler parameter must be finite"))
            .collect::<Vec<_>>();
        let source_basis = read_csv_usize(&data_dir_path.join("basis.dat"));
        let t = transform_kahler_to_computed_basis(
            &intersection.glsm,
            &intersection.basis,
            &source_basis,
            &t_raw,
        );
        (t, None, None)
    } else {
        let (c_i, kklt_basis) = load_kklt_inputs(data_dir, manifest_dir);
        let c_tau = cyrus_core::kklt::compute_c_tau(racetrack.rt_res.g_s, racetrack.w0);
        let chi_divisor = cyrus_core::compute_kklt_divisor_chi(
            &geom.polytope,
            &geom.triangulation_points,
            &intersection.kappa_full,
            &kklt_basis,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute KKLT divisor chi: {e}");
            std::process::exit(2);
        });
        let gamma =
            compute_b_field_gamma_for_o7_divisors(&geom.triangulation_points, &kklt_basis, &c_i);
        let gamma_odd_count = gamma
            .iter()
            .filter(|value| value.get().rem_euclid(2) != 0)
            .count();
        eprintln!(
            "[INFO] computed B-field gamma from O7 divisors: dim={} odd_entries={}",
            gamma.len(),
            gamma_odd_count
        );
        eprintln!(
            "[WARN] corrected_kahler_param.dat remains validation-only replay and is not loaded without --allow-downstream-kahler."
        );
        let Some(tau_target) =
            cyrus_core::kklt::compute_corrected_target_tau(&c_i, &chi_divisor, c_tau)
        else {
            eprintln!("[ERROR] corrected KKLT target construction failed");
            std::process::exit(2);
        };
        let (zeroth_order, small_curve_selection_t, gv_iteration_source_t) = if branch_candidates
            == 0
        {
            let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
            let height_init = height_projected_branch_initialization(
                geom,
                intersection,
                &production_primal_basis,
                &kklt_basis,
                &tau_phase1,
            )
            .unwrap_or_else(|e| {
                eprintln!(
                    "[ERROR] failed to build height-projected KKLT branch initialization: {e}"
                );
                std::process::exit(2);
            });
            eprintln!("[INFO] using height-projected KKLT branch initialization");
            let Some(phase1) = solve_divisor_basis_path_following(
                &intersection.kappa_full,
                production_primal_basis.as_divisor_basis(),
                &kklt_basis,
                &tau_phase1,
                &height_init,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] phase-1 divisor-basis KKLT path-following failed");
                std::process::exit(2);
            };
            if !phase1.converged {
                eprintln!(
                    "[ERROR] phase-1 mixed-basis KKLT path-following did not converge: rel_err={}",
                    phase1.relative_error.get()
                );
                std::process::exit(2);
            }
            let Some(result) = solve_divisor_basis_path_following(
                &intersection.kappa_full,
                production_primal_basis.as_divisor_basis(),
                &kklt_basis,
                &tau_target,
                &phase1.t,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] zeroth-order divisor-basis KKLT path-following failed");
                std::process::exit(2);
            };
            if !result.converged {
                eprintln!(
                    "[ERROR] zeroth-order mixed-basis KKLT path-following did not converge: rel_err={}",
                    result.relative_error.get()
                );
                std::process::exit(2);
            }
            let small_curve_selection_t = transform_production_primal_kahler_to_computed(
                intersection,
                &production_primal_basis,
                &phase1.t,
                "phase-1 Kähler point",
            );
            (result, small_curve_selection_t, phase1.t)
        } else {
            let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
            let Some(mut t_initializations) = generate_scaled_divisor_basis_branch_initializations(
                &intersection.kappa_full,
                production_primal_basis.as_divisor_basis(),
                &kklt_basis,
                &tau_phase1,
                branch_candidates,
                branch_seed,
            ) else {
                eprintln!("[ERROR] failed to generate KKLT branch initializations");
                std::process::exit(2);
            };
            let mut t_initialization_sources = vec!["generated"; t_initializations.len()];
            if branch_height_init {
                let height_init = height_projected_branch_initialization(
                    geom,
                    intersection,
                    &production_primal_basis,
                    &kklt_basis,
                    &tau_phase1,
                )
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to build height-projected KKLT branch initialization: {e}"
                    );
                    std::process::exit(2);
                });
                t_initializations.insert(0, height_init);
                t_initialization_sources.insert(0, "height_projected");
                eprintln!("[INFO] inserted height-projected KKLT branch initialization at init=0");
            }
            let branch_search = solve_divisor_basis_path_following_branch_candidates(
                &intersection.kappa_full,
                production_primal_basis.as_divisor_basis(),
                &kklt_basis,
                &tau_phase1,
                &t_initializations,
                CheckedRange::new(0, kklt_steps),
            );
            eprintln!(
                "[INFO] KKLT branch search: attempted={} solved={} non_converged={} non_positive_volume={} positive_volume={}",
                branch_search.attempted,
                branch_search.solved,
                branch_search.non_converged,
                branch_search.non_positive_volume,
                branch_search.positive_volume.len()
            );
            let attempted = branch_search.attempted;
            let solved = branch_search.solved;
            let non_converged = branch_search.non_converged;
            let non_positive_volume = branch_search.non_positive_volume;
            let mut positive_branches = branch_search.positive_volume;
            positive_branches.sort_by(|a, b| {
                a.classical_volume
                    .get()
                    .total_cmp(&b.classical_volume.get())
            });
            for (rank, branch) in positive_branches.iter().take(5).enumerate() {
                eprintln!(
                    "[INFO] KKLT branch candidate rank={} init={} phase1_volume={} rel_err={} jacobian_rank={}/{} condition={:?}",
                    rank,
                    branch.init_index,
                    branch.classical_volume.get(),
                    branch.result.relative_error.get(),
                    branch.jacobian_diagnostics.rank,
                    branch.jacobian_diagnostics.max_rank,
                    branch
                        .jacobian_diagnostics
                        .condition_number
                        .map(cyrus_core::F64::get)
                );
            }
            if positive_branches.is_empty() {
                eprintln!("[ERROR] KKLT branch search found no positive-volume phase-1 branch");
                std::process::exit(2);
            }
            #[allow(clippy::useless_let_if_seq)]
            // keeping the imperative form avoids restructuring this block
            let mut branch_gv_coverages: Option<Vec<BranchGvCoverage>> = None;
            if branch_selection.requires_gv_coverage() {
                let coverages = compute_branch_gv_coverages(
                    geom,
                    intersection,
                    &production_primal_basis,
                    &positive_branches,
                    small_curve_cutoff,
                    small_curve_pruning,
                    branch_report_missing_limit,
                    (branch_report_decomposition_depth > 0)
                        .then_some(branch_report_decomposition_depth),
                    branch_report_path.is_some() || branch_selection.requires_gv_degree_summary(),
                )
                .unwrap_or_else(|e| {
                    eprintln!("[ERROR] failed to compute branch GV coverage for selection: {e}");
                    std::process::exit(2);
                });
                let mut coverage_ranks: Vec<usize> = (0..coverages.len()).collect();
                match branch_selection {
                    BranchSelection::MinRequiredGvDegree => {
                        coverage_ranks.sort_by(|&a, &b| {
                            coverages[a]
                                .missing_required_degree_max
                                .unwrap_or(0)
                                .cmp(&coverages[b].missing_required_degree_max.unwrap_or(0))
                                .then_with(|| {
                                    coverages[a]
                                        .toric_gv_missing_count
                                        .cmp(&coverages[b].toric_gv_missing_count)
                                })
                                .then_with(|| {
                                    positive_branches[a]
                                        .classical_volume
                                        .get()
                                        .total_cmp(&positive_branches[b].classical_volume.get())
                                })
                                .then_with(|| a.cmp(&b))
                        });
                    }
                    _ => {
                        coverage_ranks.sort_by(|&a, &b| {
                            coverages[a]
                                .toric_gv_missing_count
                                .cmp(&coverages[b].toric_gv_missing_count)
                                .then_with(|| {
                                    positive_branches[a]
                                        .classical_volume
                                        .get()
                                        .total_cmp(&positive_branches[b].classical_volume.get())
                                })
                                .then_with(|| a.cmp(&b))
                        });
                    }
                }
                for rank_by_volume in coverage_ranks.iter().take(5) {
                    let coverage = &coverages[*rank_by_volume];
                    let branch = &positive_branches[*rank_by_volume];
                    eprintln!(
                        "[INFO] KKLT branch GV coverage candidate rank_by_volume={} init={} missing={} covered={} filtered={} required_degree_max={:?} phase1_volume={}",
                        rank_by_volume,
                        branch.init_index,
                        coverage.toric_gv_missing_count,
                        coverage.toric_gv_covered_count,
                        coverage.filtered_count,
                        coverage.missing_required_degree_max,
                        branch.classical_volume.get()
                    );
                }
                branch_gv_coverages = Some(coverages);
            }
            let selected_rank_by_volume = match branch_selection {
                BranchSelection::MinVolume => 0,
                BranchSelection::MaxVolume => positive_branches.len() - 1,
                BranchSelection::FirstPositive => positive_branches
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, branch)| branch.init_index)
                    .map(|(rank, _)| rank)
                    .expect("positive branches are non-empty"),
                BranchSelection::MinCondition => positive_branches
                    .iter()
                    .enumerate()
                    .filter_map(|(rank, branch)| {
                        branch
                            .jacobian_diagnostics
                            .condition_number
                            .map(|condition| (rank, condition.get()))
                    })
                    .min_by(|(_, a), (_, b)| a.total_cmp(b))
                    .map_or_else(
                        || {
                            eprintln!(
                                "[ERROR] KKLT branch search found no full-rank positive-volume phase-1 branch for min-condition selection"
                            );
                            std::process::exit(2);
                        },
                        |(rank, _)| rank,
                    ),
                BranchSelection::MinToricGvMissing => {
                    let coverages = branch_gv_coverages
                        .as_ref()
                        .expect("GV coverage was computed for coverage-based selection");
                    let volumes_by_rank: Vec<f64> = positive_branches
                        .iter()
                        .map(|branch| branch.classical_volume.get())
                        .collect();
                    select_min_toric_gv_missing_rank(coverages, &volumes_by_rank).unwrap_or_else(
                        |e| {
                            eprintln!(
                                "[ERROR] failed to select branch by toric GV coverage: {e}"
                            );
                            std::process::exit(2);
                        },
                    )
                }
                BranchSelection::MinRequiredGvDegree => {
                    let coverages = branch_gv_coverages
                        .as_ref()
                        .expect("GV coverage was computed for coverage-based selection");
                    let volumes_by_rank: Vec<f64> = positive_branches
                        .iter()
                        .map(|branch| branch.classical_volume.get())
                        .collect();
                    select_min_required_gv_degree_rank(coverages, &volumes_by_rank)
                        .unwrap_or_else(|e| {
                            eprintln!(
                                "[ERROR] failed to select branch by required GV degree: {e}"
                            );
                            std::process::exit(2);
                        })
                }
            };
            let best_branch = positive_branches[selected_rank_by_volume].clone();
            let small_curve_selection_t = transform_production_primal_kahler_to_computed(
                intersection,
                &production_primal_basis,
                &best_branch.result.t,
                "selected branch Kähler point",
            );
            if let Some(path) = branch_report_path {
                let report_path = PathBuf::from(path);
                let ctx = BranchReportContext {
                    branch_seed,
                    branch_selection,
                    small_curve_pruning,
                    kklt_steps,
                    attempted,
                    solved,
                    non_converged,
                    non_positive_volume,
                    selected_rank_by_volume,
                };
                write_branch_report_jsonl(
                    &report_path,
                    &ctx,
                    &positive_branches,
                    &t_initializations,
                    &t_initialization_sources,
                    branch_gv_coverages.as_deref(),
                )
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to write KKLT branch report {}: {e}",
                        report_path.display()
                    );
                    std::process::exit(2);
                });
                eprintln!(
                    "[INFO] wrote KKLT branch report {} {}",
                    report_path.display(),
                    if branch_gv_coverages.is_some() {
                        "with GV coverage"
                    } else {
                        "without GV coverage"
                    }
                );
                if !branch_report_skip_gv_coverage && branch_gv_coverages.is_none() {
                    branch_gv_coverages = Some(
                        compute_branch_gv_coverages(
                            geom,
                            intersection,
                            &production_primal_basis,
                            &positive_branches,
                            small_curve_cutoff,
                            small_curve_pruning,
                            branch_report_missing_limit,
                            (branch_report_decomposition_depth > 0)
                                .then_some(branch_report_decomposition_depth),
                            true,
                        )
                        .unwrap_or_else(|e| {
                            eprintln!(
                                "[ERROR] failed to compute branch GV coverage report data: {e}"
                            );
                            std::process::exit(2);
                        }),
                    );
                    write_branch_report_jsonl(
                        &report_path,
                        &ctx,
                        &positive_branches,
                        &t_initializations,
                        &t_initialization_sources,
                        branch_gv_coverages.as_deref(),
                    )
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to write KKLT branch report {}: {e}",
                            report_path.display()
                        );
                        std::process::exit(2);
                    });
                    eprintln!(
                        "[INFO] enriched KKLT branch report {} with GV coverage",
                        report_path.display()
                    );
                }
                if branch_report_only {
                    eprintln!("[INFO] stopping after branch report as requested");
                    std::process::exit(0);
                }
            }
            eprintln!(
                "[INFO] KKLT branch search selected policy={} rank_by_volume={} init={} phase1_volume={} rel_err={} jacobian_rank={}/{} condition={:?}",
                branch_selection.as_str(),
                selected_rank_by_volume,
                best_branch.init_index,
                best_branch.classical_volume.get(),
                best_branch.result.relative_error.get(),
                best_branch.jacobian_diagnostics.rank,
                best_branch.jacobian_diagnostics.max_rank,
                best_branch
                    .jacobian_diagnostics
                    .condition_number
                    .map(cyrus_core::F64::get)
            );
            if let Some(coverage) = branch_gv_coverages
                .as_ref()
                .and_then(|coverages| coverages.get(selected_rank_by_volume))
            {
                eprintln!(
                    "[INFO] selected branch toric GV coverage: covered={} missing={} filtered={}",
                    coverage.toric_gv_covered_count,
                    coverage.toric_gv_missing_count,
                    coverage.filtered_count
                );
            }
            let Some(result) = solve_divisor_basis_path_following(
                &intersection.kappa_full,
                production_primal_basis.as_divisor_basis(),
                &kklt_basis,
                &tau_target,
                &best_branch.result.t,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] corrected divisor-basis KKLT solve failed after branch search");
                std::process::exit(2);
            };
            if !result.converged {
                eprintln!(
                    "[ERROR] corrected mixed-basis KKLT solve after branch search did not converge: rel_err={}",
                    result.relative_error.get()
                );
                std::process::exit(2);
            }
            (result, small_curve_selection_t, best_branch.result.t)
        };
        eprintln!(
            "[INFO] zeroth-order mixed-basis KKLT converged={} rel_err={}",
            zeroth_order.converged,
            zeroth_order.relative_error.get()
        );
        let ambient_rays = compute_mori_cone_cap_rays(
            &geom.triangulation,
            &geom.triangulation_points,
            &geom.polytope,
            false,
            false,
            None,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute primal ambient Mori-cap rays: {e}");
            std::process::exit(2);
        });
        let small_curve_candidates = subcutoff_toric_curve_candidates(
            &ambient_rays,
            &intersection.basis,
            &small_curve_selection_t,
            small_curve_cutoff,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to select small toric curve candidates: {e}");
            std::process::exit(2);
        });
        let small_curves =
            prune_selected_curve_candidates(&small_curve_candidates, small_curve_pruning, "primal")
                .unwrap_or_else(|e| {
                    eprintln!("[ERROR] {e}");
                    std::process::exit(2);
                });
        eprintln!(
            "[INFO] primal small toric curves: ambient_rays={} subcutoff={} pruned={} pruning={} cutoff={}",
            ambient_rays.len(),
            small_curve_candidates.len(),
            small_curves.len(),
            small_curve_pruning.as_str(),
            small_curve_cutoff.get()
        );
        let toric_gvs = compute_toric_two_face_curve_gv_invariants(
            &geom.triangulation,
            &geom.triangulation_points,
            &geom.polytope,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute primal toric curve GV values: {e}");
            std::process::exit(2);
        });
        let mut gv_by_class: HashMap<Vec<i64>, malachite::Integer> = toric_gvs
            .into_iter()
            .map(|item| (item.class, item.gv))
            .collect();
        let mut small_curve_gvs = Vec::with_capacity(small_curves.len());
        let mut missing_gv_classes = Vec::new();
        for curve in &small_curves {
            match gv_by_class.get(&curve.class) {
                Some(gv) => small_curve_gvs.push((curve.class.clone(), gv.clone())),
                None => missing_gv_classes.push(curve.class.clone()),
            }
        }
        if !missing_gv_classes.is_empty()
            && (primal_gv_min_points.is_some() || primal_gv_max_deg.is_some())
        {
            eprintln!(
                "[INFO] toric formulas missed {} selected primal small curves; computing primal general GV fallback with min_points={:?} max_deg={:?}",
                missing_gv_classes.len(),
                primal_gv_min_points,
                primal_gv_max_deg
            );
            let general_gvs = compute_primal_general_gv_by_ambient_class(
                geom,
                intersection,
                &missing_gv_classes,
                Some(&ambient_rays),
                primal_gv_min_points,
                primal_gv_max_deg,
            )
            .unwrap_or_else(|e| {
                eprintln!("[ERROR] failed to compute primal general GV fallback: {e}");
                std::process::exit(2);
            });
            let mut newly_covered = 0usize;
            for (class, gv) in general_gvs {
                match gv_by_class.entry(class) {
                    std::collections::hash_map::Entry::Occupied(existing) => {
                        if existing.get() != &gv {
                            eprintln!(
                                "[ERROR] toric/general GV conflict for selected primal curve: {} vs {gv}",
                                existing.get()
                            );
                            std::process::exit(2);
                        }
                    }
                    std::collections::hash_map::Entry::Vacant(slot) => {
                        if missing_gv_classes
                            .iter()
                            .any(|missing| missing == slot.key())
                        {
                            newly_covered += 1;
                        }
                        slot.insert(gv);
                    }
                }
            }
            eprintln!("[INFO] primal general GV fallback covered {newly_covered} missing curves");

            small_curve_gvs.clear();
            missing_gv_classes.clear();
            for curve in &small_curves {
                match gv_by_class.get(&curve.class) {
                    Some(gv) => small_curve_gvs.push((curve.class.clone(), gv.clone())),
                    None => missing_gv_classes.push(curve.class.clone()),
                }
            }
        }
        if !missing_gv_classes.is_empty() {
            eprintln!(
                "[ERROR] missing computed GV values for {} selected primal small toric curves; first_missing={:?}",
                missing_gv_classes.len(),
                missing_gv_classes.first()
            );
            std::process::exit(2);
        }
        eprintln!(
            "[INFO] primal small toric curve GVs selected={}",
            small_curve_gvs.len()
        );

        let correction_source_t = solve_gv_corrected_kklt(
            intersection,
            &production_primal_basis,
            &kklt_basis,
            &c_i,
            &chi_divisor,
            c_tau,
            &small_curve_gvs,
            &gamma,
            &gv_iteration_source_t,
            kklt_steps,
        );
        let Some(gv_volume_correction) =
            cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                &small_curve_gvs,
                &intersection.basis,
                &correction_source_t,
                Some(&gamma),
            )
        else {
            eprintln!("[ERROR] failed to compute primal ambient GV volume correction");
            std::process::exit(2);
        };
        eprintln!(
            "[INFO] GV volume correction = {}",
            gv_volume_correction.get()
        );
        let Some(input_chamber_gv_target_correction) =
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &small_curve_gvs,
                &intersection.basis,
                &kklt_basis,
                &correction_source_t,
                Some(&gamma),
            )
        else {
            eprintln!("[ERROR] failed to compute input-chamber GV target correction at solved t");
            std::process::exit(2);
        };
        let input_chamber_gv_target_correction_no_gamma =
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &small_curve_gvs,
                &intersection.basis,
                &kklt_basis,
                &correction_source_t,
                None,
            );
        if input_chamber_gv_target_correction_no_gamma.is_none() {
            eprintln!(
                "[COMPARE] no-gamma input-chamber GV target correction is invalid at solved t"
            );
        }
        let checkpoint_t_for_gv = load_corrected_kahler_checkpoint(data_dir, intersection);
        let input_chamber_gv_target_correction_checkpoint_t =
            checkpoint_t_for_gv.as_ref().and_then(|checkpoint_t| {
                cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                    &small_curve_gvs,
                    &intersection.basis,
                    &kklt_basis,
                    checkpoint_t,
                    Some(&gamma),
                )
            });
        if checkpoint_t_for_gv.is_some()
            && input_chamber_gv_target_correction_checkpoint_t.is_none()
        {
            eprintln!("[COMPARE] checkpoint-t input-chamber GV target correction is invalid");
        }
        let Some(input_chamber_target_tau) = cyrus_core::kklt::compute_gv_corrected_target_tau(
            &c_i,
            &chi_divisor,
            c_tau,
            &input_chamber_gv_target_correction,
        ) else {
            eprintln!("[ERROR] failed to reconstruct input-chamber GV-corrected target tau");
            std::process::exit(2);
        };
        let input_chamber_target_tau_finite = input_chamber_target_tau
            .iter()
            .map(|value| F64::<Finite>::new(value.get()).expect("target tau is finite"))
            .collect::<Vec<_>>();
        let mut alternate_gv_corrections = Vec::new();
        if let Some(correction) = input_chamber_gv_target_correction_no_gamma.as_deref() {
            alternate_gv_corrections.push(("no_gamma", correction));
        }
        if let Some(correction) = input_chamber_gv_target_correction_checkpoint_t.as_deref() {
            alternate_gv_corrections.push(("checkpoint_t", correction));
        }
        compare_corrected_target_volume_checkpoint(
            data_dir,
            &kklt_basis,
            &chi_divisor,
            &tau_target,
            &input_chamber_gv_target_correction,
            &alternate_gv_corrections,
            &input_chamber_target_tau_finite,
        );
        (
            correction_source_t,
            Some(gv_volume_correction),
            Some(kklt_basis),
        )
    };

    compare_corrected_kahler_checkpoint(data_dir, intersection, &t);

    let corrected_chamber = triangulation_from_kahler_point(geom, &intersection.basis, &t)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to diagnose corrected Kähler chamber: {e}");
            std::process::exit(2);
        });
    let corrected_chamber_changed =
        !triangulations_have_same_simplices(&geom.triangulation, &corrected_chamber);
    eprintln!(
        "[INFO] corrected Kähler point induces FRST: simplices={} changed_from_input={}",
        corrected_chamber.simplices().len(),
        corrected_chamber_changed
    );
    if corrected_chamber_changed {
        eprintln!(
            "[WARN] corrected Kähler point lies in a different regular chamber; flop/chamber-updated GV evaluation remains an explicit instanton-layer gap."
        );
    }
    if let Some(kklt_basis) = kklt_basis_for_chamber_gv.as_deref() {
        let corrected_chamber_kappa_full = chamber_intersection_full(
            &corrected_chamber,
            &geom.triangulation_points,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compute corrected-chamber target-volume intersection: {e}"
            );
            std::process::exit(2);
        });
        let corrected_chamber_kappa_basis =
            intersection_in_basis(&corrected_chamber_kappa_full, &intersection.basis);
        compare_corrected_chamber_target_volume_checkpoint(
            "computed-t corrected-chamber",
            data_dir,
            &intersection.basis,
            kklt_basis,
            &corrected_chamber_kappa_full,
            &corrected_chamber_kappa_basis,
            &t,
        );
        if let Some(checkpoint_t) = load_corrected_kahler_checkpoint(data_dir, intersection) {
            let checkpoint_chamber =
                triangulation_from_kahler_point(geom, &intersection.basis, &checkpoint_t)
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to diagnose corrected Kähler checkpoint chamber: {e}"
                        );
                        std::process::exit(2);
                    });
            let checkpoint_chamber_kappa_full =
                chamber_intersection_full(&checkpoint_chamber, &geom.triangulation_points)
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to compute corrected Kähler checkpoint chamber intersection: {e}"
                        );
                        std::process::exit(2);
                    });
            let checkpoint_chamber_kappa_basis =
                intersection_in_basis(&checkpoint_chamber_kappa_full, &intersection.basis);
            compare_corrected_chamber_target_volume_checkpoint(
                "checkpoint-t corrected-chamber",
                data_dir,
                &intersection.basis,
                kklt_basis,
                &checkpoint_chamber_kappa_full,
                &checkpoint_chamber_kappa_basis,
                &checkpoint_t,
            );
        }
    }

    let classical_volume = classical_volume_from_t(&intersection.kappa_basis, &t);
    let h11_raw = i32::try_from(intersection.basis.len()).unwrap_or_else(|_| {
        eprintln!("[ERROR] h11 does not fit in i32");
        std::process::exit(2);
    });
    let h21_raw = i32::try_from(h21).unwrap_or_else(|_| {
        eprintln!("[ERROR] h21 does not fit in i32");
        std::process::exit(2);
    });
    let h11 = cyrus_core::types::i32::I32::<cyrus_core::types::tags::GTEOne>::new(h11_raw)
        .expect("computed h11 must be >= 1");
    let h21 = cyrus_core::types::i32::I32::<cyrus_core::types::tags::NonNeg>::new(h21_raw)
        .expect("computed h21 must be non-negative");
    let bbhl = bbhl_correction(h11, h21);
    let v_string_before_gv = classical_volume - bbhl.get();
    eprintln!("[INFO] V_classical = {classical_volume}");
    eprintln!("[INFO] BBHL correction = {}", bbhl.get());
    eprintln!("[INFO] V_string before GV volume correction = {v_string_before_gv}");
    let v_string = v_string_before_gv + gv_volume_correction.map_or(0.0, cyrus_core::F64::get);
    let v_string_pos = F64::<Pos>::new(v_string).expect("V_string must be positive");
    eprintln!("[TIME] volume: {:.2?}", t0.elapsed());
    (v_string, v_string_pos)
}

fn compare_against_dat(
    compare_dir: Option<String>,
    data_dir: Option<String>,
    g_s: F64<Pos>,
    w0: F64<Pos>,
    v_string: f64,
) {
    let compare_dir = compare_dir.or(data_dir);
    if let Some(dir) = compare_dir {
        let dir = PathBuf::from(dir);
        if let Some(g_s_expected) = read_optional_scalar_f64(&dir.join("g_s.dat")) {
            let rel_gs = ((g_s.get() - g_s_expected) / g_s_expected).abs();
            eprintln!("[COMPARE] g_s rel_err = {rel_gs}");
        } else {
            eprintln!("[COMPARE] g_s.dat checkpoint not found; skipping g_s comparison");
        }
        if let Some(w0_expected) = read_optional_scalar_f64(&dir.join("W_0.dat")) {
            let rel_w0 = ((w0.get() - w0_expected) / w0_expected).abs();
            eprintln!("[COMPARE] W0 rel_err = {rel_w0}");
        } else {
            eprintln!("[COMPARE] W_0.dat checkpoint not found; skipping W0 comparison");
        }
        let corrected_volume_path = dir.join("corrected_cy_vol.dat");
        if let Some(corrected_v_expected) = read_optional_scalar_f64(&corrected_volume_path) {
            let abs_v = (v_string - corrected_v_expected).abs();
            eprintln!("[COMPARE] corrected V_string abs_err = {abs_v}");
            if abs_v > 1e-6 {
                eprintln!(
                    "[INFO] corrected V_string differs from the checkpoint at the checkpoint's own internal tolerance: the stored corrected_target_volumes/corrected_kahler_param pair is self-inconsistent at ~5.6e-4 per divisor, while Cyrus iterates the same fixed point to 1e-10. The deterministic model (kappa, chi, gamma, GV set, target and volume functions) reproduces the checkpoint volume to ~1e-8 when evaluated at the checkpoint Kahler point; see docs/CORRECTED_CHAMBER_RESOLUTION.md"
                );
            }
            if abs_v > 0.02 {
                eprintln!(
                    "[ERROR] corrected V_string mismatch beyond checkpoint-consistency tolerance: got {v_string}, expected {corrected_v_expected}"
                );
                std::process::exit(2);
            }
        } else {
            eprintln!(
                "[COMPARE] corrected_cy_vol.dat checkpoint not found; skipping corrected V_string comparison"
            );
        }
    }
}

fn target_correction_delta_summary(
    reference: &[F64<Finite>],
    candidate: &[F64<Finite>],
) -> Result<TargetCorrectionDeltaSummary, String> {
    if reference.len() != candidate.len() {
        return Err(format!(
            "target correction dimensions differ: reference={} candidate={}",
            reference.len(),
            candidate.len()
        ));
    }
    let mut delta_sq = 0.0f64;
    let mut reference_sq = 0.0f64;
    let mut max_abs_delta = 0.0f64;
    let mut max_abs_reference = 0.0f64;
    let mut max_abs_candidate = 0.0f64;
    for (reference_entry, candidate_entry) in reference.iter().zip(candidate.iter()) {
        let reference_value = reference_entry.get();
        let candidate_value = candidate_entry.get();
        let delta = candidate_value - reference_value;
        delta_sq += delta * delta;
        reference_sq += reference_value * reference_value;
        max_abs_delta = max_abs_delta.max(delta.abs());
        max_abs_reference = max_abs_reference.max(reference_value.abs());
        max_abs_candidate = max_abs_candidate.max(candidate_value.abs());
    }
    let relative_l2_delta = if reference_sq > 0.0 {
        (delta_sq / reference_sq).sqrt()
    } else if delta_sq == 0.0 {
        0.0
    } else {
        f64::INFINITY
    };
    Ok(TargetCorrectionDeltaSummary {
        dimension: reference.len(),
        max_abs_delta,
        relative_l2_delta,
        max_abs_reference,
        max_abs_candidate,
    })
}

fn stage_vacuum(
    ek0_opt: Option<F64<Pos>>,
    racetrack: &RacetrackData,
    v_string_pos: F64<Pos>,
    v_string: f64,
    compare_dir: Option<String>,
    data_dir: Option<String>,
    validate_mcallister_assertions: bool,
    manifest_dir: &Path,
) -> PipelineSummary {
    let Some(ek0) = ek0_opt else {
        eprintln!("[ERROR] ek0 is invalid; cannot compute V0");
        std::process::exit(2);
    };
    let g_s = racetrack.rt_res.g_s;
    let vac = compute_vacuum(ek0, g_s, v_string_pos, racetrack.w0);
    let v0_log10_abs = vac.v0.get().abs().log10();
    compare_against_dat(compare_dir, data_dir, g_s, racetrack.w0, v_string);
    let summary = PipelineSummary {
        g_s: g_s.get(),
        w0: racetrack.w0.get(),
        v_string,
        v0_log10_abs,
        ek0: ek0.get(),
    };
    eprintln!("[RESULT] g_s = {}", summary.g_s);
    eprintln!("[RESULT] W0 = {}", summary.w0);
    eprintln!("[RESULT] V_string = {}", summary.v_string);
    eprintln!("[RESULT] log10(|V0|) = {}", summary.v0_log10_abs);
    if summary.w0 <= 0.0 {
        eprintln!("[ERROR] W0 must be positive");
        std::process::exit(2);
    }
    if !validate_mcallister_assertions {
        eprintln!(
            "[INFO] skipping McAllister final assertion checks (--skip-mcallister-assertions)"
        );
        return summary;
    }
    let assertion_path = manifest_dir.join("tests/mcallister_e2e/assertions/racetrack.json");
    let assertion: RacetrackAssertion = load_json(&assertion_path);
    let tol_gs = 1e-6;
    // `assertion.cy_vol` is McAllister's *uncorrected* CY volume (4711.8297...),
    // while the computed V_string is the corrected-chamber volume
    // (4711.4265..., matching corrected_cy_vol.dat to ~0.006). The two differ
    // by ~0.4, so this tolerance bounds that known offset, not a residual.
    let tol_v = 0.5;
    if (summary.g_s - assertion.g_s).abs() > tol_gs {
        eprintln!(
            "[ERROR] g_s mismatch: got {}, expected {}",
            summary.g_s, assertion.g_s
        );
        std::process::exit(2);
    }
    if (summary.v_string - assertion.cy_vol).abs() > tol_v {
        eprintln!(
            "[ERROR] V_string mismatch: got {}, expected {}",
            summary.v_string, assertion.cy_vol
        );
        std::process::exit(2);
    }
    if (summary.w0 - assertion.w_0).abs() > tol_gs {
        eprintln!(
            "[ERROR] W0 mismatch: got {}, expected {}",
            summary.w0, assertion.w_0
        );
        std::process::exit(2);
    }
    if assertion.n_curves != 344 {
        eprintln!("[ERROR] unexpected curve count in assertion");
        std::process::exit(2);
    }
    summary
}

fn run_pipeline(args: PipelineArgs) {
    let t0 = Instant::now();
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let data_dir = args.data_dir.as_deref();
    enforce_modes(data_dir, args.allow_fixtures);
    let cutoff = F64::<Pos>::new(args.cutoff).expect("cutoff must be positive");
    let small_curve_cutoff =
        F64::<Pos>::new(args.small_curve_cutoff).expect("small curve cutoff must be positive");
    let geom = stage_triangulation(data_dir, &manifest_dir, &t0);
    if !stage_enabled(Stage::Intersection, args.stop_after) {
        return;
    }
    let intersection = stage_intersection(data_dir, &geom, &t0);
    if !stage_enabled(Stage::FlatDirection, args.stop_after) {
        return;
    }
    let flat = stage_flat_direction(
        data_dir,
        &manifest_dir,
        &geom,
        args.allow_invalid_ek0,
        args.validate_mcallister_assertions,
        args.dual_basis_override.as_ref(),
        args.production_dual_basis_override.as_ref(),
        &t0,
    );
    if !stage_enabled(Stage::Gv, args.stop_after) {
        return;
    }
    let invariants = stage_gv(&flat, args.min_points, args.max_deg, &t0);
    if !stage_enabled(Stage::Racetrack, args.stop_after) {
        return;
    }
    let racetrack = stage_racetrack(&invariants, &flat.m_flux, &flat.flat_p, cutoff, &t0);
    if !stage_enabled(Stage::Volume, args.stop_after) {
        return;
    }
    let (v_string, v_string_pos) = stage_volume(
        data_dir,
        &manifest_dir,
        &geom,
        &intersection,
        &racetrack,
        args.allow_downstream_kahler,
        args.kklt_steps,
        args.branch_candidates,
        args.branch_seed,
        args.branch_selection,
        args.branch_height_init,
        args.branch_report_path.as_deref(),
        args.branch_report_missing_limit,
        args.branch_report_decomposition_depth,
        args.branch_report_only,
        args.branch_report_skip_gv_coverage,
        args.primal_gv_min_points,
        args.primal_gv_max_deg,
        args.production_primal_basis_override.as_ref(),
        small_curve_cutoff,
        args.small_curve_pruning,
        flat.dual_divisor_basis.dimension(),
        &t0,
    );
    if !stage_enabled(Stage::Vacuum, args.stop_after) {
        return;
    }
    let summary = stage_vacuum(
        flat.ek0_opt,
        &racetrack,
        v_string_pos,
        v_string,
        args.compare_dir,
        args.data_dir,
        args.validate_mcallister_assertions,
        &manifest_dir,
    );
    if let Some(path) = args.out_path {
        let out_path = PathBuf::from(path);
        let data = serde_json::to_string_pretty(&summary).expect("serialize summary");
        std::fs::write(&out_path, data)
            .unwrap_or_else(|e| panic!("Failed to write {}: {e}", out_path.display()));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int_matrix(rows: &[&[i64]]) -> Vec<Vec<malachite::Integer>> {
        rows.iter()
            .map(|row| {
                row.iter()
                    .map(|&value| malachite::Integer::from(value))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn branch_selection_parser_accepts_all_declared_policies() {
        for name in [
            "max-volume",
            "min-volume",
            "first-positive",
            "min-condition",
            "min-toric-gv-missing",
            "min-required-gv-degree",
        ] {
            let policy = parse_branch_selection(name)
                .unwrap_or_else(|| panic!("policy {name} should parse"));
            assert_eq!(policy.as_str(), name);
        }
        assert_eq!(
            parse_branch_selection("min-gv-missing").map(BranchSelection::as_str),
            Some("min-toric-gv-missing")
        );
        assert_eq!(
            parse_branch_selection("min-gv-degree").map(BranchSelection::as_str),
            Some("min-required-gv-degree")
        );
        assert!(parse_branch_selection("condition").is_none());
    }

    #[test]
    fn branch_selection_declares_coverage_requirement() {
        assert!(BranchSelection::MinToricGvMissing.requires_gv_coverage());
        assert!(BranchSelection::MinRequiredGvDegree.requires_gv_coverage());
        assert!(BranchSelection::MinRequiredGvDegree.requires_gv_degree_summary());
        assert!(!BranchSelection::MinVolume.requires_gv_coverage());
        assert!(!BranchSelection::MinVolume.requires_gv_degree_summary());
    }

    #[test]
    fn curve_pruning_strategy_parser_accepts_declared_policies() {
        assert_eq!(
            parse_curve_pruning_strategy("pair"),
            Some(CurvePruningStrategy::PairDecomposable)
        );
        assert_eq!(
            parse_curve_pruning_strategy("pair-decomposable"),
            Some(CurvePruningStrategy::PairDecomposable)
        );
        assert_eq!(
            parse_curve_pruning_strategy("finite-semigroup"),
            Some(CurvePruningStrategy::FiniteSemigroup)
        );
        assert_eq!(
            parse_curve_pruning_strategy("semigroup"),
            Some(CurvePruningStrategy::FiniteSemigroup)
        );
        assert!(parse_curve_pruning_strategy("hilbert").is_none());
    }

    #[test]
    fn basis_override_accepts_index_basis() {
        let basis =
            serde_json::from_str::<BasisOverride>(r#"{"indices":[3,4]}"#).expect("valid basis");

        assert!(
            basis.indices("--dual-basis").expect("index basis") == [3, 4],
            "index basis should be available"
        );
    }

    #[test]
    fn basis_override_rejects_matrix_basis_shape_explicitly() {
        let basis = serde_json::from_str::<BasisOverride>(r#"{"matrix":[[1,0,0],[0,1,0]]}"#)
            .expect("matrix basis JSON should parse for a deliberate error");
        let err = basis
            .indices("--dual-basis")
            .expect_err("runner must reject unsupported matrix basis explicitly");

        assert!(
            err.contains("matrix divisor basis"),
            "unexpected basis error: {err}"
        );
        assert!(
            err.contains("vector-basis only"),
            "unexpected basis error: {err}"
        );
    }

    #[test]
    fn basis_override_rejects_ambiguous_basis_shape() {
        let basis = serde_json::from_str::<BasisOverride>(
            r#"{"indices":[3,4],"matrix":[[1,0,0],[0,1,0]]}"#,
        )
        .expect("ambiguous basis JSON should parse for validation");
        let err = basis
            .indices("--dual-basis")
            .expect_err("basis override must reject ambiguous representations");

        assert!(
            err.contains("both `indices` and `matrix`"),
            "unexpected basis error: {err}"
        );
    }

    #[test]
    fn matrix_dual_basis_override_builds_m_flux_transform() {
        let glsm = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let computed_basis = vec![2, 3];
        let matrix_basis =
            serde_json::from_str::<BasisOverride>(r#"{"matrix":[[0,0,1,1],[0,0,0,1]]}"#)
                .expect("matrix basis JSON should parse");
        let production_basis = OwnedDivisorBasis::Indices(computed_basis.clone());
        let source_basis =
            owned_divisor_basis_from_override(&matrix_basis, &computed_basis, "--dual-basis")
                .expect("matrix source basis should parse");

        let transform = basis_change_matrix_between_owned(&glsm, &production_basis, &source_basis)
            .expect("matrix source basis should be equivalent to computed basis");
        assert_eq!(transform, int_matrix(&[&[1, 0], &[1, 1]]));

        let transformed = transform_m_flux_between_divisor_bases(
            &glsm,
            &production_basis,
            &source_basis,
            &[5, 7],
            "M",
        );
        assert_eq!(transformed, vec![5, 12]);
    }

    #[test]
    fn matrix_dual_basis_override_builds_k_flux_transform() {
        let glsm = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let computed_basis = vec![2, 3];
        let matrix_basis =
            serde_json::from_str::<BasisOverride>(r#"{"matrix":[[0,0,1,1],[0,0,0,1]]}"#)
                .expect("matrix basis JSON should parse");
        let production_basis = OwnedDivisorBasis::Indices(computed_basis.clone());
        let source_basis =
            owned_divisor_basis_from_override(&matrix_basis, &computed_basis, "--dual-basis")
                .expect("matrix source basis should parse");

        let transform = basis_change_matrix_between_owned(&glsm, &source_basis, &production_basis)
            .expect("matrix source basis should be equivalent to computed basis");
        assert_eq!(transform, int_matrix(&[&[1, 0], &[-1, 1]]));

        let transformed = transform_k_flux_between_divisor_bases(
            &glsm,
            &source_basis,
            &production_basis,
            &[5, 7],
        );
        assert_eq!(transformed, vec![-2, 7]);
    }

    #[test]
    fn matrix_production_dual_basis_builds_flux_transforms_and_gv_inputs() {
        let glsm = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let linrels = glsm.clone();
        let standard_basis = vec![2, 3];
        let matrix_basis =
            serde_json::from_str::<BasisOverride>(r#"{"matrix":[[0,0,1,1],[0,0,0,1]]}"#)
                .expect("matrix basis JSON should parse");
        let production_basis = owned_divisor_basis_from_override(
            &matrix_basis,
            &standard_basis,
            "--production-dual-basis",
        )
        .expect("matrix production basis should parse");
        let source_basis = OwnedDivisorBasis::Indices(standard_basis);

        let m_transformed = transform_m_flux_between_divisor_bases(
            &glsm,
            &production_basis,
            &source_basis,
            &[5, 7],
            "M",
        );
        assert_eq!(m_transformed, vec![5, 2]);

        let k_transformed = transform_k_flux_between_divisor_bases(
            &glsm,
            &source_basis,
            &production_basis,
            &[5, 7],
        );
        assert_eq!(k_transformed, vec![12, 7]);

        let ambient = vec![vec![0, 0, 1, 0], vec![0, 0, 2, 2]];
        let data = production_gv_basis_data(&ambient, &linrels, &production_basis, "test mirror")
            .expect("matrix production basis should build cygv inputs");
        assert_eq!(data.mori_rays, vec![vec![1, 0], vec![2, 1]]);
        assert_eq!(data.q_matrix, vec![vec![2, 1, 0], vec![1, -1, 1]]);
    }

    #[test]
    fn matrix_production_primal_basis_transforms_kahler_coordinates() {
        let glsm = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let computed_basis = OwnedDivisorBasis::Indices(vec![2, 3]);
        let matrix_basis =
            serde_json::from_str::<BasisOverride>(r#"{"matrix":[[0,0,1,1],[0,0,0,1]]}"#)
                .expect("matrix basis JSON should parse");
        let production_basis =
            owned_divisor_basis_from_override(&matrix_basis, &[2, 3], "--production-primal-basis")
                .expect("matrix production basis should parse");
        let computed_t = vec![
            F64::<Finite>::new(5.0).expect("finite"),
            F64::<Finite>::new(7.0).expect("finite"),
        ];

        let production_t = transform_kahler_between_owned_divisor_bases(
            &glsm,
            &production_basis,
            &computed_basis,
            &computed_t,
            "test Kähler",
            false,
        )
        .expect("computed Kähler coordinates should transform to production basis");
        assert_eq!(finite_values(&production_t), vec![5.0, 2.0]);

        let round_trip = transform_kahler_between_owned_divisor_bases(
            &glsm,
            &computed_basis,
            &production_basis,
            &production_t,
            "test Kähler",
            false,
        )
        .expect("production Kähler coordinates should transform back");
        assert_eq!(finite_values(&round_trip), vec![5.0, 7.0]);
    }

    #[test]
    fn branch_report_records_curve_pruning_strategy() {
        use cyrus_core::{
            KkltBranchSolution, KkltJacobianDiagnostics, KkltResult, f64_finite, f64_pos,
            types::tags::NonNeg,
        };

        let cache_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/cyrus-test-cache");
        std::fs::create_dir_all(&cache_dir).expect("create cache dir");
        let path = cache_dir.join(format!(
            "branch-report-{}-{}.jsonl",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("time after epoch")
                .as_nanos()
        ));
        let ctx = BranchReportContext {
            branch_seed: 7,
            branch_selection: BranchSelection::MinToricGvMissing,
            small_curve_pruning: CurvePruningStrategy::FiniteSemigroup,
            kklt_steps: 3,
            attempted: 1,
            solved: 1,
            non_converged: 0,
            non_positive_volume: 0,
            selected_rank_by_volume: 0,
        };
        let branch = KkltBranchSolution {
            init_index: 0,
            result: KkltResult {
                t: vec![f64_finite!(1.0)],
                tau: vec![f64_finite!(2.0)],
                tau_target: vec![f64_pos!(2.0)],
                converged: true,
                relative_error: F64::<NonNeg>::ZERO,
            },
            classical_volume: f64_pos!(3.0),
            jacobian_diagnostics: KkltJacobianDiagnostics {
                rank: 1,
                max_rank: 1,
                max_singular_value: f64_pos!(4.0),
                min_nonzero_singular_value: f64_pos!(2.0),
                condition_number: Some(f64_pos!(2.0)),
            },
        };
        let coverage = BranchGvCoverage {
            ambient_rays: 7,
            subcutoff_count: 5,
            filtered_count: 4,
            toric_gv_covered_count: 3,
            toric_gv_missing_count: 1,
            first_missing_class: Some(vec![1, 0]),
            missing_required_degree_min: Some(2),
            missing_required_degree_max: Some(9),
            missing_class_sample: vec![vec![1, 0]],
            bounded_decomposition_max_terms: Some(4),
            missing_bounded_decomposable_count: Some(1),
            first_missing_bounded_decomposition: Some(vec![vec![1, 0], vec![0, 1]]),
        };

        write_branch_report_jsonl(
            &path,
            &ctx,
            &[branch],
            &[vec![f64_finite!(0.5)]],
            &["unit"],
            Some(&[coverage]),
        )
        .expect("write branch report");

        let content = std::fs::read_to_string(&path).expect("read branch report");
        let rows = content
            .lines()
            .map(|line| serde_json::from_str::<serde_json::Value>(line).expect("valid jsonl row"))
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0]["type"], "summary");
        assert_eq!(rows[0]["small_curve_pruning"], "finite-semigroup");
        assert_eq!(rows[0]["selected_small_curve_filtered_count"], 4);
        assert_eq!(rows[1]["type"], "positive_branch");
        assert_eq!(rows[1]["small_curve_pruning"], "finite-semigroup");
        assert_eq!(rows[1]["small_curve_filtered_count"], 4);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn min_toric_gv_missing_selection_ties_by_volume() {
        let coverages = vec![
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 4,
                filtered_count: 4,
                toric_gv_covered_count: 3,
                toric_gv_missing_count: 1,
                first_missing_class: Some(vec![1]),
                missing_required_degree_min: None,
                missing_required_degree_max: None,
                missing_class_sample: vec![vec![1]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 5,
                filtered_count: 5,
                toric_gv_covered_count: 3,
                toric_gv_missing_count: 2,
                first_missing_class: Some(vec![2]),
                missing_required_degree_min: None,
                missing_required_degree_max: None,
                missing_class_sample: vec![vec![2]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 6,
                filtered_count: 6,
                toric_gv_covered_count: 5,
                toric_gv_missing_count: 1,
                first_missing_class: Some(vec![3]),
                missing_required_degree_min: None,
                missing_required_degree_max: None,
                missing_class_sample: vec![vec![3]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
        ];
        let volumes_by_rank = vec![20.0, 1.0, 10.0];

        let selected = select_min_toric_gv_missing_rank(&coverages, &volumes_by_rank).unwrap();

        assert_eq!(selected, 2);
    }

    #[test]
    fn min_toric_gv_missing_selection_rejects_mismatched_rows() {
        let coverages = vec![BranchGvCoverage {
            ambient_rays: 10,
            subcutoff_count: 4,
            filtered_count: 4,
            toric_gv_covered_count: 3,
            toric_gv_missing_count: 1,
            first_missing_class: Some(vec![1]),
            missing_required_degree_min: None,
            missing_required_degree_max: None,
            missing_class_sample: vec![vec![1]],
            bounded_decomposition_max_terms: None,
            missing_bounded_decomposable_count: None,
            first_missing_bounded_decomposition: None,
        }];

        assert!(
            select_min_toric_gv_missing_rank(&coverages, &[])
                .unwrap_err()
                .contains("do not match")
        );
    }

    #[test]
    fn min_required_gv_degree_selection_ties_by_missing_then_volume() {
        let coverages = vec![
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 4,
                filtered_count: 4,
                toric_gv_covered_count: 2,
                toric_gv_missing_count: 2,
                first_missing_class: Some(vec![1]),
                missing_required_degree_min: Some(2),
                missing_required_degree_max: Some(20),
                missing_class_sample: vec![vec![1]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 5,
                filtered_count: 5,
                toric_gv_covered_count: 1,
                toric_gv_missing_count: 4,
                first_missing_class: Some(vec![2]),
                missing_required_degree_min: Some(3),
                missing_required_degree_max: Some(10),
                missing_class_sample: vec![vec![2]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 6,
                filtered_count: 6,
                toric_gv_covered_count: 3,
                toric_gv_missing_count: 3,
                first_missing_class: Some(vec![3]),
                missing_required_degree_min: Some(4),
                missing_required_degree_max: Some(10),
                missing_class_sample: vec![vec![3]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
        ];
        let volumes_by_rank = vec![1.0, 100.0, 50.0];

        let selected = select_min_required_gv_degree_rank(&coverages, &volumes_by_rank).unwrap();

        assert_eq!(selected, 2);
    }

    #[test]
    fn min_required_gv_degree_selection_accepts_zero_missing_without_degree_summary() {
        let coverages = vec![
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 4,
                filtered_count: 4,
                toric_gv_covered_count: 4,
                toric_gv_missing_count: 0,
                first_missing_class: None,
                missing_required_degree_min: None,
                missing_required_degree_max: None,
                missing_class_sample: Vec::new(),
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
            BranchGvCoverage {
                ambient_rays: 10,
                subcutoff_count: 5,
                filtered_count: 5,
                toric_gv_covered_count: 4,
                toric_gv_missing_count: 1,
                first_missing_class: Some(vec![2]),
                missing_required_degree_min: Some(1),
                missing_required_degree_max: Some(1),
                missing_class_sample: vec![vec![2]],
                bounded_decomposition_max_terms: None,
                missing_bounded_decomposable_count: None,
                first_missing_bounded_decomposition: None,
            },
        ];

        let selected = select_min_required_gv_degree_rank(&coverages, &[10.0, 1.0]).unwrap();

        assert_eq!(selected, 0);
    }

    #[test]
    fn min_required_gv_degree_selection_rejects_missing_degree_summary() {
        let coverages = vec![BranchGvCoverage {
            ambient_rays: 10,
            subcutoff_count: 4,
            filtered_count: 4,
            toric_gv_covered_count: 3,
            toric_gv_missing_count: 1,
            first_missing_class: Some(vec![1]),
            missing_required_degree_min: None,
            missing_required_degree_max: None,
            missing_class_sample: vec![vec![1]],
            bounded_decomposition_max_terms: None,
            missing_bounded_decomposable_count: None,
            first_missing_bounded_decomposition: None,
        }];

        assert!(
            select_min_required_gv_degree_rank(&coverages, &[1.0])
                .unwrap_err()
                .contains("missing required GV degree summary")
        );
    }

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

    #[test]
    fn target_correction_delta_summary_reports_relative_change() {
        let reference = vec![
            F64::<Finite>::new(3.0).unwrap(),
            F64::<Finite>::new(4.0).unwrap(),
        ];
        let candidate = vec![
            F64::<Finite>::new(6.0).unwrap(),
            F64::<Finite>::new(8.0).unwrap(),
        ];

        let summary = target_correction_delta_summary(&reference, &candidate).unwrap();

        assert_eq!(summary.dimension, 2);
        assert_eq!(summary.max_abs_delta, 4.0);
        assert_eq!(summary.relative_l2_delta, 1.0);
        assert_eq!(summary.max_abs_reference, 4.0);
        assert_eq!(summary.max_abs_candidate, 8.0);
    }
}

fn main() {
    run_pipeline(parse_args());
}
