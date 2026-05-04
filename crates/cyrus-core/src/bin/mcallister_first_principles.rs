//! First-principles McAllister pipeline runner (opt-in).
//!
//! Runs the pipeline stages explicitly and fails fast on any mismatch or
//! invalid physics. This is intentionally a binary, not a test.
//!
//! Optional:
//! - `--dual-basis path/to/dual_basis.json` to supply the source coordinate
//!   basis for K/M flux vectors. Without this, McAllister validation mode uses
//!   the paper basis `[3, 4, 5, 8]`; generic mode uses Cyrus' computed dual basis.
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
//!   may be at most 3.
//! - `--branch-report-only` to stop after writing that report.
//! - `--branch-report-skip-gv-coverage` to omit the expensive per-branch
//!   small-curve/GV coverage enrichment.
//! - `--primal-gv-max-deg N` or `--primal-gv-min-points N` to compute general
//!   primal GV invariants if toric formulas do not cover selected small curves.
//! - `--diagnose-corrected-chamber-gv` to recompute toric small-curve/GV data
//!   in the FRST induced by the solved corrected Kähler point.
//!   When combined with primal GV bounds, also tries bounded general GV fallback
//!   for corrected-chamber small curves not covered by toric formulas.
//! - `--diagnose-corrected-chamber-provided-generators-gv` to also try the
//!   CYTools `mcap_generators` override path for corrected-chamber misses,
//!   explicitly without Mori-cone lattice augmentation.
//! - `--diagnose-corrected-chamber-ray-gv` to try one-dimensional HKTY ray
//!   diagnostics for missing corrected-chamber primitive Mori generators.
//! - `--diagnose-corrected-chamber-lp-face-gv` to try small HKTY face
//!   diagnostics generated from LP decomposition witnesses of missing
//!   corrected-chamber primitive Mori generators.
//! - `--diagnose-chamber-updated-kklt` to run a diagnostic-only KKLT
//!   fixed-point loop that recomputes the FRST chamber, intersections, divisor
//!   χ, and toric-covered small-curve GV target correction at each iteration.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::path::PathBuf;
use std::time::Instant;

use good_lp::{
    Expression, ProblemVariables, ResolutionError, Solution, SolverModel, default_solver, variable,
};

use cyrus_core::flat_direction::{compute_flat_direction, compute_flat_direction_full};
use cyrus_core::gv::{BoundedCurveDecompositionIndex, project_mori_cone_cap_rays_to_basis};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::range::CheckedRange;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::vacuum::compute_vacuum;
use cyrus_core::volume::bbhl_correction;
use cyrus_core::{
    Point, Polytope, Triangulation, basis_change_matrix, build_racetrack_terms,
    compute_curve_basis_matrix, compute_glsm_and_linrels, compute_grading_vector,
    compute_intersection_cytools, compute_linear_relations_no_origin, compute_mori_cone_cap_rays,
    compute_origin_circuit_curve_diagnostics, compute_regular_triangulation,
    compute_toric_two_face_curve_gv_invariants, compute_w0_from_terms,
    effective_prime_divisors_from_curve_basis, generate_scaled_kklt_branch_initializations,
    heights_to_kahler, intersection_in_basis, is_unimodular, kahler_to_heights,
    map_basis_gv_invariants_to_ambient, remove_pair_decomposable_curve_candidates,
    scale_mixed_basis_kklt_branch_initialization_to_target, solve_mixed_basis_path_following,
    solve_mixed_basis_path_following_branch_candidates, solve_racetrack,
    subcutoff_toric_curve_candidates,
};

const DEFAULT_MCALLISTER_GV_MIN_POINTS: u32 = 20_000;
const DEFAULT_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT: usize = 100_000;
const DEFAULT_CORRECTED_CHAMBER_PROVIDED_GV_GENERATOR_LIMIT: usize = 2_000;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_SPAN_GENERATOR_LIMIT: usize = 64;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_LATTICE_GENERATOR_LIMIT: usize = 64;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_LATTICE_ELEMENT_LIMIT: usize = 50_000;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_TERMS: usize = 3;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_WITNESSES: usize = 8;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_DECOMPOSITION_CLOSURE_ELEMENT_LIMIT: usize = 20_000;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_RAY_LIMIT: usize = 1_000_000;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_ANCHOR_ATTEMPTS: usize = 16;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_CUTTING_ROUNDS: usize = 64;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_SCALE_LIMIT: i64 = 100_000;
const DEFAULT_CORRECTED_CHAMBER_COVERED_GV_DIVISOR_REPRESENTATION_CLASS_LIMIT: usize = 24;

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
    fn as_str(self) -> &'static str {
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
struct BasisOverride {
    indices: Vec<usize>,
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

struct ChamberGvDiagnostic {
    ambient_rays: usize,
    subcutoff_count: usize,
    filtered_count: usize,
    toric_gv_covered_count: usize,
    toric_gv_missing_count: usize,
    basis_mori_ray_count: Option<usize>,
    degree_bounded_basis_mori_ray_count: Option<usize>,
    basis_mori_ray_degree_min: Option<i128>,
    basis_mori_ray_degree_max: Option<i128>,
    general_gv_covered_count: Option<usize>,
    ray_gv_covered_count: Option<usize>,
    ray_gv_volume_correction: Option<F64<Finite>>,
    ray_gv_sample: Vec<RayGvDiagnosticSample>,
    lp_face_gv_covered_count: Option<usize>,
    lp_face_gv_failed_count: Option<usize>,
    lp_face_gv_certified_count: Option<usize>,
    lp_face_gv_uncertified_count: Option<usize>,
    lp_face_gv_volume_correction: Option<F64<Finite>>,
    lp_face_gv_sample: Vec<FaceGvDiagnosticSample>,
    combined_diagnostic_gv_covered_count: Option<usize>,
    combined_diagnostic_gv_missing_count: Option<usize>,
    combined_diagnostic_gv_zero_count: Option<usize>,
    combined_diagnostic_gv_nonzero_count: Option<usize>,
    combined_diagnostic_gv_volume_correction: Option<F64<Finite>>,
    combined_diagnostic_gv_target_correction: Option<Vec<F64<Finite>>>,
    remaining_gv_missing_count: usize,
    first_missing_class: Option<Vec<i64>>,
    missing_required_degree_min: Option<i128>,
    missing_required_degree_max: Option<i128>,
    missing_target_stats: Option<MissingGvTargetStats>,
    covered_toric_gv_divisor_representation_baseline:
        Option<CoveredToricGvDivisorRepresentationBaseline>,
    covered_gv_target_correction: Option<Vec<F64<Finite>>>,
    covered_gv_volume_correction: Option<F64<Finite>>,
    gv_volume_correction: Option<F64<Finite>>,
}

struct ChamberToricGvSelection {
    ambient_rays: usize,
    subcutoff_count: usize,
    filtered_count: usize,
    toric_gv_covered_count: usize,
    toric_gv_missing_count: usize,
    first_missing_class: Option<Vec<i64>>,
    small_curve_gvs: Vec<(Vec<i64>, malachite::Integer)>,
}

struct ChamberUpdatedKkltDiagnostic {
    iterations: usize,
    converged: bool,
    final_t: Vec<F64<Finite>>,
    final_classical_volume: f64,
    final_gv_volume_correction: F64<Finite>,
    final_toric_missing_count: usize,
    final_first_missing_class: Option<Vec<i64>>,
}

#[derive(Clone, Debug, PartialEq)]
struct TargetCorrectionDeltaSummary {
    dimension: usize,
    max_abs_delta: f64,
    relative_l2_delta: f64,
    max_abs_reference: f64,
    max_abs_candidate: f64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MissingGvTargetStats {
    target_count: usize,
    targets_that_are_mori_generators: usize,
    targets_that_are_origin_circuits: usize,
    targets_real_cone_decomposable_by_other_generators: usize,
    targets_that_are_lp_extremal_mori_generators: usize,
    real_cone_decomposition_active_generator_min: Option<usize>,
    real_cone_decomposition_active_generator_max: Option<usize>,
    origin_circuit_resolved_conifold_count: usize,
    min_generators_le_target_degree: usize,
    max_generators_le_target_degree: usize,
    origin_coefficient_counts: BTreeMap<i64, usize>,
    origin_circuit_pattern_counts: BTreeMap<String, usize>,
    sample: Vec<MissingGvTargetSample>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MissingGvTargetSample {
    degree: i128,
    generators_le_degree: usize,
    is_mori_generator: bool,
    origin_circuit_pattern: Option<String>,
    origin_circuit_witness_count: Option<usize>,
    origin_circuit_first_witness: Option<OriginCircuitWitnessSample>,
    cms_general_divisor_shape_candidates: Option<Vec<CmsGeneralDivisorShapeCandidate>>,
    cms_general_divisor_intersection_checks: Option<Vec<CmsGeneralDivisorIntersectionCheck>>,
    real_cone_decomposable_by_other_generators: bool,
    real_cone_decomposition_active_generators: Option<usize>,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct CmsGeneralDivisorShapeCandidate {
    shrinking_divisor_index: usize,
    shrinking_divisor_coefficient: i64,
    shrinking_divisor_coordinates: Vec<i64>,
    inferred_other_normal_degree: i64,
    toric_gv1_formula_value: Option<i64>,
    all_non_origin_relation_points_are_two_face: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct CmsGeneralDivisorIntersectionCheck {
    shrinking_divisor_index: usize,
    has_rational_divisor_solution: bool,
    solution_basis_support_len: Option<usize>,
    solution_is_integral: Option<bool>,
    computed_other_normal_degree: Option<String>,
    matches_inferred_other_normal_degree: Option<bool>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CoveredToricGvDivisorRepresentationBaseline {
    checked_class_count: usize,
    class_limit: usize,
    support_index_checks: usize,
    classes_with_support_divisor_solution: usize,
    first_without_support_divisor_solution: Option<Vec<(usize, i64)>>,
    sample: Vec<CoveredToricGvDivisorRepresentationSample>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CoveredToricGvDivisorRepresentationSample {
    ambient_nonzero: Vec<(usize, i64)>,
    support_indices_with_solution: Vec<usize>,
    support_indices_without_solution: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct DivisorRepresentationSolution {
    coefficients: Vec<malachite::Rational>,
    other_normal_degree: malachite::Rational,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OriginCircuitWitnessSample {
    first_facet_exclusive_point: usize,
    second_facet_exclusive_point: usize,
    shared_two_simplex: Vec<usize>,
    first_facet_size: usize,
    second_facet_size: usize,
    sparse_relation: Vec<(usize, i64)>,
    relation_points: Vec<OriginCircuitRelationPointSample>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OriginCircuitRelationPointSample {
    point_index: usize,
    coefficient: i64,
    coordinates: Vec<i64>,
    face_dimension: Option<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RealConeDecompositionWitness {
    active_generator_indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct OneDimensionalRayGvTargets {
    candidates: Vec<Vec<i64>>,
    skipped_non_generators: usize,
    skipped_decomposable_generators: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RayGvDiagnosticSample {
    degree: i128,
    gv: malachite::Integer,
    ambient_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct FaceGvDiagnosticSample {
    degree: i128,
    generator_count: usize,
    active_generator_count: usize,
    span_generator_count: usize,
    used_span_expansion: bool,
    used_lattice_saturation: bool,
    used_integer_decomposition: bool,
    used_decomposition_diamond: bool,
    used_decomposition_closure: bool,
    integer_decomposition_term_count: Option<usize>,
    lattice_semigroup_element_count: Option<usize>,
    supporting_face_certificate: Option<SupportingFaceCertificateSummary>,
    gv: Option<malachite::Integer>,
    error: Option<String>,
    ambient_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SupportingFaceCertificateSummary {
    zero_generator_count: usize,
    positive_generator_count: usize,
    normal_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug)]
struct SupportingFaceCertificate {
    zero_generator_count: usize,
    positive_generator_count: usize,
    normal: Vec<i64>,
}

#[derive(Clone, Debug)]
struct FaceGvDiagnosticResult {
    ambient_class: Vec<i64>,
    degree: i128,
    generator_count: usize,
    active_generator_count: usize,
    span_generator_count: usize,
    used_span_expansion: bool,
    used_lattice_saturation: bool,
    used_integer_decomposition: bool,
    used_decomposition_diamond: bool,
    used_decomposition_closure: bool,
    integer_decomposition_term_count: Option<usize>,
    lattice_semigroup_element_count: Option<usize>,
    supporting_face_certificate: Option<SupportingFaceCertificate>,
    gv: Option<malachite::Integer>,
    error: Option<String>,
}

#[derive(Clone, Debug)]
struct FaceGvAttempt {
    generator_count: usize,
    active_generator_count: usize,
    span_generator_count: usize,
    used_span_expansion: bool,
    used_lattice_saturation: bool,
    used_integer_decomposition: bool,
    used_decomposition_diamond: bool,
    used_decomposition_closure: bool,
    integer_decomposition_term_count: Option<usize>,
    lattice_semigroup_element_count: Option<usize>,
    supporting_face_certificate: Option<SupportingFaceCertificate>,
    gv: Option<malachite::Integer>,
    error: Option<String>,
}

struct BranchReportContext {
    branch_seed: u64,
    branch_selection: BranchSelection,
    kklt_steps: usize,
    attempted: usize,
    solved: usize,
    non_converged: usize,
    non_positive_volume: usize,
    selected_rank_by_volume: usize,
}

struct PipelineArgs {
    stop_after: Stage,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    cutoff: f64,
    small_curve_cutoff: f64,
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
    diagnose_corrected_chamber_gv: bool,
    diagnose_corrected_chamber_provided_generators_gv: bool,
    diagnose_corrected_chamber_ray_gv: bool,
    diagnose_corrected_chamber_lp_face_gv: bool,
    diagnose_chamber_updated_kklt: bool,
    diagnose_chamber_updated_kklt_iterations: usize,
    dual_basis_override: Option<BasisOverride>,
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
    dual_basis: Vec<usize>,
    dual_kappa: cyrus_core::Intersection,
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
    let small_curve_cutoff = parse_arg_value::<f64>("--small-curve-cutoff").unwrap_or(1.0);
    let out_path = parse_arg_value::<String>("--out");
    let compare_dir = parse_arg_value::<String>("--compare-dir");
    let data_dir = parse_arg_value::<String>("--data-dir")
        .or_else(|| std::env::var("CYRUS_MCALLISTER_DATA_DIR").ok());
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
    let diagnose_corrected_chamber_gv = parse_flag("--diagnose-corrected-chamber-gv");
    let diagnose_corrected_chamber_provided_generators_gv =
        parse_flag("--diagnose-corrected-chamber-provided-generators-gv");
    let diagnose_corrected_chamber_ray_gv = parse_flag("--diagnose-corrected-chamber-ray-gv");
    let diagnose_corrected_chamber_lp_face_gv =
        parse_flag("--diagnose-corrected-chamber-lp-face-gv");
    let diagnose_chamber_updated_kklt = parse_flag("--diagnose-chamber-updated-kklt");
    let diagnose_chamber_updated_kklt_iterations =
        parse_arg_value::<usize>("--diagnose-chamber-updated-kklt-iterations").unwrap_or(6);
    let dual_basis_override = parse_arg_value::<String>("--dual-basis")
        .map(|path| load_json::<BasisOverride>(&PathBuf::from(path)));
    PipelineArgs {
        stop_after,
        max_deg,
        min_points,
        cutoff,
        small_curve_cutoff,
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
        diagnose_corrected_chamber_gv,
        diagnose_corrected_chamber_provided_generators_gv,
        diagnose_corrected_chamber_ray_gv,
        diagnose_corrected_chamber_lp_face_gv,
        diagnose_chamber_updated_kklt,
        diagnose_chamber_updated_kklt_iterations,
        dual_basis_override,
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

fn load_primal_inputs(data_dir: Option<&str>, manifest_dir: &PathBuf) -> (Vec<Vec<i64>>, Vec<f64>) {
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

fn load_flux_vectors(data_dir: Option<&str>, manifest_dir: &PathBuf) -> (Vec<i64>, Vec<i64>) {
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

fn transform_i64_coordinates(
    transform: &[Vec<malachite::Integer>],
    values: &[i64],
    label: &str,
) -> Vec<i64> {
    if transform.len() != values.len() || transform.iter().any(|row| row.len() != values.len()) {
        eprintln!("[ERROR] {label} basis transform shape does not match vector length");
        std::process::exit(2);
    }

    transform
        .iter()
        .map(|row| {
            let mut acc = malachite::Integer::from(0);
            for (coeff, &value) in row.iter().zip(values.iter()) {
                acc += coeff * malachite::Integer::from(value);
            }
            i64::try_from(&acc).unwrap_or_else(|_| {
                eprintln!("[ERROR] transformed {label} coordinate does not fit in i64");
                std::process::exit(2);
            })
        })
        .collect()
}

fn transform_f64_coordinates(
    transform: &[Vec<malachite::Integer>],
    values: &[F64<Finite>],
    label: &str,
) -> Vec<F64<Finite>> {
    if transform.len() != values.len() || transform.iter().any(|row| row.len() != values.len()) {
        eprintln!("[ERROR] {label} basis transform shape does not match vector length");
        std::process::exit(2);
    }

    transform
        .iter()
        .map(|row| {
            let value = row
                .iter()
                .zip(values.iter())
                .map(|(coeff, value)| {
                    let coeff_f = coeff
                        .to_string()
                        .parse::<f64>()
                        .expect("basis transform coefficient fits in f64");
                    coeff_f * value.get()
                })
                .sum::<f64>();
            F64::<Finite>::new(value).expect("transformed coordinate is finite")
        })
        .collect()
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
            "[INFO] transforming Kähler parameters from source basis {:?} to computed basis {:?}",
            source_basis, computed_basis
        );
    }
    transform_f64_coordinates(&transform, values, "Kähler")
}

fn compute_b_field_gamma_for_o7_divisors(
    ambient_dim: usize,
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
) -> Vec<I64<Finite>> {
    if kklt_basis.len() != c_i.len() {
        eprintln!("[ERROR] KKLT basis and c_i length mismatch when computing B-field gamma");
        std::process::exit(2);
    }

    let mut gamma = vec![0_i64; ambient_dim];
    for (&divisor_idx, ci) in kklt_basis.iter().zip(c_i.iter()) {
        if divisor_idx >= ambient_dim {
            eprintln!(
                "[ERROR] KKLT divisor index {divisor_idx} exceeds ambient dimension {ambient_dim}"
            );
            std::process::exit(2);
        }
        if ci.get() == 6 {
            gamma[divisor_idx] += 1;
        }
    }

    gamma.into_iter().map(I64::<Finite>::new).collect()
}

fn transform_i64_coordinates_transpose(
    transform: &[Vec<malachite::Integer>],
    values: &[i64],
    label: &str,
) -> Vec<i64> {
    if transform.len() != values.len() || transform.iter().any(|row| row.len() != values.len()) {
        eprintln!("[ERROR] {label} basis transform shape does not match vector length");
        std::process::exit(2);
    }

    (0..values.len())
        .map(|col| {
            let mut acc = malachite::Integer::from(0);
            for (row, &value) in transform.iter().zip(values.iter()) {
                acc += &row[col] * malachite::Integer::from(value);
            }
            i64::try_from(&acc).unwrap_or_else(|_| {
                eprintln!("[ERROR] transformed {label} coordinate does not fit in i64");
                std::process::exit(2);
            })
        })
        .collect()
}

fn transform_m_flux_to_computed_basis(
    glsm: &[Vec<malachite::Integer>],
    computed_basis: &[usize],
    flux_basis: &[usize],
    values: &[i64],
    label: &str,
) -> Vec<i64> {
    if computed_basis == flux_basis {
        return values.to_vec();
    }
    let transform = basis_change_matrix(glsm, computed_basis, flux_basis).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute {label} basis transform: {e}");
        std::process::exit(2);
    });
    if !is_unimodular(&transform) {
        eprintln!("[ERROR] {label} basis transform is not unimodular");
        std::process::exit(2);
    }
    eprintln!(
        "[INFO] transforming {label} from flux basis {:?} to computed basis {:?}",
        flux_basis, computed_basis
    );
    transform_i64_coordinates(&transform, values, label)
}

fn transform_k_flux_to_computed_basis(
    glsm: &[Vec<malachite::Integer>],
    computed_basis: &[usize],
    flux_basis: &[usize],
    values: &[i64],
) -> Vec<i64> {
    if computed_basis == flux_basis {
        return values.to_vec();
    }
    let transform = basis_change_matrix(glsm, flux_basis, computed_basis).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute K basis transform: {e}");
        std::process::exit(2);
    });
    if !is_unimodular(&transform) {
        eprintln!("[ERROR] K basis transform is not unimodular");
        std::process::exit(2);
    }
    eprintln!(
        "[INFO] transforming K from flux basis {:?} to computed basis {:?}",
        flux_basis, computed_basis
    );
    transform_i64_coordinates_transpose(&transform, values, "K")
}

fn load_kklt_inputs(data_dir: Option<&str>, manifest_dir: &PathBuf) -> (Vec<I64<Pos>>, Vec<usize>) {
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
            (c_i, kklt_basis.indices)
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

fn stage_triangulation(data_dir: Option<&str>, manifest_dir: &PathBuf, t0: &Instant) -> PrimalGeom {
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

fn select_dual_basis(override_opt: Option<&BasisOverride>, computed: Vec<usize>) -> Vec<usize> {
    override_opt.map_or_else(
        || {
            eprintln!(
                "[INFO] using computed dual basis (len={}, basis={:?})",
                computed.len(),
                computed
            );
            computed
        },
        |explicit| {
            eprintln!("[INFO] using explicit dual basis from --dual-basis");
            explicit.indices.clone()
        },
    )
}

fn stage_flat_direction(
    data_dir: Option<&str>,
    manifest_dir: &PathBuf,
    geom: &PrimalGeom,
    allow_invalid_ek0: bool,
    use_mcallister_flux_basis_default: bool,
    dual_basis_override: Option<&BasisOverride>,
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
    let (_dual_heights, dual_triangulation) =
        cyrus_core::compute_frst_heights(&dual_points_vec, dual_origin_idx)
            .expect("Failed to compute dual FRST");
    if let Some(dir) = data_dir {
        validate_dual_checkpoint(&dual_polytope, &dual_triangulation, dir);
    }
    let (dual_glsm, dual_linrel, dual_basis_auto) =
        compute_glsm_and_linrels(&dual_points_vec).expect("Failed dual GLSM");
    let dual_basis = select_dual_basis(None, dual_basis_auto);
    let flux_basis = dual_basis_override.map_or_else(
        || {
            if use_mcallister_flux_basis_default {
                eprintln!("[INFO] using McAllister flux source basis [3, 4, 5, 8]");
                vec![3, 4, 5, 8]
            } else {
                eprintln!("[INFO] using computed dual basis as flux coordinate basis");
                dual_basis.clone()
            }
        },
        |basis| {
            eprintln!("[INFO] using explicit flux source basis from --dual-basis");
            basis.indices.clone()
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
    let dual_kappa = intersection_in_basis(&dual_kappa_full, &dual_basis);
    let (k_raw, m_raw) = load_flux_vectors(data_dir, manifest_dir);
    let k_raw = transform_k_flux_to_computed_basis(&dual_glsm, &dual_basis, &flux_basis, &k_raw);
    let m_raw =
        transform_m_flux_to_computed_basis(&dual_glsm, &dual_basis, &flux_basis, &m_raw, "M");
    let k_flux: Vec<I64<Finite>> = k_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let m_flux: Vec<I64<Finite>> = m_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let (flat_p, ek0_opt) = match compute_flat_direction_full(&dual_kappa, &k_flux, &m_flux) {
        Some(v) => (v.p, Some(v.ek0)),
        None if allow_invalid_ek0 => {
            eprintln!(
                "[WARN] invalid flat direction (ek0<=0); continuing due to --allow-invalid-ek0"
            );
            let p = compute_flat_direction(&dual_kappa, &k_flux, &m_flux).unwrap_or_else(|| {
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
        dual_basis,
        dual_kappa,
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
    let rays = compute_mori_cone_cap_rays(
        &flat.dual_triangulation,
        &flat.dual_triangulation_points,
        &flat.dual_polytope,
        true,
        false,
        Some(&flat.dual_basis),
    )
    .expect("Failed mori cone cap rays");
    let grading = compute_grading_vector(&rays).expect("No grading vector found");
    let curve_basis = compute_curve_basis_matrix(&flat.dual_linrels, &flat.dual_basis)
        .expect("Failed curve basis matrix");
    let q_matrix: Vec<Vec<i64>> = curve_basis
        .iter()
        .map(|row| {
            row.iter()
                .skip(1)
                .map(|v| i64::try_from(v).expect("q entry fits i64"))
                .collect()
        })
        .collect();
    let invariants = cyrus_core::compute_gv_invariants(
        &rays,
        &grading,
        &q_matrix,
        &flat.dual_kappa,
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
        eprintln!("[ERROR] no stable racetrack solution");
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
    branches_by_volume: &[cyrus_core::KkltBranchSolution],
    small_curve_cutoff: F64<Pos>,
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
        let basis_rays = project_mori_cone_cap_rays_to_basis(&ambient_rays, &intersection.basis)
            .map_err(|e| format!("failed to project primal Mori-cap rays to basis: {e}"))?;
        Some(
            compute_grading_vector(&basis_rays)
                .ok_or_else(|| "failed to compute branch GV degree grading vector".to_string())?,
        )
    } else {
        None
    };

    branches_by_volume
        .iter()
        .map(|branch| {
            let small_curve_candidates = subcutoff_toric_curve_candidates(
                &ambient_rays,
                &intersection.basis,
                &branch.result.t,
                small_curve_cutoff,
            )
            .map_err(|e| format!("failed to select branch small toric curve candidates: {e}"))?;
            let small_curves =
                remove_pair_decomposable_curve_candidates(&small_curve_candidates)
                    .map_err(|e| format!("failed to prune branch pair-decomposable curves: {e}"))?;
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
                    ) {
                        if let Some(decomposition) =
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
    basis_rays: Option<&[Vec<i64>]>,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<HashMap<Vec<i64>, malachite::Integer>, String> {
    if min_points.is_none() && max_deg.is_none() {
        return Err(
            "primal general GV fallback requires --primal-gv-min-points or --primal-gv-max-deg"
                .into(),
        );
    }

    let computed_rays;
    let rays = if let Some(basis_rays) = basis_rays {
        basis_rays
    } else {
        computed_rays = compute_mori_cone_cap_rays(
            &geom.triangulation,
            &geom.triangulation_points,
            &geom.polytope,
            true,
            false,
            Some(&intersection.basis),
        )
        .map_err(|e| format!("failed to compute primal basis Mori-cap rays: {e}"))?;
        &computed_rays
    };
    let grading = compute_grading_vector(&rays)
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
    let curve_basis = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
        .map_err(|e| format!("failed to compute primal curve basis matrix: {e}"))?;
    let q_matrix = curve_basis
        .iter()
        .map(|row| {
            row.iter()
                .skip(1)
                .map(|value| {
                    i64::try_from(value)
                        .map_err(|_| "primal q-matrix entry does not fit in i64".to_string())
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let general_gvs = cyrus_core::compute_gv_invariants(
        &rays,
        &grading,
        &q_matrix,
        &intersection.kappa_basis,
        min_points,
        max_deg,
    )
    .map_err(|e| format!("failed to compute primal general GV invariants: {e}"))?;
    let ambient_gvs = map_basis_gv_invariants_to_ambient(&general_gvs, &curve_basis)
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

fn non_positive_basis_generator_degrees(
    rays: &[Vec<i64>],
    grading: &[i64],
) -> Result<(usize, Option<(usize, i128, Vec<i64>)>), String> {
    let mut count = 0usize;
    let mut first = None;
    for (idx, ray) in rays.iter().enumerate() {
        if ray.len() != grading.len() {
            return Err(format!(
                "basis ray length {} does not match grading vector length {}",
                ray.len(),
                grading.len()
            ));
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 {
            count += 1;
            if first.is_none() {
                first = Some((idx, degree, ray.clone()));
            }
        }
    }
    Ok((count, first))
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct GradedRayStats {
    count: usize,
    degree_bounded_count: Option<usize>,
    min_degree: Option<i128>,
    max_degree: Option<i128>,
}

fn graded_ray_stats(
    rays: &[Vec<i64>],
    grading: &[i64],
    max_degree: Option<u32>,
) -> Result<GradedRayStats, String> {
    let mut bounded = max_degree.map(|_| 0usize);
    let max_degree_i128 = max_degree.map(i128::from);
    let mut min_degree = None::<i128>;
    let mut max_degree_seen = None::<i128>;
    for ray in rays {
        if ray.len() != grading.len() {
            return Err(format!(
                "basis ray length {} does not match grading vector length {}",
                ray.len(),
                grading.len()
            ));
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        min_degree = Some(min_degree.map_or(degree, |current| current.min(degree)));
        max_degree_seen = Some(max_degree_seen.map_or(degree, |current| current.max(degree)));
        if let (Some(limit), Some(count)) = (max_degree_i128, bounded.as_mut())
            && degree <= limit
        {
            *count += 1;
        }
    }
    Ok(GradedRayStats {
        count: rays.len(),
        degree_bounded_count: bounded,
        min_degree,
        max_degree: max_degree_seen,
    })
}

fn degree_filtered_basis_rays(
    rays: &[Vec<i64>],
    grading: &[i64],
    max_degree: u32,
) -> Result<Vec<Vec<i64>>, String> {
    let limit = i128::from(max_degree);
    let mut out = Vec::new();
    for ray in rays {
        if ray.len() != grading.len() {
            return Err(format!(
                "basis ray length {} does not match grading vector length {}",
                ray.len(),
                grading.len()
            ));
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= limit {
            out.push(ray.clone());
        }
    }
    Ok(out)
}

fn missing_gv_target_stats(
    ambient_classes: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    origin_idx: usize,
    origin_circuits_by_class: &HashMap<Vec<i64>, cyrus_core::OriginCircuitCurveDiagnostic>,
    cms_intersection_checks_by_class: &HashMap<Vec<i64>, Vec<CmsGeneralDivisorIntersectionCheck>>,
    sample_limit: usize,
) -> Result<MissingGvTargetStats, String> {
    if basis.len() != grading.len() {
        return Err(format!(
            "basis length {} does not match grading vector length {}",
            basis.len(),
            grading.len()
        ));
    }
    if ambient_classes.is_empty() {
        return Ok(MissingGvTargetStats {
            target_count: 0,
            targets_that_are_mori_generators: 0,
            targets_that_are_origin_circuits: 0,
            targets_real_cone_decomposable_by_other_generators: 0,
            targets_that_are_lp_extremal_mori_generators: 0,
            real_cone_decomposition_active_generator_min: None,
            real_cone_decomposition_active_generator_max: None,
            origin_circuit_resolved_conifold_count: 0,
            min_generators_le_target_degree: 0,
            max_generators_le_target_degree: 0,
            origin_coefficient_counts: BTreeMap::new(),
            origin_circuit_pattern_counts: BTreeMap::new(),
            sample: Vec::new(),
        });
    }

    let mut sorted_degrees = Vec::with_capacity(basis_rays.len());
    let mut sorted_basis_refs: Vec<&Vec<i64>> = Vec::with_capacity(basis_rays.len());
    for ray in basis_rays {
        if ray.len() != grading.len() {
            return Err(format!(
                "basis ray length {} does not match grading vector length {}",
                ray.len(),
                grading.len()
            ));
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        sorted_degrees.push(degree);
        sorted_basis_refs.push(ray);
    }
    sorted_degrees.sort_unstable();
    sorted_basis_refs.sort_unstable();

    let mut targets_that_are_mori_generators = 0usize;
    let mut targets_that_are_origin_circuits = 0usize;
    let mut targets_real_cone_decomposable_by_other_generators = 0usize;
    let mut targets_that_are_lp_extremal_mori_generators = 0usize;
    let mut real_cone_decomposition_active_generator_min = None::<usize>;
    let mut real_cone_decomposition_active_generator_max = None::<usize>;
    let mut origin_circuit_resolved_conifold_count = 0usize;
    let mut min_generators = usize::MAX;
    let mut max_generators = 0usize;
    let mut origin_coefficient_counts = BTreeMap::new();
    let mut origin_circuit_pattern_counts = BTreeMap::new();
    let mut sample = Vec::new();
    for ambient_class in ambient_classes {
        let Some(&origin_coefficient) = ambient_class.get(origin_idx) else {
            return Err(format!(
                "origin index {origin_idx} is out of bounds for ambient curve dimension {}",
                ambient_class.len()
            ));
        };
        *origin_coefficient_counts
            .entry(origin_coefficient)
            .or_insert(0) += 1;
        let basis_class = project_ambient_curve_to_basis(ambient_class, basis)?;
        let degree = basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 {
            return Err(format!(
                "missing GV target has non-positive grading degree {degree}: {ambient_class:?}"
            ));
        }
        let generators_le_degree =
            sorted_degrees.partition_point(|&ray_degree| ray_degree <= degree);
        let is_mori_generator = sorted_basis_refs
            .binary_search_by(|ray| ray.as_slice().cmp(basis_class.as_slice()))
            .is_ok();
        if is_mori_generator {
            targets_that_are_mori_generators += 1;
        }
        let real_cone_decomposition = real_cone_decomposition_by_other_degree_bounded_generators(
            &basis_class,
            basis_rays,
            grading,
            degree,
        )?;
        let real_cone_decomposable = real_cone_decomposition.is_some();
        if real_cone_decomposable {
            targets_real_cone_decomposable_by_other_generators += 1;
        }
        let real_cone_decomposition_active_generators =
            real_cone_decomposition.as_ref().map(|witness| {
                let active_count = witness.active_generator_indices.len();
                real_cone_decomposition_active_generator_min = Some(
                    real_cone_decomposition_active_generator_min
                        .map_or(active_count, |current| current.min(active_count)),
                );
                real_cone_decomposition_active_generator_max = Some(
                    real_cone_decomposition_active_generator_max
                        .map_or(active_count, |current| current.max(active_count)),
                );
                active_count
            });
        if is_mori_generator && !real_cone_decomposable {
            targets_that_are_lp_extremal_mori_generators += 1;
        }
        let origin_circuit_diagnostic = origin_circuits_by_class.get(ambient_class);
        let origin_circuit_pattern_label = origin_circuit_diagnostic.map(|diagnostic| {
            targets_that_are_origin_circuits += 1;
            if diagnostic.is_resolved_conifold_pattern {
                origin_circuit_resolved_conifold_count += 1;
            }
            let pattern = origin_circuit_pattern(diagnostic);
            *origin_circuit_pattern_counts
                .entry(pattern.clone())
                .or_insert(0) += 1;
            pattern
        });
        let origin_circuit_witness_count =
            origin_circuit_diagnostic.map(|diagnostic| diagnostic.witnesses.len());
        let origin_circuit_first_witness = origin_circuit_diagnostic
            .and_then(|diagnostic| diagnostic.witnesses.first())
            .map(origin_circuit_witness_sample);
        let cms_general_divisor_shape_candidates = origin_circuit_diagnostic
            .map(cms_general_divisor_shape_candidates)
            .filter(|candidates| !candidates.is_empty());
        let cms_general_divisor_intersection_checks = cms_intersection_checks_by_class
            .get(ambient_class)
            .filter(|checks| !checks.is_empty())
            .cloned();
        min_generators = min_generators.min(generators_le_degree);
        max_generators = max_generators.max(generators_le_degree);
        if sample.len() < sample_limit {
            sample.push(MissingGvTargetSample {
                degree,
                generators_le_degree,
                is_mori_generator,
                origin_circuit_pattern: origin_circuit_pattern_label,
                origin_circuit_witness_count,
                origin_circuit_first_witness,
                cms_general_divisor_shape_candidates,
                cms_general_divisor_intersection_checks,
                real_cone_decomposable_by_other_generators: real_cone_decomposable,
                real_cone_decomposition_active_generators,
                ambient_nonzero: ambient_class
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
                    .collect(),
                basis_nonzero: basis_class
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
                    .collect(),
            });
        }
    }

    Ok(MissingGvTargetStats {
        target_count: ambient_classes.len(),
        targets_that_are_mori_generators,
        targets_that_are_origin_circuits,
        targets_real_cone_decomposable_by_other_generators,
        targets_that_are_lp_extremal_mori_generators,
        real_cone_decomposition_active_generator_min,
        real_cone_decomposition_active_generator_max,
        origin_circuit_resolved_conifold_count,
        min_generators_le_target_degree: min_generators,
        max_generators_le_target_degree: max_generators,
        origin_coefficient_counts,
        origin_circuit_pattern_counts,
        sample,
    })
}

fn real_cone_decomposable_by_other_degree_bounded_generators(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    target_degree: i128,
) -> Result<bool, String> {
    Ok(real_cone_decomposition_by_other_degree_bounded_generators(
        target,
        basis_rays,
        grading,
        target_degree,
    )?
    .is_some())
}

fn real_cone_decomposition_by_other_degree_bounded_generators(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    target_degree: i128,
) -> Result<Option<RealConeDecompositionWitness>, String> {
    Ok(
        real_cone_decomposition_witnesses_by_other_degree_bounded_generators(
            target,
            basis_rays,
            grading,
            target_degree,
            1,
        )?
        .into_iter()
        .next(),
    )
}

fn real_cone_decomposition_witnesses_by_other_degree_bounded_generators(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    target_degree: i128,
    max_witnesses: usize,
) -> Result<Vec<RealConeDecompositionWitness>, String> {
    if target.len() != grading.len() {
        return Err(format!(
            "target length {} does not match grading vector length {}",
            target.len(),
            grading.len()
        ));
    }

    let mut candidate_indices = Vec::new();
    for (idx, ray) in basis_rays.iter().enumerate() {
        if ray.len() != target.len() {
            return Err(format!(
                "basis ray length {} does not match target length {}",
                ray.len(),
                target.len()
            ));
        }
        if ray.as_slice() == target {
            continue;
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree > 0 && degree <= target_degree {
            candidate_indices.push(idx);
        }
    }
    if candidate_indices.is_empty() {
        return Ok(Vec::new());
    }

    let objective_attempts = max_witnesses.saturating_mul(4).max(1);
    let mut witnesses = Vec::new();
    let mut seen_active_sets = HashSet::new();
    for objective_salt in 0..objective_attempts {
        let Some(witness) = solve_real_cone_decomposition_lp(
            target,
            basis_rays,
            &candidate_indices,
            target_degree,
            objective_salt as u64,
        )?
        else {
            continue;
        };
        if seen_active_sets.insert(witness.active_generator_indices.clone()) {
            witnesses.push(witness);
            if witnesses.len() >= max_witnesses {
                break;
            }
        }
    }
    Ok(witnesses)
}

fn solve_real_cone_decomposition_lp(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    candidate_indices: &[usize],
    target_degree: i128,
    objective_salt: u64,
) -> Result<Option<RealConeDecompositionWitness>, String> {
    let mut vars = ProblemVariables::new();
    let coefficients: Vec<_> = candidate_indices
        .iter()
        .map(|_| vars.add(variable().min(0.0)))
        .collect();
    let mut objective = Expression::from(0.0);
    for (coefficient, &ray_idx) in coefficients.iter().zip(candidate_indices.iter()) {
        objective.add_mul(
            decomposition_objective_weight(ray_idx, &basis_rays[ray_idx], objective_salt),
            *coefficient,
        );
    }
    let mut model = vars.minimise(objective).using(default_solver);
    for coord in 0..target.len() {
        let mut expr = Expression::from(0.0);
        for (coefficient, &ray_idx) in coefficients.iter().zip(candidate_indices.iter()) {
            let ray_coefficient = basis_rays[ray_idx][coord];
            if ray_coefficient != 0 {
                expr.add_mul(ray_coefficient as f64, *coefficient);
            }
        }
        model = model.with(expr.eq(target[coord] as f64));
    }

    let solution = match model.solve() {
        Ok(solution) => solution,
        Err(ResolutionError::Infeasible) => return Ok(None),
        Err(err) => {
            return Err(format!(
                "failed real-cone decomposition LP for target degree {target_degree}: {err}"
            ));
        }
    };

    let mut max_residual = 0.0f64;
    for coord in 0..target.len() {
        let reconstructed = coefficients
            .iter()
            .zip(candidate_indices.iter())
            .map(|(coefficient, &ray_idx)| {
                solution.value(*coefficient) * basis_rays[ray_idx][coord] as f64
            })
            .sum::<f64>();
        max_residual = max_residual.max((reconstructed - target[coord] as f64).abs());
    }
    if max_residual > 1.0e-6 {
        return Err(format!(
            "real-cone decomposition LP returned residual {max_residual}"
        ));
    }

    let mut active_generator_indices = Vec::new();
    for (coefficient, &ray_idx) in coefficients.iter().zip(candidate_indices.iter()) {
        if solution.value(*coefficient).abs() > 1.0e-8 {
            active_generator_indices.push(ray_idx);
        }
    }
    active_generator_indices.sort_unstable();
    if active_generator_indices.is_empty() {
        return Err(
            "real-cone decomposition LP returned no active generators despite zero residual"
                .to_string(),
        );
    }
    Ok(Some(RealConeDecompositionWitness {
        active_generator_indices,
    }))
}

fn decomposition_objective_weight(ray_idx: usize, ray: &[i64], objective_salt: u64) -> f64 {
    if objective_salt == 0 {
        return 1.0;
    }
    let mut hash = 0x9e37_79b9_7f4a_7c15u64 ^ objective_salt.rotate_left(17) ^ ray_idx as u64;
    for (idx, &value) in ray.iter().enumerate().filter(|(_, value)| **value != 0) {
        hash ^= (idx as u64).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        hash = hash.rotate_left(11) ^ (value as u64).wrapping_mul(0x94d0_49bb_1331_11ebu64);
    }
    1.0 + (hash % 1_000_003) as f64 / 1_000_003.0
}

fn one_dimensional_ray_gv_targets(
    ambient_classes: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
) -> Result<OneDimensionalRayGvTargets, String> {
    let basis_ray_set: HashSet<Vec<i64>> = basis_rays.iter().cloned().collect();
    let mut candidates = Vec::new();
    let mut skipped_non_generators = 0usize;
    let mut skipped_decomposable_generators = 0usize;
    for ambient_class in ambient_classes {
        let basis_class = project_ambient_curve_to_basis(ambient_class, basis)?;
        let degree = basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 {
            return Err(format!(
                "one-dimensional ray GV target has non-positive grading degree {degree}: {ambient_class:?}"
            ));
        }
        if !basis_ray_set.contains(&basis_class) {
            skipped_non_generators += 1;
            continue;
        }
        if real_cone_decomposable_by_other_degree_bounded_generators(
            &basis_class,
            basis_rays,
            grading,
            degree,
        )? {
            skipped_decomposable_generators += 1;
            continue;
        }
        candidates.push(ambient_class.clone());
    }

    Ok(OneDimensionalRayGvTargets {
        candidates,
        skipped_non_generators,
        skipped_decomposable_generators,
    })
}

fn origin_circuit_pattern(diagnostic: &cyrus_core::OriginCircuitCurveDiagnostic) -> String {
    format!(
        "origin={};neg={:?};pos={:?};resolved_conifold={}",
        diagnostic.origin_coefficient,
        diagnostic.negative_coefficient_counts,
        diagnostic.positive_coefficient_counts,
        diagnostic.is_resolved_conifold_pattern
    )
}

fn origin_circuit_witness_sample(
    witness: &cyrus_core::OriginCircuitCurveWitness,
) -> OriginCircuitWitnessSample {
    OriginCircuitWitnessSample {
        first_facet_exclusive_point: witness.first_facet_exclusive_point,
        second_facet_exclusive_point: witness.second_facet_exclusive_point,
        shared_two_simplex: witness.shared_two_simplex.clone(),
        first_facet_size: witness.first_facet.len(),
        second_facet_size: witness.second_facet.len(),
        sparse_relation: witness.sparse_relation.clone(),
        relation_points: witness
            .relation_points
            .iter()
            .map(|point| OriginCircuitRelationPointSample {
                point_index: point.point_index,
                coefficient: point.coefficient,
                coordinates: point.coordinates.clone(),
                face_dimension: point.face_dimension,
            })
            .collect(),
    }
}

fn cms_general_divisor_shape_candidates(
    diagnostic: &cyrus_core::OriginCircuitCurveDiagnostic,
) -> Vec<CmsGeneralDivisorShapeCandidate> {
    let mut candidates: Vec<_> = diagnostic
        .witnesses
        .iter()
        .flat_map(cms_general_divisor_shape_candidates_for_witness)
        .collect();
    candidates.sort();
    candidates.dedup();
    candidates
}

fn cms_general_divisor_intersection_checks_by_class(
    ambient_classes: &[Vec<i64>],
    origin_circuits_by_class: &HashMap<Vec<i64>, cyrus_core::OriginCircuitCurveDiagnostic>,
    kappa_full: &cyrus_core::Intersection,
    basis: &[usize],
) -> Result<HashMap<Vec<i64>, Vec<CmsGeneralDivisorIntersectionCheck>>, String> {
    let mut out = HashMap::new();
    for ambient_class in ambient_classes {
        let Some(diagnostic) = origin_circuits_by_class.get(ambient_class) else {
            continue;
        };
        let mut checks = Vec::new();
        for candidate in cms_general_divisor_shape_candidates(diagnostic) {
            checks.push(cms_general_divisor_intersection_check(
                ambient_class,
                &candidate,
                kappa_full,
                basis,
            )?);
        }
        checks.sort();
        checks.dedup();
        if !checks.is_empty() {
            out.insert(ambient_class.clone(), checks);
        }
    }
    Ok(out)
}

fn cms_general_divisor_intersection_check(
    ambient_class: &[i64],
    candidate: &CmsGeneralDivisorShapeCandidate,
    kappa_full: &cyrus_core::Intersection,
    basis: &[usize],
) -> Result<CmsGeneralDivisorIntersectionCheck, String> {
    let i = candidate.shrinking_divisor_index;
    let Some(solution) = divisor_representation_solution(ambient_class, i, kappa_full, basis)?
    else {
        return Ok(CmsGeneralDivisorIntersectionCheck {
            shrinking_divisor_index: i,
            has_rational_divisor_solution: false,
            solution_basis_support_len: None,
            solution_is_integral: None,
            computed_other_normal_degree: None,
            matches_inferred_other_normal_degree: None,
        });
    };
    let support_len = solution
        .coefficients
        .iter()
        .filter(|value| **value != 0)
        .count();
    let solution_is_integral = solution.coefficients.iter().all(rational_is_integer);
    let inferred_m = malachite::Rational::from(candidate.inferred_other_normal_degree);
    Ok(CmsGeneralDivisorIntersectionCheck {
        shrinking_divisor_index: i,
        has_rational_divisor_solution: true,
        solution_basis_support_len: Some(support_len),
        solution_is_integral: Some(solution_is_integral),
        computed_other_normal_degree: Some(solution.other_normal_degree.to_string()),
        matches_inferred_other_normal_degree: Some(solution.other_normal_degree == inferred_m),
    })
}

fn compute_covered_toric_gv_divisor_representation_baseline(
    covered_classes: &[(Vec<i64>, malachite::Integer)],
    origin_idx: usize,
    kappa_full: &cyrus_core::Intersection,
    basis: &[usize],
    class_limit: usize,
) -> Result<Option<CoveredToricGvDivisorRepresentationBaseline>, String> {
    if covered_classes.is_empty() {
        return Ok(None);
    }
    if class_limit == 0 {
        return Ok(None);
    }

    let mut support_index_checks = 0usize;
    let mut classes_with_support_divisor_solution = 0usize;
    let mut first_without_support_divisor_solution = None;
    let mut sample = Vec::new();
    let mut checked_class_count = 0usize;
    for (ambient_class, _) in covered_classes.iter().take(class_limit) {
        checked_class_count += 1;
        let mut support_indices_with_solution = Vec::new();
        let mut support_indices_without_solution = Vec::new();
        for (point_idx, &coefficient) in ambient_class.iter().enumerate() {
            if coefficient == 0 || point_idx == origin_idx {
                continue;
            }
            support_index_checks += 1;
            if divisor_representation_solution(ambient_class, point_idx, kappa_full, basis)?
                .is_some()
            {
                support_indices_with_solution.push(point_idx);
            } else {
                support_indices_without_solution.push(point_idx);
            }
        }
        if support_indices_with_solution.is_empty() {
            first_without_support_divisor_solution.get_or_insert_with(|| sparse_i64(ambient_class));
        } else {
            classes_with_support_divisor_solution += 1;
        }
        if sample.len() < 8 {
            sample.push(CoveredToricGvDivisorRepresentationSample {
                ambient_nonzero: sparse_i64(ambient_class),
                support_indices_with_solution,
                support_indices_without_solution,
            });
        }
    }

    Ok(Some(CoveredToricGvDivisorRepresentationBaseline {
        checked_class_count,
        class_limit,
        support_index_checks,
        classes_with_support_divisor_solution,
        first_without_support_divisor_solution,
        sample,
    }))
}

fn covered_toric_gv_divisor_representation_class_limit() -> usize {
    std::env::var("CYRUS_CORRECTED_CHAMBER_COVERED_GV_DIVISOR_REPRESENTATION_CLASS_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_CORRECTED_CHAMBER_COVERED_GV_DIVISOR_REPRESENTATION_CLASS_LIMIT)
}

fn divisor_representation_solution(
    ambient_class: &[i64],
    point_idx: usize,
    kappa_full: &cyrus_core::Intersection,
    basis: &[usize],
) -> Result<Option<DivisorRepresentationSolution>, String> {
    if ambient_class.len() != kappa_full.dim() {
        return Err(format!(
            "divisor representation curve dimension {} does not match intersection dimension {}",
            ambient_class.len(),
            kappa_full.dim()
        ));
    }
    if point_idx >= kappa_full.dim() {
        return Err(format!(
            "divisor representation point index {point_idx} is out of bounds for dimension {}",
            kappa_full.dim()
        ));
    }
    for &basis_idx in basis {
        if basis_idx >= kappa_full.dim() {
            return Err(format!(
                "divisor representation basis index {basis_idx} is out of bounds for dimension {}",
                kappa_full.dim()
            ));
        }
    }

    let matrix: Vec<Vec<malachite::Rational>> = (0..kappa_full.dim())
        .map(|row_idx| {
            basis
                .iter()
                .map(|&basis_idx| kappa_full.get(row_idx, point_idx, basis_idx).get().clone())
                .collect()
        })
        .collect();
    let rhs: Vec<malachite::Rational> = ambient_class
        .iter()
        .map(|&coefficient| malachite::Rational::from(coefficient))
        .collect();
    let Some(coefficients) = solve_rational_linear_system(&matrix, &rhs)? else {
        return Ok(None);
    };
    let other_normal_degree =
        divisor_square_against_point(kappa_full, point_idx, basis, &coefficients)?;
    Ok(Some(DivisorRepresentationSolution {
        coefficients,
        other_normal_degree,
    }))
}

fn divisor_square_against_point(
    kappa_full: &cyrus_core::Intersection,
    point_idx: usize,
    basis: &[usize],
    divisor_coefficients: &[malachite::Rational],
) -> Result<malachite::Rational, String> {
    if basis.len() != divisor_coefficients.len() {
        return Err(format!(
            "CMS divisor square basis length {} does not match coefficient length {}",
            basis.len(),
            divisor_coefficients.len()
        ));
    }
    let mut value = malachite::Rational::from(0);
    for (a_pos, &a_idx) in basis.iter().enumerate() {
        for (b_pos, &b_idx) in basis.iter().enumerate() {
            let term = &divisor_coefficients[a_pos]
                * &divisor_coefficients[b_pos]
                * kappa_full.get(point_idx, a_idx, b_idx).get();
            value += term;
        }
    }
    Ok(value)
}

fn solve_rational_linear_system(
    matrix: &[Vec<malachite::Rational>],
    rhs: &[malachite::Rational],
) -> Result<Option<Vec<malachite::Rational>>, String> {
    if matrix.len() != rhs.len() {
        return Err(format!(
            "linear system row count {} does not match rhs length {}",
            matrix.len(),
            rhs.len()
        ));
    }
    let Some(first_row) = matrix.first() else {
        return Ok(Some(Vec::new()));
    };
    let n_cols = first_row.len();
    if matrix.iter().any(|row| row.len() != n_cols) {
        return Err("linear system matrix rows have inconsistent lengths".to_string());
    }

    let mut work: Vec<Vec<malachite::Rational>> = matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, rhs_value)| {
            let mut augmented = row.clone();
            augmented.push(rhs_value.clone());
            augmented
        })
        .collect();
    let mut pivot_cols = Vec::new();
    let mut pivot_row = 0usize;
    for col in 0..n_cols {
        let Some(found_row) = (pivot_row..work.len()).find(|&row| work[row][col] != 0) else {
            continue;
        };
        work.swap(pivot_row, found_row);
        let pivot = work[pivot_row][col].clone();
        for entry in &mut work[pivot_row][col..=n_cols] {
            *entry /= &pivot;
        }
        let normalized_pivot_row = work[pivot_row].clone();
        for row in 0..work.len() {
            if row == pivot_row || work[row][col] == 0 {
                continue;
            }
            let factor = work[row][col].clone();
            for c in col..=n_cols {
                let sub = &factor * &normalized_pivot_row[c];
                work[row][c] -= sub;
            }
        }
        pivot_cols.push(col);
        pivot_row += 1;
        if pivot_row == work.len() {
            break;
        }
    }

    for row in pivot_row..work.len() {
        let all_zero = work[row][..n_cols].iter().all(|value| *value == 0);
        if all_zero && work[row][n_cols] != 0 {
            return Ok(None);
        }
    }

    let mut solution = vec![malachite::Rational::from(0); n_cols];
    for (row, &col) in pivot_cols.iter().enumerate() {
        solution[col] = work[row][n_cols].clone();
    }
    Ok(Some(solution))
}

fn rational_is_integer(value: &malachite::Rational) -> bool {
    value.denominator_ref() == &1u32
}

fn cms_general_divisor_shape_candidates_for_witness(
    witness: &cyrus_core::OriginCircuitCurveWitness,
) -> Vec<CmsGeneralDivisorShapeCandidate> {
    let all_non_origin_relation_points_are_two_face = witness
        .relation_points
        .iter()
        .filter(|point| point.point_index != 0)
        .all(|point| point.face_dimension == Some(2));
    witness
        .relation_points
        .iter()
        .filter(|point| point.point_index != 0 && point.coefficient < 0)
        .map(|point| {
            let inferred_other_normal_degree = -2 - point.coefficient;
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: point.point_index,
                shrinking_divisor_coefficient: point.coefficient,
                shrinking_divisor_coordinates: point.coordinates.clone(),
                inferred_other_normal_degree,
                toric_gv1_formula_value: toric_gv1_formula_value(inferred_other_normal_degree),
                all_non_origin_relation_points_are_two_face,
            }
        })
        .collect()
}

fn toric_gv1_formula_value(m: i64) -> Option<i64> {
    let magnitude = m + 2;
    if magnitude < 0 {
        return None;
    }
    let sign = if (m + 1).rem_euclid(2) == 0 { 1 } else { -1 };
    Some(sign * magnitude)
}

fn project_ambient_curve_to_basis(
    ambient_class: &[i64],
    basis: &[usize],
) -> Result<Vec<i64>, String> {
    basis
        .iter()
        .map(|&idx| {
            ambient_class.get(idx).copied().ok_or_else(|| {
                format!(
                    "basis index {idx} is out of bounds for ambient curve dimension {}",
                    ambient_class.len()
                )
            })
        })
        .collect()
}

fn compute_missing_one_dimensional_ray_gvs(
    ambient_classes: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
) -> Result<Vec<(Vec<i64>, malachite::Integer, i128)>, String> {
    let mut out = Vec::with_capacity(ambient_classes.len());
    for ambient_class in ambient_classes {
        let basis_class = project_ambient_curve_to_basis(ambient_class, basis)?;
        let degree = basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 {
            return Err(format!(
                "one-dimensional ray GV target has non-positive grading degree {degree}: {ambient_class:?}"
            ));
        }
        let max_deg = u32::try_from(degree).map_err(|_| {
            format!("one-dimensional ray GV target degree {degree} does not fit in u32")
        })?;
        let target_i32 = basis_class
            .iter()
            .map(|&value| {
                i32::try_from(value).map_err(|_| {
                    "one-dimensional ray GV target coordinate does not fit in i32".to_string()
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let previous_panic_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            cyrus_core::compute_gv_invariants_with_provided_generators(
                std::slice::from_ref(&basis_class),
                grading,
                q_matrix,
                intnums,
                None,
                Some(max_deg),
            )
        }));
        std::panic::set_hook(previous_panic_hook);
        let gvs = match gvs_result {
            Ok(Ok(gvs)) => gvs,
            Ok(Err(e)) => return Err(format!("failed one-dimensional ray GV computation: {e}")),
            Err(payload) => {
                return Err(format!(
                    "one-dimensional ray GV computation panicked for target degree {degree}, ambient_nonzero={:?}; this indicates the single-ray truncation is inconsistent for this target ({})",
                    sparse_i64(ambient_class),
                    panic_payload_message(payload.as_ref())
                ));
            }
        };
        let gv = gvs
            .into_iter()
            .find_map(|(curve, gv)| (curve == target_i32).then_some(gv))
            .unwrap_or_else(|| malachite::Integer::from(0));
        out.push((ambient_class.clone(), gv, degree));
    }
    Ok(out)
}

fn compute_missing_lp_witness_face_gvs(
    ambient_classes: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
) -> Result<Vec<FaceGvDiagnosticResult>, String> {
    let mut out = Vec::new();
    for ambient_class in ambient_classes {
        let basis_class = project_ambient_curve_to_basis(ambient_class, basis)?;
        let degree = basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 {
            return Err(format!(
                "LP-witness face GV target has non-positive grading degree {degree}: {ambient_class:?}"
            ));
        }
        let witnesses = real_cone_decomposition_witnesses_by_other_degree_bounded_generators(
            &basis_class,
            basis_rays,
            grading,
            degree,
            8,
        )?;
        if witnesses.is_empty() {
            continue;
        }
        let max_deg = u32::try_from(degree).map_err(|_| {
            format!("LP-witness face GV target degree {degree} does not fit in u32")
        })?;
        let mut selected_attempt = None;
        let mut errors = Vec::new();
        for witness in &witnesses {
            let attempt = compute_lp_witness_face_attempt(
                witness,
                &basis_class,
                basis_rays,
                grading,
                q_matrix,
                intnums,
                degree,
                max_deg,
                ambient_class,
            )?;
            if attempt.gv.is_some() {
                selected_attempt = Some(attempt);
                break;
            }
            if let Some(error) = attempt.error.as_ref() {
                errors.push(error.clone());
            }
            if selected_attempt.is_none() {
                selected_attempt = Some(attempt);
            }
        }
        let mut attempt = selected_attempt.expect("at least one LP witness was available");
        if attempt.gv.is_none() && !errors.is_empty() {
            attempt.error = Some(format!(
                "{} LP witness attempts failed: {}",
                witnesses.len(),
                errors.join(" | ")
            ));
        }

        out.push(FaceGvDiagnosticResult {
            ambient_class: ambient_class.clone(),
            degree,
            generator_count: attempt.generator_count,
            active_generator_count: attempt.active_generator_count,
            span_generator_count: attempt.span_generator_count,
            used_span_expansion: attempt.used_span_expansion,
            used_lattice_saturation: attempt.used_lattice_saturation,
            used_integer_decomposition: attempt.used_integer_decomposition,
            used_decomposition_diamond: attempt.used_decomposition_diamond,
            used_decomposition_closure: attempt.used_decomposition_closure,
            integer_decomposition_term_count: attempt.integer_decomposition_term_count,
            lattice_semigroup_element_count: attempt.lattice_semigroup_element_count,
            supporting_face_certificate: attempt.supporting_face_certificate,
            gv: attempt.gv,
            error: attempt.error,
        });
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn compute_lp_witness_face_attempt(
    witness: &RealConeDecompositionWitness,
    basis_class: &[i64],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    degree: i128,
    max_deg: u32,
    ambient_class: &[i64],
) -> Result<FaceGvAttempt, String> {
    let mut provided_generators: Vec<Vec<i64>> = witness
        .active_generator_indices
        .iter()
        .map(|&idx| basis_rays[idx].clone())
        .collect();
    provided_generators.push(basis_class.to_vec());
    provided_generators.sort();
    provided_generators.dedup();
    let span_generators =
        degree_bounded_span_generators(&provided_generators, basis_rays, grading, degree)?;
    let span_generator_count = span_generators.len();
    let mut used_span_expansion = false;
    let mut used_lattice_saturation = false;
    let mut used_integer_decomposition = false;
    let mut used_decomposition_diamond = false;
    let mut used_decomposition_closure = false;
    let mut integer_decomposition_term_count = None;
    let mut lattice_semigroup_element_count = None;
    let mut supporting_face_certificate = None;
    let (gv, error) = match compute_provided_generator_target_gv(
        &provided_generators,
        basis_class,
        grading,
        q_matrix,
        intnums,
        max_deg,
        "LP-witness face",
        ambient_class,
    ) {
        Ok(gv) => (Some(gv), None),
        Err(first_error) => {
            let span_limit = std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_SPAN_GENERATOR_LIMIT")
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_SPAN_GENERATOR_LIMIT);
            let (lattice_seed_generators, accumulated_error) = if span_generators.len()
                > provided_generators.len()
                && span_generators.len() <= span_limit
            {
                used_span_expansion = true;
                match compute_provided_generator_target_gv(
                    &span_generators,
                    basis_class,
                    grading,
                    q_matrix,
                    intnums,
                    max_deg,
                    "span-expanded LP-witness face",
                    ambient_class,
                ) {
                    Ok(gv) => {
                        let (supporting_face_certificate, error) =
                            match certify_supporting_mori_face_if_requested(
                                &span_generators,
                                basis_rays,
                            ) {
                                Ok(certificate) => (certificate, None),
                                Err(certificate_error) => (
                                    None,
                                    Some(format!(
                                        "supporting-face certificate failed: {certificate_error}"
                                    )),
                                ),
                            };
                        return Ok(FaceGvAttempt {
                            generator_count: provided_generators.len(),
                            active_generator_count: witness.active_generator_indices.len(),
                            span_generator_count,
                            used_span_expansion,
                            used_lattice_saturation,
                            used_integer_decomposition,
                            used_decomposition_diamond,
                            used_decomposition_closure,
                            integer_decomposition_term_count,
                            lattice_semigroup_element_count,
                            supporting_face_certificate,
                            gv: Some(gv),
                            error,
                        });
                    }
                    Err(second_error) => (
                        span_generators.as_slice(),
                        format!("{first_error}; span-expanded retry also failed: {second_error}"),
                    ),
                }
            } else if span_generators.len() > span_limit {
                (
                    provided_generators.as_slice(),
                    format!(
                        "{first_error}; span-expanded retry skipped because {} generators exceed limit {span_limit}",
                        span_generators.len()
                    ),
                )
            } else {
                (provided_generators.as_slice(), first_error)
            };

            let lattice_generator_limit =
                std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_LATTICE_GENERATOR_LIMIT")
                    .ok()
                    .and_then(|value| value.parse::<usize>().ok())
                    .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_LATTICE_GENERATOR_LIMIT);
            if lattice_seed_generators.len() <= lattice_generator_limit {
                match compute_lattice_saturated_face_target_gv(
                    lattice_seed_generators,
                    basis_class,
                    grading,
                    q_matrix,
                    intnums,
                    max_deg,
                    ambient_class,
                ) {
                    Ok((gv, element_count)) => {
                        used_lattice_saturation = true;
                        lattice_semigroup_element_count = Some(element_count);
                        (Some(gv), None)
                    }
                    Err(lattice_error) => {
                        let exact_decomposition_max_terms = std::env::var(
                            "CYRUS_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_TERMS",
                        )
                        .ok()
                        .and_then(|value| value.parse::<usize>().ok())
                        .unwrap_or(
                            DEFAULT_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_TERMS,
                        );
                        match compute_integer_decomposition_face_target_gv(
                            basis_rays,
                            basis_class,
                            grading,
                            q_matrix,
                            intnums,
                            max_deg,
                            ambient_class,
                            exact_decomposition_max_terms,
                            std::env::var(
                                "CYRUS_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_WITNESSES",
                            )
                            .ok()
                            .and_then(|value| value.parse::<usize>().ok())
                            .unwrap_or(
                                DEFAULT_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_WITNESSES,
                            ),
                        ) {
                            Ok((
                                gv,
                                element_count,
                                term_count,
                                used_diamond,
                                used_closure,
                            )) => {
                                used_lattice_saturation = !used_diamond && !used_closure;
                                used_integer_decomposition = true;
                                used_decomposition_diamond = used_diamond;
                                used_decomposition_closure = used_closure;
                                integer_decomposition_term_count = Some(term_count);
                                lattice_semigroup_element_count = Some(element_count);
                                (Some(gv), None)
                            }
                            Err(exact_error) => (
                                None,
                                Some(format!(
                                    "{accumulated_error}; lattice-saturated retry also failed: {lattice_error}; exact-integer-decomposition retry also failed: {exact_error}"
                                )),
                            ),
                        }
                    }
                }
            } else {
                (
                    None,
                    Some(format!(
                        "{accumulated_error}; lattice-saturated retry skipped because {} seed generators exceed limit {lattice_generator_limit}",
                        lattice_seed_generators.len()
                    )),
                )
            }
        }
    };
    let error = if gv.is_some() {
        match certify_supporting_mori_face_if_requested(&span_generators, basis_rays) {
            Ok(certificate) => {
                supporting_face_certificate = certificate;
                error
            }
            Err(certificate_error) => match error {
                Some(existing_error) => Some(format!(
                    "{existing_error}; supporting-face certificate failed: {certificate_error}"
                )),
                None => Some(format!(
                    "supporting-face certificate failed: {certificate_error}"
                )),
            },
        }
    } else {
        error
    };

    Ok(FaceGvAttempt {
        generator_count: provided_generators.len(),
        active_generator_count: witness.active_generator_indices.len(),
        span_generator_count,
        used_span_expansion,
        used_lattice_saturation,
        used_integer_decomposition,
        used_decomposition_diamond,
        used_decomposition_closure,
        integer_decomposition_term_count,
        lattice_semigroup_element_count,
        supporting_face_certificate,
        gv,
        error,
    })
}

fn certify_supporting_mori_face_if_requested(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
) -> Result<Option<SupportingFaceCertificate>, String> {
    if !supporting_face_certificate_requested() {
        return Ok(None);
    }
    certify_supporting_mori_face(face_generators, basis_rays)
}

fn supporting_face_certificate_requested() -> bool {
    std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE")
        .map(|value| value != "0")
        .unwrap_or(false)
}

fn certify_supporting_mori_face(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
) -> Result<Option<SupportingFaceCertificate>, String> {
    if face_generators.is_empty() {
        return Ok(None);
    }
    if basis_rays.is_empty() {
        return Err("supporting-face certificate requires Mori generators".to_string());
    }
    let dim = face_generators[0].len();
    if dim == 0 {
        return Err("supporting-face certificate dimension is zero".to_string());
    }
    if face_generators.iter().any(|row| row.len() != dim)
        || basis_rays.iter().any(|row| row.len() != dim)
    {
        return Err("supporting-face certificate row dimensions are inconsistent".to_string());
    }

    let ray_limit = std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_RAY_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_RAY_LIMIT);
    if basis_rays.len() > ray_limit {
        return Err(format!(
            "supporting-face certificate skipped because {} Mori generators exceed limit {ray_limit}",
            basis_rays.len()
        ));
    }

    let anchors = supporting_face_anchor_candidates(face_generators, basis_rays)?;
    for anchor in anchors {
        let Some(lp_normal) =
            solve_supporting_face_normal_lp(face_generators, basis_rays, &anchor)?
        else {
            continue;
        };
        if let Some(certificate) =
            integer_supporting_face_certificate_from_lp(&lp_normal, face_generators, basis_rays)?
        {
            return Ok(Some(certificate));
        }
    }
    Ok(None)
}

fn supporting_face_anchor_candidates(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
) -> Result<Vec<Vec<i64>>, String> {
    let dim = face_generators[0].len();
    let face_rank = rational_rank_i64_rows(face_generators);
    if face_rank >= dim {
        return Ok(Vec::new());
    }

    let attempt_limit =
        std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_ANCHOR_ATTEMPTS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_ANCHOR_ATTEMPTS);
    let mut anchors = Vec::new();
    let mut seen = HashSet::new();
    for ray in basis_rays {
        if face_generators.iter().any(|generator| generator == ray) || !seen.insert(ray.clone()) {
            continue;
        }
        let mut rows = face_generators.to_vec();
        rows.push(ray.clone());
        if rational_rank_i64_rows(&rows) > face_rank {
            anchors.push(ray.clone());
            if anchors.len() >= attempt_limit {
                break;
            }
        }
    }
    Ok(anchors)
}

fn solve_supporting_face_normal_lp(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    anchor: &[i64],
) -> Result<Option<Vec<f64>>, String> {
    let max_rounds = std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_CUTTING_ROUNDS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_CUTTING_ROUNDS);
    let mut enforced_ray_indices = Vec::new();
    let mut enforced_ray_set = HashSet::new();
    for _ in 0..max_rounds {
        let Some(normal) = solve_supporting_face_normal_lp_with_enforced_rays(
            face_generators,
            basis_rays,
            anchor,
            &enforced_ray_indices,
        )?
        else {
            return Ok(None);
        };
        let Some(violating_idx) = most_negative_lp_normal_violation(&normal, basis_rays) else {
            return Ok(Some(normal));
        };
        if !enforced_ray_set.insert(violating_idx) {
            return Ok(None);
        }
        enforced_ray_indices.push(violating_idx);
    }
    Ok(None)
}

fn solve_supporting_face_normal_lp_with_enforced_rays(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    anchor: &[i64],
    enforced_ray_indices: &[usize],
) -> Result<Option<Vec<f64>>, String> {
    let dim = anchor.len();
    let mut vars = ProblemVariables::new();
    let bound = 1.0e9;
    let normal_vars: Vec<_> = (0..dim)
        .map(|_| vars.add(variable().min(-bound).max(bound)))
        .collect();

    let mut objective = Expression::from(0.0);
    objective.add_mul(0.0, normal_vars[0]);
    let mut model = vars.minimise(objective).using(default_solver);

    for generator in face_generators {
        let mut expr = Expression::from(0.0);
        for (var, &coefficient) in normal_vars.iter().zip(generator) {
            if coefficient != 0 {
                expr.add_mul(coefficient as f64, *var);
            }
        }
        model = model.with(expr.eq(0.0));
    }

    let mut anchor_expr = Expression::from(0.0);
    for (var, &coefficient) in normal_vars.iter().zip(anchor) {
        if coefficient != 0 {
            anchor_expr.add_mul(coefficient as f64, *var);
        }
    }
    model = model.with(anchor_expr.eq(1.0));

    for &ray_idx in enforced_ray_indices {
        let ray = basis_rays.get(ray_idx).ok_or_else(|| {
            format!("supporting-face enforced ray index {ray_idx} is out of bounds")
        })?;
        let mut expr = Expression::from(0.0);
        for (var, &coefficient) in normal_vars.iter().zip(ray) {
            if coefficient != 0 {
                expr.add_mul(coefficient as f64, *var);
            }
        }
        model = model.with(expr.geq(0.0));
    }

    let solution = match model.solve() {
        Ok(solution) => solution,
        Err(ResolutionError::Infeasible) => return Ok(None),
        Err(err) => {
            return Err(format!("supporting-face normal LP failed: {err}"));
        }
    };
    let normal = normal_vars
        .iter()
        .map(|var| solution.value(*var))
        .collect::<Vec<_>>();
    if normal.iter().all(|value| value.is_finite()) {
        Ok(Some(normal))
    } else {
        Err("supporting-face normal LP returned a non-finite value".to_string())
    }
}

fn most_negative_lp_normal_violation(lp_normal: &[f64], basis_rays: &[Vec<i64>]) -> Option<usize> {
    let mut worst_idx = None;
    let mut worst_dot = -1.0e-7;
    for (idx, ray) in basis_rays.iter().enumerate() {
        let dot = lp_normal
            .iter()
            .zip(ray)
            .map(|(&normal_coeff, &ray_coeff)| normal_coeff * ray_coeff as f64)
            .sum::<f64>();
        if dot < worst_dot {
            worst_dot = dot;
            worst_idx = Some(idx);
        }
    }
    worst_idx
}

fn integer_supporting_face_certificate_from_lp(
    lp_normal: &[f64],
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
) -> Result<Option<SupportingFaceCertificate>, String> {
    let scale_limit = std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_SCALE_LIMIT")
        .ok()
        .and_then(|value| value.parse::<i64>().ok())
        .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_SCALE_LIMIT);
    if scale_limit <= 0 {
        return Err("supporting-face certificate scale limit must be positive".to_string());
    }

    let mut seen_normals = HashSet::new();
    for scale in 1..=scale_limit {
        let Some(normal) = rounded_reduced_i64_normal(lp_normal, scale)? else {
            continue;
        };
        if !seen_normals.insert(normal.clone()) {
            continue;
        }
        if !normal_vanishes_on_generators(&normal, face_generators) {
            continue;
        }
        if let Some(certificate) = check_integer_supporting_face_normal(&normal, basis_rays)? {
            return Ok(Some(certificate));
        }
    }
    Ok(None)
}

fn rounded_reduced_i64_normal(lp_normal: &[f64], scale: i64) -> Result<Option<Vec<i64>>, String> {
    let mut normal = Vec::with_capacity(lp_normal.len());
    for &value in lp_normal {
        let scaled = value * scale as f64;
        if !scaled.is_finite() || scaled < i64::MIN as f64 || scaled > i64::MAX as f64 {
            return Err("supporting-face LP normal does not fit in i64 after scaling".to_string());
        }
        normal.push(scaled.round() as i64);
    }
    reduce_i64_vector_preserve_sign(&normal)
}

fn reduce_i64_vector_preserve_sign(values: &[i64]) -> Result<Option<Vec<i64>>, String> {
    let mut gcd = 0i64;
    for &value in values {
        if value == i64::MIN {
            return Err("supporting-face normal coefficient is i64::MIN".to_string());
        }
        gcd = gcd_i64(gcd, value.abs());
    }
    if gcd == 0 {
        return Ok(None);
    }
    values
        .iter()
        .map(|&value| {
            value
                .checked_div(gcd)
                .ok_or_else(|| "supporting-face normal reduction divided by zero".to_string())
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

fn normal_vanishes_on_generators(normal: &[i64], generators: &[Vec<i64>]) -> bool {
    generators
        .iter()
        .all(|generator| exact_i64_dot(normal, generator) == 0)
}

fn check_integer_supporting_face_normal(
    normal: &[i64],
    basis_rays: &[Vec<i64>],
) -> Result<Option<SupportingFaceCertificate>, String> {
    if normal.iter().all(|&value| value == 0) {
        return Ok(None);
    }
    let mut zero_generator_count = 0usize;
    let mut positive_generator_count = 0usize;
    for ray in basis_rays {
        if ray.len() != normal.len() {
            return Err(format!(
                "supporting-face normal dimension {} does not match Mori ray dimension {}",
                normal.len(),
                ray.len()
            ));
        }
        let dot = exact_i64_dot(normal, ray);
        if dot < 0 {
            return Ok(None);
        }
        if dot == 0 {
            zero_generator_count += 1;
        } else {
            positive_generator_count += 1;
        }
    }
    if positive_generator_count == 0 {
        return Ok(None);
    }
    Ok(Some(SupportingFaceCertificate {
        zero_generator_count,
        positive_generator_count,
        normal: normal.to_vec(),
    }))
}

fn exact_i64_dot(lhs: &[i64], rhs: &[i64]) -> i128 {
    lhs.iter()
        .zip(rhs)
        .map(|(&left, &right)| i128::from(left) * i128::from(right))
        .sum()
}

fn gcd_i64(a: i64, b: i64) -> i64 {
    if b == 0 { a.abs() } else { gcd_i64(b, a % b) }
}

#[allow(clippy::too_many_arguments)]
fn compute_provided_generator_target_gv(
    provided_generators: &[Vec<i64>],
    target_class: &[i64],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    max_deg: u32,
    label: &str,
    ambient_class: &[i64],
) -> Result<malachite::Integer, String> {
    let target_i32 = target_class
        .iter()
        .map(|&value| {
            i32::try_from(value)
                .map_err(|_| format!("{label} GV target coordinate does not fit in i32"))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        cyrus_core::compute_gv_invariants_with_provided_generators(
            provided_generators,
            grading,
            q_matrix,
            intnums,
            None,
            Some(max_deg),
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    match gvs_result {
        Ok(Ok(gvs)) => Ok(gvs
            .into_iter()
            .find_map(|(curve, gv)| (curve == target_i32).then_some(gv))
            .unwrap_or_else(|| malachite::Integer::from(0))),
        Ok(Err(e)) => Err(format!("failed {label} GV computation: {e}")),
        Err(payload) => Err(format!(
            "{label} GV computation panicked for ambient_nonzero={:?}: {}",
            sparse_i64(ambient_class),
            panic_payload_message(payload.as_ref())
        )),
    }
}

fn compute_lattice_saturated_face_target_gv(
    seed_generators: &[Vec<i64>],
    target_class: &[i64],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    max_deg: u32,
    ambient_class: &[i64],
) -> Result<(malachite::Integer, usize), String> {
    let semigroup_elements = degree_bounded_face_lattice_points(seed_generators, grading, max_deg)?;
    if !semigroup_elements
        .iter()
        .any(|element| element.as_slice() == target_class)
    {
        return Err(format!(
            "lattice-saturated face enumeration did not include target ambient_nonzero={:?}",
            sparse_i64(ambient_class)
        ));
    }

    let element_limit = std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_LATTICE_ELEMENT_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_LATTICE_ELEMENT_LIMIT);
    if semigroup_elements.len() > element_limit {
        return Err(format!(
            "lattice-saturated face enumeration produced {} semigroup elements, exceeding limit {element_limit}",
            semigroup_elements.len()
        ));
    }

    let gv = compute_explicit_semigroup_target_gv(
        &semigroup_elements,
        target_class,
        grading,
        q_matrix,
        intnums,
        "lattice-saturated face",
        ambient_class,
    )?;
    Ok((gv, semigroup_elements.len()))
}

#[allow(clippy::too_many_arguments)]
fn compute_explicit_semigroup_target_gv(
    semigroup_elements: &[Vec<i64>],
    target_class: &[i64],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    label: &str,
    ambient_class: &[i64],
) -> Result<malachite::Integer, String> {
    let target_i32 = target_class
        .iter()
        .map(|&value| {
            i32::try_from(value)
                .map_err(|_| format!("{label} GV target coordinate does not fit in i32"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let gvs = cyrus_core::compute_gv_invariants_with_explicit_semigroup(
        semigroup_elements,
        grading,
        q_matrix,
        intnums,
    )
    .map_err(|e| {
        format!(
            "failed {label} GV computation with {} explicit semigroup elements for ambient_nonzero={:?}: {e}",
            semigroup_elements.len(),
            sparse_i64(ambient_class),
        )
    })?;
    let gv = gvs
        .into_iter()
        .find_map(|(curve, gv)| (curve == target_i32).then_some(gv))
        .unwrap_or_else(|| malachite::Integer::from(0));
    Ok(gv)
}

#[allow(clippy::too_many_arguments)]
fn compute_integer_decomposition_face_target_gv(
    basis_rays: &[Vec<i64>],
    target_class: &[i64],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    max_deg: u32,
    ambient_class: &[i64],
    max_terms: usize,
    max_witnesses: usize,
) -> Result<(malachite::Integer, usize, usize, bool, bool), String> {
    let decompositions = exact_generator_decompositions(
        target_class,
        basis_rays,
        grading,
        i128::from(max_deg),
        max_terms,
        max_witnesses,
    )?
    .ok_or_else(|| {
        format!("no exact integer generator decomposition found with up to {max_terms} terms")
    })?;
    let mut errors = Vec::new();
    for decomposition in &decompositions {
        let diamond_elements = decomposition_diamond_elements(decomposition, target_class)?;
        match compute_explicit_semigroup_target_gv(
            &diamond_elements,
            target_class,
            grading,
            q_matrix,
            intnums,
            "decomposition-diamond face",
            ambient_class,
        ) {
            Ok(gv) => {
                return Ok((gv, diamond_elements.len(), decomposition.len(), true, false));
            }
            Err(error) => {
                let decomposition_sparse = decomposition
                    .iter()
                    .map(|ray| sparse_i64(ray))
                    .collect::<Vec<_>>();
                errors.push(format!(
                    "{}-term decomposition diamond {:?} with {} elements failed: {error}",
                    decomposition.len(),
                    decomposition_sparse,
                    diamond_elements.len()
                ));
            }
        }
    }

    let closure_limit =
        std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_DECOMPOSITION_CLOSURE_ELEMENT_LIMIT")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_DECOMPOSITION_CLOSURE_ELEMENT_LIMIT);
    match compute_decomposition_closure_face_target_gv(
        basis_rays,
        target_class,
        grading,
        q_matrix,
        intnums,
        ambient_class,
        max_terms,
        max_witnesses,
        closure_limit,
    ) {
        Ok((gv, element_count)) => return Ok((gv, element_count, max_terms, false, true)),
        Err(error) => errors.push(format!(
            "recursive decomposition closure with element limit {closure_limit} failed: {error}"
        )),
    }

    for decomposition in &decompositions {
        let mut seed_generators = decomposition.clone();
        seed_generators.push(target_class.to_vec());
        seed_generators.sort();
        seed_generators.dedup();
        match compute_lattice_saturated_face_target_gv(
            &seed_generators,
            target_class,
            grading,
            q_matrix,
            intnums,
            max_deg,
            ambient_class,
        ) {
            Ok((gv, element_count)) => {
                return Ok((gv, element_count, decomposition.len(), false, false));
            }
            Err(error) => {
                let decomposition_sparse = decomposition
                    .iter()
                    .map(|ray| sparse_i64(ray))
                    .collect::<Vec<_>>();
                errors.push(format!(
                    "{}-term decomposition {:?} with {} seed generators failed: {error}",
                    decomposition.len(),
                    decomposition_sparse,
                    seed_generators.len()
                ));
            }
        }
    }
    Err(format!(
        "{} exact integer decomposition witnesses failed: {}",
        decompositions.len(),
        errors.join(" | ")
    ))
}

#[allow(clippy::too_many_arguments)]
fn compute_decomposition_closure_face_target_gv(
    basis_rays: &[Vec<i64>],
    target_class: &[i64],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    ambient_class: &[i64],
    max_terms: usize,
    max_witnesses: usize,
    element_limit: usize,
) -> Result<(malachite::Integer, usize), String> {
    let semigroup_elements = exact_decomposition_closure_elements(
        target_class,
        basis_rays,
        grading,
        max_terms,
        max_witnesses,
        element_limit,
    )?;
    let gv = compute_explicit_semigroup_target_gv(
        &semigroup_elements,
        target_class,
        grading,
        q_matrix,
        intnums,
        "decomposition-closure face",
        ambient_class,
    )?;
    Ok((gv, semigroup_elements.len()))
}

fn exact_decomposition_closure_elements(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    max_terms: usize,
    max_witnesses: usize,
    element_limit: usize,
) -> Result<Vec<Vec<i64>>, String> {
    if element_limit == 0 {
        return Err("decomposition closure element limit is zero".to_string());
    }
    if target.len() != grading.len() {
        return Err(format!(
            "target dimension {} does not match grading dimension {}",
            target.len(),
            grading.len()
        ));
    }

    let zero = vec![0i64; target.len()];
    let mut elements = vec![zero.clone(), target.to_vec()];
    let mut seen = HashSet::from([zero, target.to_vec()]);
    let mut processed = HashSet::new();
    let mut queue = VecDeque::from([target.to_vec()]);
    while let Some(element) = queue.pop_front() {
        if !processed.insert(element.clone()) {
            continue;
        }
        if element.iter().all(|&value| value == 0) {
            continue;
        }
        let degree = element
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 {
            return Err(format!(
                "decomposition closure reached non-positive degree {degree} for element {:?}",
                sparse_i64(&element)
            ));
        }
        let Some(decompositions) = exact_generator_decompositions(
            &element,
            basis_rays,
            grading,
            degree,
            max_terms,
            max_witnesses,
        )?
        else {
            continue;
        };
        for decomposition in decompositions {
            for diamond_element in decomposition_diamond_elements(&decomposition, &element)? {
                if seen.insert(diamond_element.clone()) {
                    if seen.len() > element_limit {
                        return Err(format!(
                            "decomposition closure exceeded element limit {element_limit}"
                        ));
                    }
                    if diamond_element.iter().any(|&value| value != 0) {
                        queue.push_back(diamond_element.clone());
                    }
                    elements.push(diamond_element);
                }
            }
        }
    }
    elements.sort();
    Ok(elements)
}

fn decomposition_diamond_elements(
    decomposition: &[Vec<i64>],
    target: &[i64],
) -> Result<Vec<Vec<i64>>, String> {
    let dim = target.len();
    let zero = vec![0i64; dim];
    let mut elements = vec![zero.clone()];
    let mut seen = HashSet::from([zero]);
    for term in decomposition {
        if term.len() != dim {
            return Err(format!(
                "decomposition term dimension {} does not match target dimension {dim}",
                term.len()
            ));
        }
        let existing = elements.clone();
        for element in existing {
            let mut sum = Vec::with_capacity(dim);
            for (&lhs, &rhs) in element.iter().zip(term.iter()) {
                sum.push(lhs.checked_add(rhs).ok_or_else(|| {
                    "decomposition diamond element coordinate overflowed i64".to_string()
                })?);
            }
            if seen.insert(sum.clone()) {
                elements.push(sum);
            }
        }
    }
    if !seen.contains(target) {
        return Err(format!(
            "decomposition diamond does not contain target {:?}",
            sparse_i64(target)
        ));
    }
    elements.sort();
    Ok(elements)
}

fn exact_generator_decompositions(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    target_degree: i128,
    max_terms: usize,
    max_witnesses: usize,
) -> Result<Option<Vec<Vec<Vec<i64>>>>, String> {
    if max_terms < 2 {
        return Ok(None);
    }
    if target.len() != grading.len() {
        return Err(format!(
            "target dimension {} does not match grading vector length {}",
            target.len(),
            grading.len()
        ));
    }

    let mut candidates = Vec::new();
    let mut candidate_set = HashSet::new();
    for ray in basis_rays {
        if ray.len() != target.len() {
            return Err(format!(
                "basis ray dimension {} does not match target dimension {}",
                ray.len(),
                target.len()
            ));
        }
        if ray.as_slice() == target {
            continue;
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree > 0 && degree < target_degree && candidate_set.insert(ray.clone()) {
            candidates.push(ray.clone());
        }
    }
    candidates.sort();
    if candidates.is_empty() {
        return Ok(None);
    }

    let mut decompositions = Vec::new();
    let mut seen_decompositions = HashSet::new();
    let mut remainder = vec![0i64; target.len()];
    for first in &candidates {
        subtract_two_into(target, first, &mut remainder)?;
        if let Some(second) = candidate_set.get(remainder.as_slice()) {
            let mut decomposition = vec![first.clone(), second.clone()];
            decomposition.sort();
            if seen_decompositions.insert(decomposition.clone()) {
                decompositions.push(decomposition);
                if decompositions.len() >= max_witnesses {
                    return Ok(Some(decompositions));
                }
            }
        }
    }
    if max_terms < 3 || decompositions.len() >= max_witnesses {
        return Ok((!decompositions.is_empty()).then_some(decompositions));
    }

    for (i, first) in candidates.iter().enumerate() {
        for second in candidates.iter().skip(i) {
            subtract_three_into(target, first, second, &mut remainder)?;
            if let Some(third) = candidate_set.get(remainder.as_slice()) {
                let mut decomposition = vec![first.clone(), second.clone(), third.clone()];
                decomposition.sort();
                if seen_decompositions.insert(decomposition.clone()) {
                    decompositions.push(decomposition);
                    if decompositions.len() >= max_witnesses {
                        return Ok(Some(decompositions));
                    }
                }
            }
        }
    }
    Ok((!decompositions.is_empty()).then_some(decompositions))
}

fn subtract_two_into(target: &[i64], first: &[i64], out: &mut [i64]) -> Result<(), String> {
    if target.len() != first.len() || target.len() != out.len() {
        return Err("subtract_two_into dimensions are inconsistent".to_string());
    }
    for ((slot, &target_value), &first_value) in out.iter_mut().zip(target).zip(first) {
        *slot = target_value - first_value;
    }
    Ok(())
}

fn subtract_three_into(
    target: &[i64],
    first: &[i64],
    second: &[i64],
    out: &mut [i64],
) -> Result<(), String> {
    if target.len() != first.len() || target.len() != second.len() || target.len() != out.len() {
        return Err("subtract_three_into dimensions are inconsistent".to_string());
    }
    for (((slot, &target_value), &first_value), &second_value) in
        out.iter_mut().zip(target).zip(first).zip(second)
    {
        *slot = target_value - first_value - second_value;
    }
    Ok(())
}

fn degree_bounded_face_lattice_points(
    seed_generators: &[Vec<i64>],
    grading: &[i64],
    max_deg: u32,
) -> Result<Vec<Vec<i64>>, String> {
    if seed_generators.is_empty() {
        return Err("face lattice saturation requires at least one seed generator".to_string());
    }
    let dim = seed_generators[0].len();
    if dim != grading.len() {
        return Err(format!(
            "seed generator dimension {dim} does not match grading vector length {}",
            grading.len()
        ));
    }
    if seed_generators.iter().any(|row| row.len() != dim) {
        return Err("seed generator dimensions are inconsistent".to_string());
    }

    let support: Vec<usize> = (0..dim)
        .filter(|&coord| seed_generators.iter().any(|row| row[coord] != 0))
        .collect();
    if support.is_empty() {
        return Err("face lattice saturation seed generators are all zero".to_string());
    }

    let reduced_rays = seed_generators
        .iter()
        .map(|row| {
            support
                .iter()
                .map(|&coord| i128::from(row[coord]))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let reduced_grading = support
        .iter()
        .map(|&coord| grading[coord])
        .collect::<Vec<_>>();
    let mut cone = cyrus_core::Cone::from_rays(reduced_rays);
    let reduced_points = cone
        .find_lattice_points_ortools(None, Some(i64::from(max_deg)), &reduced_grading, 1000, 0)
        .map_err(|e| format!("failed to enumerate lattice-saturated face points: {e}"))?;

    let mut lifted_points = Vec::with_capacity(reduced_points.len());
    for point in reduced_points {
        if point.len() != support.len() {
            return Err(format!(
                "face lattice point dimension {} does not match support dimension {}",
                point.len(),
                support.len()
            ));
        }
        let mut lifted = vec![0i64; dim];
        for (&coord, value) in support.iter().zip(point) {
            lifted[coord] = value;
        }
        lifted_points.push(lifted);
    }
    lifted_points.sort();
    lifted_points.dedup();
    Ok(lifted_points)
}

fn degree_bounded_span_generators(
    seed_rows: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    grading: &[i64],
    target_degree: i128,
) -> Result<Vec<Vec<i64>>, String> {
    if seed_rows.is_empty() {
        return Ok(Vec::new());
    }
    let dim = seed_rows[0].len();
    if dim != grading.len() {
        return Err(format!(
            "seed row dimension {dim} does not match grading vector length {}",
            grading.len()
        ));
    }
    if seed_rows.iter().any(|row| row.len() != dim) {
        return Err("seed row dimensions are inconsistent".to_string());
    }
    let seed_rank = rational_rank_i64_rows(seed_rows);
    let mut out = Vec::new();
    for ray in basis_rays {
        if ray.len() != dim {
            return Err(format!(
                "basis ray dimension {} does not match seed dimension {dim}",
                ray.len()
            ));
        }
        let degree = ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= 0 || degree > target_degree {
            continue;
        }
        let mut rows = seed_rows.to_vec();
        rows.push(ray.clone());
        if rational_rank_i64_rows(&rows) == seed_rank {
            out.push(ray.clone());
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
}

fn rational_rank_i64_rows(rows: &[Vec<i64>]) -> usize {
    let rational_rows = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|&value| malachite::Rational::from(value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    cyrus_core::integer_math::matrix_rank(&rational_rows)
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

fn sparse_i64(values: &[i64]) -> Vec<(usize, i64)> {
    values
        .iter()
        .enumerate()
        .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
        .collect()
}

fn insert_missing_diagnostic_gv(
    diagnostic_gvs: &mut HashMap<Vec<i64>, malachite::Integer>,
    ambient_class: &[i64],
    gv: &malachite::Integer,
    source: &str,
) -> Result<(), String> {
    let sparse_class = sparse_i64(ambient_class);
    match diagnostic_gvs.entry(ambient_class.to_vec()) {
        std::collections::hash_map::Entry::Occupied(existing) => {
            if existing.get() != gv {
                return Err(format!(
                    "conflicting corrected-chamber diagnostic GV values for {source} class {sparse_class:?}: existing={} new={gv}",
                    existing.get()
                ));
            }
        }
        std::collections::hash_map::Entry::Vacant(slot) => {
            slot.insert(gv.clone());
        }
    }
    Ok(())
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
            .map(|value| value.get()),
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
                .map(|value| value.get()),
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
    scale_mixed_basis_kklt_branch_initialization_to_target(
        &intersection.kappa_basis,
        &intersection.kappa_full,
        &intersection.basis,
        kklt_basis,
        tau_phase1,
        &raw,
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

fn chamber_intersection_in_basis(
    tri: &Triangulation,
    points: &[Point],
    basis: &[usize],
) -> Result<cyrus_core::Intersection, String> {
    let kappa_full = chamber_intersection_full(tri, points)?;
    Ok(intersection_in_basis(&kappa_full, basis))
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

fn diagnose_chamber_gv_volume_correction(
    tri: &Triangulation,
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    kklt_basis: &[usize],
    kahler: &[F64<Finite>],
    gamma: &[I64<Finite>],
    cutoff: F64<Pos>,
    general_min_points: Option<u32>,
    general_max_deg: Option<u32>,
    provided_generators_only: bool,
    ray_gv_requested: bool,
    lp_face_gv_requested: bool,
) -> Result<ChamberGvDiagnostic, String> {
    let ambient_rays = compute_mori_cone_cap_rays(
        tri,
        &geom.triangulation_points,
        &geom.polytope,
        false,
        false,
        None,
    )
    .map_err(|e| format!("failed to compute corrected-chamber ambient Mori-cap rays: {e}"))?;
    let small_curve_candidates =
        subcutoff_toric_curve_candidates(&ambient_rays, &intersection.basis, kahler, cutoff)
            .map_err(|e| format!("failed to select corrected-chamber small curves: {e}"))?;
    let small_curves = remove_pair_decomposable_curve_candidates(&small_curve_candidates)
        .map_err(|e| format!("failed to prune corrected-chamber small curves: {e}"))?;
    let toric_gvs =
        compute_toric_two_face_curve_gv_invariants(tri, &geom.triangulation_points, &geom.polytope)
            .map_err(|e| format!("failed to compute corrected-chamber toric GV values: {e}"))?;

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
    let toric_gv_covered_count = small_curve_gvs.len();
    let toric_gv_missing_count = missing_gv_classes.len();
    let toric_small_curve_gvs = small_curve_gvs.clone();
    let toric_missing_gv_classes = missing_gv_classes.clone();

    let mut basis_ray_stats = None;
    let mut basis_rays_for_missing = None;
    let mut grading_for_missing = None;
    let mut degree_summary = None;
    let mut missing_target_stats = None;
    let mut covered_toric_gv_divisor_representation_baseline = None;
    if !missing_gv_classes.is_empty() {
        let origin_idx = geom
            .triangulation_points
            .iter()
            .position(|point| point.coords().iter().all(|&coord| coord == 0))
            .ok_or_else(|| "failed to find origin in corrected-chamber points".to_string())?;
        let basis_rays = project_mori_cone_cap_rays_to_basis(&ambient_rays, &intersection.basis)
            .map_err(|e| {
                format!("failed to project corrected-chamber Mori-cap rays to basis: {e}")
            })?;
        let origin_circuits = compute_origin_circuit_curve_diagnostics(
            tri,
            &geom.triangulation_points,
            &geom.polytope,
        )
        .map_err(|e| {
            format!("failed to compute corrected-chamber origin-circuit diagnostics: {e}")
        })?;
        let origin_circuits_by_class: HashMap<Vec<i64>, cyrus_core::OriginCircuitCurveDiagnostic> =
            origin_circuits
                .into_iter()
                .map(|diagnostic| (diagnostic.class.clone(), diagnostic))
                .collect();
        let corrected_kappa_full = chamber_intersection_full(tri, &geom.triangulation_points)?;
        let cms_intersection_checks_by_class = cms_general_divisor_intersection_checks_by_class(
            &missing_gv_classes,
            &origin_circuits_by_class,
            &corrected_kappa_full,
            &intersection.basis,
        )?;
        covered_toric_gv_divisor_representation_baseline =
            compute_covered_toric_gv_divisor_representation_baseline(
                &small_curve_gvs,
                origin_idx,
                &corrected_kappa_full,
                &intersection.basis,
                covered_toric_gv_divisor_representation_class_limit(),
            )?;
        let grading = compute_grading_vector(&basis_rays).ok_or_else(|| {
            "failed to compute corrected-chamber GV degree grading vector".to_string()
        })?;
        let summary = summarize_required_gv_degrees(
            &missing_gv_classes,
            &intersection.basis,
            &grading,
            general_max_deg,
        )?;
        let ray_stats = graded_ray_stats(&basis_rays, &grading, general_max_deg)?;
        let target_stats = missing_gv_target_stats(
            &missing_gv_classes,
            &basis_rays,
            &intersection.basis,
            &grading,
            origin_idx,
            &origin_circuits_by_class,
            &cms_intersection_checks_by_class,
            10,
        )?;
        basis_rays_for_missing = Some(basis_rays);
        grading_for_missing = Some(grading);
        degree_summary = Some(summary);
        basis_ray_stats = Some(ray_stats);
        missing_target_stats = Some(target_stats);
    }

    let general_gv_requested = general_min_points.is_some() || general_max_deg.is_some();
    let mut general_gv_covered_count = None;
    let mut ray_gv_covered_count = None;
    let mut ray_gv_sample = Vec::new();
    let mut ray_gv_volume_correction = None;
    let mut lp_face_gv_covered_count = None;
    let mut lp_face_gv_failed_count = None;
    let mut lp_face_gv_certified_count = None;
    let mut lp_face_gv_uncertified_count = None;
    let mut lp_face_gv_sample = Vec::new();
    let mut lp_face_gv_volume_correction = None;
    let mut diagnostic_missing_gvs: HashMap<Vec<i64>, malachite::Integer> = HashMap::new();
    if !missing_gv_classes.is_empty() && ray_gv_requested {
        if cfg!(panic = "abort") {
            return Err(
                "one-dimensional ray GV diagnostic is disabled in panic=abort builds because cygv reports inconsistent ray truncations by panicking; run a panic=unwind/debug build or port cygv errors to Result before using this diagnostic"
                    .to_string(),
            );
        }
        let basis_rays = basis_rays_for_missing
            .as_ref()
            .expect("basis rays computed for corrected-chamber missing curves");
        let grading = grading_for_missing
            .as_ref()
            .expect("grading computed for corrected-chamber missing curves");
        let missing_target_stats = missing_target_stats
            .as_ref()
            .expect("missing target stats computed for corrected-chamber missing curves");
        let (non_positive_count, first_non_positive) =
            non_positive_basis_generator_degrees(basis_rays, grading)?;
        if let Some((idx, degree, ray)) = first_non_positive {
            return Err(format!(
                "one-dimensional ray GV diagnostic requires a grading positive on all Mori generators; found {non_positive_count}/{} non-positive generator degrees, first index={idx} degree={degree} ray={ray:?}",
                basis_rays.len()
            ));
        }
        let ray_targets = one_dimensional_ray_gv_targets(
            &missing_gv_classes,
            basis_rays,
            &intersection.basis,
            grading,
        )?;
        if ray_targets.candidates.is_empty() {
            return Err(format!(
                "one-dimensional ray GV diagnostic found no LP-extremal primitive Mori-generator targets among {} missing corrected-chamber curves; skipped_non_generators={} skipped_decomposable_generators={}",
                missing_target_stats.target_count,
                ray_targets.skipped_non_generators,
                ray_targets.skipped_decomposable_generators
            ));
        }
        if ray_targets.candidates.len() != missing_target_stats.target_count {
            eprintln!(
                "[WARN] corrected-chamber one-dimensional ray GV diagnostic will try {}/{} missing curves; skipped_non_generators={} skipped_decomposable_generators={}",
                ray_targets.candidates.len(),
                missing_target_stats.target_count,
                ray_targets.skipped_non_generators,
                ray_targets.skipped_decomposable_generators
            );
        }
        let curve_basis = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
            .map_err(|e| format!("failed to compute corrected-chamber curve basis matrix: {e}"))?;
        let q_matrix = curve_basis
            .iter()
            .map(|row| {
                row.iter()
                    .skip(1)
                    .map(|value| {
                        i64::try_from(value).map_err(|_| {
                            "corrected-chamber q-matrix entry does not fit in i64".to_string()
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let corrected_kappa_basis =
            chamber_intersection_in_basis(tri, &geom.triangulation_points, &intersection.basis)?;
        eprintln!(
            "[WARN] corrected-chamber one-dimensional ray GV diagnostic assumes each LP-extremal primitive target spans a valid Mori-cone face; this is not yet promoted to the exact corrected-chamber GV fallback."
        );
        let ray_gvs = compute_missing_one_dimensional_ray_gvs(
            &ray_targets.candidates,
            &intersection.basis,
            grading,
            &q_matrix,
            &corrected_kappa_basis,
        )?;
        for (ambient_class, gv, _) in &ray_gvs {
            insert_missing_diagnostic_gv(
                &mut diagnostic_missing_gvs,
                ambient_class,
                gv,
                "one-dimensional ray",
            )?;
        }
        ray_gv_covered_count = Some(ray_gvs.len());
        ray_gv_sample = ray_gvs
            .iter()
            .take(10)
            .map(|(ambient_class, gv, degree)| RayGvDiagnosticSample {
                degree: *degree,
                gv: gv.clone(),
                ambient_nonzero: ambient_class
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
                    .collect(),
            })
            .collect();
        let mut ray_small_curve_gvs = small_curve_gvs.clone();
        ray_small_curve_gvs.extend(
            ray_gvs
                .iter()
                .map(|(ambient_class, gv, _)| (ambient_class.clone(), gv.clone())),
        );
        ray_gv_volume_correction = Some(
            cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                &ray_small_curve_gvs,
                &intersection.basis,
                kahler,
                Some(gamma),
            )
            .ok_or_else(|| {
                "failed to compute one-dimensional ray GV volume correction".to_string()
            })?,
        );
    }
    if !missing_gv_classes.is_empty() && lp_face_gv_requested {
        if cfg!(panic = "abort") {
            return Err(
                "LP-witness face GV diagnostic is disabled in panic=abort builds because cygv reports inconsistent face truncations by panicking; run with RUSTFLAGS='-C panic=unwind' or port cygv errors to Result before using this diagnostic"
                    .to_string(),
            );
        }
        let basis_rays = basis_rays_for_missing
            .as_ref()
            .expect("basis rays computed for corrected-chamber missing curves");
        let grading = grading_for_missing
            .as_ref()
            .expect("grading computed for corrected-chamber missing curves");
        let (non_positive_count, first_non_positive) =
            non_positive_basis_generator_degrees(basis_rays, grading)?;
        if let Some((idx, degree, ray)) = first_non_positive {
            return Err(format!(
                "LP-witness face GV diagnostic requires a grading positive on all Mori generators; found {non_positive_count}/{} non-positive generator degrees, first index={idx} degree={degree} ray={ray:?}",
                basis_rays.len()
            ));
        }
        let curve_basis = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
            .map_err(|e| format!("failed to compute corrected-chamber curve basis matrix: {e}"))?;
        let q_matrix = curve_basis
            .iter()
            .map(|row| {
                row.iter()
                    .skip(1)
                    .map(|value| {
                        i64::try_from(value).map_err(|_| {
                            "corrected-chamber q-matrix entry does not fit in i64".to_string()
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let corrected_kappa_basis =
            chamber_intersection_in_basis(tri, &geom.triangulation_points, &intersection.basis)?;
        eprintln!(
            "[WARN] corrected-chamber LP-witness face GV diagnostic uses small provided-generator cones from floating LP witnesses, with decomposition-diamond, recursive decomposition-closure, and lattice-saturated local retries on failures; this is not yet promoted to the exact corrected-chamber GV fallback."
        );
        let face_gvs = compute_missing_lp_witness_face_gvs(
            &missing_gv_classes,
            basis_rays,
            &intersection.basis,
            grading,
            &q_matrix,
            &corrected_kappa_basis,
        )?;
        let covered_count = face_gvs.iter().filter(|result| result.gv.is_some()).count();
        let failed_count = face_gvs.len().saturating_sub(covered_count);
        if supporting_face_certificate_requested() {
            let certified_count = face_gvs
                .iter()
                .filter(|result| {
                    result.gv.is_some() && result.supporting_face_certificate.is_some()
                })
                .count();
            lp_face_gv_certified_count = Some(certified_count);
            lp_face_gv_uncertified_count = Some(covered_count.saturating_sub(certified_count));
        }
        for result in &face_gvs {
            if let Some(gv) = result.gv.as_ref() {
                insert_missing_diagnostic_gv(
                    &mut diagnostic_missing_gvs,
                    &result.ambient_class,
                    gv,
                    "LP-witness face",
                )?;
            }
        }
        lp_face_gv_covered_count = Some(covered_count);
        lp_face_gv_failed_count = Some(failed_count);
        lp_face_gv_sample = face_gvs
            .iter()
            .take(10)
            .map(|result| FaceGvDiagnosticSample {
                degree: result.degree,
                generator_count: result.generator_count,
                active_generator_count: result.active_generator_count,
                span_generator_count: result.span_generator_count,
                used_span_expansion: result.used_span_expansion,
                used_lattice_saturation: result.used_lattice_saturation,
                used_integer_decomposition: result.used_integer_decomposition,
                used_decomposition_diamond: result.used_decomposition_diamond,
                used_decomposition_closure: result.used_decomposition_closure,
                integer_decomposition_term_count: result.integer_decomposition_term_count,
                lattice_semigroup_element_count: result.lattice_semigroup_element_count,
                supporting_face_certificate: result.supporting_face_certificate.as_ref().map(
                    |certificate| SupportingFaceCertificateSummary {
                        zero_generator_count: certificate.zero_generator_count,
                        positive_generator_count: certificate.positive_generator_count,
                        normal_nonzero: sparse_i64(&certificate.normal),
                    },
                ),
                gv: result.gv.clone(),
                error: result.error.clone(),
                ambient_nonzero: result
                    .ambient_class
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
                    .collect(),
            })
            .collect();
        if covered_count > 0 {
            let mut face_small_curve_gvs = small_curve_gvs.clone();
            face_small_curve_gvs.extend(face_gvs.iter().filter_map(|result| {
                result
                    .gv
                    .as_ref()
                    .map(|gv| (result.ambient_class.clone(), gv.clone()))
            }));
            lp_face_gv_volume_correction = Some(
                cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                    &face_small_curve_gvs,
                    &intersection.basis,
                    kahler,
                    Some(gamma),
                )
                .ok_or_else(|| {
                    "failed to compute LP-witness face GV volume correction".to_string()
                })?,
            );
        }
    }
    if !missing_gv_classes.is_empty() && general_gv_requested {
        let summary = degree_summary
            .as_ref()
            .expect("degree summary computed for missing corrected-chamber curves");
        if let Some(max_deg) = general_max_deg
            && summary.max_degree > i128::from(max_deg)
        {
            return Err(format!(
                "corrected-chamber general GV max_deg={max_deg} cannot cover all missing curves: required degree range {}..{} ({} curves), first_over_max={:?}",
                summary.min_degree, summary.max_degree, summary.count, summary.first_over_max
            ));
        }
        let basis_rays = basis_rays_for_missing
            .as_ref()
            .expect("basis rays computed for corrected-chamber missing curves");
        let grading = grading_for_missing
            .as_ref()
            .expect("grading computed for corrected-chamber missing curves");
        let direct_ray_limit = std::env::var("CYRUS_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(DEFAULT_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT);
        if basis_rays.len() > direct_ray_limit && !provided_generators_only {
            let bounded_note = basis_ray_stats
                .as_ref()
                .and_then(|stats| stats.degree_bounded_count)
                .map_or_else(
                    || "no max-degree ray count available".to_string(),
                    |count| format!("{count} rays have degree <= requested max_deg"),
                );
            return Err(format!(
                "corrected-chamber general GV direct fallback would dualize {} Mori generators, exceeding limit {direct_ray_limit}; {bounded_note}. Current backend computes cone hyperplanes by DDM before bounded lattice enumeration. Filtering to degree-bounded rays before dualization is not exact in general because higher-degree rays can still affect the low-degree lattice slice. Need a reduced corrected-chamber cone/lattice formulation, or set CYRUS_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT to force the direct attempt.",
                basis_rays.len(),
            ));
        }
        if basis_rays.len() > direct_ray_limit
            && provided_generators_only
            && general_max_deg.is_none()
        {
            return Err(format!(
                "corrected-chamber provided-generator GV diagnostic would pass {} Mori generators without a max-degree bound; supply --primal-gv-max-deg or increase CYRUS_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT for an explicit direct attempt",
                basis_rays.len()
            ));
        }
        let (non_positive_count, first_non_positive) =
            non_positive_basis_generator_degrees(basis_rays, grading)?;
        if let Some((idx, degree, ray)) = first_non_positive {
            return Err(format!(
                "corrected-chamber general GV fallback requires a grading positive on all Mori generators; found {non_positive_count}/{} non-positive generator degrees, first index={idx} degree={degree} ray={ray:?}",
                basis_rays.len()
            ));
        }
        eprintln!(
            "[INFO] toric formulas missed {} corrected-chamber small curves; computing corrected-chamber general GV fallback with min_points={:?} max_deg={:?}",
            missing_gv_classes.len(),
            general_min_points,
            general_max_deg
        );
        let curve_basis = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
            .map_err(|e| format!("failed to compute corrected-chamber curve basis matrix: {e}"))?;
        let q_matrix = curve_basis
            .iter()
            .map(|row| {
                row.iter()
                    .skip(1)
                    .map(|value| {
                        i64::try_from(value).map_err(|_| {
                            "corrected-chamber q-matrix entry does not fit in i64".to_string()
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let corrected_kappa_basis =
            chamber_intersection_in_basis(tri, &geom.triangulation_points, &intersection.basis)?;
        let general_gvs = if provided_generators_only {
            let provided_rays = if let Some(max_deg) = general_max_deg {
                degree_filtered_basis_rays(basis_rays, grading, max_deg)?
            } else {
                basis_rays.clone()
            };
            let provided_generator_limit =
                std::env::var("CYRUS_CORRECTED_CHAMBER_PROVIDED_GV_GENERATOR_LIMIT")
                    .ok()
                    .and_then(|value| value.parse::<usize>().ok())
                    .unwrap_or(DEFAULT_CORRECTED_CHAMBER_PROVIDED_GV_GENERATOR_LIMIT);
            if provided_rays.len() > provided_generator_limit {
                return Err(format!(
                    "corrected-chamber provided-generator GV diagnostic would pass {} generators to cygv, exceeding limit {provided_generator_limit}; increase CYRUS_CORRECTED_CHAMBER_PROVIDED_GV_GENERATOR_LIMIT for an explicit long attempt, or use a lower --primal-gv-max-deg to inspect a smaller degree window.",
                    provided_rays.len()
                ));
            }
            eprintln!(
                "[WARN] corrected-chamber provided-generator GV diagnostic is using {} caller-provided generators without Mori-cone lattice augmentation; this is not the exact corrected-chamber GV fallback.",
                provided_rays.len()
            );
            cyrus_core::compute_gv_invariants_with_provided_generators(
                &provided_rays,
                grading,
                &q_matrix,
                &corrected_kappa_basis,
                general_min_points,
                general_max_deg,
            )
        } else {
            cyrus_core::compute_gv_invariants(
                basis_rays,
                grading,
                &q_matrix,
                &corrected_kappa_basis,
                general_min_points,
                general_max_deg,
            )
        }
        .map_err(|e| format!("failed to compute corrected-chamber general GV invariants: {e}"))?;
        let ambient_gvs =
            map_basis_gv_invariants_to_ambient(&general_gvs, &curve_basis).map_err(|e| {
                format!("failed to map corrected-chamber general GV invariants to ambient: {e}")
            })?;
        let missing_set: HashSet<Vec<i64>> = missing_gv_classes.iter().cloned().collect();
        let mut newly_covered = 0usize;
        for (class, gv) in ambient_gvs {
            match gv_by_class.entry(class) {
                std::collections::hash_map::Entry::Occupied(existing) => {
                    if existing.get() != &gv {
                        return Err(format!(
                            "toric/general GV conflict for corrected-chamber curve: {} vs {gv}",
                            existing.get()
                        ));
                    }
                }
                std::collections::hash_map::Entry::Vacant(slot) => {
                    if missing_set.contains(slot.key()) {
                        newly_covered += 1;
                    }
                    slot.insert(gv);
                }
            }
        }
        general_gv_covered_count = Some(newly_covered);

        small_curve_gvs.clear();
        missing_gv_classes.clear();
        for curve in &small_curves {
            match gv_by_class.get(&curve.class) {
                Some(gv) => small_curve_gvs.push((curve.class.clone(), gv.clone())),
                None => missing_gv_classes.push(curve.class.clone()),
            }
        }
    }

    let mut combined_diagnostic_gv_covered_count = None;
    let mut combined_diagnostic_gv_missing_count = None;
    let mut combined_diagnostic_gv_zero_count = None;
    let mut combined_diagnostic_gv_nonzero_count = None;
    let mut combined_diagnostic_gv_volume_correction = None;
    let mut combined_diagnostic_gv_target_correction = None;
    if !diagnostic_missing_gvs.is_empty() {
        if let Some(unexpected) = diagnostic_missing_gvs.keys().find(|class| {
            !toric_missing_gv_classes
                .iter()
                .any(|missing| missing == *class)
        }) {
            return Err(format!(
                "corrected-chamber diagnostic GV produced a class that was not a toric miss: {:?}",
                sparse_i64(unexpected)
            ));
        }
        let (zero_count, nonzero_count) = diagnostic_gv_value_counts(&diagnostic_missing_gvs);
        combined_diagnostic_gv_zero_count = Some(zero_count);
        combined_diagnostic_gv_nonzero_count = Some(nonzero_count);
        let mut combined_small_curve_gvs = toric_small_curve_gvs.clone();
        let mut covered_count = 0usize;
        for missing_class in &toric_missing_gv_classes {
            if let Some(gv) = diagnostic_missing_gvs.get(missing_class) {
                covered_count += 1;
                combined_small_curve_gvs.push((missing_class.clone(), gv.clone()));
            }
        }
        let missing_count = toric_missing_gv_classes.len().saturating_sub(covered_count);
        combined_diagnostic_gv_covered_count = Some(covered_count);
        combined_diagnostic_gv_missing_count = Some(missing_count);
        if covered_count > 0 {
            combined_diagnostic_gv_volume_correction = Some(
                cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                    &combined_small_curve_gvs,
                    &intersection.basis,
                    kahler,
                    Some(gamma),
                )
                .ok_or_else(|| {
                    "failed to compute combined diagnostic corrected-chamber GV volume correction"
                        .to_string()
                })?,
            );
            combined_diagnostic_gv_target_correction = Some(
                cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                    &combined_small_curve_gvs,
                    &intersection.basis,
                    kklt_basis,
                    kahler,
                    Some(gamma),
                )
                .ok_or_else(|| {
                    "failed to compute combined diagnostic corrected-chamber GV target correction"
                        .to_string()
                })?,
            );
        }
    }

    let covered_gv_target_correction = if small_curve_gvs.is_empty() {
        None
    } else {
        Some(
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &small_curve_gvs,
                &intersection.basis,
                kklt_basis,
                kahler,
                Some(gamma),
            )
            .ok_or_else(|| {
                "failed to compute corrected-chamber ambient GV target correction".to_string()
            })?,
        )
    };

    let covered_gv_volume_correction = if small_curve_gvs.is_empty() {
        None
    } else {
        Some(
            cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                &small_curve_gvs,
                &intersection.basis,
                kahler,
                Some(gamma),
            )
            .ok_or_else(|| {
                "failed to compute corrected-chamber ambient GV volume correction".to_string()
            })?,
        )
    };

    let gv_volume_correction = if missing_gv_classes.is_empty() {
        covered_gv_volume_correction
    } else {
        None
    };
    let (missing_required_degree_min, missing_required_degree_max) = degree_summary
        .map_or((None, None), |summary| {
            (Some(summary.min_degree), Some(summary.max_degree))
        });

    Ok(ChamberGvDiagnostic {
        ambient_rays: ambient_rays.len(),
        subcutoff_count: small_curve_candidates.len(),
        filtered_count: small_curves.len(),
        toric_gv_covered_count,
        toric_gv_missing_count,
        basis_mori_ray_count: basis_ray_stats.as_ref().map(|stats| stats.count),
        degree_bounded_basis_mori_ray_count: basis_ray_stats
            .as_ref()
            .and_then(|stats| stats.degree_bounded_count),
        basis_mori_ray_degree_min: basis_ray_stats.as_ref().and_then(|stats| stats.min_degree),
        basis_mori_ray_degree_max: basis_ray_stats.as_ref().and_then(|stats| stats.max_degree),
        general_gv_covered_count,
        ray_gv_covered_count,
        ray_gv_volume_correction,
        ray_gv_sample,
        lp_face_gv_covered_count,
        lp_face_gv_failed_count,
        lp_face_gv_certified_count,
        lp_face_gv_uncertified_count,
        lp_face_gv_volume_correction,
        lp_face_gv_sample,
        combined_diagnostic_gv_covered_count,
        combined_diagnostic_gv_missing_count,
        combined_diagnostic_gv_zero_count,
        combined_diagnostic_gv_nonzero_count,
        combined_diagnostic_gv_volume_correction,
        combined_diagnostic_gv_target_correction,
        remaining_gv_missing_count: missing_gv_classes.len(),
        first_missing_class: missing_gv_classes.first().cloned(),
        missing_required_degree_min,
        missing_required_degree_max,
        missing_target_stats,
        covered_toric_gv_divisor_representation_baseline,
        covered_gv_target_correction,
        covered_gv_volume_correction,
        gv_volume_correction,
    })
}

fn compute_chamber_toric_gv_selection(
    tri: &Triangulation,
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    kahler: &[F64<Finite>],
    cutoff: F64<Pos>,
) -> Result<ChamberToricGvSelection, String> {
    let ambient_rays = compute_mori_cone_cap_rays(
        tri,
        &geom.triangulation_points,
        &geom.polytope,
        false,
        false,
        None,
    )
    .map_err(|e| format!("failed to compute chamber ambient Mori-cap rays: {e}"))?;
    let small_curve_candidates =
        subcutoff_toric_curve_candidates(&ambient_rays, &intersection.basis, kahler, cutoff)
            .map_err(|e| format!("failed to select chamber small curves: {e}"))?;
    let small_curves = remove_pair_decomposable_curve_candidates(&small_curve_candidates)
        .map_err(|e| format!("failed to prune chamber small curves: {e}"))?;
    let toric_gvs =
        compute_toric_two_face_curve_gv_invariants(tri, &geom.triangulation_points, &geom.polytope)
            .map_err(|e| format!("failed to compute chamber toric GV values: {e}"))?;
    let gv_by_class: HashMap<Vec<i64>, malachite::Integer> = toric_gvs
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

    Ok(ChamberToricGvSelection {
        ambient_rays: ambient_rays.len(),
        subcutoff_count: small_curve_candidates.len(),
        filtered_count: small_curves.len(),
        toric_gv_covered_count: small_curve_gvs.len(),
        toric_gv_missing_count: missing_gv_classes.len(),
        first_missing_class: missing_gv_classes.first().cloned(),
        small_curve_gvs,
    })
}

#[allow(clippy::too_many_arguments)]
fn diagnose_chamber_updated_kklt_toric_only(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    c_tau: cyrus_core::types::physics::CTau,
    gamma: &[I64<Finite>],
    start_t: &[F64<Finite>],
    cutoff: F64<Pos>,
    kklt_steps: usize,
    max_iterations: usize,
) -> Result<ChamberUpdatedKkltDiagnostic, String> {
    if max_iterations == 0 {
        return Err("chamber-updated KKLT diagnostic iterations must be positive".to_string());
    }
    let mut current_t = start_t.to_vec();
    let mut final_classical_volume = None;
    let mut final_gv_volume_correction = None;
    let mut final_toric_missing_count = 0usize;
    let mut final_first_missing_class = None;
    let mut converged = false;
    let gv_tolerance = 1e-10f64;
    let mut completed_iterations = 0usize;

    for iter in 0..max_iterations {
        let tri = triangulation_from_kahler_point(geom, &intersection.basis, &current_t)?;
        let chamber_changed_from_input =
            !triangulations_have_same_simplices(&geom.triangulation, &tri);
        let kappa_full = chamber_intersection_full(&tri, &geom.triangulation_points)?;
        let kappa_basis = intersection_in_basis(&kappa_full, &intersection.basis);
        let chi_divisor = cyrus_core::compute_kklt_divisor_chi(
            &geom.polytope,
            &geom.triangulation_points,
            &kappa_full,
            kklt_basis,
        )
        .map_err(|e| format!("failed to compute chamber-updated divisor chi: {e}"))?;
        let selection =
            compute_chamber_toric_gv_selection(&tri, geom, intersection, &current_t, cutoff)?;
        if selection.small_curve_gvs.is_empty() {
            return Err(
                "chamber-updated KKLT diagnostic found no toric-covered small-curve GV values"
                    .to_string(),
            );
        }
        if selection.toric_gv_missing_count > 0 {
            let first_missing_sparse = selection
                .first_missing_class
                .as_deref()
                .map(sparse_i64)
                .unwrap_or_default();
            eprintln!(
                "[WARN] chamber-updated KKLT diagnostic iteration {iter} uses toric-covered GV values only; missing={} first_missing_sparse={:?}",
                selection.toric_gv_missing_count, first_missing_sparse
            );
        }
        let Some(gv_correction) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &selection.small_curve_gvs,
            &intersection.basis,
            kklt_basis,
            &current_t,
            Some(gamma),
        ) else {
            return Err(format!(
                "failed to compute chamber-updated GV target correction at iteration {iter}"
            ));
        };
        let Some(tau_target) = cyrus_core::kklt::compute_gv_corrected_target_tau(
            c_i,
            &chi_divisor,
            c_tau,
            &gv_correction,
        ) else {
            return Err(format!(
                "failed to build chamber-updated GV target tau at iteration {iter}"
            ));
        };
        let Some(next) = solve_mixed_basis_path_following(
            &kappa_basis,
            &kappa_full,
            &intersection.basis,
            kklt_basis,
            &tau_target,
            &current_t,
            CheckedRange::new(0, kklt_steps),
        ) else {
            return Err(format!(
                "chamber-updated mixed-basis KKLT solve failed at iteration {iter}"
            ));
        };
        if !next.converged {
            return Err(format!(
                "chamber-updated mixed-basis KKLT solve did not converge at iteration {iter}: rel_err={}",
                next.relative_error.get()
            ));
        }
        let max_relative_step = next
            .t
            .iter()
            .zip(current_t.iter())
            .map(|(new, old)| (new.get() - old.get()).abs() / (old.get().abs() + 1e-12))
            .fold(0.0f64, f64::max);
        let Some(gv_volume_correction) =
            cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                &selection.small_curve_gvs,
                &intersection.basis,
                &next.t,
                Some(gamma),
            )
        else {
            return Err(format!(
                "failed to compute chamber-updated GV volume correction at iteration {iter}"
            ));
        };
        let classical_volume = classical_volume_from_t(&kappa_basis, &next.t);
        eprintln!(
            "[INFO] chamber-updated KKLT diagnostic iteration {iter}: changed_from_input={} ambient_rays={} subcutoff={} pair_pruned={} toric_covered={} toric_missing={} max_relative_step={} rel_err={} V_classical={} GV_volume_correction={}",
            chamber_changed_from_input,
            selection.ambient_rays,
            selection.subcutoff_count,
            selection.filtered_count,
            selection.toric_gv_covered_count,
            selection.toric_gv_missing_count,
            max_relative_step,
            next.relative_error.get(),
            classical_volume,
            gv_volume_correction.get()
        );

        final_classical_volume = Some(classical_volume);
        final_gv_volume_correction = Some(gv_volume_correction);
        final_toric_missing_count = selection.toric_gv_missing_count;
        final_first_missing_class = selection.first_missing_class;
        completed_iterations = iter + 1;
        current_t = next.t;
        if max_relative_step <= gv_tolerance {
            converged = true;
            break;
        }
    }

    let final_classical_volume = final_classical_volume.ok_or_else(|| {
        "chamber-updated KKLT diagnostic did not complete any iterations".to_string()
    })?;
    let final_gv_volume_correction = final_gv_volume_correction.ok_or_else(|| {
        "chamber-updated KKLT diagnostic did not compute a GV volume correction".to_string()
    })?;
    Ok(ChamberUpdatedKkltDiagnostic {
        iterations: completed_iterations,
        converged,
        final_t: current_t,
        final_classical_volume,
        final_gv_volume_correction,
        final_toric_missing_count,
        final_first_missing_class,
    })
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
    let Some(dir) = data_dir.map(PathBuf::from) else {
        return None;
    };
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
    alternate_gv_correction: Option<(&str, &[F64<Finite>])>,
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
    let alternate_gv_summary = alternate_gv_correction.map(|(label, correction)| {
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
            label,
            target_correction_delta_summary(&implied_gv_correction, correction).unwrap_or_else(
                |e| {
                    eprintln!(
                        "[ERROR] failed to compare alternate implied GV target correction ({label}): {e}"
                    );
                    std::process::exit(2);
                },
            ),
        )
    });
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
    if let Some((label, summary)) = alternate_gv_summary {
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

fn stage_volume(
    data_dir: Option<&str>,
    manifest_dir: &PathBuf,
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
    diagnose_corrected_chamber_gv: bool,
    diagnose_corrected_chamber_provided_generators_gv: bool,
    diagnose_corrected_chamber_ray_gv: bool,
    diagnose_corrected_chamber_lp_face_gv: bool,
    diagnose_chamber_updated_kklt: bool,
    diagnose_chamber_updated_kklt_iterations: usize,
    small_curve_cutoff: F64<Pos>,
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
    if branch_report_decomposition_depth > 3 {
        eprintln!("[ERROR] --branch-report-decomposition-depth currently supports values 0..3");
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
    if (diagnose_corrected_chamber_gv
        || diagnose_corrected_chamber_provided_generators_gv
        || diagnose_corrected_chamber_ray_gv
        || diagnose_corrected_chamber_lp_face_gv
        || diagnose_chamber_updated_kklt)
        && allow_downstream_kahler
    {
        eprintln!(
            "[ERROR] corrected-chamber/chamber-updated diagnostics are only valid for first-principles runs, not downstream Kähler replay"
        );
        std::process::exit(2);
    }
    if diagnose_chamber_updated_kklt && diagnose_chamber_updated_kklt_iterations == 0 {
        eprintln!(
            "[ERROR] --diagnose-chamber-updated-kklt-iterations must be positive when --diagnose-chamber-updated-kklt is set"
        );
        std::process::exit(2);
    }
    let (
        t,
        gv_volume_correction,
        gamma_for_chamber_gv,
        kklt_basis_for_chamber_gv,
        input_chamber_gv_target_correction,
    ) = if allow_downstream_kahler {
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
        (t, None, None, None, None)
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
            compute_b_field_gamma_for_o7_divisors(intersection.kappa_full.dim(), &kklt_basis, &c_i);
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
        let (zeroth_order, small_curve_selection_t) = if branch_candidates == 0 {
            let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
            let height_init = height_projected_branch_initialization(
                geom,
                intersection,
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
            let Some(phase1) = cyrus_core::kklt::solve_mixed_basis_path_following(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
                &kklt_basis,
                &tau_phase1,
                &height_init,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] phase-1 mixed-basis KKLT path-following failed");
                std::process::exit(2);
            };
            if !phase1.converged {
                eprintln!(
                    "[ERROR] phase-1 mixed-basis KKLT path-following did not converge: rel_err={}",
                    phase1.relative_error.get()
                );
                std::process::exit(2);
            }
            let Some(result) = cyrus_core::kklt::solve_mixed_basis_path_following(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
                &kklt_basis,
                &tau_target,
                &phase1.t,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] zeroth-order mixed-basis KKLT path-following failed");
                std::process::exit(2);
            };
            if !result.converged {
                eprintln!(
                    "[ERROR] zeroth-order mixed-basis KKLT path-following did not converge: rel_err={}",
                    result.relative_error.get()
                );
                std::process::exit(2);
            }
            let small_curve_selection_t = phase1.t.clone();
            (result, small_curve_selection_t)
        } else {
            let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
            let Some(mut t_initializations) = generate_scaled_kklt_branch_initializations(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
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
            let branch_search = solve_mixed_basis_path_following_branch_candidates(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
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
                        .map(|value| value.get())
                );
            }
            if positive_branches.is_empty() {
                eprintln!("[ERROR] KKLT branch search found no positive-volume phase-1 branch");
                std::process::exit(2);
            }
            let mut branch_gv_coverages: Option<Vec<BranchGvCoverage>> = None;
            if branch_selection.requires_gv_coverage() {
                let coverages = compute_branch_gv_coverages(
                    geom,
                    intersection,
                    &positive_branches,
                    small_curve_cutoff,
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
                    .map(|(rank, _)| rank)
                    .unwrap_or_else(|| {
                        eprintln!(
                            "[ERROR] KKLT branch search found no full-rank positive-volume phase-1 branch for min-condition selection"
                        );
                        std::process::exit(2);
                    }),
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
            let small_curve_selection_t = best_branch.result.t.clone();
            if let Some(path) = branch_report_path {
                let report_path = PathBuf::from(path);
                let ctx = BranchReportContext {
                    branch_seed,
                    branch_selection,
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
                            &positive_branches,
                            small_curve_cutoff,
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
                    .map(|value| value.get())
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
            let Some(result) = solve_mixed_basis_path_following(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
                &kklt_basis,
                &tau_target,
                &best_branch.result.t,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] corrected mixed-basis KKLT solve failed after branch search");
                std::process::exit(2);
            };
            if !result.converged {
                eprintln!(
                    "[ERROR] corrected mixed-basis KKLT solve after branch search did not converge: rel_err={}",
                    result.relative_error.get()
                );
                std::process::exit(2);
            }
            (result, small_curve_selection_t)
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
        let small_curves = remove_pair_decomposable_curve_candidates(&small_curve_candidates)
            .unwrap_or_else(|e| {
                eprintln!("[ERROR] failed to prune pair-decomposable small curves: {e}");
                std::process::exit(2);
            });
        eprintln!(
            "[INFO] primal small toric curves: ambient_rays={} subcutoff={} pair_pruned={} cutoff={}",
            ambient_rays.len(),
            small_curve_candidates.len(),
            small_curves.len(),
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
            let basis_rays =
                project_mori_cone_cap_rays_to_basis(&ambient_rays, &intersection.basis)
                    .unwrap_or_else(|e| {
                        eprintln!("[ERROR] failed to project primal Mori-cap rays to basis: {e}");
                        std::process::exit(2);
                    });
            let general_gvs = compute_primal_general_gv_by_ambient_class(
                geom,
                intersection,
                &missing_gv_classes,
                Some(&basis_rays),
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

        let mut corrected = zeroth_order;
        let mut correction_source_t = small_curve_selection_t;
        let max_gv_iterations = 20usize;
        let gv_tolerance = 1e-10f64;
        let mut gv_converged = false;
        for iter in 0..max_gv_iterations {
            let previous_t = correction_source_t.clone();
            let Some(gv_correction) =
                cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                    &small_curve_gvs,
                    &intersection.basis,
                    &kklt_basis,
                    &previous_t,
                    Some(&gamma),
                )
            else {
                eprintln!(
                    "[ERROR] failed to compute primal ambient GV target correction at iteration {iter}"
                );
                std::process::exit(2);
            };
            let Some(tau_target) = cyrus_core::kklt::compute_gv_corrected_target_tau(
                &c_i,
                &chi_divisor,
                c_tau,
                &gv_correction,
            ) else {
                eprintln!(
                    "[ERROR] GV-corrected KKLT target construction failed at iteration {iter}"
                );
                std::process::exit(2);
            };
            let Some(next) = solve_mixed_basis_path_following(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
                &kklt_basis,
                &tau_target,
                &previous_t,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] GV-corrected mixed-basis KKLT solve failed at iteration {iter}");
                std::process::exit(2);
            };
            if !next.converged {
                eprintln!(
                    "[ERROR] GV-corrected mixed-basis KKLT solve did not converge at iteration {iter}: rel_err={}",
                    next.relative_error.get()
                );
                std::process::exit(2);
            }
            let max_relative_step = next
                .t
                .iter()
                .zip(previous_t.iter())
                .map(|(new, old)| (new.get() - old.get()).abs() / (old.get().abs() + 1e-12))
                .fold(0.0f64, f64::max);
            eprintln!(
                "[INFO] GV-corrected KKLT iteration {iter}: max_relative_step={} rel_err={}",
                max_relative_step,
                next.relative_error.get()
            );
            corrected = next;
            correction_source_t = corrected.t.clone();
            if max_relative_step <= gv_tolerance {
                gv_converged = true;
                break;
            }
        }
        if !gv_converged {
            eprintln!(
                "[ERROR] GV-corrected KKLT fixed-point iteration did not converge in {max_gv_iterations} iterations"
            );
            std::process::exit(2);
        }
        let Some(gv_volume_correction) =
            cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                &small_curve_gvs,
                &intersection.basis,
                &corrected.t,
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
                &corrected.t,
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
                &corrected.t,
                None,
            );
        if input_chamber_gv_target_correction_no_gamma.is_none() {
            eprintln!(
                "[COMPARE] no-gamma input-chamber GV target correction is invalid at solved t"
            );
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
        compare_corrected_target_volume_checkpoint(
            data_dir,
            &kklt_basis,
            &chi_divisor,
            &tau_target,
            &input_chamber_gv_target_correction,
            input_chamber_gv_target_correction_no_gamma
                .as_deref()
                .map(|correction| ("no_gamma", correction)),
            &input_chamber_target_tau_finite,
        );
        if diagnose_chamber_updated_kklt {
            eprintln!(
                "[WARN] chamber-updated KKLT diagnostic is toric-covered-only when a chamber has missing toric GV values; it is not promoted to the production volume."
            );
            let diagnostic = diagnose_chamber_updated_kklt_toric_only(
                geom,
                intersection,
                &kklt_basis,
                &c_i,
                c_tau,
                &gamma,
                &corrected.t,
                small_curve_cutoff,
                kklt_steps,
                diagnose_chamber_updated_kklt_iterations,
            )
            .unwrap_or_else(|e| {
                eprintln!("[ERROR] chamber-updated KKLT diagnostic failed: {e}");
                std::process::exit(2);
            });
            let max_relative_t_delta = diagnostic
                .final_t
                .iter()
                .zip(corrected.t.iter())
                .map(|(new, old)| (new.get() - old.get()).abs() / (old.get().abs() + 1e-12))
                .fold(0.0f64, f64::max);
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
            let diagnostic_v_string = diagnostic.final_classical_volume - bbhl.get()
                + diagnostic.final_gv_volume_correction.get();
            let production_v_string =
                classical_volume_from_t(&intersection.kappa_basis, &corrected.t) - bbhl.get()
                    + gv_volume_correction.get();
            let final_first_missing_sparse = diagnostic
                .final_first_missing_class
                .as_deref()
                .map(sparse_i64)
                .unwrap_or_default();
            eprintln!(
                "[INFO] chamber-updated KKLT diagnostic final: iterations={} converged={} max_relative_t_delta_vs_production={} toric_missing={} first_missing_sparse={:?} V_classical={} GV_volume_correction={} V_string_partial={} delta_vs_production_partial={}",
                diagnostic.iterations,
                diagnostic.converged,
                max_relative_t_delta,
                diagnostic.final_toric_missing_count,
                final_first_missing_sparse,
                diagnostic.final_classical_volume,
                diagnostic.final_gv_volume_correction.get(),
                diagnostic_v_string,
                diagnostic_v_string - production_v_string
            );
        }
        (
            corrected.t,
            Some(gv_volume_correction),
            Some(gamma),
            Some(kklt_basis),
            Some(input_chamber_gv_target_correction),
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
    if diagnose_corrected_chamber_gv
        || diagnose_corrected_chamber_provided_generators_gv
        || diagnose_corrected_chamber_ray_gv
        || diagnose_corrected_chamber_lp_face_gv
    {
        let gamma = gamma_for_chamber_gv.as_deref().unwrap_or_else(|| {
            eprintln!("[ERROR] corrected-chamber GV diagnostic requires first-principles gamma");
            std::process::exit(2);
        });
        let kklt_basis = kklt_basis_for_chamber_gv.as_deref().unwrap_or_else(|| {
            eprintln!(
                "[ERROR] corrected-chamber GV diagnostic requires first-principles KKLT basis"
            );
            std::process::exit(2);
        });
        let diag = diagnose_chamber_gv_volume_correction(
            &corrected_chamber,
            geom,
            intersection,
            kklt_basis,
            &t,
            gamma,
            small_curve_cutoff,
            primal_gv_min_points,
            primal_gv_max_deg,
            diagnose_corrected_chamber_provided_generators_gv,
            diagnose_corrected_chamber_ray_gv,
            diagnose_corrected_chamber_lp_face_gv,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] corrected-chamber GV diagnostic failed: {e}");
            std::process::exit(2);
        });
        eprintln!(
            "[INFO] corrected-chamber small toric curves: ambient_rays={} subcutoff={} pair_pruned={} toric_covered={} toric_missing={} cutoff={}",
            diag.ambient_rays,
            diag.subcutoff_count,
            diag.filtered_count,
            diag.toric_gv_covered_count,
            diag.toric_gv_missing_count,
            small_curve_cutoff.get()
        );
        if let Some(general_covered) = diag.general_gv_covered_count {
            eprintln!(
                "[INFO] corrected-chamber general GV fallback covered {} missing curves; remaining_missing={}",
                general_covered, diag.remaining_gv_missing_count
            );
        }
        if let Some(ray_covered) = diag.ray_gv_covered_count {
            eprintln!(
                "[INFO] corrected-chamber one-dimensional ray GV diagnostic covered {} missing primitive generators; sample={:?}",
                ray_covered, diag.ray_gv_sample
            );
            if let Some(ray_correction) = diag.ray_gv_volume_correction.as_ref() {
                let ray_scope = if ray_covered < diag.toric_gv_missing_count {
                    "partial "
                } else {
                    ""
                };
                eprintln!(
                    "[INFO] corrected-chamber one-dimensional ray GV {ray_scope}volume correction (diagnostic) = {}",
                    ray_correction.get()
                );
                if let Some(input_chamber_correction) = gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber one-dimensional ray GV volume correction delta_vs_input_chamber (diagnostic) = {}",
                        ray_correction.get() - input_chamber_correction.get()
                    );
                }
            }
        }
        if let Some(face_covered) = diag.lp_face_gv_covered_count {
            eprintln!(
                "[INFO] corrected-chamber LP-witness face GV diagnostic covered {} decomposable missing curves; failed={:?}; sample={:?}",
                face_covered, diag.lp_face_gv_failed_count, diag.lp_face_gv_sample
            );
            if let Some(certified_count) = diag.lp_face_gv_certified_count {
                eprintln!(
                    "[INFO] corrected-chamber LP-witness face GV supporting-face certificates: certified={} uncertified={:?}",
                    certified_count, diag.lp_face_gv_uncertified_count
                );
            }
            if let Some(face_correction) = diag.lp_face_gv_volume_correction.as_ref() {
                eprintln!(
                    "[INFO] corrected-chamber LP-witness face GV partial volume correction (diagnostic) = {}",
                    face_correction.get()
                );
                if let Some(input_chamber_correction) = gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber LP-witness face GV volume correction delta_vs_input_chamber (diagnostic) = {}",
                        face_correction.get() - input_chamber_correction.get()
                    );
                }
            }
        }
        if let Some(combined_covered) = diag.combined_diagnostic_gv_covered_count {
            let combined_missing = diag.combined_diagnostic_gv_missing_count.unwrap_or(0);
            eprintln!(
                "[INFO] corrected-chamber combined diagnostic GV covered {}/{} toric-missing curves; remaining={}",
                combined_covered, diag.toric_gv_missing_count, combined_missing
            );
            if let (Some(zero_count), Some(nonzero_count)) = (
                diag.combined_diagnostic_gv_zero_count,
                diag.combined_diagnostic_gv_nonzero_count,
            ) {
                eprintln!(
                    "[INFO] corrected-chamber combined diagnostic GV values: zero={} nonzero={}",
                    zero_count, nonzero_count
                );
            }
            if let Some(combined_correction) =
                diag.combined_diagnostic_gv_volume_correction.as_ref()
            {
                let scope = if combined_missing == 0 {
                    ""
                } else {
                    "partial "
                };
                eprintln!(
                    "[INFO] corrected-chamber combined diagnostic GV {scope}volume correction (diagnostic, not promoted) = {}",
                    combined_correction.get()
                );
                if let Some(input_chamber_correction) = gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber combined diagnostic GV volume correction delta_vs_input_chamber (diagnostic, not promoted) = {}",
                        combined_correction.get() - input_chamber_correction.get()
                    );
                }
                if let Some(covered_correction) = diag.covered_gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber combined diagnostic GV volume correction delta_vs_toric_covered (diagnostic, not promoted) = {}",
                        combined_correction.get() - covered_correction.get()
                    );
                }
                if let (Some(input_target_correction), Some(combined_target_correction)) = (
                    input_chamber_gv_target_correction.as_ref(),
                    diag.combined_diagnostic_gv_target_correction.as_ref(),
                ) {
                    let summary = target_correction_delta_summary(
                        input_target_correction,
                        combined_target_correction,
                    )
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to compare combined diagnostic GV target correction: {e}"
                        );
                        std::process::exit(2);
                    });
                    eprintln!(
                        "[INFO] corrected-chamber combined diagnostic GV target correction delta_vs_input_chamber (diagnostic, not promoted) = {:?}",
                        summary
                    );
                }
            }
        }
        if let Some(ray_count) = diag.basis_mori_ray_count {
            let degree_window = match (
                diag.basis_mori_ray_degree_min,
                diag.basis_mori_ray_degree_max,
            ) {
                (Some(min_degree), Some(max_degree)) => {
                    format!("{min_degree}..{max_degree}")
                }
                _ => "unknown".to_string(),
            };
            eprintln!(
                "[INFO] corrected-chamber basis Mori rays for general GV: total={} degree_le_max={:?} degree_range={}",
                ray_count, diag.degree_bounded_basis_mori_ray_count, degree_window
            );
        }
        if let Some(stats) = diag.missing_target_stats.as_ref() {
            let active_generator_range = match (
                stats.real_cone_decomposition_active_generator_min,
                stats.real_cone_decomposition_active_generator_max,
            ) {
                (Some(min_active), Some(max_active)) => {
                    format!("{min_active}..{max_active}")
                }
                _ => "none".to_string(),
            };
            eprintln!(
                "[INFO] corrected-chamber missing GV target reduction: targets={} targets_as_mori_generators={} targets_as_origin_circuits={} real_cone_decomposable_by_other_generators={} lp_extremal_mori_generators={} real_cone_active_generator_range={} origin_circuit_resolved_conifold={} generators_le_target_degree_range={}..{} origin_coefficients={:?}",
                stats.target_count,
                stats.targets_that_are_mori_generators,
                stats.targets_that_are_origin_circuits,
                stats.targets_real_cone_decomposable_by_other_generators,
                stats.targets_that_are_lp_extremal_mori_generators,
                active_generator_range,
                stats.origin_circuit_resolved_conifold_count,
                stats.min_generators_le_target_degree,
                stats.max_generators_le_target_degree,
                stats.origin_coefficient_counts
            );
            eprintln!(
                "[INFO] corrected-chamber missing GV origin-circuit patterns: {:?}",
                stats.origin_circuit_pattern_counts
            );
            eprintln!(
                "[INFO] corrected-chamber missing GV target sample: {:?}",
                stats.sample
            );
        }
        if let Some(baseline) = diag
            .covered_toric_gv_divisor_representation_baseline
            .as_ref()
        {
            eprintln!(
                "[INFO] corrected-chamber covered toric GV divisor-representation baseline: checked_classes={} class_limit={} support_index_checks={} classes_with_support_solution={} first_without_support_solution={:?} sample={:?}",
                baseline.checked_class_count,
                baseline.class_limit,
                baseline.support_index_checks,
                baseline.classes_with_support_divisor_solution,
                baseline.first_without_support_divisor_solution,
                baseline.sample
            );
        }
        if let Some(first_missing) = &diag.first_missing_class {
            if let (Some(min_degree), Some(max_degree)) = (
                diag.missing_required_degree_min,
                diag.missing_required_degree_max,
            ) {
                eprintln!(
                    "[INFO] corrected-chamber missing toric GV required degree range: min={} max={} count={}",
                    min_degree, max_degree, diag.toric_gv_missing_count
                );
            }
            if let Some(chamber_correction) = diag.covered_gv_volume_correction.as_ref() {
                eprintln!(
                    "[INFO] corrected-chamber covered GV volume correction (partial) = {}",
                    chamber_correction.get()
                );
                if let Some(input_chamber_correction) = gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber covered GV volume correction delta_vs_input_chamber (partial) = {}",
                        chamber_correction.get() - input_chamber_correction.get()
                    );
                }
            }
            if let (Some(input_target_correction), Some(covered_target_correction)) = (
                input_chamber_gv_target_correction.as_ref(),
                diag.covered_gv_target_correction.as_ref(),
            ) {
                let summary = target_correction_delta_summary(
                    input_target_correction,
                    covered_target_correction,
                )
                .unwrap_or_else(|e| {
                    eprintln!("[ERROR] failed to compare covered GV target correction: {e}");
                    std::process::exit(2);
                });
                eprintln!(
                    "[INFO] corrected-chamber covered GV target correction delta_vs_input_chamber (partial) = {:?}",
                    summary
                );
            }
            eprintln!(
                "[WARN] corrected-chamber GV volume correction unavailable: missing GV values, first_missing={first_missing:?}"
            );
        } else if let Some(chamber_correction) = diag.gv_volume_correction.as_ref() {
            eprintln!(
                "[INFO] corrected-chamber GV volume correction = {}",
                chamber_correction.get()
            );
            if let Some(input_chamber_correction) = gv_volume_correction.as_ref() {
                eprintln!(
                    "[INFO] corrected-chamber GV volume correction delta_vs_input_chamber = {}",
                    chamber_correction.get() - input_chamber_correction.get()
                );
            }
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
    let v_string = v_string_before_gv + gv_volume_correction.map_or(0.0, |value| value.get());
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
                    "[WARN] corrected V_string comparison is not exact; this residual remains an unresolved instanton/chamber discrepancy, not a reproduced result"
                );
            }
            if abs_v > 0.1 {
                eprintln!(
                    "[ERROR] corrected V_string mismatch: got {v_string}, expected {corrected_v_expected}"
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

fn diagnostic_gv_value_counts(
    diagnostic_gvs: &HashMap<Vec<i64>, malachite::Integer>,
) -> (usize, usize) {
    let zero = malachite::Integer::from(0);
    let zero_count = diagnostic_gvs.values().filter(|gv| *gv == &zero).count();
    let nonzero_count = diagnostic_gvs.len().saturating_sub(zero_count);
    (zero_count, nonzero_count)
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
    manifest_dir: &PathBuf,
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
    let tol_v = 50.0;
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
        args.diagnose_corrected_chamber_gv,
        args.diagnose_corrected_chamber_provided_generators_gv,
        args.diagnose_corrected_chamber_ray_gv,
        args.diagnose_corrected_chamber_lp_face_gv,
        args.diagnose_chamber_updated_kklt,
        args.diagnose_chamber_updated_kklt_iterations,
        small_curve_cutoff,
        flat.dual_basis.len(),
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
    fn non_positive_basis_generator_degree_reports_first_bad_ray() {
        let rays = vec![vec![1, 0], vec![-2, 1], vec![0, 1]];
        let (count, first) = non_positive_basis_generator_degrees(&rays, &[1, 1]).unwrap();

        assert_eq!(count, 1);
        assert_eq!(first, Some((1, -1, vec![-2, 1])));
    }

    #[test]
    fn non_positive_basis_generator_degree_rejects_shape_mismatch() {
        assert!(
            non_positive_basis_generator_degrees(&[vec![1, 2, 3]], &[1, 1])
                .unwrap_err()
                .contains("basis ray length")
        );
    }

    #[test]
    fn graded_ray_stats_counts_degree_bounded_rays() {
        let stats =
            graded_ray_stats(&[vec![1, 0], vec![2, 1], vec![0, 5]], &[2, 3], Some(7)).unwrap();

        assert_eq!(
            stats,
            GradedRayStats {
                count: 3,
                degree_bounded_count: Some(2),
                min_degree: Some(2),
                max_degree: Some(15),
            }
        );
    }

    #[test]
    fn graded_ray_stats_rejects_shape_mismatch() {
        assert!(
            graded_ray_stats(&[vec![1, 2, 3]], &[1, 1], Some(3))
                .unwrap_err()
                .contains("basis ray length")
        );
    }

    #[test]
    fn degree_filtered_basis_rays_keeps_only_requested_degree_window() {
        let filtered =
            degree_filtered_basis_rays(&[vec![1, 0], vec![2, 1], vec![0, 5]], &[2, 3], 7).unwrap();

        assert_eq!(filtered, vec![vec![1, 0], vec![2, 1]]);
    }

    #[test]
    fn degree_filtering_rays_is_not_exact_for_cone_lattice_slice() {
        let rays = [vec![1, 0], vec![1, 100]];
        let grading = [1, 1];
        let filtered = degree_filtered_basis_rays(&rays, &grading, 2).unwrap();

        assert_eq!(filtered, vec![vec![1, 0]]);

        let low_degree_lattice_point = [1, 1];
        let degree = low_degree_lattice_point
            .iter()
            .zip(grading)
            .map(|(&coefficient, weight)| coefficient * weight)
            .sum::<i64>();
        assert_eq!(degree, 2);

        // The point is in cone((1,0),(1,100)):
        // (1,1) = 99/100 * (1,0) + 1/100 * (1,100).
        let full_cone_numerator = [99 * rays[0][0] + rays[1][0], 99 * rays[0][1] + rays[1][1]];
        assert_eq!(full_cone_numerator, [100, 100]);

        // But it is not in cone((1,0)), the degree-filtered cone. This is why
        // the exact corrected-chamber GV path cannot simply drop high-degree
        // Mori rays before computing the full cone inequalities.
        assert_ne!(low_degree_lattice_point[1], 0);
    }

    #[test]
    fn missing_gv_target_stats_reports_generator_membership_and_degree_window() {
        let ambient_classes = vec![vec![0, 1, 1], vec![0, 2, 0]];
        let basis_rays = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        let stats = missing_gv_target_stats(
            &ambient_classes,
            &basis_rays,
            &[1, 2],
            &[2, 3],
            0,
            &HashMap::new(),
            &HashMap::new(),
            4,
        )
        .unwrap();

        assert_eq!(stats.target_count, 2);
        assert_eq!(stats.targets_that_are_mori_generators, 1);
        assert_eq!(stats.targets_that_are_origin_circuits, 0);
        assert_eq!(stats.targets_real_cone_decomposable_by_other_generators, 2);
        assert_eq!(stats.targets_that_are_lp_extremal_mori_generators, 0);
        assert_eq!(stats.real_cone_decomposition_active_generator_min, Some(1));
        assert_eq!(stats.real_cone_decomposition_active_generator_max, Some(2));
        assert_eq!(stats.origin_circuit_resolved_conifold_count, 0);
        assert_eq!(stats.min_generators_le_target_degree, 2);
        assert_eq!(stats.max_generators_le_target_degree, 3);
        assert_eq!(stats.origin_coefficient_counts, BTreeMap::from([(0, 2)]));
        assert_eq!(stats.origin_circuit_pattern_counts, BTreeMap::new());
        assert_eq!(
            stats.sample,
            vec![
                MissingGvTargetSample {
                    degree: 5,
                    generators_le_degree: 3,
                    is_mori_generator: true,
                    origin_circuit_pattern: None,
                    origin_circuit_witness_count: None,
                    origin_circuit_first_witness: None,
                    cms_general_divisor_shape_candidates: None,
                    cms_general_divisor_intersection_checks: None,
                    real_cone_decomposable_by_other_generators: true,
                    real_cone_decomposition_active_generators: Some(2),
                    ambient_nonzero: vec![(1, 1), (2, 1)],
                    basis_nonzero: vec![(0, 1), (1, 1)]
                },
                MissingGvTargetSample {
                    degree: 4,
                    generators_le_degree: 2,
                    is_mori_generator: false,
                    origin_circuit_pattern: None,
                    origin_circuit_witness_count: None,
                    origin_circuit_first_witness: None,
                    cms_general_divisor_shape_candidates: None,
                    cms_general_divisor_intersection_checks: None,
                    real_cone_decomposable_by_other_generators: true,
                    real_cone_decomposition_active_generators: Some(1),
                    ambient_nonzero: vec![(1, 2)],
                    basis_nonzero: vec![(0, 2)]
                }
            ]
        );
    }

    #[test]
    fn cms_general_divisor_shape_candidates_report_negative_two_face_divisors() {
        let witness = cyrus_core::OriginCircuitCurveWitness {
            class: vec![-1, 0, 2, 1, 1, -3],
            first_facet_exclusive_point: 3,
            second_facet_exclusive_point: 4,
            shared_two_simplex: vec![2, 5, 1],
            first_facet: vec![1, 2, 3, 5],
            second_facet: vec![1, 2, 4, 5],
            sparse_relation: vec![(0, -1), (2, 2), (3, 1), (4, 1), (5, -3)],
            relation_points: vec![
                cyrus_core::OriginCircuitRelationPoint {
                    point_index: 0,
                    coefficient: -1,
                    coordinates: vec![0, 0, 0, 0],
                    face_dimension: Some(4),
                },
                cyrus_core::OriginCircuitRelationPoint {
                    point_index: 2,
                    coefficient: 2,
                    coordinates: vec![1, 2, 2, 2],
                    face_dimension: Some(2),
                },
                cyrus_core::OriginCircuitRelationPoint {
                    point_index: 3,
                    coefficient: 1,
                    coordinates: vec![2, 2, 1, 2],
                    face_dimension: Some(2),
                },
                cyrus_core::OriginCircuitRelationPoint {
                    point_index: 4,
                    coefficient: 1,
                    coordinates: vec![2, 3, 1, 3],
                    face_dimension: Some(2),
                },
                cyrus_core::OriginCircuitRelationPoint {
                    point_index: 5,
                    coefficient: -3,
                    coordinates: vec![2, 3, 2, 3],
                    face_dimension: Some(2),
                },
            ],
        };

        let diagnostic = cyrus_core::OriginCircuitCurveDiagnostic {
            class: witness.class.clone(),
            origin_coefficient: -1,
            negative_coefficient_counts: BTreeMap::from([(-3, 1)]),
            positive_coefficient_counts: BTreeMap::from([(1, 2), (2, 1)]),
            is_resolved_conifold_pattern: false,
            witnesses: vec![witness.clone(), witness],
        };
        let candidates = cms_general_divisor_shape_candidates(&diagnostic);

        assert_eq!(
            candidates,
            vec![CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 5,
                shrinking_divisor_coefficient: -3,
                shrinking_divisor_coordinates: vec![2, 3, 2, 3],
                inferred_other_normal_degree: 1,
                toric_gv1_formula_value: Some(3),
                all_non_origin_relation_points_are_two_face: true,
            }]
        );
    }

    #[test]
    fn rational_linear_system_finds_exact_overdetermined_solution() {
        let matrix = vec![
            vec![malachite::Rational::from(1), malachite::Rational::from(0)],
            vec![malachite::Rational::from(0), malachite::Rational::from(1)],
            vec![malachite::Rational::from(1), malachite::Rational::from(1)],
        ];
        let rhs = vec![
            malachite::Rational::from(2),
            malachite::Rational::from(3),
            malachite::Rational::from(5),
        ];

        let solution = solve_rational_linear_system(&matrix, &rhs).unwrap();

        assert_eq!(
            solution,
            Some(vec![
                malachite::Rational::from(2),
                malachite::Rational::from(3)
            ])
        );

        let inconsistent_rhs = vec![
            malachite::Rational::from(2),
            malachite::Rational::from(3),
            malachite::Rational::from(6),
        ];
        assert_eq!(
            solve_rational_linear_system(&matrix, &inconsistent_rhs).unwrap(),
            None
        );
    }

    #[test]
    fn cms_general_divisor_intersection_check_solves_divisor_class_and_normal_degree() {
        let mut kappa = cyrus_core::Intersection::new(3);
        kappa.set(
            2,
            1,
            1,
            cyrus_core::types::rational::Rational::<Finite>::from_raw(malachite::Rational::from(4)),
        );
        kappa.set(
            2,
            2,
            2,
            cyrus_core::types::rational::Rational::<Finite>::from_raw(malachite::Rational::from(
                -3,
            )),
        );
        let candidate = CmsGeneralDivisorShapeCandidate {
            shrinking_divisor_index: 2,
            shrinking_divisor_coefficient: -3,
            shrinking_divisor_coordinates: vec![2, 3, 2, 3],
            inferred_other_normal_degree: 1,
            toric_gv1_formula_value: Some(3),
            all_non_origin_relation_points_are_two_face: true,
        };

        let check =
            cms_general_divisor_intersection_check(&[0, 4, -3], &candidate, &kappa, &[1, 2])
                .unwrap();

        assert_eq!(
            check,
            CmsGeneralDivisorIntersectionCheck {
                shrinking_divisor_index: 2,
                has_rational_divisor_solution: true,
                solution_basis_support_len: Some(2),
                solution_is_integral: Some(true),
                computed_other_normal_degree: Some("1".to_string()),
                matches_inferred_other_normal_degree: Some(true),
            }
        );
    }

    #[test]
    fn covered_toric_gv_divisor_representation_baseline_checks_support_divisors() {
        let mut kappa = cyrus_core::Intersection::new(3);
        kappa.set(
            2,
            1,
            1,
            cyrus_core::types::rational::Rational::<Finite>::from_raw(malachite::Rational::from(4)),
        );
        kappa.set(
            2,
            2,
            2,
            cyrus_core::types::rational::Rational::<Finite>::from_raw(malachite::Rational::from(
                -3,
            )),
        );
        let covered_classes = vec![(vec![0, 4, -3], malachite::Integer::from(3))];

        let baseline = compute_covered_toric_gv_divisor_representation_baseline(
            &covered_classes,
            0,
            &kappa,
            &[1, 2],
            1,
        )
        .unwrap()
        .unwrap();

        assert_eq!(
            baseline,
            CoveredToricGvDivisorRepresentationBaseline {
                checked_class_count: 1,
                class_limit: 1,
                support_index_checks: 2,
                classes_with_support_divisor_solution: 1,
                first_without_support_divisor_solution: None,
                sample: vec![CoveredToricGvDivisorRepresentationSample {
                    ambient_nonzero: vec![(1, 4), (2, -3)],
                    support_indices_with_solution: vec![1, 2],
                    support_indices_without_solution: vec![],
                }],
            }
        );
    }

    #[test]
    fn one_dimensional_ray_gv_targets_keep_only_lp_extremal_generators() {
        let targets = one_dimensional_ray_gv_targets(
            &[vec![0, 1, 0], vec![0, 1, 1], vec![0, 2, 0]],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &[1, 2],
            &[2, 3],
        )
        .unwrap();

        assert_eq!(
            targets,
            OneDimensionalRayGvTargets {
                candidates: vec![vec![0, 1, 0]],
                skipped_non_generators: 1,
                skipped_decomposable_generators: 1,
            }
        );
    }

    #[test]
    fn degree_bounded_span_generator_count_uses_exact_rank() {
        let count = degree_bounded_span_generators(
            &[vec![1, 0]],
            &[vec![1, 0], vec![2, 0], vec![0, 1], vec![1, 1]],
            &[1, 1],
            2,
        )
        .unwrap()
        .len();

        assert_eq!(count, 2);
    }

    #[test]
    fn degree_bounded_face_lattice_points_saturates_non_unimodular_cone() {
        let points =
            degree_bounded_face_lattice_points(&[vec![1, 0], vec![1, 2]], &[1, 1], 3).unwrap();

        assert!(points.contains(&vec![1, 1]));
        assert!(points.contains(&vec![1, 0]));
        assert!(points.contains(&vec![1, 2]));
    }

    #[test]
    fn supporting_mori_face_certificate_finds_exact_integer_normal() {
        let certificate =
            certify_supporting_mori_face(&[vec![1, 0]], &[vec![1, 0], vec![0, 1], vec![1, 1]])
                .unwrap()
                .unwrap();

        assert_eq!(exact_i64_dot(&certificate.normal, &[1, 0]), 0);
        assert!(exact_i64_dot(&certificate.normal, &[0, 1]) > 0);
        assert_eq!(certificate.zero_generator_count, 1);
        assert_eq!(certificate.positive_generator_count, 2);
    }

    #[test]
    fn supporting_mori_face_certificate_rejects_interior_ray() {
        let certificate =
            certify_supporting_mori_face(&[vec![1, 1]], &[vec![1, 0], vec![0, 1], vec![1, 1]])
                .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn exact_generator_decomposition_finds_integer_sums() {
        let mut decompositions = exact_generator_decompositions(
            &[2, 1],
            &[vec![1, 0], vec![1, 1], vec![0, 1], vec![2, 1]],
            &[1, 1],
            3,
            2,
            1,
        )
        .unwrap()
        .unwrap();
        let decomposition = decompositions.pop().unwrap();

        assert_eq!(decomposition, vec![vec![1, 0], vec![1, 1]]);
    }

    #[test]
    fn decomposition_diamond_elements_are_sub_sums_of_exact_decomposition() {
        let elements =
            decomposition_diamond_elements(&[vec![1, 0], vec![1, 0], vec![0, 1]], &[2, 1]).unwrap();

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
    fn insert_missing_diagnostic_gv_rejects_conflicting_duplicates() {
        let mut diagnostic_gvs = HashMap::new();
        insert_missing_diagnostic_gv(
            &mut diagnostic_gvs,
            &[0, 2, -1],
            &malachite::Integer::from(3),
            "test-a",
        )
        .unwrap();
        insert_missing_diagnostic_gv(
            &mut diagnostic_gvs,
            &[0, 2, -1],
            &malachite::Integer::from(3),
            "test-b",
        )
        .unwrap();

        let err = insert_missing_diagnostic_gv(
            &mut diagnostic_gvs,
            &[0, 2, -1],
            &malachite::Integer::from(4),
            "test-c",
        )
        .unwrap_err();

        assert!(err.contains("conflicting corrected-chamber diagnostic GV values"));
        assert_eq!(
            diagnostic_gvs.get(&vec![0, 2, -1]),
            Some(&malachite::Integer::from(3))
        );
    }

    #[test]
    fn diagnostic_gv_value_counts_separates_zero_and_nonzero_values() {
        let diagnostic_gvs = HashMap::from([
            (vec![1, 0], malachite::Integer::from(0)),
            (vec![0, 1], malachite::Integer::from(4)),
            (vec![1, 1], malachite::Integer::from(-2)),
        ]);

        assert_eq!(diagnostic_gv_value_counts(&diagnostic_gvs), (1, 2));
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

    #[test]
    fn exact_decomposition_closure_recurses_into_intermediate_elements() {
        let elements = exact_decomposition_closure_elements(
            &[2, 1],
            &[vec![1, 0], vec![0, 1], vec![1, 1], vec![2, 1]],
            &[1, 1],
            2,
            8,
            16,
        )
        .unwrap();

        assert_eq!(
            elements,
            vec![vec![0, 0], vec![0, 1], vec![1, 0], vec![1, 1], vec![2, 1],]
        );
    }
}

fn main() {
    run_pipeline(parse_args());
}
