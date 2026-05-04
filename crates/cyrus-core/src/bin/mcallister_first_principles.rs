//! First-principles McAllister pipeline runner (opt-in).
//!
//! Runs the pipeline stages explicitly and fails fast on any mismatch or
//! invalid physics. This is intentionally a binary, not a test.
//!
//! Optional:
//! - `--dual-basis path/to/dual_basis.json` to supply the flux coordinate basis.
//! - `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.
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

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
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

fn validate_basis_checkpoint(
    glsm: &[Vec<malachite::Integer>],
    computed: &[usize],
    data_dir: &str,
    label: &str,
) {
    let basis_path = PathBuf::from(data_dir).join("basis.dat");
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
    let expected_dual_points = read_points(&dir.join("dual_points.dat"));
    let expected_dual_set = {
        let mut points = expected_dual_points;
        points.sort();
        points
    };
    let actual_dual_set = sorted_point_coords(dual_polytope.vertices());
    if actual_dual_set != expected_dual_set {
        eprintln!("[ERROR] computed dual polytope points differ from dual_points.dat checkpoint");
        std::process::exit(2);
    }

    let expected_simplices_i64 = read_points(&dir.join("dual_simplices.dat"));
    let mut expected_simplices: Vec<Vec<usize>> = expected_simplices_i64
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|value| {
                    usize::try_from(value).expect("dual_simplices.dat index must be non-negative")
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
    eprintln!(
        "[INFO] computed dual polytope/FRST match dual_points.dat and dual_simplices.dat checkpoints"
    );
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
    lp_face_gv_volume_correction: Option<F64<Finite>>,
    lp_face_gv_sample: Vec<FaceGvDiagnosticSample>,
    remaining_gv_missing_count: usize,
    first_missing_class: Option<Vec<i64>>,
    missing_required_degree_min: Option<i128>,
    missing_required_degree_max: Option<i128>,
    missing_target_stats: Option<MissingGvTargetStats>,
    covered_gv_volume_correction: Option<F64<Finite>>,
    gv_volume_correction: Option<F64<Finite>>,
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
    real_cone_decomposable_by_other_generators: bool,
    real_cone_decomposition_active_generators: Option<usize>,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
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
    gv: Option<malachite::Integer>,
    error: Option<String>,
    ambient_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug)]
struct FaceGvDiagnosticResult {
    ambient_class: Vec<i64>,
    degree: i128,
    generator_count: usize,
    active_generator_count: usize,
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
    eprintln!(
        "[INFO] transforming Kähler parameters from source basis {:?} to computed basis {:?}",
        source_basis, computed_basis
    );
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
    let flux_basis =
        dual_basis_override.map_or_else(|| vec![3, 4, 5, 8], |basis| basis.indices.clone());
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
        let origin_circuit_pattern_label =
            origin_circuits_by_class
                .get(ambient_class)
                .map(|diagnostic| {
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
        min_generators = min_generators.min(generators_le_degree);
        max_generators = max_generators.max(generators_le_degree);
        if sample.len() < sample_limit {
            sample.push(MissingGvTargetSample {
                degree,
                generators_le_degree,
                is_mori_generator,
                origin_circuit_pattern: origin_circuit_pattern_label,
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
        return Ok(None);
    }

    let mut vars = ProblemVariables::new();
    let coefficients: Vec<_> = candidate_indices
        .iter()
        .map(|_| vars.add(variable().min(0.0)))
        .collect();
    let mut objective = Expression::from(0.0);
    for coefficient in &coefficients {
        objective.add_mul(1.0, *coefficient);
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
        let Some(witness) = real_cone_decomposition_by_other_degree_bounded_generators(
            &basis_class,
            basis_rays,
            grading,
            degree,
        )?
        else {
            continue;
        };
        let max_deg = u32::try_from(degree).map_err(|_| {
            format!("LP-witness face GV target degree {degree} does not fit in u32")
        })?;
        let mut provided_generators: Vec<Vec<i64>> = witness
            .active_generator_indices
            .iter()
            .map(|&idx| basis_rays[idx].clone())
            .collect();
        provided_generators.push(basis_class.clone());
        provided_generators.sort();
        provided_generators.dedup();

        let previous_panic_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            cyrus_core::compute_gv_invariants_with_provided_generators(
                &provided_generators,
                grading,
                q_matrix,
                intnums,
                None,
                Some(max_deg),
            )
        }));
        std::panic::set_hook(previous_panic_hook);

        let (gv, error) = match gvs_result {
            Ok(Ok(gvs)) => {
                let target_i32 = basis_class
                    .iter()
                    .map(|&value| {
                        i32::try_from(value).map_err(|_| {
                            "LP-witness face GV target coordinate does not fit in i32".to_string()
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let gv = gvs
                    .into_iter()
                    .find_map(|(curve, gv)| (curve == target_i32).then_some(gv))
                    .unwrap_or_else(|| malachite::Integer::from(0));
                (Some(gv), None)
            }
            Ok(Err(e)) => (
                None,
                Some(format!("failed LP-witness face GV computation: {e}")),
            ),
            Err(payload) => (
                None,
                Some(format!(
                    "LP-witness face GV computation panicked for target degree {degree}, ambient_nonzero={:?}: {}",
                    sparse_i64(ambient_class),
                    panic_payload_message(payload.as_ref())
                )),
            ),
        };

        out.push(FaceGvDiagnosticResult {
            ambient_class: ambient_class.clone(),
            degree,
            generator_count: provided_generators.len(),
            active_generator_count: witness.active_generator_indices.len(),
            gv,
            error,
        });
    }
    Ok(out)
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
    Ok(intersection_in_basis(&kappa_full, basis))
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

    let mut basis_ray_stats = None;
    let mut basis_rays_for_missing = None;
    let mut grading_for_missing = None;
    let mut degree_summary = None;
    let mut missing_target_stats = None;
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
    let mut lp_face_gv_sample = Vec::new();
    let mut lp_face_gv_volume_correction = None;
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
            "[WARN] corrected-chamber LP-witness face GV diagnostic uses small provided-generator cones from floating LP witnesses; this is not yet promoted to the exact corrected-chamber GV fallback."
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
        lp_face_gv_covered_count = Some(covered_count);
        lp_face_gv_failed_count = Some(failed_count);
        lp_face_gv_sample = face_gvs
            .iter()
            .take(10)
            .map(|result| FaceGvDiagnosticSample {
                degree: result.degree,
                generator_count: result.generator_count,
                active_generator_count: result.active_generator_count,
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
                "corrected-chamber general GV direct fallback would dualize {} Mori generators, exceeding limit {direct_ray_limit}; {bounded_note}. Current backend computes cone hyperplanes by DDM before bounded lattice enumeration. Need a reduced corrected-chamber cone/lattice formulation, or set CYRUS_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT to force the direct attempt.",
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
        lp_face_gv_volume_correction,
        lp_face_gv_sample,
        remaining_gv_missing_count: missing_gv_classes.len(),
        first_missing_class: missing_gv_classes.first().cloned(),
        missing_required_degree_min,
        missing_required_degree_max,
        missing_target_stats,
        covered_gv_volume_correction,
        gv_volume_correction,
    })
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
        || diagnose_corrected_chamber_lp_face_gv)
        && allow_downstream_kahler
    {
        eprintln!(
            "[ERROR] --diagnose-corrected-chamber-gv is only valid for first-principles runs, not downstream Kähler replay"
        );
        std::process::exit(2);
    }
    let (t, gv_volume_correction, gamma_for_chamber_gv) = if allow_downstream_kahler {
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
        (corrected.t, Some(gv_volume_correction), Some(gamma))
    };

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
    if diagnose_corrected_chamber_gv
        || diagnose_corrected_chamber_provided_generators_gv
        || diagnose_corrected_chamber_ray_gv
        || diagnose_corrected_chamber_lp_face_gv
    {
        let gamma = gamma_for_chamber_gv.as_deref().unwrap_or_else(|| {
            eprintln!("[ERROR] corrected-chamber GV diagnostic requires first-principles gamma");
            std::process::exit(2);
        });
        let diag = diagnose_chamber_gv_volume_correction(
            &corrected_chamber,
            geom,
            intersection,
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
        let g_s_expected: f64 = std::fs::read_to_string(dir.join("g_s.dat"))
            .expect("Failed to read g_s.dat")
            .trim()
            .parse()
            .expect("Invalid g_s.dat");
        let w0_expected: f64 = std::fs::read_to_string(dir.join("W_0.dat"))
            .expect("Failed to read W_0.dat")
            .trim()
            .parse()
            .expect("Invalid W_0.dat");
        let rel_gs = ((g_s.get() - g_s_expected) / g_s_expected).abs();
        let rel_w0 = ((w0.get() - w0_expected) / w0_expected).abs();
        eprintln!("[COMPARE] g_s rel_err = {rel_gs}");
        eprintln!("[COMPARE] W0 rel_err = {rel_w0}");
        let corrected_volume_path = dir.join("corrected_cy_vol.dat");
        if corrected_volume_path.exists() {
            let corrected_v_expected: f64 = std::fs::read_to_string(&corrected_volume_path)
                .expect("Failed to read corrected_cy_vol.dat")
                .trim()
                .parse()
                .expect("Invalid corrected_cy_vol.dat");
            let abs_v = (v_string - corrected_v_expected).abs();
            eprintln!("[COMPARE] corrected V_string abs_err = {abs_v}");
            if abs_v > 0.1 {
                eprintln!(
                    "[ERROR] corrected V_string mismatch: got {v_string}, expected {corrected_v_expected}"
                );
                std::process::exit(2);
            }
        }
    }
}

fn stage_vacuum(
    ek0_opt: Option<F64<Pos>>,
    racetrack: &RacetrackData,
    v_string_pos: F64<Pos>,
    v_string: f64,
    compare_dir: Option<String>,
    data_dir: Option<String>,
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
    let assertion_path = manifest_dir.join("tests/mcallister_e2e/assertions/racetrack.json");
    let assertion: RacetrackAssertion = load_json(&assertion_path);
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
    if summary.w0 <= 0.0 {
        eprintln!("[ERROR] W0 must be positive");
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
                    real_cone_decomposable_by_other_generators: true,
                    real_cone_decomposition_active_generators: Some(1),
                    ambient_nonzero: vec![(1, 2)],
                    basis_nonzero: vec![(0, 2)]
                }
            ]
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
}

fn main() {
    run_pipeline(parse_args());
}
