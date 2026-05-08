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
//! - `--dump-corrected-chamber-gv-context path` to export the corrected-chamber
//!   missing-GV context as JSON for source-level CYTools/cygv reconstruction.
//! - `--dump-corrected-chamber-secondary-certificate path` to export only the
//!   corrected-chamber secondary-cone height certificate as JSON without
//!   building the full missing-GV context.
//! - `--dump-corrected-chamber-2face-secondary-certificate path` to export the
//!   CYTools-style 4D two-face-restricted secondary-cone certificate.
//! - `--dump-corrected-chamber-expanded-secondary-fan-certificate path` to
//!   export the CYTools-style expanded-secondary support certificate for the
//!   corrected-chamber height vector.
//! - `--dump-corrected-chamber-face-triangulation-choice-summary path` to
//!   export exact/sampled 4D two-face FRT choice counts used by the expanded
//!   chamber search machinery. Sampling controls are
//!   `--face-triangulation-max-exact-points`,
//!   `--face-triangulation-samples-per-large-face`,
//!   `--face-triangulation-max-sampling-attempts-per-face`, and
//!   `--face-triangulation-seed`.
//! - `--diagnose-chamber-updated-kklt` to run a diagnostic-only KKLT
//!   fixed-point loop that recomputes the FRST chamber, intersections, divisor
//!   χ, and toric-covered small-curve GV target correction at each iteration.
//! - `--small-curve-pruning <pair|finite-semigroup>` to choose the selected
//!   toric curve pruning rule. The default `pair` reproduces McAllister
//!   checkpoints; `finite-semigroup` is stricter for GA/search diagnostics.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::path::{Path, PathBuf};
use std::time::Instant;

use good_lp::{
    Expression, ProblemVariables, ResolutionError, Solution, SolverModel, default_solver, variable,
};

use cyrus_core::flat_direction::{compute_flat_direction, compute_flat_direction_full};
use cyrus_core::gv::{
    BoundedCurveDecompositionIndex, certify_supporting_mori_face_by_exact_kernel,
    check_supporting_mori_face_normal,
};
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
    compute_origin_circuit_curve_diagnostics, compute_regular_triangulation,
    compute_toric_curve_gv_diagnostics, compute_toric_two_face_curve_gv_invariants,
    compute_w0_from_terms, divisor_basis_change_matrix, effective_prime_divisors_from_curve_basis,
    expanded_secondary_fan_hyperplanes_on_polytope_2faces_4d,
    fine_regular_triangulation_choices_on_polytope_2faces_4d_with_sampling,
    generate_scaled_divisor_basis_branch_initializations, gv_divisor_basis_data, heights_to_kahler,
    intersection_in_basis, intersection_in_divisor_basis, is_unimodular, kahler_to_heights,
    map_basis_gv_invariants_to_ambient, project_ambient_curve_to_basis,
    prune_decomposable_curve_candidates, scale_divisor_basis_kklt_branch_initialization_to_target,
    secondary_cone_height_pairings, secondary_cone_hyperplanes_native,
    secondary_cone_hyperplanes_native_on_polytope_2faces_4d, solve_divisor_basis_path_following,
    solve_divisor_basis_path_following_branch_candidates, solve_mixed_basis_path_following,
    solve_racetrack, subcutoff_toric_curve_candidates,
};

const DEFAULT_MCALLISTER_GV_MIN_POINTS: u32 = 20_000;
const DEFAULT_CORRECTED_CHAMBER_GENERAL_GV_DIRECT_RAY_LIMIT: usize = 100_000;
const DEFAULT_CORRECTED_CHAMBER_PROVIDED_GV_GENERATOR_LIMIT: usize = 2_000;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_SPAN_GENERATOR_LIMIT: usize = 64;
const DEFAULT_FACE_TRIANGULATION_MAX_EXACT_POINTS: usize = 17;
const DEFAULT_FACE_TRIANGULATION_SAMPLES_PER_LARGE_FACE: usize = 1_000;
const DEFAULT_FACE_TRIANGULATION_MAX_SAMPLING_ATTEMPTS_PER_FACE: usize = 100_000;
const DEFAULT_FACE_TRIANGULATION_SEED: u64 = 0;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_LATTICE_GENERATOR_LIMIT: usize = 64;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_LATTICE_ELEMENT_LIMIT: usize = 50_000;
const DEFAULT_CORRECTED_CHAMBER_LP_FACE_INTEGER_DECOMPOSITION_MAX_TERMS: usize = 4;
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

    fn dimension(&self) -> usize {
        match self {
            Self::Indices(indices) => indices.len(),
            Self::Matrix { basis_matrix, .. } => basis_matrix.len(),
        }
    }

    fn description(&self) -> &'static str {
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
    formula_sum_diagnostic_gv_covered_count: Option<usize>,
    formula_sum_diagnostic_gv_missing_count: Option<usize>,
    formula_sum_diagnostic_gv_zero_count: Option<usize>,
    formula_sum_diagnostic_gv_nonzero_count: Option<usize>,
    formula_sum_diagnostic_gv_volume_correction: Option<F64<Finite>>,
    formula_sum_diagnostic_gv_target_correction: Option<Vec<F64<Finite>>>,
    remaining_gv_missing_count: usize,
    first_missing_class: Option<Vec<i64>>,
    missing_required_degree_min: Option<i128>,
    missing_required_degree_max: Option<i128>,
    missing_target_stats: Option<MissingGvTargetStats>,
    uncovered_source_ray_stats_degree_bound_for_missing: Option<i128>,
    uncovered_source_ray_stats_for_missing: Option<MissingGvTargetStats>,
    shared_facet_unresolved_source_ray_stats_for_missing: Option<MissingGvTargetStats>,
    uncovered_source_ray_toric_diagnostic_sample: Option<Vec<ToricGvDiagnosticContextSample>>,
    degree_bounded_toric_gv_diagnostic_context_for_missing:
        Option<Vec<ToricGvDiagnosticContextSample>>,
    secondary_cone_height_certificate: Option<SecondaryConeHeightCertificate>,
    secondary_cone_2face_height_certificate: Option<SecondaryConeHeightCertificate>,
    expanded_secondary_fan_height_certificate: Option<SecondaryConeHeightCertificate>,
    secondary_cone_heights_for_missing: Option<Vec<f64>>,
    basis_mori_rays_for_missing_degree_bound: Option<i128>,
    basis_mori_rays_for_missing_degree_bounded: Option<Vec<Vec<i64>>>,
    degree_bounded_mori_ray_context_for_missing: Option<Vec<DegreeBoundedMoriRayContextSample>>,
    covered_toric_gv_context_for_missing: Option<Vec<CoveredToricGvContextSample>>,
    gv_q_matrix_for_missing: Option<Vec<Vec<i64>>>,
    gv_curve_basis_matrix_for_missing: Option<Vec<Vec<String>>>,
    grading_for_missing: Option<Vec<i64>>,
    corrected_kappa_basis_for_missing: Option<Vec<SparseIntersectionEntry>>,
    covered_toric_gv_divisor_representation_baseline:
        Option<CoveredToricGvDivisorRepresentationBaseline>,
    covered_gv_target_correction: Option<Vec<F64<Finite>>>,
    covered_gv_volume_correction: Option<F64<Finite>>,
    gv_volume_correction: Option<F64<Finite>>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
struct SecondaryConeHeightCertificate {
    status: String,
    epsilon: f64,
    hyperplane_count: usize,
    pairing_count: usize,
    min_pairing: Option<f64>,
    max_pairing: Option<f64>,
    strictly_inside: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct FaceTriangulationChoiceSummary {
    status: String,
    max_exact_face_points: usize,
    samples_per_large_face: usize,
    max_sampling_attempts_per_face: usize,
    seed: u64,
    face_count: usize,
    exact_face_count: usize,
    sampled_face_count: usize,
    empty_choice_face_count: usize,
    min_face_points: Option<usize>,
    max_face_points: Option<usize>,
    min_choice_count: Option<usize>,
    max_choice_count: Option<usize>,
    total_choice_count: String,
    face_point_counts: Vec<usize>,
    choice_counts: Vec<usize>,
    sampled_face_indices: Vec<usize>,
    sampled_face_point_counts: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct SparseIntersectionEntry {
    indices: [usize; 3],
    value: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct DegreeBoundedMoriRayContextSample {
    degree: i128,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct CoveredToricGvContextSample {
    degree: i128,
    gv: String,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct ToricGvDiagnosticContextSample {
    degree: i128,
    gv: String,
    source_bucket: String,
    source_summary: String,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Serialize)]
struct CorrectedChamberGvContextExport<'a> {
    schema_version: u32,
    source: &'static str,
    small_curve_cutoff: f64,
    small_curve_pruning: &'static str,
    ambient_rays: usize,
    subcutoff_count: usize,
    filtered_count: usize,
    toric_gv_covered_count: usize,
    toric_gv_missing_count: usize,
    remaining_gv_missing_count: usize,
    first_missing_class: Option<&'a Vec<i64>>,
    missing_required_degree_min: Option<i128>,
    missing_required_degree_max: Option<i128>,
    basis_mori_ray_count: Option<usize>,
    degree_bounded_basis_mori_ray_count: Option<usize>,
    basis_mori_ray_degree_min: Option<i128>,
    basis_mori_ray_degree_max: Option<i128>,
    basis_mori_rays_for_missing_degree_bound: Option<i128>,
    basis_mori_rays_for_missing_degree_bounded: Option<&'a Vec<Vec<i64>>>,
    degree_bounded_mori_ray_context_for_missing: Option<&'a Vec<DegreeBoundedMoriRayContextSample>>,
    covered_toric_gv_context_for_missing: Option<&'a Vec<CoveredToricGvContextSample>>,
    degree_bounded_toric_gv_diagnostic_context_for_missing:
        Option<&'a Vec<ToricGvDiagnosticContextSample>>,
    secondary_cone_height_certificate: Option<&'a SecondaryConeHeightCertificate>,
    secondary_cone_2face_height_certificate: Option<&'a SecondaryConeHeightCertificate>,
    expanded_secondary_fan_height_certificate: Option<&'a SecondaryConeHeightCertificate>,
    secondary_cone_heights_for_missing: Option<&'a Vec<f64>>,
    gv_q_matrix_for_missing: Option<&'a Vec<Vec<i64>>>,
    gv_curve_basis_matrix_for_missing: Option<&'a Vec<Vec<String>>>,
    grading_for_missing: Option<&'a Vec<i64>>,
    corrected_kappa_basis_for_missing: Option<&'a Vec<SparseIntersectionEntry>>,
    missing_target_sample_limit: usize,
    missing_target_sample_is_complete: Option<bool>,
    missing_target_stats: Option<&'a MissingGvTargetStats>,
    formula_sum_diagnostic_gv_covered_count: Option<usize>,
    formula_sum_diagnostic_gv_missing_count: Option<usize>,
    formula_sum_diagnostic_gv_zero_count: Option<usize>,
    formula_sum_diagnostic_gv_nonzero_count: Option<usize>,
    formula_sum_diagnostic_gv_volume_correction: Option<f64>,
    formula_sum_diagnostic_gv_target_correction: Option<Vec<f64>>,
    uncovered_source_ray_stats_degree_bound_for_missing: Option<i128>,
    uncovered_source_ray_stats_for_missing: Option<&'a MissingGvTargetStats>,
    shared_facet_unresolved_source_ray_stats_for_missing: Option<&'a MissingGvTargetStats>,
    uncovered_source_ray_toric_diagnostic_sample: Option<&'a Vec<ToricGvDiagnosticContextSample>>,
}

struct ChamberToricGvSelection {
    ambient_rays: usize,
    subcutoff_count: usize,
    filtered_count: usize,
    subcutoff_toric_gv_covered_count: usize,
    subcutoff_toric_gv_missing_count: usize,
    toric_gv_covered_count: usize,
    toric_gv_missing_count: usize,
    first_missing_class: Option<Vec<i64>>,
    subcutoff_missing_gv_classes: Vec<Vec<i64>>,
    missing_gv_classes: Vec<Vec<i64>>,
    small_curve_candidates: Vec<ToricCurveCandidate>,
    small_curves: Vec<ToricCurveCandidate>,
    subcutoff_curve_gvs: Vec<(Vec<i64>, malachite::Integer)>,
    small_curve_gvs: Vec<(Vec<i64>, malachite::Integer)>,
}

struct AmbientTargetContributionRow {
    class: Vec<i64>,
    gv: malachite::Integer,
    q_divisor: i64,
    q_dot_t: f64,
    parity: i128,
    contribution: f64,
}

#[derive(Clone, Debug)]
struct TopToricLocalGvTarget {
    kklt_idx: usize,
    point_idx: usize,
    class: Vec<i64>,
    toric_gv: malachite::Integer,
}

struct GvBranchBucketAccumulator {
    count: usize,
    q_dot_min: f64,
    q_dot_max: f64,
    correction: Vec<f64>,
}

struct GvBranchBucketReport {
    parity_mod2: i128,
    q_dot_bucket: &'static str,
    count: usize,
    q_dot_min: f64,
    q_dot_max: f64,
    max_abs_correction: f64,
    l2_correction: f64,
    correction: Vec<f64>,
}

struct GvSourceBucketAccumulator {
    count: usize,
    source_count_min: usize,
    source_count_max: usize,
    q_dot_min: f64,
    q_dot_max: f64,
    correction: Vec<f64>,
}

struct GvSourceBucketReport {
    label: String,
    count: usize,
    source_count_min: usize,
    source_count_max: usize,
    q_dot_min: f64,
    q_dot_max: f64,
    max_abs_correction: f64,
    l2_correction: f64,
    correction: Vec<f64>,
}

struct GvWeightLpReport {
    label: &'static str,
    lower_bound: f64,
    upper_bound: f64,
    objective_epsilon: f64,
    max_abs_delta: f64,
    relative_l2_delta: f64,
    weight_min: f64,
    weight_max: f64,
    weight_mean: f64,
    near_lower_count: usize,
    near_upper_count: usize,
    interior_count: usize,
    max_abs_delta_from_one: f64,
    weights: Vec<f64>,
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

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
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
    origin_circuit_affine_rank_counts: BTreeMap<usize, usize>,
    branch_status_counts: BTreeMap<String, usize>,
    branch_bucket_counts: BTreeMap<String, usize>,
    real_cone_decomposition_exact_kind_counts: BTreeMap<String, usize>,
    sample: Vec<MissingGvTargetSample>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct MissingGvTargetSample {
    degree: i128,
    generators_le_degree: usize,
    is_mori_generator: bool,
    origin_circuit_pattern: Option<String>,
    origin_circuit_witness_count: Option<usize>,
    origin_circuit_first_witness: Option<OriginCircuitWitnessSample>,
    origin_circuit_witnesses: Option<Vec<OriginCircuitWitnessSample>>,
    origin_circuit_affine_support: Option<OriginCircuitAffineSupportSample>,
    cms_general_divisor_shape_candidates: Option<Vec<CmsGeneralDivisorShapeCandidate>>,
    cms_general_divisor_intersection_checks: Option<Vec<CmsGeneralDivisorIntersectionCheck>>,
    branch_diagnostic: Option<MissingGvBranchDiagnostic>,
    real_cone_decomposable_by_other_generators: bool,
    real_cone_decomposition_active_generators: Option<usize>,
    real_cone_decomposition_active_generator_basis_nonzero: Option<Vec<Vec<(usize, i64)>>>,
    real_cone_decomposition_exact_coefficients: Option<Vec<String>>,
    real_cone_decomposition_exact_kind: Option<&'static str>,
    ambient_nonzero: Vec<(usize, i64)>,
    basis_nonzero: Vec<(usize, i64)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct MissingGvBranchDiagnostic {
    q_dot_t: String,
    parity: i128,
    parity_mod2: i128,
    q_dot_bucket: &'static str,
    dilog_status: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct OriginCircuitAffineSupportSample {
    affine_rank: usize,
    coefficient_counts: BTreeMap<i64, usize>,
    local_charge_basis: Vec<Vec<i64>>,
    local_coordinates: Vec<OriginCircuitLocalCoordinateSample>,
    local_coordinates_2d: Option<Vec<OriginCircuitLocalCoordinate2DSample>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct OriginCircuitLocalCoordinateSample {
    point_index: usize,
    coordinates: Vec<i64>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct OriginCircuitLocalCoordinate2DSample {
    point_index: usize,
    coordinates: [i64; 2],
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct CmsGeneralDivisorShapeCandidate {
    shrinking_divisor_index: usize,
    shrinking_divisor_coefficient: i64,
    shrinking_divisor_coordinates: Vec<i64>,
    inferred_other_normal_degree: i64,
    toric_gv1_formula_value: Option<i64>,
    all_non_origin_relation_points_are_two_face: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct CmsGeneralDivisorIntersectionCheck {
    shrinking_divisor_index: usize,
    has_rational_divisor_solution: bool,
    linear_system: Option<CmsGeneralDivisorLinearSystemDiagnostic>,
    solution_basis_support_len: Option<usize>,
    solution_basis_nonzero: Option<Vec<(usize, String)>>,
    solution_ambient_basis_nonzero: Option<Vec<(usize, String)>>,
    solution_is_integral: Option<bool>,
    computed_other_normal_degree: Option<String>,
    matches_inferred_other_normal_degree: Option<bool>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct CmsGeneralDivisorLinearSystemDiagnostic {
    row_count: usize,
    column_count: usize,
    rank: usize,
    augmented_rank: usize,
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

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct OriginCircuitWitnessSample {
    first_facet_exclusive_point: usize,
    second_facet_exclusive_point: usize,
    shared_two_simplex: Vec<usize>,
    shared_two_simplex_points: Vec<OriginCircuitRelationPointSample>,
    shared_two_simplex_star_simplices: Vec<Vec<usize>>,
    shared_two_simplex_star_extra_point_samples: Vec<Vec<OriginCircuitRelationPointSample>>,
    first_facet: Vec<usize>,
    second_facet: Vec<usize>,
    first_facet_size: usize,
    second_facet_size: usize,
    sparse_relation: Vec<(usize, i64)>,
    relation_points: Vec<OriginCircuitRelationPointSample>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
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
    small_curve_pruning: CurvePruningStrategy,
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
    diagnose_corrected_chamber_gv: bool,
    diagnose_corrected_chamber_provided_generators_gv: bool,
    diagnose_corrected_chamber_ray_gv: bool,
    diagnose_corrected_chamber_lp_face_gv: bool,
    dump_corrected_chamber_gv_context_path: Option<String>,
    dump_corrected_chamber_secondary_certificate_path: Option<String>,
    dump_corrected_chamber_2face_secondary_certificate_path: Option<String>,
    dump_corrected_chamber_expanded_secondary_fan_certificate_path: Option<String>,
    dump_corrected_chamber_face_triangulation_choice_summary_path: Option<String>,
    face_triangulation_max_exact_points: usize,
    face_triangulation_samples_per_large_face: usize,
    face_triangulation_max_sampling_attempts_per_face: usize,
    face_triangulation_seed: u64,
    diagnose_chamber_updated_kklt: bool,
    diagnose_chamber_updated_kklt_iterations: usize,
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
    let small_curve_cutoff = parse_arg_value::<f64>("--small-curve-cutoff").unwrap_or(1.0);
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
    let dump_corrected_chamber_gv_context_path =
        parse_arg_value::<String>("--dump-corrected-chamber-gv-context");
    let dump_corrected_chamber_secondary_certificate_path =
        parse_arg_value::<String>("--dump-corrected-chamber-secondary-certificate");
    let dump_corrected_chamber_2face_secondary_certificate_path =
        parse_arg_value::<String>("--dump-corrected-chamber-2face-secondary-certificate");
    let dump_corrected_chamber_expanded_secondary_fan_certificate_path =
        parse_arg_value::<String>("--dump-corrected-chamber-expanded-secondary-fan-certificate");
    let dump_corrected_chamber_face_triangulation_choice_summary_path =
        parse_arg_value::<String>("--dump-corrected-chamber-face-triangulation-choice-summary");
    let face_triangulation_max_exact_points =
        parse_arg_value::<usize>("--face-triangulation-max-exact-points")
            .unwrap_or(DEFAULT_FACE_TRIANGULATION_MAX_EXACT_POINTS);
    let face_triangulation_samples_per_large_face =
        parse_arg_value::<usize>("--face-triangulation-samples-per-large-face")
            .unwrap_or(DEFAULT_FACE_TRIANGULATION_SAMPLES_PER_LARGE_FACE);
    let face_triangulation_max_sampling_attempts_per_face =
        parse_arg_value::<usize>("--face-triangulation-max-sampling-attempts-per-face")
            .unwrap_or(DEFAULT_FACE_TRIANGULATION_MAX_SAMPLING_ATTEMPTS_PER_FACE);
    let face_triangulation_seed = parse_arg_value::<u64>("--face-triangulation-seed")
        .unwrap_or(DEFAULT_FACE_TRIANGULATION_SEED);
    let diagnose_chamber_updated_kklt = parse_flag("--diagnose-chamber-updated-kklt");
    let diagnose_chamber_updated_kklt_iterations =
        parse_arg_value::<usize>("--diagnose-chamber-updated-kklt-iterations").unwrap_or(6);
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
        diagnose_corrected_chamber_gv,
        diagnose_corrected_chamber_provided_generators_gv,
        diagnose_corrected_chamber_ray_gv,
        diagnose_corrected_chamber_lp_face_gv,
        dump_corrected_chamber_gv_context_path,
        dump_corrected_chamber_secondary_certificate_path,
        dump_corrected_chamber_2face_secondary_certificate_path,
        dump_corrected_chamber_expanded_secondary_fan_certificate_path,
        dump_corrected_chamber_face_triangulation_choice_summary_path,
        face_triangulation_max_exact_points,
        face_triangulation_samples_per_large_face,
        face_triangulation_max_sampling_attempts_per_face,
        face_triangulation_seed,
        diagnose_chamber_updated_kklt,
        diagnose_chamber_updated_kklt_iterations,
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
    manifest_dir: &PathBuf,
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
    let (_dual_heights, dual_triangulation) =
        cyrus_core::compute_frst_heights(&dual_points_vec, dual_origin_idx)
            .expect("Failed to compute dual FRST");
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

fn sparse_intersection_entries(kappa: &cyrus_core::Intersection) -> Vec<SparseIntersectionEntry> {
    kappa
        .iter()
        .map(|(&(i, j, k), value)| SparseIntersectionEntry {
            indices: [i, j, k],
            value: value.get().to_string(),
        })
        .collect()
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
    let general_gvs = cyrus_core::compute_gv_invariants(
        &rays,
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
    triangulation_points: &[Point],
    triangulation: Option<&Triangulation>,
    basis_rays: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    kahler: Option<&[F64<Finite>]>,
    gamma: Option<&[I64<Finite>]>,
    origin_idx: usize,
    origin_circuits_by_class: &HashMap<Vec<i64>, cyrus_core::OriginCircuitCurveDiagnostic>,
    cms_intersection_checks_by_class: &HashMap<Vec<i64>, Vec<CmsGeneralDivisorIntersectionCheck>>,
    check_real_cone_decomposition: bool,
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
            origin_circuit_affine_rank_counts: BTreeMap::new(),
            branch_status_counts: BTreeMap::new(),
            branch_bucket_counts: BTreeMap::new(),
            real_cone_decomposition_exact_kind_counts: BTreeMap::new(),
            sample: Vec::new(),
        });
    }
    if let Some(kahler) = kahler
        && kahler.len() != basis.len()
    {
        return Err(format!(
            "Kähler vector length {} does not match basis length {}",
            kahler.len(),
            basis.len()
        ));
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
    let mut origin_circuit_affine_rank_counts = BTreeMap::new();
    let mut branch_status_counts = BTreeMap::new();
    let mut branch_bucket_counts = BTreeMap::new();
    let mut real_cone_decomposition_exact_kind_counts = BTreeMap::new();
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
        let basis_class =
            project_ambient_curve_to_basis(ambient_class, basis).map_err(|e| e.to_string())?;
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
        let real_cone_decomposition = if check_real_cone_decomposition {
            real_cone_decomposition_by_other_degree_bounded_generators(
                &basis_class,
                basis_rays,
                grading,
                degree,
            )?
        } else {
            None
        };
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
        let real_cone_decomposition_active_generator_basis_nonzero =
            real_cone_decomposition.as_ref().map(|witness| {
                witness
                    .active_generator_indices
                    .iter()
                    .map(|&idx| {
                        basis_rays[idx]
                            .iter()
                            .enumerate()
                            .filter_map(|(coord, &value)| (value != 0).then_some((coord, value)))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            });
        let real_cone_decomposition_exact_coefficients = real_cone_decomposition
            .as_ref()
            .map(|witness| {
                exact_active_generator_coefficients(
                    &basis_class,
                    basis_rays,
                    &witness.active_generator_indices,
                )
            })
            .transpose()?;
        let real_cone_decomposition_exact_kind = real_cone_decomposition_exact_coefficients
            .as_ref()
            .map(|coefficients| classify_exact_active_generator_coefficients(coefficients));
        if let Some(kind) = real_cone_decomposition_exact_kind {
            *real_cone_decomposition_exact_kind_counts
                .entry(kind.to_string())
                .or_insert(0) += 1;
        }
        let real_cone_decomposition_exact_coefficients = real_cone_decomposition_exact_coefficients
            .map(|coefficients| coefficients.iter().map(ToString::to_string).collect());
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
            .map(|witness| {
                origin_circuit_witness_sample(witness, triangulation_points, triangulation)
            });
        let origin_circuit_witnesses = origin_circuit_diagnostic
            .map(|diagnostic| {
                diagnostic
                    .witnesses
                    .iter()
                    .map(|witness| {
                        origin_circuit_witness_sample(witness, triangulation_points, triangulation)
                    })
                    .collect::<Vec<_>>()
            })
            .filter(|witnesses| !witnesses.is_empty());
        let origin_circuit_affine_support = if origin_circuit_diagnostic.is_some() {
            let affine_diagnostic =
                cyrus_core::diagnose_affine_toric_circuit(ambient_class, triangulation_points)
                    .map_err(|e| {
                        format!(
                            "failed to diagnose origin-circuit affine support for {:?}: {e}",
                            sparse_i64(ambient_class)
                        )
                    })?
                    .ok_or_else(|| {
                        format!(
                            "origin-circuit class is not an affine toric circuit: {:?}",
                            sparse_i64(ambient_class)
                        )
                    })?;
            *origin_circuit_affine_rank_counts
                .entry(affine_diagnostic.affine_rank)
                .or_insert(0) += 1;
            Some(origin_circuit_affine_support_sample(&affine_diagnostic))
        } else {
            None
        };
        let cms_general_divisor_shape_candidates = origin_circuit_diagnostic
            .map(cms_general_divisor_shape_candidates)
            .filter(|candidates| !candidates.is_empty());
        let cms_general_divisor_intersection_checks = cms_intersection_checks_by_class
            .get(ambient_class)
            .filter(|checks| !checks.is_empty())
            .cloned();
        let branch_diagnostic = match (kahler, gamma) {
            (Some(kahler), Some(gamma)) => {
                let diagnostic = missing_gv_branch_diagnostic(ambient_class, basis, kahler, gamma)?;
                let status_label = diagnostic.dilog_status.clone();
                let bucket_label = format!(
                    "parity_mod2={};qdot_bin={};status={}",
                    diagnostic.parity_mod2, diagnostic.q_dot_bucket, diagnostic.dilog_status
                );
                *branch_status_counts.entry(status_label).or_insert(0) += 1;
                *branch_bucket_counts.entry(bucket_label).or_insert(0) += 1;
                Some(diagnostic)
            }
            _ => None,
        };
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
                origin_circuit_witnesses,
                origin_circuit_affine_support,
                cms_general_divisor_shape_candidates,
                cms_general_divisor_intersection_checks,
                branch_diagnostic,
                real_cone_decomposable_by_other_generators: real_cone_decomposable,
                real_cone_decomposition_active_generators,
                real_cone_decomposition_active_generator_basis_nonzero,
                real_cone_decomposition_exact_coefficients,
                real_cone_decomposition_exact_kind,
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
        origin_circuit_affine_rank_counts,
        branch_status_counts,
        branch_bucket_counts,
        real_cone_decomposition_exact_kind_counts,
        sample,
    })
}

fn missing_gv_branch_diagnostic(
    ambient_class: &[i64],
    basis: &[usize],
    kahler: &[F64<Finite>],
    gamma: &[I64<Finite>],
) -> Result<MissingGvBranchDiagnostic, String> {
    if basis.len() != kahler.len() {
        return Err(format!(
            "Kähler vector length {} does not match basis length {}",
            kahler.len(),
            basis.len()
        ));
    }
    if basis.iter().any(|&idx| idx >= ambient_class.len()) {
        return Err(format!(
            "basis index out of bounds for missing GV class of dimension {}",
            ambient_class.len()
        ));
    }
    let q_dot_t = basis
        .iter()
        .zip(kahler.iter())
        .map(|(&idx, ti)| ambient_class[idx] as f64 * ti.get())
        .sum::<f64>();
    let parity = ambient_curve_b_field_parity_diagnostic(ambient_class, basis, gamma)
        .ok_or_else(|| "failed to compute missing GV B-field parity".to_string())?;
    let parity_mod2 = parity.rem_euclid(2);
    let dilog_status = match cyrus_core::kklt::gv_dilog_from_curve_volume_checked(q_dot_t, parity) {
        Ok(_) => "real_ok".to_string(),
        Err(failure) => format!("{failure:?}"),
    };
    Ok(MissingGvBranchDiagnostic {
        q_dot_t: format!("{q_dot_t:.17}"),
        parity,
        parity_mod2,
        q_dot_bucket: gv_branch_q_dot_bucket(q_dot_t),
        dilog_status,
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
        let basis_class =
            project_ambient_curve_to_basis(ambient_class, basis).map_err(|e| e.to_string())?;
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
    triangulation_points: &[Point],
    triangulation: Option<&Triangulation>,
) -> OriginCircuitWitnessSample {
    OriginCircuitWitnessSample {
        first_facet_exclusive_point: witness.first_facet_exclusive_point,
        second_facet_exclusive_point: witness.second_facet_exclusive_point,
        shared_two_simplex: witness.shared_two_simplex.clone(),
        shared_two_simplex_points: origin_circuit_shared_two_simplex_point_samples(
            witness,
            triangulation_points,
        ),
        shared_two_simplex_star_simplices: origin_circuit_shared_two_simplex_star_simplices(
            witness,
            triangulation,
        ),
        shared_two_simplex_star_extra_point_samples:
            origin_circuit_shared_two_simplex_star_extra_point_samples(
                witness,
                triangulation_points,
                triangulation,
            ),
        first_facet: witness.first_facet.clone(),
        second_facet: witness.second_facet.clone(),
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

fn origin_circuit_shared_two_simplex_point_samples(
    witness: &cyrus_core::OriginCircuitCurveWitness,
    triangulation_points: &[Point],
) -> Vec<OriginCircuitRelationPointSample> {
    let relation_by_point = witness
        .relation_points
        .iter()
        .map(|point| (point.point_index, point))
        .collect::<HashMap<_, _>>();
    witness
        .shared_two_simplex
        .iter()
        .filter_map(|&point_index| {
            if let Some(relation_point) = relation_by_point.get(&point_index) {
                return Some(OriginCircuitRelationPointSample {
                    point_index,
                    coefficient: relation_point.coefficient,
                    coordinates: relation_point.coordinates.clone(),
                    face_dimension: relation_point.face_dimension,
                });
            }
            let coordinates = triangulation_points.get(point_index)?.coords().to_vec();
            Some(OriginCircuitRelationPointSample {
                point_index,
                coefficient: 0,
                coordinates,
                face_dimension: None,
            })
        })
        .collect()
}

fn origin_circuit_point_sample(
    point_index: usize,
    relation_by_point: &HashMap<usize, &cyrus_core::OriginCircuitRelationPoint>,
    triangulation_points: &[Point],
) -> Option<OriginCircuitRelationPointSample> {
    if let Some(relation_point) = relation_by_point.get(&point_index) {
        return Some(OriginCircuitRelationPointSample {
            point_index,
            coefficient: relation_point.coefficient,
            coordinates: relation_point.coordinates.clone(),
            face_dimension: relation_point.face_dimension,
        });
    }
    let coordinates = triangulation_points.get(point_index)?.coords().to_vec();
    Some(OriginCircuitRelationPointSample {
        point_index,
        coefficient: 0,
        coordinates,
        face_dimension: None,
    })
}

fn origin_circuit_shared_two_simplex_star_simplices(
    witness: &cyrus_core::OriginCircuitCurveWitness,
    triangulation: Option<&Triangulation>,
) -> Vec<Vec<usize>> {
    let Some(triangulation) = triangulation else {
        return Vec::new();
    };
    let shared = witness
        .shared_two_simplex
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    let mut simplices = triangulation
        .simplices()
        .iter()
        .filter(|simplex| {
            shared
                .iter()
                .all(|point_index| simplex.contains(point_index))
        })
        .map(|simplex| {
            let mut simplex = simplex.clone();
            simplex.sort_unstable();
            simplex
        })
        .collect::<Vec<_>>();
    simplices.sort();
    simplices.dedup();
    simplices
}

fn origin_circuit_shared_two_simplex_star_extra_point_samples(
    witness: &cyrus_core::OriginCircuitCurveWitness,
    triangulation_points: &[Point],
    triangulation: Option<&Triangulation>,
) -> Vec<Vec<OriginCircuitRelationPointSample>> {
    let relation_by_point = witness
        .relation_points
        .iter()
        .map(|point| (point.point_index, point))
        .collect::<HashMap<_, _>>();
    let shared = witness
        .shared_two_simplex
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    origin_circuit_shared_two_simplex_star_simplices(witness, triangulation)
        .iter()
        .map(|simplex| {
            simplex
                .iter()
                .copied()
                .filter(|point_index| *point_index != 0 && !shared.contains(point_index))
                .filter_map(|point_index| {
                    origin_circuit_point_sample(
                        point_index,
                        &relation_by_point,
                        triangulation_points,
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn origin_circuit_affine_support_sample(
    diagnostic: &cyrus_core::AffineToricCircuitDiagnostic,
) -> OriginCircuitAffineSupportSample {
    let mut coefficient_counts = BTreeMap::new();
    for point in &diagnostic.relation_points {
        *coefficient_counts.entry(point.coefficient).or_insert(0) += 1;
    }
    OriginCircuitAffineSupportSample {
        affine_rank: diagnostic.affine_rank,
        coefficient_counts,
        local_charge_basis: diagnostic.local_charge_basis.clone(),
        local_coordinates: diagnostic
            .local_coordinates
            .iter()
            .map(|point| OriginCircuitLocalCoordinateSample {
                point_index: point.point_index,
                coordinates: point.coordinates.clone(),
            })
            .collect(),
        local_coordinates_2d: diagnostic.local_coordinates_2d.as_ref().map(|coordinates| {
            coordinates
                .iter()
                .map(|point| OriginCircuitLocalCoordinate2DSample {
                    point_index: point.point_index,
                    coordinates: point.coordinates,
                })
                .collect()
        }),
    }
}

fn toric_curve_gv_source_summary(source: &cyrus_core::ToricCurveGvSource) -> String {
    match source {
        cyrus_core::ToricCurveGvSource::TwoFace {
            edge,
            two_face_points,
            two_face_genus,
            edge_coefficients,
            edge_face_dimensions,
            edge_one_face_genera,
        } => format!(
            "two_face(edge={edge:?};coeffs={edge_coefficients:?};face_dims={edge_face_dimensions:?};one_face_genera={edge_one_face_genera:?};two_face_genus={two_face_genus};two_face_points={two_face_points:?})"
        ),
        cyrus_core::ToricCurveGvSource::ResolvedConifoldOriginCircuit {
            origin_index,
            witness,
        } => format!(
            "resolved_conifold_origin(origin={origin_index};shared_two_simplex={:?};relation={:?})",
            witness.shared_two_simplex, witness.sparse_relation
        ),
    }
}

fn toric_curve_gv_source_bucket_label(diagnostic: &cyrus_core::ToricCurveGvDiagnostic) -> String {
    let mut labels = diagnostic
        .sources
        .iter()
        .map(|source| match source {
            cyrus_core::ToricCurveGvSource::TwoFace {
                two_face_genus,
                edge_coefficients,
                edge_face_dimensions,
                edge_one_face_genera,
                ..
            } => format!(
                "two_face:g={two_face_genus}:coeffs={edge_coefficients:?}:edge_dims={edge_face_dimensions:?}:edge_one_face_genera={edge_one_face_genera:?}"
            ),
            cyrus_core::ToricCurveGvSource::ResolvedConifoldOriginCircuit { .. } => {
                "resolved_conifold_origin".to_string()
            }
        })
        .collect::<Vec<_>>();
    labels.sort();
    labels.dedup();
    if labels.len() == 1 {
        labels.pop().expect("one label is present")
    } else {
        format!("multi:{}", labels.join("|"))
    }
}

fn toric_curve_gv_diagnostic_summary(diagnostic: &cyrus_core::ToricCurveGvDiagnostic) -> String {
    let mut sources = diagnostic
        .sources
        .iter()
        .map(toric_curve_gv_source_summary)
        .collect::<Vec<_>>();
    sources.sort();
    format!(
        "gv={};source_count={};sources=[{}]",
        diagnostic.gv,
        diagnostic.sources.len(),
        sources.join("|")
    )
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
    let linear_system = Some(divisor_representation_linear_system_diagnostic(
        ambient_class,
        i,
        kappa_full,
        basis,
    )?);
    let Some(solution) = divisor_representation_solution(ambient_class, i, kappa_full, basis)?
    else {
        return Ok(CmsGeneralDivisorIntersectionCheck {
            shrinking_divisor_index: i,
            has_rational_divisor_solution: false,
            linear_system,
            solution_basis_support_len: None,
            solution_basis_nonzero: None,
            solution_ambient_basis_nonzero: None,
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
    let solution_basis_nonzero = sparse_rational_strings_by_index(&solution.coefficients);
    let solution_ambient_basis_nonzero = solution_basis_nonzero
        .iter()
        .map(|(basis_pos, value)| {
            basis
                .get(*basis_pos)
                .map(|&ambient_idx| (ambient_idx, value.clone()))
                .ok_or_else(|| {
                    format!(
                        "CMS divisor solution basis position {basis_pos} is out of bounds for basis length {}",
                        basis.len()
                    )
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let solution_is_integral = solution.coefficients.iter().all(rational_is_integer);
    let inferred_m = malachite::Rational::from(candidate.inferred_other_normal_degree);
    Ok(CmsGeneralDivisorIntersectionCheck {
        shrinking_divisor_index: i,
        has_rational_divisor_solution: true,
        linear_system,
        solution_basis_support_len: Some(support_len),
        solution_basis_nonzero: Some(solution_basis_nonzero),
        solution_ambient_basis_nonzero: Some(solution_ambient_basis_nonzero),
        solution_is_integral: Some(solution_is_integral),
        computed_other_normal_degree: Some(solution.other_normal_degree.to_string()),
        matches_inferred_other_normal_degree: Some(solution.other_normal_degree == inferred_m),
    })
}

fn divisor_representation_linear_system_diagnostic(
    ambient_class: &[i64],
    point_idx: usize,
    kappa_full: &cyrus_core::Intersection,
    basis: &[usize],
) -> Result<CmsGeneralDivisorLinearSystemDiagnostic, String> {
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
    let augmented = matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, rhs_value)| {
            let mut augmented_row = row.clone();
            augmented_row.push(rhs_value.clone());
            augmented_row
        })
        .collect::<Vec<_>>();

    Ok(CmsGeneralDivisorLinearSystemDiagnostic {
        row_count: matrix.len(),
        column_count: basis.len(),
        rank: cyrus_core::integer_math::matrix_rank(&matrix),
        augmented_rank: cyrus_core::integer_math::matrix_rank(&augmented),
    })
}

fn sparse_rational_strings_by_index(values: &[malachite::Rational]) -> Vec<(usize, String)> {
    values
        .iter()
        .enumerate()
        .filter(|(_, value)| **value != 0)
        .map(|(idx, value)| (idx, value.to_string()))
        .collect()
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

fn exact_active_generator_coefficients(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    active_generator_indices: &[usize],
) -> Result<Vec<malachite::Rational>, String> {
    if active_generator_indices.is_empty() {
        return Err("exact active-generator decomposition has no active generators".to_string());
    }
    let dim = target.len();
    let mut matrix = vec![Vec::with_capacity(active_generator_indices.len()); dim];
    for &ray_idx in active_generator_indices {
        let Some(ray) = basis_rays.get(ray_idx) else {
            return Err(format!(
                "active generator index {ray_idx} is out of bounds for {} basis rays",
                basis_rays.len()
            ));
        };
        if ray.len() != dim {
            return Err(format!(
                "active generator dimension {} does not match target dimension {dim}",
                ray.len()
            ));
        }
        for (coord, &value) in ray.iter().enumerate() {
            matrix[coord].push(malachite::Rational::from(value));
        }
    }
    let rhs = target
        .iter()
        .map(|&value| malachite::Rational::from(value))
        .collect::<Vec<_>>();
    let Some(coefficients) = solve_rational_linear_system(&matrix, &rhs)? else {
        return Err("LP active generators have no exact rational decomposition".to_string());
    };
    verify_exact_active_generator_coefficients(
        target,
        basis_rays,
        active_generator_indices,
        &coefficients,
    )?;
    Ok(coefficients)
}

fn verify_exact_active_generator_coefficients(
    target: &[i64],
    basis_rays: &[Vec<i64>],
    active_generator_indices: &[usize],
    coefficients: &[malachite::Rational],
) -> Result<(), String> {
    if active_generator_indices.len() != coefficients.len() {
        return Err(format!(
            "active generator count {} does not match coefficient count {}",
            active_generator_indices.len(),
            coefficients.len()
        ));
    }
    for coord in 0..target.len() {
        let mut reconstructed = malachite::Rational::from(0);
        for (&ray_idx, coefficient) in active_generator_indices.iter().zip(coefficients.iter()) {
            reconstructed += coefficient * malachite::Rational::from(basis_rays[ray_idx][coord]);
        }
        if reconstructed != malachite::Rational::from(target[coord]) {
            return Err(format!(
                "exact active-generator decomposition failed at coordinate {coord}: reconstructed {reconstructed}, target {}",
                target[coord]
            ));
        }
    }
    Ok(())
}

fn classify_exact_active_generator_coefficients(
    coefficients: &[malachite::Rational],
) -> &'static str {
    let zero = malachite::Rational::from(0);
    if coefficients.iter().all(|value| value >= &zero) {
        if coefficients.iter().all(rational_is_integer) {
            "integer_semigroup"
        } else {
            "rational_cone"
        }
    } else {
        "signed_rational"
    }
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

fn compute_missing_one_dimensional_ray_gvs(
    ambient_classes: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
) -> Result<Vec<(Vec<i64>, malachite::Integer, i128)>, String> {
    let mut out = Vec::with_capacity(ambient_classes.len());
    for ambient_class in ambient_classes {
        let basis_class =
            project_ambient_curve_to_basis(ambient_class, basis).map_err(|e| e.to_string())?;
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
        let series = cyrus_core::compute_ambient_one_dimensional_ray_gv_series(
            ambient_class,
            basis,
            grading,
            q_matrix,
            intnums,
            1,
        )
        .map_err(|e| {
            format!(
                "failed one-dimensional ray GV computation for target degree {degree}, ambient_nonzero={:?}: {e}",
                sparse_i64(ambient_class)
            )
        })?;
        let gv = series
            .values
            .into_iter()
            .next()
            .expect("max_multiple=1 returns one GV value");
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
        let basis_class =
            project_ambient_curve_to_basis(ambient_class, basis).map_err(|e| e.to_string())?;
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

fn corrected_chamber_top_toric_local_gv_requested() -> bool {
    std::env::var("CYRUS_CORRECTED_CHAMBER_TOP_TORIC_LOCAL_GV")
        .map(|value| value != "0")
        .unwrap_or(false)
}

fn corrected_chamber_gv_weight_lp_requested() -> bool {
    std::env::var("CYRUS_CORRECTED_CHAMBER_GV_WEIGHT_LP")
        .map(|value| value != "0")
        .unwrap_or(false)
}

fn corrected_chamber_top_toric_local_gv_limit() -> usize {
    std::env::var("CYRUS_CORRECTED_CHAMBER_TOP_TORIC_LOCAL_GV_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(8)
}

fn corrected_chamber_top_toric_local_gv_witness_limit() -> usize {
    std::env::var("CYRUS_CORRECTED_CHAMBER_TOP_TORIC_LOCAL_GV_WITNESS_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(4)
}

fn corrected_chamber_top_toric_local_gv_delta_limit() -> usize {
    std::env::var("CYRUS_CORRECTED_CHAMBER_TOP_TORIC_LOCAL_GV_DELTA_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(4)
}

fn corrected_chamber_top_toric_local_gv_contribution_limit() -> usize {
    std::env::var("CYRUS_CORRECTED_CHAMBER_TOP_TORIC_LOCAL_GV_CONTRIBUTION_LIMIT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(6)
}

fn report_corrected_chamber_top_toric_local_gv_diagnostics(
    tri: &Triangulation,
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    targets: &[TopToricLocalGvTarget],
) {
    if targets.is_empty() {
        return;
    }
    if cfg!(panic = "abort") {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic skipped: panic=abort build cannot contain cygv truncation panics"
        );
        return;
    }

    let limit = corrected_chamber_top_toric_local_gv_limit();
    let witness_limit = corrected_chamber_top_toric_local_gv_witness_limit();
    let ambient_rays = match compute_mori_cone_cap_rays(
        tri,
        &geom.triangulation_points,
        &geom.polytope,
        false,
        false,
        None,
    ) {
        Ok(rays) => rays,
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic failed to compute Mori-cap rays: {e}"
            );
            return;
        }
    };
    let gv_basis_data = match vector_gv_basis_data(
        &ambient_rays,
        &intersection.linrels,
        &intersection.basis,
        "checkpoint-t corrected-chamber top-toric local GV diagnostic",
    ) {
        Ok(data) => data,
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic failed to build GV basis data: {e}"
            );
            return;
        }
    };
    let basis_rays = &gv_basis_data.mori_rays;
    let Some(grading) = compute_grading_vector(basis_rays) else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic failed to compute grading vector"
        );
        return;
    };
    match non_positive_basis_generator_degrees(basis_rays, &grading) {
        Ok((non_positive_count, Some((idx, degree, ray)))) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic skipped: grading has {non_positive_count}/{} non-positive Mori generator degrees; first index={idx} degree={degree} ray={ray:?}",
                basis_rays.len()
            );
            return;
        }
        Ok((_non_positive_count, None)) => {}
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic failed to inspect grading degrees: {e}"
            );
            return;
        }
    }
    let corrected_kappa_basis = match chamber_intersection_in_basis(
        tri,
        &geom.triangulation_points,
        &intersection.basis,
    ) {
        Ok(kappa) => kappa,
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV diagnostic failed to compute basis intersections: {e}"
            );
            return;
        }
    };

    let mut seen = HashSet::new();
    let mut attempted = 0usize;
    for target in targets {
        if attempted >= limit {
            break;
        }
        if !seen.insert(target.class.clone()) {
            continue;
        }
        attempted += 1;
        report_single_top_toric_local_gv_diagnostic(
            target,
            basis_rays,
            &intersection.basis,
            &grading,
            &gv_basis_data.q_matrix,
            &corrected_kappa_basis,
            witness_limit,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn report_single_top_toric_local_gv_diagnostic(
    target: &TopToricLocalGvTarget,
    basis_rays: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &cyrus_core::Intersection,
    witness_limit: usize,
) {
    let basis_class = match project_ambient_curve_to_basis(&target.class, basis) {
        Ok(class) => class,
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} toric_gv={} error=projection_failed:{e} class={:?}",
                target.kklt_idx,
                target.point_idx,
                target.toric_gv,
                sparse_i64(&target.class)
            );
            return;
        }
    };
    let degree = basis_class
        .iter()
        .zip(grading.iter())
        .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
        .sum::<i128>();
    if degree <= 0 {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} toric_gv={} degree={} error=non_positive_degree class={:?}",
            target.kklt_idx,
            target.point_idx,
            target.toric_gv,
            degree,
            sparse_i64(&target.class)
        );
        return;
    }
    let Ok(max_deg) = u32::try_from(degree) else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} toric_gv={} degree={} error=degree_too_large class={:?}",
            target.kklt_idx,
            target.point_idx,
            target.toric_gv,
            degree,
            sparse_i64(&target.class)
        );
        return;
    };

    let witnesses = match real_cone_decomposition_witnesses_by_other_degree_bounded_generators(
        &basis_class,
        basis_rays,
        grading,
        degree,
        witness_limit,
    ) {
        Ok(witnesses) => witnesses,
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} toric_gv={} degree={} error=witness_failed:{e} class={:?}",
                target.kklt_idx,
                target.point_idx,
                target.toric_gv,
                degree,
                sparse_i64(&target.class)
            );
            return;
        }
    };
    for witness in &witnesses {
        match compute_lp_witness_face_attempt(
            witness,
            &basis_class,
            basis_rays,
            grading,
            q_matrix,
            intnums,
            degree,
            max_deg,
            &target.class,
        ) {
            Ok(attempt) => {
                if let Some(local_gv) = attempt.gv.as_ref() {
                    eprintln!(
                        "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=lp_witness degree={} toric_gv={} local_gv={} matches={} generator_count={} active_generator_count={} span_generator_count={} used_span_expansion={} used_lattice_saturation={} used_integer_decomposition={} used_decomposition_diamond={} used_decomposition_closure={} integer_decomposition_terms={:?} lattice_elements={:?} error={:?} class={:?}",
                        target.kklt_idx,
                        target.point_idx,
                        degree,
                        target.toric_gv,
                        local_gv,
                        local_gv == &target.toric_gv,
                        attempt.generator_count,
                        attempt.active_generator_count,
                        attempt.span_generator_count,
                        attempt.used_span_expansion,
                        attempt.used_lattice_saturation,
                        attempt.used_integer_decomposition,
                        attempt.used_decomposition_diamond,
                        attempt.used_decomposition_closure,
                        attempt.integer_decomposition_term_count,
                        attempt.lattice_semigroup_element_count,
                        attempt.error,
                        sparse_i64(&target.class)
                    );
                    return;
                }
                if let Some(error) = attempt.error {
                    eprintln!(
                        "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=lp_witness degree={} toric_gv={} local_gv=None error={} class={:?}",
                        target.kklt_idx,
                        target.point_idx,
                        degree,
                        target.toric_gv,
                        error,
                        sparse_i64(&target.class)
                    );
                }
            }
            Err(e) => {
                eprintln!(
                    "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=lp_witness degree={} toric_gv={} local_gv=None error={} class={:?}",
                    target.kklt_idx,
                    target.point_idx,
                    degree,
                    target.toric_gv,
                    e,
                    sparse_i64(&target.class)
                );
            }
        }
    }

    match one_dimensional_ray_gv_targets(
        std::slice::from_ref(&target.class),
        basis_rays,
        basis,
        grading,
    ) {
        Ok(ray_targets) if !ray_targets.candidates.is_empty() => {
            match compute_missing_one_dimensional_ray_gvs(
                &ray_targets.candidates,
                basis,
                grading,
                q_matrix,
                intnums,
            ) {
                Ok(ray_gvs) => {
                    for (ambient_class, local_gv, ray_degree) in ray_gvs {
                        eprintln!(
                            "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=one_dimensional_ray degree={} toric_gv={} local_gv={} matches={} class={:?}",
                            target.kklt_idx,
                            target.point_idx,
                            ray_degree,
                            target.toric_gv,
                            local_gv,
                            local_gv == target.toric_gv,
                            sparse_i64(&ambient_class)
                        );
                    }
                }
                Err(e) => {
                    eprintln!(
                        "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=one_dimensional_ray degree={} toric_gv={} local_gv=None error={} class={:?}",
                        target.kklt_idx,
                        target.point_idx,
                        degree,
                        target.toric_gv,
                        e,
                        sparse_i64(&target.class)
                    );
                }
            }
        }
        Ok(ray_targets) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=none degree={} toric_gv={} witness_count={} skipped_non_generators={} skipped_decomposable_generators={} class={:?}",
                target.kklt_idx,
                target.point_idx,
                degree,
                target.toric_gv,
                witnesses.len(),
                ray_targets.skipped_non_generators,
                ray_targets.skipped_decomposable_generators,
                sparse_i64(&target.class)
            );
        }
        Err(e) => {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber top-toric local GV kklt_idx={} point_idx={} method=none degree={} toric_gv={} witness_count={} error={} class={:?}",
                target.kklt_idx,
                target.point_idx,
                degree,
                target.toric_gv,
                witnesses.len(),
                e,
                sparse_i64(&target.class)
            );
        }
    }
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

    if let Some(certificate) =
        certify_supporting_mori_face_by_exact_kernel(face_generators, basis_rays)
            .map_err(|err| err.to_string())?
    {
        return Ok(Some(SupportingFaceCertificate {
            zero_generator_count: certificate.zero_generator_count,
            positive_generator_count: certificate.positive_generator_count,
            normal: certificate.normal,
        }));
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

    if let Some(lp_normal) = solve_supporting_face_normal_aggregate_lp(face_generators, basis_rays)?
        && let Some(certificate) =
            integer_supporting_face_certificate_from_lp(&lp_normal, face_generators, basis_rays)?
    {
        return Ok(Some(certificate));
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

fn solve_supporting_face_normal_aggregate_lp(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
) -> Result<Option<Vec<f64>>, String> {
    let max_rounds = std::env::var("CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_CUTTING_ROUNDS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE_CUTTING_ROUNDS);
    let aggregate = aggregate_mori_ray_coefficients(basis_rays)?;
    let mut enforced_ray_indices = Vec::new();
    let mut enforced_ray_set = HashSet::new();
    for _ in 0..max_rounds {
        let Some(normal) = solve_supporting_face_normal_aggregate_lp_with_enforced_rays(
            face_generators,
            basis_rays,
            &aggregate,
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

fn aggregate_mori_ray_coefficients(basis_rays: &[Vec<i64>]) -> Result<Vec<i128>, String> {
    let Some(first_ray) = basis_rays.first() else {
        return Err("supporting-face aggregate normal LP requires Mori generators".to_string());
    };
    let dim = first_ray.len();
    let mut aggregate = vec![0_i128; dim];
    for ray in basis_rays {
        if ray.len() != dim {
            return Err(
                "supporting-face aggregate normal LP ray dimensions are inconsistent".to_string(),
            );
        }
        for (slot, &coefficient) in aggregate.iter_mut().zip(ray) {
            *slot = (*slot)
                .checked_add(i128::from(coefficient))
                .ok_or_else(|| {
                    "supporting-face aggregate normal LP coefficient overflowed i128".to_string()
                })?;
        }
    }
    Ok(aggregate)
}

fn solve_supporting_face_normal_aggregate_lp_with_enforced_rays(
    face_generators: &[Vec<i64>],
    basis_rays: &[Vec<i64>],
    aggregate: &[i128],
    enforced_ray_indices: &[usize],
) -> Result<Option<Vec<f64>>, String> {
    let dim = aggregate.len();
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

    let mut aggregate_expr = Expression::from(0.0);
    for (var, &coefficient) in normal_vars.iter().zip(aggregate) {
        if coefficient != 0 {
            aggregate_expr.add_mul(coefficient as f64, *var);
        }
    }
    model = model.with(aggregate_expr.geq(1.0));

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
            return Err(format!("supporting-face aggregate normal LP failed: {err}"));
        }
    };
    let normal = normal_vars
        .iter()
        .map(|var| solution.value(*var))
        .collect::<Vec<_>>();
    if normal.iter().all(|value| value.is_finite()) {
        Ok(Some(normal))
    } else {
        Err("supporting-face aggregate normal LP returned a non-finite value".to_string())
    }
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
        if let Some(certificate) =
            check_supporting_mori_face_normal(&normal, face_generators, basis_rays)
                .map_err(|err| err.to_string())?
        {
            return Ok(Some(SupportingFaceCertificate {
                zero_generator_count: certificate.zero_generator_count,
                positive_generator_count: certificate.positive_generator_count,
                normal: certificate.normal,
            }));
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
    let candidates = candidates
        .into_iter()
        .map(|ray| {
            let degree = ray
                .iter()
                .zip(grading.iter())
                .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
                .sum::<i128>();
            DecompositionCandidate { ray, degree }
        })
        .collect::<Vec<_>>();
    let mut partial = Vec::new();
    let mut partial_sum = vec![0i64; target.len()];
    for term_count in 2..=max_terms {
        search_generator_decompositions(
            target,
            &candidates,
            &candidate_set,
            target_degree,
            term_count,
            0,
            0,
            &mut partial,
            &mut partial_sum,
            &mut seen_decompositions,
            &mut decompositions,
            max_witnesses,
        )?;
        if decompositions.len() >= max_witnesses {
            return Ok(Some(decompositions));
        }
    }
    Ok((!decompositions.is_empty()).then_some(decompositions))
}

struct DecompositionCandidate {
    ray: Vec<i64>,
    degree: i128,
}

#[allow(clippy::too_many_arguments)]
fn search_generator_decompositions(
    target: &[i64],
    candidates: &[DecompositionCandidate],
    candidate_set: &HashSet<Vec<i64>>,
    target_degree: i128,
    desired_terms: usize,
    start_index: usize,
    partial_degree: i128,
    partial: &mut Vec<Vec<i64>>,
    partial_sum: &mut [i64],
    seen_decompositions: &mut HashSet<Vec<Vec<i64>>>,
    decompositions: &mut Vec<Vec<Vec<i64>>>,
    max_witnesses: usize,
) -> Result<(), String> {
    if decompositions.len() >= max_witnesses {
        return Ok(());
    }
    if partial.len() + 1 == desired_terms {
        let remaining_degree = target_degree - partial_degree;
        if remaining_degree <= 0 {
            return Ok(());
        }
        let mut remainder = Vec::with_capacity(target.len());
        for (&target_value, &partial_value) in target.iter().zip(partial_sum.iter()) {
            remainder.push(target_value.checked_sub(partial_value).ok_or_else(|| {
                "generator decomposition remainder coordinate overflowed i64".to_string()
            })?);
        }
        if !candidate_set.contains(&remainder) {
            return Ok(());
        }
        let Some(candidate) = candidates
            .iter()
            .skip(start_index)
            .find(|candidate| candidate.degree == remaining_degree && candidate.ray == remainder)
        else {
            return Ok(());
        };
        partial.push(candidate.ray.clone());
        let mut decomposition = partial.clone();
        decomposition.sort();
        if seen_decompositions.insert(decomposition.clone()) {
            decompositions.push(decomposition);
        }
        partial.pop();
        return Ok(());
    }

    for (idx, candidate) in candidates.iter().enumerate().skip(start_index) {
        let next_degree = partial_degree + candidate.degree;
        if next_degree >= target_degree {
            continue;
        }
        for (slot, &value) in partial_sum.iter_mut().zip(candidate.ray.iter()) {
            *slot = slot.checked_add(value).ok_or_else(|| {
                "generator decomposition partial sum coordinate overflowed i64".to_string()
            })?;
        }
        partial.push(candidate.ray.clone());
        search_generator_decompositions(
            target,
            candidates,
            candidate_set,
            target_degree,
            desired_terms,
            idx,
            next_degree,
            partial,
            partial_sum,
            seen_decompositions,
            decompositions,
            max_witnesses,
        )?;
        partial.pop();
        for (slot, &value) in partial_sum.iter_mut().zip(candidate.ray.iter()) {
            *slot = slot.checked_sub(value).ok_or_else(|| {
                "generator decomposition partial sum coordinate underflowed i64".to_string()
            })?;
        }
        if decompositions.len() >= max_witnesses {
            return Ok(());
        }
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

fn dense_i64_from_sparse(
    sparse: &[(usize, i64)],
    dim: usize,
    label: &str,
) -> Result<Vec<i64>, String> {
    let mut dense = vec![0; dim];
    let mut seen = HashSet::new();
    for &(idx, value) in sparse {
        if idx >= dim {
            return Err(format!(
                "{label} sparse coordinate {idx} is out of bounds for dimension {dim}"
            ));
        }
        if value == 0 {
            return Err(format!("{label} sparse coordinate {idx} has zero value"));
        }
        if !seen.insert(idx) {
            return Err(format!("{label} sparse coordinate {idx} is duplicated"));
        }
        dense[idx] = value;
    }
    Ok(dense)
}

fn degree_bounded_mori_ray_context_samples(
    ambient_rays: &[Vec<i64>],
    basis: &[usize],
    grading: &[i64],
    max_degree: i128,
) -> Result<Vec<DegreeBoundedMoriRayContextSample>, String> {
    let mut samples = Vec::new();
    for ambient_ray in ambient_rays {
        let basis_ray =
            project_ambient_curve_to_basis(ambient_ray, basis).map_err(|e| e.to_string())?;
        if basis_ray.len() != grading.len() {
            return Err(format!(
                "projected Mori ray dimension {} does not match grading dimension {}",
                basis_ray.len(),
                grading.len()
            ));
        }
        let degree = basis_ray
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= max_degree {
            samples.push(DegreeBoundedMoriRayContextSample {
                degree,
                ambient_nonzero: sparse_i64(ambient_ray),
                basis_nonzero: sparse_i64(&basis_ray),
            });
        }
    }
    Ok(samples)
}

fn uncovered_degree_bounded_source_ray_classes(
    ambient_rays: &[Vec<i64>],
    covered_gvs: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    grading: &[i64],
    max_degree: i128,
) -> Result<Vec<Vec<i64>>, String> {
    let covered_classes = covered_gvs
        .iter()
        .map(|(class, _)| class.clone())
        .collect::<HashSet<_>>();
    let mut out = Vec::new();
    for ambient_ray in ambient_rays {
        let degree = ambient_curve_grading_degree(ambient_ray, basis, grading)?;
        if degree <= max_degree && !covered_classes.contains(ambient_ray) {
            out.push(ambient_ray.clone());
        }
    }
    Ok(out)
}

fn active_noncovered_dependency_source_ray_classes(
    target_stats: &MissingGvTargetStats,
    degree_bounded_ray_context: &[DegreeBoundedMoriRayContextSample],
    covered_toric_context: &[CoveredToricGvContextSample],
    ambient_dim: usize,
) -> Result<Vec<Vec<i64>>, String> {
    let covered_basis_supports = covered_toric_context
        .iter()
        .map(|sample| sample.basis_nonzero.clone())
        .collect::<HashSet<_>>();
    let missing_target_basis_supports = target_stats
        .sample
        .iter()
        .map(|sample| sample.basis_nonzero.clone())
        .collect::<HashSet<_>>();
    let mut out = HashSet::new();
    for (target_idx, target) in target_stats.sample.iter().enumerate() {
        let Some(active_generators) = target
            .real_cone_decomposition_active_generator_basis_nonzero
            .as_ref()
        else {
            continue;
        };
        for active_generator in active_generators {
            if covered_basis_supports.contains(active_generator)
                || missing_target_basis_supports.contains(active_generator)
            {
                continue;
            }
            let Some(ray_context) = degree_bounded_ray_context
                .iter()
                .find(|sample| sample.basis_nonzero == *active_generator)
            else {
                return Err(format!(
                    "missing target {target_idx} active dependency {:?} is not in the degree-bounded source-ray context",
                    active_generator
                ));
            };
            out.insert(dense_i64_from_sparse(
                &ray_context.ambient_nonzero,
                ambient_dim,
                "active dependency source ray",
            )?);
        }
    }
    let mut out = out.into_iter().collect::<Vec<_>>();
    out.sort();
    Ok(out)
}

fn shared_facet_unresolved_source_ray_classes(
    target_stats: &MissingGvTargetStats,
    degree_bounded_ray_context: &[DegreeBoundedMoriRayContextSample],
    covered_toric_context: &[CoveredToricGvContextSample],
    source_derived_toric_context: &[ToricGvDiagnosticContextSample],
    already_diagnosed_source_stats: &MissingGvTargetStats,
    ambient_dim: usize,
) -> Result<Vec<Vec<i64>>, String> {
    let mut known_basis_supports = covered_toric_context
        .iter()
        .map(|sample| sample.basis_nonzero.clone())
        .collect::<HashSet<_>>();
    known_basis_supports.extend(
        source_derived_toric_context
            .iter()
            .map(|sample| sample.basis_nonzero.clone()),
    );
    known_basis_supports.extend(
        target_stats
            .sample
            .iter()
            .map(|sample| sample.basis_nonzero.clone()),
    );
    known_basis_supports.extend(
        already_diagnosed_source_stats
            .sample
            .iter()
            .map(|sample| sample.basis_nonzero.clone()),
    );

    let mut out = HashSet::new();
    for (target_idx, target) in target_stats.sample.iter().enumerate() {
        for witness in missing_target_origin_circuit_witnesses(target) {
            let shared_facet = origin_circuit_witness_shared_facet_support(witness);
            for ray in degree_bounded_ray_context {
                if ray.degree <= 0 || ray.degree > target.degree {
                    continue;
                }
                if known_basis_supports.contains(&ray.basis_nonzero) {
                    continue;
                }
                if !ray
                    .ambient_nonzero
                    .iter()
                    .all(|(idx, _)| shared_facet.contains(idx))
                {
                    continue;
                }
                out.insert(dense_i64_from_sparse(
                    &ray.ambient_nonzero,
                    ambient_dim,
                    &format!("shared-facet source ray for target {target_idx}"),
                )?);
            }
        }
    }
    let mut out = out.into_iter().collect::<Vec<_>>();
    out.sort();
    Ok(out)
}

fn missing_target_origin_circuit_witnesses(
    sample: &MissingGvTargetSample,
) -> Vec<&OriginCircuitWitnessSample> {
    sample
        .origin_circuit_witnesses
        .as_ref()
        .filter(|witnesses| !witnesses.is_empty())
        .map(|witnesses| witnesses.iter().collect())
        .or_else(|| {
            sample
                .origin_circuit_first_witness
                .as_ref()
                .map(|witness| vec![witness])
        })
        .unwrap_or_default()
}

fn origin_circuit_witness_shared_facet_support(
    witness: &OriginCircuitWitnessSample,
) -> HashSet<usize> {
    let first_facet = witness.first_facet.iter().copied().collect::<HashSet<_>>();
    let second_facet = witness.second_facet.iter().copied().collect::<HashSet<_>>();
    let mut shared_facet = first_facet
        .intersection(&second_facet)
        .copied()
        .collect::<HashSet<_>>();
    shared_facet.insert(0);
    shared_facet.insert(witness.first_facet_exclusive_point);
    shared_facet.insert(witness.second_facet_exclusive_point);
    shared_facet
}

fn merge_unique_ambient_classes(mut lhs: Vec<Vec<i64>>, rhs: Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    let mut seen = lhs.iter().cloned().collect::<HashSet<_>>();
    for class in rhs {
        if seen.insert(class.clone()) {
            lhs.push(class);
        }
    }
    lhs.sort();
    lhs
}

fn covered_toric_gv_context_samples(
    covered_gvs: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    grading: &[i64],
    max_degree: i128,
) -> Result<Vec<CoveredToricGvContextSample>, String> {
    let mut samples = Vec::new();
    for (ambient_class, gv) in covered_gvs {
        let basis_class =
            project_ambient_curve_to_basis(ambient_class, basis).map_err(|e| e.to_string())?;
        if basis_class.len() != grading.len() {
            return Err(format!(
                "covered toric GV basis class dimension {} does not match grading dimension {}",
                basis_class.len(),
                grading.len()
            ));
        }
        let degree = basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        if degree <= max_degree {
            samples.push(CoveredToricGvContextSample {
                degree,
                gv: gv.to_string(),
                ambient_nonzero: sparse_i64(ambient_class),
                basis_nonzero: sparse_i64(&basis_class),
            });
        }
    }
    samples.sort_by(|lhs, rhs| {
        lhs.degree
            .cmp(&rhs.degree)
            .then_with(|| lhs.basis_nonzero.cmp(&rhs.basis_nonzero))
            .then_with(|| lhs.ambient_nonzero.cmp(&rhs.ambient_nonzero))
    });
    Ok(samples)
}

fn toric_gv_diagnostic_context_samples(
    ambient_classes: &[Vec<i64>],
    diagnostics_by_class: &HashMap<Vec<i64>, cyrus_core::ToricCurveGvDiagnostic>,
    basis: &[usize],
    grading: &[i64],
    sample_limit: usize,
) -> Result<Vec<ToricGvDiagnosticContextSample>, String> {
    let mut samples = Vec::new();
    for ambient_class in ambient_classes {
        if samples.len() >= sample_limit {
            break;
        }
        let Some(diagnostic) = diagnostics_by_class.get(ambient_class) else {
            continue;
        };
        let basis_class =
            project_ambient_curve_to_basis(ambient_class, basis).map_err(|e| e.to_string())?;
        if basis_class.len() != grading.len() {
            return Err(format!(
                "toric GV diagnostic basis class dimension {} does not match grading dimension {}",
                basis_class.len(),
                grading.len()
            ));
        }
        let degree = basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
            .sum::<i128>();
        samples.push(ToricGvDiagnosticContextSample {
            degree,
            gv: diagnostic.gv.to_string(),
            source_bucket: toric_curve_gv_source_bucket_label(diagnostic),
            source_summary: toric_curve_gv_diagnostic_summary(diagnostic),
            ambient_nonzero: sparse_i64(ambient_class),
            basis_nonzero: sparse_i64(&basis_class),
        });
    }
    samples.sort_by(|lhs, rhs| {
        lhs.degree
            .cmp(&rhs.degree)
            .then_with(|| lhs.basis_nonzero.cmp(&rhs.basis_nonzero))
            .then_with(|| lhs.ambient_nonzero.cmp(&rhs.ambient_nonzero))
    });
    Ok(samples)
}

fn ambient_curve_b_field_parity_diagnostic(
    curve: &[i64],
    basis: &[usize],
    gamma: &[I64<Finite>],
) -> Option<i128> {
    if gamma.len() == basis.len() {
        if basis.iter().any(|&idx| idx >= curve.len()) {
            return None;
        }
        return Some(
            basis
                .iter()
                .zip(gamma.iter())
                .map(|(&idx, gi)| i128::from(curve[idx]) * i128::from(gi.get()))
                .sum(),
        );
    }
    if gamma.len() == curve.len() {
        return Some(
            curve
                .iter()
                .zip(gamma.iter())
                .map(|(&qi, gi)| i128::from(qi) * i128::from(gi.get()))
                .sum(),
        );
    }
    None
}

fn shifted_ambient_gamma(gamma: &[I64<Finite>], offset: isize) -> Option<Vec<I64<Finite>>> {
    let mut shifted = vec![0_i64; gamma.len()];
    for (idx, value) in gamma.iter().enumerate() {
        let value = value.get();
        if value == 0 {
            continue;
        }
        let shifted_idx = isize::try_from(idx).ok()?.checked_add(offset)?;
        let shifted_idx = usize::try_from(shifted_idx).ok()?;
        let slot = shifted.get_mut(shifted_idx)?;
        *slot = slot.checked_add(value)?;
    }
    Some(shifted.into_iter().map(I64::<Finite>::new).collect())
}

fn gamma_with_toggled_index(gamma: &[I64<Finite>], index: usize) -> Option<Vec<I64<Finite>>> {
    let mut toggled = gamma.to_vec();
    let entry = toggled.get_mut(index)?;
    *entry = I64::<Finite>::new(entry.get().checked_add(1)?);
    Some(toggled)
}

fn project_ambient_gamma_to_curve_basis(
    curve_basis_matrix: &[Vec<malachite::Integer>],
    ambient_gamma: &[I64<Finite>],
) -> Option<Vec<I64<Finite>>> {
    let ambient_dim = curve_basis_matrix.first()?.len();
    if ambient_dim == 0
        || ambient_gamma.len() != ambient_dim
        || curve_basis_matrix
            .iter()
            .any(|row| row.len() != ambient_dim)
    {
        return None;
    }

    curve_basis_matrix
        .iter()
        .map(|row| {
            let mut acc = malachite::Integer::from(0);
            for (coefficient, gamma_entry) in row.iter().zip(ambient_gamma.iter()) {
                acc += coefficient * malachite::Integer::from(gamma_entry.get());
            }
            i64::try_from(&acc).ok().map(I64::<Finite>::new)
        })
        .collect()
}

fn count_gamma_parity_mismatches(
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    ambient_gamma: &[I64<Finite>],
    basis_gamma: &[I64<Finite>],
) -> Option<usize> {
    let mut mismatches = 0usize;
    for (curve, _invariant) in gv_invariants {
        let ambient_parity = ambient_curve_b_field_parity_diagnostic(curve, basis, ambient_gamma)?;
        let basis_parity = ambient_curve_b_field_parity_diagnostic(curve, basis, basis_gamma)?;
        if ambient_parity.rem_euclid(2) != basis_parity.rem_euclid(2) {
            mismatches += 1;
        }
    }
    Some(mismatches)
}

fn ambient_target_contribution_rows(
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    divisor_idx: usize,
    t: &[F64<Finite>],
    gamma: &[I64<Finite>],
) -> Option<Vec<AmbientTargetContributionRow>> {
    if basis.len() != t.len() {
        return None;
    }
    let mut rows = Vec::with_capacity(gv_invariants.len());
    for (curve, invariant) in gv_invariants {
        if divisor_idx >= curve.len() || basis.iter().any(|&idx| idx >= curve.len()) {
            return None;
        }
        let q_dot_t = basis
            .iter()
            .zip(t.iter())
            .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
            .sum::<f64>();
        let parity = ambient_curve_b_field_parity_diagnostic(curve, basis, gamma)?;
        let single = [(curve.clone(), invariant.clone())];
        let contribution = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &single,
            basis,
            &[divisor_idx],
            t,
            Some(gamma),
        )?
        .first()?
        .get();
        rows.push(AmbientTargetContributionRow {
            class: curve.clone(),
            gv: invariant.clone(),
            q_divisor: curve[divisor_idx],
            q_dot_t,
            parity,
            contribution,
        });
    }
    Some(rows)
}

fn gv_branch_q_dot_bucket(q_dot_t: f64) -> &'static str {
    if q_dot_t <= 0.0 {
        "nonpositive"
    } else if q_dot_t < 0.005 {
        "0_0.005"
    } else if q_dot_t < 0.01 {
        "0.005_0.01"
    } else if q_dot_t < 0.02 {
        "0.01_0.02"
    } else if q_dot_t < 0.05 {
        "0.02_0.05"
    } else if q_dot_t < 0.1 {
        "0.05_0.1"
    } else {
        "gte_0.1"
    }
}

fn f64_target_vector(values: impl IntoIterator<Item = f64>) -> Option<Vec<F64<Finite>>> {
    values.into_iter().map(F64::<Finite>::new).collect()
}

fn max_abs_and_l2(values: &[f64]) -> (f64, f64) {
    let mut max_abs = 0.0f64;
    let mut l2_sq = 0.0f64;
    for &value in values {
        max_abs = max_abs.max(value.abs());
        l2_sq += value * value;
    }
    (max_abs, l2_sq.sqrt())
}

fn report_corrected_chamber_gv_branch_buckets(
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: &[I64<Finite>],
    input_chi_reference: &[F64<Finite>],
    corrected_chi_reference: &[F64<Finite>],
    covered_target: &[F64<Finite>],
) {
    if basis.len() != t.len() || kklt_basis.len() != covered_target.len() {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV branch buckets unavailable: dimension mismatch"
        );
        return;
    }

    let mut buckets: BTreeMap<(i128, &'static str), GvBranchBucketAccumulator> = BTreeMap::new();
    for (curve, invariant) in gv_invariants {
        if kklt_basis.iter().any(|&idx| idx >= curve.len())
            || basis.iter().any(|&idx| idx >= curve.len())
        {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV branch buckets unavailable: curve dimension mismatch"
            );
            return;
        }
        let q_dot_t = basis
            .iter()
            .zip(t.iter())
            .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
            .sum::<f64>();
        let Some(parity) = ambient_curve_b_field_parity_diagnostic(curve, basis, gamma) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV branch buckets unavailable: gamma dimension mismatch"
            );
            return;
        };
        let single = [(curve.clone(), invariant.clone())];
        let Some(contribution) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &single,
            basis,
            kklt_basis,
            t,
            Some(gamma),
        ) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV branch buckets unavailable: invalid single-curve contribution q_dot_t={q_dot_t} parity_mod2={}",
                parity.rem_euclid(2)
            );
            return;
        };
        let bucket_key = (parity.rem_euclid(2), gv_branch_q_dot_bucket(q_dot_t));
        let bucket = buckets
            .entry(bucket_key)
            .or_insert_with(|| GvBranchBucketAccumulator {
                count: 0,
                q_dot_min: f64::INFINITY,
                q_dot_max: f64::NEG_INFINITY,
                correction: vec![0.0; kklt_basis.len()],
            });
        bucket.count += 1;
        bucket.q_dot_min = bucket.q_dot_min.min(q_dot_t);
        bucket.q_dot_max = bucket.q_dot_max.max(q_dot_t);
        for (accumulated, value) in bucket.correction.iter_mut().zip(contribution.iter()) {
            *accumulated += value.get();
        }
    }

    let mut reports = buckets
        .into_iter()
        .map(|((parity_mod2, q_dot_bucket), bucket)| {
            let (max_abs_correction, l2_correction) = max_abs_and_l2(&bucket.correction);
            GvBranchBucketReport {
                parity_mod2,
                q_dot_bucket,
                count: bucket.count,
                q_dot_min: bucket.q_dot_min,
                q_dot_max: bucket.q_dot_max,
                max_abs_correction,
                l2_correction,
                correction: bucket.correction,
            }
        })
        .collect::<Vec<_>>();
    reports.sort_unstable_by(|lhs, rhs| {
        rhs.l2_correction
            .total_cmp(&lhs.l2_correction)
            .then_with(|| lhs.parity_mod2.cmp(&rhs.parity_mod2))
            .then_with(|| lhs.q_dot_bucket.cmp(rhs.q_dot_bucket))
    });

    for report in reports {
        let Some(target_without_bucket) = f64_target_vector(
            covered_target
                .iter()
                .zip(report.correction.iter())
                .map(|(covered, bucket)| covered.get() - bucket),
        ) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV branch bucket parity_mod2={} qdot_bin={} unavailable: non-finite drop vector",
                report.parity_mod2, report.q_dot_bucket
            );
            continue;
        };
        let input_drop_summary = target_correction_delta_summary(
            input_chi_reference,
            &target_without_bucket,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare input-chi corrected-chamber GV branch bucket drop: {e}"
            );
            std::process::exit(2);
        });
        let corrected_drop_summary =
            target_correction_delta_summary(corrected_chi_reference, &target_without_bucket)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare corrected-chi corrected-chamber GV branch bucket drop: {e}"
                    );
                    std::process::exit(2);
                });
        let Some(target_with_flipped_bucket) = f64_target_vector(
            covered_target
                .iter()
                .zip(report.correction.iter())
                .map(|(covered, bucket)| covered.get() - 2.0 * bucket),
        ) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV branch bucket parity_mod2={} qdot_bin={} unavailable: non-finite flip vector",
                report.parity_mod2, report.q_dot_bucket
            );
            continue;
        };
        let input_flip_summary = target_correction_delta_summary(
            input_chi_reference,
            &target_with_flipped_bucket,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare input-chi corrected-chamber GV branch bucket flip: {e}"
            );
            std::process::exit(2);
        });
        let corrected_flip_summary = target_correction_delta_summary(
            corrected_chi_reference,
            &target_with_flipped_bucket,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare corrected-chi corrected-chamber GV branch bucket flip: {e}"
            );
            std::process::exit(2);
        });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV branch bucket parity_mod2={} qdot_bin={} count={} qdot_min={} qdot_max={} bucket_max_abs={} bucket_l2={} drop_input_chi_max_abs={} drop_input_chi_relative_l2={} drop_corrected_chi_max_abs={} drop_corrected_chi_relative_l2={} flip_input_chi_max_abs={} flip_input_chi_relative_l2={} flip_corrected_chi_max_abs={} flip_corrected_chi_relative_l2={}",
            report.parity_mod2,
            report.q_dot_bucket,
            report.count,
            report.q_dot_min,
            report.q_dot_max,
            report.max_abs_correction,
            report.l2_correction,
            input_drop_summary.max_abs_delta,
            input_drop_summary.relative_l2_delta,
            corrected_drop_summary.max_abs_delta,
            corrected_drop_summary.relative_l2_delta,
            input_flip_summary.max_abs_delta,
            input_flip_summary.relative_l2_delta,
            corrected_flip_summary.max_abs_delta,
            corrected_flip_summary.relative_l2_delta
        );
    }
}

fn report_corrected_chamber_toric_gv_source_buckets(
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    toric_gv_diagnostic_by_class: &HashMap<Vec<i64>, cyrus_core::ToricCurveGvDiagnostic>,
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: &[I64<Finite>],
    input_chi_reference: &[F64<Finite>],
    corrected_chi_reference: &[F64<Finite>],
    covered_target: &[F64<Finite>],
) {
    if basis.len() != t.len() || kklt_basis.len() != covered_target.len() {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV source buckets unavailable: dimension mismatch"
        );
        return;
    }

    let mut buckets: BTreeMap<String, GvSourceBucketAccumulator> = BTreeMap::new();
    for (curve, invariant) in gv_invariants {
        if kklt_basis.iter().any(|&idx| idx >= curve.len())
            || basis.iter().any(|&idx| idx >= curve.len())
        {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV source buckets unavailable: curve dimension mismatch"
            );
            return;
        }
        let q_dot_t = basis
            .iter()
            .zip(t.iter())
            .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
            .sum::<f64>();
        let Some(parity) = ambient_curve_b_field_parity_diagnostic(curve, basis, gamma) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV source buckets unavailable: gamma dimension mismatch"
            );
            return;
        };
        let single = [(curve.clone(), invariant.clone())];
        let Some(contribution) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &single,
            basis,
            kklt_basis,
            t,
            Some(gamma),
        ) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV source buckets unavailable: invalid single-curve contribution q_dot_t={q_dot_t} parity_mod2={}",
                parity.rem_euclid(2)
            );
            return;
        };
        let (bucket_key, source_count) =
            if let Some(diagnostic) = toric_gv_diagnostic_by_class.get(curve) {
                (
                    toric_curve_gv_source_bucket_label(diagnostic),
                    diagnostic.sources.len(),
                )
            } else {
                ("missing_diagnostic".to_string(), 0)
            };
        let bucket = buckets
            .entry(bucket_key)
            .or_insert_with(|| GvSourceBucketAccumulator {
                count: 0,
                source_count_min: usize::MAX,
                source_count_max: 0,
                q_dot_min: f64::INFINITY,
                q_dot_max: f64::NEG_INFINITY,
                correction: vec![0.0; kklt_basis.len()],
            });
        bucket.count += 1;
        bucket.source_count_min = bucket.source_count_min.min(source_count);
        bucket.source_count_max = bucket.source_count_max.max(source_count);
        bucket.q_dot_min = bucket.q_dot_min.min(q_dot_t);
        bucket.q_dot_max = bucket.q_dot_max.max(q_dot_t);
        for (accumulated, value) in bucket.correction.iter_mut().zip(contribution.iter()) {
            *accumulated += value.get();
        }
    }

    let mut reports = buckets
        .into_iter()
        .map(|(label, bucket)| {
            let (max_abs_correction, l2_correction) = max_abs_and_l2(&bucket.correction);
            GvSourceBucketReport {
                label,
                count: bucket.count,
                source_count_min: bucket.source_count_min,
                source_count_max: bucket.source_count_max,
                q_dot_min: bucket.q_dot_min,
                q_dot_max: bucket.q_dot_max,
                max_abs_correction,
                l2_correction,
                correction: bucket.correction,
            }
        })
        .collect::<Vec<_>>();
    reports.sort_unstable_by(|lhs, rhs| {
        rhs.l2_correction
            .total_cmp(&lhs.l2_correction)
            .then_with(|| lhs.label.cmp(&rhs.label))
    });

    for report in reports.into_iter().take(16) {
        let Some(target_without_bucket) = f64_target_vector(
            covered_target
                .iter()
                .zip(report.correction.iter())
                .map(|(covered, bucket)| covered.get() - bucket),
        ) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV source bucket label={} unavailable: non-finite drop vector",
                report.label
            );
            continue;
        };
        let input_drop_summary = target_correction_delta_summary(
            input_chi_reference,
            &target_without_bucket,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare input-chi corrected-chamber GV source bucket drop: {e}"
            );
            std::process::exit(2);
        });
        let corrected_drop_summary =
            target_correction_delta_summary(corrected_chi_reference, &target_without_bucket)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare corrected-chi corrected-chamber GV source bucket drop: {e}"
                    );
                    std::process::exit(2);
                });
        let Some(target_with_flipped_bucket) = f64_target_vector(
            covered_target
                .iter()
                .zip(report.correction.iter())
                .map(|(covered, bucket)| covered.get() - 2.0 * bucket),
        ) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV source bucket label={} unavailable: non-finite flip vector",
                report.label
            );
            continue;
        };
        let input_flip_summary = target_correction_delta_summary(
            input_chi_reference,
            &target_with_flipped_bucket,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare input-chi corrected-chamber GV source bucket flip: {e}"
            );
            std::process::exit(2);
        });
        let corrected_flip_summary =
            target_correction_delta_summary(corrected_chi_reference, &target_with_flipped_bucket)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare corrected-chi corrected-chamber GV source bucket flip: {e}"
                    );
                    std::process::exit(2);
                });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV source bucket label={} count={} source_count_min={} source_count_max={} qdot_min={} qdot_max={} bucket_max_abs={} bucket_l2={} drop_input_chi_max_abs={} drop_input_chi_relative_l2={} drop_corrected_chi_max_abs={} drop_corrected_chi_relative_l2={} flip_input_chi_max_abs={} flip_input_chi_relative_l2={} flip_corrected_chi_max_abs={} flip_corrected_chi_relative_l2={}",
            report.label,
            report.count,
            report.source_count_min,
            report.source_count_max,
            report.q_dot_min,
            report.q_dot_max,
            report.max_abs_correction,
            report.l2_correction,
            input_drop_summary.max_abs_delta,
            input_drop_summary.relative_l2_delta,
            corrected_drop_summary.max_abs_delta,
            corrected_drop_summary.relative_l2_delta,
            input_flip_summary.max_abs_delta,
            input_flip_summary.relative_l2_delta,
            corrected_flip_summary.max_abs_delta,
            corrected_flip_summary.relative_l2_delta
        );
    }
}

fn ambient_target_contribution_vectors(
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: &[I64<Finite>],
) -> Result<Vec<Vec<f64>>, String> {
    let mut contributions = Vec::with_capacity(gv_invariants.len());
    for (curve_index, (curve, invariant)) in gv_invariants.iter().enumerate() {
        let single = [(curve.clone(), invariant.clone())];
        let Some(contribution) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &single,
            basis,
            kklt_basis,
            t,
            Some(gamma),
        ) else {
            return Err(format!(
                "failed to compute single-curve contribution vector for curve_index={curve_index} class={:?}",
                sparse_i64(curve)
            ));
        };
        contributions.push(contribution.iter().map(|value| value.get()).collect());
    }
    Ok(contributions)
}

fn solve_gv_weight_lp(
    label: &'static str,
    contributions: &[Vec<f64>],
    target: &[F64<Finite>],
    lower_bound: f64,
    upper_bound: f64,
) -> Result<GvWeightLpReport, String> {
    if contributions.is_empty() {
        return Err("GV weight LP requires at least one contribution vector".to_string());
    }
    if !(lower_bound.is_finite() && upper_bound.is_finite() && lower_bound <= upper_bound) {
        return Err(format!(
            "invalid GV weight LP bounds: lower={lower_bound} upper={upper_bound}"
        ));
    }
    let dim = target.len();
    if dim == 0 || contributions.iter().any(|row| row.len() != dim) {
        return Err(format!(
            "GV weight LP dimension mismatch: target_dim={} contribution_dims={:?}",
            dim,
            contributions
                .iter()
                .take(8)
                .map(Vec::len)
                .collect::<Vec<_>>()
        ));
    }

    let mut vars = ProblemVariables::new();
    let weights = contributions
        .iter()
        .map(|_| vars.add(variable().min(lower_bound).max(upper_bound)))
        .collect::<Vec<_>>();
    let epsilon = vars.add(variable().min(0.0));
    let mut objective = Expression::from(0.0);
    objective.add_mul(1.0, epsilon);
    let mut model = vars.minimise(objective).using(default_solver);

    for coord in 0..dim {
        let target_value = target[coord].get();

        let mut upper = Expression::from(-target_value);
        for (weight, contribution) in weights.iter().zip(contributions.iter()) {
            let coefficient = contribution[coord];
            if coefficient != 0.0 {
                upper.add_mul(coefficient, *weight);
            }
        }
        upper.add_mul(-1.0, epsilon);
        model = model.with(upper.leq(0.0));

        let mut lower = Expression::from(target_value);
        for (weight, contribution) in weights.iter().zip(contributions.iter()) {
            let coefficient = contribution[coord];
            if coefficient != 0.0 {
                lower.add_mul(-coefficient, *weight);
            }
        }
        lower.add_mul(-1.0, epsilon);
        model = model.with(lower.leq(0.0));
    }

    let solution = match model.solve() {
        Ok(solution) => solution,
        Err(ResolutionError::Infeasible) => {
            return Err("GV weight LP is unexpectedly infeasible".to_string());
        }
        Err(err) => return Err(format!("GV weight LP solve failed for {label}: {err}")),
    };

    let weights = weights
        .iter()
        .map(|weight| solution.value(*weight))
        .collect::<Vec<_>>();
    if weights.iter().any(|value| !value.is_finite()) {
        return Err("GV weight LP produced a non-finite weight".to_string());
    }

    let mut fitted = vec![0.0f64; dim];
    for (weight, contribution) in weights.iter().zip(contributions.iter()) {
        for (entry, value) in fitted.iter_mut().zip(contribution.iter()) {
            *entry += *weight * *value;
        }
    }

    let residuals = target
        .iter()
        .zip(fitted.iter())
        .map(|(reference, candidate)| candidate - reference.get())
        .collect::<Vec<_>>();
    let (max_abs_delta, residual_l2) = max_abs_and_l2(&residuals);
    let (_, target_l2) =
        max_abs_and_l2(&target.iter().map(|value| value.get()).collect::<Vec<_>>());
    let relative_l2_delta = if target_l2 == 0.0 {
        residual_l2
    } else {
        residual_l2 / target_l2
    };

    let weight_min = weights.iter().copied().fold(f64::INFINITY, f64::min);
    let weight_max = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let weight_mean = weights.iter().sum::<f64>() / weights.len() as f64;
    let tolerance = 1e-7;
    let near_lower_count = weights
        .iter()
        .filter(|&&weight| weight <= lower_bound + tolerance)
        .count();
    let near_upper_count = weights
        .iter()
        .filter(|&&weight| weight >= upper_bound - tolerance)
        .count();
    let interior_count = weights
        .len()
        .saturating_sub(near_lower_count + near_upper_count);
    let max_abs_delta_from_one = weights
        .iter()
        .map(|weight| (weight - 1.0).abs())
        .fold(0.0f64, f64::max);

    Ok(GvWeightLpReport {
        label,
        lower_bound,
        upper_bound,
        objective_epsilon: solution.value(epsilon),
        max_abs_delta,
        relative_l2_delta,
        weight_min,
        weight_max,
        weight_mean,
        near_lower_count,
        near_upper_count,
        interior_count,
        max_abs_delta_from_one,
        weights,
    })
}

fn solve_gf2_system(rows: &[Vec<u8>], rhs: &[u8], n_vars: usize) -> Option<Vec<u8>> {
    if rows.len() != rhs.len() || rows.iter().any(|row| row.len() != n_vars) {
        return None;
    }

    let mut augmented = rows
        .iter()
        .zip(rhs.iter())
        .map(|(row, &rhs_value)| {
            let mut augmented_row = row.clone();
            augmented_row.push(rhs_value & 1);
            augmented_row
        })
        .collect::<Vec<_>>();
    let mut pivot_columns = Vec::new();
    let mut pivot_row = 0usize;
    for col in 0..n_vars {
        let Some(found_row) = (pivot_row..augmented.len()).find(|&row| augmented[row][col] == 1)
        else {
            continue;
        };
        augmented.swap(pivot_row, found_row);
        for row in 0..augmented.len() {
            if row != pivot_row && augmented[row][col] == 1 {
                for entry in col..=n_vars {
                    augmented[row][entry] ^= augmented[pivot_row][entry];
                }
            }
        }
        pivot_columns.push(col);
        pivot_row += 1;
        if pivot_row == augmented.len() {
            break;
        }
    }

    if augmented[pivot_row..]
        .iter()
        .any(|row| row[..n_vars].iter().all(|&value| value == 0) && row[n_vars] == 1)
    {
        return None;
    }

    let mut solution = vec![0u8; n_vars];
    for (row_idx, &col) in pivot_columns.iter().enumerate() {
        solution[col] = augmented[row_idx][n_vars] & 1;
    }
    Some(solution)
}

fn report_bounded_sign_lp_parity_fit(
    report: &GvWeightLpReport,
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: &[I64<Finite>],
    checkpoint_implied_gv: &[F64<Finite>],
) {
    let ambient_dim = gamma.len();
    let tolerance = 1e-7;
    let mut rows = Vec::new();
    let mut rhs = Vec::new();
    let mut negative_constraints = 0usize;
    let mut positive_constraints = 0usize;
    for ((curve, _invariant), &weight) in gv_invariants.iter().zip(report.weights.iter()) {
        let desired_flip = if weight <= -1.0 + tolerance {
            negative_constraints += 1;
            Some(1u8)
        } else if weight >= 1.0 - tolerance {
            positive_constraints += 1;
            Some(0u8)
        } else {
            None
        };
        let Some(desired_flip) = desired_flip else {
            continue;
        };
        if curve.len() != ambient_dim {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV weight LP parity-fit unavailable: curve dimension {} != gamma dimension {}",
                curve.len(),
                ambient_dim
            );
            return;
        }
        rows.push(
            curve
                .iter()
                .map(|coefficient| coefficient.rem_euclid(2) as u8)
                .collect::<Vec<_>>(),
        );
        rhs.push(desired_flip);
    }

    let constraints = rows.len();
    let Some(delta_gamma_mod2) = solve_gf2_system(&rows, &rhs, ambient_dim) else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV weight LP parity-fit label={} compatible=false constraints={} negative_constraints={} positive_constraints={}",
            report.label, constraints, negative_constraints, positive_constraints
        );
        return;
    };

    let mut candidate_gamma = gamma.to_vec();
    for (entry, &delta) in candidate_gamma.iter_mut().zip(delta_gamma_mod2.iter()) {
        if delta == 1 {
            let Some(new_value) = entry.get().checked_add(1) else {
                eprintln!(
                    "[COMPARE] checkpoint-t corrected-chamber GV weight LP parity-fit unavailable: gamma overflow"
                );
                return;
            };
            *entry = I64::<Finite>::new(new_value);
        }
    }
    let support = delta_gamma_mod2.iter().filter(|&&value| value == 1).count();
    let Some(candidate_target) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
        gv_invariants,
        basis,
        kklt_basis,
        t,
        Some(&candidate_gamma),
    ) else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV weight LP parity-fit label={} compatible=true support={} but candidate target correction is invalid",
            report.label, support
        );
        return;
    };
    let summary = target_correction_delta_summary(checkpoint_implied_gv, &candidate_target)
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare bounded-sign LP parity-fit GV target correction: {e}"
            );
            std::process::exit(2);
        });
    eprintln!(
        "[COMPARE] checkpoint-t corrected-chamber GV weight LP parity-fit label={} compatible=true constraints={} negative_constraints={} positive_constraints={} delta_gamma_support={} max_abs={} relative_l2={} max_abs_candidate={}",
        report.label,
        constraints,
        negative_constraints,
        positive_constraints,
        support,
        summary.max_abs_delta,
        summary.relative_l2_delta,
        summary.max_abs_candidate
    );
}

fn report_corrected_chamber_gv_weight_lp(
    gv_invariants: &[(Vec<i64>, malachite::Integer)],
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: &[I64<Finite>],
    checkpoint_implied_gv: &[F64<Finite>],
) {
    let contributions =
        match ambient_target_contribution_vectors(gv_invariants, basis, kklt_basis, t, gamma) {
            Ok(contributions) => contributions,
            Err(err) => {
                eprintln!(
                    "[COMPARE] checkpoint-t corrected-chamber GV weight LP unavailable: {err}"
                );
                return;
            }
        };
    for (label, lower_bound, upper_bound) in
        [("continuous_drop", 0.0, 1.0), ("bounded_sign", -1.0, 1.0)]
    {
        let report = solve_gv_weight_lp(
            label,
            &contributions,
            checkpoint_implied_gv,
            lower_bound,
            upper_bound,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed corrected-chamber GV weight LP diagnostic: {e}");
            std::process::exit(2);
        });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV weight LP label={} bounds=[{},{}] eps={} max_abs={} relative_l2={} weight_min={} weight_max={} weight_mean={} near_lower={} near_upper={} interior={} max_abs_delta_from_one={}",
            report.label,
            report.lower_bound,
            report.upper_bound,
            report.objective_epsilon,
            report.max_abs_delta,
            report.relative_l2_delta,
            report.weight_min,
            report.weight_max,
            report.weight_mean,
            report.near_lower_count,
            report.near_upper_count,
            report.interior_count,
            report.max_abs_delta_from_one
        );

        let mut weight_deltas = report
            .weights
            .iter()
            .enumerate()
            .map(|(idx, &weight)| (idx, (weight - 1.0).abs(), weight))
            .collect::<Vec<_>>();
        weight_deltas.sort_unstable_by(|lhs, rhs| rhs.1.total_cmp(&lhs.1));
        for (curve_index, _delta, weight) in weight_deltas.into_iter().take(6) {
            let (curve, invariant) = &gv_invariants[curve_index];
            let q_dot_t = basis
                .iter()
                .zip(t.iter())
                .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
                .sum::<f64>();
            let parity = ambient_curve_b_field_parity_diagnostic(curve, basis, gamma)
                .map(|value| value.rem_euclid(2))
                .unwrap_or(-1);
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV weight LP top_weight label={} curve_index={} weight={} q_dot_t={} parity_mod2={} gv={} class={:?}",
                report.label,
                curve_index,
                weight,
                q_dot_t,
                parity,
                invariant,
                sparse_i64(curve)
            );
        }
        if report.label == "bounded_sign" {
            report_bounded_sign_lp_parity_fit(
                &report,
                gv_invariants,
                basis,
                kklt_basis,
                t,
                gamma,
                checkpoint_implied_gv,
            );
        }
    }
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

fn secondary_cone_height_certificate_for_kahler(
    tri: &Triangulation,
    geom: &PrimalGeom,
    basis: &[usize],
    kahler: &[F64<Finite>],
) -> Result<SecondaryConeHeightCertificate, String> {
    let heights = secondary_cone_typed_heights_for_kahler(geom, basis, kahler)?;
    let hyperplanes = secondary_cone_hyperplanes_native(&geom.triangulation_points, tri)
        .map_err(|e| format!("failed to compute secondary-cone hyperplanes: {e}"))?;
    secondary_cone_height_certificate_from_hyperplanes(&hyperplanes, &heights)
}

fn secondary_cone_2face_height_certificate_for_kahler(
    tri: &Triangulation,
    geom: &PrimalGeom,
    basis: &[usize],
    kahler: &[F64<Finite>],
) -> Result<SecondaryConeHeightCertificate, String> {
    let heights = secondary_cone_typed_heights_for_kahler(geom, basis, kahler)?;
    let hyperplanes = secondary_cone_hyperplanes_native_on_polytope_2faces_4d(
        &geom.triangulation_points,
        tri,
        &geom.polytope,
    )
    .map_err(|e| format!("failed to compute 2-face secondary-cone hyperplanes: {e}"))?;
    secondary_cone_height_certificate_from_hyperplanes(&hyperplanes, &heights)
}

fn expanded_secondary_fan_height_certificate_for_kahler(
    geom: &PrimalGeom,
    basis: &[usize],
    kahler: &[F64<Finite>],
) -> Result<SecondaryConeHeightCertificate, String> {
    let heights = secondary_cone_typed_heights_for_kahler(geom, basis, kahler)?;
    let hyperplanes = expanded_secondary_fan_hyperplanes_on_polytope_2faces_4d(
        &geom.triangulation_points,
        &geom.polytope,
    )
    .map_err(|e| format!("failed to compute expanded-secondary fan hyperplanes: {e}"))?;
    secondary_cone_height_certificate_from_hyperplanes(&hyperplanes, &heights)
}

fn corrected_chamber_face_triangulation_choice_summary(
    geom: &PrimalGeom,
    max_exact_face_points: usize,
    samples_per_large_face: usize,
    max_sampling_attempts_per_face: usize,
    seed: u64,
) -> Result<FaceTriangulationChoiceSummary, String> {
    let faces = geom
        .polytope
        .faces_4d_for_points(&geom.triangulation_points)
        .map_err(|e| format!("failed to compute 4D two-faces for face choices: {e}"))?;
    let face_point_counts = faces.twofaces.iter().map(Vec::len).collect::<Vec<_>>();
    let choices = fine_regular_triangulation_choices_on_polytope_2faces_4d_with_sampling(
        &geom.triangulation_points,
        &geom.polytope,
        max_exact_face_points,
        samples_per_large_face,
        max_sampling_attempts_per_face,
        seed,
    )
    .map_err(|e| format!("failed to compute 4D two-face triangulation choices: {e}"))?;
    let choice_counts = choices.iter().map(Vec::len).collect::<Vec<_>>();
    if choice_counts.len() != face_point_counts.len() {
        return Err(format!(
            "face choice count mismatch: faces={} choice_blocks={}",
            face_point_counts.len(),
            choice_counts.len()
        ));
    }

    let sampled_face_indices = face_point_counts
        .iter()
        .enumerate()
        .filter_map(|(idx, &point_count)| (point_count > max_exact_face_points).then_some(idx))
        .collect::<Vec<_>>();
    let sampled_face_point_counts = sampled_face_indices
        .iter()
        .map(|&idx| face_point_counts[idx])
        .collect::<Vec<_>>();
    let exact_face_count = face_point_counts
        .len()
        .checked_sub(sampled_face_indices.len())
        .ok_or_else(|| "sampled face count exceeds face count".to_string())?;
    let empty_choice_face_count = choice_counts
        .iter()
        .filter(|&&choice_count| choice_count == 0)
        .count();
    let total_choice_count = choice_counts
        .iter()
        .fold(malachite::Integer::from(1), |acc, &choice_count| {
            acc * malachite::Integer::from(choice_count)
        });
    let status = if empty_choice_face_count > 0 {
        "face_triangulation_choices_have_empty_face_blocks"
    } else if sampled_face_indices.is_empty() {
        "face_triangulation_choices_exact"
    } else {
        "face_triangulation_choices_exact_and_sampled"
    };

    Ok(FaceTriangulationChoiceSummary {
        status: status.to_string(),
        max_exact_face_points,
        samples_per_large_face,
        max_sampling_attempts_per_face,
        seed,
        face_count: face_point_counts.len(),
        exact_face_count,
        sampled_face_count: sampled_face_indices.len(),
        empty_choice_face_count,
        min_face_points: face_point_counts.iter().copied().min(),
        max_face_points: face_point_counts.iter().copied().max(),
        min_choice_count: choice_counts.iter().copied().min(),
        max_choice_count: choice_counts.iter().copied().max(),
        total_choice_count: total_choice_count.to_string(),
        face_point_counts,
        choice_counts,
        sampled_face_indices,
        sampled_face_point_counts,
    })
}

fn secondary_cone_height_certificate_from_hyperplanes(
    hyperplanes: &[Vec<i64>],
    heights: &[F64<Finite>],
) -> Result<SecondaryConeHeightCertificate, String> {
    let pairings = secondary_cone_height_pairings(&hyperplanes, &heights)
        .map_err(|e| format!("failed to pair secondary-cone hyperplanes with heights: {e}"))?;
    let epsilon = F64::<Pos>::new(1e-6).expect("CYTools secondary-cone epsilon is positive");
    let strictly_inside = pairings.iter().all(|pairing| pairing.get() > epsilon.get());
    let min_pairing = pairings
        .iter()
        .map(|pairing| pairing.get())
        .reduce(f64::min);
    let max_pairing = pairings
        .iter()
        .map(|pairing| pairing.get())
        .reduce(f64::max);
    let status = if hyperplanes.is_empty() {
        "no_secondary_cone_hyperplanes"
    } else if strictly_inside {
        "strictly_inside_secondary_cone"
    } else {
        "height_vector_on_or_outside_secondary_cone"
    };

    Ok(SecondaryConeHeightCertificate {
        status: status.to_string(),
        epsilon: epsilon.get(),
        hyperplane_count: hyperplanes.len(),
        pairing_count: pairings.len(),
        min_pairing,
        max_pairing,
        strictly_inside,
    })
}

fn secondary_cone_heights_for_kahler(
    geom: &PrimalGeom,
    basis: &[usize],
    kahler: &[F64<Finite>],
) -> Result<Vec<f64>, String> {
    Ok(
        secondary_cone_typed_heights_for_kahler(geom, basis, kahler)?
            .iter()
            .map(|height| height.get())
            .collect(),
    )
}

fn secondary_cone_typed_heights_for_kahler(
    geom: &PrimalGeom,
    basis: &[usize],
    kahler: &[F64<Finite>],
) -> Result<Vec<F64<Finite>>, String> {
    let basis_non_origin = basis
        .iter()
        .map(|&idx| {
            idx.checked_sub(1).ok_or_else(|| {
                "secondary-cone Kähler basis unexpectedly contains origin".to_string()
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let non_origin_count = geom
        .triangulation_points
        .len()
        .checked_sub(1)
        .ok_or_else(|| "secondary-cone point set is empty".to_string())?;
    let heights = kahler_to_heights(kahler, &basis_non_origin, non_origin_count)
        .ok_or_else(|| "failed to embed Kähler point into secondary-cone heights".to_string())?;
    Ok(heights)
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
    small_curve_pruning: CurvePruningStrategy,
    general_min_points: Option<u32>,
    general_max_deg: Option<u32>,
    provided_generators_only: bool,
    ray_gv_requested: bool,
    lp_face_gv_requested: bool,
    missing_target_sample_limit: usize,
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
    let small_curves = prune_selected_curve_candidates(
        &small_curve_candidates,
        small_curve_pruning,
        "corrected-chamber",
    )?;
    let toric_gv_diagnostic_by_class: HashMap<Vec<i64>, cyrus_core::ToricCurveGvDiagnostic> =
        compute_toric_curve_gv_diagnostics(tri, &geom.triangulation_points, &geom.polytope)
            .map_err(|e| format!("failed to compute corrected-chamber toric GV values: {e}"))?
            .into_iter()
            .map(|item| (item.class.clone(), item))
            .collect();
    let mut gv_by_class: HashMap<Vec<i64>, malachite::Integer> = toric_gv_diagnostic_by_class
        .iter()
        .map(|(class, diagnostic)| (class.clone(), diagnostic.gv.clone()))
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
    let secondary_cone_height_certificate = Some(secondary_cone_height_certificate_for_kahler(
        tri,
        geom,
        &intersection.basis,
        kahler,
    )?);
    let secondary_cone_2face_height_certificate = Some(
        secondary_cone_2face_height_certificate_for_kahler(tri, geom, &intersection.basis, kahler)?,
    );
    let expanded_secondary_fan_height_certificate = Some(
        expanded_secondary_fan_height_certificate_for_kahler(geom, &intersection.basis, kahler)?,
    );
    let secondary_cone_heights_for_missing = Some(secondary_cone_heights_for_kahler(
        geom,
        &intersection.basis,
        kahler,
    )?);

    let mut basis_ray_stats = None;
    let mut basis_rays_for_missing = None;
    let mut basis_rays_for_missing_degree_bound = None;
    let mut basis_rays_for_missing_degree_bounded = None;
    let mut degree_bounded_mori_ray_context_for_missing = None;
    let mut covered_toric_gv_context_for_missing = None;
    let mut gv_basis_data_for_missing = None;
    let mut grading_for_missing = None;
    let mut corrected_kappa_basis_for_missing = None;
    let mut degree_summary = None;
    let mut missing_target_stats = None;
    let mut uncovered_source_ray_stats_degree_bound_for_missing = None;
    let mut uncovered_source_ray_stats_for_missing = None;
    let mut shared_facet_unresolved_source_ray_stats_for_missing = None;
    let mut uncovered_source_ray_toric_diagnostic_sample = None;
    let mut degree_bounded_toric_gv_diagnostic_context_for_missing = None;
    let mut covered_toric_gv_divisor_representation_baseline = None;
    if !missing_gv_classes.is_empty() {
        let origin_idx = geom
            .triangulation_points
            .iter()
            .position(|point| point.coords().iter().all(|&coord| coord == 0))
            .ok_or_else(|| "failed to find origin in corrected-chamber points".to_string())?;
        let gv_basis_data = vector_gv_basis_data(
            &ambient_rays,
            &intersection.linrels,
            &intersection.basis,
            "corrected-chamber",
        )?;
        let basis_rays = gv_basis_data.mori_rays.clone();
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
        let corrected_kappa_basis =
            intersection_in_basis(&corrected_kappa_full, &intersection.basis);
        corrected_kappa_basis_for_missing =
            Some(sparse_intersection_entries(&corrected_kappa_basis));
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
        let degree_bounded_basis_rays = basis_rays
            .iter()
            .filter(|ray| {
                ray.iter()
                    .zip(grading.iter())
                    .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
                    .sum::<i128>()
                    <= summary.max_degree
            })
            .cloned()
            .collect::<Vec<_>>();
        let degree_bounded_ambient_ray_context = degree_bounded_mori_ray_context_samples(
            &ambient_rays,
            &intersection.basis,
            &grading,
            summary.max_degree,
        )?;
        let ambient_dimension = ambient_rays
            .first()
            .map(Vec::len)
            .ok_or_else(|| "corrected-chamber Mori-cap ambient ray context is empty".to_string())?;
        let degree_bounded_ambient_classes = degree_bounded_ambient_ray_context
            .iter()
            .enumerate()
            .map(|(idx, sample)| {
                dense_i64_from_sparse(
                    &sample.ambient_nonzero,
                    ambient_dimension,
                    &format!("degree-bounded Mori ray context sample {idx} ambient class"),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        degree_bounded_toric_gv_diagnostic_context_for_missing =
            Some(toric_gv_diagnostic_context_samples(
                &degree_bounded_ambient_classes,
                &toric_gv_diagnostic_by_class,
                &intersection.basis,
                &grading,
                usize::MAX,
            )?);
        let lower_source_ray_degree_bound = summary.min_degree.saturating_sub(1);
        let mut uncovered_source_ray_classes = uncovered_degree_bounded_source_ray_classes(
            &ambient_rays,
            &small_curve_gvs,
            &intersection.basis,
            &grading,
            lower_source_ray_degree_bound,
        )?;
        let covered_toric_gv_context = covered_toric_gv_context_samples(
            &small_curve_gvs,
            &intersection.basis,
            &grading,
            summary.max_degree,
        )?;
        let ray_stats = graded_ray_stats(&basis_rays, &grading, general_max_deg)?;
        let target_stats = missing_gv_target_stats(
            &missing_gv_classes,
            &geom.triangulation_points,
            Some(tri),
            &basis_rays,
            &intersection.basis,
            &grading,
            Some(kahler),
            Some(gamma),
            origin_idx,
            &origin_circuits_by_class,
            &cms_intersection_checks_by_class,
            true,
            missing_target_sample_limit,
        )?;
        let active_dependency_source_ray_classes = active_noncovered_dependency_source_ray_classes(
            &target_stats,
            &degree_bounded_ambient_ray_context,
            &covered_toric_gv_context,
            ambient_rays.first().map(Vec::len).ok_or_else(|| {
                "corrected-chamber Mori-cap ambient ray context is empty".to_string()
            })?,
        )?;
        uncovered_source_ray_classes = merge_unique_ambient_classes(
            uncovered_source_ray_classes,
            active_dependency_source_ray_classes,
        );
        uncovered_source_ray_toric_diagnostic_sample = Some(toric_gv_diagnostic_context_samples(
            &uncovered_source_ray_classes,
            &toric_gv_diagnostic_by_class,
            &intersection.basis,
            &grading,
            missing_target_sample_limit,
        )?);
        let uncovered_source_ray_cms_intersection_checks_by_class =
            cms_general_divisor_intersection_checks_by_class(
                &uncovered_source_ray_classes,
                &origin_circuits_by_class,
                &corrected_kappa_full,
                &intersection.basis,
            )?;
        let uncovered_source_ray_stats = missing_gv_target_stats(
            &uncovered_source_ray_classes,
            &geom.triangulation_points,
            Some(tri),
            &basis_rays,
            &intersection.basis,
            &grading,
            Some(kahler),
            Some(gamma),
            origin_idx,
            &origin_circuits_by_class,
            &uncovered_source_ray_cms_intersection_checks_by_class,
            false,
            missing_target_sample_limit,
        )?;
        let shared_facet_unresolved_source_ray_classes =
            shared_facet_unresolved_source_ray_classes(
                &target_stats,
                &degree_bounded_ambient_ray_context,
                &covered_toric_gv_context,
                degree_bounded_toric_gv_diagnostic_context_for_missing
                    .as_deref()
                    .unwrap_or(&[]),
                &uncovered_source_ray_stats,
                ambient_dimension,
            )?;
        let shared_facet_unresolved_source_ray_cms_intersection_checks_by_class =
            cms_general_divisor_intersection_checks_by_class(
                &shared_facet_unresolved_source_ray_classes,
                &origin_circuits_by_class,
                &corrected_kappa_full,
                &intersection.basis,
            )?;
        let shared_facet_unresolved_source_ray_stats = missing_gv_target_stats(
            &shared_facet_unresolved_source_ray_classes,
            &geom.triangulation_points,
            Some(tri),
            &basis_rays,
            &intersection.basis,
            &grading,
            Some(kahler),
            Some(gamma),
            origin_idx,
            &origin_circuits_by_class,
            &shared_facet_unresolved_source_ray_cms_intersection_checks_by_class,
            false,
            missing_target_sample_limit,
        )?;
        basis_rays_for_missing_degree_bound = Some(summary.max_degree);
        basis_rays_for_missing_degree_bounded = Some(degree_bounded_basis_rays);
        degree_bounded_mori_ray_context_for_missing = Some(degree_bounded_ambient_ray_context);
        covered_toric_gv_context_for_missing = Some(covered_toric_gv_context);
        basis_rays_for_missing = Some(basis_rays);
        gv_basis_data_for_missing = Some(gv_basis_data);
        grading_for_missing = Some(grading);
        degree_summary = Some(summary);
        basis_ray_stats = Some(ray_stats);
        missing_target_stats = Some(target_stats);
        uncovered_source_ray_stats_degree_bound_for_missing = Some(lower_source_ray_degree_bound);
        uncovered_source_ray_stats_for_missing = Some(uncovered_source_ray_stats);
        shared_facet_unresolved_source_ray_stats_for_missing =
            Some(shared_facet_unresolved_source_ray_stats);
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
        let gv_basis_data = gv_basis_data_for_missing
            .as_ref()
            .expect("GV basis data computed for corrected-chamber missing curves");
        let corrected_kappa_basis =
            chamber_intersection_in_basis(tri, &geom.triangulation_points, &intersection.basis)?;
        eprintln!(
            "[WARN] corrected-chamber one-dimensional ray GV diagnostic assumes each LP-extremal primitive target spans a valid Mori-cone face; this is not yet promoted to the exact corrected-chamber GV fallback."
        );
        let ray_gvs = compute_missing_one_dimensional_ray_gvs(
            &ray_targets.candidates,
            &intersection.basis,
            grading,
            &gv_basis_data.q_matrix,
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
        let gv_basis_data = gv_basis_data_for_missing
            .as_ref()
            .expect("GV basis data computed for corrected-chamber missing curves");
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
            &gv_basis_data.q_matrix,
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
        let gv_basis_data = gv_basis_data_for_missing
            .as_ref()
            .expect("GV basis data computed for corrected-chamber missing curves");
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
                &gv_basis_data.q_matrix,
                &corrected_kappa_basis,
                general_min_points,
                general_max_deg,
            )
        } else {
            cyrus_core::compute_gv_invariants(
                basis_rays,
                grading,
                &gv_basis_data.q_matrix,
                &corrected_kappa_basis,
                general_min_points,
                general_max_deg,
            )
        }
        .map_err(|e| format!("failed to compute corrected-chamber general GV invariants: {e}"))?;
        let ambient_gvs =
            map_basis_gv_invariants_to_ambient(&general_gvs, &gv_basis_data.curve_basis_matrix)
                .map_err(|e| {
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

    let ambient_dim = geom.triangulation_points.len();
    let formula_sum_diagnostic_gvs =
        formula_sum_diagnostic_missing_gvs(missing_target_stats.as_ref(), ambient_dim)?;
    let mut formula_sum_diagnostic_gv_covered_count = None;
    let mut formula_sum_diagnostic_gv_missing_count = None;
    let mut formula_sum_diagnostic_gv_zero_count = None;
    let mut formula_sum_diagnostic_gv_nonzero_count = None;
    let mut formula_sum_diagnostic_gv_volume_correction = None;
    let mut formula_sum_diagnostic_gv_target_correction = None;
    if !formula_sum_diagnostic_gvs.is_empty() {
        if let Some(unexpected) = formula_sum_diagnostic_gvs.keys().find(|class| {
            !toric_missing_gv_classes
                .iter()
                .any(|missing| missing == *class)
        }) {
            return Err(format!(
                "corrected-chamber formula-sum diagnostic GV produced a class that was not a toric miss: {:?}",
                sparse_i64(unexpected)
            ));
        }
        let (zero_count, nonzero_count) = diagnostic_gv_value_counts(&formula_sum_diagnostic_gvs);
        formula_sum_diagnostic_gv_zero_count = Some(zero_count);
        formula_sum_diagnostic_gv_nonzero_count = Some(nonzero_count);
        let mut formula_sum_small_curve_gvs = toric_small_curve_gvs.clone();
        let mut covered_count = 0usize;
        for missing_class in &toric_missing_gv_classes {
            if let Some(gv) = formula_sum_diagnostic_gvs.get(missing_class) {
                covered_count += 1;
                formula_sum_small_curve_gvs.push((missing_class.clone(), gv.clone()));
            }
        }
        let missing_count = toric_missing_gv_classes.len().saturating_sub(covered_count);
        formula_sum_diagnostic_gv_covered_count = Some(covered_count);
        formula_sum_diagnostic_gv_missing_count = Some(missing_count);
        if covered_count > 0 {
            formula_sum_diagnostic_gv_volume_correction = Some(
                cyrus_core::kklt::compute_gv_volume_correction_for_ambient_curves(
                    &formula_sum_small_curve_gvs,
                    &intersection.basis,
                    kahler,
                    Some(gamma),
                )
                .ok_or_else(|| {
                    "failed to compute formula-sum diagnostic corrected-chamber GV volume correction"
                        .to_string()
                })?,
            );
            formula_sum_diagnostic_gv_target_correction = Some(
                cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                    &formula_sum_small_curve_gvs,
                    &intersection.basis,
                    kklt_basis,
                    kahler,
                    Some(gamma),
                )
                .ok_or_else(|| {
                    "failed to compute formula-sum diagnostic corrected-chamber GV target correction"
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
        formula_sum_diagnostic_gv_covered_count,
        formula_sum_diagnostic_gv_missing_count,
        formula_sum_diagnostic_gv_zero_count,
        formula_sum_diagnostic_gv_nonzero_count,
        formula_sum_diagnostic_gv_volume_correction,
        formula_sum_diagnostic_gv_target_correction,
        remaining_gv_missing_count: missing_gv_classes.len(),
        first_missing_class: missing_gv_classes.first().cloned(),
        missing_required_degree_min,
        missing_required_degree_max,
        missing_target_stats,
        uncovered_source_ray_stats_degree_bound_for_missing,
        uncovered_source_ray_stats_for_missing,
        shared_facet_unresolved_source_ray_stats_for_missing,
        uncovered_source_ray_toric_diagnostic_sample,
        basis_mori_rays_for_missing_degree_bound: basis_rays_for_missing_degree_bound,
        basis_mori_rays_for_missing_degree_bounded: basis_rays_for_missing_degree_bounded,
        degree_bounded_mori_ray_context_for_missing,
        covered_toric_gv_context_for_missing,
        degree_bounded_toric_gv_diagnostic_context_for_missing,
        secondary_cone_height_certificate,
        secondary_cone_2face_height_certificate,
        expanded_secondary_fan_height_certificate,
        secondary_cone_heights_for_missing,
        gv_q_matrix_for_missing: gv_basis_data_for_missing
            .as_ref()
            .map(|data| data.q_matrix.clone()),
        gv_curve_basis_matrix_for_missing: gv_basis_data_for_missing.as_ref().map(|data| {
            data.curve_basis_matrix
                .iter()
                .map(|row| row.iter().map(ToString::to_string).collect())
                .collect()
        }),
        grading_for_missing,
        corrected_kappa_basis_for_missing,
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
    small_curve_pruning: CurvePruningStrategy,
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
    let small_curves =
        prune_selected_curve_candidates(&small_curve_candidates, small_curve_pruning, "chamber")?;
    let toric_gvs =
        compute_toric_two_face_curve_gv_invariants(tri, &geom.triangulation_points, &geom.polytope)
            .map_err(|e| format!("failed to compute chamber toric GV values: {e}"))?;
    let gv_by_class: HashMap<Vec<i64>, malachite::Integer> = toric_gvs
        .into_iter()
        .map(|item| (item.class, item.gv))
        .collect();

    let mut subcutoff_curve_gvs = Vec::with_capacity(small_curve_candidates.len());
    let mut subcutoff_missing_gv_classes = Vec::new();
    for curve in &small_curve_candidates {
        match gv_by_class.get(&curve.class) {
            Some(gv) => subcutoff_curve_gvs.push((curve.class.clone(), gv.clone())),
            None => subcutoff_missing_gv_classes.push(curve.class.clone()),
        }
    }

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
        subcutoff_toric_gv_covered_count: subcutoff_curve_gvs.len(),
        subcutoff_toric_gv_missing_count: subcutoff_missing_gv_classes.len(),
        toric_gv_covered_count: small_curve_gvs.len(),
        toric_gv_missing_count: missing_gv_classes.len(),
        first_missing_class: missing_gv_classes.first().cloned(),
        subcutoff_missing_gv_classes,
        missing_gv_classes,
        small_curve_candidates,
        small_curves,
        subcutoff_curve_gvs,
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
    small_curve_pruning: CurvePruningStrategy,
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
        let selection = compute_chamber_toric_gv_selection(
            &tri,
            geom,
            intersection,
            &current_t,
            cutoff,
            small_curve_pruning,
        )?;
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
            "[INFO] chamber-updated KKLT diagnostic iteration {iter}: changed_from_input={} ambient_rays={} subcutoff={} pruned={} pruning={} toric_covered={} toric_missing={} max_relative_step={} rel_err={} V_classical={} GV_volume_correction={}",
            chamber_changed_from_input,
            selection.ambient_rays,
            selection.subcutoff_count,
            selection.filtered_count,
            small_curve_pruning.as_str(),
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

fn compare_gv_target_correction_to_corrected_target_checkpoint(
    label: &str,
    data_dir: Option<&str>,
    kklt_basis: &[usize],
    base_target_tau: &[F64<Pos>],
    candidate_gv_correction: &[F64<Finite>],
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
    if checkpoint.len() != base_target_tau.len()
        || checkpoint.len() != candidate_gv_correction.len()
    {
        eprintln!(
            "[COMPARE] corrected target-volume diagnostic GV correction length mismatch ({label}): checkpoint={} base={} candidate={}",
            checkpoint.len(),
            base_target_tau.len(),
            candidate_gv_correction.len()
        );
        return;
    }
    let implied_gv_correction = base_target_tau
        .iter()
        .zip(checkpoint.iter())
        .map(|(base, target)| {
            F64::<Finite>::new(base.get() - target.get()).expect("implied GV correction is finite")
        })
        .collect::<Vec<_>>();
    let summary = target_correction_delta_summary(&implied_gv_correction, candidate_gv_correction)
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare corrected target-volume diagnostic GV correction ({label}): {e}"
            );
            std::process::exit(2);
        });
    eprintln!(
        "[COMPARE] corrected target-volume implied GV correction delta ({label}): max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_candidate={}",
        summary.max_abs_delta,
        summary.relative_l2_delta,
        summary.max_abs_reference,
        summary.max_abs_candidate
    );
    let mut deltas = implied_gv_correction
        .iter()
        .zip(candidate_gv_correction.iter())
        .enumerate()
        .map(|(idx, (implied, candidate))| {
            let delta = candidate.get() - implied.get();
            (idx, delta.abs(), delta, implied.get(), candidate.get())
        })
        .collect::<Vec<_>>();
    deltas.sort_unstable_by(|lhs, rhs| rhs.1.total_cmp(&lhs.1));
    for (idx, _abs_delta, delta, implied, candidate) in deltas.into_iter().take(8) {
        let point_idx = kklt_basis.get(idx).copied();
        eprintln!(
            "[COMPARE] corrected target-volume implied GV correction top_delta ({label}) kklt_idx={} point_idx={:?} delta={} checkpoint_implied={} candidate={}",
            idx, point_idx, delta, implied, candidate
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

fn corrected_divisor_values_from_parts(
    classical_tau: &[F64<Finite>],
    chi_divisor: &[I64<Finite>],
    gv_correction: &[F64<Finite>],
) -> Option<Vec<F64<Finite>>> {
    if classical_tau.len() != chi_divisor.len() || classical_tau.len() != gv_correction.len() {
        return None;
    }
    let twenty_four = F64::<Pos>::new(24.0).expect("24 is positive");
    Some(
        classical_tau
            .iter()
            .zip(chi_divisor.iter())
            .zip(gv_correction.iter())
            .map(|((classical, chi), gv)| *classical - chi.to_f64() / twenty_four + *gv)
            .collect(),
    )
}

fn corrected_chamber_gv_trace_json_path() -> Option<PathBuf> {
    std::env::var_os("CYRUS_CORRECTED_CHAMBER_GV_TRACE_JSON").map(PathBuf::from)
}

#[allow(clippy::too_many_arguments)]
fn write_corrected_chamber_gv_trace_json(
    path: &Path,
    basis: &[usize],
    kklt_basis: &[usize],
    checkpoint_t: &[F64<Finite>],
    gamma: &[I64<Finite>],
    c_i: &[I64<Pos>],
    c_tau: F64<Pos>,
    corrected_target_volumes: &[F64<Finite>],
    raw_kklt_target: &[F64<Finite>],
    checkpoint_chamber_chi: &[I64<Finite>],
    checkpoint_chamber_classical_tau: &[F64<Finite>],
    checkpoint_implied_gv: &[F64<Finite>],
    covered_gv_target: &[F64<Finite>],
    selection: &ChamberToricGvSelection,
    small_curve_pruning: CurvePruningStrategy,
) -> Result<(), String> {
    #[derive(Serialize)]
    struct TraceCurve {
        class: Vec<i64>,
        gv: String,
        q_dot_t: f64,
        parity_mod2: i8,
        target_correction_nonzero: Vec<(usize, f64)>,
        target_correction_l2: f64,
    }

    #[derive(Serialize)]
    struct TraceMissingCurve {
        class: Vec<i64>,
        q_dot_t: f64,
        parity_mod2: i8,
    }

    #[derive(Serialize)]
    struct Trace {
        basis: Vec<usize>,
        kklt_basis: Vec<usize>,
        checkpoint_t: Vec<f64>,
        gamma: Vec<i64>,
        c_i: Vec<i64>,
        c_tau: f64,
        corrected_target_volumes: Vec<f64>,
        raw_kklt_target: Vec<f64>,
        checkpoint_chamber_chi: Vec<i64>,
        checkpoint_chamber_classical_tau: Vec<f64>,
        checkpoint_implied_gv: Vec<f64>,
        toric_covered_gv_target: Vec<f64>,
        ambient_rays: usize,
        subcutoff_count: usize,
        small_curve_pruning: &'static str,
        pruned_count: usize,
        subcutoff_toric_covered_count: usize,
        subcutoff_toric_missing_count: usize,
        pruned_toric_covered_count: usize,
        pruned_toric_missing_count: usize,
        subcutoff_toric_curves: Vec<TraceCurve>,
        pruned_toric_curves: Vec<TraceCurve>,
        subcutoff_toric_missing_curves: Vec<TraceMissingCurve>,
        pruned_toric_missing_curves: Vec<TraceMissingCurve>,
    }

    fn trace_curves(
        curves: &[(Vec<i64>, malachite::Integer)],
        basis: &[usize],
        kklt_basis: &[usize],
        t: &[F64<Finite>],
        gamma: &[I64<Finite>],
    ) -> Result<Vec<TraceCurve>, String> {
        curves
            .iter()
            .map(|(class, gv)| {
                let q_dot_t = basis
                    .iter()
                    .zip(t.iter())
                    .map(|(&idx, ti)| class[idx] as f64 * ti.get())
                    .sum::<f64>();
                let parity = ambient_curve_b_field_parity_diagnostic(class, basis, gamma)
                    .ok_or_else(|| "failed to compute trace curve B-field parity".to_string())?;
                let single = [(class.clone(), gv.clone())];
                let contribution =
                    cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                        &single,
                        basis,
                        kklt_basis,
                        t,
                        Some(gamma),
                    )
                    .ok_or_else(|| {
                        format!(
                            "failed to compute trace target correction for q_dot_t={q_dot_t} parity_mod2={}",
                            parity.rem_euclid(2)
                        )
                    })?;
                let target_correction_nonzero = contribution
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, value)| {
                        let raw = value.get();
                        (raw != 0.0).then_some((idx, raw))
                    })
                    .collect::<Vec<_>>();
                let target_correction_l2 = contribution
                    .iter()
                    .map(|value| value.get() * value.get())
                    .sum::<f64>()
                    .sqrt();
                Ok(TraceCurve {
                    class: class.clone(),
                    gv: gv.to_string(),
                    q_dot_t,
                    parity_mod2: i8::try_from(parity.rem_euclid(2)).expect("mod-2 parity fits i8"),
                    target_correction_nonzero,
                    target_correction_l2,
                })
            })
            .collect()
    }

    fn trace_missing_curves(
        curves: &[Vec<i64>],
        basis: &[usize],
        t: &[F64<Finite>],
        gamma: &[I64<Finite>],
    ) -> Result<Vec<TraceMissingCurve>, String> {
        curves
            .iter()
            .map(|class| {
                let q_dot_t = basis
                    .iter()
                    .zip(t.iter())
                    .map(|(&idx, ti)| class[idx] as f64 * ti.get())
                    .sum::<f64>();
                let parity = ambient_curve_b_field_parity_diagnostic(class, basis, gamma)
                    .ok_or_else(|| "failed to compute trace curve B-field parity".to_string())?;
                Ok(TraceMissingCurve {
                    class: class.clone(),
                    q_dot_t,
                    parity_mod2: i8::try_from(parity.rem_euclid(2)).expect("mod-2 parity fits i8"),
                })
            })
            .collect()
    }

    let trace = Trace {
        basis: basis.to_vec(),
        kklt_basis: kklt_basis.to_vec(),
        checkpoint_t: finite_values(checkpoint_t),
        gamma: gamma.iter().map(|entry| entry.get()).collect(),
        c_i: c_i.iter().map(|entry| entry.get()).collect(),
        c_tau: c_tau.get(),
        corrected_target_volumes: finite_values(corrected_target_volumes),
        raw_kklt_target: finite_values(raw_kklt_target),
        checkpoint_chamber_chi: checkpoint_chamber_chi
            .iter()
            .map(|entry| entry.get())
            .collect(),
        checkpoint_chamber_classical_tau: finite_values(checkpoint_chamber_classical_tau),
        checkpoint_implied_gv: finite_values(checkpoint_implied_gv),
        toric_covered_gv_target: finite_values(covered_gv_target),
        ambient_rays: selection.ambient_rays,
        subcutoff_count: selection.subcutoff_count,
        small_curve_pruning: small_curve_pruning.as_str(),
        pruned_count: selection.filtered_count,
        subcutoff_toric_covered_count: selection.subcutoff_toric_gv_covered_count,
        subcutoff_toric_missing_count: selection.subcutoff_toric_gv_missing_count,
        pruned_toric_covered_count: selection.toric_gv_covered_count,
        pruned_toric_missing_count: selection.toric_gv_missing_count,
        subcutoff_toric_curves: trace_curves(
            &selection.subcutoff_curve_gvs,
            basis,
            kklt_basis,
            checkpoint_t,
            gamma,
        )?,
        pruned_toric_curves: trace_curves(
            &selection.small_curve_gvs,
            basis,
            kklt_basis,
            checkpoint_t,
            gamma,
        )?,
        subcutoff_toric_missing_curves: trace_missing_curves(
            &selection.subcutoff_missing_gv_classes,
            basis,
            checkpoint_t,
            gamma,
        )?,
        pruned_toric_missing_curves: trace_missing_curves(
            &selection.missing_gv_classes,
            basis,
            checkpoint_t,
            gamma,
        )?,
    };

    let content = serde_json::to_string_pretty(&trace)
        .map_err(|e| format!("failed to serialize corrected-chamber GV trace: {e}"))?;
    std::fs::write(path, content)
        .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
    Ok(())
}

fn write_corrected_chamber_gv_context_export(
    path: &Path,
    diag: &ChamberGvDiagnostic,
    small_curve_cutoff: F64<Pos>,
    small_curve_pruning: CurvePruningStrategy,
    missing_target_sample_limit: usize,
) -> Result<(), String> {
    let missing_target_sample_is_complete = diag
        .missing_target_stats
        .as_ref()
        .map(|stats| stats.sample.len() == stats.target_count);
    let export = CorrectedChamberGvContextExport {
        schema_version: 4,
        source: "mcallister_first_principles corrected-chamber GV diagnostic",
        small_curve_cutoff: small_curve_cutoff.get(),
        small_curve_pruning: small_curve_pruning.as_str(),
        ambient_rays: diag.ambient_rays,
        subcutoff_count: diag.subcutoff_count,
        filtered_count: diag.filtered_count,
        toric_gv_covered_count: diag.toric_gv_covered_count,
        toric_gv_missing_count: diag.toric_gv_missing_count,
        remaining_gv_missing_count: diag.remaining_gv_missing_count,
        first_missing_class: diag.first_missing_class.as_ref(),
        missing_required_degree_min: diag.missing_required_degree_min,
        missing_required_degree_max: diag.missing_required_degree_max,
        basis_mori_ray_count: diag.basis_mori_ray_count,
        degree_bounded_basis_mori_ray_count: diag.degree_bounded_basis_mori_ray_count,
        basis_mori_ray_degree_min: diag.basis_mori_ray_degree_min,
        basis_mori_ray_degree_max: diag.basis_mori_ray_degree_max,
        basis_mori_rays_for_missing_degree_bound: diag.basis_mori_rays_for_missing_degree_bound,
        basis_mori_rays_for_missing_degree_bounded: diag
            .basis_mori_rays_for_missing_degree_bounded
            .as_ref(),
        degree_bounded_mori_ray_context_for_missing: diag
            .degree_bounded_mori_ray_context_for_missing
            .as_ref(),
        covered_toric_gv_context_for_missing: diag.covered_toric_gv_context_for_missing.as_ref(),
        degree_bounded_toric_gv_diagnostic_context_for_missing: diag
            .degree_bounded_toric_gv_diagnostic_context_for_missing
            .as_ref(),
        secondary_cone_height_certificate: diag.secondary_cone_height_certificate.as_ref(),
        secondary_cone_2face_height_certificate: diag
            .secondary_cone_2face_height_certificate
            .as_ref(),
        expanded_secondary_fan_height_certificate: diag
            .expanded_secondary_fan_height_certificate
            .as_ref(),
        secondary_cone_heights_for_missing: diag.secondary_cone_heights_for_missing.as_ref(),
        gv_q_matrix_for_missing: diag.gv_q_matrix_for_missing.as_ref(),
        gv_curve_basis_matrix_for_missing: diag.gv_curve_basis_matrix_for_missing.as_ref(),
        grading_for_missing: diag.grading_for_missing.as_ref(),
        corrected_kappa_basis_for_missing: diag.corrected_kappa_basis_for_missing.as_ref(),
        missing_target_sample_limit,
        missing_target_sample_is_complete,
        missing_target_stats: diag.missing_target_stats.as_ref(),
        formula_sum_diagnostic_gv_covered_count: diag.formula_sum_diagnostic_gv_covered_count,
        formula_sum_diagnostic_gv_missing_count: diag.formula_sum_diagnostic_gv_missing_count,
        formula_sum_diagnostic_gv_zero_count: diag.formula_sum_diagnostic_gv_zero_count,
        formula_sum_diagnostic_gv_nonzero_count: diag.formula_sum_diagnostic_gv_nonzero_count,
        formula_sum_diagnostic_gv_volume_correction: diag
            .formula_sum_diagnostic_gv_volume_correction
            .as_ref()
            .map(|value| value.get()),
        formula_sum_diagnostic_gv_target_correction: diag
            .formula_sum_diagnostic_gv_target_correction
            .as_ref()
            .map(|correction| correction.iter().map(|value| value.get()).collect()),
        uncovered_source_ray_stats_degree_bound_for_missing: diag
            .uncovered_source_ray_stats_degree_bound_for_missing,
        uncovered_source_ray_stats_for_missing: diag
            .uncovered_source_ray_stats_for_missing
            .as_ref(),
        shared_facet_unresolved_source_ray_stats_for_missing: diag
            .shared_facet_unresolved_source_ray_stats_for_missing
            .as_ref(),
        uncovered_source_ray_toric_diagnostic_sample: diag
            .uncovered_source_ray_toric_diagnostic_sample
            .as_ref(),
    };
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create {}: {e}", parent.display()))?;
    }
    let content = serde_json::to_string_pretty(&export)
        .map_err(|e| format!("failed to serialize corrected-chamber GV context: {e}"))?;
    std::fs::write(path, content)
        .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
    Ok(())
}

fn write_secondary_cone_height_certificate(
    path: &Path,
    certificate: &SecondaryConeHeightCertificate,
) -> Result<(), String> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create {}: {e}", parent.display()))?;
    }
    let content = serde_json::to_string_pretty(certificate)
        .map_err(|e| format!("failed to serialize secondary-cone height certificate: {e}"))?;
    std::fs::write(path, content)
        .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
    Ok(())
}

fn write_corrected_chamber_face_triangulation_choice_summary(
    path: &Path,
    geom: &PrimalGeom,
    max_exact_face_points: usize,
    samples_per_large_face: usize,
    max_sampling_attempts_per_face: usize,
    seed: u64,
) -> Result<FaceTriangulationChoiceSummary, String> {
    let summary = corrected_chamber_face_triangulation_choice_summary(
        geom,
        max_exact_face_points,
        samples_per_large_face,
        max_sampling_attempts_per_face,
        seed,
    )?;
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create {}: {e}", parent.display()))?;
    }
    let content = serde_json::to_string_pretty(&summary).map_err(|e| {
        format!("failed to serialize corrected-chamber face triangulation choice summary: {e}")
    })?;
    std::fs::write(path, content)
        .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
    Ok(summary)
}

#[allow(clippy::too_many_arguments)]
fn compare_checkpoint_t_corrected_chamber_gv_target(
    data_dir: Option<&str>,
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    c_tau: F64<Pos>,
    base_target_tau: &[F64<Pos>],
    checkpoint_t: &[F64<Finite>],
    gamma: &[I64<Finite>],
    cutoff: F64<Pos>,
    small_curve_pruning: CurvePruningStrategy,
) {
    let Some(dir) = data_dir.map(PathBuf::from) else {
        return;
    };
    let target_path = dir.join("corrected_target_volumes.dat");
    if !target_path.exists() {
        return;
    }
    let checkpoint_target = read_csv_f64(&target_path)
        .into_iter()
        .map(|value| {
            F64::<Finite>::new(value).expect("corrected target-volume checkpoint is finite")
        })
        .collect::<Vec<_>>();
    if checkpoint_target.len() != base_target_tau.len()
        || checkpoint_target.len() != kklt_basis.len()
    {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target length mismatch: checkpoint={} base={} kklt={}",
            checkpoint_target.len(),
            base_target_tau.len(),
            kklt_basis.len()
        );
        return;
    }
    let checkpoint_implied_gv = base_target_tau
        .iter()
        .zip(checkpoint_target.iter())
        .map(|(base, target)| {
            F64::<Finite>::new(base.get() - target.get())
                .expect("checkpoint implied GV correction is finite")
        })
        .collect::<Vec<_>>();
    let checkpoint_chamber =
        triangulation_from_kahler_point(geom, &intersection.basis, checkpoint_t).unwrap_or_else(
            |e| {
                eprintln!("[ERROR] failed to build checkpoint-t corrected chamber: {e}");
                std::process::exit(2);
            },
        );
    let checkpoint_chamber_kappa_full = chamber_intersection_full(
        &checkpoint_chamber,
        &geom.triangulation_points,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute checkpoint-t corrected-chamber intersections: {e}");
        std::process::exit(2);
    });
    let checkpoint_chamber_kappa_basis =
        intersection_in_basis(&checkpoint_chamber_kappa_full, &intersection.basis);
    let Some(checkpoint_chamber_classical_tau) = cyrus_core::kklt::compute_kklt_divisor_volumes(
        &checkpoint_chamber_kappa_basis,
        &checkpoint_chamber_kappa_full,
        &intersection.basis,
        kklt_basis,
        checkpoint_t,
    ) else {
        eprintln!("[ERROR] failed to compute checkpoint-t corrected-chamber classical tau");
        std::process::exit(2);
    };
    let checkpoint_chamber_chi = cyrus_core::compute_kklt_divisor_chi(
        &geom.polytope,
        &geom.triangulation_points,
        &checkpoint_chamber_kappa_full,
        kklt_basis,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute checkpoint-t corrected-chamber divisor chi: {e}");
        std::process::exit(2);
    });
    let Some(checkpoint_chamber_base_target) =
        cyrus_core::kklt::compute_corrected_target_tau(c_i, &checkpoint_chamber_chi, c_tau)
    else {
        eprintln!("[ERROR] failed to compute checkpoint-t corrected-chamber base target tau");
        std::process::exit(2);
    };
    let checkpoint_chamber_implied_gv = checkpoint_chamber_base_target
        .iter()
        .zip(checkpoint_target.iter())
        .map(|(base, target)| {
            F64::<Finite>::new(base.get() - target.get())
                .expect("checkpoint corrected-chamber implied GV correction is finite")
        })
        .collect::<Vec<_>>();
    let chi_shift_summary = target_correction_delta_summary(
        &checkpoint_implied_gv,
        &checkpoint_chamber_implied_gv,
    )
    .unwrap_or_else(|e| {
        eprintln!(
            "[ERROR] failed to compare input-vs-corrected-chamber implied GV target correction: {e}"
        );
        std::process::exit(2);
    });
    eprintln!(
        "[COMPARE] checkpoint-t corrected-chamber implied GV correction chi-shift delta: max_abs={} relative_l2={} max_abs_input_chi_implied={} max_abs_corrected_chi_implied={}",
        chi_shift_summary.max_abs_delta,
        chi_shift_summary.relative_l2_delta,
        chi_shift_summary.max_abs_reference,
        chi_shift_summary.max_abs_candidate
    );
    let corrected_heights_path = dir.join("corrected_heights.dat");
    let corrected_heights_chamber = if corrected_heights_path.exists() {
        let corrected_heights = read_csv_f64(&corrected_heights_path);
        if corrected_heights.len() != geom.triangulation_points.len() {
            eprintln!(
                "[COMPARE] corrected_heights.dat length mismatch: heights={} points={}",
                corrected_heights.len(),
                geom.triangulation_points.len()
            );
            None
        } else {
            let height_chamber =
                compute_regular_triangulation(&geom.triangulation_points, &corrected_heights)
                    .unwrap_or_else(|e| {
                        eprintln!("[ERROR] failed to build corrected_heights.dat chamber: {e}");
                        std::process::exit(2);
                    });
            let same_as_checkpoint_t =
                triangulations_have_same_simplices(&checkpoint_chamber, &height_chamber);
            eprintln!(
                "[COMPARE] corrected_heights.dat chamber: simplices={} same_as_checkpoint_t={}",
                height_chamber.simplices().len(),
                same_as_checkpoint_t
            );
            (!same_as_checkpoint_t).then_some(height_chamber)
        }
    } else {
        None
    };
    let selection = compute_chamber_toric_gv_selection(
        &checkpoint_chamber,
        geom,
        intersection,
        checkpoint_t,
        cutoff,
        small_curve_pruning,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compute checkpoint-t corrected-chamber GV selection: {e}");
        std::process::exit(2);
    });
    let Some(covered_gv_target) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
        &selection.small_curve_gvs,
        &intersection.basis,
        kklt_basis,
        checkpoint_t,
        Some(gamma),
    ) else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber toric-covered GV target correction is invalid"
        );
        return;
    };
    let summary = target_correction_delta_summary(&checkpoint_implied_gv, &covered_gv_target)
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare checkpoint-t corrected-chamber GV target correction: {e}"
            );
            std::process::exit(2);
        });
    eprintln!(
        "[COMPARE] checkpoint-t corrected-chamber GV target correction delta: max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_toric_covered={} ambient_rays={} subcutoff={} pruned={} pruning={} toric_covered={} toric_missing={}",
        summary.max_abs_delta,
        summary.relative_l2_delta,
        summary.max_abs_reference,
        summary.max_abs_candidate,
        selection.ambient_rays,
        selection.subcutoff_count,
        selection.filtered_count,
        small_curve_pruning.as_str(),
        selection.toric_gv_covered_count,
        selection.toric_gv_missing_count
    );
    if let Some(subcutoff_target) =
        cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &selection.subcutoff_curve_gvs,
            &intersection.basis,
            kklt_basis,
            checkpoint_t,
            Some(gamma),
        )
    {
        let subcutoff_summary =
            target_correction_delta_summary(&checkpoint_implied_gv, &subcutoff_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare checkpoint-t corrected-chamber unpruned GV target correction: {e}"
                    );
                    std::process::exit(2);
                });
        let subcutoff_vs_pruned_summary =
            target_correction_delta_summary(&covered_gv_target, &subcutoff_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare unpruned vs pruned corrected-chamber GV target correction: {e}"
                    );
                    std::process::exit(2);
                });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (unpruned_subcutoff): max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_unpruned={} subcutoff_toric_covered={} subcutoff_toric_missing={} vs_pruned_max_abs={} vs_pruned_relative_l2={} pruning={}",
            subcutoff_summary.max_abs_delta,
            subcutoff_summary.relative_l2_delta,
            subcutoff_summary.max_abs_reference,
            subcutoff_summary.max_abs_candidate,
            selection.subcutoff_toric_gv_covered_count,
            selection.subcutoff_toric_gv_missing_count,
            subcutoff_vs_pruned_summary.max_abs_delta,
            subcutoff_vs_pruned_summary.relative_l2_delta,
            small_curve_pruning.as_str()
        );
    } else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber unpruned subcutoff GV target correction is invalid"
        );
    }
    let raw_kklt_target = c_i
        .iter()
        .map(|ci| {
            let value = (ci.to_f64() / c_tau).get();
            F64::<Finite>::new(value).expect("raw KKLT target is finite")
        })
        .collect::<Vec<_>>();
    if let Some(path) = corrected_chamber_gv_trace_json_path() {
        write_corrected_chamber_gv_trace_json(
            &path,
            &intersection.basis,
            kklt_basis,
            checkpoint_t,
            gamma,
            c_i,
            c_tau,
            &checkpoint_target,
            &raw_kklt_target,
            &checkpoint_chamber_chi,
            &checkpoint_chamber_classical_tau,
            &checkpoint_chamber_implied_gv,
            &covered_gv_target,
            &selection,
            small_curve_pruning,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] {e}");
            std::process::exit(2);
        });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV trace JSON written: {}",
            path.display()
        );
    }
    let checkpoint_equation_values = corrected_divisor_values_from_parts(
        &checkpoint_chamber_classical_tau,
        &checkpoint_chamber_chi,
        &checkpoint_chamber_implied_gv,
    )
    .unwrap_or_else(|| {
        eprintln!("[ERROR] failed to compute checkpoint-implied corrected divisor values");
        std::process::exit(2);
    });
    let checkpoint_equation_summary =
        target_correction_delta_summary(&raw_kklt_target, &checkpoint_equation_values)
            .unwrap_or_else(|e| {
                eprintln!(
                    "[ERROR] failed to compare checkpoint-implied corrected divisor equation: {e}"
                );
                std::process::exit(2);
            });
    let toric_equation_values = corrected_divisor_values_from_parts(
        &checkpoint_chamber_classical_tau,
        &checkpoint_chamber_chi,
        &covered_gv_target,
    )
    .unwrap_or_else(|| {
        eprintln!("[ERROR] failed to compute toric-covered corrected divisor values");
        std::process::exit(2);
    });
    let toric_equation_summary = target_correction_delta_summary(
        &raw_kklt_target,
        &toric_equation_values,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to compare toric-covered corrected divisor equation: {e}");
        std::process::exit(2);
    });
    eprintln!(
        "[COMPARE] checkpoint-t corrected divisor equation: checkpoint_implied_max_abs={} checkpoint_implied_relative_l2={} toric_covered_max_abs={} toric_covered_relative_l2={} max_abs_raw_target={} max_abs_toric_corrected_T={}",
        checkpoint_equation_summary.max_abs_delta,
        checkpoint_equation_summary.relative_l2_delta,
        toric_equation_summary.max_abs_delta,
        toric_equation_summary.relative_l2_delta,
        toric_equation_summary.max_abs_reference,
        toric_equation_summary.max_abs_candidate
    );
    let corrected_chi_summary =
        target_correction_delta_summary(&checkpoint_chamber_implied_gv, &covered_gv_target)
            .unwrap_or_else(|e| {
                eprintln!(
                    "[ERROR] failed to compare corrected-chi checkpoint-t corrected-chamber GV target correction: {e}"
                );
                std::process::exit(2);
            });
    eprintln!(
        "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (corrected_chamber_chi): max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_toric_covered={}",
        corrected_chi_summary.max_abs_delta,
        corrected_chi_summary.relative_l2_delta,
        corrected_chi_summary.max_abs_reference,
        corrected_chi_summary.max_abs_candidate
    );
    if let Some(no_gamma_target) = cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
        &selection.small_curve_gvs,
        &intersection.basis,
        kklt_basis,
        checkpoint_t,
        None,
    ) {
        let no_gamma_input_chi_summary =
            target_correction_delta_summary(&checkpoint_implied_gv, &no_gamma_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare checkpoint-t corrected-chamber no-gamma GV target correction: {e}"
                    );
                    std::process::exit(2);
                });
        let no_gamma_corrected_chi_summary =
            target_correction_delta_summary(&checkpoint_chamber_implied_gv, &no_gamma_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare corrected-chi checkpoint-t corrected-chamber no-gamma GV target correction: {e}"
                    );
                    std::process::exit(2);
                });
        let no_gamma_vs_gamma_summary =
            target_correction_delta_summary(&covered_gv_target, &no_gamma_target).unwrap_or_else(
                |e| {
                    eprintln!(
                        "[ERROR] failed to compare checkpoint-t corrected-chamber no-gamma-vs-gamma GV target correction: {e}"
                    );
                    std::process::exit(2);
                },
            );
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (no_gamma): input_chi_max_abs={} input_chi_relative_l2={} corrected_chi_max_abs={} corrected_chi_relative_l2={} vs_gamma_max_abs={} vs_gamma_relative_l2={}",
            no_gamma_input_chi_summary.max_abs_delta,
            no_gamma_input_chi_summary.relative_l2_delta,
            no_gamma_corrected_chi_summary.max_abs_delta,
            no_gamma_corrected_chi_summary.relative_l2_delta,
            no_gamma_vs_gamma_summary.max_abs_delta,
            no_gamma_vs_gamma_summary.relative_l2_delta
        );
    } else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber no-gamma GV target correction is invalid"
        );
    }
    let curve_basis_matrix = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compute curve-basis matrix for B-field gamma projection: {e}"
            );
            std::process::exit(2);
        });
    if let Some(basis_projected_gamma) =
        project_ambient_gamma_to_curve_basis(&curve_basis_matrix, gamma)
    {
        let parity_mismatches = count_gamma_parity_mismatches(
            &selection.small_curve_gvs,
            &intersection.basis,
            gamma,
            &basis_projected_gamma,
        )
        .unwrap_or_else(|| {
            eprintln!("[ERROR] failed to compare ambient and basis-projected gamma parities");
            std::process::exit(2);
        });
        let basis_projected_odd_count = basis_projected_gamma
            .iter()
            .filter(|entry| entry.get().rem_euclid(2) != 0)
            .count();
        let Some(basis_projected_target) =
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &selection.small_curve_gvs,
                &intersection.basis,
                kklt_basis,
                checkpoint_t,
                Some(&basis_projected_gamma),
            )
        else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (basis_projected_gamma) unavailable: correction is invalid"
            );
            return;
        };
        let basis_projected_summary =
            target_correction_delta_summary(&checkpoint_implied_gv, &basis_projected_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare checkpoint-t corrected-chamber basis-projected gamma correction: {e}"
                    );
                    std::process::exit(2);
                });
        let basis_projected_vs_ambient_summary = target_correction_delta_summary(
            &covered_gv_target,
            &basis_projected_target,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compare ambient vs basis-projected gamma corrections: {e}"
            );
            std::process::exit(2);
        });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (basis_projected_gamma): gamma_len={} odd_entries={} parity_mismatches_vs_ambient={} max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_basis_projected={} vs_ambient_max_abs={} vs_ambient_relative_l2={}",
            basis_projected_gamma.len(),
            basis_projected_odd_count,
            parity_mismatches,
            basis_projected_summary.max_abs_delta,
            basis_projected_summary.relative_l2_delta,
            basis_projected_summary.max_abs_reference,
            basis_projected_summary.max_abs_candidate,
            basis_projected_vs_ambient_summary.max_abs_delta,
            basis_projected_vs_ambient_summary.relative_l2_delta
        );
    } else {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (basis_projected_gamma) unavailable: failed to project ambient gamma"
        );
    }
    report_corrected_chamber_gv_branch_buckets(
        &selection.small_curve_gvs,
        &intersection.basis,
        kklt_basis,
        checkpoint_t,
        gamma,
        &checkpoint_implied_gv,
        &checkpoint_chamber_implied_gv,
        &covered_gv_target,
    );
    for threshold in [0.005_f64, 0.01, 0.02, 0.05, 0.1] {
        let thresholded_gvs = selection
            .small_curve_gvs
            .iter()
            .filter(|(curve, _invariant)| {
                intersection
                    .basis
                    .iter()
                    .zip(checkpoint_t.iter())
                    .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
                    .sum::<f64>()
                    >= threshold
            })
            .cloned()
            .collect::<Vec<_>>();
        let Some(thresholded_target) =
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &thresholded_gvs,
                &intersection.basis,
                kklt_basis,
                checkpoint_t,
                Some(gamma),
            )
        else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (drop_qdot_lt_{threshold}) unavailable: thresholded correction is invalid"
            );
            continue;
        };
        let input_chi_threshold_summary =
            target_correction_delta_summary(&checkpoint_implied_gv, &thresholded_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare input-chi thresholded corrected-chamber GV target correction: {e}"
                    );
                    std::process::exit(2);
                });
        let corrected_chi_threshold_summary =
            target_correction_delta_summary(&checkpoint_chamber_implied_gv, &thresholded_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare corrected-chi thresholded corrected-chamber GV target correction: {e}"
                    );
                    std::process::exit(2);
                });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (drop_qdot_lt_{threshold}): kept={} input_chi_max_abs={} input_chi_relative_l2={} corrected_chi_max_abs={} corrected_chi_relative_l2={}",
            thresholded_gvs.len(),
            input_chi_threshold_summary.max_abs_delta,
            input_chi_threshold_summary.relative_l2_delta,
            corrected_chi_threshold_summary.max_abs_delta,
            corrected_chi_threshold_summary.relative_l2_delta
        );
    }
    if let Some(height_chamber) = corrected_heights_chamber.as_ref() {
        let height_selection = compute_chamber_toric_gv_selection(
            height_chamber,
            geom,
            intersection,
            checkpoint_t,
            cutoff,
            small_curve_pruning,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute corrected_heights.dat GV selection: {e}");
            std::process::exit(2);
        });
        if let Some(height_target) =
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &height_selection.small_curve_gvs,
                &intersection.basis,
                kklt_basis,
                checkpoint_t,
                Some(gamma),
            )
        {
            let height_summary = target_correction_delta_summary(
                &checkpoint_implied_gv,
                &height_target,
            )
            .unwrap_or_else(|e| {
                eprintln!(
                    "[ERROR] failed to compare corrected_heights.dat GV target correction: {e}"
                );
                std::process::exit(2);
            });
            eprintln!(
                "[COMPARE] corrected_heights.dat GV target correction delta: max_abs={} relative_l2={} ambient_rays={} subcutoff={} pruned={} pruning={} toric_covered={} toric_missing={}",
                height_summary.max_abs_delta,
                height_summary.relative_l2_delta,
                height_selection.ambient_rays,
                height_selection.subcutoff_count,
                height_selection.filtered_count,
                small_curve_pruning.as_str(),
                height_selection.toric_gv_covered_count,
                height_selection.toric_gv_missing_count
            );
        } else {
            eprintln!(
                "[COMPARE] corrected_heights.dat GV target correction is invalid at checkpoint t"
            );
        }
    }
    for (label, offset) in [
        ("ambient_gamma_shift_minus_1", -1_isize),
        ("ambient_gamma_shift_plus_1", 1_isize),
    ] {
        let Some(shifted_gamma) = shifted_ambient_gamma(gamma, offset) else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV target correction delta ({label}) unavailable: shifted gamma would leave ambient range"
            );
            continue;
        };
        let Some(shifted_target) =
            cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                &selection.small_curve_gvs,
                &intersection.basis,
                kklt_basis,
                checkpoint_t,
                Some(&shifted_gamma),
            )
        else {
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV target correction delta ({label}) unavailable: shifted correction is invalid"
            );
            continue;
        };
        let shifted_summary =
            target_correction_delta_summary(&checkpoint_implied_gv, &shifted_target)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compare checkpoint-t corrected-chamber shifted GV target correction ({label}): {e}"
                    );
                    std::process::exit(2);
                });
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target correction delta ({label}): max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_shifted={}",
            shifted_summary.max_abs_delta,
            shifted_summary.relative_l2_delta,
            shifted_summary.max_abs_reference,
            shifted_summary.max_abs_candidate
        );
    }
    let origin_idx = geom
        .triangulation_points
        .iter()
        .position(|point| point.coords().iter().all(|&coord| coord == 0));
    if let Some(origin_idx) = origin_idx {
        if let Some(origin_toggled_gamma) = gamma_with_toggled_index(gamma, origin_idx) {
            let Some(origin_toggled_target) =
                cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                    &selection.small_curve_gvs,
                    &intersection.basis,
                    kklt_basis,
                    checkpoint_t,
                    Some(&origin_toggled_gamma),
                )
            else {
                eprintln!(
                    "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (ambient_gamma_toggle_origin) unavailable: toggled correction is invalid"
                );
                return;
            };
            let origin_toggled_summary =
                target_correction_delta_summary(&checkpoint_implied_gv, &origin_toggled_target)
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to compare checkpoint-t corrected-chamber origin-toggled GV target correction: {e}"
                        );
                        std::process::exit(2);
                    });
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV target correction delta (ambient_gamma_toggle_origin): origin_idx={} max_abs={} relative_l2={} max_abs_checkpoint_implied={} max_abs_toggled={}",
                origin_idx,
                origin_toggled_summary.max_abs_delta,
                origin_toggled_summary.relative_l2_delta,
                origin_toggled_summary.max_abs_reference,
                origin_toggled_summary.max_abs_candidate
            );
        }
    }
    let mut gv_deltas = checkpoint_implied_gv
        .iter()
        .zip(covered_gv_target.iter())
        .enumerate()
        .map(|(idx, (checkpoint_implied, toric_covered))| {
            let delta = toric_covered.get() - checkpoint_implied.get();
            (
                idx,
                delta.abs(),
                delta,
                checkpoint_implied.get(),
                toric_covered.get(),
            )
        })
        .collect::<Vec<_>>();
    gv_deltas.sort_unstable_by(|lhs, rhs| rhs.1.total_cmp(&lhs.1));
    for &(idx, _abs_delta, delta, checkpoint_implied, toric_covered) in gv_deltas.iter().take(8) {
        eprintln!(
            "[COMPARE] checkpoint-t corrected-chamber GV target top_delta kklt_idx={} point_idx={} delta={} checkpoint_implied={} toric_covered={} base_tau={} checkpoint_tau={}",
            idx,
            kklt_basis[idx],
            delta,
            checkpoint_implied,
            toric_covered,
            base_target_tau[idx].get(),
            checkpoint_target[idx].get()
        );
    }
    let decomposition_index = BoundedCurveDecompositionIndex::new(
        &selection.small_curve_candidates,
    )
    .unwrap_or_else(|e| {
        eprintln!(
            "[ERROR] failed to build checkpoint-t corrected-chamber decomposition index: {e}"
        );
        std::process::exit(2);
    });
    let small_curve_by_class = selection
        .small_curves
        .iter()
        .map(|candidate| (candidate.class.clone(), candidate))
        .collect::<HashMap<_, _>>();
    let toric_gv_diagnostic_by_class = compute_toric_curve_gv_diagnostics(
        &checkpoint_chamber,
        &geom.triangulation_points,
        &geom.polytope,
    )
    .unwrap_or_else(|e| {
        eprintln!(
            "[ERROR] failed to compute checkpoint-t corrected-chamber toric GV diagnostics: {e}"
        );
        std::process::exit(2);
    })
    .into_iter()
    .map(|diagnostic| (diagnostic.class.clone(), diagnostic))
    .collect::<HashMap<_, _>>();
    report_corrected_chamber_toric_gv_source_buckets(
        &selection.small_curve_gvs,
        &toric_gv_diagnostic_by_class,
        &intersection.basis,
        kklt_basis,
        checkpoint_t,
        gamma,
        &checkpoint_implied_gv,
        &checkpoint_chamber_implied_gv,
        &covered_gv_target,
    );
    if corrected_chamber_gv_weight_lp_requested() {
        report_corrected_chamber_gv_weight_lp(
            &selection.small_curve_gvs,
            &intersection.basis,
            kklt_basis,
            checkpoint_t,
            gamma,
            &checkpoint_implied_gv,
        );
    }
    let collect_top_toric_local_gv = corrected_chamber_top_toric_local_gv_requested();
    let mut top_toric_local_gv_targets = Vec::new();
    let mut top_toric_local_gv_seen = HashSet::new();
    let top_delta_limit = corrected_chamber_top_toric_local_gv_delta_limit();
    let top_contribution_limit = corrected_chamber_top_toric_local_gv_contribution_limit();
    for &(idx, _abs_delta, _delta, _checkpoint_implied, _toric_covered) in
        gv_deltas.iter().take(top_delta_limit)
    {
        let divisor_idx = kklt_basis[idx];
        let Some(mut rows) = ambient_target_contribution_rows(
            &selection.small_curve_gvs,
            &intersection.basis,
            divisor_idx,
            checkpoint_t,
            gamma,
        ) else {
            eprintln!(
                "[ERROR] failed to compute checkpoint-t corrected-chamber GV target contribution rows for kklt_idx={idx}"
            );
            continue;
        };
        rows.sort_unstable_by(|lhs, rhs| rhs.contribution.abs().total_cmp(&lhs.contribution.abs()));
        for row in rows.into_iter().take(top_contribution_limit) {
            let decomp_terms = small_curve_by_class
                .get(&row.class)
                .and_then(|candidate| {
                    decomposition_index
                        .find_decomposition(candidate, 4)
                        .ok()
                        .flatten()
                })
                .map(|decomposition| decomposition.len());
            let gv_source = toric_gv_diagnostic_by_class
                .get(&row.class)
                .map_or_else(|| "missing".to_string(), toric_curve_gv_diagnostic_summary);
            if collect_top_toric_local_gv && top_toric_local_gv_seen.insert(row.class.clone()) {
                top_toric_local_gv_targets.push(TopToricLocalGvTarget {
                    kklt_idx: idx,
                    point_idx: divisor_idx,
                    class: row.class.clone(),
                    toric_gv: row.gv.clone(),
                });
            }
            eprintln!(
                "[COMPARE] checkpoint-t corrected-chamber GV target contribution kklt_idx={} point_idx={} contribution={} q_i={} q_dot_t={} parity_mod2={} gv={} decomp_terms_le4={:?} gv_source={} class={:?}",
                idx,
                divisor_idx,
                row.contribution,
                row.q_divisor,
                row.q_dot_t,
                row.parity.rem_euclid(2),
                row.gv,
                decomp_terms,
                gv_source,
                sparse_i64(&row.class)
            );
        }
    }
    if collect_top_toric_local_gv {
        report_corrected_chamber_top_toric_local_gv_diagnostics(
            &checkpoint_chamber,
            geom,
            intersection,
            &top_toric_local_gv_targets,
        );
    }
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
    dump_corrected_chamber_gv_context_path: Option<&str>,
    dump_corrected_chamber_secondary_certificate_path: Option<&str>,
    dump_corrected_chamber_2face_secondary_certificate_path: Option<&str>,
    dump_corrected_chamber_expanded_secondary_fan_certificate_path: Option<&str>,
    dump_corrected_chamber_face_triangulation_choice_summary_path: Option<&str>,
    face_triangulation_max_exact_points: usize,
    face_triangulation_samples_per_large_face: usize,
    face_triangulation_max_sampling_attempts_per_face: usize,
    face_triangulation_seed: u64,
    diagnose_chamber_updated_kklt: bool,
    diagnose_chamber_updated_kklt_iterations: usize,
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
    if (diagnose_corrected_chamber_gv
        || diagnose_corrected_chamber_provided_generators_gv
        || diagnose_corrected_chamber_ray_gv
        || diagnose_corrected_chamber_lp_face_gv
        || dump_corrected_chamber_face_triangulation_choice_summary_path.is_some()
        || diagnose_chamber_updated_kklt)
        && allow_downstream_kahler
    {
        eprintln!(
            "[ERROR] corrected-chamber/chamber-updated diagnostics are only valid for first-principles runs, not downstream Kähler replay"
        );
        std::process::exit(2);
    }
    if dump_corrected_chamber_face_triangulation_choice_summary_path.is_some()
        && face_triangulation_max_sampling_attempts_per_face == 0
    {
        eprintln!(
            "[ERROR] --face-triangulation-max-sampling-attempts-per-face must be positive when dumping face triangulation choices"
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
        base_target_tau_for_chamber_gv,
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
        (t, None, None, None, None, None)
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
            (result, small_curve_selection_t, phase1.t.clone())
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
            (
                result,
                small_curve_selection_t,
                best_branch.result.t.clone(),
            )
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

        let mut correction_source_t = small_curve_selection_t;
        let mut correction_source_t_production = gv_iteration_source_t;
        let max_gv_iterations = 20usize;
        let gv_tolerance = 1e-10f64;
        let mut gv_converged = false;
        for iter in 0..max_gv_iterations {
            let previous_t = correction_source_t.clone();
            let previous_t_production = correction_source_t_production.clone();
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
            let Some(next) = solve_divisor_basis_path_following(
                &intersection.kappa_full,
                production_primal_basis.as_divisor_basis(),
                &kklt_basis,
                &tau_target,
                &previous_t_production,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!(
                    "[ERROR] GV-corrected divisor-basis KKLT solve failed at iteration {iter}"
                );
                std::process::exit(2);
            };
            if !next.converged {
                eprintln!(
                    "[ERROR] GV-corrected mixed-basis KKLT solve did not converge at iteration {iter}: rel_err={}",
                    next.relative_error.get()
                );
                std::process::exit(2);
            }
            let next_computed_t = transform_production_primal_kahler_to_computed(
                intersection,
                &production_primal_basis,
                &next.t,
                "GV-corrected Kähler point",
            );
            let max_relative_step = next_computed_t
                .iter()
                .zip(previous_t.iter())
                .map(|(new, old)| (new.get() - old.get()).abs() / (old.get().abs() + 1e-12))
                .fold(0.0f64, f64::max);
            eprintln!(
                "[INFO] GV-corrected KKLT iteration {iter}: max_relative_step={} rel_err={}",
                max_relative_step,
                next.relative_error.get()
            );
            correction_source_t_production = next.t.clone();
            correction_source_t = next_computed_t;
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
        if let Some(checkpoint_t) = checkpoint_t_for_gv.as_ref() {
            compare_checkpoint_t_corrected_chamber_gv_target(
                data_dir,
                geom,
                intersection,
                &kklt_basis,
                &c_i,
                c_tau,
                &tau_target,
                checkpoint_t,
                &gamma,
                small_curve_cutoff,
                small_curve_pruning,
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
                &correction_source_t,
                small_curve_cutoff,
                small_curve_pruning,
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
                .zip(correction_source_t.iter())
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
                classical_volume_from_t(&intersection.kappa_basis, &correction_source_t)
                    - bbhl.get()
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
            correction_source_t,
            Some(gv_volume_correction),
            Some(gamma),
            Some(kklt_basis),
            Some(tau_target),
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
    if let Some(path) = dump_corrected_chamber_secondary_certificate_path {
        let certificate = secondary_cone_height_certificate_for_kahler(
            &corrected_chamber,
            geom,
            &intersection.basis,
            &t,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to compute corrected-chamber secondary certificate: {e}");
            std::process::exit(2);
        });
        let path = PathBuf::from(path);
        write_secondary_cone_height_certificate(&path, &certificate).unwrap_or_else(|e| {
            eprintln!("[ERROR] failed to write corrected-chamber secondary certificate: {e}");
            std::process::exit(2);
        });
        eprintln!(
            "[INFO] corrected-chamber secondary certificate JSON written: {}",
            path.display()
        );
    }
    if let Some(path) = dump_corrected_chamber_2face_secondary_certificate_path {
        let certificate = secondary_cone_2face_height_certificate_for_kahler(
            &corrected_chamber,
            geom,
            &intersection.basis,
            &t,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to compute corrected-chamber 2-face secondary certificate: {e}"
            );
            std::process::exit(2);
        });
        let path = PathBuf::from(path);
        write_secondary_cone_height_certificate(&path, &certificate).unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to write corrected-chamber 2-face secondary certificate: {e}"
            );
            std::process::exit(2);
        });
        eprintln!(
            "[INFO] corrected-chamber 2-face secondary certificate JSON written: {}",
            path.display()
        );
    }
    if let Some(path) = dump_corrected_chamber_expanded_secondary_fan_certificate_path {
        let certificate =
            expanded_secondary_fan_height_certificate_for_kahler(geom, &intersection.basis, &t)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to compute corrected-chamber expanded-secondary fan certificate: {e}"
                    );
                    std::process::exit(2);
                });
        let path = PathBuf::from(path);
        write_secondary_cone_height_certificate(&path, &certificate).unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to write corrected-chamber expanded-secondary fan certificate: {e}"
            );
            std::process::exit(2);
        });
        eprintln!(
            "[INFO] corrected-chamber expanded-secondary fan certificate JSON written: {}",
            path.display()
        );
    }
    if let Some(path) = dump_corrected_chamber_face_triangulation_choice_summary_path {
        let path = PathBuf::from(path);
        let summary = write_corrected_chamber_face_triangulation_choice_summary(
            &path,
            geom,
            face_triangulation_max_exact_points,
            face_triangulation_samples_per_large_face,
            face_triangulation_max_sampling_attempts_per_face,
            face_triangulation_seed,
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "[ERROR] failed to write corrected-chamber face triangulation choice summary: {e}"
            );
            std::process::exit(2);
        });
        eprintln!(
            "[INFO] corrected-chamber face triangulation choice summary JSON written: {} status={} faces={} sampled_faces={} empty_choice_faces={} total_choices={}",
            path.display(),
            summary.status,
            summary.face_count,
            summary.sampled_face_count,
            summary.empty_choice_face_count,
            summary.total_choice_count
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
    let dump_corrected_chamber_gv_context = dump_corrected_chamber_gv_context_path.is_some();
    if diagnose_corrected_chamber_gv
        || diagnose_corrected_chamber_provided_generators_gv
        || diagnose_corrected_chamber_ray_gv
        || diagnose_corrected_chamber_lp_face_gv
        || dump_corrected_chamber_gv_context
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
        let missing_target_sample_limit = if dump_corrected_chamber_gv_context {
            usize::MAX
        } else {
            10
        };
        let diag = diagnose_chamber_gv_volume_correction(
            &corrected_chamber,
            geom,
            intersection,
            kklt_basis,
            &t,
            gamma,
            small_curve_cutoff,
            small_curve_pruning,
            primal_gv_min_points,
            primal_gv_max_deg,
            diagnose_corrected_chamber_provided_generators_gv,
            diagnose_corrected_chamber_ray_gv,
            diagnose_corrected_chamber_lp_face_gv,
            missing_target_sample_limit,
        )
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] corrected-chamber GV diagnostic failed: {e}");
            std::process::exit(2);
        });
        if let Some(path) = dump_corrected_chamber_gv_context_path {
            let path = PathBuf::from(path);
            write_corrected_chamber_gv_context_export(
                &path,
                &diag,
                small_curve_cutoff,
                small_curve_pruning,
                missing_target_sample_limit,
            )
            .unwrap_or_else(|e| {
                eprintln!("[ERROR] failed to write corrected-chamber GV context export: {e}");
                std::process::exit(2);
            });
            eprintln!(
                "[INFO] corrected-chamber GV context JSON written: {}",
                path.display()
            );
        }
        eprintln!(
            "[INFO] corrected-chamber small toric curves: ambient_rays={} subcutoff={} pruned={} pruning={} toric_covered={} toric_missing={} cutoff={}",
            diag.ambient_rays,
            diag.subcutoff_count,
            diag.filtered_count,
            small_curve_pruning.as_str(),
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
        if let Some(formula_sum_covered) = diag.formula_sum_diagnostic_gv_covered_count {
            let formula_sum_missing = diag.formula_sum_diagnostic_gv_missing_count.unwrap_or(0);
            eprintln!(
                "[INFO] corrected-chamber local formula-sum diagnostic GV covered {}/{} toric-missing curves; remaining={}",
                formula_sum_covered, diag.toric_gv_missing_count, formula_sum_missing
            );
            if let (Some(zero_count), Some(nonzero_count)) = (
                diag.formula_sum_diagnostic_gv_zero_count,
                diag.formula_sum_diagnostic_gv_nonzero_count,
            ) {
                eprintln!(
                    "[INFO] corrected-chamber local formula-sum diagnostic GV values: zero={} nonzero={}",
                    zero_count, nonzero_count
                );
            }
            if let Some(formula_sum_correction) =
                diag.formula_sum_diagnostic_gv_volume_correction.as_ref()
            {
                let scope = if formula_sum_missing == 0 {
                    ""
                } else {
                    "partial "
                };
                eprintln!(
                    "[INFO] corrected-chamber local formula-sum diagnostic GV {scope}volume correction (diagnostic, not promoted) = {}",
                    formula_sum_correction.get()
                );
                if let Some(input_chamber_correction) = gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber local formula-sum diagnostic GV volume correction delta_vs_input_chamber (diagnostic, not promoted) = {}",
                        formula_sum_correction.get() - input_chamber_correction.get()
                    );
                }
                if let Some(covered_correction) = diag.covered_gv_volume_correction.as_ref() {
                    eprintln!(
                        "[INFO] corrected-chamber local formula-sum diagnostic GV volume correction delta_vs_toric_covered (diagnostic, not promoted) = {}",
                        formula_sum_correction.get() - covered_correction.get()
                    );
                }
                if let (Some(input_target_correction), Some(formula_sum_target_correction)) = (
                    input_chamber_gv_target_correction.as_ref(),
                    diag.formula_sum_diagnostic_gv_target_correction.as_ref(),
                ) {
                    let summary = target_correction_delta_summary(
                        input_target_correction,
                        formula_sum_target_correction,
                    )
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "[ERROR] failed to compare formula-sum diagnostic GV target correction: {e}"
                        );
                        std::process::exit(2);
                    });
                    eprintln!(
                        "[INFO] corrected-chamber local formula-sum diagnostic GV target correction delta_vs_input_chamber (diagnostic, not promoted) = {:?}",
                        summary
                    );
                }
                if let (Some(base_target_tau), Some(formula_sum_target_correction)) = (
                    base_target_tau_for_chamber_gv.as_deref(),
                    diag.formula_sum_diagnostic_gv_target_correction.as_deref(),
                ) {
                    compare_gv_target_correction_to_corrected_target_checkpoint(
                        "corrected_chamber_formula_sum_diagnostic",
                        data_dir,
                        kklt_basis,
                        base_target_tau,
                        formula_sum_target_correction,
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
                "[INFO] corrected-chamber missing GV target reduction: targets={} targets_as_mori_generators={} targets_as_origin_circuits={} real_cone_decomposable_by_other_generators={} lp_extremal_mori_generators={} real_cone_active_generator_range={} exact_decomposition_kinds={:?} origin_circuit_resolved_conifold={} origin_circuit_affine_ranks={:?} generators_le_target_degree_range={}..{} origin_coefficients={:?}",
                stats.target_count,
                stats.targets_that_are_mori_generators,
                stats.targets_that_are_origin_circuits,
                stats.targets_real_cone_decomposable_by_other_generators,
                stats.targets_that_are_lp_extremal_mori_generators,
                active_generator_range,
                stats.real_cone_decomposition_exact_kind_counts,
                stats.origin_circuit_resolved_conifold_count,
                stats.origin_circuit_affine_rank_counts,
                stats.min_generators_le_target_degree,
                stats.max_generators_le_target_degree,
                stats.origin_coefficient_counts
            );
            eprintln!(
                "[INFO] corrected-chamber missing GV origin-circuit patterns: {:?}",
                stats.origin_circuit_pattern_counts
            );
            eprintln!(
                "[INFO] corrected-chamber missing GV branch statuses: {:?}",
                stats.branch_status_counts
            );
            eprintln!(
                "[INFO] corrected-chamber missing GV branch buckets: {:?}",
                stats.branch_bucket_counts
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

fn cms_general_divisor_formula_sum(
    candidates: Option<&[CmsGeneralDivisorShapeCandidate]>,
) -> Option<malachite::Integer> {
    let mut found = false;
    let mut sum = malachite::Integer::from(0);
    for value in candidates
        .unwrap_or_default()
        .iter()
        .filter_map(|candidate| candidate.toric_gv1_formula_value)
    {
        found = true;
        sum += malachite::Integer::from(value);
    }
    found.then_some(sum)
}

fn formula_sum_diagnostic_missing_gvs(
    missing_target_stats: Option<&MissingGvTargetStats>,
    ambient_dim: usize,
) -> Result<HashMap<Vec<i64>, malachite::Integer>, String> {
    let mut diagnostic_gvs = HashMap::new();
    let Some(stats) = missing_target_stats else {
        return Ok(diagnostic_gvs);
    };
    for (idx, sample) in stats.sample.iter().enumerate() {
        let Some(gv) =
            cms_general_divisor_formula_sum(sample.cms_general_divisor_shape_candidates.as_deref())
        else {
            continue;
        };
        let ambient_class = dense_i64_from_sparse(
            &sample.ambient_nonzero,
            ambient_dim,
            &format!("formula-sum diagnostic missing target {idx}"),
        )?;
        insert_missing_diagnostic_gv(
            &mut diagnostic_gvs,
            &ambient_class,
            &gv,
            "local formula-sum",
        )?;
    }
    Ok(diagnostic_gvs)
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
        args.diagnose_corrected_chamber_gv,
        args.diagnose_corrected_chamber_provided_generators_gv,
        args.diagnose_corrected_chamber_ray_gv,
        args.diagnose_corrected_chamber_lp_face_gv,
        args.dump_corrected_chamber_gv_context_path.as_deref(),
        args.dump_corrected_chamber_secondary_certificate_path
            .as_deref(),
        args.dump_corrected_chamber_2face_secondary_certificate_path
            .as_deref(),
        args.dump_corrected_chamber_expanded_secondary_fan_certificate_path
            .as_deref(),
        args.dump_corrected_chamber_face_triangulation_choice_summary_path
            .as_deref(),
        args.face_triangulation_max_exact_points,
        args.face_triangulation_samples_per_large_face,
        args.face_triangulation_max_sampling_attempts_per_face,
        args.face_triangulation_seed,
        args.diagnose_chamber_updated_kklt,
        args.diagnose_chamber_updated_kklt_iterations,
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
    fn corrected_chamber_trace_records_curve_pruning_strategy() {
        use cyrus_core::{f64_finite, f64_pos, i64_pos};

        let cache_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/cyrus-test-cache");
        std::fs::create_dir_all(&cache_dir).expect("create cache dir");
        let path = cache_dir.join(format!(
            "corrected-chamber-trace-{}-{}.json",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("time after epoch")
                .as_nanos()
        ));
        let selection = ChamberToricGvSelection {
            ambient_rays: 3,
            subcutoff_count: 2,
            filtered_count: 1,
            subcutoff_toric_gv_covered_count: 1,
            subcutoff_toric_gv_missing_count: 1,
            toric_gv_covered_count: 1,
            toric_gv_missing_count: 0,
            first_missing_class: None,
            subcutoff_missing_gv_classes: vec![vec![2]],
            missing_gv_classes: Vec::new(),
            small_curve_candidates: Vec::new(),
            small_curves: Vec::new(),
            subcutoff_curve_gvs: vec![(vec![1], malachite::Integer::from(1))],
            small_curve_gvs: vec![(vec![1], malachite::Integer::from(1))],
        };

        write_corrected_chamber_gv_trace_json(
            &path,
            &[0],
            &[0],
            &[f64_finite!(0.5)],
            &[I64::<Finite>::new(1)],
            &[i64_pos!(2)],
            f64_pos!(3.0),
            &[f64_finite!(4.0)],
            &[f64_finite!(5.0)],
            &[I64::<Finite>::new(6)],
            &[f64_finite!(7.0)],
            &[f64_finite!(8.0)],
            &[f64_finite!(9.0)],
            &selection,
            CurvePruningStrategy::FiniteSemigroup,
        )
        .expect("write corrected chamber trace");

        let content = std::fs::read_to_string(&path).expect("read corrected chamber trace");
        let value = serde_json::from_str::<serde_json::Value>(&content).expect("valid trace JSON");
        assert_eq!(value["small_curve_pruning"], "finite-semigroup");
        assert_eq!(value["pruned_count"], 1);
        assert_eq!(value["pruned_toric_covered_count"], 1);
        assert_eq!(value["pruned_toric_missing_count"], 0);
        assert!(value.get("pair_pruned_count").is_none());
        let pruned_curves = value["pruned_toric_curves"]
            .as_array()
            .expect("pruned toric curves should be exported");
        assert_eq!(pruned_curves.len(), 1);
        assert_eq!(pruned_curves[0]["target_correction_nonzero"][0][0], 0);
        assert!(
            pruned_curves[0]["target_correction_l2"]
                .as_f64()
                .expect("l2 is finite")
                > 0.0
        );
        assert!(value["pruned_toric_missing_curves"].is_array());

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn corrected_chamber_context_export_includes_uncovered_source_toric_diagnostics() {
        let cache_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/cyrus-test-cache");
        std::fs::create_dir_all(&cache_dir).expect("create cache dir");
        let path = cache_dir.join(format!(
            "corrected-chamber-context-{}-{}.json",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("time after epoch")
                .as_nanos()
        ));
        let diag = ChamberGvDiagnostic {
            ambient_rays: 3,
            subcutoff_count: 2,
            filtered_count: 1,
            toric_gv_covered_count: 0,
            toric_gv_missing_count: 1,
            basis_mori_ray_count: None,
            degree_bounded_basis_mori_ray_count: None,
            basis_mori_ray_degree_min: None,
            basis_mori_ray_degree_max: None,
            general_gv_covered_count: None,
            ray_gv_covered_count: None,
            ray_gv_volume_correction: None,
            ray_gv_sample: Vec::new(),
            lp_face_gv_covered_count: None,
            lp_face_gv_failed_count: None,
            lp_face_gv_certified_count: None,
            lp_face_gv_uncertified_count: None,
            lp_face_gv_volume_correction: None,
            lp_face_gv_sample: Vec::new(),
            combined_diagnostic_gv_covered_count: None,
            combined_diagnostic_gv_missing_count: None,
            combined_diagnostic_gv_zero_count: None,
            combined_diagnostic_gv_nonzero_count: None,
            combined_diagnostic_gv_volume_correction: None,
            combined_diagnostic_gv_target_correction: None,
            formula_sum_diagnostic_gv_covered_count: None,
            formula_sum_diagnostic_gv_missing_count: None,
            formula_sum_diagnostic_gv_zero_count: None,
            formula_sum_diagnostic_gv_nonzero_count: None,
            formula_sum_diagnostic_gv_volume_correction: None,
            formula_sum_diagnostic_gv_target_correction: None,
            remaining_gv_missing_count: 1,
            first_missing_class: Some(vec![0, 1, -1]),
            missing_required_degree_min: Some(4),
            missing_required_degree_max: Some(4),
            missing_target_stats: None,
            uncovered_source_ray_stats_degree_bound_for_missing: Some(4),
            uncovered_source_ray_stats_for_missing: None,
            shared_facet_unresolved_source_ray_stats_for_missing: None,
            uncovered_source_ray_toric_diagnostic_sample: Some(vec![
                ToricGvDiagnosticContextSample {
                    degree: 4,
                    gv: "-2".to_string(),
                    source_bucket: "origin_circuit".to_string(),
                    source_summary: "diagnostic".to_string(),
                    ambient_nonzero: vec![(1, 1), (2, -1)],
                    basis_nonzero: vec![(0, 1), (1, -1)],
                },
            ]),
            degree_bounded_toric_gv_diagnostic_context_for_missing: Some(vec![
                ToricGvDiagnosticContextSample {
                    degree: 2,
                    gv: "1".to_string(),
                    source_bucket: "two_face".to_string(),
                    source_summary: "degree-bounded diagnostic".to_string(),
                    ambient_nonzero: vec![(0, 1), (2, -1)],
                    basis_nonzero: vec![(0, 1)],
                },
            ]),
            secondary_cone_height_certificate: Some(SecondaryConeHeightCertificate {
                status: "strictly_inside_secondary_cone".to_string(),
                epsilon: 1e-6,
                hyperplane_count: 1,
                pairing_count: 1,
                min_pairing: Some(1.0),
                max_pairing: Some(1.0),
                strictly_inside: true,
            }),
            secondary_cone_2face_height_certificate: Some(SecondaryConeHeightCertificate {
                status: "strictly_inside_secondary_cone".to_string(),
                epsilon: 1e-6,
                hyperplane_count: 2,
                pairing_count: 2,
                min_pairing: Some(0.5),
                max_pairing: Some(1.5),
                strictly_inside: true,
            }),
            expanded_secondary_fan_height_certificate: None,
            secondary_cone_heights_for_missing: Some(vec![0.0, 0.0, 1.5]),
            basis_mori_rays_for_missing_degree_bound: None,
            basis_mori_rays_for_missing_degree_bounded: None,
            degree_bounded_mori_ray_context_for_missing: None,
            covered_toric_gv_context_for_missing: None,
            gv_q_matrix_for_missing: None,
            gv_curve_basis_matrix_for_missing: None,
            grading_for_missing: None,
            corrected_kappa_basis_for_missing: None,
            covered_toric_gv_divisor_representation_baseline: None,
            covered_gv_target_correction: None,
            covered_gv_volume_correction: None,
            gv_volume_correction: None,
        };

        write_corrected_chamber_gv_context_export(
            &path,
            &diag,
            F64::<Pos>::new(1.0).expect("positive cutoff"),
            CurvePruningStrategy::PairDecomposable,
            8,
        )
        .expect("write corrected chamber context");

        let content = std::fs::read_to_string(&path).expect("read corrected chamber context");
        let value =
            serde_json::from_str::<serde_json::Value>(&content).expect("valid context JSON");
        let sample = value["uncovered_source_ray_toric_diagnostic_sample"]
            .as_array()
            .expect("source toric diagnostic sample should be exported");
        assert_eq!(sample.len(), 1);
        assert_eq!(sample[0]["degree"], 4);
        assert_eq!(sample[0]["gv"], "-2");
        assert_eq!(sample[0]["source_bucket"], "origin_circuit");
        let degree_bounded_sample = value["degree_bounded_toric_gv_diagnostic_context_for_missing"]
            .as_array()
            .expect("degree-bounded toric diagnostic context should be exported");
        assert_eq!(degree_bounded_sample.len(), 1);
        assert_eq!(degree_bounded_sample[0]["gv"], "1");
        assert_eq!(degree_bounded_sample[0]["source_bucket"], "two_face");
        let chamber_certificate = &value["secondary_cone_height_certificate"];
        assert_eq!(
            chamber_certificate["status"],
            "strictly_inside_secondary_cone"
        );
        assert_eq!(chamber_certificate["hyperplane_count"], 1);
        assert_eq!(chamber_certificate["strictly_inside"], true);
        let twoface_chamber_certificate = &value["secondary_cone_2face_height_certificate"];
        assert_eq!(twoface_chamber_certificate["hyperplane_count"], 2);
        assert_eq!(twoface_chamber_certificate["strictly_inside"], true);
        assert_eq!(
            value["secondary_cone_heights_for_missing"]
                .as_array()
                .expect("secondary-cone heights should be exported")
                .iter()
                .map(|entry| entry.as_f64().expect("height should be numeric"))
                .collect::<Vec<_>>(),
            vec![0.0, 0.0, 1.5]
        );

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn secondary_cone_height_certificate_writer_serializes_certificate() {
        let cache_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/cyrus-test-cache");
        std::fs::create_dir_all(&cache_dir).expect("create cache dir");
        let path = cache_dir.join(format!(
            "secondary-cone-certificate-{}-{}.json",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("time after epoch")
                .as_nanos()
        ));
        let certificate = SecondaryConeHeightCertificate {
            status: "strictly_inside_secondary_cone".to_string(),
            epsilon: 1e-6,
            hyperplane_count: 2,
            pairing_count: 2,
            min_pairing: Some(0.25),
            max_pairing: Some(3.0),
            strictly_inside: true,
        };

        write_secondary_cone_height_certificate(&path, &certificate)
            .expect("write secondary-cone certificate");

        let content = std::fs::read_to_string(&path).expect("read secondary-cone certificate");
        let value =
            serde_json::from_str::<serde_json::Value>(&content).expect("valid certificate JSON");
        assert_eq!(value["status"], "strictly_inside_secondary_cone");
        assert_eq!(value["hyperplane_count"], 2);
        assert_eq!(value["strictly_inside"], true);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn corrected_chamber_face_triangulation_choice_summary_tracks_sampled_faces() {
        let vertices = vec![
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![0, 0, 0, 1]),
            Point::new(vec![-1, -1, -1, -1]),
        ];
        let mut triangulation_points = vec![Point::new(vec![0, 0, 0, 0])];
        triangulation_points.extend(vertices.clone());
        let heights = triangulation_points
            .iter()
            .map(|_| F64::<Finite>::new(0.0).expect("height is finite"))
            .collect::<Vec<_>>();
        let geom = PrimalGeom {
            polytope: Polytope::from_vertices(vertices).expect("build simplex polytope"),
            heights,
            triangulation_points,
            triangulation: Triangulation::new(vec![vec![1, 2, 3, 4, 5]]),
        };

        let summary = corrected_chamber_face_triangulation_choice_summary(&geom, 2, 1, 16, 0)
            .expect("summarize sampled simplex face choices");

        assert_eq!(
            summary.status,
            "face_triangulation_choices_exact_and_sampled"
        );
        assert_eq!(summary.face_count, 10);
        assert_eq!(summary.exact_face_count, 0);
        assert_eq!(summary.sampled_face_count, 10);
        assert_eq!(summary.empty_choice_face_count, 0);
        assert_eq!(summary.min_face_points, Some(3));
        assert_eq!(summary.max_face_points, Some(3));
        assert_eq!(summary.choice_counts, vec![1; 10]);
        assert_eq!(summary.total_choice_count, "1");
        assert_eq!(summary.sampled_face_point_counts, vec![3; 10]);
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
            &[],
            None,
            &basis_rays,
            &[1, 2],
            &[2, 3],
            None,
            None,
            0,
            &HashMap::new(),
            &HashMap::new(),
            true,
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
        assert_eq!(stats.origin_circuit_affine_rank_counts, BTreeMap::new());
        assert_eq!(stats.branch_status_counts, BTreeMap::new());
        assert_eq!(stats.branch_bucket_counts, BTreeMap::new());
        assert_eq!(
            stats.real_cone_decomposition_exact_kind_counts,
            BTreeMap::from([("integer_semigroup".to_string(), 2)])
        );
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
                    origin_circuit_witnesses: None,
                    origin_circuit_affine_support: None,
                    cms_general_divisor_shape_candidates: None,
                    cms_general_divisor_intersection_checks: None,
                    branch_diagnostic: None,
                    real_cone_decomposable_by_other_generators: true,
                    real_cone_decomposition_active_generators: Some(2),
                    real_cone_decomposition_active_generator_basis_nonzero: Some(vec![
                        vec![(0, 1)],
                        vec![(1, 1)]
                    ]),
                    real_cone_decomposition_exact_coefficients: Some(vec![
                        "1".to_string(),
                        "1".to_string()
                    ]),
                    real_cone_decomposition_exact_kind: Some("integer_semigroup"),
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
                    origin_circuit_witnesses: None,
                    origin_circuit_affine_support: None,
                    cms_general_divisor_shape_candidates: None,
                    cms_general_divisor_intersection_checks: None,
                    branch_diagnostic: None,
                    real_cone_decomposable_by_other_generators: true,
                    real_cone_decomposition_active_generators: Some(1),
                    real_cone_decomposition_active_generator_basis_nonzero: Some(vec![vec![(
                        0, 1
                    )]]),
                    real_cone_decomposition_exact_coefficients: Some(vec!["2".to_string()]),
                    real_cone_decomposition_exact_kind: Some("integer_semigroup"),
                    ambient_nonzero: vec![(1, 2)],
                    basis_nonzero: vec![(0, 2)]
                }
            ]
        );
    }

    #[test]
    fn missing_gv_target_stats_can_skip_real_cone_decomposition_for_source_rays() {
        let stats = missing_gv_target_stats(
            &[vec![0, 2, 0]],
            &[],
            None,
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &[1, 2],
            &[2, 3],
            None,
            None,
            0,
            &HashMap::new(),
            &HashMap::new(),
            false,
            1,
        )
        .unwrap();

        assert_eq!(stats.targets_real_cone_decomposable_by_other_generators, 0);
        assert_eq!(
            stats.sample[0].real_cone_decomposable_by_other_generators,
            false
        );
        assert_eq!(stats.sample[0].real_cone_decomposition_exact_kind, None);
    }

    #[test]
    fn active_dependency_source_rays_extend_uncovered_source_context() {
        let target_stats = MissingGvTargetStats {
            target_count: 2,
            targets_that_are_mori_generators: 0,
            targets_that_are_origin_circuits: 0,
            targets_real_cone_decomposable_by_other_generators: 1,
            targets_that_are_lp_extremal_mori_generators: 0,
            real_cone_decomposition_active_generator_min: Some(4),
            real_cone_decomposition_active_generator_max: Some(4),
            origin_circuit_resolved_conifold_count: 0,
            min_generators_le_target_degree: 0,
            max_generators_le_target_degree: 0,
            origin_coefficient_counts: BTreeMap::new(),
            origin_circuit_pattern_counts: BTreeMap::new(),
            origin_circuit_affine_rank_counts: BTreeMap::new(),
            branch_status_counts: BTreeMap::new(),
            branch_bucket_counts: BTreeMap::new(),
            real_cone_decomposition_exact_kind_counts: BTreeMap::new(),
            sample: vec![
                MissingGvTargetSample {
                    degree: 5,
                    generators_le_degree: 4,
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
                    real_cone_decomposition_active_generators: Some(4),
                    real_cone_decomposition_active_generator_basis_nonzero: Some(vec![
                        vec![(0, 1)],
                        vec![(1, 1)],
                        vec![(2, 1)],
                        vec![(3, 1)],
                    ]),
                    real_cone_decomposition_exact_coefficients: Some(vec![
                        "1".to_string(),
                        "1".to_string(),
                        "1".to_string(),
                        "1".to_string(),
                    ]),
                    real_cone_decomposition_exact_kind: Some("integer_semigroup"),
                    ambient_nonzero: vec![(4, 1)],
                    basis_nonzero: vec![(4, 1)],
                },
                MissingGvTargetSample {
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
                    ambient_nonzero: vec![(1, 1)],
                    basis_nonzero: vec![(1, 1)],
                },
            ],
        };
        let degree_bounded_context = vec![
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, 1)],
                basis_nonzero: vec![(0, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(1, 1)],
                basis_nonzero: vec![(1, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, -1), (2, 1)],
                basis_nonzero: vec![(2, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, -1), (3, 1)],
                basis_nonzero: vec![(3, 1)],
            },
        ];
        let covered_context = vec![CoveredToricGvContextSample {
            degree: 1,
            gv: "7".to_string(),
            ambient_nonzero: vec![(0, 1)],
            basis_nonzero: vec![(0, 1)],
        }];

        let active_classes = active_noncovered_dependency_source_ray_classes(
            &target_stats,
            &degree_bounded_context,
            &covered_context,
            4,
        )
        .unwrap();
        assert_eq!(active_classes, vec![vec![-1, 0, 0, 1], vec![-1, 0, 1, 0]]);
        assert_eq!(
            merge_unique_ambient_classes(vec![vec![-1, 0, 1, 0]], active_classes),
            vec![vec![-1, 0, 0, 1], vec![-1, 0, 1, 0]]
        );
    }

    #[test]
    fn shared_facet_unresolved_source_rays_skip_existing_source_context() {
        let witness = OriginCircuitWitnessSample {
            first_facet_exclusive_point: 1,
            second_facet_exclusive_point: 3,
            shared_two_simplex: vec![2],
            shared_two_simplex_points: Vec::new(),
            shared_two_simplex_star_simplices: Vec::new(),
            shared_two_simplex_star_extra_point_samples: Vec::new(),
            first_facet: vec![1, 2],
            second_facet: vec![2, 3],
            first_facet_size: 2,
            second_facet_size: 2,
            sparse_relation: vec![(0, -1), (1, 1), (3, 1)],
            relation_points: Vec::new(),
        };
        let target_stats = MissingGvTargetStats {
            target_count: 1,
            targets_that_are_mori_generators: 0,
            targets_that_are_origin_circuits: 1,
            targets_real_cone_decomposable_by_other_generators: 0,
            targets_that_are_lp_extremal_mori_generators: 0,
            real_cone_decomposition_active_generator_min: None,
            real_cone_decomposition_active_generator_max: None,
            origin_circuit_resolved_conifold_count: 0,
            min_generators_le_target_degree: 0,
            max_generators_le_target_degree: 0,
            origin_coefficient_counts: BTreeMap::new(),
            origin_circuit_pattern_counts: BTreeMap::new(),
            origin_circuit_affine_rank_counts: BTreeMap::new(),
            branch_status_counts: BTreeMap::new(),
            branch_bucket_counts: BTreeMap::new(),
            real_cone_decomposition_exact_kind_counts: BTreeMap::new(),
            sample: vec![MissingGvTargetSample {
                degree: 4,
                generators_le_degree: 0,
                is_mori_generator: false,
                origin_circuit_pattern: None,
                origin_circuit_witness_count: Some(1),
                origin_circuit_first_witness: Some(witness),
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
                ambient_nonzero: vec![(0, -1), (1, 1), (3, 1)],
                basis_nonzero: vec![(2, 1)],
            }],
        };
        let degree_bounded_context = vec![
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, -1), (1, 1)],
                basis_nonzero: vec![(0, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, -1), (2, 1)],
                basis_nonzero: vec![(1, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 1,
                ambient_nonzero: vec![(0, -1), (3, 1)],
                basis_nonzero: vec![(3, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 2,
                ambient_nonzero: vec![(0, -1), (2, 1), (3, 1)],
                basis_nonzero: vec![(4, 1)],
            },
            DegreeBoundedMoriRayContextSample {
                degree: 2,
                ambient_nonzero: vec![(0, -1), (4, 1)],
                basis_nonzero: vec![(5, 1)],
            },
        ];
        let covered_context = vec![CoveredToricGvContextSample {
            degree: 1,
            gv: "7".to_string(),
            ambient_nonzero: vec![(0, -1), (1, 1)],
            basis_nonzero: vec![(0, 1)],
        }];
        let source_derived_toric_context = vec![ToricGvDiagnosticContextSample {
            degree: 1,
            gv: "-2".to_string(),
            source_bucket: "two_face".to_string(),
            source_summary: "diagnostic".to_string(),
            ambient_nonzero: vec![(0, -1), (2, 1)],
            basis_nonzero: vec![(1, 1)],
        }];
        let already_diagnosed_source_stats = MissingGvTargetStats {
            target_count: 1,
            targets_that_are_mori_generators: 1,
            targets_that_are_origin_circuits: 0,
            targets_real_cone_decomposable_by_other_generators: 0,
            targets_that_are_lp_extremal_mori_generators: 1,
            real_cone_decomposition_active_generator_min: None,
            real_cone_decomposition_active_generator_max: None,
            origin_circuit_resolved_conifold_count: 0,
            min_generators_le_target_degree: 0,
            max_generators_le_target_degree: 0,
            origin_coefficient_counts: BTreeMap::new(),
            origin_circuit_pattern_counts: BTreeMap::new(),
            origin_circuit_affine_rank_counts: BTreeMap::new(),
            branch_status_counts: BTreeMap::new(),
            branch_bucket_counts: BTreeMap::new(),
            real_cone_decomposition_exact_kind_counts: BTreeMap::new(),
            sample: vec![MissingGvTargetSample {
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
                ambient_nonzero: vec![(0, -1), (3, 1)],
                basis_nonzero: vec![(3, 1)],
            }],
        };

        let classes = shared_facet_unresolved_source_ray_classes(
            &target_stats,
            &degree_bounded_context,
            &covered_context,
            &source_derived_toric_context,
            &already_diagnosed_source_stats,
            5,
        )
        .unwrap();
        assert_eq!(classes, vec![vec![-1, 0, 1, 1, 0]]);
    }

    #[test]
    fn missing_gv_target_stats_records_branch_cut_status() {
        let kahler = vec![F64::<Finite>::new(-0.1).expect("finite")];
        let gamma = vec![I64::<Finite>::new(0)];
        let stats = missing_gv_target_stats(
            &[vec![1, 0]],
            &[],
            None,
            &[vec![1]],
            &[0],
            &[1],
            Some(&kahler),
            Some(&gamma),
            1,
            &HashMap::new(),
            &HashMap::new(),
            true,
            1,
        )
        .unwrap();

        assert_eq!(
            stats.branch_status_counts,
            BTreeMap::from([("RealBranchCut".to_string(), 1)])
        );
        assert_eq!(
            stats.branch_bucket_counts,
            BTreeMap::from([(
                "parity_mod2=0;qdot_bin=nonpositive;status=RealBranchCut".to_string(),
                1
            )])
        );
        assert_eq!(
            stats.sample[0].branch_diagnostic,
            Some(MissingGvBranchDiagnostic {
                q_dot_t: "-0.10000000000000001".to_string(),
                parity: 0,
                parity_mod2: 0,
                q_dot_bucket: "nonpositive",
                dilog_status: "RealBranchCut".to_string(),
            })
        );
    }

    #[test]
    fn missing_gv_target_stats_records_origin_circuit_affine_support() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 1]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
        ];
        let class = vec![-1, -1, 1, 1];
        let first_witness = cyrus_core::OriginCircuitCurveWitness {
            class: class.clone(),
            first_facet_exclusive_point: 2,
            second_facet_exclusive_point: 3,
            shared_two_simplex: vec![1],
            first_facet: vec![1, 2],
            second_facet: vec![1, 3],
            sparse_relation: vec![(0, -1), (1, -1), (2, 1), (3, 1)],
            relation_points: Vec::new(),
        };
        let second_witness = cyrus_core::OriginCircuitCurveWitness {
            class: class.clone(),
            first_facet_exclusive_point: 3,
            second_facet_exclusive_point: 2,
            shared_two_simplex: vec![1],
            first_facet: vec![1, 3],
            second_facet: vec![1, 2],
            sparse_relation: vec![(0, -1), (1, -1), (2, 1), (3, 1)],
            relation_points: Vec::new(),
        };
        let mut origin_circuits_by_class = HashMap::new();
        origin_circuits_by_class.insert(
            class.clone(),
            cyrus_core::OriginCircuitCurveDiagnostic {
                class: class.clone(),
                origin_coefficient: -1,
                negative_coefficient_counts: BTreeMap::from([(-1, 1)]),
                positive_coefficient_counts: BTreeMap::from([(1, 2)]),
                is_resolved_conifold_pattern: true,
                witnesses: vec![first_witness, second_witness],
            },
        );
        let triangulation = Triangulation::new(vec![vec![0, 1, 2], vec![0, 1, 3]]);

        let stats = missing_gv_target_stats(
            std::slice::from_ref(&class),
            &points,
            Some(&triangulation),
            &[vec![-1, 1, 1]],
            &[1, 2, 3],
            &[3, 2, 2],
            None,
            None,
            0,
            &origin_circuits_by_class,
            &HashMap::new(),
            true,
            1,
        )
        .unwrap();

        assert_eq!(stats.targets_that_are_origin_circuits, 1);
        assert_eq!(
            stats.origin_circuit_affine_rank_counts,
            BTreeMap::from([(2, 1)])
        );
        let support = stats.sample[0]
            .origin_circuit_affine_support
            .as_ref()
            .expect("origin circuit should include affine support diagnostics");
        assert_eq!(support.affine_rank, 2);
        assert_eq!(
            support.coefficient_counts,
            BTreeMap::from([(-1, 2), (1, 2)])
        );
        assert_eq!(support.local_charge_basis, vec![vec![1, 1, -1, -1]]);
        assert_eq!(support.local_coordinates.len(), 4);
        for point in &support.local_coordinates {
            assert_eq!(point.coordinates.len(), 2);
        }
        assert!(support.local_coordinates_2d.is_some());
        assert_eq!(stats.sample[0].origin_circuit_witness_count, Some(2));
        assert_eq!(
            stats.sample[0]
                .origin_circuit_witnesses
                .as_ref()
                .map(Vec::len),
            Some(2)
        );
        let first_witness = stats.sample[0]
            .origin_circuit_first_witness
            .as_ref()
            .expect("origin circuit should serialize the first witness");
        assert_eq!(first_witness.shared_two_simplex_points.len(), 1);
        assert_eq!(first_witness.shared_two_simplex_points[0].point_index, 1);
        assert_eq!(first_witness.shared_two_simplex_points[0].coefficient, 0);
        assert_eq!(
            first_witness.shared_two_simplex_points[0].coordinates,
            vec![1, 1]
        );
        assert_eq!(
            first_witness.shared_two_simplex_star_simplices,
            vec![vec![0, 1, 2], vec![0, 1, 3]]
        );
        let star_extra_points = first_witness
            .shared_two_simplex_star_extra_point_samples
            .iter()
            .map(|simplex| {
                simplex
                    .iter()
                    .map(|point| {
                        (
                            point.point_index,
                            point.coefficient,
                            point.coordinates.clone(),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            star_extra_points,
            vec![vec![(2, 0, vec![1, 0])], vec![(3, 0, vec![0, 1])]]
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
    fn exact_active_generator_coefficients_classify_semigroup_vs_rational_cone() {
        let integer_coefficients =
            exact_active_generator_coefficients(&[1, 1], &[vec![1, 0], vec![0, 1]], &[0, 1])
                .unwrap();
        assert_eq!(
            integer_coefficients,
            vec![malachite::Rational::from(1), malachite::Rational::from(1)]
        );
        assert_eq!(
            classify_exact_active_generator_coefficients(&integer_coefficients),
            "integer_semigroup"
        );

        let rational_coefficients =
            exact_active_generator_coefficients(&[1, 1], &[vec![2, 0], vec![0, 2]], &[0, 1])
                .unwrap();
        assert_eq!(
            rational_coefficients,
            vec![
                malachite::Rational::from(1) / malachite::Rational::from(2),
                malachite::Rational::from(1) / malachite::Rational::from(2)
            ]
        );
        assert_eq!(
            classify_exact_active_generator_coefficients(&rational_coefficients),
            "rational_cone"
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
                linear_system: Some(CmsGeneralDivisorLinearSystemDiagnostic {
                    row_count: 3,
                    column_count: 2,
                    rank: 2,
                    augmented_rank: 2,
                }),
                solution_basis_support_len: Some(2),
                solution_basis_nonzero: Some(vec![(0, "1".to_string()), (1, "1".to_string()),]),
                solution_ambient_basis_nonzero: Some(vec![
                    (1, "1".to_string()),
                    (2, "1".to_string()),
                ]),
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

        assert_eq!(integer_dot_for_test(&certificate.normal, &[1, 0]), 0);
        assert!(integer_dot_for_test(&certificate.normal, &[0, 1]) > 0);
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
    fn supporting_mori_face_aggregate_lp_finds_higher_codimension_normal() {
        let lp_normal = solve_supporting_face_normal_aggregate_lp(
            &[vec![1, 0, 0]],
            &[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1], vec![1, 1, 1]],
        )
        .unwrap()
        .unwrap();

        let certificate = integer_supporting_face_certificate_from_lp(
            &lp_normal,
            &[vec![1, 0, 0]],
            &[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1], vec![1, 1, 1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(integer_dot_for_test(&certificate.normal, &[1, 0, 0]), 0);
        assert!(certificate.zero_generator_count >= 1);
        assert!(certificate.positive_generator_count >= 1);
    }

    fn integer_dot_for_test(lhs: &[i64], rhs: &[i64]) -> i128 {
        lhs.iter()
            .zip(rhs)
            .map(|(&left, &right)| i128::from(left) * i128::from(right))
            .sum()
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
    fn exact_generator_decomposition_honors_four_term_limit() {
        assert!(
            exact_generator_decompositions(&[4], &[vec![1], vec![4]], &[1], 4, 3, 4)
                .unwrap()
                .is_none()
        );

        let mut decompositions =
            exact_generator_decompositions(&[4], &[vec![1], vec![4]], &[1], 4, 4, 4)
                .unwrap()
                .unwrap();
        let decomposition = decompositions.pop().unwrap();

        assert_eq!(decomposition, vec![vec![1], vec![1], vec![1], vec![1]]);
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
    fn cms_general_divisor_formula_sum_preserves_candidate_multiplicity() {
        let candidates = vec![
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 2,
                shrinking_divisor_coefficient: -3,
                shrinking_divisor_coordinates: vec![1, 0],
                inferred_other_normal_degree: 1,
                toric_gv1_formula_value: Some(-2),
                all_non_origin_relation_points_are_two_face: true,
            },
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 5,
                shrinking_divisor_coefficient: -1,
                shrinking_divisor_coordinates: vec![0, 1],
                inferred_other_normal_degree: -1,
                toric_gv1_formula_value: Some(1),
                all_non_origin_relation_points_are_two_face: true,
            },
            CmsGeneralDivisorShapeCandidate {
                shrinking_divisor_index: 7,
                shrinking_divisor_coefficient: -1,
                shrinking_divisor_coordinates: vec![1, 1],
                inferred_other_normal_degree: -1,
                toric_gv1_formula_value: Some(1),
                all_non_origin_relation_points_are_two_face: true,
            },
        ];

        assert_eq!(
            cms_general_divisor_formula_sum(Some(&candidates)),
            Some(malachite::Integer::from(0))
        );
        assert_eq!(cms_general_divisor_formula_sum(None), None);
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
    fn gv_weight_lp_recovers_continuous_drop_solution() {
        let contributions = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let target = vec![
            F64::<Finite>::new(0.25).unwrap(),
            F64::<Finite>::new(1.0).unwrap(),
        ];

        let report = solve_gv_weight_lp("test", &contributions, &target, 0.0, 1.0).unwrap();

        assert!(report.max_abs_delta < 1e-8);
        assert!(report.relative_l2_delta < 1e-8);
        assert!((report.weights[0] - 0.25).abs() < 1e-8);
        assert!((report.weights[1] - 1.0).abs() < 1e-8);
    }

    #[test]
    fn gf2_solver_finds_binary_solution_and_rejects_inconsistency() {
        let rows = vec![vec![1, 0, 1], vec![0, 1, 1]];
        let rhs = vec![1, 0];

        let solution = solve_gf2_system(&rows, &rhs, 3).unwrap();
        assert_eq!((solution[0] ^ solution[2]) & 1, 1);
        assert_eq!((solution[1] ^ solution[2]) & 1, 0);

        let inconsistent_rows = vec![vec![1, 0], vec![1, 0]];
        let inconsistent_rhs = vec![0, 1];
        assert!(solve_gf2_system(&inconsistent_rows, &inconsistent_rhs, 2).is_none());
    }

    #[test]
    fn projected_ambient_gamma_preserves_curve_parity() {
        let curve_basis = vec![
            vec![
                malachite::Integer::from(1),
                malachite::Integer::from(0),
                malachite::Integer::from(4),
            ],
            vec![
                malachite::Integer::from(0),
                malachite::Integer::from(1),
                malachite::Integer::from(-1),
            ],
        ];
        let ambient_gamma = vec![
            I64::<Finite>::new(1),
            I64::<Finite>::new(0),
            I64::<Finite>::new(1),
        ];

        let basis_gamma =
            project_ambient_gamma_to_curve_basis(&curve_basis, &ambient_gamma).unwrap();

        assert_eq!(
            basis_gamma
                .iter()
                .map(|entry| entry.get())
                .collect::<Vec<_>>(),
            vec![5, -1]
        );
        let gv_invariants = vec![(vec![2, 3, 5], malachite::Integer::from(1))];
        assert_eq!(
            count_gamma_parity_mismatches(&gv_invariants, &[0, 1], &ambient_gamma, &basis_gamma),
            Some(0)
        );
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
