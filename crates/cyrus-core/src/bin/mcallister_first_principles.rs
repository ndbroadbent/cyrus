//! First-principles McAllister pipeline runner (opt-in).
//!
//! Runs the pipeline stages explicitly and fails fast on any mismatch or
//! invalid physics. This is intentionally a binary, not a test.
//!
//! Optional:
//! - `--dual-basis path/to/dual_basis.json` to supply the flux coordinate basis.
//! - `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.
//! - `--branch-candidates N --branch-selection <max-volume|min-volume|first-positive|min-condition>`
//!   to run deterministic KKLT branch search without loading Kähler checkpoints.
//! - `--branch-report-jsonl path` to export positive phase-1 branch candidates
//!   discovered by that search for GA/debug ranking.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

use cyrus_core::flat_direction::{compute_flat_direction, compute_flat_direction_full};
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
    compute_regular_triangulation, compute_toric_two_face_curve_gv_invariants,
    compute_w0_from_terms, generate_scaled_kklt_branch_initializations, intersection_in_basis,
    is_unimodular, remove_pair_decomposable_curve_candidates, solve_mixed_basis_path_following,
    solve_mixed_basis_path_following_branch_candidates, solve_racetrack,
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
}

impl BranchSelection {
    fn as_str(self) -> &'static str {
        match self {
            Self::MaxVolume => "max-volume",
            Self::MinVolume => "min-volume",
            Self::FirstPositive => "first-positive",
            Self::MinCondition => "min-condition",
        }
    }
}

fn parse_branch_selection(name: &str) -> Option<BranchSelection> {
    match name {
        "max-volume" => Some(BranchSelection::MaxVolume),
        "min-volume" => Some(BranchSelection::MinVolume),
        "first-positive" => Some(BranchSelection::FirstPositive),
        "min-condition" => Some(BranchSelection::MinCondition),
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
struct DualPointsFixture {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct DualSimplicesFixture {
    simplices: Vec<Vec<usize>>,
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
    selected_phase1_volume: f64,
    selected_phase1_rel_err: f64,
    selected_jacobian_rank: usize,
    selected_jacobian_max_rank: usize,
    selected_jacobian_condition_number: Option<f64>,
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
    branch_report_path: Option<String>,
    dual_basis_override: Option<BasisOverride>,
}

struct PrimalGeom {
    polytope: Polytope,
    triangulation_points: Vec<Point>,
    triangulation: Triangulation,
}

struct PrimalIntersection {
    glsm: Vec<Vec<malachite::Integer>>,
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
    let branch_report_path = parse_arg_value::<String>("--branch-report-jsonl");
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
        branch_report_path,
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

fn load_dual_inputs(
    data_dir: Option<&str>,
    manifest_dir: &PathBuf,
) -> (Vec<Vec<i64>>, Vec<Vec<usize>>) {
    data_dir.map_or_else(
        || {
            eprintln!("[WARN] using dual points/simplices JSON fixtures");
            let dual_points_path =
                manifest_dir.join("tests/mcallister_e2e/overrides/dual_points.json");
            let dual_points: DualPointsFixture = load_json(&dual_points_path);
            let dual_simplices_path =
                manifest_dir.join("tests/mcallister_e2e/overrides/dual_simplices.json");
            let dual_simplices: DualSimplicesFixture = load_json(&dual_simplices_path);
            (dual_points.points, dual_simplices.simplices)
        },
        |dir| {
            let dir = PathBuf::from(dir);
            let points = read_points(&dir.join("dual_points.dat"));
            let simplices_i64 = read_points(&dir.join("dual_simplices.dat"));
            let simplices = simplices_i64
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|v| usize::try_from(v).expect("simplex index must be non-negative"))
                        .collect::<Vec<usize>>()
                })
                .collect();
            (points, simplices)
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

fn simplex_point_count(simplices: &[Vec<usize>]) -> usize {
    simplices
        .iter()
        .flatten()
        .copied()
        .max()
        .map_or(0, |idx| idx + 1)
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
        triangulation_points,
        triangulation,
    }
}

fn stage_intersection(
    data_dir: Option<&str>,
    geom: &PrimalGeom,
    t0: &Instant,
) -> PrimalIntersection {
    let (glsm, _linrels, basis) =
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
    allow_invalid_ek0: bool,
    dual_basis_override: Option<&BasisOverride>,
    t0: &Instant,
) -> FlatDirectionData {
    let (dual_points_raw, dual_simplices_raw) = load_dual_inputs(data_dir, manifest_dir);
    let dual_all_points: Vec<Point> = dual_points_raw
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();
    let dual_polytope = Polytope::from_vertices(dual_all_points).expect("Failed dual polytope");
    let triangulation_point_count = simplex_point_count(&dual_simplices_raw);
    let dual_points_vec: Vec<Point> = dual_points_raw
        .iter()
        .take(triangulation_point_count)
        .map(|coords| Point::new(coords.clone()))
        .collect();
    let dual_triangulation = Triangulation::new(dual_simplices_raw);
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

fn write_branch_report_jsonl(
    path: &PathBuf,
    ctx: &BranchReportContext,
    branches_by_volume: &[cyrus_core::KkltBranchSolution],
    t_initializations: &[Vec<F64<Finite>>],
) -> Result<(), String> {
    let Some(selected) = branches_by_volume.get(ctx.selected_rank_by_volume) else {
        return Err(format!(
            "selected branch rank {} is outside {} branch rows",
            ctx.selected_rank_by_volume,
            branches_by_volume.len()
        ));
    };

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
        selected_phase1_volume: selected.classical_volume.get(),
        selected_phase1_rel_err: selected.result.relative_error.get(),
        selected_jacobian_rank: selected.jacobian_diagnostics.rank,
        selected_jacobian_max_rank: selected.jacobian_diagnostics.max_rank,
        selected_jacobian_condition_number: selected
            .jacobian_diagnostics
            .condition_number
            .map(|value| value.get()),
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
        let row = BranchReportBranch {
            record_type: "positive_branch",
            branch_seed: ctx.branch_seed,
            branch_selection: ctx.branch_selection.as_str(),
            rank_by_volume,
            selected: rank_by_volume == ctx.selected_rank_by_volume,
            init_index: branch.init_index,
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
        };
        lines.push_str(&serde_json::to_string(&row).map_err(|e| e.to_string())?);
        lines.push('\n');
    }

    std::fs::write(path, lines).map_err(|e| e.to_string())
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
    branch_report_path: Option<&str>,
    small_curve_cutoff: F64<Pos>,
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

    let t = if allow_downstream_kahler {
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
        transform_kahler_to_computed_basis(
            &intersection.glsm,
            &intersection.basis,
            &source_basis,
            &t_raw,
        )
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
        eprintln!(
            "[WARN] corrected_kahler_param.dat remains validation-only replay and is not loaded without --allow-downstream-kahler."
        );
        let Some(tau_target) =
            cyrus_core::kklt::compute_corrected_target_tau(&c_i, &chi_divisor, c_tau)
        else {
            eprintln!("[ERROR] corrected KKLT target construction failed");
            std::process::exit(2);
        };
        let zeroth_order = if branch_candidates == 0 {
            let Some(result) = cyrus_core::kklt::solve_two_phase_mixed_basis_path_following(
                &intersection.kappa_basis,
                &intersection.kappa_full,
                &intersection.basis,
                &kklt_basis,
                &tau_target,
                &c_i,
                CheckedRange::new(0, kklt_steps),
            ) else {
                eprintln!("[ERROR] zeroth-order mixed-basis KKLT path-following failed");
                std::process::exit(2);
            };
            result
        } else {
            let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
            let Some(t_initializations) = generate_scaled_kklt_branch_initializations(
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
            };
            let best_branch = positive_branches[selected_rank_by_volume].clone();
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
                )
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[ERROR] failed to write KKLT branch report {}: {e}",
                        report_path.display()
                    );
                    std::process::exit(2);
                });
                eprintln!("[INFO] wrote KKLT branch report {}", report_path.display());
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
            result
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
            &zeroth_order.t,
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
            "[INFO] primal small toric curves: ambient_rays={} subcutoff={} filtered_hilbert_candidates={} cutoff={}",
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
        let max_gv_iterations = 20usize;
        let gv_tolerance = 1e-10f64;
        let mut gv_converged = false;
        for iter in 0..max_gv_iterations {
            let previous_t = corrected.t.clone();
            let Some(gv_correction) =
                cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
                    &small_curve_gvs,
                    &intersection.basis,
                    &kklt_basis,
                    &previous_t,
                    None,
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
        corrected.t
    };

    let classical_volume = classical_volume_from_t(&intersection.kappa_basis, &t);
    let h11 = cyrus_core::types::i32::I32::<cyrus_core::types::tags::GTEOne>::new(214).unwrap();
    let h21 = cyrus_core::types::i32::I32::<cyrus_core::types::tags::NonNeg>::new(4).unwrap();
    let bbhl = bbhl_correction(h11, h21);
    let v_string = classical_volume - bbhl.get();
    let v_string_pos = F64::<Pos>::new(v_string).expect("V_string must be positive");
    eprintln!("[TIME] volume: {:.2?}", t0.elapsed());
    (v_string, v_string_pos)
}

fn compare_against_dat(
    compare_dir: Option<String>,
    data_dir: Option<String>,
    g_s: F64<Pos>,
    w0: F64<Pos>,
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
    compare_against_dat(compare_dir, data_dir, g_s, racetrack.w0);
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
        args.branch_report_path.as_deref(),
        small_curve_cutoff,
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
        ] {
            let policy = parse_branch_selection(name)
                .unwrap_or_else(|| panic!("policy {name} should parse"));
            assert_eq!(policy.as_str(), name);
        }
        assert!(parse_branch_selection("condition").is_none());
    }
}

fn main() {
    run_pipeline(parse_args());
}
