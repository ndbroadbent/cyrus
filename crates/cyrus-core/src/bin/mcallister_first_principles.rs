//! First-principles McAllister pipeline runner (opt-in).
//!
//! Runs the pipeline stages explicitly and fails fast on any mismatch or
//! invalid physics. This is intentionally a binary, not a test.
//!
//! Optional:
//! - `--dual-basis path/to/dual_basis.json` to supply the flux coordinate basis.
//! - `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.

use serde::{Deserialize, Serialize};
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
    compute_regular_triangulation, compute_w0_from_terms, intersection_in_basis, is_unimodular,
    solve_racetrack,
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

struct PipelineArgs {
    stop_after: Stage,
    max_deg: Option<u32>,
    min_points: Option<u32>,
    cutoff: f64,
    out_path: Option<String>,
    compare_dir: Option<String>,
    data_dir: Option<String>,
    allow_invalid_ek0: bool,
    allow_fixtures: bool,
    allow_downstream_kahler: bool,
    kklt_steps: usize,
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
    let out_path = parse_arg_value::<String>("--out");
    let compare_dir = parse_arg_value::<String>("--compare-dir");
    let data_dir = parse_arg_value::<String>("--data-dir")
        .or_else(|| std::env::var("CYRUS_MCALLISTER_DATA_DIR").ok());
    let allow_invalid_ek0 = parse_flag("--allow-invalid-ek0");
    let allow_fixtures = parse_flag("--allow-fixtures");
    let allow_downstream_kahler = parse_flag("--allow-downstream-kahler");
    let kklt_steps = parse_arg_value::<usize>("--kklt-steps").unwrap_or(200);
    let dual_basis_override = parse_arg_value::<String>("--dual-basis")
        .map(|path| load_json::<BasisOverride>(&PathBuf::from(path)));
    PipelineArgs {
        stop_after,
        max_deg,
        min_points,
        cutoff,
        out_path,
        compare_dir,
        data_dir,
        allow_invalid_ek0,
        allow_fixtures,
        allow_downstream_kahler,
        kklt_steps,
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

fn stage_volume(
    data_dir: Option<&str>,
    manifest_dir: &PathBuf,
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    racetrack: &RacetrackData,
    allow_downstream_kahler: bool,
    kklt_steps: usize,
    t0: &Instant,
) -> (f64, F64<Pos>) {
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
        let Some(zeroth_order) = cyrus_core::kklt::solve_two_phase_mixed_basis_path_following(
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
        eprintln!(
            "[INFO] zeroth-order mixed-basis KKLT converged={} rel_err={}",
            zeroth_order.converged,
            zeroth_order.relative_error.get()
        );
        eprintln!(
            "[ERROR] full primal KKLT GV corrections are not available in this runner yet; refusing to reuse dual racetrack GV data"
        );
        std::process::exit(2);
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

fn main() {
    run_pipeline(parse_args());
}
