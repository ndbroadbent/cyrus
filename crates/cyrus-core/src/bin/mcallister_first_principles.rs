//! First-principles McAllister pipeline runner (opt-in).
//!
//! Runs the pipeline stages explicitly and fails fast on any mismatch or
//! invalid physics. This is intentionally a binary, not a test.
//!
//! Optional:
//! - `--dual-basis path/to/dual_basis.json` to supply explicit dual basis indices.
//! - `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.

use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::time::Instant;

use cyrus_core::flat_direction::{compute_flat_direction, compute_flat_direction_full};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::vacuum::compute_vacuum;
use cyrus_core::volume::bbhl_correction;
use cyrus_core::{
    build_racetrack_terms, compute_curve_basis_matrix, compute_glsm_and_linrels,
    compute_grading_vector, compute_intersection_cytools, compute_linear_relations_no_origin,
    compute_mori_cone_cap_rays, compute_regular_triangulation, compute_w0, solve_racetrack,
    intersection_in_basis, Point, Polytope, Triangulation,
};

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
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>().unwrap_or_else(|e| panic!("Invalid float {s} in {}: {e}", path.display())))
        .collect()
}

fn read_csv_i64(path: &PathBuf) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .flat_map(|line| line.split(','))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<i64>().unwrap_or_else(|e| panic!("Invalid int {s} in {}: {e}", path.display())))
        .collect()
}

fn read_csv_usize(path: &PathBuf) -> Vec<usize> {
    read_csv_i64(path)
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index must be non-negative"))
        .collect()
}

fn assert_basis_matches_dat(computed: &[usize], data_dir: &str, label: &str) {
    let basis_path = PathBuf::from(data_dir).join("basis.dat");
    let expected = read_csv_usize(&basis_path);
    if computed != expected {
        eprintln!(
            "[ERROR] {label} basis mismatch: computed len={} expected len={}",
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
            eprintln!("[ERROR]   idx {i}: computed={a}, expected={b}");
        }
        eprintln!("[ERROR]   total mismatches: {}", diff.len());
        std::process::exit(2);
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
                .map(|s| s.trim().parse::<i64>().unwrap_or_else(|e| panic!("Invalid point entry {s} in {}: {e}", path.display())))
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

fn main() {
    let t0 = Instant::now();
    let stop_after = parse_arg_value::<String>("--stop-after")
        .and_then(|s| parse_stage(&s))
        .unwrap_or(Stage::Vacuum);
    let max_deg = parse_arg_value::<u32>("--max-deg");
    let min_points = if max_deg.is_some() {
        None
    } else {
        parse_arg_value::<u32>("--min-points").or(Some(100))
    };
    let cutoff = parse_arg_value::<f64>("--cutoff").unwrap_or(1.0);
    let out_path = parse_arg_value::<String>("--out");
    let compare_dir = parse_arg_value::<String>("--compare-dir");
    let data_dir = parse_arg_value::<String>("--data-dir")
        .or_else(|| std::env::var("CYRUS_MCALLISTER_DATA_DIR").ok());
    let allow_invalid_ek0 = parse_flag("--allow-invalid-ek0");
    let allow_fixtures = parse_flag("--allow-fixtures");
    let first_principles_env = std::env::var_os("CYRUS_FIRST_PRINCIPLES").is_some();

    if allow_fixtures && first_principles_env {
        eprintln!("[ERROR] --allow-fixtures cannot be used with CYRUS_FIRST_PRINCIPLES");
        std::process::exit(2);
    }
    let dual_basis_override =
        parse_arg_value::<String>("--dual-basis").map(|path| load_json::<BasisOverride>(&PathBuf::from(path)));

    if let Some(dir) = &data_dir {
        eprintln!("[INFO] using McAllister data dir: {}", dir);
        eprintln!("[MODE] first-principles (.dat)");
    } else if !allow_fixtures {
        eprintln!(
            "[ERROR] No McAllister data dir set. Refusing to fall back to JSON fixtures."
        );
        eprintln!("[ERROR] Set CYRUS_MCALLISTER_DATA_DIR or pass --allow-fixtures.");
        std::process::exit(2);
    } else {
        eprintln!("[MODE] fixtures (JSON)");
        eprintln!("[WARN] using JSON fixtures (not a first-principles run)");
    }

    let cutoff = F64::<Pos>::new(cutoff).expect("cutoff must be positive");

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // === Stage 1: Load polytope ===
    let (points_raw, heights_raw) = if let Some(dir) = &data_dir {
        let dir = PathBuf::from(dir);
        let points = read_points(&dir.join("points.dat"));
        let heights = read_csv_f64(&dir.join("heights.dat"));
        (points, heights)
    } else {
        let poly_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
        let input: PolytopeInput = load_json(&poly_path);
        let heights_path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
        let heights: HeightsInput = load_json(&heights_path);
        (input.points, heights.values)
    };

    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("Failed to build polytope");

    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter triangulation points");
    let _origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found");

    // === Stage 2: Triangulation (McAllister heights) ===
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights_raw)
        .expect("Failed to compute triangulation");
    eprintln!("[TIME] triangulation: {:.2?}", t0.elapsed());
    if !stage_enabled(Stage::Intersection, stop_after) {
        return;
    }

    // === Stage 3: GLSM + intersection numbers (primal) ===
    let (_glsm, linrels, basis) =
        compute_glsm_and_linrels(&triangulation_points).expect("Failed GLSM/linrels");
    if let Some(dir) = &data_dir {
        assert_basis_matches_dat(&basis, dir, "primal");
    }
    let points_i64: Vec<Vec<i64>> = triangulation_points
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
        &triangulation,
        &triangulation_points,
        &linrels_i64,
    )
    .expect("Failed intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);
    eprintln!("[TIME] intersection: {:.2?}", t0.elapsed());

    if !stage_enabled(Stage::FlatDirection, stop_after) {
        return;
    }

    // === Stage 4: Flat direction (dual data from .dat or fixtures fallback) ===
    let (dual_points_raw, dual_simplices_raw): (Vec<Vec<i64>>, Vec<Vec<usize>>) = if let Some(dir) = &data_dir {
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
    } else {
        eprintln!("[WARN] using dual points/simplices JSON fixtures");
        let dual_points_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_points.json");
        let dual_points: DualPointsFixture = load_json(&dual_points_path);
        let dual_simplices_path =
            manifest_dir.join("tests/mcallister_e2e/overrides/dual_simplices.json");
        let dual_simplices: DualSimplicesFixture = load_json(&dual_simplices_path);
        (dual_points.points, dual_simplices.simplices)
    };

    let dual_points_vec: Vec<Point> = dual_points_raw
        .iter()
        .take(9)
        .map(|coords| Point::new(coords.clone()))
        .collect();
    let dual_triangulation = Triangulation::new(dual_simplices_raw);
    let (_dual_glsm, _dual_linrel, dual_basis_auto) =
        compute_glsm_and_linrels(&dual_points_vec).expect("Failed dual GLSM");
    let dual_basis = if let Some(explicit_basis) = &dual_basis_override {
        eprintln!("[INFO] using explicit dual basis from --dual-basis");
        explicit_basis.indices.clone()
    } else {
        eprintln!(
            "[INFO] using computed dual basis (len={}, basis={:?})",
            dual_basis_auto.len(),
            dual_basis_auto
        );
        dual_basis_auto
    };
    if dual_basis_override.is_none() && dual_basis != vec![3, 4, 5, 8] {
        eprintln!(
            "[ERROR] computed dual basis {:?} does not match McAllister [3,4,5,8]",
            dual_basis
        );
        std::process::exit(2);
    }

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

    let (k_raw, m_raw) = if let Some(dir) = &data_dir {
        let dir = PathBuf::from(dir);
        let k = read_csv_i64(&dir.join("K_vec.dat"));
        let m = read_csv_i64(&dir.join("M_vec.dat"));
        (k, m)
    } else {
        let flux_path = manifest_dir.join("tests/mcallister_e2e/inputs/flux.json");
        let flux: FluxInput = load_json(&flux_path);
        (flux.k, flux.m)
    };
    let k_flux: Vec<I64<Finite>> = k_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let m_flux: Vec<I64<Finite>> = m_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();

    let (flat_p, ek0_opt) = match compute_flat_direction_full(&dual_kappa, &k_flux, &m_flux) {
        Some(v) => (v.p, Some(v.ek0)),
        None if allow_invalid_ek0 => {
            eprintln!("[WARN] invalid flat direction (ek0<=0); continuing due to --allow-invalid-ek0");
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

    if !stage_enabled(Stage::Gv, stop_after) {
        return;
    }

    // === Stage 5: GV invariants ===
    let rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        true,
        false,
        Some(&basis),
    )
    .expect("Failed mori cone cap rays");
    let grading = compute_grading_vector(&rays).expect("No grading vector found");

    let curve_basis =
        compute_curve_basis_matrix(&linrels, &basis).expect("Failed curve basis matrix");
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
        &kappa_basis,
        min_points,
        max_deg,
    )
    .expect("GV computation failed");
    eprintln!("[TIME] gv_invariants: {:.2?}", t0.elapsed());

    if !stage_enabled(Stage::Racetrack, stop_after) {
        return;
    }

    // === Stage 6-8: Racetrack + W0 ===
    let gv: Vec<cyrus_core::racetrack::GvInvariant> = invariants
        .into_iter()
        .map(|(curve, value)| {
            let curve = curve
                .into_iter()
                .map(|v| I64::<Finite>::new(i64::from(v)))
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

    let terms = build_racetrack_terms(&gv, &m_flux, &flat_p, cutoff);
    if terms.len() < 2 {
        eprintln!("[ERROR] not enough racetrack terms");
        std::process::exit(2);
    }
    let Some(rt_res) = solve_racetrack(&terms) else {
        eprintln!("[ERROR] no stable racetrack solution");
        std::process::exit(2);
    };
    let w0 = compute_w0(&rt_res, &terms[0], &terms[1]);
    eprintln!("[TIME] racetrack: {:.2?}", t0.elapsed());

    if !stage_enabled(Stage::Volume, stop_after) {
        return;
    }

    // === Stage 9-10: V_string (primal) ===
    let kappa_primal = intersection_in_basis(&kappa_full, &basis);

    let Some(data_dir_path) = data_dir.clone().map(PathBuf::from) else {
        eprintln!(
            "[ERROR] Volume stage requires McAllister data dir for corrected_kahler_param.dat"
        );
        std::process::exit(2);
    };
    let corrected_t = read_csv_f64(&data_dir_path.join("corrected_kahler_param.dat"));

    let mut volume_sum = 0.0f64;
    for (&(i, j, k), val) in kappa_primal.iter() {
        let (kappa_val, _): (f64, _) =
            malachite::num::conversion::traits::RoundingFrom::rounding_from(
                val.get(),
                malachite::rounding_modes::RoundingMode::Nearest,
            );
        let t_product = corrected_t[i] * corrected_t[j] * corrected_t[k];
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        volume_sum += mult * kappa_val * t_product;
    }
    let classical_volume = volume_sum / 6.0;
    let h11 = cyrus_core::types::i32::I32::<cyrus_core::types::tags::GTEOne>::new(214).unwrap();
    let h21 = cyrus_core::types::i32::I32::<cyrus_core::types::tags::NonNeg>::new(4).unwrap();
    let bbhl = bbhl_correction(h11, h21);
    let v_string = classical_volume - bbhl.get();
    let v_string_pos = F64::<Pos>::new(v_string).expect("V_string must be positive");
    eprintln!("[TIME] volume: {:.2?}", t0.elapsed());

    if !stage_enabled(Stage::Vacuum, stop_after) {
        return;
    }

    // === Stage 11: V0 ===
    let ek0 = match ek0_opt {
        Some(v) => v,
        None => {
            eprintln!("[ERROR] ek0 is invalid; cannot compute V0");
            std::process::exit(2);
        }
    };
    let g_s = rt_res.g_s;
    let w0_pos = F64::<Pos>::new(w0.get()).expect("W0 must be positive");
    let vac = compute_vacuum(ek0, g_s, v_string_pos, w0_pos);
    let v0_log10_abs = vac.v0.get().abs().log10();

    // === Comparison (optional) ===
    let compare_dir = compare_dir.or_else(|| data_dir.clone());
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
        eprintln!("[COMPARE] g_s rel_err = {}", rel_gs);
        eprintln!("[COMPARE] W0 rel_err = {}", rel_w0);
    }

    let assertion_path = manifest_dir.join("tests/mcallister_e2e/assertions/racetrack.json");
    let assertion: RacetrackAssertion = load_json(&assertion_path);

    let summary = PipelineSummary {
        g_s: g_s.get(),
        w0: w0.get(),
        v_string,
        v0_log10_abs,
        ek0: ek0.get(),
    };

    eprintln!("[RESULT] g_s = {}", summary.g_s);
    eprintln!("[RESULT] W0 = {}", summary.w0);
    eprintln!("[RESULT] V_string = {}", summary.v_string);
    eprintln!("[RESULT] log10(|V0|) = {}", summary.v0_log10_abs);

    // Strict assertions (fail fast, no silent fallback)
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

    if let Some(path) = out_path {
        let out_path = PathBuf::from(path);
        let data = serde_json::to_string_pretty(&summary).expect("serialize summary");
        std::fs::write(&out_path, data)
            .unwrap_or_else(|e| panic!("Failed to write {}: {e}", out_path.display()));
    }
}
