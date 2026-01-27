//! One-off McAllister racetrack computation (opt-in).
//!
//! Computes flat direction, GV invariants, and the racetrack solution from
//! first principles. This is intentionally a binary, not a test.
//!
//! Use `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.

use serde::Deserialize;
use serde::Serialize;
use std::path::PathBuf;
use std::time::Instant;

use cyrus_core::{
    Point, Polytope, build_racetrack_terms, compute_curve_basis_matrix,
    compute_glsm_and_linrels, compute_grading_vector, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_mori_cone_cap_rays, compute_gv_invariants,
    compute_regular_triangulation, compute_w0, intersection_in_basis, solve_racetrack,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::flat_direction::{compute_n_matrix, solve_linear_system_faer};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct FluxInput {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

fn load_flux() -> FluxInput {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/flux.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()))
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &PathBuf) -> T {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
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

fn read_scalar(path: &std::path::Path) -> f64 {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .trim()
        .parse::<f64>()
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
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<usize>().expect("invalid usize"))
        .collect()
}

fn assert_basis_matches_dat(computed: &[usize], data_dir: &str) {
    let basis_path = PathBuf::from(data_dir).join("basis.dat");
    let expected = read_csv_usize(&basis_path);
    if computed != expected {
        eprintln!(
            "[ERROR] basis mismatch: computed len={} expected len={}",
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

#[derive(Debug, Serialize)]
struct RacetrackTermOut {
    exponent: f64,
    coefficient: f64,
}

#[derive(Debug, Serialize)]
struct RacetrackOut {
    g_s: f64,
    im_tau: f64,
    w0: f64,
    terms: Vec<RacetrackTermOut>,
    compare: Option<CompareOut>,
}

#[derive(Debug, Serialize)]
struct CompareOut {
    g_s_expected: f64,
    w0_expected: f64,
    g_s_rel_err: f64,
    w0_rel_err: f64,
}

fn main() {
    let t0 = Instant::now();
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
    let allow_fixtures = parse_flag("--allow-fixtures");
    let first_principles_env = std::env::var_os("CYRUS_FIRST_PRINCIPLES").is_some();

    if allow_fixtures && first_principles_env {
        eprintln!("[ERROR] --allow-fixtures cannot be used with CYRUS_FIRST_PRINCIPLES");
        std::process::exit(2);
    }

    if data_dir.is_none() && !allow_fixtures {
        eprintln!(
            "[ERROR] No McAllister data dir set. Refusing to fall back to JSON fixtures."
        );
        eprintln!("[ERROR] Set CYRUS_MCALLISTER_DATA_DIR or pass --allow-fixtures.");
        std::process::exit(2);
    }
    if data_dir.is_none() && allow_fixtures {
        eprintln!("[MODE] fixtures (JSON)");
        eprintln!("[WARN] using JSON fixtures (not a first-principles run)");
    }
    if data_dir.is_some() {
        eprintln!("[MODE] first-principles (.dat)");
    }
    let cutoff = F64::<Pos>::new(cutoff).expect("cutoff must be positive");
    if let Some(dir) = &data_dir {
        eprintln!("[INFO] using McAllister data dir: {}", dir);
    }

    let (k_raw, m_raw) = if let Some(dir) = &data_dir {
        let dir = PathBuf::from(dir);
        let k = read_csv_i64(&dir.join("K_vec.dat"));
        let m = read_csv_i64(&dir.join("M_vec.dat"));
        (k, m)
    } else {
        let flux = load_flux();
        (flux.k, flux.m)
    };
    let k_typed: Vec<I64<Finite>> = k_raw.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = m_raw.iter().map(|&x| I64::<Finite>::new(x)).collect();

    let (points_raw, heights_raw) = if let Some(dir) = &data_dir {
        let dir = PathBuf::from(dir);
        let points = read_points(&dir.join("points.dat"));
        let heights = read_csv_f64(&dir.join("heights.dat"));
        (points, heights)
    } else {
        let poly_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/mcallister_e2e/inputs/polytope.json");
        let poly: PolytopeInput = load_json(&poly_path);
        let heights_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/mcallister_e2e/inputs/heights.json");
        let heights: HeightsInput = load_json(&heights_path);
        (poly.points, heights.values)
    };

    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("Failed to build polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter triangulation points");
    let triangulation =
        compute_regular_triangulation(&triangulation_points, &heights_raw)
            .expect("Failed to compute triangulation from heights");
    eprintln!("[TIME] triangulation: {:.2?}", t0.elapsed());

    let (_glsm, linrels, basis) =
        compute_glsm_and_linrels(&triangulation_points).expect("Failed GLSM/linrels");
    if let Some(dir) = &data_dir {
        assert_basis_matches_dat(&basis, dir);
    }
    eprintln!("[TIME] glsm/linrels: {:.2?}", t0.elapsed());
    let curve_basis =
        compute_curve_basis_matrix(&linrels, &basis).expect("Failed curve basis matrix");
    eprintln!("[TIME] curve_basis: {:.2?}", t0.elapsed());

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

    let intersection_full = compute_intersection_cytools(
        &triangulation,
        &triangulation_points,
        &linrels_i64,
    )
    .expect("Failed intersection numbers");
    eprintln!("[TIME] intersection: {:.2?}", t0.elapsed());
    let intersection = intersection_in_basis(&intersection_full, &basis);
    eprintln!("[TIME] intersection_in_basis: {:.2?}", t0.elapsed());

    let rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        true,
        false,
        Some(&basis),
    )
    .expect("Failed mori cone cap rays");
    eprintln!(
        "[DEBUG] mori_cone_cap rays: count={}, dim={}",
        rays.len(),
        rays.first().map(Vec::len).unwrap_or(0)
    );
    eprintln!("[TIME] mori_cone_cap: {:.2?}", t0.elapsed());
    let grading = compute_grading_vector(&rays).expect("No grading vector found");
    eprintln!("[TIME] grading_vector: {:.2?}", t0.elapsed());

    let q_matrix: Vec<Vec<i64>> = curve_basis
        .iter()
        .map(|row| {
            row.iter()
                .skip(1)
                .map(|v| i64::try_from(v).expect("q entry fits i64"))
                .collect()
        })
        .collect();

    let invariants = compute_gv_invariants(
        &rays,
        &grading,
        &q_matrix,
        &intersection,
        min_points,
        max_deg,
    )
    .expect("GV computation failed");
    eprintln!("[TIME] gv_invariants: {:.2?}", t0.elapsed());
    eprintln!("[INFO] gv invariants count: {}", invariants.len());

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

    let n_mat = compute_n_matrix(&intersection, &m_typed);
    let Some(p) = solve_linear_system_faer(&n_mat, &k_typed) else {
        eprintln!("[ERROR] N matrix is singular");
        std::process::exit(2);
    };

    let terms = build_racetrack_terms(&gv, &m_typed, &p, cutoff);
    eprintln!("[INFO] racetrack terms: {}", terms.len());
    if terms.len() < 2 {
        eprintln!("[ERROR] not enough racetrack terms");
        std::process::exit(2);
    }

    let Some(rt_res) = solve_racetrack(&terms) else {
        eprintln!("[ERROR] no stable racetrack solution");
        std::process::exit(2);
    };

    let w0 = compute_w0(&rt_res, &terms[0], &terms[1]);
    eprintln!("[RESULT] g_s = {}", rt_res.g_s.get());
    eprintln!("[RESULT] Im(tau) = {}", rt_res.im_tau.get());
    eprintln!("[RESULT] W0 = {}", w0.get());

    let compare_dir = compare_dir.or_else(|| data_dir.clone());
    let compare = if let Some(dir) = compare_dir {
        let dir = PathBuf::from(dir);
        let g_s_path = dir.join("g_s.dat");
        let w0_path = dir.join("W_0.dat");
        let g_s_expected = read_scalar(&g_s_path);
        let w0_expected = read_scalar(&w0_path);
        let g_s_rel_err = ((rt_res.g_s.get() - g_s_expected) / g_s_expected).abs();
        let w0_rel_err = ((w0.get() - w0_expected) / w0_expected).abs();
        eprintln!(
            "[COMPARE] g_s expected={}, rel_err={}",
            g_s_expected, g_s_rel_err
        );
        eprintln!(
            "[COMPARE] W0 expected={}, rel_err={}",
            w0_expected, w0_rel_err
        );
        Some(CompareOut {
            g_s_expected,
            w0_expected,
            g_s_rel_err,
            w0_rel_err,
        })
    } else {
        None
    };

    if let Some(path) = out_path {
        let payload = RacetrackOut {
            g_s: rt_res.g_s.get(),
            im_tau: rt_res.im_tau.get(),
            w0: w0.get(),
            terms: terms
                .iter()
                .map(|t| RacetrackTermOut {
                    exponent: t.exponent.get(),
                    coefficient: t.coefficient.get(),
                })
                .collect(),
            compare,
        };
        let json = serde_json::to_string_pretty(&payload)
            .expect("serialize racetrack output");
        std::fs::write(&path, json).unwrap_or_else(|e| {
            panic!("Failed to write {}: {e}", path);
        });
        eprintln!("[INFO] wrote racetrack output: {}", path);
    }
    eprintln!("[TIME] total: {:.2?}", t0.elapsed());
}
