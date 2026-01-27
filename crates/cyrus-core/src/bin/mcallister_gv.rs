//! One-off McAllister GV computation (opt-in).
//!
//! This is intentionally a binary, not a test. It is expensive and should be
//! run explicitly when you want to compute/cache GV invariants.
//!
//! Use `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.

use serde::Deserialize;
use std::path::PathBuf;
use std::time::Instant;

use cyrus_core::{
    Point, Polytope, compute_curve_basis_matrix, compute_glsm_and_linrels,
    compute_grading_vector, compute_intersection_cytools, compute_linear_relations_no_origin,
    compute_mori_cone_cap_rays, compute_gv_invariants, compute_regular_triangulation,
    intersection_in_basis,
};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
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

fn main() {
    let t0 = Instant::now();
    let max_deg = parse_arg_value::<u32>("--max-deg");
    let min_points = if max_deg.is_some() {
        None
    } else {
        parse_arg_value::<u32>("--min-points").or(Some(100))
    };
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
    if let Some(dir) = &data_dir {
        eprintln!("[INFO] using McAllister data dir: {}", dir);
    }

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
}
