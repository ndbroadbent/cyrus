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
    Point, Polytope, Triangulation, compute_curve_basis_matrix, compute_frst_heights,
    compute_glsm_and_linrels, compute_grading_vector, compute_gv_invariants,
    compute_intersection_cytools, compute_linear_relations_no_origin, compute_mori_cone_cap_rays,
    curve_basis_matrix_without_origin_i64, intersection_in_basis,
};

const DEFAULT_MCALLISTER_GV_MIN_POINTS: u32 = 20_000;

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
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

fn read_simplices(path: &PathBuf) -> Vec<Vec<usize>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| {
                    s.trim().parse::<usize>().unwrap_or_else(|e| {
                        panic!("Invalid simplex entry {s} in {}: {e}", path.display())
                    })
                })
                .collect::<Vec<usize>>()
        })
        .collect()
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
    let expected_dual_points = {
        let mut points = read_points(&dir.join("dual_points.dat"));
        points.sort();
        points
    };
    let actual_dual_points = sorted_point_coords(dual_polytope.vertices());
    if actual_dual_points != expected_dual_points {
        eprintln!("[ERROR] computed dual polytope differs from dual_points.dat checkpoint");
        std::process::exit(2);
    }

    let mut expected_simplices = read_simplices(&dir.join("dual_simplices.dat"));
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

fn main() {
    let t0 = Instant::now();
    let max_deg = parse_arg_value::<u32>("--max-deg");
    let min_points = if max_deg.is_some() {
        None
    } else {
        parse_arg_value::<u32>("--min-points").or(Some(DEFAULT_MCALLISTER_GV_MIN_POINTS))
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
        eprintln!("[ERROR] No McAllister data dir set. Refusing to fall back to JSON fixtures.");
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
        eprintln!("[INFO] using McAllister data dir: {dir}");
    }

    let primal_points_raw = data_dir.as_ref().map_or_else(
        || {
            let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
            let poly_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
            let poly: PolytopeInput = load_json(&poly_path);
            poly.points
        },
        |dir| {
            let dir = PathBuf::from(dir);
            read_points(&dir.join("points.dat"))
        },
    );

    let primal_points: Vec<Point> = primal_points_raw
        .iter()
        .map(|p| Point::new(p.clone()))
        .collect::<Vec<_>>();
    let primal_polytope = Polytope::from_vertices(primal_points).expect("Failed primal polytope");
    let polytope = primal_polytope
        .compute_dual()
        .expect("Failed dual polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter dual triangulation points");
    let origin_idx = triangulation_points
        .iter()
        .position(|point| point.coords().iter().all(|&coord| coord == 0))
        .expect("dual origin not found");
    let (_heights, triangulation) =
        compute_frst_heights(&triangulation_points, origin_idx).expect("Failed dual FRST");
    if let Some(dir) = &data_dir {
        validate_dual_checkpoint(&polytope, &triangulation, dir);
    }
    eprintln!("[TIME] triangulation: {:.2?}", t0.elapsed());

    let (_glsm, linrels, basis_auto) =
        compute_glsm_and_linrels(&triangulation_points).expect("Failed GLSM/linrels");
    let basis = basis_auto;
    eprintln!("[DEBUG] gv divisor basis: {basis:?}");
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

    let intersection_full =
        compute_intersection_cytools(&triangulation, &triangulation_points, &linrels_i64)
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
        rays.first().map_or(0, Vec::len)
    );
    eprintln!("[TIME] mori_cone_cap: {:.2?}", t0.elapsed());
    let grading = compute_grading_vector(&rays).expect("No grading vector found");
    eprintln!("[TIME] grading_vector: {:.2?}", t0.elapsed());

    let q_matrix = curve_basis_matrix_without_origin_i64(&curve_basis).expect("q entries fit i64");

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
