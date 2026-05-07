//! One-off McAllister GV computation (opt-in).
//!
//! This is intentionally a binary, not a test. It is expensive and should be
//! run explicitly when you want to compute/cache GV invariants.
//!
//! Use `--allow-fixtures` to permit JSON fixture fallback when no data dir is set.

use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

use malachite::Integer;

use cyrus_core::{
    DivisorBasis, Point, Polytope, Triangulation, compute_frst_heights, compute_glsm_and_linrels,
    compute_grading_vector, compute_gv_invariants, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_mori_cone_cap_rays, gv_divisor_basis_data,
    intersection_in_divisor_basis, map_basis_gv_invariants_to_ambient,
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

fn read_integer_vector(path: &PathBuf) -> Vec<Integer> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .replace('\n', ",")
        .split(',')
        .filter(|entry| !entry.trim().is_empty())
        .map(|entry| {
            entry
                .trim()
                .parse::<Integer>()
                .unwrap_or_else(|()| panic!("Invalid integer entry {entry} in {}", path.display()))
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

fn validate_dual_gv_checkpoint(
    computed_basis_gv: &[(Vec<i32>, Integer)],
    curve_basis_matrix: &[Vec<Integer>],
    data_dir: &str,
) {
    let dir = PathBuf::from(data_dir);
    let expected_curves = read_points(&dir.join("dual_curves.dat"));
    let expected_gvs = read_integer_vector(&dir.join("dual_curves_gv.dat"));
    if expected_curves.len() != expected_gvs.len() {
        eprintln!(
            "[ERROR] dual_curves.dat row count {} differs from dual_curves_gv.dat count {}",
            expected_curves.len(),
            expected_gvs.len()
        );
        std::process::exit(2);
    }

    let computed_ambient =
        map_basis_gv_invariants_to_ambient(computed_basis_gv, curve_basis_matrix)
            .unwrap_or_else(|e| panic!("Failed to map GV invariants to ambient basis: {e}"));
    let mut computed = BTreeMap::new();
    for (curve, gv) in computed_ambient {
        if let Some(previous) = computed.insert(curve.clone(), gv.clone()) {
            if previous != gv {
                eprintln!(
                    "[ERROR] computed GV conflict for ambient curve {curve:?}: {previous} vs {gv}"
                );
                std::process::exit(2);
            }
        }
    }

    let mut expected = BTreeMap::new();
    for (curve, gv) in expected_curves.into_iter().zip(expected_gvs) {
        if let Some(previous) = expected.insert(curve.clone(), gv.clone()) {
            if previous != gv {
                eprintln!(
                    "[ERROR] checkpoint GV conflict for ambient curve {curve:?}: {previous} vs {gv}"
                );
                std::process::exit(2);
            }
        }
    }

    let mut matches = 0usize;
    let mut missing = Vec::new();
    let mut mismatches = Vec::new();
    for (curve, expected_gv) in &expected {
        match computed.get(curve) {
            Some(computed_gv) if computed_gv == expected_gv => matches += 1,
            Some(computed_gv) => {
                mismatches.push((curve.clone(), expected_gv.clone(), computed_gv.clone()))
            }
            None => missing.push((curve.clone(), expected_gv.clone())),
        }
    }

    eprintln!(
        "[INFO] dual GV checkpoint: matches={}/{} computed_ambient={}",
        matches,
        expected.len(),
        computed.len()
    );
    if !missing.is_empty() || !mismatches.is_empty() {
        eprintln!(
            "[ERROR] dual GV checkpoint mismatch: missing={} mismatches={}",
            missing.len(),
            mismatches.len()
        );
        if let Some((curve, expected_gv)) = missing.first() {
            eprintln!(
                "[ERROR] first missing checkpoint GV curve {curve:?}: expected {expected_gv}"
            );
        }
        if let Some((curve, expected_gv, computed_gv)) = mismatches.first() {
            eprintln!(
                "[ERROR] first mismatched checkpoint GV curve {curve:?}: expected {expected_gv}, computed {computed_gv}"
            );
        }
        std::process::exit(2);
    }
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
    let intersection =
        intersection_in_divisor_basis(&intersection_full, DivisorBasis::Indices(&basis))
            .expect("Failed intersection basis handoff");
    eprintln!("[TIME] intersection_in_basis: {:.2?}", t0.elapsed());

    let ambient_rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        false,
        false,
        None,
    )
    .expect("Failed ambient mori cone cap rays");
    let gv_basis = gv_divisor_basis_data(&ambient_rays, &linrels, DivisorBasis::Indices(&basis))
        .expect("Failed GV divisor-basis data");
    eprintln!(
        "[DEBUG] mori_cone_cap rays: count={}, dim={}",
        gv_basis.mori_rays.len(),
        gv_basis.mori_rays.first().map_or(0, Vec::len)
    );
    eprintln!("[TIME] mori_cone_cap: {:.2?}", t0.elapsed());
    let grading = compute_grading_vector(&gv_basis.mori_rays).expect("No grading vector found");
    eprintln!("[TIME] grading_vector: {:.2?}", t0.elapsed());

    let invariants = compute_gv_invariants(
        &gv_basis.mori_rays,
        &grading,
        &gv_basis.q_matrix,
        &intersection,
        min_points,
        max_deg,
    )
    .expect("GV computation failed");
    eprintln!("[TIME] gv_invariants: {:.2?}", t0.elapsed());
    eprintln!("[INFO] gv invariants count: {}", invariants.len());
    if let Some(dir) = &data_dir {
        validate_dual_gv_checkpoint(&invariants, &gv_basis.curve_basis_matrix, dir);
    }
}
