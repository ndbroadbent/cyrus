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

use cyrus_core::flat_direction::{compute_n_matrix, solve_linear_system_faer};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::{
    DivisorBasis, Point, Polytope, Triangulation, apply_integer_basis_transform,
    apply_integer_basis_transform_transpose, basis_change_matrix, build_racetrack_terms,
    compute_frst_heights, compute_glsm_and_linrels, compute_grading_vector, compute_gv_invariants,
    compute_intersection_cytools, compute_linear_relations_no_origin, compute_mori_cone_cap_rays,
    compute_w0_from_terms, gv_divisor_basis_data, intersection_in_divisor_basis, is_unimodular,
    solve_racetrack,
};

const DEFAULT_MCALLISTER_GV_MIN_POINTS: u32 = 20_000;

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
    apply_integer_basis_transform(&transform, values, label).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to apply {label} basis transform: {e}");
        std::process::exit(2);
    })
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
    apply_integer_basis_transform_transpose(&transform, values, "K").unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to apply K basis transform: {e}");
        std::process::exit(2);
    })
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

fn mirror_geometry_from_primal_points(
    primal_points_raw: Vec<Vec<i64>>,
) -> (Polytope, Vec<Point>, Triangulation) {
    let primal_points = primal_points_raw
        .into_iter()
        .map(Point::new)
        .collect::<Vec<_>>();
    let primal_polytope =
        Polytope::from_vertices(primal_points).expect("Failed to build primal polytope");
    let dual_polytope = primal_polytope
        .compute_dual()
        .expect("Failed to compute dual polytope");
    let dual_points = dual_polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter dual triangulation points");
    let dual_origin_idx = dual_points
        .iter()
        .position(|point| point.coords().iter().all(|&coord| coord == 0))
        .expect("dual origin not found in triangulation points");
    let (_dual_heights, dual_triangulation) =
        compute_frst_heights(&dual_points, dual_origin_idx).expect("Failed to compute dual FRST");
    (dual_polytope, dual_points, dual_triangulation)
}

fn sorted_point_coords(points: &[Point]) -> Vec<Vec<i64>> {
    let mut coords = points
        .iter()
        .map(|point| point.coords().to_vec())
        .collect::<Vec<_>>();
    coords.sort();
    coords
}

fn sorted_simplices(triangulation: &Triangulation) -> Vec<Vec<usize>> {
    let mut simplices = triangulation.simplices().to_vec();
    simplices.sort();
    simplices
}

fn validate_dual_geometry_checkpoint(
    dual_polytope: &Polytope,
    dual_triangulation: &Triangulation,
    data_dir: &str,
) {
    let dir = PathBuf::from(data_dir);
    let dual_points_path = dir.join("dual_points.dat");
    if dual_points_path.exists() {
        let mut expected_points = read_points(&dual_points_path);
        expected_points.sort();
        let actual_points = sorted_point_coords(dual_polytope.vertices());
        if actual_points != expected_points {
            eprintln!("[ERROR] computed dual polytope differs from dual_points.dat checkpoint");
            std::process::exit(2);
        }
        eprintln!("[INFO] computed dual polytope matches dual_points.dat checkpoint");
    }

    let dual_simplices_path = dir.join("dual_simplices.dat");
    if dual_simplices_path.exists() {
        let mut expected_simplices = read_points(&dual_simplices_path)
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|value| {
                        usize::try_from(value)
                            .expect("dual_simplices.dat index must be non-negative")
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        expected_simplices.sort();
        let actual_simplices = sorted_simplices(dual_triangulation);
        if actual_simplices != expected_simplices {
            eprintln!("[ERROR] computed dual FRST differs from dual_simplices.dat checkpoint");
            std::process::exit(2);
        }
        eprintln!("[INFO] computed dual FRST matches dual_simplices.dat checkpoint");
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
    let cutoff = F64::<Pos>::new(cutoff).expect("cutoff must be positive");
    if let Some(dir) = &data_dir {
        eprintln!("[INFO] using McAllister data dir: {dir}");
    }

    let (k_raw, m_raw) = data_dir.as_ref().map_or_else(
        || {
            let flux = load_flux();
            (flux.k, flux.m)
        },
        |dir| {
            let dir = PathBuf::from(dir);
            let k = read_csv_i64(&dir.join("K_vec.dat"));
            let m = read_csv_i64(&dir.join("M_vec.dat"));
            (k, m)
        },
    );

    let primal_points_raw = data_dir.as_ref().map_or_else(
        || {
            let poly_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests/mcallister_e2e/inputs/polytope.json");
            let poly: PolytopeInput = load_json(&poly_path);
            poly.points
        },
        |dir| {
            let dir = PathBuf::from(dir);
            read_points(&dir.join("points.dat"))
        },
    );
    let (polytope, triangulation_points, triangulation) =
        mirror_geometry_from_primal_points(primal_points_raw);
    if let Some(dir) = &data_dir {
        validate_dual_geometry_checkpoint(&polytope, &triangulation, dir);
    }
    eprintln!("[TIME] mirror geometry: {:.2?}", t0.elapsed());

    let (glsm, linrels, basis_auto) =
        compute_glsm_and_linrels(&triangulation_points).expect("Failed GLSM/linrels");
    let basis = basis_auto;
    eprintln!("[DEBUG] racetrack divisor basis: {basis:?}");
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
    let flux_basis = vec![3, 4, 5, 8];
    let k_raw = transform_k_flux_to_computed_basis(&glsm, &basis, &flux_basis, &k_raw);
    let m_raw = transform_m_flux_to_computed_basis(&glsm, &basis, &flux_basis, &m_raw, "M");
    let k_typed: Vec<I64<Finite>> = k_raw.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = m_raw.iter().map(|&x| I64::<Finite>::new(x)).collect();

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

    let Some(w0) = compute_w0_from_terms(&rt_res, &terms) else {
        eprintln!("[ERROR] racetrack W0 computation failed or cancelled exactly");
        std::process::exit(2);
    };
    eprintln!("[RESULT] g_s = {}", rt_res.g_s.get());
    eprintln!("[RESULT] Im(tau) = {}", rt_res.im_tau.get());
    eprintln!("[RESULT] W0 = {}", w0.get());

    let compare_dir = compare_dir.or_else(|| data_dir.clone());
    let compare = compare_dir.map_or(None, |dir| {
        let dir = PathBuf::from(dir);
        let g_s_path = dir.join("g_s.dat");
        let w0_path = dir.join("W_0.dat");
        let g_s_expected = read_scalar(&g_s_path);
        let w0_expected = read_scalar(&w0_path);
        let g_s_rel_err = ((rt_res.g_s.get() - g_s_expected) / g_s_expected).abs();
        let w0_rel_err = ((w0.get() - w0_expected) / w0_expected).abs();
        eprintln!("[COMPARE] g_s expected={g_s_expected}, rel_err={g_s_rel_err}");
        eprintln!("[COMPARE] W0 expected={w0_expected}, rel_err={w0_rel_err}");
        Some(CompareOut {
            g_s_expected,
            w0_expected,
            g_s_rel_err,
            w0_rel_err,
        })
    });

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
        let json = serde_json::to_string_pretty(&payload).expect("serialize racetrack output");
        std::fs::write(&path, json).unwrap_or_else(|e| {
            panic!("Failed to write {path}: {e}");
        });
        eprintln!("[INFO] wrote racetrack output: {path}");
    }
    eprintln!("[TIME] total: {:.2?}", t0.elapsed());
}
