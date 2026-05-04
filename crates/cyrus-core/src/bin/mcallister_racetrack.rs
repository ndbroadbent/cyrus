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
    Point, Polytope, Triangulation, basis_change_matrix, build_racetrack_terms,
    compute_curve_basis_matrix, compute_glsm_and_linrels, compute_grading_vector,
    compute_gv_invariants, compute_intersection_cytools, compute_linear_relations_no_origin,
    compute_mori_cone_cap_rays, compute_w0_from_terms, intersection_in_basis, is_unimodular,
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

#[derive(Debug, Deserialize)]
struct SimplicesInput {
    simplices: Vec<Vec<usize>>,
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

    let (points_raw, simplices_raw) = data_dir.as_ref().map_or_else(
        || {
            let poly_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests/mcallister_e2e/overrides/dual_points.json");
            let poly: PolytopeInput = load_json(&poly_path);
            let simplices_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests/mcallister_e2e/overrides/dual_simplices.json");
            let simplices: SimplicesInput = load_json(&simplices_path);
            (poly.points, simplices.simplices)
        },
        |dir| {
            let dir = PathBuf::from(dir);
            let points = read_points(&dir.join("dual_points.dat"));
            let simplices = read_simplices(&dir.join("dual_simplices.dat"));
            (points, simplices)
        },
    );

    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("Failed to build polytope");
    let triangulation_point_count = simplex_point_count(&simplices_raw);
    let triangulation_points = points_raw
        .iter()
        .take(triangulation_point_count)
        .map(|p| Point::new(p.clone()))
        .collect::<Vec<_>>();
    let triangulation = Triangulation::new(simplices_raw);
    eprintln!("[TIME] triangulation: {:.2?}", t0.elapsed());

    let (glsm, linrels, basis_auto) =
        compute_glsm_and_linrels(&triangulation_points).expect("Failed GLSM/linrels");
    let basis = basis_auto;
    eprintln!("[DEBUG] racetrack divisor basis: {basis:?}");
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
    let flux_basis = vec![3, 4, 5, 8];
    let k_raw = transform_k_flux_to_computed_basis(&glsm, &basis, &flux_basis, &k_raw);
    let m_raw = transform_m_flux_to_computed_basis(&glsm, &basis, &flux_basis, &m_raw, "M");
    let k_typed: Vec<I64<Finite>> = k_raw.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = m_raw.iter().map(|&x| I64::<Finite>::new(x)).collect();

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
