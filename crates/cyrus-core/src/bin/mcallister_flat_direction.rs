#![allow(dead_code)] // diagnostic helpers are retained for ad-hoc investigation runs

//! Flat-direction diagnostics for McAllister data.
//!
//! Computes the dual intersection numbers, builds the N matrix,
//! solves for p, and reports the sign and magnitude of κ p^3.
//! This is a focused diagnostic to isolate ek0 sign issues.
//!
//! Usage (examples):
//! - Basic run (requires CYRUS_MCALLISTER_DATA_DIR unless --allow-fixtures is set):
//!   cargo run -p cyrus-core --bin mcallister_flat_direction --release
//! - Compute basis automatically from dual points (default unless --basis-* provided)
//! - Compare two bases:
//!   --basis-indices "3,4,5,8" --compare-basis-indices "1,2,3,4"
//! - Deterministic swap:
//!   --sweep-pairs "3:1;4:2"
//! - Two-swap bounded sweep + report:
//!   --sweep-two-swaps --sweep-two-max-attempts 5000 --sweep-two-report /tmp/flat_sweep.json

use serde::Deserialize;
use serde_json::json;
use std::path::{Path, PathBuf};
use std::time::Instant;

use cyrus_core::flat_direction::{
    compute_flat_direction, compute_n_matrix, solve_linear_system_faer,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    Intersection, Point, Polytope, Triangulation, compute_frst_heights, compute_glsm_and_linrels,
    compute_intersection_cytools, compute_linear_relations_no_origin, intersection_in_basis,
};

#[derive(Debug, Deserialize)]
struct FluxInput {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
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

fn load_json<T: for<'de> Deserialize<'de>>(path: &PathBuf) -> T {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&content)
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

fn read_csv_usize(path: &PathBuf) -> Vec<usize> {
    read_csv_i64(path)
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index must be non-negative"))
        .collect()
}

fn parse_basis_indices(arg: &str) -> Vec<usize> {
    arg.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<usize>()
                .unwrap_or_else(|e| panic!("Invalid basis index {s}: {e}"))
        })
        .collect()
}

fn parse_usize_list(arg: &str) -> Vec<usize> {
    arg.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<usize>()
                .unwrap_or_else(|e| panic!("Invalid value {s}: {e}"))
        })
        .collect()
}

fn parse_pair_list(arg: &str) -> Vec<(usize, usize)> {
    arg.split(';')
        .map(str::trim)
        .filter(|pair| !pair.is_empty())
        .map(|pair| {
            let mut parts = pair.split(':');
            let a = parts
                .next()
                .unwrap_or("0")
                .trim()
                .parse::<usize>()
                .unwrap_or_else(|e| panic!("Invalid pair {pair}: {e}"));
            let b = parts
                .next()
                .unwrap_or("0")
                .trim()
                .parse::<usize>()
                .unwrap_or_else(|e| panic!("Invalid pair {pair}: {e}"));
            (a, b)
        })
        .collect()
}

fn apply_swaps(mut basis: Vec<usize>, swaps: &[(usize, usize)]) -> Option<Vec<usize>> {
    for &(out_idx, in_idx) in swaps {
        if basis.contains(&in_idx) {
            return None;
        }
        if let Some(pos) = basis.iter().position(|&v| v == out_idx) {
            basis[pos] = in_idx;
        } else {
            return None;
        }
    }
    Some(basis)
}

#[derive(Debug)]
struct TermSummary {
    pos_sum: f64,
    neg_sum: f64,
    total: f64,
    terms: Vec<(f64, (usize, usize, usize))>,
}

fn summarize_terms(kappa: &cyrus_core::Intersection, p: &[F64<Finite>]) -> TermSummary {
    let mut terms: Vec<(f64, (usize, usize, usize))> = Vec::new();
    let mut pos_sum = 0.0f64;
    let mut neg_sum = 0.0f64;
    for (&(i, j, k), val) in kappa.iter() {
        let (kappa_val, _): (f64, _) =
            malachite::num::conversion::traits::RoundingFrom::rounding_from(
                val.get(),
                malachite::rounding_modes::RoundingMode::Nearest,
            );
        let t_product = p[i].get() * p[j].get() * p[k].get();
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        let contrib = mult * kappa_val * t_product;
        if contrib >= 0.0 {
            pos_sum += contrib;
        } else {
            neg_sum += contrib;
        }
        terms.push((contrib, (i, j, k)));
    }
    terms.sort_by(|a, b| b.0.abs().partial_cmp(&a.0.abs()).unwrap());
    TermSummary {
        pos_sum,
        neg_sum,
        total: pos_sum + neg_sum,
        terms,
    }
}

fn print_terms(label: &str, summary: &TermSummary, top_terms: usize) {
    eprintln!(
        "[INFO] {label} contributions: pos_sum={:.6e}, neg_sum={:.6e}, total={:.6e}",
        summary.pos_sum, summary.neg_sum, summary.total
    );
    eprintln!("[INFO] {label} top κ p^3 contributions:");
    for (idx, (contrib, key)) in summary.terms.iter().take(top_terms).enumerate() {
        eprintln!("  {:>2}. κ_{:?} = {:>14.6e}", idx + 1, key, contrib);
    }
}

fn summary_json(summary: &TermSummary, top_terms: usize) -> serde_json::Value {
    let top: Vec<serde_json::Value> = summary
        .terms
        .iter()
        .take(top_terms)
        .map(|(contrib, key)| {
            json!({
                "index": [key.0, key.1, key.2],
                "contrib": contrib,
            })
        })
        .collect();
    json!({
        "pos_sum": summary.pos_sum,
        "neg_sum": summary.neg_sum,
        "total": summary.total,
        "top_terms": top,
    })
}

struct FlatArgs {
    data_dir: Option<String>,
    allow_fixtures: bool,
    top_terms: usize,
    basis_indices_arg: Option<String>,
    basis_file_arg: Option<String>,
    compare_basis_indices_arg: Option<String>,
    compare_basis_file_arg: Option<String>,
    allow_nonpositive: bool,
    sweep_swaps: usize,
    sweep_replace_from: Option<String>,
    sweep_replace_to: Option<String>,
    sweep_pairs_arg: Option<String>,
    sweep_max_attempts: usize,
    sweep_report_path: Option<String>,
    sweep_two: bool,
    sweep_two_max_attempts: usize,
    sweep_two_report_path: Option<String>,
    out_path: Option<String>,
}

fn parse_args() -> FlatArgs {
    let data_dir = parse_arg_value::<String>("--data-dir")
        .or_else(|| std::env::var("CYRUS_MCALLISTER_DATA_DIR").ok());
    FlatArgs {
        data_dir,
        allow_fixtures: parse_flag("--allow-fixtures"),
        top_terms: parse_arg_value::<usize>("--top-terms").unwrap_or(12),
        basis_indices_arg: parse_arg_value::<String>("--basis-indices"),
        basis_file_arg: parse_arg_value::<String>("--basis-file"),
        compare_basis_indices_arg: parse_arg_value::<String>("--compare-basis-indices"),
        compare_basis_file_arg: parse_arg_value::<String>("--compare-basis-file"),
        allow_nonpositive: parse_flag("--allow-nonpositive"),
        sweep_swaps: parse_arg_value::<usize>("--sweep-swaps").unwrap_or(0),
        sweep_replace_from: parse_arg_value::<String>("--sweep-replace-from"),
        sweep_replace_to: parse_arg_value::<String>("--sweep-replace-to"),
        sweep_pairs_arg: parse_arg_value::<String>("--sweep-pairs"),
        sweep_max_attempts: parse_arg_value::<usize>("--sweep-max-attempts").unwrap_or(5000),
        sweep_report_path: parse_arg_value::<String>("--sweep-report"),
        sweep_two: parse_flag("--sweep-two-swaps"),
        sweep_two_max_attempts: parse_arg_value::<usize>("--sweep-two-max-attempts")
            .unwrap_or(10000),
        sweep_two_report_path: parse_arg_value::<String>("--sweep-two-report"),
        out_path: parse_arg_value::<String>("--out"),
    }
}

fn enforce_modes(data_dir: Option<&str>, allow_fixtures: bool) {
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
}

fn load_primal_points(data_dir: Option<&str>, manifest_dir: &Path) -> Vec<Vec<i64>> {
    data_dir.map_or_else(
        || {
            let poly_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
            let poly: PolytopeInput = load_json(&poly_path);
            poly.points
        },
        |dir| {
            let dir = PathBuf::from(dir);
            read_points(&dir.join("points.dat"))
        },
    )
}

fn compute_dual_inputs_from_primal_points(
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
    let dual_points_vec = dual_polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter dual triangulation points");
    let dual_origin_idx = dual_points_vec
        .iter()
        .position(|point| point.coords().iter().all(|&coord| coord == 0))
        .expect("dual origin not found in triangulation points");
    let (_dual_heights, dual_triangulation) =
        compute_frst_heights(&dual_points_vec, dual_origin_idx).expect("Failed dual FRST");
    (dual_polytope, dual_points_vec, dual_triangulation)
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

fn validate_dual_checkpoint(
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

fn compute_dual_kappa_full(
    dual_points_vec: &[Point],
    dual_triangulation: &Triangulation,
) -> Intersection {
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
    compute_intersection_cytools(dual_triangulation, dual_points_vec, &dual_linrels_i64)
        .expect("Failed dual intersection numbers")
}

fn resolve_basis_indices(
    basis_indices_arg: Option<&String>,
    basis_file_arg: Option<&String>,
    dual_points_vec: &[Point],
) -> Vec<usize> {
    basis_indices_arg.map_or_else(
        || {
            if let Some(file) = basis_file_arg {
                read_csv_usize(&PathBuf::from(file))
            } else {
                let (_glsm, _linrel, basis) =
                    compute_glsm_and_linrels(dual_points_vec).expect("Failed dual GLSM");
                basis
            }
        },
        |indices| parse_basis_indices(indices),
    )
}

fn load_flux_vectors(
    data_dir: Option<&str>,
    manifest_dir: &Path,
) -> (Vec<I64<Finite>>, Vec<I64<Finite>>) {
    let (k_raw, m_raw) = data_dir.map_or_else(
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
    );
    let k_flux: Vec<I64<Finite>> = k_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let m_flux: Vec<I64<Finite>> = m_raw.iter().map(|&v| I64::<Finite>::new(v)).collect();
    (k_flux, m_flux)
}

fn contract_ppp(kappa: &Intersection, p: &[F64<Finite>]) -> f64 {
    kappa.contract_triple_finite(p).map_or(0.0, F64::get)
}

fn main() {
    let t0 = Instant::now();
    let data_dir = parse_arg_value::<String>("--data-dir")
        .or_else(|| std::env::var("CYRUS_MCALLISTER_DATA_DIR").ok());
    let allow_fixtures = parse_flag("--allow-fixtures");
    let first_principles_env = std::env::var_os("CYRUS_FIRST_PRINCIPLES").is_some();
    let top_terms = parse_arg_value::<usize>("--top-terms").unwrap_or(12);
    let basis_indices_arg = parse_arg_value::<String>("--basis-indices");
    let basis_file_arg = parse_arg_value::<String>("--basis-file");
    let compare_basis_indices_arg = parse_arg_value::<String>("--compare-basis-indices");
    let compare_basis_file_arg = parse_arg_value::<String>("--compare-basis-file");
    let allow_nonpositive = parse_flag("--allow-nonpositive");
    let sweep_swaps = parse_arg_value::<usize>("--sweep-swaps").unwrap_or(0);
    let sweep_replace_from = parse_arg_value::<String>("--sweep-replace-from");
    let sweep_replace_to = parse_arg_value::<String>("--sweep-replace-to");
    let sweep_pairs_arg = parse_arg_value::<String>("--sweep-pairs");
    let sweep_max_attempts = parse_arg_value::<usize>("--sweep-max-attempts").unwrap_or(5000);
    let sweep_report_path = parse_arg_value::<String>("--sweep-report");
    let sweep_two = parse_flag("--sweep-two-swaps");
    let sweep_two_max_attempts =
        parse_arg_value::<usize>("--sweep-two-max-attempts").unwrap_or(10000);

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
    let sweep_two_report_path = parse_arg_value::<String>("--sweep-two-report");
    let out_path = parse_arg_value::<String>("--out");
    let mut compare_report: Option<serde_json::Value> = None;
    let mut sweep_report: Option<serde_json::Value> = None;
    let mut sweep_two_report: Option<serde_json::Value> = None;
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let primal_points_raw = load_primal_points(data_dir.as_deref(), &manifest_dir);
    let (dual_polytope, dual_points_vec, dual_triangulation) =
        compute_dual_inputs_from_primal_points(primal_points_raw);
    if let Some(dir) = &data_dir {
        validate_dual_checkpoint(&dual_polytope, &dual_triangulation, dir);
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

    let basis_indices = if let Some(indices) = basis_indices_arg.as_ref() {
        parse_basis_indices(indices)
    } else if let Some(file) = basis_file_arg.as_ref() {
        read_csv_usize(&PathBuf::from(file))
    } else {
        let (_glsm, _linrel, basis) =
            compute_glsm_and_linrels(&dual_points_vec).expect("Failed dual GLSM");
        basis
    };

    let dual_kappa = intersection_in_basis(&dual_kappa_full, &basis_indices);

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

    let compute_ppp = |basis: &Vec<usize>| -> Option<f64> {
        let kappa_cmp = intersection_in_basis(&dual_kappa_full, basis);
        let n_cmp = compute_n_matrix(&kappa_cmp, &m_flux);
        let p_cmp = solve_linear_system_faer(&n_cmp, &k_flux)?;
        kappa_cmp
            .contract_triple_finite(&p_cmp)
            .map(cyrus_core::F64::get)
    };

    let n_mat = compute_n_matrix(&dual_kappa, &m_flux);
    let p = solve_linear_system_faer(&n_mat, &k_flux).expect("N matrix solve failed");
    let p_alt =
        compute_flat_direction(&dual_kappa, &k_flux, &m_flux).expect("flat direction failed");

    let kappa_ppp = dual_kappa
        .contract_triple_finite(&p)
        .map_or(0.0, cyrus_core::F64::get);

    let kappa_ppp_alt = dual_kappa
        .contract_triple_finite(&p_alt)
        .map_or(0.0, cyrus_core::F64::get);

    eprintln!("[TIME] flat-direction diagnostics: {:.2?}", t0.elapsed());
    eprintln!("[INFO] basis indices: {basis_indices:?}");
    eprintln!("[INFO] kappa dim: {}", dual_kappa.dim());
    eprintln!("[INFO] kappa nonzero: {}", dual_kappa.num_nonzero());
    eprintln!(
        "[INFO] p (first 4): {:?}",
        &p.iter().take(4).map(|v| v.get()).collect::<Vec<_>>()
    );
    eprintln!("[INFO] kappa ppp: {kappa_ppp}");
    eprintln!("[INFO] kappa ppp alt: {kappa_ppp_alt}");

    let summary_p = summarize_terms(&dual_kappa, &p);
    let summary_p_alt = summarize_terms(&dual_kappa, &p_alt);
    print_terms("p", &summary_p, top_terms);
    print_terms("p_alt", &summary_p_alt, top_terms);

    if kappa_ppp <= 0.0 && !allow_nonpositive {
        eprintln!("[ERROR] kappa p^3 is non-positive");
        std::process::exit(2);
    }

    let ek0 = F64::<Finite>::new(1.0 / (4.0 / 3.0 * kappa_ppp)).expect("finite ek0");
    let ek0_expected = 0.234_393_f64;
    eprintln!("[RESULT] ek0 = {}", ek0.get());
    eprintln!(
        "[COMPARE] ek0 expected={} rel_err={}",
        ek0_expected,
        ((ek0.get() - ek0_expected) / ek0_expected).abs()
    );

    // Optional basis comparison
    let compare_basis = if let Some(indices) = compare_basis_indices_arg.as_ref() {
        Some(parse_basis_indices(indices))
    } else {
        compare_basis_file_arg
            .as_ref()
            .map(|file| read_csv_usize(&PathBuf::from(file)))
    };

    if let Some(compare_indices) = compare_basis {
        if compare_indices.len() != basis_indices.len() {
            eprintln!(
                "[ERROR] compare basis length {} does not match primary basis length {}",
                compare_indices.len(),
                basis_indices.len()
            );
            std::process::exit(2);
        }
        let kappa_cmp = intersection_in_basis(&dual_kappa_full, &compare_indices);
        let n_cmp = compute_n_matrix(&kappa_cmp, &m_flux);
        let p_cmp =
            solve_linear_system_faer(&n_cmp, &k_flux).expect("N matrix solve failed (compare)");
        let kappa_ppp_cmp = kappa_cmp
            .contract_triple_finite(&p_cmp)
            .map_or(0.0, cyrus_core::F64::get);
        eprintln!("[INFO] compare basis indices: {compare_indices:?}");
        eprintln!("[INFO] compare kappa ppp: {kappa_ppp_cmp}");

        let summary_cmp = summarize_terms(&kappa_cmp, &p_cmp);
        print_terms("compare", &summary_cmp, top_terms);

        let mut diffs: Vec<(f64, (usize, usize, usize))> = Vec::new();
        for (contrib, key) in summary_p.terms.iter().take(top_terms) {
            let cmp_val = summary_cmp
                .terms
                .iter()
                .find(|(_, k)| k == key)
                .map_or(0.0, |(v, _)| *v);
            diffs.push((contrib - cmp_val, *key));
        }
        diffs.sort_by(|a, b| b.0.abs().partial_cmp(&a.0.abs()).unwrap());
        eprintln!("[INFO] top term diffs (primary - compare):");
        for (idx, (diff, key)) in diffs.iter().take(top_terms).enumerate() {
            eprintln!("  {:>2}. Δκ_{:?} = {:>14.6e}", idx + 1, key, diff);
        }

        compare_report = Some(json!({
            "basis": compare_indices,
            "kappa_ppp": kappa_ppp_cmp,
            "summary": summary_json(&summary_cmp, top_terms),
            "diffs": diffs.iter().take(top_terms).map(|(d, k)| {
                json!({
                    "index": [k.0, k.1, k.2],
                    "diff": d,
                })
            }).collect::<Vec<_>>(),
        }));
    }

    // Optional sweep: explicit swap pairs (deterministic)
    if let Some(pairs_arg) = sweep_pairs_arg.as_ref() {
        let pairs = parse_pair_list(pairs_arg);
        if pairs.is_empty() {
            eprintln!("[ERROR] sweep pairs list is empty");
            std::process::exit(2);
        }
        if let Some(candidate) = apply_swaps(basis_indices.clone(), &pairs) {
            let Some(kappa_ppp_cmp) = compute_ppp(&candidate) else {
                eprintln!("[ERROR] sweep pairs: N matrix solve failed");
                std::process::exit(2);
            };
            eprintln!("[INFO] sweep pairs: kappa ppp = {kappa_ppp_cmp}");
            if kappa_ppp_cmp > 0.0 {
                eprintln!("[FOUND] sweep pairs yield positive kappa p^3");
            }
        } else {
            eprintln!("[ERROR] sweep pairs invalid (duplicate or missing indices)");
            std::process::exit(2);
        }
    }

    // Optional sweep: try swapping a small number of indices from basis to candidate pool.
    if sweep_swaps > 0 {
        let replace_from = sweep_replace_from
            .as_ref()
            .map_or_else(|| basis_indices.clone(), |s| parse_usize_list(s));
        let replace_to = sweep_replace_to
            .as_ref()
            .map_or_else(|| basis_indices.clone(), |s| parse_usize_list(s));

        if replace_from.is_empty() || replace_to.is_empty() {
            eprintln!("[ERROR] sweep replace lists cannot be empty");
            std::process::exit(2);
        }

        eprintln!(
            "[INFO] sweep: swaps={}, from_pool={}, to_pool={}",
            sweep_swaps,
            replace_from.len(),
            replace_to.len()
        );

        let mut attempts = 0usize;
        let mut found_swap: Option<(usize, usize, f64, usize)> = None;
        let mut found = false;
        for &out_idx in replace_from.iter().take(sweep_swaps) {
            for &in_idx in replace_to.iter().take(sweep_swaps) {
                if out_idx == in_idx {
                    continue;
                }
                let mut candidate = basis_indices.clone();
                if !candidate.contains(&in_idx) {
                    if let Some(pos) = candidate.iter().position(|&v| v == out_idx) {
                        candidate[pos] = in_idx;
                    } else {
                        continue;
                    }
                }

                attempts += 1;
                let Some(kappa_ppp_cmp) = compute_ppp(&candidate) else {
                    continue;
                };
                if kappa_ppp_cmp > 0.0 {
                    eprintln!(
                        "[FOUND] swap {out_idx} -> {in_idx} yields kappa ppp = {kappa_ppp_cmp} (attempt {attempts})"
                    );
                    found_swap = Some((out_idx, in_idx, kappa_ppp_cmp, attempts));
                    if let Some(path) = sweep_report_path.as_ref() {
                        let report = json!({
                            "swap": [out_idx, in_idx],
                            "kappa_ppp": kappa_ppp_cmp,
                            "attempts": attempts,
                        });
                        std::fs::write(path, report.to_string())
                            .unwrap_or_else(|e| panic!("Failed to write {path}: {e}"));
                    }
                    found = true;
                    break;
                }
                if attempts >= sweep_max_attempts {
                    eprintln!("[INFO] sweep max attempts reached ({sweep_max_attempts})");
                    break;
                }
            }
            if found {
                break;
            }
            if attempts >= sweep_max_attempts {
                break;
            }
        }
        if !found {
            eprintln!("[INFO] sweep completed: no positive kappa p^3 found (attempts={attempts})");
        }

        sweep_report = Some(json!({
            "attempts": attempts,
            "found": found_swap.map(|(out_idx, in_idx, ppp, tries)| {
                json!({
                    "swap": [out_idx, in_idx],
                    "kappa_ppp": ppp,
                    "attempts": tries,
                })
            }),
        }));
    }

    // Optional sweep: two-swap exhaustive (bounded)
    if sweep_two {
        let replace_from = sweep_replace_from
            .as_ref()
            .map_or_else(|| basis_indices.clone(), |s| parse_usize_list(s));
        let replace_to = sweep_replace_to
            .as_ref()
            .map_or_else(|| basis_indices.clone(), |s| parse_usize_list(s));
        if replace_from.is_empty() || replace_to.is_empty() {
            eprintln!("[ERROR] sweep-two replace lists cannot be empty");
            std::process::exit(2);
        }

        #[derive(Clone, Debug)]
        struct SweepHit {
            swap1: (usize, usize),
            swap2: (usize, usize),
            kappa_ppp: f64,
        }

        let mut attempts = 0usize;
        let mut hits: Vec<SweepHit> = Vec::new();
        for &out1 in &replace_from {
            for &in1 in &replace_to {
                if out1 == in1 {
                    continue;
                }
                for &out2 in &replace_from {
                    if out2 == out1 {
                        continue;
                    }
                    for &in2 in &replace_to {
                        if in2 == in1 || in2 == out2 {
                            continue;
                        }
                        let swaps = [(out1, in1), (out2, in2)];
                        let Some(candidate) = apply_swaps(basis_indices.clone(), &swaps) else {
                            continue;
                        };
                        attempts += 1;
                        let Some(kappa_ppp_cmp) = compute_ppp(&candidate) else {
                            continue;
                        };
                        hits.push(SweepHit {
                            swap1: (out1, in1),
                            swap2: (out2, in2),
                            kappa_ppp: kappa_ppp_cmp,
                        });
                        if attempts >= sweep_two_max_attempts {
                            break;
                        }
                    }
                    if attempts >= sweep_two_max_attempts {
                        break;
                    }
                }
                if attempts >= sweep_two_max_attempts {
                    break;
                }
            }
            if attempts >= sweep_two_max_attempts {
                break;
            }
        }

        hits.sort_by(|a, b| b.kappa_ppp.partial_cmp(&a.kappa_ppp).unwrap());
        if let Some(best) = hits.first() {
            eprintln!(
                "[INFO] sweep-two best: swap {:?} + {:?} => kappa ppp {} (attempts={})",
                best.swap1, best.swap2, best.kappa_ppp, attempts
            );
        } else {
            eprintln!("[INFO] sweep-two: no viable candidates (attempts={attempts})");
        }
        if let Some(path) = sweep_two_report_path.as_ref() {
            let report = json!({
                "attempts": attempts,
                "best": hits.first().map(|h| json!({
                    "swap1": [h.swap1.0, h.swap1.1],
                    "swap2": [h.swap2.0, h.swap2.1],
                    "kappa_ppp": h.kappa_ppp,
                })),
                "top": hits.iter().take(5).map(|h| json!({
                    "swap1": [h.swap1.0, h.swap1.1],
                    "swap2": [h.swap2.0, h.swap2.1],
                    "kappa_ppp": h.kappa_ppp,
                })).collect::<Vec<_>>(),
            });
            std::fs::write(path, report.to_string())
                .unwrap_or_else(|e| panic!("Failed to write {path}: {e}"));
        }
        sweep_two_report = Some(json!({
            "attempts": attempts,
            "best": hits.first().map(|h| json!({
                "swap1": [h.swap1.0, h.swap1.1],
                "swap2": [h.swap2.0, h.swap2.1],
                "kappa_ppp": h.kappa_ppp,
            })),
            "top": hits.iter().take(5).map(|h| json!({
                "swap1": [h.swap1.0, h.swap1.1],
                "swap2": [h.swap2.0, h.swap2.1],
                "kappa_ppp": h.kappa_ppp,
            })).collect::<Vec<_>>(),
        }));
    }

    if let Some(path) = out_path.as_ref() {
        let report = json!({
            "basis": basis_indices,
            "kappa_ppp": kappa_ppp,
            "kappa_ppp_alt": kappa_ppp_alt,
            "ek0": ek0.get(),
            "summary": {
                "p": summary_json(&summary_p, top_terms),
                "p_alt": summary_json(&summary_p_alt, top_terms),
            },
            "compare": compare_report,
            "sweep": sweep_report,
            "sweep_two": sweep_two_report,
        });
        std::fs::write(path, report.to_string())
            .unwrap_or_else(|e| panic!("Failed to write {path}: {e}"));
    }
}
