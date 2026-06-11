//! Landscape smoke test worker: compute chamber-fixed geometry objects for
//! one polytope and dump them as JSON for comparison against CYTools.
//!
//! Input: `--points <csv>` (full CYTools-ordered lattice point list) and
//! optionally `--heights <csv>` (heights over the not-facet points, fixing
//! the chamber; typically CYTools' default triangulation heights so both
//! sides compute in the same phase). Output: `--out <json>`.
//!
//! This is deliberately a thin driver over the same cyrus-core entry points
//! the McAllister runner uses, so a comparison failure here indicts the
//! production pipeline, not a parallel implementation.

use std::collections::BTreeMap;
use std::path::PathBuf;

use cyrus_core::{
    Point, Polytope, compute_frst_heights, compute_glsm_and_linrels, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_regular_triangulation, intersection_in_basis,
};

fn parse_arg_value(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|arg| arg == name)
        .and_then(|idx| args.get(idx + 1).cloned())
}

fn read_csv_rows(path: &str) -> Vec<Vec<i64>> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    text.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|cell| {
                    cell.trim()
                        .parse::<i64>()
                        .unwrap_or_else(|e| panic!("bad integer {cell:?} in {path}: {e}"))
                })
                .collect()
        })
        .collect()
}

fn read_csv_f64(path: &str) -> Vec<f64> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    text.split([',', '\n'])
        .filter(|cell| !cell.trim().is_empty())
        .map(|cell| {
            cell.trim()
                .parse::<f64>()
                .unwrap_or_else(|e| panic!("bad float {cell:?} in {path}: {e}"))
        })
        .collect()
}

fn sorted_simplices(simplices: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = simplices
        .iter()
        .map(|simplex| {
            let mut simplex = simplex.clone();
            simplex.sort_unstable();
            simplex
        })
        .collect();
    out.sort();
    out
}

fn main() {
    let points_path = parse_arg_value("--points").expect("--points <csv> is required");
    let out_path = parse_arg_value("--out").expect("--out <json> is required");
    let heights_path = parse_arg_value("--heights");

    let points_raw = read_csv_rows(&points_path);
    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("failed to build polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter triangulation points");
    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("origin must be among triangulation points");

    let dual = polytope.compute_dual().expect("failed to compute dual");
    let mut dual_not_facet: Vec<Vec<i64>> = dual
        .points_not_interior_to_facets()
        .expect("failed to filter dual points")
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    dual_not_facet.sort();

    // Cyrus's own default chamber (chamber-agreement metric).
    let (_, default_triangulation) = compute_frst_heights(&triangulation_points, origin_idx)
        .expect("failed to compute default FRST");

    // The comparison chamber: provided heights (CYTools default) if given.
    let triangulation = heights_path.map_or_else(
        || default_triangulation.clone(),
        |path| {
            let heights = read_csv_f64(&path);
            assert_eq!(
                heights.len(),
                triangulation_points.len(),
                "heights length must match triangulation point count"
            );
            compute_regular_triangulation(&triangulation_points, &heights)
                .expect("failed to triangulate with provided heights")
        },
    );

    let (_, _, basis) =
        compute_glsm_and_linrels(&triangulation_points).expect("failed GLSM/linrels");
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
    let kappa_full =
        compute_intersection_cytools(&triangulation, &triangulation_points, &linrels_i64)
            .expect("failed intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);

    // Exact rational entries of the in-basis tensor, keyed "i,j,k" sorted.
    let mut kappa_entries: BTreeMap<String, String> = BTreeMap::new();
    let dim = basis.len();
    for i in 0..dim {
        for j in i..dim {
            for k in j..dim {
                let value = kappa_basis.get(i, j, k).to_string();
                if value != "0" {
                    kappa_entries.insert(format!("{i},{j},{k}"), value);
                }
            }
        }
    }

    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&polytope, &triangulation_points)
            .expect("failed Batyrev h11 correction");

    let report = serde_json::json!({
        "n_points": points_raw.len(),
        "n_triangulation_points": triangulation_points.len(),
        "dual_not_facet_points": dual_not_facet,
        "simplices": sorted_simplices(triangulation.simplices()),
        "default_chamber_simplices": sorted_simplices(default_triangulation.simplices()),
        "basis": basis,
        "h11_toric": basis.len(),
        "h11_extra": h11_extra,
        "kappa_basis": kappa_entries,
    });
    std::fs::write(
        PathBuf::from(&out_path),
        serde_json::to_string(&report).expect("serialize report"),
    )
    .unwrap_or_else(|e| panic!("failed to write {out_path}: {e}"));
}
