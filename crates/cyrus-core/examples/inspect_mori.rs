#![allow(missing_docs)]

use std::path::PathBuf;

use cyrus_core::{
    compute_frst_heights, compute_glsm_and_linrels, compute_grading_vector,
    compute_mori_cone_cap_rays, Cone, Point, Polytope,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

fn load_polytope() -> Polytope {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    let points: Vec<Point> = input.points.into_iter().map(Point::new).collect();
    Polytope::from_vertices(points).expect("Failed to create polytope")
}

fn main() {
    let polytope = load_polytope();
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter triangulation points");
    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found");

    let (_heights, triangulation) =
        compute_frst_heights(&triangulation_points, origin_idx)
            .expect("Failed to compute FRST triangulation");

    let (_glsm, _linrels, basis) =
        compute_glsm_and_linrels(&triangulation_points).expect("Failed GLSM/linrels");

    let rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        true,
        false,
        Some(&basis),
    )
    .expect("Failed mori cone cap rays");

    let mut min_ray = i64::MAX;
    let mut max_ray = i64::MIN;
    let mut neg_entries = 0usize;
    let dim = rays[0].len();
    let mut unit_hits = vec![false; dim];
    for r in &rays {
        for &v in r {
            min_ray = min_ray.min(v);
            max_ray = max_ray.max(v);
            if v < 0 {
                neg_entries += 1;
            }
        }
        let mut idx = None;
        let mut ok = true;
        for (i, &v) in r.iter().enumerate() {
            if v == 1 {
                if idx.is_some() {
                    ok = false;
                    break;
                }
                idx = Some(i);
            } else if v != 0 {
                ok = false;
                break;
            }
        }
        if ok {
            if let Some(i) = idx {
                unit_hits[i] = true;
            }
        }
    }
    let unit_count = unit_hits.iter().filter(|&&v| v).count();
    println!(
        "mori rays: count={}, dim={}, min_entry={}, max_entry={}, neg_entries={}, unit_rays={}",
        rays.len(),
        dim,
        min_ray,
        max_ray,
        neg_entries,
        unit_count
    );
    let rays_i128: Vec<Vec<i128>> = rays
        .iter()
        .map(|r| r.iter().map(|&x| x as i128).collect())
        .collect();
    let mut cone = Cone::from_rays(rays_i128);
    let h = cone.hyperplanes().to_vec();
    println!(
        "hyperplanes: count={}, dim={}",
        h.len(),
        h.first().map(Vec::len).unwrap_or(0)
    );
    let grading = compute_grading_vector(&rays).expect("No grading vector found");
    let mut min_g = i64::MAX;
    let mut max_g = i64::MIN;
    let mut neg_g = 0usize;
    let mut g1 = 0usize;
    let mut g2 = 0usize;
    let mut g3 = 0usize;
    let mut g4 = 0usize;
    for &v in &grading {
        min_g = min_g.min(v);
        max_g = max_g.max(v);
        if v < 0 {
            neg_g += 1;
        }
        match v {
            1 => g1 += 1,
            2 => g2 += 1,
            3 => g3 += 1,
            4 => g4 += 1,
            _ => {}
        }
    }
    println!(
        "grading: len={}, min={}, max={}, neg_count={}, count(1)={}, count(2)={}, count(3)={}, count(4)={}",
        grading.len(),
        min_g,
        max_g,
        neg_g
        ,g1,g2,g3,g4
    );
}
