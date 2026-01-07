#![allow(missing_docs)]
//! GLSM and intersection number benchmarks.
//!
//! Target: <1 second for intersection numbers on 4-214-647.

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use serde::Deserialize;
use std::path::PathBuf;

use cyrus_core::{
    Point, Polytope, compute_glsm_charge_matrix, compute_intersection_numbers,
    compute_regular_triangulation,
};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

struct Fixture {
    triangulation: cyrus_core::Triangulation,
    triangulation_points: Vec<Point>,
    glsm: Vec<Vec<malachite::Integer>>,
}

fn load_fixture() -> Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    let all_points: Vec<Point> = input
        .points
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(all_points).expect("Failed to create polytope");

    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter points");

    let heights_path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
    let heights_content = std::fs::read_to_string(&heights_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", heights_path.display()));
    let heights_input: HeightsInput = serde_json::from_str(&heights_content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", heights_path.display()));

    let triangulation = compute_regular_triangulation(&triangulation_points, &heights_input.values)
        .expect("Failed to compute triangulation");

    let glsm =
        compute_glsm_charge_matrix(&triangulation_points, true).expect("Failed to compute GLSM");

    Fixture {
        triangulation,
        triangulation_points,
        glsm,
    }
}

fn bench_glsm_charge_matrix(c: &mut Criterion) {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path).unwrap();
    let input: PolytopeInput = serde_json::from_str(&content).unwrap();

    let all_points: Vec<Point> = input
        .points
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(all_points).unwrap();
    let triangulation_points = polytope.points_not_interior_to_facets().unwrap();

    c.bench_function("glsm/charge_matrix/219pts", |b| {
        b.iter(|| {
            compute_glsm_charge_matrix(black_box(&triangulation_points), true)
                .expect("GLSM computation failed")
        });
    });
}

fn bench_intersection_numbers(c: &mut Criterion) {
    let fixture = load_fixture();

    let n_simplices = fixture.triangulation.simplices().len();

    let mut group = c.benchmark_group("intersection");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(30));

    group.bench_function(format!("numbers/219pts/{n_simplices}simplices"), |b| {
        b.iter(|| {
            compute_intersection_numbers(
                black_box(&fixture.triangulation),
                black_box(&fixture.triangulation_points),
                black_box(&fixture.glsm),
            )
            .expect("Intersection computation failed")
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_glsm_charge_matrix,
    bench_intersection_numbers,
);
criterion_main!(benches);
