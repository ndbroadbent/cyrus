//! Triangulation benchmarks.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use serde::Deserialize;
use std::path::PathBuf;

use cyrus_core::{Point, Polytope, compute_frst_heights, compute_regular_triangulation};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

struct Fixture {
    triangulation_points: Vec<Point>,
    origin_idx: usize,
    heights: Vec<f64>,
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

    let polytope = Polytope::from_vertices(all_points)
        .expect("Failed to create polytope");

    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter points");

    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found");

    let heights_path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
    let heights_content = std::fs::read_to_string(&heights_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", heights_path.display()));
    let heights_input: HeightsInput = serde_json::from_str(&heights_content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", heights_path.display()));

    Fixture {
        triangulation_points,
        origin_idx,
        heights: heights_input.values,
    }
}

fn bench_regular_triangulation(c: &mut Criterion) {
    let fixture = load_fixture();

    c.bench_function("triangulation/regular/219pts", |b| {
        b.iter(|| {
            compute_regular_triangulation(
                black_box(&fixture.triangulation_points),
                black_box(&fixture.heights),
            ).expect("Triangulation failed")
        });
    });
}

fn bench_frst_heights(c: &mut Criterion) {
    let fixture = load_fixture();

    let mut group = c.benchmark_group("triangulation");
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(30));

    group.bench_function("frst_heights/219pts", |b| {
        b.iter(|| {
            compute_frst_heights(
                black_box(&fixture.triangulation_points),
                black_box(fixture.origin_idx),
            ).expect("FRST failed")
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_regular_triangulation,
    bench_frst_heights,
);
criterion_main!(benches);
