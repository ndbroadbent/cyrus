//! Polytope operation benchmarks.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use serde::Deserialize;
use std::path::PathBuf;

use cyrus_core::{Point, Polytope};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

struct Fixture {
    all_points: Vec<Point>,
    polytope: Polytope,
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

    let polytope = Polytope::from_vertices(all_points.clone())
        .expect("Failed to create polytope");

    Fixture { all_points, polytope }
}

fn bench_from_vertices(c: &mut Criterion) {
    let fixture = load_fixture();

    c.bench_function("polytope/from_vertices/294pts", |b| {
        b.iter(|| {
            Polytope::from_vertices(black_box(fixture.all_points.clone()))
                .expect("Failed to create polytope")
        });
    });
}

fn bench_compute_dual(c: &mut Criterion) {
    let fixture = load_fixture();

    c.bench_function("polytope/compute_dual/294pts", |b| {
        b.iter(|| {
            black_box(&fixture.polytope)
                .compute_dual()
                .expect("Failed to compute dual")
        });
    });
}

fn bench_points_not_interior_to_facets(c: &mut Criterion) {
    let fixture = load_fixture();

    c.bench_function("polytope/points_not_interior_to_facets/294pts", |b| {
        b.iter(|| {
            black_box(&fixture.polytope)
                .points_not_interior_to_facets()
                .expect("Failed to filter points")
        });
    });
}

criterion_group!(
    benches,
    bench_from_vertices,
    bench_compute_dual,
    bench_points_not_interior_to_facets,
);
criterion_main!(benches);
