//! Dump triangulation points (points not interior to facets) with indices.

use cyrus_core::{Point, Polytope};
use std::path::PathBuf;

fn read_points(path: &PathBuf) -> Vec<Point> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let coords = line
                .split(',')
                .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
                .collect::<Vec<_>>();
            Point::new(coords)
        })
        .collect()
}

fn main() {
    let data_dir = std::env::args()
        .nth(1)
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("CYRUS_MCALLISTER_DATA_DIR").map(PathBuf::from))
        .unwrap_or_else(|| {
            eprintln!(
                "Usage: dump_triangulation_points <data_dir> (or set CYRUS_MCALLISTER_DATA_DIR)"
            );
            std::process::exit(2);
        });

    let points = read_points(&data_dir.join("points.dat"));
    let polytope = Polytope::from_vertices(points).expect("Failed to build polytope");
    let tri_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to compute triangulation points");

    for (i, p) in tri_points.iter().enumerate() {
        let coords = p
            .coords()
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>()
            .join(",");
        println!("{i},{coords}");
    }
}
