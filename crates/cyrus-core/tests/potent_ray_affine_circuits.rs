//! First-principles validation of McAllister potent-ray affine circuit structure.

use std::collections::BTreeMap;
use std::path::Path;

use cyrus_core::{LocalToricCircuitKind, Point, Polytope, diagnose_affine_toric_circuit};

fn first_principles_enabled() -> bool {
    std::env::var_os("CYRUS_FIRST_PRINCIPLES").is_some()
}

fn mcallister_data_dir() -> Option<std::path::PathBuf> {
    std::env::var_os("CYRUS_MCALLISTER_DATA_DIR").map(std::path::PathBuf::from)
}

fn read_csv_rows_i64(path: &Path) -> Vec<Vec<i64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|value| {
                    value
                        .trim()
                        .parse::<i64>()
                        .unwrap_or_else(|err| panic!("invalid integer {value}: {err}"))
                })
                .collect()
        })
        .collect()
}

#[test]
fn mcallister_potent_rays_are_affine_toric_circuits() {
    if !first_principles_enabled() {
        return;
    }
    let Some(data_dir) = mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let potent_rays = read_csv_rows_i64(&data_dir.join("potent_rays.dat"));
    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter triangulation points");

    let mut affine_count = 0usize;
    let mut local_p2_count = 0usize;
    let mut affine_rank_counts = BTreeMap::new();
    for ray in &potent_rays {
        let diagnostic = diagnose_affine_toric_circuit(ray, &triangulation_points)
            .expect("affine circuit diagnostic should accept McAllister dimensions")
            .expect("saved potent ray should be an affine toric circuit");
        affine_count += 1;
        *affine_rank_counts
            .entry(diagnostic.affine_rank)
            .or_insert(0usize) += 1;
        if matches!(
            diagnostic.kind,
            Some(LocalToricCircuitKind::LocalP2Triangle { .. })
        ) {
            local_p2_count += 1;
        }
    }

    assert_eq!(potent_rays.len(), 411);
    assert_eq!(affine_count, 411);
    assert_eq!(affine_rank_counts, BTreeMap::from([(2, 395), (4, 16)]));
    assert_eq!(local_p2_count, 56);

    let first = diagnose_affine_toric_circuit(&potent_rays[0], &triangulation_points)
        .expect("first potent-ray diagnostic should succeed")
        .expect("first potent ray should be affine");
    assert_eq!(first.affine_rank, 2);
    assert_eq!(first.coefficient_counts, BTreeMap::from([(-3, 1), (1, 3)]));
    assert_eq!(
        first.kind,
        Some(LocalToricCircuitKind::LocalP2Triangle {
            interior_point: 168,
            vertex_points: vec![43, 155, 169],
            interior_coefficient: -3,
            vertex_coefficient: 1,
        })
    );
}
