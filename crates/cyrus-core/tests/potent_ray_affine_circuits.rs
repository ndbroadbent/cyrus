//! First-principles validation of McAllister potent-ray affine circuit structure.
#![allow(missing_docs, clippy::too_many_lines)]

use std::collections::BTreeMap;
use std::path::Path;

use cyrus_core::{
    AffineToricCircuitDiagnostic, LocalToricCircuitKind, LocalToricCoordinate2D, Point, Polytope,
    RankTwoLocalSupportSignature, compute_local_toric_circuit_gv_series,
    curve_in_rational_row_span, diagnose_affine_toric_circuit, rank_two_local_charge_model,
    rank_two_local_support_signature,
};
use malachite::Integer;

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

fn read_csv_rows_integer(path: &Path) -> Vec<Vec<Integer>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|value| {
                    value
                        .trim()
                        .parse::<Integer>()
                        .unwrap_or_else(|err| panic!("invalid integer {value}: {err:?}"))
                })
                .collect()
        })
        .collect()
}

fn assert_rank_two_local_coordinates_satisfy_relation(diagnostic: &AffineToricCircuitDiagnostic) {
    let coordinates = diagnostic
        .local_coordinates_2d
        .as_ref()
        .expect("rank-two affine support should have local coordinates");
    assert_eq!(
        coordinates.len(),
        diagnostic.relation_points.len(),
        "local coordinate count should match affine support size"
    );

    let coordinates_by_point: BTreeMap<usize, [i64; 2]> = coordinates
        .iter()
        .map(|point| (point.point_index, point.coordinates))
        .collect();
    let mut coordinate_sum = [0i128; 2];
    for point in &diagnostic.relation_points {
        let coordinates = coordinates_by_point
            .get(&point.point_index)
            .expect("every relation point should have local coordinates");
        coordinate_sum[0] += i128::from(point.coefficient) * i128::from(coordinates[0]);
        coordinate_sum[1] += i128::from(point.coefficient) * i128::from(coordinates[1]);
    }
    assert_eq!(
        coordinate_sum,
        [0, 0],
        "affine relation should remain zero in reconstructed local coordinates"
    );
    assert!(
        rank_two_local_support_signature(diagnostic).is_some(),
        "rank-two affine support should have a normalized local signature"
    );
    let signature = rank_two_local_support_signature(diagnostic)
        .expect("rank-two affine support should have a normalized local signature");
    let model = rank_two_local_charge_model(&signature)
        .expect("rank-two affine support should produce a local charge model");
    assert_eq!(
        model.points.len(),
        diagnostic.relation_points.len(),
        "local charge-model point count should match affine support size"
    );
    assert_eq!(
        model.charge_basis.len(),
        diagnostic.relation_points.len() - 3,
        "local charge-model rank should equal support size minus rank-two affine dimension plus one"
    );
    assert!(
        curve_in_rational_row_span(&model.target_relation, &model.charge_basis)
            .expect("local charge-model span check should be exact"),
        "target potent-ray relation should lie in the local charge lattice"
    );
    assert_eq!(
        model.target_relation_in_charge_basis.len(),
        model.charge_basis.len(),
        "target charge coordinates should have one entry per local charge-basis row"
    );
    let mut reconstructed_target = vec![0i64; model.target_relation.len()];
    for (coordinate, charge) in model
        .target_relation_in_charge_basis
        .iter()
        .zip(model.charge_basis.iter())
    {
        for (out, charge_value) in reconstructed_target.iter_mut().zip(charge.iter()) {
            *out += coordinate * charge_value;
        }
    }
    assert_eq!(
        reconstructed_target, model.target_relation,
        "target charge coordinates should reconstruct the potent-ray relation exactly"
    );
    for charge in &model.charge_basis {
        let mut local_coordinate_sum = [0i128; 2];
        let mut coefficient_sum = 0i128;
        for (coefficient, point) in charge.iter().zip(model.points.iter()) {
            coefficient_sum += i128::from(*coefficient);
            local_coordinate_sum[0] += i128::from(*coefficient) * i128::from(point.coordinates[0]);
            local_coordinate_sum[1] += i128::from(*coefficient) * i128::from(point.coordinates[1]);
        }
        assert_eq!(
            coefficient_sum, 0,
            "each local charge-model basis row should be affine"
        );
        assert_eq!(
            local_coordinate_sum,
            [0, 0],
            "each local charge-model basis row should annihilate local coordinates"
        );
    }
}

fn rank_two_signature_coefficient_pattern(signature: &RankTwoLocalSupportSignature) -> Vec<i64> {
    let mut coefficients = signature
        .entries
        .iter()
        .map(|entry| entry.coefficient)
        .collect::<Vec<_>>();
    coefficients.sort_unstable();
    coefficients
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
    let mut rank_two_local_coordinate_count = 0usize;
    let mut affine_rank_counts = BTreeMap::new();
    let mut signature_counts: BTreeMap<RankTwoLocalSupportSignature, usize> = BTreeMap::new();
    let mut coefficient_pattern_counts: BTreeMap<Vec<i64>, usize> = BTreeMap::new();
    let mut local_p2_signature: Option<RankTwoLocalSupportSignature> = None;
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
        if diagnostic.affine_rank == 2 {
            rank_two_local_coordinate_count += 1;
            assert_rank_two_local_coordinates_satisfy_relation(&diagnostic);
            let signature = rank_two_local_support_signature(&diagnostic)
                .expect("rank-two diagnostic should have a local signature");
            *signature_counts.entry(signature.clone()).or_insert(0) += 1;
            *coefficient_pattern_counts
                .entry(rank_two_signature_coefficient_pattern(&signature))
                .or_insert(0) += 1;
            if matches!(
                diagnostic.kind,
                Some(LocalToricCircuitKind::LocalP2Triangle { .. })
            ) {
                match &local_p2_signature {
                    Some(existing) => assert_eq!(
                        existing, &signature,
                        "all recognized local P2 supports should share one normalized signature"
                    ),
                    None => local_p2_signature = Some(signature),
                }
            }
        } else {
            assert!(
                diagnostic.local_coordinates_2d.is_none(),
                "non-rank-two supports should not expose rank-two local coordinates"
            );
        }
    }

    assert_eq!(potent_rays.len(), 411);
    assert_eq!(affine_count, 411);
    assert_eq!(affine_rank_counts, BTreeMap::from([(2, 395), (4, 16)]));
    assert_eq!(rank_two_local_coordinate_count, 395);
    assert_eq!(local_p2_count, 56);
    assert_eq!(
        signature_counts.len(),
        16,
        "rank-two local supports should collapse to the audited signature inventory"
    );
    assert_eq!(
        coefficient_pattern_counts,
        BTreeMap::from([
            (vec![-14, 1, 4, 4, 5], 7),
            (vec![-12, 1, 1, 5, 5], 31),
            (vec![-11, 1, 1, 4, 5], 19),
            (vec![-11, 1, 3, 3, 4], 9),
            (vec![-10, 1, 1, 4, 4], 32),
            (vec![-10, 2, 2, 3, 3], 34),
            (vec![-9, 1, 1, 2, 2, 3], 1),
            (vec![-9, 1, 1, 3, 4], 25),
            (vec![-8, 1, 1, 3, 3], 34),
            (vec![-8, 1, 2, 2, 3], 32),
            (vec![-7, 1, 1, 1, 2, 2], 1),
            (vec![-7, 1, 1, 2, 3], 32),
            (vec![-6, 1, 1, 2, 2], 31),
            (vec![-5, 1, 1, 1, 2], 34),
            (vec![-4, 1, 1, 1, 1], 17),
            (vec![-3, 1, 1, 1], 56),
        ]),
        "rank-two local support families should be inventoried without using GV values"
    );
    assert_eq!(
        signature_counts
            .get(
                local_p2_signature
                    .as_ref()
                    .expect("at least one local P2 signature should be present")
            )
            .copied(),
        Some(56),
        "the normalized local P2 signature should account for every recognized local P2 ray"
    );
    assert!(
        signature_counts.len() > 1,
        "rank-two local support signatures should distinguish non-P2 families"
    );

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
            local_coordinates: vec![
                LocalToricCoordinate2D {
                    point_index: 43,
                    coordinates: [1, 0],
                },
                LocalToricCoordinate2D {
                    point_index: 155,
                    coordinates: [0, 1],
                },
                LocalToricCoordinate2D {
                    point_index: 168,
                    coordinates: [0, 0],
                },
                LocalToricCoordinate2D {
                    point_index: 169,
                    coordinates: [-1, -1],
                },
            ],
        })
    );
    assert_rank_two_local_coordinates_satisfy_relation(&first);
}

#[test]
fn mcallister_potent_rays_have_local_affine_charge_contexts() {
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

    for ray in &potent_rays {
        let diagnostic = diagnose_affine_toric_circuit(ray, &triangulation_points)
            .expect("affine circuit diagnostic should accept McAllister dimensions")
            .expect("saved potent ray should be an affine toric circuit");
        let expected_relation_rank = diagnostic.relation_points.len() - diagnostic.affine_rank - 1;
        assert_eq!(
            diagnostic.local_charge_basis.len(),
            expected_relation_rank,
            "local charge-basis rank should match support affine rank"
        );
        if diagnostic.affine_rank == 2 {
            assert_rank_two_local_coordinates_satisfy_relation(&diagnostic);
        } else {
            assert!(
                diagnostic.local_coordinates_2d.is_none(),
                "non-rank-two supports should not expose rank-two local coordinates"
            );
        }

        let support_relation: Vec<i64> = diagnostic
            .relation_points
            .iter()
            .map(|point| point.coefficient)
            .collect();
        assert!(
            curve_in_rational_row_span(&support_relation, &diagnostic.local_charge_basis)
                .expect("local charge-basis span check should be exact"),
            "saved potent-ray relation must lie in the reconstructed local affine charge lattice"
        );
    }
}

#[test]
fn first_mcallister_local_p2_potent_ray_gvs_are_reconstructed() {
    if !first_principles_enabled() {
        return;
    }
    let Some(data_dir) = mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let potent_rays = read_csv_rows_i64(&data_dir.join("potent_rays.dat"));
    let expected_gv_rows = read_csv_rows_integer(&data_dir.join("potent_rays_gv.dat"));
    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter triangulation points");

    let first_ray = potent_rays
        .first()
        .expect("McAllister potent-ray checkpoint is nonempty");
    let first_expected_gvs = expected_gv_rows
        .first()
        .expect("McAllister potent-ray GV checkpoint is nonempty");
    let diagnostic = diagnose_affine_toric_circuit(first_ray, &triangulation_points)
        .expect("first potent ray diagnostic should accept McAllister dimensions")
        .expect("first potent ray should be an affine circuit");
    let kind = diagnostic
        .kind
        .as_ref()
        .expect("first potent ray should be the local P2 triangle");

    let computed = compute_local_toric_circuit_gv_series(kind, first_expected_gvs.len())
        .expect("local P2 GV reconstruction should succeed")
        .expect("local P2 circuit should have a GV series");

    assert_eq!(
        computed, *first_expected_gvs,
        "first saved potent-ray GV row must be reproduced from the reconstructed local P2 model"
    );
}

#[test]
fn mcallister_rank_two_local_charge_models_are_inventoried() {
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

    let mut inventory = BTreeMap::new();
    for ray in &potent_rays {
        let diagnostic = diagnose_affine_toric_circuit(ray, &triangulation_points)
            .expect("affine circuit diagnostic should accept McAllister dimensions")
            .expect("saved potent ray should be an affine toric circuit");
        if diagnostic.affine_rank != 2 {
            continue;
        }
        let signature = rank_two_local_support_signature(&diagnostic)
            .expect("rank-two diagnostic should have a local signature");
        let coefficient_pattern = rank_two_signature_coefficient_pattern(&signature);
        let model =
            rank_two_local_charge_model(&signature).expect("rank-two support should have a model");
        *inventory
            .entry((
                coefficient_pattern,
                model.charge_basis,
                model.target_relation_in_charge_basis,
            ))
            .or_insert(0usize) += 1;
    }

    assert_eq!(
        inventory,
        BTreeMap::from([
            (
                (
                    vec![-14, 1, 4, 4, 5],
                    vec![vec![1, 1, -1, -1, 0], vec![2, -1, 0, 0, -1]],
                    vec![-4, -5],
                ),
                7,
            ),
            (
                (
                    vec![-12, 1, 1, 5, 5],
                    vec![vec![2, -1, -1, 0, 0], vec![2, 0, 0, -1, -1]],
                    vec![-1, -5],
                ),
                31,
            ),
            (
                (
                    vec![-11, 1, 1, 4, 5],
                    vec![vec![1, -1, -1, 1, 0], vec![3, -1, -1, 0, -1]],
                    vec![4, -5],
                ),
                19,
            ),
            (
                (
                    vec![-11, 1, 3, 3, 4],
                    vec![vec![1, 1, -1, -1, 0], vec![2, -1, 0, 0, -1]],
                    vec![-3, -4],
                ),
                9,
            ),
            (
                (
                    vec![-10, 1, 1, 4, 4],
                    vec![vec![2, -1, -1, 0, 0], vec![2, 0, 0, -1, -1]],
                    vec![-1, -4],
                ),
                32,
            ),
            (
                (
                    vec![-10, 2, 2, 3, 3],
                    vec![vec![2, -1, -1, 0, 0], vec![2, 0, 0, -1, -1]],
                    vec![-2, -3],
                ),
                34,
            ),
            (
                (
                    vec![-9, 1, 1, 2, 2, 3],
                    vec![
                        vec![1, -1, 1, 0, -1, 0],
                        vec![1, 1, -1, 0, 0, -1],
                        vec![2, -1, 0, -1, 0, 0],
                    ],
                    vec![-2, -3, -2],
                ),
                1,
            ),
            (
                (
                    vec![-9, 1, 1, 3, 4],
                    vec![vec![1, -1, -1, 1, 0], vec![3, -1, -1, 0, -1]],
                    vec![3, -4],
                ),
                25,
            ),
            (
                (
                    vec![-8, 1, 1, 3, 3],
                    vec![vec![2, -1, -1, 0, 0], vec![2, 0, 0, -1, -1]],
                    vec![-1, -3],
                ),
                34,
            ),
            (
                (
                    vec![-8, 1, 2, 2, 3],
                    vec![vec![1, 1, -1, -1, 0], vec![2, -1, 0, 0, -1]],
                    vec![-2, -3],
                ),
                32,
            ),
            (
                (
                    vec![-7, 1, 1, 1, 2, 2],
                    vec![
                        vec![1, -1, -1, 1, 0, 0],
                        vec![2, -1, 0, 0, 0, -1],
                        vec![2, 0, -1, 0, -1, 0],
                    ],
                    vec![1, -2, -2],
                ),
                1,
            ),
            (
                (
                    vec![-7, 1, 1, 2, 3],
                    vec![vec![1, -1, -1, 1, 0], vec![3, -1, -1, 0, -1]],
                    vec![2, -3],
                ),
                32,
            ),
            (
                (
                    vec![-6, 1, 1, 2, 2],
                    vec![vec![2, -1, -1, 0, 0], vec![2, 0, 0, -1, -1]],
                    vec![-1, -2],
                ),
                31,
            ),
            (
                (
                    vec![-5, 1, 1, 1, 2],
                    vec![vec![1, -1, 1, -1, 0], vec![2, 0, -1, 0, -1]],
                    vec![-1, -2],
                ),
                34,
            ),
            (
                (
                    vec![-4, 1, 1, 1, 1],
                    vec![vec![2, -1, 0, 0, -1], vec![2, 0, -1, -1, 0]],
                    vec![-1, -1],
                ),
                17,
            ),
            (
                (vec![-3, 1, 1, 1], vec![vec![3, -1, -1, -1]], vec![-1],),
                56,
            ),
        ]),
        "rank-two local charge models should be derived from supports before GV assignment"
    );
}
