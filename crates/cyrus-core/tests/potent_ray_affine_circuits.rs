//! First-principles validation of McAllister potent-ray affine circuit structure.
#![allow(missing_docs, clippy::too_many_lines)]

use std::collections::BTreeMap;
use std::path::Path;

use cyrus_core::{
    AffineToricCircuitDiagnostic, CkyzLocalSurfaceKind, LocalToricCircuitKind,
    LocalToricCoordinate2D, Point, Polytope, RankTwoLocalChargeModel, RankTwoLocalSupportSignature,
    ckyz_local_surface_cover_weight_coefficients, compute_ckyz_local_gv_invariants_for_degrees,
    compute_local_toric_circuit_gv_series, curve_in_rational_row_span,
    diagnose_affine_toric_circuit, identify_ckyz_local_surface, rank_two_local_charge_model,
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

fn ckyz_kind_index(kind: &CkyzLocalSurfaceKind) -> usize {
    match kind {
        CkyzLocalSurfaceKind::LocalP2 => 0,
        CkyzLocalSurfaceKind::HirzebruchF0 => 1,
        CkyzLocalSurfaceKind::HirzebruchF1 => 2,
        CkyzLocalSurfaceKind::Polygon5 => 3,
    }
}

fn scale_ckyz_degree(direction: &[usize], multiple: usize) -> Vec<usize> {
    direction
        .iter()
        .map(|entry| entry * multiple)
        .collect::<Vec<_>>()
}

fn gcd_abs(mut a: i64, mut b: i64) -> i64 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let rem = a % b;
        a = b;
        b = rem;
    }
    a
}

fn cross(origin: [i64; 2], lhs: [i64; 2], rhs: [i64; 2]) -> i64 {
    (lhs[0] - origin[0]) * (rhs[1] - origin[1]) - (lhs[1] - origin[1]) * (rhs[0] - origin[0])
}

fn convex_hull_2d(points: &[[i64; 2]]) -> Vec<[i64; 2]> {
    let mut sorted = points.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    if sorted.len() <= 1 {
        return sorted;
    }

    let mut lower = Vec::new();
    for point in &sorted {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], *point) <= 0
        {
            lower.pop();
        }
        lower.push(*point);
    }

    let mut upper = Vec::new();
    for point in sorted.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], *point) <= 0
        {
            upper.pop();
        }
        upper.push(*point);
    }

    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

fn point_on_segment(point: [i64; 2], start: [i64; 2], end: [i64; 2]) -> bool {
    cross(start, end, point) == 0
        && point[0] >= start[0].min(end[0])
        && point[0] <= start[0].max(end[0])
        && point[1] >= start[1].min(end[1])
        && point[1] <= start[1].max(end[1])
}

fn point_in_convex_polygon(point: [i64; 2], hull: &[[i64; 2]]) -> bool {
    hull.iter()
        .zip(hull.iter().cycle().skip(1))
        .all(|(&start, &end)| cross(start, end, point) >= 0)
}

fn point_on_polygon_boundary(point: [i64; 2], hull: &[[i64; 2]]) -> bool {
    hull.iter()
        .zip(hull.iter().cycle().skip(1))
        .any(|(&start, &end)| point_on_segment(point, start, end))
}

fn assert_rank_two_signature_is_reflexive_polygon(signature: &RankTwoLocalSupportSignature) {
    let negative_points: Vec<[i64; 2]> = signature
        .entries
        .iter()
        .filter(|entry| entry.coefficient < 0)
        .map(|entry| entry.coordinates)
        .collect();
    assert_eq!(
        negative_points.len(),
        1,
        "rank-two support should have one compact interior point"
    );
    let interior = negative_points[0];
    let shifted: Vec<[i64; 2]> = signature
        .entries
        .iter()
        .map(|entry| {
            [
                entry.coordinates[0] - interior[0],
                entry.coordinates[1] - interior[1],
            ]
        })
        .collect();
    let hull = convex_hull_2d(&shifted);
    assert!(
        hull.len() >= 3,
        "rank-two support should have a two-dimensional polygon hull"
    );

    let min_x = shifted.iter().map(|point| point[0]).min().unwrap();
    let max_x = shifted.iter().map(|point| point[0]).max().unwrap();
    let min_y = shifted.iter().map(|point| point[1]).min().unwrap();
    let max_y = shifted.iter().map(|point| point[1]).max().unwrap();
    let mut interior_lattice_points = Vec::new();
    for x in min_x..=max_x {
        for y in min_y..=max_y {
            let point = [x, y];
            if point_in_convex_polygon(point, &hull) && !point_on_polygon_boundary(point, &hull) {
                interior_lattice_points.push(point);
            }
        }
    }
    assert_eq!(
        interior_lattice_points,
        vec![[0, 0]],
        "rank-two local polygon should have the compact point as its unique interior lattice point"
    );

    for (&start, &end) in hull.iter().zip(hull.iter().cycle().skip(1)) {
        let dx = end[0] - start[0];
        let dy = end[1] - start[1];
        let edge_gcd = gcd_abs(dx, dy);
        assert!(edge_gcd > 0, "polygon edge should be nonzero");
        let primitive_normal = [-dy / edge_gcd, dx / edge_gcd];
        let distance = primitive_normal[0] * start[0] + primitive_normal[1] * start[1];
        assert_eq!(
            distance.abs(),
            1,
            "rank-two local polygon edge should have lattice distance one from the compact point"
        );
    }
}

fn permute_columns(matrix: &[Vec<i64>], permutation: &[usize]) -> Vec<Vec<i64>> {
    matrix
        .iter()
        .map(|row| permutation.iter().map(|&idx| row[idx]).collect())
        .collect()
}

fn permute_vector(vector: &[i64], permutation: &[usize]) -> Vec<i64> {
    permutation.iter().map(|&idx| vector[idx]).collect()
}

fn transform_rows(transform: &[Vec<i64>], matrix: &[Vec<i64>]) -> Vec<Vec<i64>> {
    transform
        .iter()
        .map(|row| {
            (0..matrix[0].len())
                .map(|col| {
                    row.iter()
                        .zip(matrix.iter())
                        .map(|(coefficient, source_row)| coefficient * source_row[col])
                        .sum()
                })
                .collect()
        })
        .collect()
}

fn combine_source_charges(coordinates: &[i64], source_relations: &[Vec<i64>]) -> Vec<i64> {
    let mut out = vec![0; source_relations[0].len()];
    for (coordinate, relation) in coordinates.iter().zip(source_relations.iter()) {
        for (value, relation_value) in out.iter_mut().zip(relation.iter()) {
            *value += coordinate * relation_value;
        }
    }
    out
}

fn rank_two_models_by_coefficient_pattern(
    data_dir: &Path,
) -> BTreeMap<Vec<i64>, RankTwoLocalChargeModel> {
    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let potent_rays = read_csv_rows_i64(&data_dir.join("potent_rays.dat"));
    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter triangulation points");

    let mut models_by_pattern = BTreeMap::new();
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
        models_by_pattern
            .entry(coefficient_pattern)
            .and_modify(|existing| assert_eq!(existing, &model))
            .or_insert(model);
    }
    models_by_pattern
}

#[test]
fn mcallister_rank_four_local_affine_supports_are_inventoried() {
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

    let mut inventory: BTreeMap<(Vec<usize>, Vec<Vec<i64>>), BTreeMap<Vec<i64>, usize>> =
        BTreeMap::new();
    for ray in &potent_rays {
        let diagnostic = diagnose_affine_toric_circuit(ray, &triangulation_points)
            .expect("affine circuit diagnostic should accept McAllister dimensions")
            .expect("saved potent ray should be an affine toric circuit");
        if diagnostic.affine_rank != 4 {
            continue;
        }
        let support_indices = diagnostic
            .relation_points
            .iter()
            .map(|point| point.point_index)
            .collect::<Vec<_>>();
        let relation = diagnostic
            .relation_points
            .iter()
            .map(|point| point.coefficient)
            .collect::<Vec<_>>();
        *inventory
            .entry((support_indices, diagnostic.local_charge_basis))
            .or_default()
            .entry(relation)
            .or_insert(0usize) += 1;
    }

    assert_eq!(
        inventory,
        BTreeMap::from([
            (
                (vec![0, 3, 8, 9, 17, 23], vec![vec![1, 3, -1, -1, -1, -1]],),
                BTreeMap::from([(vec![-1, -3, 1, 1, 1, 1], 1)]),
            ),
            (
                (
                    vec![0, 3, 8, 9, 17, 23, 60],
                    vec![vec![1, 3, -1, -1, -1, -1, 0], vec![0, 1, -1, -1, 0, 0, 1],],
                ),
                BTreeMap::from([
                    (vec![-9, -26, 8, 8, 9, 9, 1], 1),
                    (vec![-8, -23, 7, 7, 8, 8, 1], 1),
                    (vec![-7, -20, 6, 6, 7, 7, 1], 1),
                    (vec![-7, -19, 5, 5, 7, 7, 2], 1),
                    (vec![-7, -18, 4, 4, 7, 7, 3], 1),
                    (vec![-6, -17, 5, 5, 6, 6, 1], 1),
                    (vec![-5, -14, 4, 4, 5, 5, 1], 1),
                    (vec![-5, -13, 3, 3, 5, 5, 2], 1),
                    (vec![-5, -12, 2, 2, 5, 5, 3], 1),
                    (vec![-5, -11, 1, 1, 5, 5, 4], 1),
                    (vec![-4, -11, 3, 3, 4, 4, 1], 1),
                    (vec![-4, -9, 1, 1, 4, 4, 3], 1),
                    (vec![-3, -8, 2, 2, 3, 3, 1], 1),
                    (vec![-3, -7, 1, 1, 3, 3, 2], 1),
                    (vec![-2, -5, 1, 1, 2, 2, 1], 1),
                ]),
            ),
        ]),
        "rank-four potent-ray supports should remain explicit local affine charge contexts"
    );
}

#[test]
fn mcallister_rank_two_models_identify_ckyz_sources() {
    if !first_principles_enabled() {
        return;
    }
    let Some(data_dir) = mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let models_by_pattern = rank_two_models_by_coefficient_pattern(&data_dir);
    let expected = BTreeMap::from([
        (
            vec![-14, 1, 4, 4, 5],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![5, 4]),
        ),
        (
            vec![-12, 1, 1, 5, 5],
            (CkyzLocalSurfaceKind::HirzebruchF0, vec![1, 5]),
        ),
        (
            vec![-11, 1, 1, 4, 5],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![5, 1]),
        ),
        (
            vec![-11, 1, 3, 3, 4],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![4, 3]),
        ),
        (
            vec![-10, 1, 1, 4, 4],
            (CkyzLocalSurfaceKind::HirzebruchF0, vec![1, 4]),
        ),
        (
            vec![-10, 2, 2, 3, 3],
            (CkyzLocalSurfaceKind::HirzebruchF0, vec![2, 3]),
        ),
        (
            vec![-9, 1, 1, 2, 2, 3],
            (CkyzLocalSurfaceKind::Polygon5, vec![4, 3, 2]),
        ),
        (
            vec![-9, 1, 1, 3, 4],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![4, 1]),
        ),
        (
            vec![-8, 1, 1, 3, 3],
            (CkyzLocalSurfaceKind::HirzebruchF0, vec![1, 3]),
        ),
        (
            vec![-8, 1, 2, 2, 3],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![3, 2]),
        ),
        (
            vec![-7, 1, 1, 1, 2, 2],
            (CkyzLocalSurfaceKind::Polygon5, vec![3, 2, 2]),
        ),
        (
            vec![-7, 1, 1, 2, 3],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![3, 1]),
        ),
        (
            vec![-6, 1, 1, 2, 2],
            (CkyzLocalSurfaceKind::HirzebruchF0, vec![1, 2]),
        ),
        (
            vec![-5, 1, 1, 1, 2],
            (CkyzLocalSurfaceKind::HirzebruchF1, vec![2, 1]),
        ),
        (
            vec![-4, 1, 1, 1, 1],
            (CkyzLocalSurfaceKind::HirzebruchF0, vec![1, 1]),
        ),
        (vec![-3, 1, 1, 1], (CkyzLocalSurfaceKind::LocalP2, vec![1])),
    ]);

    assert_eq!(models_by_pattern.len(), expected.len());
    for (coefficient_pattern, (expected_kind, expected_target)) in expected {
        let model = models_by_pattern
            .get(&coefficient_pattern)
            .expect("expected rank-two local charge model should be present");
        let identification = identify_ckyz_local_surface(model)
            .expect("CKYZ identification should run exactly")
            .expect("rank-two McAllister model should match a CKYZ source");
        assert_eq!(
            identification.kind, expected_kind,
            "CKYZ source kind mismatch for {coefficient_pattern:?}"
        );
        assert_eq!(
            identification.source_target_direction, expected_target,
            "CKYZ target direction mismatch for {coefficient_pattern:?}"
        );
    }
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
fn mcallister_rank_two_local_supports_are_reflexive_polygons() {
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

    let mut signature_counts: BTreeMap<RankTwoLocalSupportSignature, usize> = BTreeMap::new();
    for ray in &potent_rays {
        let diagnostic = diagnose_affine_toric_circuit(ray, &triangulation_points)
            .expect("affine circuit diagnostic should accept McAllister dimensions")
            .expect("saved potent ray should be an affine toric circuit");
        if diagnostic.affine_rank != 2 {
            continue;
        }
        let signature = rank_two_local_support_signature(&diagnostic)
            .expect("rank-two diagnostic should have a local signature");
        *signature_counts.entry(signature).or_insert(0) += 1;
    }

    assert_eq!(
        signature_counts.len(),
        16,
        "rank-two supports should match the 16 reflexive-polygon local families"
    );
    for signature in signature_counts.keys() {
        assert_rank_two_signature_is_reflexive_polygon(signature);
    }
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
fn mcallister_rank_two_ckyz_potent_ray_gvs_are_reconstructed() {
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

    let mut gv_cache: BTreeMap<(usize, Vec<usize>, usize), BTreeMap<Vec<usize>, Integer>> =
        BTreeMap::new();
    let mut checked_rows = 0usize;
    for (row_index, (ray, expected_gvs)) in
        potent_rays.iter().zip(expected_gv_rows.iter()).enumerate()
    {
        let diagnostic = diagnose_affine_toric_circuit(ray, &triangulation_points)
            .expect("potent-ray diagnostic should accept McAllister dimensions")
            .expect("saved potent ray should be an affine toric circuit");
        if diagnostic.affine_rank != 2 {
            continue;
        }
        let signature = rank_two_local_support_signature(&diagnostic)
            .expect("rank-two diagnostic should have a local signature");
        let model = rank_two_local_charge_model(&signature)
            .expect("rank-two signature should produce a local charge model");
        let identification = identify_ckyz_local_surface(&model)
            .expect("CKYZ identification should run exactly")
            .expect("rank-two McAllister model should match a CKYZ source");
        let cover_weights = ckyz_local_surface_cover_weight_coefficients(&identification.kind);
        let source_target_direction = identification
            .source_target_direction
            .iter()
            .map(|&entry| {
                usize::try_from(entry).expect("CKYZ source target direction should be nonnegative")
            })
            .collect::<Vec<_>>();
        let multiples_to_check = expected_gvs.len().min(2);
        let target_degrees = (1..=multiples_to_check)
            .map(|multiple| scale_ckyz_degree(&source_target_direction, multiple))
            .collect::<Vec<_>>();
        let cache_key = (
            ckyz_kind_index(&identification.kind),
            source_target_direction.clone(),
            multiples_to_check,
        );
        let source_gvs = gv_cache.entry(cache_key).or_insert_with(|| {
            compute_ckyz_local_gv_invariants_for_degrees(
                &identification.source_relations,
                &identification.local_intersection_terms,
                cover_weights,
                &target_degrees,
            )
            .expect("source-derived CKYZ local GV extraction should succeed")
        });

        for (multiple, expected_gv) in (1..=multiples_to_check).zip(expected_gvs.iter()) {
            let degree = scale_ckyz_degree(&source_target_direction, multiple);
            let actual = source_gvs
                .get(&degree)
                .cloned()
                .unwrap_or_else(|| Integer::from(0));
            assert_eq!(
                &actual, expected_gv,
                "CKYZ source GV mismatch for potent-ray row {row_index}, multiple {multiple}, degree {degree:?}"
            );
        }
        checked_rows += 1;
    }

    assert_eq!(
        checked_rows, 395,
        "all rank-two CKYZ potent-ray rows should now have CKYZ-reconstructed first-two GV checks"
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

#[test]
fn mcallister_five_point_rank_two_models_match_ckyz_hirzebruch_data() {
    if !first_principles_enabled() {
        return;
    }
    let Some(data_dir) = mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let models_by_pattern = rank_two_models_by_coefficient_pattern(&data_dir);

    let ckyz_f0 = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
    let ckyz_f1 = vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]];
    let minus_identity = vec![vec![-1, 0], vec![0, -1]];
    let swap_and_negate = vec![vec![0, -1], vec![-1, 0]];
    let f1_shear = vec![vec![1, -1], vec![-1, 0]];

    let expected = BTreeMap::from([
        (
            vec![-14, 1, 4, 4, 5],
            (&ckyz_f1, vec![0, 4, 2, 1, 3], &swap_and_negate, vec![5, 4]),
        ),
        (
            vec![-12, 1, 1, 5, 5],
            (&ckyz_f0, vec![0, 1, 3, 2, 4], &minus_identity, vec![1, 5]),
        ),
        (
            vec![-11, 1, 1, 4, 5],
            (&ckyz_f1, vec![0, 4, 1, 3, 2], &f1_shear, vec![5, 1]),
        ),
        (
            vec![-11, 1, 3, 3, 4],
            (&ckyz_f1, vec![0, 4, 2, 1, 3], &swap_and_negate, vec![4, 3]),
        ),
        (
            vec![-10, 1, 1, 4, 4],
            (&ckyz_f0, vec![0, 1, 3, 2, 4], &minus_identity, vec![1, 4]),
        ),
        (
            vec![-10, 2, 2, 3, 3],
            (&ckyz_f0, vec![0, 1, 3, 2, 4], &minus_identity, vec![2, 3]),
        ),
        (
            vec![-9, 1, 1, 3, 4],
            (&ckyz_f1, vec![0, 4, 1, 3, 2], &f1_shear, vec![4, 1]),
        ),
        (
            vec![-8, 1, 1, 3, 3],
            (&ckyz_f0, vec![0, 1, 3, 2, 4], &minus_identity, vec![1, 3]),
        ),
        (
            vec![-8, 1, 2, 2, 3],
            (&ckyz_f1, vec![0, 4, 2, 1, 3], &swap_and_negate, vec![3, 2]),
        ),
        (
            vec![-7, 1, 1, 2, 3],
            (&ckyz_f1, vec![0, 4, 1, 3, 2], &f1_shear, vec![3, 1]),
        ),
        (
            vec![-6, 1, 1, 2, 2],
            (&ckyz_f0, vec![0, 1, 3, 2, 4], &minus_identity, vec![1, 2]),
        ),
        (
            vec![-5, 1, 1, 1, 2],
            (&ckyz_f1, vec![0, 4, 1, 2, 3], &swap_and_negate, vec![2, 1]),
        ),
        (
            vec![-4, 1, 1, 1, 1],
            (&ckyz_f0, vec![0, 1, 2, 4, 3], &minus_identity, vec![1, 1]),
        ),
    ]);

    for (coefficient_pattern, (source_relations, permutation, transform, source_target)) in expected
    {
        let model = models_by_pattern
            .get(&coefficient_pattern)
            .expect("expected CKYZ-matched local charge model should be present");
        assert_eq!(
            model.points.len(),
            5,
            "CKYZ Hirzebruch comparison only covers five-point local polygons"
        );
        let permuted_charge_basis = permute_columns(&model.charge_basis, &permutation);
        assert_eq!(
            transform_rows(transform, &permuted_charge_basis),
            *source_relations,
            "local charge basis should match the CKYZ source relations for {coefficient_pattern:?}"
        );
        assert_eq!(
            combine_source_charges(&source_target, source_relations),
            permute_vector(&model.target_relation, &permutation),
            "CKYZ target coordinates should reconstruct the potent-ray relation for {coefficient_pattern:?}"
        );
    }
}

#[test]
fn mcallister_six_point_rank_two_models_match_ckyz_polygon5_data() {
    if !first_principles_enabled() {
        return;
    }
    let Some(data_dir) = mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let models_by_pattern = rank_two_models_by_coefficient_pattern(&data_dir);
    let ckyz_polygon5 = vec![
        vec![-1, 1, -1, 1, 0, 0],
        vec![-1, -1, 1, 0, 0, 1],
        vec![-1, 0, 1, -1, 1, 0],
    ];
    let expected = BTreeMap::from([
        (
            vec![-9, 1, 1, 2, 2, 3],
            (
                vec![0, 1, 2, 4, 3, 5],
                vec![vec![-1, 0, 0], vec![0, -1, 0], vec![1, 0, -1]],
                vec![4, 3, 2],
            ),
        ),
        (
            vec![-7, 1, 1, 1, 2, 2],
            (
                vec![0, 1, 3, 2, 5, 4],
                vec![vec![-1, 0, 0], vec![1, 0, -1], vec![1, -1, 0]],
                vec![3, 2, 2],
            ),
        ),
    ]);

    for (coefficient_pattern, (permutation, transform, source_target)) in expected {
        let model = models_by_pattern
            .get(&coefficient_pattern)
            .expect("expected CKYZ polygon-5 local charge model should be present");
        assert_eq!(
            model.points.len(),
            6,
            "CKYZ polygon-5 comparison only covers six-point local polygons"
        );
        let permuted_charge_basis = permute_columns(&model.charge_basis, &permutation);
        assert_eq!(
            transform_rows(&transform, &permuted_charge_basis),
            ckyz_polygon5,
            "local charge basis should match the CKYZ polygon-5 source relations for {coefficient_pattern:?}"
        );
        assert_eq!(
            combine_source_charges(&source_target, &ckyz_polygon5),
            permute_vector(&model.target_relation, &permutation),
            "CKYZ polygon-5 target coordinates should reconstruct the potent-ray relation for {coefficient_pattern:?}"
        );
    }
}
