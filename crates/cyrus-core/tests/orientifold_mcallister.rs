//! Validation: the automated orientifold layer must reproduce the declared
//! orientifold data of all five published McAllister examples - the
//! involution parity, the O7/so(8) assignment pattern (c_i = 6 vs 1), the
//! O7 disjointness, the rigidity census, and the published D3 tadpoles
//! (110 / 44 / 60 / 60 / 30).
//!
//! Requires CYRUS_MCALLISTER_PAPER_DATA pointing at the paper_data root
//! (skipped otherwise).

use cyrus_core::orientifold::enumerate_involutions;
use cyrus_core::{
    Point, Polytope, compute_glsm_and_linrels, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_regular_triangulation,
};

struct Example {
    dir: &'static str,
    expected_tadpole: f64,
    /// Published count of rigid prime toric divisors (None = not stated).
    expected_rigid: Option<usize>,
}

const EXAMPLES: &[Example] = &[
    Example {
        dir: "4-214-647",
        expected_tadpole: 110.0,
        expected_rigid: Some(216),
    },
    Example {
        dir: "5-81-3213",
        expected_tadpole: 44.0,
        expected_rigid: Some(84),
    },
    Example {
        dir: "5-113-4627-main",
        expected_tadpole: 60.0,
        expected_rigid: None,
    },
    Example {
        dir: "5-113-4627-alternative",
        expected_tadpole: 60.0,
        expected_rigid: None,
    },
    Example {
        dir: "7-51-13590",
        expected_tadpole: 30.0,
        expected_rigid: Some(53),
    },
];

fn read_rows(path: &std::path::Path) -> Vec<Vec<i64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(',')
                .map(|c| c.trim().parse().expect("integer"))
                .collect()
        })
        .collect()
}

fn read_floats(path: &std::path::Path) -> Vec<f64> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .split([',', '\n'])
        .filter(|c| !c.trim().is_empty())
        .map(|c| c.trim().parse().expect("float"))
        .collect()
}

#[test]
fn orientifold_layer_reproduces_all_five_published_examples() {
    let Ok(root) = std::env::var("CYRUS_MCALLISTER_PAPER_DATA") else {
        eprintln!("CYRUS_MCALLISTER_PAPER_DATA not set; skipping");
        return;
    };
    let root = std::path::PathBuf::from(root);

    for example in EXAMPLES {
        let dir = root.join(example.dir);
        assert!(
            dir.join("points.dat").exists(),
            "{}: points.dat missing",
            example.dir
        );
        let points_raw = read_rows(&dir.join("points.dat"));
        let heights = read_floats(&dir.join("heights.dat"));
        let kklt_basis: Vec<usize> = read_rows(&dir.join("kklt_basis.dat"))[0]
            .iter()
            .map(|&x| usize::try_from(x).expect("index"))
            .collect();
        let c_i = &read_rows(&dir.join("target_volumes.dat"))[0];

        let point_objs: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
        let polytope = Polytope::from_vertices(point_objs).expect("polytope");
        let tri_points = polytope
            .points_not_interior_to_facets()
            .expect("triangulation points");
        let triangulation =
            compute_regular_triangulation(&tri_points, &heights).expect("declared FRST");
        let (_, _, basis) = compute_glsm_and_linrels(&tri_points).expect("glsm");
        let points_i64: Vec<Vec<i64>> = tri_points.iter().map(|p| p.coords().to_vec()).collect();
        let linrels = compute_linear_relations_no_origin(&points_i64);
        let linrels_i64: Vec<Vec<i64>> = linrels
            .iter()
            .map(|r| r.iter().map(|x| i64::try_from(x).expect("i64")).collect())
            .collect();
        let kappa_full = compute_intersection_cytools(&triangulation, &tri_points, &linrels_i64)
            .expect("intersection numbers");

        let h11_extra = cyrus_core::divisor::batyrev_h11_extra_classes(&polytope, &tri_points)
            .expect("batyrev");
        let h11 = basis.len() + h11_extra;
        let h21: usize = example.dir.split('-').next().unwrap().parse().unwrap();

        let models = enumerate_involutions(&polytope, &tri_points, &kappa_full, &basis, h11, h21)
            .expect("enumeration");
        assert_eq!(models.len(), 15, "{}: 15 involution classes", example.dir);

        // The declared involution: c_i = 6 divisors fix sigma.
        let declared_o7: Vec<usize> = kklt_basis
            .iter()
            .zip(c_i.iter())
            .filter(|&(_, &c)| c == 6)
            .map(|(&idx, _)| idx)
            .collect();
        assert!(!declared_o7.is_empty(), "{}: no declared O7s", example.dir);
        let sigma: Vec<i64> = tri_points[declared_o7[0]]
            .coords()
            .iter()
            .map(|c| c.rem_euclid(2))
            .collect();

        let model = models
            .iter()
            .find(|m| m.sigma == sigma)
            .unwrap_or_else(|| panic!("{}: declared sigma {sigma:?} not enumerated", example.dir));

        // Every declared c_i must agree with the automated assignment:
        // so(8) on O7s, one Euclidean D3-brane per irreducible component.
        for (&idx, &declared_c) in kklt_basis.iter().zip(c_i.iter()) {
            let automated_c = if model.gamma[idx] == 1 {
                6
            } else {
                i64::try_from(model.components[idx]).expect("components")
            };
            assert_eq!(
                automated_c, declared_c,
                "{}: c_i mismatch at divisor {idx}",
                example.dir
            );
        }

        // Published tadpole.
        assert!(
            (model.q_d3 - example.expected_tadpole).abs() < 1e-9,
            "{}: tadpole {} vs published {}",
            example.dir,
            model.q_d3,
            example.expected_tadpole
        );

        // O7s must be mutually disjoint for the so(8) story.
        assert!(
            model.o7_pairwise_disjoint,
            "{}: O7 divisors intersect",
            example.dir
        );

        // Rigidity census vs the paper, when stated. The paper counts PURE
        // rigid irreducible components (purity excludes a divisor when
        // h^{2,1}(D-hat) != 0, which this module does not compute), so the
        // published number must be bracketed by our rigid-point count from
        // below and our rigid-component count from above. For the favorable
        // examples the three counts coincide exactly.
        let rigid_points = model
            .divisor_classes
            .iter()
            .filter(|(_, c)| c.is_rigid())
            .count();
        let rigid_components: usize = model
            .divisor_classes
            .iter()
            .filter(|(_, c)| c.is_rigid())
            .map(|(idx, _)| model.components[*idx])
            .sum();
        if let Some(expected) = example.expected_rigid {
            assert!(
                rigid_points <= expected && expected <= rigid_components,
                "{}: paper's {} pure rigid divisors not bracketed by [{} points, {} components]",
                example.dir,
                expected,
                rigid_points,
                rigid_components
            );
        }
        let rigid_count = rigid_points;
        // Every divisor McAllister stabilized must be rigid in our census.
        for &idx in &kklt_basis {
            let rigid = model
                .divisor_classes
                .iter()
                .find(|(i, _)| *i == idx)
                .is_some_and(|(_, c)| c.is_rigid());
            assert!(
                rigid,
                "{}: declared KKLT divisor {idx} classified non-rigid",
                example.dir
            );
        }
        eprintln!(
            "{}: sigma={sigma:?} O7s={} rigid={}/{} tadpole={} OK",
            example.dir,
            model.o7_points.len(),
            rigid_count,
            model.divisor_classes.len(),
            model.q_d3
        );
    }
}
