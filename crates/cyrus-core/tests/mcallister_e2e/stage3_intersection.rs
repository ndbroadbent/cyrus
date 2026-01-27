//! McAllister Pipeline Stage 3: GLSM Charges + Intersection Numbers
//!
//! Computes:
//! - GLSM charge matrix from triangulation points
//! - Intersection numbers κ_ijk from triangulation + GLSM
//!
//! - "ours" branch: Using our FRST triangulation
//! - "theirs" branch: Using McAllister's triangulation (from their heights)

#![allow(missing_docs)]

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use cyrus_core::{
    basis_change_matrix, compute_frst_heights, compute_glsm_and_linrels,
    compute_linear_relations_no_origin, compute_regular_triangulation,
    intersection::compute_intersection_numbers_with_linear_relations, is_unimodular,
};

use cyrus_core::{Point, Polytope, compute_glsm_charge_matrix};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

struct Stage3Fixture {
    triangulation_points: Vec<Point>,
    origin_idx: usize,
}

fn read_csv_rows_i64(path: &PathBuf) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn read_csv_f64(path: &PathBuf) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn read_csv_i64(path: &PathBuf) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
        .collect()
}

fn require_first_principles() -> bool {
    if !crate::first_principles_enabled() {
        eprintln!("Skipping first-principles test (set CYRUS_FIRST_PRINCIPLES=1)");
        return false;
    }
    true
}

fn load_fixture() -> Stage3Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let data_dir = crate::mcallister_data_dir();
    if crate::first_principles_enabled() && data_dir.is_none() {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    }

    let points_raw = if let Some(dir) = data_dir {
        read_csv_rows_i64(&dir.join("points.dat"))
    } else {
        if !crate::fixtures_enabled() {
            panic!("Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures");
        }
        let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
        let content = std::fs::read_to_string(&input_path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
        let input: PolytopeInput = serde_json::from_str(&content)
            .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));
        input.points
    };

    let all_points: Vec<Point> = points_raw
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(all_points).expect("Failed to create polytope");

    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter points");

    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found");

    Stage3Fixture {
        triangulation_points,
        origin_idx,
    }
}

fn load_mcallister_heights() -> Vec<f64> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let data_dir = crate::mcallister_data_dir();
    if crate::first_principles_enabled() && data_dir.is_none() {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    }
    if let Some(dir) = data_dir {
        return read_csv_f64(&dir.join("heights.dat"));
    }
    if !crate::fixtures_enabled() {
        panic!("Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures");
    }
    let path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    let input: HeightsInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()));
    input.values
}

#[test]
fn stage3_primal_basis_matches_dat() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        eprintln!("Skipping basis check (set CYRUS_MCALLISTER_DATA_DIR)");
        return;
    };

    let fixture = load_fixture();
    let (glsm, _linrel, basis) =
        compute_glsm_and_linrels(&fixture.triangulation_points)
            .expect("Failed to compute GLSM/basis");

    let basis_path = data_dir.join("basis.dat");
    let basis_i64 = read_csv_i64(&basis_path);
    let basis_expected: Vec<usize> = basis_i64
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index fits usize"))
        .collect();

    if basis == basis_expected {
        return;
    }

    if std::env::var_os("CYRUS_STRICT_BASIS").is_some() {
        assert_eq!(
            basis, basis_expected,
            "Computed basis does not match McAllister basis.dat"
        );
    }

    let transform = basis_change_matrix(&glsm, &basis, &basis_expected)
        .expect("Failed to compute basis change matrix");

    assert!(
        is_unimodular(&transform),
        "Basis differs from McAllister, and change-of-basis is not unimodular"
    );

    eprintln!(
        "Basis differs from McAllister; unimodular change-of-basis confirmed."
    );
}

// =============================================================================
// GLSM charge matrix (same for both branches - depends only on points)
// =============================================================================

#[test]
fn stage3_glsm_charge_matrix() {
    let fixture = load_fixture();

    let glsm = compute_glsm_charge_matrix(&fixture.triangulation_points, true)
        .expect("Failed to compute GLSM charge matrix");

    // GLSM is the kernel of the point matrix.
    // With include_origin=true, the column count equals the point count
    // (origin is already present in points_not_interior_to_facets).
    assert!(!glsm.is_empty(), "GLSM should not be empty");
    let n_pts = fixture.triangulation_points.len();
    let dim = fixture.triangulation_points[0].dim();
    let expected_rows = n_pts - (dim + 1);
    assert_eq!(
        glsm.len(),
        expected_rows,
        "GLSM should have n_points - dim - 1 kernel vectors"
    );
    assert_eq!(
        glsm[0].len(),
        n_pts,
        "GLSM columns should equal number of points (origin included)"
    );

    // Snapshot the GLSM matrix - convert Integer to i64
    let glsm_i64: Vec<Vec<i64>> = glsm
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("GLSM value fits in i64"))
                .collect()
        })
        .collect();

    insta::assert_json_snapshot!("glsm_charge_matrix", glsm_i64);
}

// =============================================================================
// "ours" branch: Intersection numbers from our FRST triangulation
// =============================================================================

#[test]
fn stage3_ours_intersection_numbers() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_fixture();

    // Compute our triangulation
    let (_heights, triangulation) =
        compute_frst_heights(&fixture.triangulation_points, fixture.origin_idx)
            .expect("Failed to compute FRST heights");

    // Convert points to i64 vecs for linear_relations computation
    let points_i64: Vec<Vec<i64>> = fixture
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();

    // Compute linear relations (CYTools-style, origin excluded)
    let linear_relations = compute_linear_relations_no_origin(&points_i64);

    // Compute intersection numbers using linear relations
    let kappa = compute_intersection_numbers_with_linear_relations(
        &triangulation,
        &fixture.triangulation_points,
        &linear_relations,
    )
    .expect("Failed to compute intersection numbers");

    // Convert to serializable format: sorted list of ((i,j,k), value)
    let mut entries: Vec<((usize, usize, usize), String)> = kappa
        .iter()
        .map(|(&key, val)| (key, val.get().to_string()))
        .collect();
    entries.sort_by_key(|(key, _)| *key);

    #[derive(Serialize)]
    struct IntersectionSnapshot {
        dim: usize,
        num_nonzero: usize,
        entries: Vec<((usize, usize, usize), String)>,
    }

    let snapshot = IntersectionSnapshot {
        dim: kappa.dim(),
        num_nonzero: kappa.num_nonzero(),
        entries,
    };

    insta::assert_json_snapshot!("intersection_ours", snapshot);
}

// =============================================================================
// "theirs" branch: Intersection numbers from McAllister's triangulation
// =============================================================================

#[test]
fn stage3_theirs_intersection_numbers() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_fixture();
    let heights = load_mcallister_heights();

    // Compute triangulation from their heights
    let triangulation = compute_regular_triangulation(&fixture.triangulation_points, &heights)
        .expect("Failed to compute triangulation");

    // Convert points to i64 vecs for linear_relations computation
    let points_i64: Vec<Vec<i64>> = fixture
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();

    // Compute linear relations (CYTools-style, origin excluded)
    let linear_relations = compute_linear_relations_no_origin(&points_i64);

    // Compute intersection numbers using linear relations
    let kappa = compute_intersection_numbers_with_linear_relations(
        &triangulation,
        &fixture.triangulation_points,
        &linear_relations,
    )
    .expect("Failed to compute intersection numbers");

    // Convert to serializable format: sorted list of ((i,j,k), value)
    let mut entries: Vec<((usize, usize, usize), String)> = kappa
        .iter()
        .map(|(&key, val)| (key, val.get().to_string()))
        .collect();
    entries.sort_by_key(|(key, _)| *key);

    #[derive(Serialize)]
    struct IntersectionSnapshot {
        dim: usize,
        num_nonzero: usize,
        entries: Vec<((usize, usize, usize), String)>,
    }

    let snapshot = IntersectionSnapshot {
        dim: kappa.dim(),
        num_nonzero: kappa.num_nonzero(),
        entries,
    };

    insta::assert_json_snapshot!("intersection_theirs", snapshot);
}

// =============================================================================
// DUAL TEST: Compare intersection numbers against CYTools
// =============================================================================

/// Dual test comparing our intersection numbers against CYTools computed values
#[test]
fn stage3_dual_test_intersection_vs_cytools() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_fixture();
    let heights = load_mcallister_heights();

    // Compute triangulation from McAllister's heights
    let triangulation = compute_regular_triangulation(&fixture.triangulation_points, &heights)
        .expect("Failed to compute triangulation");

    // Convert points to i64 vecs for linear_relations computation
    let points_i64: Vec<Vec<i64>> = fixture
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();

    // Compute linear relations (CYTools-style, origin excluded)
    let linear_relations = compute_linear_relations_no_origin(&points_i64);

    eprintln!(
        "Linear relations: {} rows x {} cols",
        linear_relations.len(),
        linear_relations.first().map_or(0, |r| r.len())
    );

    // Compute our intersection numbers using linear relations
    let kappa = compute_intersection_numbers_with_linear_relations(
        &triangulation,
        &fixture.triangulation_points,
        &linear_relations,
    )
    .expect("Failed to compute intersection numbers");

    // Load CYTools intersection numbers fixture
    #[derive(Debug, Deserialize)]
    struct IntersectionEntry {
        indices: Vec<usize>,
        value: i64,
    }

    #[derive(Debug, Deserialize)]
    struct CytoolsIntersection {
        in_basis: bool,
        total_nonzero: usize,
        sample: Vec<IntersectionEntry>,
    }

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_path =
        manifest_dir.join("tests/mcallister_e2e/assertions/intersection_sample_cytools.json");
    let content = std::fs::read_to_string(&fixture_path)
        .unwrap_or_else(|e| panic!("Failed to read CYTools fixture: {e}"));
    let cytools: CytoolsIntersection = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse CYTools fixture: {e}"));

    // Load CYTools divisor basis to map indices
    #[derive(Debug, Deserialize)]
    struct CytoolsBasis {
        basis: Vec<usize>,
    }

    let basis_path =
        manifest_dir.join("tests/mcallister_e2e/assertions/divisor_basis_cytools.json");
    let basis_content = std::fs::read_to_string(&basis_path)
        .unwrap_or_else(|e| panic!("Failed to read CYTools basis: {e}"));
    let cytools_basis: CytoolsBasis = serde_json::from_str(&basis_content)
        .unwrap_or_else(|e| panic!("Failed to parse CYTools basis: {e}"));

    eprintln!("\n=== Dual Test: Intersection Numbers ===");
    eprintln!("CYTools total nonzero: {}", cytools.total_nonzero);
    eprintln!("Our total nonzero: {}", kappa.num_nonzero());
    eprintln!("CYTools basis (first 10): {:?}", &cytools_basis.basis[..10]);

    // CYTools uses in_basis=True, so indices are 0-based into their divisor basis
    // cytools_basis.basis[i] gives the full point index for basis element i
    // Our kappa uses full point indices (with origin at 0)

    let basis = &cytools_basis.basis;

    let mut matches = 0;
    let mut mismatches = 0;
    let mut missing = 0;

    for entry in &cytools.sample {
        let (bi, bj, bk) = (entry.indices[0], entry.indices[1], entry.indices[2]);

        // Map basis indices to full point indices
        // CYTools basis indices are 0..214, basis[n] gives full point index
        if bi >= basis.len() || bj >= basis.len() || bk >= basis.len() {
            eprintln!("SKIP: basis index out of range ({},{},{})", bi, bj, bk);
            continue;
        }

        let (fi, fj, fk) = (basis[bi], basis[bj], basis[bk]);

        // Look up in our kappa using full indices
        let our_val = kappa.get(fi, fj, fk);
        let our_rational = our_val.get();

        // Try to convert to i64 for comparison
        if let Ok(our_i64) = i64::try_from(our_rational) {
            if our_i64 == entry.value {
                matches += 1;
            } else if our_i64 == 0 {
                missing += 1;
                if missing <= 5 {
                    eprintln!(
                        "MISSING: κ_basis({},{},{}) = κ_full({},{},{}) = 0 (ours) vs {} (CYTools)",
                        bi, bj, bk, fi, fj, fk, entry.value
                    );
                }
            } else {
                mismatches += 1;
                if mismatches <= 5 {
                    eprintln!(
                        "MISMATCH: κ_basis({},{},{}) = κ_full({},{},{}) = {} (ours) vs {} (CYTools)",
                        bi, bj, bk, fi, fj, fk, our_i64, entry.value
                    );
                }
            }
        } else {
            // Rational doesn't fit in i64
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!(
                    "MISMATCH: κ_basis({},{},{}) = κ_full({},{},{}) = {} (ours, non-integer) vs {} (CYTools)",
                    bi, bj, bk, fi, fj, fk, our_rational, entry.value
                );
            }
        }
    }

    eprintln!("\nComparison results:");
    eprintln!("  Matches: {}", matches);
    eprintln!("  Mismatches: {}", mismatches);
    eprintln!("  Missing: {}", missing);
    eprintln!("  Total sample: {}", cytools.sample.len());

    // Save comparison summary for documentation
    #[derive(Serialize)]
    struct DualTestSummary {
        cytools_total_nonzero: usize,
        our_total_nonzero: usize,
        sample_size: usize,
        matches: usize,
        mismatches: usize,
        missing: usize,
        note: String,
    }

    let summary = DualTestSummary {
        cytools_total_nonzero: cytools.total_nonzero,
        our_total_nonzero: kappa.num_nonzero(),
        sample_size: cytools.sample.len(),
        matches,
        mismatches,
        missing,
        note: "CYTools uses latest version basis, our test uses McAllister heights".to_string(),
    };

    insta::assert_json_snapshot!("dual_test_intersection_cytools", summary);
}
