//! McAllister Pipeline Stage 4: Flat Direction and e^{K₀}
//!
//! Computes the flat direction p and Kähler potential factor e^{K₀}.
//!
//! **Inputs:**
//! - Dual intersection numbers κ̃_abc (from Stage 3, dim=h21=4)
//! - Flux vectors K, M (h21=4 dimensional)
//!
//! **Outputs:**
//! - N matrix: N_ab = κ̃_abc M^c (4×4 matrix)
//! - Flat direction: p = N⁻¹ K (4-vector)
//! - e^{K₀} = (4/3 × κ̃_abc p^a p^b p^c)^{-1}
//!
//! **Expected values (McAllister Table 6.4):**
//! - e^{K₀} ≈ 0.234393
//!
//! **Branches:**
//! - "ours": Uses our computed divisor basis
//! - "theirs": Uses McAllister's dual basis override [3,4,5,8] - must EXACTLY reproduce e^K₀

#![allow(missing_docs)]

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Round a float to N decimal places to avoid floating point non-determinism in snapshots
fn round_to_decimals(value: f64, decimals: u32) -> f64 {
    let multiplier = 10f64.powi(decimals as i32);
    (value * multiplier).round() / multiplier
}

use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    Point, Polytope, Triangulation, compute_divisor_basis, compute_flat_direction_full,
    compute_frst_heights, compute_glsm_charge_matrix, intersection::compute_intersection_cytools,
    intersection_in_basis,
};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct FluxInput {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct DualBasisOverride {
    indices: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct DualPointsOverride {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct DualSimplicesOverride {
    simplices: Vec<Vec<usize>>,
}

/// Fixture containing dual intersection tensor and flux vectors
struct Stage4Fixture {
    /// Dual intersection numbers in basis (h21=4)
    kappa_dual: cyrus_core::Intersection,
    /// K flux vector (h21=4)
    k_flux: Vec<I64<Finite>>,
    /// M flux vector (h21=4)
    m_flux: Vec<I64<Finite>>,
    /// Which basis was used (for documentation)
    basis_used: Vec<usize>,
}

/// Load fixture using OUR computed divisor basis and triangulation
fn load_stage4_fixture_ours() -> Stage4Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // === Compute dual intersection numbers using our triangulation ===

    // Load primal polytope
    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    // Create primal polytope (exclude origin for vertices)
    let primal_verts: Vec<Point> = input
        .points
        .iter()
        .filter(|p| !p.iter().all(|&x| x == 0))
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let primal = Polytope::from_vertices(primal_verts).expect("Failed to create primal polytope");

    // Compute dual polytope
    let dual = primal
        .compute_dual()
        .expect("Failed to compute dual polytope");

    let dual_points: Vec<Point> = dual.vertices().to_vec();

    // Create dual polytope and get triangulation points
    let dual_polytope =
        Polytope::from_vertices(dual_points).expect("Failed to create dual polytope");

    let triangulation_points = dual_polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter dual points");

    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found in dual triangulation points");

    // Compute FRST triangulation
    let (_heights, triangulation) = compute_frst_heights(&triangulation_points, origin_idx)
        .expect("Failed to compute dual FRST heights");

    // Compute linear_relations from points (used by CYTools algorithm)
    let points_as_vecs: Vec<Vec<i64>> = triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linear_relations_int =
        cyrus_core::integer_math::compute_linear_relations_no_origin(&points_as_vecs);

    // Convert to i64
    let linear_relations_no_origin: Vec<Vec<i64>> = linear_relations_int
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("fits in i64"))
                .collect()
        })
        .collect();

    // Compute full intersection tensor using CYTools algorithm
    let kappa_dual_full = compute_intersection_cytools(
        &triangulation,
        &triangulation_points,
        &linear_relations_no_origin,
    )
    .expect("Failed to compute dual intersection numbers");

    // Compute GLSM for divisor basis computation
    let glsm = compute_glsm_charge_matrix(&triangulation_points, true)
        .expect("Failed to compute dual GLSM");

    // Compute our own basis
    let basis = compute_divisor_basis(&triangulation_points, &glsm)
        .expect("Failed to compute dual divisor basis");

    let kappa_dual = intersection_in_basis(&kappa_dual_full, &basis);

    // === Load flux vectors ===

    let flux_path = manifest_dir.join("tests/mcallister_e2e/inputs/flux.json");
    let flux_content = std::fs::read_to_string(&flux_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", flux_path.display()));
    let flux: FluxInput = serde_json::from_str(&flux_content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", flux_path.display()));

    let k_flux: Vec<I64<Finite>> = flux.k.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let m_flux: Vec<I64<Finite>> = flux.m.iter().map(|&v| I64::<Finite>::new(v)).collect();

    Stage4Fixture {
        kappa_dual,
        k_flux,
        m_flux,
        basis_used: basis,
    }
}

/// Load fixture using McAllister's triangulation, basis [3,4,5,8], and flux vectors
fn load_stage4_fixture_theirs() -> Stage4Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // === Load McAllister's dual triangulation data ===

    // Load dual points (McAllister's point ordering)
    let dual_points_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_points.json");
    let content = std::fs::read_to_string(&dual_points_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", dual_points_path.display()));
    let dual_points_data: DualPointsOverride = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", dual_points_path.display()));

    // McAllister's triangulation uses points 0-8 only
    // Point 0 is the origin, points 1-8 are the actual rays (divisors)
    let triangulation_points: Vec<Point> = dual_points_data
        .points
        .iter()
        .take(9)
        .map(|coords| Point::new(coords.clone()))
        .collect();

    // Load dual simplices (McAllister's triangulation)
    let dual_simplices_path =
        manifest_dir.join("tests/mcallister_e2e/overrides/dual_simplices.json");
    let content = std::fs::read_to_string(&dual_simplices_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", dual_simplices_path.display()));
    let dual_simplices_data: DualSimplicesOverride = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", dual_simplices_path.display()));

    let triangulation = Triangulation::new(dual_simplices_data.simplices);

    // Load linear_relations_no_origin from override file (to match McAllister's exact setup)
    #[derive(Debug, Deserialize)]
    struct LinearRelationsOverride {
        linear_relations_no_origin: Vec<Vec<i64>>,
    }
    let lr_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_linear_relations.json");
    let content = std::fs::read_to_string(&lr_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", lr_path.display()));
    let lr_data: LinearRelationsOverride = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", lr_path.display()));
    let linear_relations_no_origin = lr_data.linear_relations_no_origin;

    // Compute intersection numbers using CYTools algorithm
    let kappa_dual_full = compute_intersection_cytools(
        &triangulation,
        &triangulation_points,
        &linear_relations_no_origin,
    )
    .expect("Failed to compute dual intersection numbers");

    // Load McAllister's basis override
    let dual_basis_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_basis.json");
    let content = std::fs::read_to_string(&dual_basis_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", dual_basis_path.display()));
    let basis_override: DualBasisOverride = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", dual_basis_path.display()));

    let basis = basis_override.indices;
    let kappa_dual = intersection_in_basis(&kappa_dual_full, &basis);

    // === Load flux vectors ===

    let flux_path = manifest_dir.join("tests/mcallister_e2e/inputs/flux.json");
    let flux_content = std::fs::read_to_string(&flux_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", flux_path.display()));
    let flux: FluxInput = serde_json::from_str(&flux_content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", flux_path.display()));

    let k_flux: Vec<I64<Finite>> = flux.k.iter().map(|&v| I64::<Finite>::new(v)).collect();
    let m_flux: Vec<I64<Finite>> = flux.m.iter().map(|&v| I64::<Finite>::new(v)).collect();

    Stage4Fixture {
        kappa_dual,
        k_flux,
        m_flux,
        basis_used: basis,
    }
}

// =============================================================================
// BRANCH "OURS": Using our computed divisor basis
//
// NOTE: The "ours" branch computes N matrix and p vector but does NOT require
// a valid e^K₀. This is because McAllister's flux vectors K, M are tuned for
// their specific basis [3,4,5,8]. With our different basis, κ_abc p^a p^b p^c
// may be negative, making e^K₀ invalid.
//
// The "theirs" branch tests the physics correctness with McAllister's exact setup.
// The "ours" branch tests that our computation pipeline works correctly.
// =============================================================================

#[test]
fn stage4_ours_n_matrix() {
    let fixture = load_stage4_fixture_ours();

    // Compute N matrix directly (not through compute_flat_direction_full)
    let n_matrix =
        cyrus_core::flat_direction::compute_n_matrix(&fixture.kappa_dual, &fixture.m_flux);

    // N matrix should be 4x4
    let n_rows = 4;
    let n_cols = 4;

    // Convert to Vec<Vec<f64>> for snapshot
    let n_matrix_f64: Vec<Vec<f64>> = (0..n_rows)
        .map(|i| {
            (0..n_cols)
                .map(|j| round_to_decimals(n_matrix[(i, j)], 10))
                .collect()
        })
        .collect();

    insta::assert_json_snapshot!("ours_n_matrix", n_matrix_f64);
}

#[test]
fn stage4_ours_flat_direction_p() {
    let fixture = load_stage4_fixture_ours();

    // Compute N matrix and solve for p
    let n_matrix =
        cyrus_core::flat_direction::compute_n_matrix(&fixture.kappa_dual, &fixture.m_flux);
    let p = cyrus_core::flat_direction::solve_linear_system_faer(&n_matrix, &fixture.k_flux)
        .expect("Linear system solve failed");

    // p should be 4-dimensional
    assert_eq!(p.len(), 4, "Flat direction p should have 4 components");

    // Snapshot the flat direction (rounded to avoid floating point non-determinism)
    let p_f64: Vec<f64> = p.iter().map(|x| round_to_decimals(x.get(), 10)).collect();

    insta::assert_json_snapshot!("ours_flat_direction_p", p_f64);
}

#[test]
fn stage4_ours_ek0() {
    let fixture = load_stage4_fixture_ours();

    // Compute N matrix and solve for p
    let n_matrix =
        cyrus_core::flat_direction::compute_n_matrix(&fixture.kappa_dual, &fixture.m_flux);
    let p = cyrus_core::flat_direction::solve_linear_system_faer(&n_matrix, &fixture.k_flux)
        .expect("Linear system solve failed");

    // Compute κ_abc p^a p^b p^c
    let mut kappa_ppp = 0.0f64;
    for (&(i, j, k), val) in fixture.kappa_dual.iter() {
        let (val_f, _): (f64, _) = malachite::num::conversion::traits::RoundingFrom::rounding_from(
            val.get(),
            malachite::rounding_modes::RoundingMode::Nearest,
        );
        let contrib = val_f * p[i].get() * p[j].get() * p[k].get();
        // Account for symmetry multiplicity
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        kappa_ppp += mult * contrib;
    }

    // NOTE: With our computed basis and McAllister's flux vectors,
    // κ_abc p^a p^b p^c may be negative. This is expected because
    // the flux vectors were tuned for their specific basis [3,4,5,8].
    // The "theirs" branch verifies correct physics with their exact setup.
    #[derive(Serialize)]
    struct Ek0Snapshot {
        kappa_ppp: f64,
        kappa_ppp_is_positive: bool,
        ek0: Option<f64>,
        basis_used: Vec<usize>,
        note: &'static str,
    }

    let ek0 = if kappa_ppp > 0.0 {
        Some(round_to_decimals(1.0 / (4.0 / 3.0 * kappa_ppp), 10))
    } else {
        None
    };

    let snapshot = Ek0Snapshot {
        kappa_ppp: round_to_decimals(kappa_ppp, 10),
        kappa_ppp_is_positive: kappa_ppp > 0.0,
        ek0,
        basis_used: fixture.basis_used,
        note: "Computed with OUR basis. May be invalid with McAllister's flux vectors. See theirs branch for McAllister reproduction.",
    };

    insta::assert_json_snapshot!("ours_ek0", snapshot);
}

#[test]
fn stage4_ours_full_output() {
    let fixture = load_stage4_fixture_ours();

    // Verify dimensions match
    assert_eq!(fixture.kappa_dual.dim(), 4, "κ̃ should be h21=4 dimensional");
    assert_eq!(fixture.k_flux.len(), 4, "K should be h21=4 dimensional");
    assert_eq!(fixture.m_flux.len(), 4, "M should be h21=4 dimensional");

    // Compute N matrix and p directly
    let n_matrix =
        cyrus_core::flat_direction::compute_n_matrix(&fixture.kappa_dual, &fixture.m_flux);
    let p = cyrus_core::flat_direction::solve_linear_system_faer(&n_matrix, &fixture.k_flux)
        .expect("Linear system solve failed");

    // Compute κ_abc p^a p^b p^c
    let mut kappa_ppp = 0.0f64;
    for (&(i, j, k), val) in fixture.kappa_dual.iter() {
        let (val_f, _): (f64, _) = malachite::num::conversion::traits::RoundingFrom::rounding_from(
            val.get(),
            malachite::rounding_modes::RoundingMode::Nearest,
        );
        let contrib = val_f * p[i].get() * p[j].get() * p[k].get();
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        kappa_ppp += mult * contrib;
    }

    let ek0 = if kappa_ppp > 0.0 {
        Some(round_to_decimals(1.0 / (4.0 / 3.0 * kappa_ppp), 10))
    } else {
        None
    };

    #[derive(Serialize)]
    struct Stage4Output {
        h21: usize,
        basis_used: Vec<usize>,
        k_flux: Vec<i64>,
        m_flux: Vec<i64>,
        n_matrix: Vec<Vec<f64>>,
        p: Vec<f64>,
        kappa_ppp: f64,
        kappa_ppp_is_positive: bool,
        ek0: Option<f64>,
        note: &'static str,
    }

    let output = Stage4Output {
        h21: 4,
        basis_used: fixture.basis_used.clone(),
        k_flux: fixture.k_flux.iter().map(|x| x.get()).collect(),
        m_flux: fixture.m_flux.iter().map(|x| x.get()).collect(),
        n_matrix: (0..4)
            .map(|i| {
                (0..4)
                    .map(|j| round_to_decimals(n_matrix[(i, j)], 10))
                    .collect()
            })
            .collect(),
        p: p.iter().map(|x| round_to_decimals(x.get(), 10)).collect(),
        kappa_ppp: round_to_decimals(kappa_ppp, 10),
        kappa_ppp_is_positive: kappa_ppp > 0.0,
        ek0,
        note: "Computed with OUR basis. May be invalid with McAllister's flux vectors.",
    };

    insta::assert_json_snapshot!("ours_stage4_full_output", output);
}

// =============================================================================
// BRANCH "THEIRS": Using McAllister's triangulation and basis [3,4,5,8]
//
// BUG UNDER INVESTIGATION:
// κ_abc p^a p^b p^c comes out NEGATIVE when it should be positive.
// This means e^K₀ cannot be computed (would be negative/invalid).
//
// This is a CRITICAL discrepancy that must be understood and fixed.
// See CLAUDE.md: "Discrepancies Are GOLD" - we do NOT ignore this.
// =============================================================================

#[test]
fn stage4_theirs_n_matrix() {
    let fixture = load_stage4_fixture_theirs();

    // === DEBUG: Investigate κ_abc p^a p^b p^c negativity ===
    eprintln!("\n=== THEIRS BRANCH DEBUG ===");
    eprintln!("Basis used: {:?}", fixture.basis_used);
    eprintln!(
        "K flux: {:?}",
        fixture.k_flux.iter().map(|x| x.get()).collect::<Vec<_>>()
    );
    eprintln!(
        "M flux: {:?}",
        fixture.m_flux.iter().map(|x| x.get()).collect::<Vec<_>>()
    );

    eprintln!("\nDual intersection numbers (κ_abc):");
    for (&(i, j, k), val) in fixture.kappa_dual.iter() {
        eprintln!("  κ[{},{},{}] = {}", i, j, k, val.get());
    }

    // Compute N matrix and solve for p
    let n_matrix =
        cyrus_core::flat_direction::compute_n_matrix(&fixture.kappa_dual, &fixture.m_flux);
    eprintln!("\nN matrix (N_ab = κ_abc M^c):");
    for i in 0..4 {
        eprintln!(
            "  row {}: {:?}",
            i,
            (0..4).map(|j| n_matrix[(i, j)]).collect::<Vec<_>>()
        );
    }

    let p = cyrus_core::flat_direction::solve_linear_system_faer(&n_matrix, &fixture.k_flux)
        .expect("Linear system solve failed");
    eprintln!("\nFlat direction p = N^{{-1}} K:");
    eprintln!("  p = {:?}", p.iter().map(|x| x.get()).collect::<Vec<_>>());

    // Compute κ_abc p^a p^b p^c manually to see the value
    let mut kappa_ppp = 0.0f64;
    for (&(i, j, k), val) in fixture.kappa_dual.iter() {
        let (val_f, _): (f64, _) = malachite::num::conversion::traits::RoundingFrom::rounding_from(
            val.get(),
            malachite::rounding_modes::RoundingMode::Nearest,
        );
        let contrib = val_f * p[i].get() * p[j].get() * p[k].get();
        // Account for symmetry multiplicity
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        kappa_ppp += mult * contrib;
        eprintln!(
            "  κ[{},{},{}] * p^{} * p^{} * p^{} = {} * {} = {} (mult {})",
            i,
            j,
            k,
            i,
            j,
            k,
            val_f,
            p[i].get() * p[j].get() * p[k].get(),
            contrib,
            mult
        );
    }
    eprintln!("\nκ_abc p^a p^b p^c = {kappa_ppp}");
    eprintln!("Expected: POSITIVE (for valid e^K₀)");
    eprintln!("e^K₀ would be: {}", 1.0 / (4.0 / 3.0 * kappa_ppp));
    eprintln!("=== END DEBUG ===\n");

    let result = compute_flat_direction_full(&fixture.kappa_dual, &fixture.k_flux, &fixture.m_flux)
        .expect("Failed to compute flat direction - κ_abc p^a p^b p^c is negative!");

    // N matrix should be 4x4
    assert_eq!(result.n_matrix.len(), 4, "N matrix should have 4 rows");
    assert!(
        result.n_matrix.iter().all(|row| row.len() == 4),
        "N matrix should have 4 columns"
    );

    // Snapshot the N matrix (rounded to avoid floating point non-determinism)
    let n_matrix_f64: Vec<Vec<f64>> = result
        .n_matrix
        .iter()
        .map(|row| row.iter().map(|x| round_to_decimals(x.get(), 10)).collect())
        .collect();

    insta::assert_json_snapshot!("theirs_n_matrix", n_matrix_f64);
}

#[test]
fn stage4_theirs_flat_direction_p() {
    let fixture = load_stage4_fixture_theirs();

    let result = compute_flat_direction_full(&fixture.kappa_dual, &fixture.k_flux, &fixture.m_flux)
        .expect("Failed to compute flat direction");

    // p should be 4-dimensional
    assert_eq!(
        result.p.len(),
        4,
        "Flat direction p should have 4 components"
    );

    // Snapshot the flat direction (rounded to avoid floating point non-determinism)
    let p_f64: Vec<f64> = result
        .p
        .iter()
        .map(|x| round_to_decimals(x.get(), 10))
        .collect();

    insta::assert_json_snapshot!("theirs_flat_direction_p", p_f64);
}

#[test]
fn stage4_theirs_ek0() {
    let fixture = load_stage4_fixture_theirs();

    let result = compute_flat_direction_full(&fixture.kappa_dual, &fixture.k_flux, &fixture.m_flux)
        .expect("Failed to compute flat direction");

    let ek0 = result.ek0.get();

    // McAllister Table 6.4: e^{K₀} ≈ 0.234393
    // Using their basis [3,4,5,8], we MUST reproduce this exactly
    let mcallister_ek0 = 0.234393;
    let tolerance = 1e-5; // Reasonable tolerance for floating point

    assert!(
        (ek0 - mcallister_ek0).abs() < tolerance,
        "e^{{K₀}} = {ek0} does not match McAllister's {mcallister_ek0} (tolerance {tolerance})"
    );

    #[derive(Serialize)]
    struct Ek0Snapshot {
        value: f64,
        basis_used: Vec<usize>,
        mcallister_expected: f64,
        note: &'static str,
    }

    let snapshot = Ek0Snapshot {
        value: round_to_decimals(ek0, 10),
        basis_used: fixture.basis_used,
        mcallister_expected: mcallister_ek0,
        note: "Using McAllister's dual basis override - must match exactly",
    };

    insta::assert_json_snapshot!("theirs_ek0", snapshot);
}

#[test]
fn stage4_theirs_full_output() {
    let fixture = load_stage4_fixture_theirs();

    // Verify dimensions match
    assert_eq!(fixture.kappa_dual.dim(), 4, "κ̃ should be h21=4 dimensional");
    assert_eq!(fixture.k_flux.len(), 4, "K should be h21=4 dimensional");
    assert_eq!(fixture.m_flux.len(), 4, "M should be h21=4 dimensional");
    assert_eq!(
        fixture.basis_used,
        vec![3, 4, 5, 8],
        "Must use McAllister's basis"
    );

    let result = compute_flat_direction_full(&fixture.kappa_dual, &fixture.k_flux, &fixture.m_flux)
        .expect("Failed to compute flat direction");

    #[derive(Serialize)]
    struct Stage4Output {
        /// Dimensions
        h21: usize,
        /// Which basis was used (must be [3,4,5,8])
        basis_used: Vec<usize>,
        /// Input flux vectors
        k_flux: Vec<i64>,
        m_flux: Vec<i64>,
        /// N matrix (4x4)
        n_matrix: Vec<Vec<f64>>,
        /// Flat direction p (4-vector)
        p: Vec<f64>,
        /// e^{K₀} - must match McAllister's 0.234393
        ek0: f64,
    }

    let output = Stage4Output {
        h21: 4,
        basis_used: fixture.basis_used.clone(),
        k_flux: fixture.k_flux.iter().map(|x| x.get()).collect(),
        m_flux: fixture.m_flux.iter().map(|x| x.get()).collect(),
        n_matrix: result
            .n_matrix
            .iter()
            .map(|row| row.iter().map(|x| round_to_decimals(x.get(), 10)).collect())
            .collect(),
        p: result
            .p
            .iter()
            .map(|x| round_to_decimals(x.get(), 10))
            .collect(),
        ek0: round_to_decimals(result.ek0.get(), 10),
    };

    insta::assert_json_snapshot!("theirs_stage4_full_output", output);
}
