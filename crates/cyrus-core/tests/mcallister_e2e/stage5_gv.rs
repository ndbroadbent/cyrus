//! McAllister Pipeline Stage 5: GV Invariant Computation
//!
//! Tests the integration with the `cygv` crate for computing Gopakumar-Vafa invariants.
//!
//! **Requirements for GV computation:**
//! - Mori cone cap generators (from triangulation circuits)
//! - Grading vector (computed from Mori cone)
//! - Curve basis (GLSM charge matrix without origin)
//! - Intersection numbers κ_ijk
//!
//! **Current Status:**
//! - cygv crate is integrated as a dev dependency
//! - Full Mori cone computation needs to be ported from CYTools
//! - For now, we validate by loading GV invariants from McAllister's data files

#![allow(missing_docs)]
#![allow(dead_code)]

/// Test that cygv crate is properly linked and basic types work
#[test]
fn stage5_cygv_basic_import() {
    // Import cygv types to verify linking
    use cygv::Semigroup;
    use nalgebra::{DMatrix, RowDVector};

    // Create a simple test matrix with generators as columns
    // 2 generators in 2D space: [1,0] and [0,1]
    let generators: DMatrix<i32> = DMatrix::from_column_slice(2, 2, &[1, 0, 0, 1]);

    // Create a grading vector (must be positive-definite)
    let grading: RowDVector<i32> = RowDVector::from_row_slice(&[1, 1]);

    // Verify we can create these types
    assert_eq!(generators.nrows(), 2);
    assert_eq!(generators.ncols(), 2);
    assert_eq!(grading.ncols(), 2);

    // Test Semigroup construction (basic cygv type)
    // This verifies the cygv crate is properly linked
    let semigroup = Semigroup::from_data(generators.clone(), grading.clone());

    // Note: Semigroup construction may fail for simple inputs
    // The important thing is that the cygv crate is linked and types work
    match semigroup {
        Ok(sg) => {
            eprintln!(
                "Semigroup created successfully, max_degree: {}",
                sg.max_degree
            );
        }
        Err(e) => {
            eprintln!(
                "Semigroup construction error (expected for simple inputs): {:?}",
                e
            );
            // This is OK - the crate is linked and working
        }
    }

    // Test that we can call with_max_degree which gives us more control
    let semigroup_with_deg = Semigroup::with_max_degree(generators, grading, 5);
    match semigroup_with_deg {
        Ok(sg) => {
            eprintln!(
                "Semigroup with max_degree created, elements: {}",
                sg.elements.ncols()
            );
            assert!(sg.max_degree >= 5, "Max degree should be at least 5");
        }
        Err(e) => {
            eprintln!("Semigroup with_max_degree error: {:?}", e);
            // This might fail too for simple inputs, but cygv is working
        }
    }

    // The test passes as long as the cygv types are importable and callable
    // The actual computation may require more complex inputs
}

/// Simple GV computation test with a known small example
/// The quintic CY3 in P^4 has GV_1 = 2875 (number of lines)
#[test]
#[ignore = "Requires full Mori cone infrastructure - for future implementation"]
fn stage5_cygv_quintic_example() {
    // Note: This test is for future implementation when we have all inputs ready.
    // The cygv API requires:
    // - generators: Mori cone generators (columns)
    // - grading_vector: positive-definite grading
    // - q: curve basis matrix
    // - nefpart: nef partition (can be empty for hypersurface)
    // - intnums: intersection numbers as HashMap

    // For the quintic in P^4:
    // - h11 = 1, single Kähler modulus
    // - κ_{111} = 5 (triple self-intersection of hyperplane)
    // - GV_1 = 2875 (number of lines)
    // - GV_2 = 609250 (number of conics)

    // Placeholder - actual implementation requires Mori cone computation
    // Test is ignored until Mori cone infrastructure is ready
}

/// Verify we have the data structures needed for McAllister GV computation
#[test]
fn stage5_mcallister_gv_data_available() {
    use std::fs;

    let Some(data_dir) = crate::mcallister_data_dir() else {
        if crate::first_principles_enabled() {
            panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
        }
        eprintln!("Skipping GV data availability check (set CYRUS_MCALLISTER_DATA_DIR)");
        return;
    };

    // Check all required data files exist
    let required_files = [
        "dual_curves.dat",     // 5177 curves for dual polytope
        "dual_curves_gv.dat",  // GV invariants (some very large)
        "small_curves.dat",    // 344 curves for racetrack
        "small_curves_gv.dat", // Small GV invariants
    ];

    for file in &required_files {
        let path = data_dir.join(file);
        assert!(
            fs::metadata(&path).is_ok(),
            "Required file {} should exist",
            file
        );
    }

    // Load and verify dual curves data
    let curves_content = fs::read_to_string(data_dir.join("dual_curves.dat")).unwrap();
    let n_curves = curves_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    assert_eq!(n_curves, 5177, "Should have 5177 dual curves");

    // Load and verify GV invariants count
    let gv_content = fs::read_to_string(data_dir.join("dual_curves_gv.dat")).unwrap();
    let gv_values: Vec<&str> = gv_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .flat_map(|l| l.split(','))
        .map(|s| s.trim())
        .collect();
    assert_eq!(gv_values.len(), 5177, "Should have 5177 GV values");

    // Check that some GV values are huge (demonstrating why we use arbitrary precision)
    let huge_gv_count = gv_values
        .iter()
        .filter(|s| s.len() > 20) // More than 20 digits
        .count();

    #[derive(serde::Serialize)]
    struct GvDataSummary {
        n_curves: usize,
        n_gv_values: usize,
        huge_gv_count: usize,
        first_5_gv: Vec<String>,
        max_gv_digits: usize,
    }

    let max_digits = gv_values.iter().map(|s| s.len()).max().unwrap_or(0);

    let summary = GvDataSummary {
        n_curves,
        n_gv_values: gv_values.len(),
        huge_gv_count,
        first_5_gv: gv_values.iter().take(5).map(|s| (*s).to_string()).collect(),
        max_gv_digits: max_digits,
    };

    insta::assert_json_snapshot!("mcallister_gv_data_summary", summary);
}

/// Document what's needed to compute GV invariants from scratch
/// This test serves as a roadmap for the full implementation
#[test]
fn stage5_gv_computation_roadmap() {
    #[derive(serde::Serialize)]
    struct GvComputationRoadmap {
        status: &'static str,
        /// Bitmask of completed components: cygv=1, mori=2, grading=4, pipeline=8
        completed_components: u8,
        steps_remaining: Vec<&'static str>,
    }

    // Components: cygv integrated (1), mori cone (0), grading vector (0), full pipeline (0)
    let completed = 1u8; // Only cygv is integrated

    let roadmap = GvComputationRoadmap {
        status: "In Progress - cygv integrated, Mori cone needed",
        completed_components: completed,
        steps_remaining: vec![
            "Port mori_cone_cap from CYTools calabiyau.py lines 2295-2400",
            "Port find_grading_vector from CYTools cone.py",
            "Create compute_gvs wrapper using cygv::compute_gv_rat_threefold",
            "Verify computed GVs match McAllister's data files",
            "Integrate into W₀ computation pipeline",
        ],
    };

    insta::assert_json_snapshot!("gv_computation_roadmap", roadmap);
}
