//! McAllister Pipeline Stages 10-11: Volume and Cosmological Constant
//!
//! Tests the final stages of the pipeline:
//! - Stage 10: Compute V_string from intersection numbers + Kähler params
//! - Stage 11: Compute V₀ = -3 e^K₀ (g_s⁷ / (4 V_string)²) W₀²
//!
//! **Expected values (McAllister arXiv:2107.09064):**
//! - V_string ≈ 4711.83
//! - V₀ ≈ -5.5 × 10⁻²⁰³ Mpl⁴

#![allow(missing_docs)]
#![allow(clippy::cast_possible_truncation)] // Exponent values from log10 are small

use serde::Deserialize;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Round a float to N decimal places
fn round_to_decimals(value: f64, decimals: u32) -> f64 {
    let multiplier = 10f64.powi(decimals as i32);
    (value * multiplier).round() / multiplier
}

fn require_data_dir() -> Option<PathBuf> {
    let Some(dir) = crate::mcallister_data_dir() else {
        assert!(
            !crate::first_principles_enabled(),
            "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
        );
        eprintln!("Skipping volume checks (set CYRUS_MCALLISTER_DATA_DIR)");
        return None;
    };
    Some(dir)
}

fn require_first_principles() -> bool {
    if !crate::first_principles_enabled() {
        eprintln!("Skipping first-principles test (set CYRUS_FIRST_PRINCIPLES=1)");
        return false;
    }
    true
}

fn require_runner_heavy() -> bool {
    if std::env::var_os("CYRUS_MCALLISTER_RUNNER_HEAVY").is_none() {
        eprintln!(
            "Skipping full first-principles runner test (set CYRUS_MCALLISTER_RUNNER_HEAVY=1)"
        );
        return false;
    }
    true
}

fn read_csv_i64(path: &PathBuf) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
        .collect()
}

fn read_scalar_f64(path: &PathBuf) -> f64 {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()))
        .trim()
        .parse::<f64>()
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn extract_result_value(stderr: &str, name: &str) -> f64 {
    let prefix = format!("[RESULT] {name} = ");
    stderr
        .lines()
        .find_map(|line| line.strip_prefix(&prefix))
        .unwrap_or_else(|| panic!("missing result line for {name} in runner stderr:\n{stderr}"))
        .parse::<f64>()
        .unwrap_or_else(|e| panic!("failed to parse runner result {name}: {e}"))
}

#[derive(Debug, Deserialize)]
struct RacetrackAssertion {
    g_s: f64,
    #[serde(rename = "W_0")]
    w_0: f64,
    cy_vol: f64,
    #[allow(dead_code)]
    n_curves: usize,
}

/// Load McAllister's ground truth values
fn load_racetrack_assertion() -> RacetrackAssertion {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let path = manifest_dir.join("tests/mcallister_e2e/assertions/racetrack.json");
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Load McAllister's Kähler parameters (t^i values)
fn load_kahler_params(data_dir: &Path) -> (Vec<f64>, Vec<f64>) {
    use std::fs;

    // Load uncorrected Kähler parameters
    let uncorrected_content = fs::read_to_string(data_dir.join("kahler_param.dat"))
        .expect("Failed to read kahler_param.dat");

    let uncorrected: Vec<f64> = uncorrected_content
        .trim()
        .split(',')
        .map(|s| s.trim().parse::<f64>().expect("Invalid Kähler parameter"))
        .collect();

    // Load corrected (with instanton corrections) Kähler parameters
    let corrected_content = fs::read_to_string(data_dir.join("corrected_kahler_param.dat"))
        .expect("Failed to read corrected_kahler_param.dat");

    let corrected: Vec<f64> = corrected_content
        .trim()
        .split(',')
        .map(|s| {
            s.trim()
                .parse::<f64>()
                .expect("Invalid corrected Kähler parameter")
        })
        .collect();

    (uncorrected, corrected)
}

#[test]
fn stage9_load_kahler_params() {
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let (uncorrected, corrected) = load_kahler_params(&data_dir);

    // Both should have 214 parameters (one per basis divisor)
    assert_eq!(uncorrected.len(), 214, "Should have 214 Kähler parameters");
    assert_eq!(
        corrected.len(),
        214,
        "Should have 214 corrected Kähler parameters"
    );

    // Check some specific values from the data
    let uncorrected_sum: f64 = uncorrected.iter().sum();
    let corrected_sum: f64 = corrected.iter().sum();

    #[derive(serde::Serialize)]
    struct KahlerParamsSummary {
        n_params: usize,
        uncorrected_first_5: Vec<f64>,
        corrected_first_5: Vec<f64>,
        uncorrected_sum: f64,
        corrected_sum: f64,
        ratio: f64,
    }

    let summary = KahlerParamsSummary {
        n_params: uncorrected.len(),
        uncorrected_first_5: uncorrected
            .iter()
            .take(5)
            .map(|&v| round_to_decimals(v, 6))
            .collect(),
        corrected_first_5: corrected
            .iter()
            .take(5)
            .map(|&v| round_to_decimals(v, 6))
            .collect(),
        uncorrected_sum: round_to_decimals(uncorrected_sum, 4),
        corrected_sum: round_to_decimals(corrected_sum, 4),
        ratio: round_to_decimals(corrected_sum / uncorrected_sum, 4),
    };

    insta::assert_json_snapshot!("kahler_params_summary", summary);
}

/// Compute V_string from primal intersection numbers and Kähler parameters
/// This is a first-principles test that computes intersection numbers from scratch
#[test]
fn stage10_compute_v_string() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    use cyrus_core::types::i32::I32;
    use cyrus_core::types::tags::{GTEOne, NonNeg};
    use cyrus_core::volume::bbhl_correction;
    use cyrus_core::{
        Point, Polytope, compute_regular_triangulation, intersection::compute_intersection_cytools,
        intersection_in_basis,
    };

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // === Load polytope and compute triangulation ===

    #[derive(Debug, Deserialize)]
    struct PolytopeInput {
        points: Vec<Vec<i64>>,
    }

    #[derive(Debug, Deserialize)]
    struct HeightsInput {
        values: Vec<f64>,
    }

    let points_raw = if let Some(dir) = crate::mcallister_data_dir() {
        let content =
            std::fs::read_to_string(dir.join("points.dat")).expect("Failed to read points.dat");
        content
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|line| {
                line.split(',')
                    .map(|s| s.trim().parse::<i64>().expect("Invalid point value"))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<Vec<i64>>>()
    } else {
        assert!(
            crate::fixtures_enabled(),
            "Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures"
        );
        let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
        let content = std::fs::read_to_string(&input_path).expect("Failed to read polytope.json");
        let input: PolytopeInput =
            serde_json::from_str(&content).expect("Failed to parse polytope.json");
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

    let heights = if let Some(dir) = crate::mcallister_data_dir() {
        let content =
            std::fs::read_to_string(dir.join("heights.dat")).expect("Failed to read heights.dat");
        content
            .split([',', '\n', '\r'])
            .filter(|s| !s.trim().is_empty())
            .map(|s| s.trim().parse::<f64>().expect("Invalid height value"))
            .collect::<Vec<f64>>()
    } else {
        assert!(
            crate::fixtures_enabled(),
            "Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures"
        );
        let heights_path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
        let heights_content =
            std::fs::read_to_string(&heights_path).expect("Failed to read heights.json");
        let heights_input: HeightsInput =
            serde_json::from_str(&heights_content).expect("Failed to parse heights.json");
        heights_input.values
    };

    // Compute triangulation from their heights
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights)
        .expect("Failed to compute triangulation");

    // === Compute intersection numbers using CYTools algorithm ===

    // Compute linear relations (needed for CYTools algorithm)
    let points_as_vecs: Vec<Vec<i64>> = triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linear_relations_int =
        cyrus_core::integer_math::compute_linear_relations_no_origin(&points_as_vecs);
    let linear_relations_no_origin: Vec<Vec<i64>> = linear_relations_int
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("fits in i64"))
                .collect()
        })
        .collect();

    let kappa_full = compute_intersection_cytools(
        &triangulation,
        &triangulation_points,
        &linear_relations_no_origin,
    )
    .expect("Failed to compute intersection numbers");

    let basis_path = data_dir.join("basis.dat");
    let basis_i64 = read_csv_i64(&basis_path);
    let basis: Vec<usize> = basis_i64
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index fits usize"))
        .collect();

    eprintln!("Using McAllister's primal basis: {} indices", basis.len());
    eprintln!("First 10 basis indices: {:?}", &basis[..10]);

    // Extract intersection numbers in basis
    let kappa = intersection_in_basis(&kappa_full, &basis);

    eprintln!(
        "Intersection numbers: dim={}, non-zero={}",
        kappa.dim(),
        kappa.num_nonzero()
    );

    // === Load Kähler parameters ===

    let (uncorrected_t, corrected_t) = load_kahler_params(&data_dir);

    // Verify dimensions match
    assert_eq!(
        kappa.dim(),
        basis.len(),
        "Intersection dim {} should match basis size {}",
        kappa.dim(),
        basis.len()
    );
    assert_eq!(
        corrected_t.len(),
        basis.len(),
        "Kähler params {} should match basis size {}",
        corrected_t.len(),
        basis.len()
    );

    // === Compute classical volumes from Kähler parameters ===

    let volume_from_t = |t: &[f64]| -> f64 {
        let mut volume_sum = 0.0f64;
        for (&(i, j, k), val) in kappa.iter() {
            // Convert rational to f64
            let (kappa_val, _): (f64, _) =
                malachite::num::conversion::traits::RoundingFrom::rounding_from(
                    val.get(),
                    malachite::rounding_modes::RoundingMode::Nearest,
                );

            let t_product = t[i] * t[j] * t[k];

            // Account for symmetry multiplicity
            let mult = if i == j && j == k {
                1.0
            } else if i == j || j == k || i == k {
                3.0
            } else {
                6.0
            };

            volume_sum += mult * kappa_val * t_product;
        }

        volume_sum / 6.0
    };

    let classical_volume_uncorrected = volume_from_t(&uncorrected_t);
    let classical_volume_corrected = volume_from_t(&corrected_t);

    // === Apply BBHL correction ===

    let h11 = I32::<GTEOne>::new(214).unwrap();
    let h21 = I32::<NonNeg>::new(4).unwrap();
    let bbhl = bbhl_correction(h11, h21);

    let v_string = classical_volume_corrected - bbhl.get();

    eprintln!(
        "Classical volume (uncorrected): {}",
        classical_volume_uncorrected
    );
    eprintln!(
        "Classical volume (corrected): {}",
        classical_volume_corrected
    );
    eprintln!("BBHL correction: {}", bbhl.get());
    eprintln!("V_string (computed): {}", v_string);

    // === Compare against McAllister's value ===

    let assertion = load_racetrack_assertion();
    let expected_v_string = assertion.cy_vol;

    eprintln!("V_string (expected): {}", expected_v_string);

    // Allow some tolerance due to floating point and possible basis differences
    let tolerance = 50.0; // ~1% of 4711
    assert!(
        (v_string - expected_v_string).abs() < tolerance,
        "V_string {} differs from expected {} by more than {}",
        v_string,
        expected_v_string,
        tolerance
    );

    // Compare against McAllister's cy_vol.dat (V_string without worldsheet instantons).
    let cy_vol_expected: f64 = std::fs::read_to_string(data_dir.join("cy_vol.dat"))
        .expect("Failed to read cy_vol.dat")
        .trim()
        .parse()
        .expect("Invalid cy_vol.dat");
    let dat_tol = 1e-6;
    assert!(
        (v_string - cy_vol_expected).abs() < dat_tol,
        "V_string mismatch: computed={}, expected={}",
        v_string,
        cy_vol_expected
    );

    #[derive(serde::Serialize)]
    struct VStringSummary {
        kappa_dim: usize,
        kappa_nonzero: usize,
        basis_size: usize,
        classical_volume_uncorrected: f64,
        classical_volume_corrected: f64,
        bbhl_correction: f64,
        v_string_computed: f64,
        v_string_expected: f64,
        difference: f64,
        difference_percent: f64,
    }

    // Avoid -0.0 in snapshot by normalizing near-zero values
    let diff = round_to_decimals(v_string - expected_v_string, 4);
    let diff_pct = round_to_decimals(
        100.0 * (v_string - expected_v_string) / expected_v_string,
        4,
    );

    let summary = VStringSummary {
        kappa_dim: kappa.dim(),
        kappa_nonzero: kappa.num_nonzero(),
        basis_size: basis.len(),
        classical_volume_uncorrected: round_to_decimals(classical_volume_uncorrected, 4),
        classical_volume_corrected: round_to_decimals(classical_volume_corrected, 4),
        bbhl_correction: round_to_decimals(bbhl.get(), 6),
        v_string_computed: round_to_decimals(v_string, 4),
        v_string_expected: round_to_decimals(expected_v_string, 4),
        difference: if diff.abs() < 1e-10 { 0.0 } else { diff },
        difference_percent: if diff_pct.abs() < 1e-10 {
            0.0
        } else {
            diff_pct
        },
    };

    insta::assert_json_snapshot!("v_string_computation", summary);
}

/// Compute V₀ from McAllister's published values
/// V₀ = -3 × e^K₀ × (g_s⁷ / (4×V_string)²) × W₀²
#[test]
fn stage11_cosmological_constant() {
    // Load all the physics values from McAllister
    let assertion = load_racetrack_assertion();

    // McAllister values (Table 6.4 of arXiv:2107.09064)
    let e_k0 = 0.234393; // From stage 4 flat direction
    let g_s = assertion.g_s; // 0.00911134
    let v_string = assertion.cy_vol; // 4711.829675...
    let w_0 = assertion.w_0; // 2.30012e-90

    eprintln!("=== McAllister Physics Values ===");
    eprintln!("e^K₀ = {}", e_k0);
    eprintln!("g_s = {}", g_s);
    eprintln!("V_string = {}", v_string);
    eprintln!("W₀ = {:e}", w_0);

    // Compute V₀ = -3 × e^K₀ × (g_s⁷ / (4×V_string)²) × W₀²
    //
    // Breaking it down:
    // - g_s^7 = (0.00911134)^7 ≈ 5.05e-15
    // - (4×V_string)² = (4×4711.83)² ≈ 3.55e8
    // - g_s^7 / (4×V_string)² ≈ 1.42e-23
    // - W₀² = (2.30e-90)² ≈ 5.29e-180
    // - e^K₀ × ratio × W₀² ≈ 0.234 × 1.42e-23 × 5.29e-180 ≈ 1.76e-203
    // - V₀ = -3 × 1.76e-203 ≈ -5.3e-203 Mpl⁴

    let g_s_7 = g_s.powi(7);
    let denominator = (4.0 * v_string).powi(2);
    let ratio = g_s_7 / denominator;
    let w_0_squared = w_0 * w_0;
    let v_0 = -3.0 * e_k0 * ratio * w_0_squared;

    eprintln!("\n=== V₀ Computation Breakdown ===");
    eprintln!("g_s^7 = {:e}", g_s_7);
    eprintln!("(4×V_string)² = {:e}", denominator);
    eprintln!("g_s^7 / (4×V_string)² = {:e}", ratio);
    eprintln!("W₀² = {:e}", w_0_squared);
    eprintln!("V₀ = {:e} Mpl⁴", v_0);

    // Expected: V₀ ≈ -5.5 × 10⁻²⁰³ Mpl⁴ (McAllister section 6.4)
    // This is an extremely small negative value - successful moduli stabilization!

    // Verify the exponent is correct (around -203)
    let v_0_exponent = v_0.abs().log10().round() as i32;
    eprintln!("\nV₀ exponent: 10^{}", v_0_exponent);

    assert!(v_0 < 0.0, "V₀ should be negative (AdS vacuum)");
    assert!(
        (-210..=-195).contains(&v_0_exponent),
        "V₀ exponent should be around -203, got {}",
        v_0_exponent
    );

    #[derive(serde::Serialize)]
    struct V0Summary {
        e_k0: f64,
        g_s: f64,
        v_string: f64,
        w_0_exponent: i32,
        g_s_7: f64,
        denominator: f64,
        ratio: f64,
        v_0_exponent: i32,
        v_0_sign: String,
        note: String,
    }

    let summary = V0Summary {
        e_k0: round_to_decimals(e_k0, 6),
        g_s: round_to_decimals(g_s, 8),
        v_string: round_to_decimals(v_string, 4),
        w_0_exponent: w_0.abs().log10().round() as i32,
        g_s_7,
        denominator: round_to_decimals(denominator, 2),
        ratio,
        v_0_exponent,
        v_0_sign: "negative (AdS)".to_string(),
        note: "V₀ ≈ -5.5 × 10⁻²⁰³ Mpl⁴ - extremely small CC achieved!".to_string(),
    };

    insta::assert_json_snapshot!("v0_cosmological_constant", summary);
}

#[test]
fn stage11_first_principles_runner_reaches_corrected_volume_and_v0() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };

    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let runner = std::env::var_os("CYRUS_MCALLISTER_RUNNER_BIN").map_or_else(
        || workspace_root.join("target/release/mcallister_first_principles"),
        PathBuf::from,
    );
    if !runner.exists() {
        eprintln!(
            "Skipping full first-principles runner test (build release runner or set CYRUS_MCALLISTER_RUNNER_BIN)"
        );
        return;
    }
    let output = Command::new(&runner)
        .current_dir(&workspace_root)
        .arg("--data-dir")
        .arg(&data_dir)
        .arg("--branch-candidates")
        .arg("1")
        .arg("--branch-height-init")
        .arg("--branch-selection")
        .arg("first-positive")
        .arg("--kklt-steps")
        .arg("64")
        .output()
        .unwrap_or_else(|e| panic!("failed to run {}: {e}", runner.display()));

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "first-principles runner failed with status {:?}\nstderr:\n{}",
        output.status.code(),
        stderr
    );
    assert!(
        stderr.contains("primal small toric curve GVs selected=344"),
        "runner did not compute the expected 344 small toric GV values:\n{stderr}"
    );
    assert!(
        stderr.contains("GV volume correction ="),
        "runner did not compute the GV volume correction:\n{stderr}"
    );
    assert!(
        !stderr.contains("transforming Kähler parameters from source basis"),
        "runner appears to have replayed downstream corrected Kähler parameters:\n{stderr}"
    );

    let v_string = extract_result_value(&stderr, "V_string");
    let v0_log10_abs = extract_result_value(&stderr, "log10(|V0|)");
    let corrected_checkpoint = read_scalar_f64(&data_dir.join("corrected_cy_vol.dat"));
    let uncorrected_checkpoint = read_scalar_f64(&data_dir.join("cy_vol.dat"));

    assert!(
        (v_string - corrected_checkpoint).abs() < 0.1,
        "runner V_string {v_string} should match corrected_cy_vol.dat {corrected_checkpoint}"
    );
    assert!(
        (v_string - corrected_checkpoint).abs() < (v_string - uncorrected_checkpoint).abs(),
        "runner V_string should be closer to corrected than uncorrected checkpoint"
    );
    assert!(
        (-203.0..-201.0).contains(&v0_log10_abs),
        "runner log10(|V0|) should be near -202, got {v0_log10_abs}"
    );
}
