//! McAllister data integrity checks.
//!
//! Verifies that extracted JSON fixtures match the original `.dat` files
//! from the McAllister dataset.
//!
//! JSON fixtures are deprecated and only checked when both
//! `CYRUS_ALLOW_FIXTURES=1` and `CYRUS_MCALLISTER_DATA_DIR` are set.

#![allow(missing_docs)]

use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct DualPointsFixture {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct DualSimplicesFixture {
    simplices: Vec<Vec<usize>>,
}

#[derive(Debug, Deserialize)]
struct BasisFixture {
    indices: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct TargetVolumesFixture {
    #[serde(rename = "c_i")]
    values: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct CurvesFixture {
    curve_classes: Vec<Vec<i64>>,
    gv_invariants: Vec<i64>,
}

#[derive(Debug, Deserialize)]
struct DualLinearRelationsFixture {
    linear_relations_no_origin: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct FluxInput {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

fn read_csv_f64(path: &str) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn read_csv_i64(path: &str) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
        .collect()
}

fn read_csv_rows_i64(path: &str) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
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

fn read_csv_rows_usize(path: &str) -> Vec<Vec<usize>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
    content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<usize>().expect("invalid integer"))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn require_data_dir() -> Option<PathBuf> {
    let Some(dir) = crate::mcallister_data_dir() else {
        if crate::first_principles_enabled() {
            panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
        }
        if crate::fixtures_enabled() {
            panic!("CYRUS_MCALLISTER_DATA_DIR must be set to validate JSON fixtures");
        }
        eprintln!("Skipping data integrity checks (set CYRUS_MCALLISTER_DATA_DIR)");
        return None;
    };
    Some(dir)
}

fn fixtures_enabled() -> bool {
    if !crate::fixtures_enabled() {
        eprintln!("Skipping JSON fixture checks (set CYRUS_ALLOW_FIXTURES=1)");
        return false;
    }
    true
}

#[test]
fn stage0_heights_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: HeightsInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("heights.dat");
    let dat = read_csv_f64(dat_path.to_str().expect("valid heights.dat path"));

    assert_eq!(json.values.len(), dat.len(), "heights length mismatch");
    for (i, (a, b)) in json.values.iter().zip(dat.iter()).enumerate() {
        let diff = (a - b).abs();
        assert!(
            diff < 1e-8,
            "heights mismatch at index {i}: json={a}, dat={b}"
        );
    }
}

#[test]
fn stage0_polytope_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("points.dat");
    let dat = read_csv_rows_i64(dat_path.to_str().expect("valid points.dat path"));

    assert_eq!(json.points.len(), dat.len(), "points length mismatch");
    for (i, (a, b)) in json.points.iter().zip(dat.iter()).enumerate() {
        assert_eq!(a, b, "point row {i} mismatch");
    }
}

#[test]
fn stage0_dual_points_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_points.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: DualPointsFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("dual_points.dat");
    let dat = read_csv_rows_i64(dat_path.to_str().expect("valid dual_points.dat path"));

    assert_eq!(json.points.len(), dat.len(), "dual points length mismatch");
    for (i, (a, b)) in json.points.iter().zip(dat.iter()).enumerate() {
        assert_eq!(a, b, "dual point row {i} mismatch");
    }
}

#[test]
fn stage0_dual_simplices_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_simplices.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: DualSimplicesFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("dual_simplices.dat");
    let dat = read_csv_rows_usize(dat_path.to_str().expect("valid dual_simplices.dat path"));

    assert_eq!(json.simplices.len(), dat.len(), "dual simplices length mismatch");
    for (i, (a, b)) in json.simplices.iter().zip(dat.iter()).enumerate() {
        assert_eq!(a, b, "dual simplex row {i} mismatch");
    }
}

#[test]
fn stage0_dual_linear_relations_fixture_matches_recompute() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/overrides/dual_linear_relations.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: DualLinearRelationsFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dual_points = read_csv_rows_i64(
        data_dir
            .join("dual_points.dat")
            .to_str()
            .expect("valid dual_points.dat path"),
    );
    let dual_points_i64: Vec<Vec<i64>> = dual_points.into_iter().take(9).collect();
    let recomputed = cyrus_core::compute_linear_relations_no_origin(&dual_points_i64)
        .into_iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("linrel fits i64"))
                .collect::<Vec<i64>>()
        })
        .collect::<Vec<Vec<i64>>>();

    assert_eq!(
        json.linear_relations_no_origin.len(),
        recomputed.len(),
        "dual linear relations row count mismatch"
    );
    for (i, (a, b)) in json
        .linear_relations_no_origin
        .iter()
        .zip(recomputed.iter())
        .enumerate()
    {
        assert_eq!(a, b, "dual linear relations row {i} mismatch");
    }
}

#[test]
fn stage0_dual_basis_matches_expected() {
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let dual_points =
        read_csv_rows_i64(data_dir.join("dual_points.dat").to_str().expect("valid dual_points.dat path"));
    let dual_points_vec: Vec<cyrus_core::Point> = dual_points
        .into_iter()
        .take(9)
        .map(cyrus_core::Point::new)
        .collect();

    let (glsm, _linrel, basis) =
        cyrus_core::compute_glsm_and_linrels(&dual_points_vec)
            .expect("Failed to compute dual GLSM/basis");

    let expected = vec![3, 4, 5, 8];
    if basis == expected {
        return;
    }

    if std::env::var_os("CYRUS_STRICT_BASIS").is_some() {
        assert_eq!(
            basis,
            expected,
            "Computed dual basis must match McAllister [3,4,5,8]"
        );
    }

    let transform = cyrus_core::basis_change_matrix(&glsm, &basis, &expected)
        .expect("Failed to compute dual basis change matrix");
    assert!(
        cyrus_core::is_unimodular(&transform),
        "Dual basis differs from McAllister, and change-of-basis is not unimodular"
    );
    eprintln!("Dual basis differs from McAllister; unimodular change-of-basis confirmed.");
}

#[test]
#[ignore = "primal_basis.json is deprecated; basis should come from basis.dat"]
fn stage0_primal_basis_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/overrides/primal_basis.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: BasisFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("basis.dat");
    let dat: Vec<i64> = read_csv_i64(dat_path.to_str().expect("valid basis.dat path"));
    let dat_usize: Vec<usize> = dat
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index fits usize"))
        .collect();

    assert_eq!(json.indices, dat_usize, "primal basis mismatch");
}

#[test]
fn stage0_kklt_basis_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/assertions/kklt_basis.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: BasisFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("kklt_basis.dat");
    let dat: Vec<i64> = read_csv_i64(dat_path.to_str().expect("valid kklt_basis.dat path"));
    let dat_usize: Vec<usize> = dat
        .into_iter()
        .map(|v| usize::try_from(v).expect("basis index fits usize"))
        .collect();

    assert_eq!(json.indices, dat_usize, "kklt basis mismatch");
}

#[test]
fn stage0_target_volumes_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/overrides/target_volumes.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: TargetVolumesFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let dat_path = data_dir.join("target_volumes.dat");
    let dat = read_csv_i64(dat_path.to_str().expect("valid target_volumes.dat path"));
    assert_eq!(json.values, dat, "target volumes mismatch");
}

#[test]
fn stage0_small_curves_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/overrides/curves.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: CurvesFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let curves_path = data_dir.join("small_curves.dat");
    let curves = read_csv_rows_i64(curves_path.to_str().expect("valid small_curves.dat path"));

    let gv_path = data_dir.join("small_curves_gv.dat");
    let gv = read_csv_i64(gv_path.to_str().expect("valid small_curves_gv.dat path"));

    assert_eq!(json.curve_classes.len(), curves.len(), "curve count mismatch");
    assert_eq!(json.gv_invariants.len(), gv.len(), "gv count mismatch");
    assert_eq!(json.curve_classes, curves, "curve classes mismatch");
    assert_eq!(json.gv_invariants, gv, "gv invariants mismatch");
}

#[test]
fn stage0_flux_fixture_matches_dat() {
    if !fixtures_enabled() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let json_path = manifest_dir.join("tests/mcallister_e2e/inputs/flux.json");
    let content = std::fs::read_to_string(&json_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", json_path.display()));
    let json: FluxInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", json_path.display()));

    let k_path = data_dir.join("K_vec.dat");
    let k = read_csv_i64(k_path.to_str().expect("valid K_vec.dat path"));
    let m_path = data_dir.join("M_vec.dat");
    let m = read_csv_i64(m_path.to_str().expect("valid M_vec.dat path"));

    assert_eq!(json.k, k, "K flux mismatch");
    assert_eq!(json.m, m, "M flux mismatch");
}
