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
    let content =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn read_csv_i64(path: &str) -> Vec<i64> {
    let content =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
    content
        .split(|c| c == ',' || c == '\n' || c == '\r')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
        .collect()
}

fn read_csv_rows_i64(path: &str) -> Vec<Vec<i64>> {
    let content =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
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
    let content =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {}: {e}", path));
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

fn runner_heavy_enabled() -> bool {
    if std::env::var_os("CYRUS_MCALLISTER_RUNNER_HEAVY").is_none() {
        eprintln!(
            "Skipping full first-principles runner test (set CYRUS_MCALLISTER_RUNNER_HEAVY=1)"
        );
        return false;
    }
    true
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ArtifactUse {
    DeclaredInput,
    ValidationCheckpoint,
    ValidationReplayOnly,
    TemporaryBridge,
}

#[derive(Clone, Copy, Debug)]
struct ArtifactPolicy {
    file: &'static str,
    usage: ArtifactUse,
    note: &'static str,
}

const ARTIFACT_POLICIES: &[ArtifactPolicy] = &[
    ArtifactPolicy {
        file: "points.dat",
        usage: ArtifactUse::DeclaredInput,
        note: "polytope vertices for the validation target",
    },
    ArtifactPolicy {
        file: "heights.dat",
        usage: ArtifactUse::DeclaredInput,
        note: "selected regular triangulation heights",
    },
    ArtifactPolicy {
        file: "K_vec.dat",
        usage: ArtifactUse::DeclaredInput,
        note: "integer K flux vector",
    },
    ArtifactPolicy {
        file: "M_vec.dat",
        usage: ArtifactUse::DeclaredInput,
        note: "integer M flux vector",
    },
    ArtifactPolicy {
        file: "kklt_basis.dat",
        usage: ArtifactUse::DeclaredInput,
        note: "declared KKLT/non-perturbative divisor set",
    },
    ArtifactPolicy {
        file: "target_volumes.dat",
        usage: ArtifactUse::DeclaredInput,
        note: "declared c_i targets from the KKLT/orientifold choice",
    },
    ArtifactPolicy {
        file: "dual_simplices.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed dual FRST checkpoint",
    },
    ArtifactPolicy {
        file: "dual_points.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "dual polytope points must be computed by Cyrus and compared here",
    },
    ArtifactPolicy {
        file: "basis.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed primal divisor basis checkpoint",
    },
    ArtifactPolicy {
        file: "dual_curves.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed dual curve classes checkpoint",
    },
    ArtifactPolicy {
        file: "dual_curves_gv.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed dual GV invariants checkpoint",
    },
    ArtifactPolicy {
        file: "small_curves.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed racetrack curve subset checkpoint",
    },
    ArtifactPolicy {
        file: "small_curves_cutoff.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "checkpoint for McAllister's small-curve cutoff; production runner takes an explicit cutoff parameter",
    },
    ArtifactPolicy {
        file: "small_curves_gv.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed racetrack GV subset checkpoint",
    },
    ArtifactPolicy {
        file: "small_curves_vols.dat",
        usage: ArtifactUse::ValidationReplayOnly,
        note: "downstream small-curve volumes; only used to pin checkpoint semantics",
    },
    ArtifactPolicy {
        file: "g_s.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed racetrack string coupling checkpoint",
    },
    ArtifactPolicy {
        file: "W_0.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "computed racetrack superpotential checkpoint",
    },
    ArtifactPolicy {
        file: "kahler_param.dat",
        usage: ArtifactUse::ValidationReplayOnly,
        note: "downstream KKLT output; never a first-principles computation input",
    },
    ArtifactPolicy {
        file: "corrected_kahler_param.dat",
        usage: ArtifactUse::ValidationReplayOnly,
        note: "downstream corrected KKLT output; gated for checkpoint replay only",
    },
    ArtifactPolicy {
        file: "corrected_heights.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "corrected-chamber FRST checkpoint for diagnostics, never a production chamber input",
    },
    ArtifactPolicy {
        file: "corrected_target_volumes.dat",
        usage: ArtifactUse::ValidationCheckpoint,
        note: "corrected-chamber classical KKLT divisor-volume checkpoint",
    },
    ArtifactPolicy {
        file: "cy_vol.dat",
        usage: ArtifactUse::ValidationReplayOnly,
        note: "downstream volume output; compare only",
    },
    ArtifactPolicy {
        file: "corrected_cy_vol.dat",
        usage: ArtifactUse::ValidationReplayOnly,
        note: "downstream corrected volume output; compare only",
    },
];

fn artifact_policy(file: &str) -> Option<ArtifactPolicy> {
    ARTIFACT_POLICIES
        .iter()
        .copied()
        .find(|policy| policy.file == file)
}

fn dat_tokens(source: &str) -> impl Iterator<Item = &str> {
    source
        .split(|c: char| !(c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.'))
        .filter(|token| {
            token.ends_with(".dat")
                && token[..token.len() - 4]
                    .bytes()
                    .any(|b| b.is_ascii_alphanumeric())
        })
}

fn source_region_before_fn<'a>(source: &'a str, start_fn: &str, end_fn: &str) -> &'a str {
    let start_marker = format!("fn {start_fn}(");
    let end_marker = format!("\nfn {end_fn}(");
    let start = source
        .find(&start_marker)
        .unwrap_or_else(|| panic!("missing function {start_fn}"));
    let end = source[start..]
        .find(&end_marker)
        .unwrap_or_else(|| panic!("missing function {end_fn} after {start_fn}"));
    &source[start..start + end]
}

#[test]
fn stage0_artifact_policy_is_explicit_and_complete() {
    use std::collections::{BTreeMap, BTreeSet};

    let mut policy_names = BTreeSet::new();
    for policy in ARTIFACT_POLICIES {
        assert!(
            policy_names.insert(policy.file),
            "duplicate artifact policy for {}",
            policy.file
        );
        assert!(
            !policy.note.is_empty(),
            "{} needs a policy note",
            policy.file
        );
    }

    assert_eq!(
        artifact_policy("points.dat").unwrap().usage,
        ArtifactUse::DeclaredInput
    );
    assert_eq!(
        artifact_policy("corrected_kahler_param.dat").unwrap().usage,
        ArtifactUse::ValidationReplayOnly
    );

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let source_roots = [
        manifest_dir.join("tests/mcallister_e2e"),
        manifest_dir.join("src/bin"),
    ];

    let mut referenced = BTreeMap::<String, BTreeSet<String>>::new();
    let mut stack = source_roots.to_vec();
    while let Some(path) = stack.pop() {
        if path.is_dir() {
            for entry in std::fs::read_dir(&path)
                .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()))
            {
                stack.push(entry.expect("valid directory entry").path());
            }
            continue;
        }
        if path.extension().and_then(|ext| ext.to_str()) != Some("rs") {
            continue;
        }
        let source = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
        for token in dat_tokens(&source) {
            referenced.entry(token.to_string()).or_default().insert(
                path.strip_prefix(&manifest_dir)
                    .unwrap_or(&path)
                    .display()
                    .to_string(),
            );
        }
    }

    for file in referenced.keys() {
        assert!(
            artifact_policy(file).is_some(),
            "{file} is referenced but missing from ARTIFACT_POLICIES"
        );
    }

    for policy in ARTIFACT_POLICIES {
        assert!(
            referenced.contains_key(policy.file),
            "{} has a policy but is not referenced by McAllister sources",
            policy.file
        );
    }
}

#[test]
fn stage0_first_principles_runner_does_not_silently_replay_downstream_outputs() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let runner_path = manifest_dir.join("src/bin/mcallister_first_principles.rs");
    let runner = std::fs::read_to_string(&runner_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", runner_path.display()));
    let stage_volume = source_region_before_fn(&runner, "stage_volume", "compare_against_dat");
    let stage_volume_tokens = dat_tokens(stage_volume).collect::<std::collections::BTreeSet<_>>();

    for forbidden in ["kahler_param.dat", "cy_vol.dat", "corrected_cy_vol.dat"] {
        assert!(
            !stage_volume_tokens.contains(forbidden),
            "{forbidden} must not be read by the first-principles volume computation"
        );
    }

    assert!(
        stage_volume.contains("corrected_kahler_param.dat"),
        "runner should name the currently unresolved corrected Kahler checkpoint"
    );
    assert!(
        stage_volume.contains("--allow-downstream-kahler"),
        "corrected_kahler_param.dat replay must require an explicit validation-only flag"
    );
    assert!(
        stage_volume.contains("if allow_downstream_kahler"),
        "corrected_kahler_param.dat replay must be isolated behind allow_downstream_kahler"
    );

    let compare = source_region_before_fn(&runner, "compare_against_dat", "stage_vacuum");
    let compare_tokens = dat_tokens(compare).collect::<std::collections::BTreeSet<_>>();
    assert!(
        compare_tokens.contains("corrected_cy_vol.dat"),
        "corrected_cy_vol.dat should only be used as a validation comparison checkpoint"
    );
    assert!(
        !compare_tokens.contains("corrected_kahler_param.dat"),
        "validation comparison must not load corrected Kähler parameters"
    );

    let stage_vacuum = source_region_before_fn(&runner, "stage_vacuum", "run_pipeline");
    assert!(
        !stage_volume.contains("compare_against_dat("),
        "volume computation must not invoke downstream checkpoint comparison"
    );
    assert!(
        stage_vacuum.contains("compare_against_dat("),
        "downstream checkpoint comparison belongs after V_string has been computed"
    );
    assert!(
        runner.contains("--skip-mcallister-assertions")
            && stage_vacuum.contains("if !validate_mcallister_assertions"),
        "McAllister final assertion checks must be explicitly skippable for non-validation candidates"
    );
    let stage_flat_direction =
        source_region_before_fn(&runner, "stage_flat_direction", "sparse_i64");
    assert!(
        stage_flat_direction.contains("use_mcallister_flux_basis_default")
            && stage_flat_direction.contains("using computed dual basis as flux coordinate basis")
            && stage_flat_direction.contains("using McAllister flux source basis [3, 4, 5, 8]"),
        "generic runs must not default to the McAllister flux source basis"
    );
    let validate_basis =
        source_region_before_fn(&runner, "validate_basis_checkpoint", "sorted_point_coords");
    assert!(
        validate_basis.contains("basis_path.exists()"),
        "basis.dat checkpoint validation must be optional after Cyrus computes the basis"
    );
    let validate_dual = source_region_before_fn(&runner, "validate_dual_checkpoint", "read_points");
    assert!(
        validate_dual.contains("dual_points_path.exists()")
            && validate_dual.contains("dual_simplices_path.exists()"),
        "dual_points.dat and dual_simplices.dat must be optional validation checkpoints after Cyrus computes dual geometry"
    );
    assert!(
        compare.contains("read_optional_scalar_f64(&dir.join(\"g_s.dat\"))")
            && compare.contains("read_optional_scalar_f64(&dir.join(\"W_0.dat\"))")
            && compare.contains("read_optional_scalar_f64(&corrected_volume_path)"),
        "g_s.dat, W_0.dat, and corrected_cy_vol.dat must be optional validation checkpoints, not required inputs"
    );
    assert!(
        compare.contains("corrected V_string comparison is not exact"),
        "non-exact corrected-volume comparisons must be logged as unresolved discrepancies"
    );

    for gv_artifact in [
        "dual_curves.dat",
        "dual_curves_gv.dat",
        "small_curves.dat",
        "small_curves_gv.dat",
    ] {
        assert!(
            !runner.contains(gv_artifact),
            "{gv_artifact} must be computed by the runner, not loaded"
        );
    }

    assert!(
        runner.contains(".compute_dual()"),
        "runner must compute the dual polytope instead of loading dual_points.dat as an input"
    );
    assert!(
        runner.contains("compute_frst_heights"),
        "runner must compute the dual FRST instead of loading dual_simplices.dat as an input"
    );
    assert!(
        !runner.contains("load_dual_inputs"),
        "runner must not keep the old dual_points/dual_simplices loader"
    );
    assert!(
        !runner.contains("overrides/dual_points.json")
            && !runner.contains("overrides/dual_simplices.json"),
        "runner must not fall back to legacy dual fixture overrides"
    );
}

#[test]
fn stage0_mcallister_binaries_do_not_use_validation_replay_artifacts() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bin_dir = manifest_dir.join("src/bin");
    let allowed_replay = (
        "src/bin/mcallister_first_principles.rs",
        "corrected_kahler_param.dat",
    );
    let forbidden_computed_artifacts = [
        "dual_curves.dat",
        "dual_curves_gv.dat",
        "small_curves.dat",
        "small_curves_gv.dat",
    ];

    for entry in std::fs::read_dir(&bin_dir)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", bin_dir.display()))
    {
        let path = entry.expect("valid directory entry").path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("rs") {
            continue;
        }
        let rel_path = path
            .strip_prefix(&manifest_dir)
            .unwrap_or(&path)
            .display()
            .to_string();
        let source = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
        let tokens = dat_tokens(&source).collect::<std::collections::BTreeSet<_>>();

        for artifact in forbidden_computed_artifacts {
            assert!(
                !tokens.contains(artifact),
                "{rel_path} must compute {artifact}; it may only be read by validation tests"
            );
        }

        for token in &tokens {
            let policy =
                artifact_policy(token).unwrap_or_else(|| panic!("{token} has no artifact policy"));
            if policy.usage != ArtifactUse::ValidationReplayOnly {
                continue;
            }
            match (rel_path.as_str(), *token) {
                path_token if path_token == allowed_replay => {
                    assert!(
                        source.contains("--allow-downstream-kahler"),
                        "{rel_path} must gate {token} behind --allow-downstream-kahler"
                    );
                }
                ("src/bin/mcallister_first_principles.rs", "corrected_cy_vol.dat") => {
                    let compare =
                        source_region_before_fn(&source, "compare_against_dat", "stage_vacuum");
                    assert!(
                        dat_tokens(compare).any(|compare_token| compare_token == *token),
                        "{rel_path} may only read {token} inside compare_against_dat"
                    );
                }
                _ => {
                    panic!(
                        "{rel_path} reads validation-only artifact {token} without an explicit exception"
                    );
                }
            }
        }
    }
}

#[test]
fn stage0_first_principles_runner_accepts_declared_inputs_only_data_dir() {
    if !crate::first_principles_enabled() || !runner_heavy_enabled() {
        return;
    }
    let Some(source_dir) = require_data_dir() else {
        return;
    };

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace_root = manifest_dir.join("../..");
    let runner = std::env::var_os("CYRUS_MCALLISTER_RUNNER_BIN")
        .map(PathBuf::from)
        .unwrap_or_else(|| workspace_root.join("target/release/mcallister_first_principles"));
    if !runner.exists() {
        eprintln!(
            "Skipping declared-input runner test (build release runner or set CYRUS_MCALLISTER_RUNNER_BIN)"
        );
        return;
    }

    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time should be after UNIX_EPOCH")
        .as_nanos();
    let temp_dir = std::env::temp_dir().join(format!(
        "cyrus-mcallister-declared-inputs-{}-{stamp}",
        std::process::id()
    ));
    std::fs::create_dir(&temp_dir)
        .unwrap_or_else(|e| panic!("Failed to create {}: {e}", temp_dir.display()));
    for file in [
        "points.dat",
        "heights.dat",
        "K_vec.dat",
        "M_vec.dat",
        "kklt_basis.dat",
        "target_volumes.dat",
    ] {
        std::fs::copy(source_dir.join(file), temp_dir.join(file))
            .unwrap_or_else(|e| panic!("Failed to copy declared input {file}: {e}"));
    }

    let output = std::process::Command::new(&runner)
        .current_dir(&workspace_root)
        .arg("--data-dir")
        .arg(&temp_dir)
        .arg("--kklt-steps")
        .arg("64")
        .output()
        .unwrap_or_else(|e| panic!("failed to run {}: {e}", runner.display()));
    let _ = std::fs::remove_dir_all(&temp_dir);

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "first-principles runner failed with status {:?}\nstderr:\n{}",
        output.status.code(),
        stderr
    );
    for checkpoint in [
        "basis.dat checkpoint not found",
        "dual_points.dat checkpoint not found",
        "dual_simplices.dat checkpoint not found",
        "g_s.dat checkpoint not found",
        "W_0.dat checkpoint not found",
        "corrected_cy_vol.dat checkpoint not found",
    ] {
        assert!(
            stderr.contains(checkpoint),
            "declared-input run should skip missing validation checkpoint {checkpoint}:\n{stderr}"
        );
    }
    assert!(
        stderr.contains("[RESULT] V_string = 4711.504666573377"),
        "declared-input run should reach the current no-replay V_string result:\n{stderr}"
    );
    assert!(
        stderr.contains("[RESULT] log10(|V0|) = -202.26279591106547"),
        "declared-input run should reach the current no-replay V0 result:\n{stderr}"
    );
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

    assert_eq!(
        json.simplices.len(),
        dat.len(),
        "dual simplices length mismatch"
    );
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
    let dual_points = read_csv_rows_i64(
        data_dir
            .join("dual_points.dat")
            .to_str()
            .expect("valid dual_points.dat path"),
    );
    let dual_points_vec: Vec<cyrus_core::Point> = dual_points
        .into_iter()
        .take(9)
        .map(cyrus_core::Point::new)
        .collect();

    let (glsm, _linrel, basis) = cyrus_core::compute_glsm_and_linrels(&dual_points_vec)
        .expect("Failed to compute dual GLSM/basis");

    let expected = vec![3, 4, 5, 8];
    if basis == expected {
        return;
    }

    if std::env::var_os("CYRUS_STRICT_BASIS").is_some() {
        assert_eq!(
            basis, expected,
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

    assert_eq!(
        json.curve_classes.len(),
        curves.len(),
        "curve count mismatch"
    );
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
