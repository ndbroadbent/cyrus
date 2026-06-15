//! Stage 11 (in-process): drive `verify_kklt_vacuum` directly.
//!
//! The heavy `stage10_volume::stage11_first_principles_runner_reaches_*` test
//! runs the whole `mcallister_first_principles` binary as a subprocess. This
//! test instead exercises the library orchestrator IN-PROCESS: it computes the
//! primal geometry and intersection numbers from first principles (the same
//! production code the binary's triangulation/intersection stages call), then
//! calls `cyrus_core::kklt_vacuum::verify_kklt_vacuum` and checks it reproduces
//! McAllister 4-214's corrected volume and `V0`.
//!
//! The orchestrator's *inputs* `g_s`, `|W0|`, `e^{K0}`, `c_i` and the KKLT
//! basis are NOT computed by the orchestrator - they are produced by the
//! upstream flat-direction/racetrack stages (validated by stage4/stage5) and
//! are passed in here from the published checkpoint. The orchestrator's own
//! scope - the corrected target, branch solve, small-curve GV coverage, the
//! instanton-corrected solve, and the BBHL/GV-corrected volume down to `V0` -
//! is computed live, so a regression in any of those fails this test.

use std::path::{Path, PathBuf};

use cyrus_core::kklt_vacuum::{
    BranchSelection, PrimalGeom, PrimalIntersection, StabilizationInputs, VacuumConfig,
    computed_primal_basis, verify_kklt_vacuum,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::{
    CurvePruningStrategy, Point, Polytope, compute_glsm_and_linrels, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_regular_triangulation, intersection_in_basis,
};

/// Published `e^{K0}` for 4-214-647 (arXiv:2107.09064 SS6.4). The checkpoint
/// ships no `e^{K0}` scalar file; `e^{K0}` is an orchestrator input from the
/// flat-direction stage (stage4) and only enters `V0`, not the volume.
const MCALLISTER_4_214_EK0: f64 = 0.234_393;
/// `h^{2,1}` for 4-214-647 (so `χ = 2(h11 - h21) = 2(214 - 4) = 420`).
const MCALLISTER_4_214_H21: usize = 4;

fn require_first_principles() -> bool {
    if !crate::first_principles_enabled() {
        eprintln!("Skipping first-principles test (set CYRUS_FIRST_PRINCIPLES=1)");
        return false;
    }
    true
}

fn require_runner_heavy() -> bool {
    if std::env::var_os("CYRUS_MCALLISTER_RUNNER_HEAVY").is_none() {
        eprintln!("Skipping heavy orchestrator test (set CYRUS_MCALLISTER_RUNNER_HEAVY=1)");
        return false;
    }
    true
}

fn require_data_dir() -> Option<PathBuf> {
    let Some(dir) = crate::mcallister_data_dir() else {
        assert!(
            !crate::first_principles_enabled(),
            "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
        );
        return None;
    };
    Some(dir)
}

fn read_floats_csv(path: &Path) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn read_ints_csv(path: &Path) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
        .collect()
}

fn read_usize_csv(path: &Path) -> Vec<usize> {
    read_ints_csv(path)
        .into_iter()
        .map(|v| usize::try_from(v).expect("index must be non-negative"))
        .collect()
}

fn read_scalar_f64(path: &Path) -> f64 {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()))
        .trim()
        .parse::<f64>()
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn read_points(path: &Path) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<i64>().expect("invalid point entry"))
                .collect()
        })
        .collect()
}

/// Build the primal geometry from `points.dat` + `heights.dat` (production
/// triangulation code).
fn build_geom(data_dir: &Path) -> PrimalGeom {
    let points_raw = read_points(&data_dir.join("points.dat"));
    let heights_raw = read_floats_csv(&data_dir.join("heights.dat"));
    let heights: Vec<F64<Finite>> = heights_raw
        .iter()
        .map(|&h| F64::<Finite>::new(h).expect("height must be finite"))
        .collect();
    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("failed to build polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter triangulation points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights_raw)
        .expect("failed to compute triangulation");
    PrimalGeom {
        polytope,
        heights,
        triangulation_points,
        triangulation,
    }
}

/// Build the primal intersection data (GLSM/linrels/basis + intersection
/// tensors) from the geometry (production intersection code).
fn build_intersection(geom: &PrimalGeom) -> PrimalIntersection {
    let (glsm, linrels, basis) =
        compute_glsm_and_linrels(&geom.triangulation_points).expect("failed GLSM/linrels");
    let points_i64: Vec<Vec<i64>> = geom
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linrels_i64: Vec<Vec<i64>> = compute_linear_relations_no_origin(&points_i64)
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("linrel fits in i64"))
                .collect()
        })
        .collect();
    let kappa_full = compute_intersection_cytools(
        &geom.triangulation,
        &geom.triangulation_points,
        &linrels_i64,
    )
    .expect("failed intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);
    PrimalIntersection {
        glsm,
        linrels,
        basis,
        kappa_full,
        kappa_basis,
    }
}

/// Cross-example characterization of McAllister's `kklt_basis` selection.
///
/// Per the paper SS2 (the F-flat algorithm), the basis is "a set of h^{1,1}
/// linearly independent RIGID prime toric divisors such that, with their
/// volumes set to unity, the remaining four prime toric divisors have strictly
/// positive volume" - i.e. tau_* = (1,...,1) lies in the dual effective cone.
/// This prints, for whichever example `CYRUS_MCALLISTER_DATA_DIR` points at,
/// the rigid/non-rigid split and the four "remaining" (non-basis) divisors, to
/// confirm `kklt_basis` is always a subset of the rigid divisors and the
/// remaining four are the non-rigid ones plus the dropped rigid ones.
#[test]
fn characterize_kklt_basis_selection() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);
    let (classes, _components) = cyrus_core::orientifold::classify_divisor_rigidity(
        &geom.polytope,
        &geom.triangulation_points,
    )
    .expect("rigidity classification");
    let rigid: std::collections::HashSet<usize> = classes
        .iter()
        .filter(|(_, c)| c.is_rigid())
        .map(|(idx, _)| *idx)
        .collect();
    let non_origin: Vec<usize> = (0..geom.triangulation_points.len())
        .filter(|&i| {
            geom.triangulation_points[i]
                .coords()
                .iter()
                .any(|&c| c != 0)
        })
        .collect();
    let kklt: std::collections::HashSet<usize> = read_usize_csv(&data_dir.join("kklt_basis.dat"))
        .into_iter()
        .collect();

    let non_rigid: Vec<usize> = non_origin
        .iter()
        .copied()
        .filter(|i| !rigid.contains(i))
        .collect();
    let remaining: Vec<usize> = non_origin
        .iter()
        .copied()
        .filter(|i| !kklt.contains(i))
        .collect();
    let dropped_rigid: Vec<usize> = remaining
        .iter()
        .copied()
        .filter(|i| rigid.contains(i))
        .collect();
    let kklt_subset_rigid = kklt.iter().all(|i| rigid.contains(i));

    eprintln!(
        "[CHAR] {}: h11_toric={} n_prime={} n_rigid={} n_nonrigid={} kklt_len={}",
        data_dir.file_name().unwrap().to_string_lossy(),
        intersection.basis.len(),
        non_origin.len(),
        rigid.len(),
        non_rigid.len(),
        kklt.len()
    );
    eprintln!("[CHAR]   kklt is subset of rigid: {kklt_subset_rigid}");
    eprintln!("[CHAR]   non-rigid divisors: {non_rigid:?}");
    eprintln!("[CHAR]   remaining (non-basis) divisors: {remaining:?}");
    eprintln!("[CHAR]   of which rigid (dropped from basis): {dropped_rigid:?}");
    assert!(
        kklt_subset_rigid,
        "kklt_basis must be a subset of rigid divisors"
    );
    assert_eq!(
        non_origin.len(),
        intersection.basis.len() + 4,
        "there are h^{{1,1}}+4 prime toric divisors (paper SS2)"
    );
    assert_eq!(
        remaining.len(),
        4,
        "exactly four prime divisors are non-basis"
    );
    // The remaining four are precisely the non-rigid divisors plus the
    // dropped-rigid ones - no other divisors are excluded.
    assert_eq!(
        non_rigid.len() + dropped_rigid.len(),
        4,
        "remaining four = non-rigid + dropped-rigid"
    );

    // The paper's effective-cone criterion must accept McAllister's own basis.
    assert!(
        cyrus_core::orientifold::complement_is_admissible(&geom.triangulation_points, &remaining),
        "McAllister's basis complement must satisfy tau_* in E^0"
    );

    // select_kklt_basis derives an admissible basis from scratch.
    let derived =
        cyrus_core::orientifold::select_kklt_basis(&geom.polytope, &geom.triangulation_points)
            .expect("an admissible KKLT basis exists");
    eprintln!(
        "[CHAR]   select_kklt_basis len={} == McAllister kklt_basis.dat: {}",
        derived.len(),
        derived
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
            == kklt
    );
    assert_eq!(
        derived.len(),
        intersection.basis.len(),
        "derived basis has h11 divisors"
    );
    assert!(
        derived.iter().all(|i| rigid.contains(i)),
        "derived basis is all-rigid"
    );
}

/// The effective-cone criterion (paper SS2) must separate the admissible basis
/// from the inadmissible ones the V0-invariance experiment tried - cheaply,
/// with no KKLT solve. 4-214 only.
#[test]
fn kklt_basis_criterion_separates_good_and_bad() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    if data_dir.file_name().unwrap().to_string_lossy() != "4-214-647" {
        return;
    }
    let geom = build_geom(&data_dir);
    let pts = &geom.triangulation_points;
    let admissible = |c: &[usize]| cyrus_core::orientifold::complement_is_admissible(pts, c);

    for c in [vec![1, 2, 46, 130], vec![1, 2, 3, 4], vec![1, 2, 217, 218]] {
        let (det, nums) = cyrus_core::orientifold::complement_remaining_volumes(pts, &c);
        let tau: Vec<f64> = nums.iter().map(|&n| n as f64 / det as f64).collect();
        eprintln!(
            "[CRIT] complement {c:?}: det={det} admissible={} tau_C={tau:?}",
            admissible(&c)
        );
    }
    // McAllister's {46,130}: all remaining volumes positive (det=2, NOT
    // unimodular - independence is what matters, not unimodularity).
    assert!(
        admissible(&[1, 2, 46, 130]),
        "McAllister complement is admissible"
    );
    // {3,4} drops a structural divisor -> divisor 3 has negative volume.
    assert!(
        !admissible(&[1, 2, 3, 4]),
        "dropping structural {{3,4}} is inadmissible (negative remaining volume)"
    );
    // {217,218}: divisor 218 has negative volume - so the criterion DOES reject
    // the basis whose solve landed in the GV-uncoverable 685-curve chamber.
    assert!(
        !admissible(&[1, 2, 217, 218]),
        "dropping {{217,218}} is inadmissible (negative remaining volume)"
    );
}

/// De-risking the GA deep-verify integration: how much of the orientifold
/// model that `verify_kklt_vacuum` needs (the involution, `c_i`, `kklt_basis`)
/// can be DERIVED from the geometry rather than read from McAllister's data
/// files?
///
/// Established here for 4-214-647:
/// 1. `enumerate_involutions` finds the correct involution `sigma=[1,0,0,0]`
///    with 51 O7 divisors (matching the runner's own gamma derivation).
/// 2. The `c_i` rule - `6` on O7-parity divisors, else the irreducible-
///    component count - reproduces `target_volumes.dat` EXACTLY on
///    McAllister's `kklt_basis`.
/// 3. Basis SELECTION is the derivable rule below (RESOLVED via paper SS2).
///    For 4-214 the rigid prime toric divisors number 216 and `kklt_basis` is
///    those minus exactly `{46, 130}` (both rigid; `{1, 2}` are non-rigid).
///
/// Resolution - the paper SS2 F-flat algorithm gives the selection rule
/// explicitly: among the h^{1,1}+4 prime toric divisors, `kklt_basis` is a set
/// of h^{1,1} linearly independent RIGID divisors such that, with their volumes
/// set to unity, the volumes of the remaining four divisors are strictly
/// positive - i.e. tau_* = (1,...,1) lies in the dual effective cone E(X)^0.
/// Confirmed across all five examples (see `characterize_kklt_basis_selection`):
/// the "remaining four" are always exactly the non-rigid divisors plus the
/// dropped-rigid ones. This is a brute-force search over rigid bases with a
/// concrete positivity criterion, NOT a hand-picked convention.
///
/// The earlier V0-invariance experiment failures are now explained by this
/// criterion: dropping `{3, 4}` (4-214) puts a remaining divisor at negative
/// volume, so tau_* leaves E(X)^0 and the F-flat solve has no positive-volume
/// branch (exactly the paper's autochthonous-negative-volume failure mode);
/// dropping `{217, 218}` happens to keep tau_* in E^0 but lands the solve in a
/// different, GV-uncoverable chamber. So the GA's deep-verify derives the basis
/// by the paper's criterion (rigid + remaining-four-positive at tau_*).
#[test]
fn derived_orientifold_model_vs_mcallister_data() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);

    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&geom.polytope, &geom.triangulation_points)
            .expect("batyrev h11 extra");
    let h11 = intersection.basis.len() + h11_extra;

    let models = cyrus_core::orientifold::enumerate_involutions(
        &geom.polytope,
        &geom.triangulation_points,
        &intersection.kappa_full,
        &intersection.basis,
        h11,
        MCALLISTER_4_214_H21,
    )
    .expect("enumerate involutions");

    // (1) The correct involution is found.
    let model = models
        .iter()
        .find(|m| m.sigma == vec![1, 0, 0, 0])
        .expect("sigma=[1,0,0,0] model present");
    assert_eq!(model.o7_points.len(), 51, "expected 51 O7 divisors");

    let data_kklt_basis: Vec<usize> = read_usize_csv(&data_dir.join("kklt_basis.dat"));
    let data_c_i: Vec<i64> = read_ints_csv(&data_dir.join("target_volumes.dat"));

    let c_i_rule = |idx: usize| -> i64 {
        if model.gamma[idx] == 1 {
            6
        } else {
            i64::try_from(model.components[idx]).unwrap()
        }
    };

    // (2) The c_i rule reproduces target_volumes.dat on McAllister's basis.
    let c_i_on_data: Vec<i64> = data_kklt_basis.iter().map(|&idx| c_i_rule(idx)).collect();
    assert_eq!(
        c_i_on_data, data_c_i,
        "derived c_i rule must reproduce target_volumes.dat on McAllister's kklt_basis"
    );

    // (3) Basis-selection gap: rigid prime divisors minus McAllister's basis
    // is exactly {46, 130}, and every McAllister divisor is rigid.
    let mut rigid: Vec<usize> = model
        .divisor_classes
        .iter()
        .filter(|(_, class)| class.is_rigid())
        .map(|(idx, _)| *idx)
        .collect();
    rigid.sort_unstable();
    let dropped: Vec<usize> = rigid
        .iter()
        .filter(|i| !data_kklt_basis.contains(i))
        .copied()
        .collect();
    let extra: Vec<usize> = data_kklt_basis
        .iter()
        .filter(|i| !rigid.contains(i))
        .copied()
        .collect();
    assert_eq!(dropped, vec![46, 130], "rigid-minus-McAllister drop set");
    assert!(extra.is_empty(), "every McAllister kklt divisor is rigid");
}

/// Which admissible drop-pair maximizes the minimum remaining divisor volume
/// (the "most interior" tau_*), and is it McAllister's {46,130}? Cheap - just
/// 4x4 determinants over all rigid drop-pairs, no KKLT solve. 4-214 only.
#[test]
fn most_interior_admissible_basis_4214() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    if data_dir.file_name().unwrap().to_string_lossy() != "4-214-647" {
        return;
    }
    let geom = build_geom(&data_dir);
    let pts = &geom.triangulation_points;
    let (classes, _) =
        cyrus_core::orientifold::classify_divisor_rigidity(&geom.polytope, pts).expect("rigidity");
    let mut rigid: Vec<usize> = classes
        .iter()
        .filter(|(_, c)| c.is_rigid())
        .map(|(i, _)| *i)
        .collect();
    rigid.sort_unstable();

    let min_remaining_vol = |a: usize, b: usize| -> Option<f64> {
        let (det, nums) = cyrus_core::orientifold::complement_remaining_volumes(pts, &[1, 2, a, b]);
        if det == 0 || !nums.iter().all(|&n| n != 0 && (n > 0) == (det > 0)) {
            return None;
        }
        Some(
            nums.iter()
                .map(|&n| n as f64 / det as f64)
                .fold(f64::INFINITY, f64::min),
        )
    };

    let mut scored: Vec<((usize, usize), f64)> = Vec::new();
    for (i, &a) in rigid.iter().enumerate() {
        for &b in &rigid[i + 1..] {
            if let Some(m) = min_remaining_vol(a, b) {
                scored.push(((a, b), m));
            }
        }
    }
    scored.sort_by(|x, y| y.1.total_cmp(&x.1));
    eprintln!(
        "[INTERIOR] admissible pairs={} of {} total",
        scored.len(),
        rigid.len() * (rigid.len() - 1) / 2
    );
    eprintln!("[INTERIOR] top 8 by min-remaining-volume:");
    for ((a, b), m) in scored.iter().take(8) {
        eprintln!("[INTERIOR]   drop {{{a},{b}}}: min_remaining_vol={m:.3}");
    }
    let mcallister = min_remaining_vol(46, 130).expect("McAllister admissible");
    let mcallister_rank = scored.iter().position(|((a, b), _)| *a == 46 && *b == 130);
    eprintln!(
        "[INTERIOR] McAllister {{46,130}}: min_remaining_vol={mcallister:.3} rank={mcallister_rank:?}"
    );

    // CYTools-style tiebreak: keep low-L1-norm divisors in the basis, drop the
    // high-norm ones. Rank admissible pairs by the dropped divisors' total L1
    // norm (descending) and see where McAllister's {46,130} lands.
    let l1 = |idx: usize| -> i64 { pts[idx].coords().iter().map(|c| c.abs()).sum() };
    let mut by_norm: Vec<((usize, usize), i64)> = scored
        .iter()
        .map(|&((a, b), _)| ((a, b), l1(a) + l1(b)))
        .collect();
    by_norm.sort_by(|x, y| y.1.cmp(&x.1));
    eprintln!("[NORM] top 8 admissible by dropped-divisor L1 norm (descending):");
    for ((a, b), n) in by_norm.iter().take(8) {
        eprintln!(
            "[NORM]   drop {{{a},{b}}}: dropped_norm={n} (norms {},{})",
            l1(*a),
            l1(*b)
        );
    }
    let mc_norm_rank = by_norm.iter().position(|((a, b), _)| *a == 46 && *b == 130);
    eprintln!(
        "[NORM] McAllister {{46,130}}: dropped_norm={} rank_by_norm={mc_norm_rank:?}",
        l1(46) + l1(130)
    );
    // Also ascending (drop lowest-norm) for comparison.
    by_norm.sort_by(|x, y| x.1.cmp(&y.1));
    let mc_norm_rank_asc = by_norm.iter().position(|((a, b), _)| *a == 46 && *b == 130);
    eprintln!("[NORM] McAllister rank ascending (drop low-norm)={mc_norm_rank_asc:?}");
}

fn require_general_gv() -> bool {
    if std::env::var_os("CYRUS_GENERAL_GV").is_none() {
        eprintln!("Skipping general-GV validation (set CYRUS_GENERAL_GV=1)");
        return false;
    }
    true
}

/// Build the orientifold involution model ONCE for a geometry: the involution
/// defining the O7 set is the one with the most O7-planes among the
/// pairwise-disjoint involutions (strongest racetrack hierarchy). Viability
/// against the GLSM basis is irrelevant - the KKLT basis is all-rigid by
/// construction; we only read gamma/components for c_i.
fn build_o7_model(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
) -> cyrus_core::orientifold::OrientifoldModel {
    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&geom.polytope, &geom.triangulation_points)
            .expect("batyrev h11 extra");
    let h11 = intersection.basis.len() + h11_extra;
    let models = cyrus_core::orientifold::enumerate_involutions(
        &geom.polytope,
        &geom.triangulation_points,
        &intersection.kappa_full,
        &intersection.basis,
        h11,
        MCALLISTER_4_214_H21,
    )
    .expect("enumerate involutions");
    models
        .into_iter()
        .filter(|m| m.o7_pairwise_disjoint && !m.o7_points.is_empty())
        .max_by_key(|m| m.o7_points.len())
        .expect("a disjoint-O7 involution exists")
}

/// Derive the c_i (dual Coxeter / instanton coefficients) for a derived KKLT
/// basis from a prebuilt involution model: 6 on O7-parity divisors (so(8)
/// gaugino condensation), else the irreducible-component count. This is the
/// rule `derived_orientifold_model_vs_mcallister_data` proves reproduces
/// `target_volumes.dat` on McAllister's own basis.
fn c_i_for_basis(
    model: &cyrus_core::orientifold::OrientifoldModel,
    kklt_basis: &[usize],
) -> Vec<I64<Pos>> {
    kklt_basis
        .iter()
        .map(|&idx| {
            let c = if model.gamma[idx] == 1 {
                6
            } else {
                i64::try_from(model.components[idx]).expect("component count fits i64")
            };
            I64::<Pos>::new(c).expect("c_i positive")
        })
        .collect()
}

/// END-TO-END BASIS-AGNOSTIC VALIDATION (the GA's deep-verify path).
///
/// `verify_kklt_vacuum_reaches_corrected_volume_and_v0` feeds McAllister's
/// hand-shipped `kklt_basis.dat` ({46,130} dropped). The GA has no such file: it
/// must DERIVE an admissible basis from the geometry (paper SS2 effective-cone
/// criterion). There are thousands of admissible bases, and they are NOT
/// physically interchangeable for the worldsheet-instanton expansion: each
/// places the phase-1 selection point in a different chamber, and only some
/// chambers select small curves the cheap toric GV formulas can cover. A
/// non-McAllister chamber can select a curve at GV grading degree ~700, where
/// the general HKTY-series fallback is intractable (verified empirically: the
/// least-interior basis needs degree 716).
///
/// The tractable, McAllister-faithful route is therefore NOT "any basis +
/// general GV" but a SCAN: walk the admissible bases (paper criterion) in
/// most-interior-first order and take the first whose chamber is fully covered
/// by the cheap toric + minimal-degree GV layers (`primal_gv_min_points=None`).
/// That chamber is McAllister's (or an equivalent low-degree one), and it
/// reproduces V0 ~ -202 with no expensive general-GV computation. This is the
/// exact procedure `GaGeometry::deep_verify` runs. It is slow (each rejected
/// basis costs one branch solve; McAllister's lands around interior-rank 436),
/// so it is gated behind `CYRUS_GENERAL_GV=1` and intended to run detached.
#[test]
fn derived_basis_auto_scan_reproduces_v0() {
    if !require_first_principles() || !require_runner_heavy() || !require_general_gv() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };

    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);
    let production_primal_basis = computed_primal_basis(&intersection);

    // Derive ALL admissible KKLT bases from first principles (NO kklt_basis.dat),
    // most-interior first, each with its own involution-derived c_i.
    let bases = cyrus_core::orientifold::admissible_kklt_bases(
        &geom.polytope,
        &geom.triangulation_points,
        0,
    )
    .expect("admissible KKLT bases exist");
    eprintln!(
        "[SCAN] {} admissible KKLT bases to scan (interior order)",
        bases.len()
    );
    let o7_model = build_o7_model(&geom, &intersection);
    eprintln!(
        "[SCAN] involution sigma={:?} o7_count={}",
        o7_model.sigma,
        o7_model.o7_points.len()
    );
    let candidates: Vec<(Vec<usize>, Vec<I64<Pos>>)> = bases
        .into_iter()
        .map(|basis| {
            let c_i = c_i_for_basis(&o7_model, &basis);
            (basis, c_i)
        })
        .collect();

    // g_s / |W0| are racetrack (upstream) outputs; e^{K0} is the flat-direction
    // factor. These are orchestrator inputs, validated by stage4/stage5.
    let g_s = F64::<Pos>::new(read_scalar_f64(&data_dir.join("g_s.dat"))).expect("g_s positive");
    let w0 = F64::<Pos>::new(read_scalar_f64(&data_dir.join("W_0.dat"))).expect("W0 positive");
    let ek0 = F64::<Pos>::new(MCALLISTER_4_214_EK0).expect("ek0 positive");
    // The GA has no small_curves_cutoff.dat; use the runner's default (1.0).
    let small_curve_cutoff = F64::<Pos>::new(1.0).expect("cutoff positive");

    let base_inputs = StabilizationInputs {
        geom: &geom,
        intersection: &intersection,
        production_primal_basis: &production_primal_basis,
        c_i: &candidates[0].1,
        kklt_basis: &candidates[0].0,
        g_s,
        w0,
        ek0,
        h21: MCALLISTER_4_214_H21,
    };
    let config = VacuumConfig {
        kklt_steps: 64,
        branch_candidates: 1,
        branch_seed: 42,
        branch_height_init: true,
        branch_selection: BranchSelection::FirstPositive,
        small_curve_cutoff,
        small_curve_pruning: CurvePruningStrategy::PairDecomposable,
        // Cheap coverage ONLY: the scan finds a chamber the toric + minimal
        // layers fully cover, so general GV is never needed.
        primal_gv_min_points: None,
        primal_gv_max_deg: None,
        max_missing_gv_impact: 0.0,
        missing_gv_abs_bound: 10,
    };

    let verdict =
        cyrus_core::kklt_vacuum::verify_kklt_vacuum_auto_basis(&base_inputs, &candidates, &config)
            .unwrap_or_else(|e| panic!("derived-basis auto-scan verify failed: {e}"));

    eprintln!(
        "[SCAN] V_string={} small_curve_gv_count={} gv_controlled={} deferred={}",
        verdict.v_string.get(),
        verdict.small_curve_gv_count,
        verdict.gv_controlled,
        verdict.deferred_missing_gv_count,
    );
    let v0_log10_abs = verdict.v0.get().abs().log10();
    eprintln!("[SCAN] log10(|V0|)={v0_log10_abs}");

    // WHICH chamber did the scan pick? Not necessarily McAllister's {46,130}:
    // any admissible basis lands on the same stabilized minimum (V0 is
    // basis-invariant); the scan just takes the first toric-coverable one. We
    // load kklt_basis.dat ONLY for this diagnostic comparison (not used in the
    // computation path).
    let mcallister: std::collections::HashSet<usize> =
        read_usize_csv(&data_dir.join("kklt_basis.dat"))
            .into_iter()
            .collect();
    let selected: std::collections::HashSet<usize> =
        verdict.selected_kklt_basis.iter().copied().collect();
    let dropped: Vec<usize> = (0..geom.triangulation_points.len())
        .filter(|i| {
            geom.triangulation_points[*i]
                .coords()
                .iter()
                .any(|&c| c != 0)
                && !selected.contains(i)
        })
        .collect();
    eprintln!(
        "[SCAN] selected basis dropped divisors={dropped:?}; same as McAllister {{46,130}}: {}",
        selected == mcallister
    );

    assert!(
        verdict.gv_controlled && verdict.deferred_missing_gv_count == 0,
        "scanned basis must be fully GV-controlled (no deferred curves)"
    );
    assert!(
        (-203.0..-201.0).contains(&v0_log10_abs),
        "derived-basis auto-scan log10(|V0|) should be near -202, got {v0_log10_abs}"
    );
}

#[test]
fn verify_kklt_vacuum_reaches_corrected_volume_and_v0() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };

    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);
    let production_primal_basis = computed_primal_basis(&intersection);

    let c_i: Vec<I64<Pos>> = read_ints_csv(&data_dir.join("target_volumes.dat"))
        .into_iter()
        .map(|v| I64::<Pos>::new(v).expect("target_volumes c_i must be positive"))
        .collect();
    let kklt_basis: Vec<usize> = read_usize_csv(&data_dir.join("kklt_basis.dat"));

    let g_s = F64::<Pos>::new(read_scalar_f64(&data_dir.join("g_s.dat"))).expect("g_s positive");
    let w0 = F64::<Pos>::new(read_scalar_f64(&data_dir.join("W_0.dat"))).expect("W0 positive");
    let ek0 = F64::<Pos>::new(MCALLISTER_4_214_EK0).expect("ek0 positive");
    let small_curve_cutoff =
        F64::<Pos>::new(read_scalar_f64(&data_dir.join("small_curves_cutoff.dat")))
            .expect("small-curve cutoff positive");

    let inputs = StabilizationInputs {
        geom: &geom,
        intersection: &intersection,
        production_primal_basis: &production_primal_basis,
        c_i: &c_i,
        kklt_basis: &kklt_basis,
        g_s,
        w0,
        ek0,
        h21: MCALLISTER_4_214_H21,
    };
    // Mirrors the validated runner invocation:
    // --branch-candidates 1 --branch-height-init --branch-selection first-positive --kklt-steps 64
    let config = VacuumConfig {
        kklt_steps: 64,
        branch_candidates: 1,
        branch_seed: 42,
        branch_height_init: true,
        branch_selection: BranchSelection::FirstPositive,
        small_curve_cutoff,
        small_curve_pruning: CurvePruningStrategy::PairDecomposable,
        primal_gv_min_points: None,
        primal_gv_max_deg: None,
        max_missing_gv_impact: 0.0,
        missing_gv_abs_bound: 10,
    };

    let verdict = verify_kklt_vacuum(&inputs, &config)
        .unwrap_or_else(|e| panic!("verify_kklt_vacuum failed: {e}"));

    assert_eq!(
        verdict.small_curve_gv_count, 344,
        "expected 344 selected small-curve GVs, got {}",
        verdict.small_curve_gv_count
    );
    assert!(
        verdict.gv_controlled && verdict.deferred_missing_gv_count == 0,
        "4-214 should be fully GV-controlled with no deferred curves"
    );

    let corrected_checkpoint = read_scalar_f64(&data_dir.join("corrected_cy_vol.dat"));
    let uncorrected_checkpoint = read_scalar_f64(&data_dir.join("cy_vol.dat"));
    let v_string = verdict.v_string.get();
    assert!(
        (v_string - corrected_checkpoint).abs() < 0.1,
        "orchestrator V_string {v_string} should match corrected_cy_vol.dat {corrected_checkpoint}"
    );
    assert!(
        (v_string - corrected_checkpoint).abs() < (v_string - uncorrected_checkpoint).abs(),
        "orchestrator V_string should be closer to corrected than uncorrected checkpoint"
    );

    let v0_log10_abs = verdict.v0.get().abs().log10();
    assert!(
        (-203.0..-201.0).contains(&v0_log10_abs),
        "orchestrator log10(|V0|) should be near -202, got {v0_log10_abs}"
    );
}
