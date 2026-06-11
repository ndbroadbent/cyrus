//! Chamber-search hook at the missing-GV failure site.
//!
//! When the selected small curves exceed GV coverage, this hook (enabled by
//! `--chamber-search-steps <n>`) walks neighbor chambers via
//! `cyrus_core::chamber_search`, gating each on phase-1 KKLT reachability,
//! and emits certified interior heights for the best chamber found. The
//! caller re-runs the verification with the emitted heights - the search is
//! advisory, the re-run is the proof.

use cyrus_core::chamber_search::search_covered_chamber;
use cyrus_core::gv::CurvePruningStrategy;
use cyrus_core::height_kahler::{effective_prime_divisors_from_curve_basis, heights_to_kahler};
use cyrus_core::kklt::{
    scale_divisor_basis_kklt_branch_initialization_to_target, solve_divisor_basis_path_following,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::range::CheckedRange;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::{Triangulation, compute_curve_basis_matrix};

use crate::{
    OwnedDivisorBasis, PrimalGeom, PrimalIntersection, chamber_intersection_full, parse_arg_value,
    transform_computed_primal_kahler_to_production,
};

/// Everything the failure-site hook needs from the volume stage.
pub struct ChamberSearchContext<'a> {
    pub missing_gv_classes: &'a [Vec<i64>],
    pub geom: &'a PrimalGeom,
    pub intersection: &'a PrimalIntersection,
    pub production_basis: &'a OwnedDivisorBasis,
    pub kklt_basis: &'a [usize],
    pub c_i: &'a [I64<Pos>],
    pub kklt_steps: usize,
    pub small_curve_selection_t: &'a [F64<Finite>],
    pub small_curve_cutoff: F64<Pos>,
    pub small_curve_pruning: CurvePruningStrategy,
    pub data_dir: Option<&'a str>,
}

/// No-op when every selected curve has a GV value; otherwise report the
/// failure, run the chamber search when enabled, and exit loudly.
pub fn enforce_gv_coverage(ctx: &ChamberSearchContext<'_>) {
    if ctx.missing_gv_classes.is_empty() {
        return;
    }
    eprintln!(
        "[ERROR] missing computed GV values for {} selected primal small toric curves; first_missing={:?}",
        ctx.missing_gv_classes.len(),
        ctx.missing_gv_classes.first()
    );
    run_chamber_search_on_missing_gv(ctx);
    std::process::exit(2);
}

/// Phase-1 KKLT reachability for a neighbor chamber: recompute its
/// intersection numbers, project its interior heights to a scaled branch
/// initialization, and require the divisor-basis path-following to
/// converge. A chamber that covers every curve but cannot reach the
/// solution is useless.
#[allow(clippy::too_many_arguments)] // closure context for the search
fn neighbor_chamber_phase1_converges(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    tau_phase1: &[F64<Pos>],
    kklt_steps: usize,
    neighbor_tri: &Triangulation,
    neighbor_heights: &[f64],
) -> bool {
    let Ok(kappa_neighbor) = chamber_intersection_full(neighbor_tri, &geom.triangulation_points)
    else {
        return false;
    };
    let Ok(basis_non_origin) = intersection
        .basis
        .iter()
        .map(|&idx| idx.checked_sub(1).ok_or(()))
        .collect::<Result<Vec<_>, ()>>()
    else {
        return false;
    };
    let Ok(curve_basis) = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
    else {
        return false;
    };
    let Some(prime_divisors) =
        effective_prime_divisors_from_curve_basis(&curve_basis, &basis_non_origin)
    else {
        return false;
    };
    let heights_typed: Vec<F64<Finite>> = neighbor_heights
        .iter()
        .map(|&h| F64::<Finite>::new(h))
        .collect::<Option<Vec<_>>>()
        .unwrap_or_default();
    if heights_typed.is_empty() {
        return false;
    }
    let Some(raw) = heights_to_kahler(&heights_typed, &basis_non_origin, &prime_divisors) else {
        return false;
    };
    let Ok(production_raw) = transform_computed_primal_kahler_to_production(
        intersection,
        production_basis,
        &raw,
        "neighbor-chamber height projection",
    ) else {
        return false;
    };
    let Some(init) = scale_divisor_basis_kklt_branch_initialization_to_target(
        &kappa_neighbor,
        production_basis.as_divisor_basis(),
        kklt_basis,
        tau_phase1,
        &production_raw,
    ) else {
        return false;
    };
    solve_divisor_basis_path_following(
        &kappa_neighbor,
        production_basis.as_divisor_basis(),
        kklt_basis,
        tau_phase1,
        &init,
        CheckedRange::new(0, kklt_steps),
    )
    .is_some_and(|result| result.converged)
}

/// Walk neighbor chambers for one whose own two-face formulas cover its
/// selected small curves, gating on phase-1 reachability, and emit the
/// winner's certified heights. Enabled by --chamber-search-steps <n>.
fn run_chamber_search_on_missing_gv(ctx: &ChamberSearchContext<'_>) {
    let &ChamberSearchContext {
        missing_gv_classes: _,
        geom,
        intersection,
        production_basis,
        kklt_basis,
        c_i,
        kklt_steps,
        small_curve_selection_t,
        small_curve_cutoff,
        small_curve_pruning,
        data_dir,
    } = ctx;
    let Some(max_expansions) = parse_arg_value::<usize>("--chamber-search-steps") else {
        eprintln!(
            "[HINT] pass --chamber-search-steps <n> to search neighbor chambers for full GV coverage"
        );
        return;
    };
    eprintln!(
        "[INFO] chamber search: walking bistellar flips (budget {max_expansions} expansions) scoring two-face GV coverage"
    );
    let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
    let reachable = |neighbor_tri: &Triangulation, neighbor_heights: &[f64]| {
        neighbor_chamber_phase1_converges(
            geom,
            intersection,
            production_basis,
            kklt_basis,
            &tau_phase1,
            kklt_steps,
            neighbor_tri,
            neighbor_heights,
        )
    };
    let outcome = search_covered_chamber(
        &geom.polytope,
        &geom.triangulation_points,
        geom.triangulation.simplices(),
        &intersection.basis,
        small_curve_selection_t,
        small_curve_cutoff,
        small_curve_pruning,
        max_expansions,
        &reachable,
    )
    .unwrap_or_else(|e| {
        eprintln!("[ERROR] chamber search failed: {e}");
        std::process::exit(2);
    });
    eprintln!(
        "[INFO] chamber search: start_uncovered={} best_uncovered={} flips={} expanded={}",
        outcome.start_uncovered, outcome.uncovered, outcome.flips_from_start, outcome.expanded
    );
    if outcome.flips_from_start == 0 {
        eprintln!("[ERROR] chamber search found no better chamber within budget");
        return;
    }
    let Some(dir) = data_dir else {
        eprintln!("[ERROR] chamber search needs --data-dir to emit heights");
        return;
    };
    let path = std::path::PathBuf::from(dir).join("heights_chamber_search.dat");
    let csv = outcome
        .heights
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    std::fs::write(&path, csv).unwrap_or_else(|e| {
        eprintln!("[ERROR] failed to write {}: {e}", path.display());
        std::process::exit(2);
    });
    eprintln!(
        "[INFO] chamber search: wrote candidate heights to {} ({} uncovered remain in that chamber)",
        path.display(),
        outcome.uncovered
    );
    eprintln!(
        "[INFO] to verify in the new chamber: cp {} {}/heights.dat and re-run",
        path.display(),
        dir
    );
}
