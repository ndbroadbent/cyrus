// Test-specific clippy allows
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::uninlined_format_args)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::redundant_closure_for_method_calls)]
#![allow(clippy::unnecessary_sort_by)]
#![allow(clippy::needless_for_each)]
#![allow(clippy::iter_cloned_collect)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::float_cmp)]
#![allow(clippy::struct_field_names)]
#![allow(clippy::suboptimal_flops)]
#![allow(clippy::redundant_clone)]
#![allow(clippy::stable_sort_primitive)]
#![allow(clippy::bool_to_int_with_if)]
#![allow(dead_code)]
#![allow(unused_variables)]

//! McAllister E2E Pipeline Tests
//!
//! End-to-end validation against arXiv:2107.09064.
//! Each stage computes ONE thing from first principles.
//!
//! Stage 1: Polytope points (raw primitive)
//! Stage 2: Triangulation (from points + heights)
//! Stage 3: GLSM charges + Intersection numbers (from triangulation)
//! Stage 4: Flat direction (from intersection + flux)
//! Stage 5: Racetrack solution (from GV + flat direction)
//! Stage 6: W₀ (from racetrack)
//! Stage 7: V_string (from moduli)
//! Stage 8: V₀ (from all above)

pub(crate) fn first_principles_enabled() -> bool {
    std::env::var_os("CYRUS_FIRST_PRINCIPLES").is_some()
}

pub(crate) fn fixtures_enabled() -> bool {
    let enabled = std::env::var_os("CYRUS_ALLOW_FIXTURES").is_some();
    assert!(
        !(enabled && first_principles_enabled()),
        "CYRUS_ALLOW_FIXTURES cannot be used with CYRUS_FIRST_PRINCIPLES"
    );
    enabled
}

pub(crate) fn mcallister_data_dir() -> Option<std::path::PathBuf> {
    if let Some(dir) = std::env::var_os("CYRUS_MCALLISTER_DATA_DIR") {
        return Some(std::path::PathBuf::from(dir));
    }

    let default = std::path::PathBuf::from(
        "/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647",
    );
    if default.exists() {
        Some(default)
    } else {
        None
    }
}

#[path = "mcallister_e2e/stage0_data_integrity.rs"]
mod stage0_data_integrity;

#[path = "mcallister_e2e/stage1_polytope.rs"]
mod stage1_polytope;

#[path = "mcallister_e2e/stage2_triangulation.rs"]
mod stage2_triangulation;

#[path = "mcallister_e2e/stage3_intersection.rs"]
mod stage3_intersection;

#[path = "mcallister_e2e/stage4_flat_direction.rs"]
mod stage4_flat_direction;

#[path = "mcallister_e2e/stage5_racetrack.rs"]
mod stage5_racetrack;

#[path = "mcallister_e2e/stage5_gv.rs"]
mod stage5_gv;

#[path = "mcallister_e2e/stage10_volume.rs"]
mod stage10_volume;

#[path = "mcallister_e2e/stage11_verify_orchestrator.rs"]
mod stage11_verify_orchestrator;

#[path = "mcallister_e2e/stage11_bench.rs"]
mod stage11_bench;
