// Test-specific clippy allows
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::uninlined_format_args)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::cognitive_complexity)]
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
