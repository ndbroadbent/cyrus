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
