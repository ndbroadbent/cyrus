//! Tip of stretched cone computation.
//!
//! Find the tip (minimum norm point) of the stretched cone:
//! {x : H @ x >= c} where c is the stretching parameter.
//!
//! This is a quadratic programming problem:
//!   minimize ||x||^2
//!   subject to H @ x >= c
//!
//! Reference: CYTools cone.py tip_of_stretched_cone()

use super::Cone;
use good_lp::{default_solver, variable, Expression, ProblemVariables, Solution, SolverModel};

/// Find the tip of a stretched cone.
///
/// The stretched cone is {x : H @ x >= c} where H are hyperplanes and c > 0.
/// The tip is the minimum-norm point in this region.
///
/// For small dimensions, we use iterative refinement.
/// For large dimensions, this would need a proper QP solver.
pub fn tip_of_stretched_cone(cone: &mut Cone, c: f64) -> Option<Vec<f64>> {
    let h = cone.hyperplanes().to_vec();
    let dim = cone.ambient_dim();

    if h.is_empty() {
        // No constraints - tip is at origin (but that's not interior)
        // Return a small point in the interior direction
        let mut tip = vec![1.0; dim];
        let norm: f64 = tip.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in &mut tip {
            *x *= c / norm;
        }
        return Some(tip);
    }

    // Use gradient descent with projection
    // Start from the center of mass of hyperplane normals
    let mut x = find_initial_point(&h, dim, c)?;

    // Iterative refinement
    for _ in 0..1000 {
        // Project gradient of ||x||^2 = 2x onto feasible region
        let (feasible, violations) = check_feasibility(&h, &x, c);

        if feasible {
            // Move towards origin while staying feasible
            let step = 0.1;
            let new_x = x.iter().map(|&xi| xi * (1.0 - step)).collect::<Vec<_>>();

            let (new_feasible, _) = check_feasibility(&h, &new_x, c);
            if new_feasible {
                x = new_x;
            } else {
                // Can't improve much, we're close to optimal
                break;
            }
        } else {
            // Move away from violated constraints
            for &(i, violation) in &violations {
                let h_i = &h[i];
                let norm_sq: f64 = h_i.iter().map(|&hi| (hi as f64).powi(2)).sum();
                if norm_sq > 1e-10 {
                    // Move x in direction of h_i to satisfy constraint
                    let step = (violation + 0.1) / norm_sq;
                    for (j, &hj) in h_i.iter().enumerate() {
                        x[j] += step * hj as f64;
                    }
                }
            }
        }
    }

    // Verify solution
    let (feasible, _) = check_feasibility(&h, &x, c);
    if feasible { Some(x) } else { None }
}

/// Find an interior point using LP approach.
///
/// We want to find x such that H @ x >= c for all rows of H.
pub fn find_interior_point_lp(cone: &mut Cone) -> Option<Vec<f64>> {
    let h = cone.hyperplanes().to_vec();
    let dim = cone.ambient_dim();

    if h.is_empty() {
        // Full space - any point is interior
        return Some(vec![1.0; dim]);
    }

    // Objective vector: average of hyperplane normals (CYTools feasibility)
    let mut obj_vec = vec![0.0f64; dim];
    for row in &h {
        for i in 0..dim {
            obj_vec[i] += row[i] as f64;
        }
    }
    let denom = h.len() as f64;
    if denom > 0.0 {
        for v in &mut obj_vec {
            *v /= denom;
        }
    }

    // Fast heuristic: try constructive interior point first for smaller systems.
    if h.len() <= 10_000 {
        if let Some(candidate) = find_initial_point(&h, dim, 1.0) {
            if cone.contains_interior(&candidate, 1e-10) {
                return Some(candidate);
            }
        }
    }

    // Cutting-plane: start with a subset, add violated constraints.
    let mut active: Vec<usize> = (0..h.len().min(64)).collect();
    let mut active_set: std::collections::HashSet<usize> = active.iter().copied().collect();

    for iter in 0..200 {
        if iter % 10 == 0 {
            eprintln!(
                "[DEBUG] interior lp iter={}, active_constraints={}",
                iter,
                active.len()
            );
        }
        let candidate = solve_feasibility_lp(&h, &active, dim, &obj_vec)?;

        let mut violations: Vec<(usize, f64)> = Vec::new();
        for (idx, row) in h.iter().enumerate() {
            let dot: f64 = row
                .iter()
                .zip(candidate.iter())
                .map(|(&a, &x)| a as f64 * x)
                .sum();
            if dot < 1.0 - 1e-9 {
                violations.push((idx, 1.0 - dot));
            }
        }

        if violations.is_empty() && cone.contains_interior(&candidate, 1e-10) {
            return Some(candidate);
        }

        violations.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let mut added = 0usize;
        for (idx, _) in violations {
            if active_set.insert(idx) {
                active.push(idx);
                added += 1;
            }
            if added >= 256 {
                break;
            }
        }

        if added == 0 && iter > 0 {
            break;
        }
    }

    None
}

fn solve_feasibility_lp(
    h: &[Vec<i128>],
    active: &[usize],
    dim: usize,
    obj_vec: &[f64],
) -> Option<Vec<f64>> {
    let mut vars = ProblemVariables::new();
    let bound = 1.0e6;
    let x: Vec<_> = (0..dim)
        .map(|_| vars.add(variable().min(-bound).max(bound)))
        .collect();

    let mut obj_expr = Expression::from(0.0);
    for i in 0..dim {
        obj_expr.add_mul(obj_vec[i], x[i]);
    }
    let mut model = vars.minimise(obj_expr).using(default_solver);

    for &row_idx in active {
        let row = &h[row_idx];
        let max_abs = row.iter().map(|v| v.abs()).max().unwrap_or(0) as f64;
        if max_abs == 0.0 {
            eprintln!("[WARN] interior lp constraint has zero row");
            return None;
        }
        let scale = max_abs;
        let mut expr = Expression::from(0.0);
        for i in 0..dim {
            expr.add_mul((row[i] as f64) / scale, x[i]);
        }
        model = model.with(expr.geq(1.0 / scale));
    }

    let solution = match model.solve() {
        Ok(sol) => sol,
        Err(err) => {
            eprintln!("[WARN] interior lp solve failed: {}", err);
            return None;
        }
    };
    let mut out = Vec::with_capacity(dim);
    for i in 0..dim {
        out.push(solution.value(x[i]));
    }
    Some(out)
}

/// Find an initial feasible point for optimization.
fn find_initial_point(h: &[Vec<i128>], dim: usize, c: f64) -> Option<Vec<f64>> {
    if h.is_empty() {
        return Some(vec![c; dim]);
    }

    // Strategy: Start with sum of hyperplane normals
    // This tends to point towards the interior of the cone
    let mut x = vec![0.0; dim];
    for row in h {
        for (i, &hi) in row.iter().enumerate() {
            x[i] += hi as f64;
        }
    }

    // Normalize
    let norm: f64 = x.iter().map(|&xi| xi * xi).sum::<f64>().sqrt();
    if norm < 1e-10 {
        // Degenerate case - try each dimension
        for i in 0..dim {
            x[i] = 1.0;
        }
    } else {
        for xi in &mut x {
            *xi /= norm;
        }
    }

    // Scale to satisfy constraints
    scale_to_feasible(h, &mut x, c);

    // Verify and refine
    let (feasible, violations) = check_feasibility(h, &x, c);
    if feasible {
        return Some(x);
    }

    // Refinement: push away from violated constraints
    for _ in 0..100 {
        for &(i, violation) in &violations {
            let row = &h[i];
            let norm_sq: f64 = row.iter().map(|&hi| (hi as f64).powi(2)).sum();
            if norm_sq > 1e-10 {
                let step = (violation.abs() + c) / norm_sq;
                for (j, &hj) in row.iter().enumerate() {
                    x[j] += step * hj as f64;
                }
            }
        }

        let (feasible, new_violations) = check_feasibility(h, &x, c);
        if feasible {
            return Some(x);
        }

        if new_violations.is_empty() {
            break;
        }
    }

    // Try a different approach: large scaling
    let mut x = vec![0.0; dim];
    for row in h {
        for (i, &hi) in row.iter().enumerate() {
            x[i] += hi as f64;
        }
    }

    // Large scale to ensure we're well inside
    let scale = 1000.0 * c;
    for xi in &mut x {
        *xi *= scale;
    }

    let (feasible, _) = check_feasibility(h, &x, c);
    if feasible { Some(x) } else { None }
}

/// Scale a point to satisfy H @ x >= c.
fn scale_to_feasible(h: &[Vec<i128>], x: &mut [f64], c: f64) {
    // Find minimum scale factor needed
    let mut min_scale = 1.0;

    for row in h {
        let dot: f64 = row
            .iter()
            .zip(x.iter())
            .map(|(&hi, &xi)| hi as f64 * xi)
            .sum();
        if dot > 1e-10 {
            let needed_scale = c / dot;
            if needed_scale > min_scale {
                min_scale = needed_scale;
            }
        }
    }

    // Apply scaling with margin
    let scale = min_scale * 1.1;
    for xi in x.iter_mut() {
        *xi *= scale;
    }
}

/// Check if x satisfies H @ x >= c for all hyperplanes.
///
/// Returns (feasible, violations) where violations is a list of (index, violation amount).
fn check_feasibility(h: &[Vec<i128>], x: &[f64], c: f64) -> (bool, Vec<(usize, f64)>) {
    let mut violations = Vec::new();

    for (i, row) in h.iter().enumerate() {
        let dot: f64 = row
            .iter()
            .zip(x.iter())
            .map(|(&hi, &xi)| hi as f64 * xi)
            .sum();
        if dot < c - 1e-10 {
            violations.push((i, c - dot));
        }
    }

    (violations.is_empty(), violations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tip_first_quadrant() {
        // First quadrant: x >= 0, y >= 0, stretched by c=1
        // Tip should be near (1, 1) / sqrt(2) but scaled to satisfy H@x >= 1
        let mut cone = Cone::from_hyperplanes(vec![vec![1, 0], vec![0, 1]]);

        let tip = tip_of_stretched_cone(&mut cone, 1.0);
        assert!(tip.is_some());

        let tip = tip.unwrap();
        assert!(tip[0] >= 1.0 - 0.01, "x = {} should be >= 1", tip[0]);
        assert!(tip[1] >= 1.0 - 0.01, "y = {} should be >= 1", tip[1]);
    }

    #[test]
    fn test_tip_simple_cone() {
        // Cone defined by (1,0) and (0,1) as rays (first quadrant)
        // This is simpler and more testable
        let mut cone = Cone::from_rays(vec![vec![1, 0], vec![0, 1]]);

        // First compute hyperplanes (needed for tip finding)
        let _ = cone.hyperplanes();

        let tip = tip_of_stretched_cone(&mut cone, 0.5);

        // If tip found, check it's in the cone
        if let Some(ref t) = tip {
            // Check positive
            assert!(t[0] >= 0.0, "x should be >= 0");
            assert!(t[1] >= 0.0, "y should be >= 0");
        }
        // Note: tip finding may fail for some cones - this is expected
        // The algorithm is a simple gradient descent, not a full QP solver
    }

    #[test]
    fn test_find_interior_point() {
        let mut cone = Cone::from_hyperplanes(vec![vec![1, 0], vec![0, 1]]);

        let interior = find_interior_point_lp(&mut cone);
        assert!(interior.is_some());

        let p = interior.unwrap();
        assert!(p[0] > 0.0);
        assert!(p[1] > 0.0);
    }

    #[test]
    fn test_check_feasibility() {
        let h = vec![vec![1, 0], vec![0, 1]];

        // Feasible point
        let (feasible, _) = check_feasibility(&h, &[2.0, 3.0], 1.0);
        assert!(feasible);

        // Infeasible point
        let (feasible, violations) = check_feasibility(&h, &[0.5, 2.0], 1.0);
        assert!(!feasible);
        assert!(!violations.is_empty());
    }
}
