//! Polyhedral cone computations.
//!
//! This module provides cone operations needed for Kähler/Mori cone computations.
//! Implements the Double Description Method (DDM) for cone dualization.
//!
//! Reference: CYTools cone.py, PPL (Parma Polyhedra Library)

mod ddm;
mod tip;

use std::collections::HashSet;
use std::env;
use std::hash::{Hash, Hasher};
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

use good_lp::{constraint, default_solver, variable, ProblemVariables, SolverModel};
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum LatticeC {
    Scalar(f64),
    List(Vec<f64>),
}

#[derive(Serialize, Deserialize)]
struct LatticeRequest {
    hyperplanes: Vec<Vec<i64>>,
    grading_vector: Vec<i64>,
    min_points: Option<usize>,
    max_deg: Option<i64>,
    max_coord: i64,
    deg_window: i64,
    max_deg_limit: Option<i64>,
    max_time_sec: Option<f64>,
    max_solutions: Option<usize>,
    num_search_workers: Option<i32>,
    strict: bool,
    c: LatticeC,
    check_grading: bool,
}

/// A rational polyhedral cone, represented by either rays or hyperplanes.
///
/// A cone can be defined in two equivalent ways:
/// - V-representation: As the positive span of rays: `{R^T @ λ : λ >= 0}`
/// - H-representation: As intersection of halfspaces: `{x : H @ x >= 0}`
///
/// The Double Description Method (DDM) converts between these representations.
#[derive(Debug, Clone)]
pub struct Cone {
    /// Dimension of the ambient space
    ambient_dim: usize,

    /// Rays (V-representation), if computed
    rays: Option<Vec<Vec<i128>>>,

    /// Hyperplane normals (H-representation), if computed
    hyperplanes: Option<Vec<Vec<i128>>>,

    /// Cached dimension
    dim: Option<usize>,
}

impl Cone {
    /// Create a cone from generating rays.
    ///
    /// The cone is the positive span: `{R^T @ λ : λ >= 0}` where R is the matrix
    /// with rays as rows.
    pub fn from_rays(rays: Vec<Vec<i128>>) -> Self {
        let ambient_dim = if rays.is_empty() { 0 } else { rays[0].len() };

        // Normalize rays by GCD
        let normalized = normalize_vectors(rays);

        Self {
            ambient_dim,
            rays: Some(normalized),
            hyperplanes: None,
            dim: None,
        }
    }

    /// Create a cone from hyperplane normals (inward-pointing).
    ///
    /// The cone is the intersection: `{x : H @ x >= 0}` where H is the matrix
    /// with hyperplane normals as rows.
    pub fn from_hyperplanes(hyperplanes: Vec<Vec<i128>>) -> Self {
        let ambient_dim = if hyperplanes.is_empty() {
            0
        } else {
            hyperplanes[0].len()
        };

        // Handle empty hyperplanes (full space)
        if hyperplanes.is_empty() {
            return Self::full_space(ambient_dim);
        }

        // Normalize hyperplanes by GCD
        let normalized = normalize_vectors(hyperplanes);

        Self {
            ambient_dim,
            rays: None,
            hyperplanes: Some(normalized),
            dim: None,
        }
    }

    /// Create the full space cone (no constraints).
    pub fn full_space(ambient_dim: usize) -> Self {
        // Full space has rays ±e_i for each coordinate
        let mut rays = Vec::with_capacity(2 * ambient_dim);
        for i in 0..ambient_dim {
            let mut pos = vec![0i128; ambient_dim];
            let mut neg = vec![0i128; ambient_dim];
            pos[i] = 1;
            neg[i] = -1;
            rays.push(pos);
            rays.push(neg);
        }

        Self {
            ambient_dim,
            rays: Some(rays),
            hyperplanes: Some(Vec::new()),
            dim: Some(ambient_dim),
        }
    }

    /// Get the ambient dimension.
    pub const fn ambient_dim(&self) -> usize {
        self.ambient_dim
    }

    /// Get the dimension of the cone.
    pub fn dim(&mut self) -> usize {
        if let Some(d) = self.dim {
            return d;
        }

        let d = if let Some(ref rays) = self.rays {
            matrix_rank(rays)
        } else {
            // Need to compute rays first
            let _ = self.rays();
            matrix_rank(self.rays.as_ref().unwrap())
        };

        self.dim = Some(d);
        d
    }

    /// Get the rays (V-representation).
    ///
    /// If only hyperplanes were provided, computes rays via DDM.
    pub fn rays(&mut self) -> &[Vec<i128>] {
        if self.rays.is_none() {
            // Compute from hyperplanes using DDM
            let h = self.hyperplanes.as_ref().unwrap();
            let rays = ddm::dualize(h, self.ambient_dim);
            self.rays = Some(rays);
        }
        self.rays.as_ref().unwrap()
    }

    /// Get the hyperplanes (H-representation).
    ///
    /// If only rays were provided, computes hyperplanes via DDM.
    pub fn hyperplanes(&mut self) -> &[Vec<i128>] {
        if self.hyperplanes.is_none() {
            // Compute from rays using DDM (with cache)
            let rays_ref = self.rays.as_ref().unwrap();
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            let mut key_rays = canonicalize_vectors(rays_ref);
            key_rays.sort();
            key_rays.hash(&mut hasher);
            let key = hasher.finish();
            let cache_dir = std::env::var("CYRUS_CACHE_DIR")
                .map(PathBuf::from)
                .unwrap_or_else(|_| PathBuf::from("target/cyrus-cache"));
            let cache_path = cache_dir.join(format!("hyperplanes_{key:x}.json"));

            if cache_path.exists() {
                if let Ok(data) = std::fs::read_to_string(&cache_path) {
                    if let Ok(h) = serde_json::from_str::<Vec<Vec<i128>>>(&data) {
                        eprintln!(
                            "[DEBUG] hyperplanes: cache hit (count={})",
                            h.len()
                        );
                        self.hyperplanes = Some(h);
                    }
                }
            }

            if self.hyperplanes.is_none() {
                if let Some(h) = load_compatible_hyperplanes(&cache_dir, rays_ref) {
                    eprintln!(
                        "[DEBUG] hyperplanes: reused compatible cache (count={})",
                        h.len()
                    );
                    self.hyperplanes = Some(h.clone());
                    if let Err(e) = std::fs::create_dir_all(&cache_dir) {
                        eprintln!(
                            "[WARN] failed to create hyperplane cache dir {}: {}",
                            cache_dir.display(),
                            e
                        );
                    } else if let Ok(json) = serde_json::to_string(&h) {
                        if let Err(e) = std::fs::write(&cache_path, json) {
                            eprintln!(
                                "[WARN] failed to write hyperplane cache {}: {}",
                                cache_path.display(),
                                e
                            );
                        }
                    }
                }
            }

            if self.hyperplanes.is_none() {
                eprintln!(
                    "[DEBUG] hyperplanes: reducing rays (input count={})",
                    rays_ref.len()
                );
                let rays = self.extremal_rays();
                eprintln!("[DEBUG] hyperplanes: extremal rays count={}", rays.len());
                let hyperplanes = ddm::dualize(&rays, self.ambient_dim);
                if let Err(e) = std::fs::create_dir_all(&cache_dir) {
                    eprintln!(
                        "[WARN] failed to create hyperplane cache dir {}: {}",
                        cache_dir.display(),
                        e
                    );
                } else if let Ok(json) = serde_json::to_string(&hyperplanes) {
                    if let Err(e) = std::fs::write(&cache_path, json) {
                        eprintln!(
                            "[WARN] failed to write hyperplane cache {}: {}",
                            cache_path.display(),
                            e
                        );
                    }
                }
                self.hyperplanes = Some(hyperplanes);
            }
        }
        self.hyperplanes.as_ref().unwrap()
    }

    /// Get the dual cone.
    ///
    /// The dual of a cone C is: `C* = {y : y·x >= 0 for all x in C}`
    ///
    /// If C is defined by rays R, then C* is defined by hyperplanes R.
    /// If C is defined by hyperplanes H, then C* is defined by rays H.
    pub fn dual(&self) -> Self {
        if let Some(ref rays) = self.rays {
            Self::from_hyperplanes(rays.clone())
        } else {
            Self::from_rays(self.hyperplanes.as_ref().unwrap().clone())
        }
    }

    /// Check if the cone contains a point.
    ///
    /// A point x is in the cone if H @ x >= 0 for all hyperplanes H.
    pub fn contains(&mut self, point: &[f64]) -> bool {
        self.contains_with_tolerance(point, 0.0)
    }

    /// Check if the cone contains a point with tolerance.
    ///
    /// A point x is in the cone if H @ x >= -eps for all hyperplanes H.
    pub fn contains_with_tolerance(&mut self, point: &[f64], eps: f64) -> bool {
        let h = self.hyperplanes();
        for row in h {
            let dot: f64 = row
                .iter()
                .zip(point)
                .map(|(&a, &b)| a as f64 * b)
                .sum();
            if dot < -eps {
                return false;
            }
        }
        true
    }

    /// Check if a point is in the strict interior.
    ///
    /// A point x is in the strict interior if H @ x > eps for all hyperplanes H.
    pub fn contains_interior(&mut self, point: &[f64], eps: f64) -> bool {
        let h = self.hyperplanes();
        for row in h {
            let dot: f64 = row
                .iter()
                .zip(point)
                .map(|(&a, &b)| a as f64 * b)
                .sum();
            if dot <= eps {
                return false;
            }
        }
        true
    }

    /// Check if the cone is pointed (strongly convex).
    ///
    /// A cone is pointed if it contains no lines, i.e., there is no x != 0
    /// such that both x and -x are in the cone.
    ///
    /// Equivalent to: the dual cone is solid (full-dimensional).
    pub fn is_pointed(&mut self) -> bool {
        self.dual().is_solid()
    }

    /// Check if the cone is solid (full-dimensional).
    ///
    /// A cone is solid if its dimension equals the ambient dimension.
    pub fn is_solid(&mut self) -> bool {
        self.dim() == self.ambient_dim
    }

    /// Find a point in the strict interior of the cone.
    ///
    /// If rays are known, returns the sum of all rays (normalized).
    /// Otherwise uses LP to find an interior point.
    pub fn find_interior_point(&mut self) -> Option<Vec<f64>> {
        // Simple case: if we have rays, sum them
        if let Some(ref rays) = self.rays {
            if rays.is_empty() {
                return None;
            }

            let mut sum = vec![0.0; self.ambient_dim];
            for ray in rays {
                for (i, &r) in ray.iter().enumerate() {
                    sum[i] += r as f64;
                }
            }

            // Check it's actually interior
            if self.contains_interior(&sum, 1e-10) {
                return Some(sum);
            }
        }

        // Fall back to LP
        tip::find_interior_point_lp(self)
    }

    /// Find the tip of the stretched cone.
    ///
    /// The stretched cone is the region at least distance `c` from each hyperplane.
    /// The tip is the point with minimum norm in this region.
    ///
    /// This is used to find a "canonical" interior point of the Kähler cone.
    pub fn tip_of_stretched_cone(&mut self, c: f64) -> Option<Vec<f64>> {
        tip::tip_of_stretched_cone(self, c)
    }

    /// Find lattice points in the cone using an OR-Tools CP-SAT helper.
    ///
    /// This mirrors CYTools `Cone.find_lattice_points` for the integer search.
    /// Either `min_points` or `max_deg` must be provided (but not both).
    pub fn find_lattice_points_ortools(
        &mut self,
        min_points: Option<usize>,
        max_deg: Option<i64>,
        grading_vector: &[i64],
        max_coord: i64,
        deg_window: i64,
    ) -> Result<Vec<Vec<i64>>> {
        if (min_points.is_some() && max_deg.is_some()) || (min_points.is_none() && max_deg.is_none())
        {
            return Err(Error::InvalidInput(
                "Either min_points or max_deg must be set (exclusively)".into(),
            ));
        }
        if grading_vector.is_empty() {
            return Err(Error::InvalidInput("grading vector is empty".into()));
        }

        eprintln!("[DEBUG] lattice_points: computing hyperplanes...");
        let hyperplanes = self.hyperplanes().to_vec();
        eprintln!(
            "[DEBUG] lattice_points: hyperplanes count={}, dim={}",
            hyperplanes.len(),
            hyperplanes.first().map(Vec::len).unwrap_or(0)
        );
        let mut hyperplanes_i64 = Vec::with_capacity(hyperplanes.len());
        for row in hyperplanes {
            let mut out = Vec::with_capacity(row.len());
            for v in row {
                let v64 = i64::try_from(v).map_err(|_| {
                    Error::InvalidInput("hyperplane coefficient does not fit in i64".into())
                })?;
                out.push(v64);
            }
            hyperplanes_i64.push(out);
        }

        let max_time_sec = env::var("CYRUS_LATTICE_MAX_TIME_SEC")
            .ok()
            .and_then(|v| v.parse::<f64>().ok());
        let max_solutions = env::var("CYRUS_LATTICE_MAX_SOLUTIONS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok());
        let max_deg_limit = env::var("CYRUS_LATTICE_MAX_DEG")
            .ok()
            .and_then(|v| v.parse::<i64>().ok());
        let num_search_workers = env::var("CYRUS_LATTICE_THREADS")
            .ok()
            .and_then(|v| v.parse::<i32>().ok());
        let max_coord = env::var("CYRUS_LATTICE_MAX_COORD")
            .ok()
            .and_then(|v| v.parse::<i64>().ok())
            .unwrap_or(max_coord);
        let deg_window = env::var("CYRUS_LATTICE_DEG_WINDOW")
            .ok()
            .and_then(|v| v.parse::<i64>().ok())
            .unwrap_or(deg_window);
        let strict = env::var("CYRUS_LATTICE_STRICT")
            .map(|v| v != "0")
            .unwrap_or(true);
        if max_time_sec.is_some()
            || max_solutions.is_some()
            || max_deg_limit.is_some()
            || num_search_workers.is_some()
            || max_coord != 1000
            || deg_window != 0
            || !strict
        {
            eprintln!(
                "[DEBUG] lattice_points overrides: max_time_sec={:?}, max_solutions={:?}, max_deg_limit={:?}, num_search_workers={:?}, max_coord={}, deg_window={}, strict={}",
                max_time_sec,
                max_solutions,
                max_deg_limit,
                num_search_workers,
                max_coord,
                deg_window,
                strict
            );
        }

        let req = LatticeRequest {
            hyperplanes: hyperplanes_i64,
            grading_vector: grading_vector.to_vec(),
            min_points,
            max_deg,
            max_coord,
            deg_window,
            max_deg_limit,
            max_time_sec,
            max_solutions,
            num_search_workers,
            strict,
            c: LatticeC::Scalar(0.0),
            check_grading: true,
        };

        let script = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("scripts")
            .join("find_lattice_points_ortools.py");
        if !script.exists() {
            return Err(Error::InvalidInput(format!(
                "Missing lattice-point helper script: {}",
                script.display()
            )));
        }

        let python = std::env::var("CYRUS_PYTHON").unwrap_or_else(|_| "python3".to_string());
        let mut child = Command::new(python)
            .arg(script)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| Error::InvalidInput(format!("Failed to spawn python3: {e}")))?;

        {
            let stdin = child
                .stdin
                .as_mut()
                .ok_or_else(|| Error::InvalidInput("Failed to open python stdin".into()))?;
            let payload = serde_json::to_vec(&req)
                .map_err(|e| Error::InvalidInput(format!("JSON encode failed: {e}")))?;
            stdin
                .write_all(&payload)
                .map_err(|e| Error::InvalidInput(format!("Write to python stdin failed: {e}")))?;
        }

        let output = child
            .wait_with_output()
            .map_err(|e| Error::InvalidInput(format!("Failed to run python: {e}")))?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::InvalidInput(format!(
                "Lattice-point helper failed: {stderr}"
            )));
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        let mut points: Vec<Vec<i64>> = Vec::new();
        for line in stdout.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let pt: Vec<i64> = serde_json::from_str(line).map_err(|e| {
                Error::InvalidInput(format!("Failed to parse lattice-point line: {e}"))
            })?;
            points.push(pt);
        }

        // Sort by degree, then lexicographically (CYTools behavior).
        let degs: Vec<i64> = points
            .iter()
            .map(|p| {
                p.iter()
                    .zip(grading_vector.iter())
                    .map(|(&x, &g)| x * g)
                    .sum()
            })
            .collect();
        let mut idx: Vec<usize> = (0..points.len()).collect();
        idx.sort_by(|&a, &b| {
            let da = degs[a];
            let db = degs[b];
            da.cmp(&db).then_with(|| points[a].cmp(&points[b]))
        });
        let mut sorted = Vec::with_capacity(points.len());
        for i in idx {
            sorted.push(std::mem::take(&mut points[i]));
        }
        Ok(sorted)
    }

    /// Compute extremal rays (minimal generating set).
    ///
    /// For pointed cones, returns the unique minimal set of rays.
    pub fn extremal_rays(&mut self) -> Vec<Vec<i128>> {
        let rays = self.rays().to_vec();

        // Remove duplicates
        let unique: HashSet<Vec<i128>> = rays.into_iter().collect();
        let rays: Vec<Vec<i128>> = unique.into_iter().collect();

        if rays.len() <= 1 {
            return rays;
        }

        // Filter to extremal rays using LP
        let mut extremal = Vec::new();
        for (i, ray) in rays.iter().enumerate() {
            if is_extremal(&rays, i) {
                extremal.push(ray.clone());
            }
        }

        extremal
    }

    /// Intersection of this cone with another cone.
    ///
    /// Returns a new cone defined by the union of hyperplanes.
    pub fn intersection(&mut self, other: &mut Self) -> Self {
        let mut h1 = self.hyperplanes().to_vec();
        let h2 = other.hyperplanes();
        h1.extend(h2.iter().cloned());
        Self::from_hyperplanes(h1)
    }
}

/// Normalize vectors by dividing by their GCD.
fn normalize_vectors(vectors: Vec<Vec<i128>>) -> Vec<Vec<i128>> {
    vectors
        .into_iter()
        .filter_map(|v| {
            let g = gcd_vec(&v);
            if g == 0 {
                None // Skip zero vectors
            } else {
                Some(v.iter().map(|&x| x / g).collect())
            }
        })
        .collect()
}

/// Canonicalize vectors for stable hashing: normalize by GCD and sign.
fn canonicalize_vectors(vectors: &[Vec<i128>]) -> Vec<Vec<i128>> {
    vectors
        .iter()
        .filter_map(|v| {
            let g = gcd_vec(v);
            if g == 0 {
                return None;
            }
            let mut out: Vec<i128> = v.iter().map(|&x| x / g).collect();
            for &x in &out {
                if x > 0 {
                    return Some(out);
                } else if x < 0 {
                    out = out.iter().map(|&y| -y).collect();
                    return Some(out);
                }
            }
            Some(out)
        })
        .collect()
}

fn load_compatible_hyperplanes(
    cache_dir: &PathBuf,
    rays: &[Vec<i128>],
) -> Option<Vec<Vec<i128>>> {
    let entries = std::fs::read_dir(cache_dir).ok()?;
    for entry in entries.flatten() {
        let path = entry.path();
        if !path
            .file_name()
            .and_then(|s| s.to_str())
            .map(|s| s.starts_with("hyperplanes_") && s.ends_with(".json"))
            .unwrap_or(false)
        {
            continue;
        }
        let data = std::fs::read_to_string(&path).ok()?;
        let h: Vec<Vec<i128>> = serde_json::from_str(&data).ok()?;
        if h.is_empty() {
            continue;
        }
        if h[0].len() != rays[0].len() {
            continue;
        }
        if hyperplanes_compatible(&h, rays) {
            return Some(h);
        }
    }
    None
}

fn hyperplanes_compatible(h: &[Vec<i128>], rays: &[Vec<i128>]) -> bool {
    for ray in rays {
        for row in h {
            let mut dot: i128 = 0;
            for (a, b) in row.iter().zip(ray.iter()) {
                dot += a * b;
            }
            if dot < 0 {
                return false;
            }
        }
    }
    true
}

/// GCD of a vector of integers.
fn gcd_vec(v: &[i128]) -> i128 {
    v.iter().fold(0, |acc, &x| gcd(acc, x.abs()))
}

/// GCD of two integers.
fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Compute matrix rank using Gaussian elimination.
fn matrix_rank(matrix: &[Vec<i128>]) -> usize {
    if matrix.is_empty() {
        return 0;
    }

    let rows = matrix.len();
    let cols = matrix[0].len();

    // Convert to f64 for numerical rank computation
    let mut mat: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| row.iter().map(|&x| x as f64).collect())
        .collect();

    let mut rank = 0;
    let mut pivot_col = 0;

    for row in 0..rows {
        if pivot_col >= cols {
            break;
        }

        // Find pivot
        let mut max_row = row;
        let mut max_val = mat[row][pivot_col].abs();
        for r in (row + 1)..rows {
            if mat[r][pivot_col].abs() > max_val {
                max_val = mat[r][pivot_col].abs();
                max_row = r;
            }
        }

        if max_val < 1e-10 {
            pivot_col += 1;
            continue;
        }

        mat.swap(row, max_row);
        rank += 1;

        // Eliminate
        for r in (row + 1)..rows {
            if mat[r][pivot_col].abs() > 1e-15 {
                let factor = mat[r][pivot_col] / mat[row][pivot_col];
                for c in pivot_col..cols {
                    mat[r][c] -= factor * mat[row][c];
                }
            }
        }

        pivot_col += 1;
    }

    rank
}

/// Check if ray i is extremal in the set of rays.
///
/// A ray r is extremal if it cannot be written as a non-negative
/// combination of the other rays.
fn is_extremal(rays: &[Vec<i128>], i: usize) -> bool {
    let target = &rays[i];
    let others: Vec<&Vec<i128>> = rays
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .map(|(_, r)| r)
        .collect();

    if others.is_empty() {
        return true;
    }

    let dim = target.len();
    let mut vars = ProblemVariables::new();
    let lambdas: Vec<_> = (0..others.len())
        .map(|_| vars.add(variable().min(0.0)))
        .collect();

    let mut model = vars.minimise(0.0).using(default_solver);

    for k in 0..dim {
        let mut expr = good_lp::Expression::from(0.0);
        for (j, ray) in others.iter().enumerate() {
            expr.add_mul(ray[k] as f64, lambdas[j]);
        }
        model = model.with(constraint!(expr == target[k] as f64));
    }

    model.solve().is_err()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cone_from_rays() {
        let rays = vec![vec![1, 0], vec![0, 1]];
        let mut cone = Cone::from_rays(rays);

        assert_eq!(cone.ambient_dim(), 2);
        assert_eq!(cone.rays().len(), 2);
    }

    #[test]
    fn test_cone_from_hyperplanes() {
        let hyperplanes = vec![vec![1, 0], vec![0, 1]];
        let mut cone = Cone::from_hyperplanes(hyperplanes);

        assert_eq!(cone.ambient_dim(), 2);
        assert_eq!(cone.hyperplanes().len(), 2);
    }

    #[test]
    fn test_cone_contains() {
        // First quadrant: x >= 0, y >= 0
        let mut cone = Cone::from_hyperplanes(vec![vec![1, 0], vec![0, 1]]);

        assert!(cone.contains(&[1.0, 1.0]));
        assert!(cone.contains(&[0.0, 0.0]));
        assert!(!cone.contains(&[-1.0, 1.0]));
        assert!(!cone.contains(&[1.0, -1.0]));
    }

    #[test]
    fn test_cone_dual() {
        let rays = vec![vec![1, 0], vec![1, 1]];
        let cone = Cone::from_rays(rays);
        let dual = cone.dual();

        // Dual is defined by hyperplanes = original rays
        assert!(dual.hyperplanes.is_some());
    }

    #[test]
    fn test_normalize_vectors() {
        let v = vec![vec![2, 4, 6], vec![0, 0, 0], vec![3, 6, 9]];
        let normalized = normalize_vectors(v);

        // Should have 2 vectors (zero vector removed)
        assert_eq!(normalized.len(), 2);
        // [2,4,6] / 2 = [1,2,3]
        assert!(normalized.contains(&vec![1, 2, 3]));
        // [3,6,9] / 3 = [1,2,3]
        assert!(normalized.contains(&vec![1, 2, 3]));
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(17, 13), 1);
        assert_eq!(gcd(0, 5), 5);
        assert_eq!(gcd(5, 0), 5);
    }

    #[test]
    fn test_matrix_rank() {
        // Full rank 2x2
        let m1 = vec![vec![1, 0], vec![0, 1]];
        assert_eq!(matrix_rank(&m1), 2);

        // Rank 1 (dependent rows)
        let m2 = vec![vec![1, 2], vec![2, 4]];
        assert_eq!(matrix_rank(&m2), 1);

        // Rank 2 for 3x2
        let m3 = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        assert_eq!(matrix_rank(&m3), 2);
    }
}
