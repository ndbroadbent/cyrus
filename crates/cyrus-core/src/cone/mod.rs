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

use good_lp::{ProblemVariables, Solution, SolverModel, constraint, default_solver, variable};
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

const HYPERPLANE_CACHE_VERSION: &str = "hyperplanes-exact-v3";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum ConeDualBackend {
    Ddm,
    PplLcdd,
}

impl ConeDualBackend {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Ddm => "ddm",
            Self::PplLcdd => "ppl_lcdd",
        }
    }
}

fn cone_dual_backend_from_env() -> ConeDualBackend {
    match env::var("CYRUS_CONE_DUAL_BACKEND").as_deref() {
        Ok("ppl_lcdd") => ConeDualBackend::PplLcdd,
        Ok("ddm") | Err(_) => ConeDualBackend::Ddm,
        Ok(other) => {
            panic!("unsupported CYRUS_CONE_DUAL_BACKEND={other}; expected ddm or ppl_lcdd")
        }
    }
}

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
            let rays_ref = self.rays.as_ref().unwrap().clone();
            let (cache_dir, cache_path) = hyperplane_cache_paths(&rays_ref);
            if let Some(h) = load_hyperplanes_cache(&cache_path) {
                eprintln!("[DEBUG] hyperplanes: cache hit (count={})", h.len());
                self.hyperplanes = Some(h);
            }

            if self.hyperplanes.is_none() {
                let hyperplanes = compute_hyperplanes_from_rays(
                    &rays_ref,
                    self.ambient_dim,
                    &cache_dir,
                    &cache_path,
                );
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
            let dot: f64 = row.iter().zip(point).map(|(&a, &b)| a as f64 * b).sum();
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
            let dot: f64 = row.iter().zip(point).map(|(&a, &b)| a as f64 * b).sum();
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
        validate_lattice_args(min_points, max_deg, grading_vector)?;

        eprintln!("[DEBUG] lattice_points: computing hyperplanes...");
        let hyperplanes = self.hyperplanes().to_vec();
        eprintln!(
            "[DEBUG] lattice_points: hyperplanes count={}, dim={}",
            hyperplanes.len(),
            hyperplanes.first().map_or(0, Vec::len)
        );
        let hyperplanes_i64 = hyperplanes_to_i64(hyperplanes)?;

        let overrides = LatticeOverrides::from_env(max_coord, deg_window);
        overrides.log_if_overridden();

        let req = LatticeRequest {
            hyperplanes: hyperplanes_i64,
            grading_vector: grading_vector.to_vec(),
            min_points,
            max_deg,
            max_coord: overrides.max_coord,
            deg_window: overrides.deg_window,
            max_deg_limit: overrides.max_deg_limit,
            max_time_sec: overrides.max_time_sec,
            max_solutions: overrides.max_solutions,
            num_search_workers: overrides.num_search_workers,
            strict: overrides.strict,
            c: LatticeC::Scalar(0.0),
            check_grading: true,
        };

        if let Some(points) = find_lattice_points_native(&req)? {
            eprintln!(
                "[DEBUG] lattice_points: native bounded enumerator returned {} points",
                points.len()
            );
            return Ok(sort_lattice_points(points, grading_vector));
        }

        let script = lattice_helper_script()?;
        let stdout = run_lattice_helper(&req, &script)?;
        let points = parse_lattice_points(&stdout)?;
        Ok(sort_lattice_points(points, grading_vector))
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

fn validate_lattice_args(
    min_points: Option<usize>,
    max_deg: Option<i64>,
    grading_vector: &[i64],
) -> Result<()> {
    if (min_points.is_some() && max_deg.is_some()) || (min_points.is_none() && max_deg.is_none()) {
        return Err(Error::InvalidInput(
            "Either min_points or max_deg must be set (exclusively)".into(),
        ));
    }
    if grading_vector.is_empty() {
        return Err(Error::InvalidInput("grading vector is empty".into()));
    }
    Ok(())
}

fn hyperplanes_to_i64(hyperplanes: Vec<Vec<i128>>) -> Result<Vec<Vec<i64>>> {
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
    Ok(hyperplanes_i64)
}

struct LatticeOverrides {
    max_time_sec: Option<f64>,
    max_solutions: Option<usize>,
    max_deg_limit: Option<i64>,
    num_search_workers: Option<i32>,
    max_coord: i64,
    deg_window: i64,
    strict: bool,
    has_overrides: bool,
}

impl LatticeOverrides {
    fn from_env(default_max_coord: i64, default_deg_window: i64) -> Self {
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
            .unwrap_or(default_max_coord);
        let deg_window = env::var("CYRUS_LATTICE_DEG_WINDOW")
            .ok()
            .and_then(|v| v.parse::<i64>().ok())
            .unwrap_or(default_deg_window);
        let strict = env::var("CYRUS_LATTICE_STRICT")
            .map(|v| v != "0")
            .unwrap_or(true);
        let has_overrides = max_time_sec.is_some()
            || max_solutions.is_some()
            || max_deg_limit.is_some()
            || num_search_workers.is_some()
            || max_coord != default_max_coord
            || deg_window != default_deg_window
            || !strict;
        Self {
            max_time_sec,
            max_solutions,
            max_deg_limit,
            num_search_workers,
            max_coord,
            deg_window,
            strict,
            has_overrides,
        }
    }

    fn log_if_overridden(&self) {
        if !self.has_overrides {
            return;
        }
        eprintln!(
            "[DEBUG] lattice_points overrides: max_time_sec={:?}, max_solutions={:?}, max_deg_limit={:?}, num_search_workers={:?}, max_coord={}, deg_window={}, strict={}",
            self.max_time_sec,
            self.max_solutions,
            self.max_deg_limit,
            self.num_search_workers,
            self.max_coord,
            self.deg_window,
            self.strict
        );
    }
}

fn lattice_helper_script() -> Result<PathBuf> {
    let script = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("scripts")
        .join("find_lattice_points_ortools.py");
    if !script.exists() {
        return Err(Error::InvalidInput(format!(
            "Missing lattice-point helper script: {}",
            script.display()
        )));
    }
    Ok(script)
}

fn run_lattice_helper(req: &LatticeRequest, script: &PathBuf) -> Result<String> {
    let python = env::var("CYRUS_PYTHON").unwrap_or_else(|_| "python3".to_string());
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
        let payload = serde_json::to_vec(req)
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
    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

fn parse_lattice_points(stdout: &str) -> Result<Vec<Vec<i64>>> {
    let mut points: Vec<Vec<i64>> = Vec::new();
    for line in stdout.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let pt: Vec<i64> = serde_json::from_str(line)
            .map_err(|e| Error::InvalidInput(format!("Failed to parse lattice-point line: {e}")))?;
        points.push(pt);
    }
    Ok(points)
}

fn find_lattice_points_native(req: &LatticeRequest) -> Result<Option<Vec<Vec<i64>>>> {
    if env::var("CYRUS_LATTICE_NATIVE")
        .map(|value| value == "0")
        .unwrap_or(false)
    {
        return Ok(None);
    }
    if req.max_time_sec.is_some() || req.num_search_workers.is_some() {
        return Ok(None);
    }
    if req.hyperplanes.is_empty() || req.grading_vector.is_empty() {
        return Ok(None);
    }
    if !matches!(req.c, LatticeC::Scalar(c) if c == 0.0) {
        return Ok(None);
    }

    let dim = req.grading_vector.len();
    if dim > 6 || req.hyperplanes.iter().any(|row| row.len() != dim) {
        return Ok(None);
    }
    if req.deg_window < 0 || req.max_coord < 0 {
        return Err(Error::InvalidInput(
            "native lattice search requires non-negative max_coord and deg_window".into(),
        ));
    }

    if req.check_grading {
        let check = enumerate_lattice_degree_window(req, None, Some(0), Some(2))?;
        if check.len() > 1 {
            return Err(Error::InvalidInput(
                "grading vector must be positive on the cone".into(),
            ));
        }
    }

    let mut points = Vec::new();
    if let Some(max_deg) = req.max_deg {
        points = enumerate_lattice_degree_window(req, None, Some(max_deg), req.max_solutions)?;
        if let Some(limit) = req.max_solutions
            && points.len() >= limit
        {
            return Err(Error::InvalidInput(
                "max_solutions reached before completion".into(),
            ));
        }
    } else if let Some(min_points) = req.min_points {
        let mut degree = 0i64;
        while points.len() < min_points {
            if let Some(max_deg_limit) = req.max_deg_limit
                && degree > max_deg_limit
            {
                return Err(Error::InvalidInput(
                    "max_deg_limit reached before min_points".into(),
                ));
            }
            let remaining_limit = req
                .max_solutions
                .map(|limit| limit.saturating_sub(points.len()));
            let slice = enumerate_lattice_degree_window(
                req,
                Some(degree),
                Some(degree + req.deg_window),
                remaining_limit,
            )?;
            points.extend(slice);
            if let Some(limit) = req.max_solutions
                && points.len() >= limit
                && points.len() < min_points
            {
                return Err(Error::InvalidInput(
                    "max_solutions reached before min_points".into(),
                ));
            }
            degree += req.deg_window + 1;
        }
    } else {
        return Ok(None);
    }

    Ok(Some(points))
}

fn enumerate_lattice_degree_window(
    req: &LatticeRequest,
    deg_low: Option<i64>,
    deg_high: Option<i64>,
    max_solutions: Option<usize>,
) -> Result<Vec<Vec<i64>>> {
    let bounds = lattice_coordinate_bounds(req, deg_low, deg_high)?;
    if bounds.iter().any(|(lo, hi)| lo > hi) {
        return Ok(Vec::new());
    }

    let mut points = Vec::new();
    let mut candidate = vec![0i64; req.grading_vector.len()];
    enumerate_lattice_recursive(
        req,
        &bounds,
        deg_low,
        deg_high,
        max_solutions,
        0,
        &mut candidate,
        &mut points,
    );
    Ok(points)
}

#[allow(clippy::too_many_arguments)]
fn enumerate_lattice_recursive(
    req: &LatticeRequest,
    bounds: &[(i64, i64)],
    deg_low: Option<i64>,
    deg_high: Option<i64>,
    max_solutions: Option<usize>,
    coord: usize,
    candidate: &mut [i64],
    points: &mut Vec<Vec<i64>>,
) {
    if max_solutions.is_some_and(|limit| points.len() >= limit) {
        return;
    }
    if coord == candidate.len() {
        if lattice_point_satisfies(req, candidate, deg_low, deg_high) {
            points.push(candidate.to_vec());
        }
        return;
    }

    let (lo, hi) = bounds[coord];
    for value in lo..=hi {
        candidate[coord] = value;
        enumerate_lattice_recursive(
            req,
            bounds,
            deg_low,
            deg_high,
            max_solutions,
            coord + 1,
            candidate,
            points,
        );
        if max_solutions.is_some_and(|limit| points.len() >= limit) {
            return;
        }
    }
}

fn lattice_coordinate_bounds(
    req: &LatticeRequest,
    deg_low: Option<i64>,
    deg_high: Option<i64>,
) -> Result<Vec<(i64, i64)>> {
    let dim = req.grading_vector.len();
    let mut bounds = Vec::with_capacity(dim);
    for coord in 0..dim {
        let min = optimize_lattice_coordinate(req, deg_low, deg_high, coord, false)?;
        let max = optimize_lattice_coordinate(req, deg_low, deg_high, coord, true)?;
        let lo = (min - 1e-9).ceil() as i64;
        let hi = (max + 1e-9).floor() as i64;
        bounds.push((lo.max(-req.max_coord), hi.min(req.max_coord)));
    }
    Ok(bounds)
}

fn optimize_lattice_coordinate(
    req: &LatticeRequest,
    deg_low: Option<i64>,
    deg_high: Option<i64>,
    coord: usize,
    maximize: bool,
) -> Result<f64> {
    let dim = req.grading_vector.len();
    let mut vars = ProblemVariables::new();
    let x: Vec<_> = (0..dim)
        .map(|_| {
            vars.add(
                variable()
                    .min(-(req.max_coord as f64))
                    .max(req.max_coord as f64),
            )
        })
        .collect();

    let mut objective = good_lp::Expression::from(0.0);
    objective.add_mul(if maximize { -1.0 } else { 1.0 }, x[coord]);
    let mut model = vars.minimise(objective).using(default_solver);

    for row in &req.hyperplanes {
        let mut expr = good_lp::Expression::from(0.0);
        for (coeff, var) in row.iter().zip(x.iter()) {
            expr.add_mul(*coeff as f64, *var);
        }
        model = model.with(constraint!(expr >= 0.0));
    }

    let mut degree_expr = good_lp::Expression::from(0.0);
    for (coeff, var) in req.grading_vector.iter().zip(x.iter()) {
        degree_expr.add_mul(*coeff as f64, *var);
    }
    if let Some(low) = deg_low {
        model = model.with(constraint!(degree_expr.clone() >= low as f64));
    }
    if let Some(high) = deg_high {
        model = model.with(constraint!(degree_expr <= high as f64));
    }

    let solution = model
        .solve()
        .map_err(|e| Error::InvalidInput(format!("native lattice LP bound failed: {e}")))?;
    Ok(solution.value(x[coord]))
}

fn lattice_point_satisfies(
    req: &LatticeRequest,
    candidate: &[i64],
    deg_low: Option<i64>,
    deg_high: Option<i64>,
) -> bool {
    if req.hyperplanes.iter().any(|row| {
        row.iter()
            .zip(candidate.iter())
            .map(|(&a, &x)| i128::from(a) * i128::from(x))
            .sum::<i128>()
            < 0
    }) {
        return false;
    }

    let degree = req
        .grading_vector
        .iter()
        .zip(candidate.iter())
        .map(|(&a, &x)| i128::from(a) * i128::from(x))
        .sum::<i128>();
    if let Some(low) = deg_low
        && degree < i128::from(low)
    {
        return false;
    }
    if let Some(high) = deg_high
        && degree > i128::from(high)
    {
        return false;
    }
    true
}

fn sort_lattice_points(mut points: Vec<Vec<i64>>, grading_vector: &[i64]) -> Vec<Vec<i64>> {
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
    sorted
}

/// Normalize vectors by dividing by their GCD.
fn normalize_vectors(vectors: Vec<Vec<i128>>) -> Vec<Vec<i128>> {
    let mut seen = HashSet::new();
    let mut normalized = Vec::new();
    for v in vectors {
        let g = gcd_vec(&v);
        if g == 0 {
            continue;
        }
        let row = v.iter().map(|&x| x / g).collect::<Vec<_>>();
        if seen.insert(row.clone()) {
            normalized.push(row);
        }
    }
    normalized
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

fn hyperplane_cache_paths(rays: &[Vec<i128>]) -> (PathBuf, PathBuf) {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    HYPERPLANE_CACHE_VERSION.hash(&mut hasher);
    cone_dual_backend_from_env().as_str().hash(&mut hasher);
    let mut key_rays = canonicalize_vectors(rays);
    key_rays.sort();
    let key_len = key_rays.len();
    key_rays.hash(&mut hasher);
    let key = hasher.finish() ^ (key_len as u64);
    let cache_dir = env::var("CYRUS_CACHE_DIR")
        .map_or_else(|_| PathBuf::from("target/cyrus-cache"), PathBuf::from);
    let cache_path = cache_dir.join(format!("hyperplanes_{key:x}.json"));
    (cache_dir, cache_path)
}

fn load_hyperplanes_cache(cache_path: &PathBuf) -> Option<Vec<Vec<i128>>> {
    if !cache_path.exists() {
        return None;
    }
    let data = std::fs::read_to_string(cache_path).ok()?;
    serde_json::from_str::<Vec<Vec<i128>>>(&data).ok()
}

fn write_hyperplanes_cache(cache_dir: &PathBuf, cache_path: &PathBuf, h: &[Vec<i128>]) {
    if let Err(e) = std::fs::create_dir_all(cache_dir) {
        eprintln!(
            "[WARN] failed to create hyperplane cache dir {}: {}",
            cache_dir.display(),
            e
        );
        return;
    }
    if let Ok(json) = serde_json::to_string(h)
        && let Err(e) = std::fs::write(cache_path, json)
    {
        eprintln!(
            "[WARN] failed to write hyperplane cache {}: {}",
            cache_path.display(),
            e
        );
    }
}

fn compute_hyperplanes_from_rays(
    rays_ref: &[Vec<i128>],
    ambient_dim: usize,
    cache_dir: &PathBuf,
    cache_path: &PathBuf,
) -> Vec<Vec<i128>> {
    eprintln!(
        "[DEBUG] hyperplanes: dualizing rays (input count={})",
        rays_ref.len()
    );
    let hyperplanes = match cone_dual_backend_from_env() {
        ConeDualBackend::PplLcdd => compute_hyperplanes_from_rays_ppl_lcdd(rays_ref, ambient_dim)
            .unwrap_or_else(|err| panic!("ppl_lcdd cone dualization failed: {err}")),
        ConeDualBackend::Ddm => ddm::dualize(rays_ref, ambient_dim),
    };
    write_hyperplanes_cache(cache_dir, cache_path, &hyperplanes);
    hyperplanes
}

fn compute_hyperplanes_from_rays_ppl_lcdd(
    rays_ref: &[Vec<i128>],
    ambient_dim: usize,
) -> Result<Vec<Vec<i128>>> {
    for row in rays_ref {
        if row.len() != ambient_dim {
            return Err(Error::InvalidInput(
                "ppl_lcdd cone dualization requires consistent ray dimensions".into(),
            ));
        }
    }

    let binary = env::var("CYRUS_PPL_LCDD").unwrap_or_else(|_| "ppl_lcdd".to_string());
    let mut command = Command::new(binary);
    if let Ok(seconds) = env::var("CYRUS_PPL_LCDD_MAX_CPU_SEC") {
        command.arg(format!("--max-cpu={seconds}"));
    }
    if let Ok(megabytes) = env::var("CYRUS_PPL_LCDD_MAX_MEMORY_MB") {
        command.arg(format!("--max-memory={megabytes}"));
    }
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|err| Error::InvalidInput(format!("failed to spawn ppl_lcdd: {err}")))?;

    {
        let stdin = child
            .stdin
            .as_mut()
            .ok_or_else(|| Error::InvalidInput("failed to open ppl_lcdd stdin".into()))?;
        stdin
            .write_all(cdd_v_representation(rays_ref, ambient_dim).as_bytes())
            .map_err(|err| Error::InvalidInput(format!("failed to write ppl_lcdd input: {err}")))?;
    }

    let output = child
        .wait_with_output()
        .map_err(|err| Error::InvalidInput(format!("failed to wait for ppl_lcdd: {err}")))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(Error::InvalidInput(format!(
            "ppl_lcdd exited with status {}: {stderr}",
            output.status
        )));
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    parse_ppl_lcdd_h_representation(&stdout, ambient_dim)
}

fn cdd_v_representation(rays_ref: &[Vec<i128>], ambient_dim: usize) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    out.push_str("V-representation\nbegin\n");
    writeln!(out, "{} {} integer", rays_ref.len(), ambient_dim + 1)
        .expect("writing to a String cannot fail");
    for ray in rays_ref {
        out.push('0');
        for value in ray {
            out.push(' ');
            out.push_str(&value.to_string());
        }
        out.push('\n');
    }
    out.push_str("end\n");
    out
}

fn parse_ppl_lcdd_h_representation(output: &str, ambient_dim: usize) -> Result<Vec<Vec<i128>>> {
    let mut linearity_rows = HashSet::new();
    let mut lines = output
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty());
    let Some(header) = lines.next() else {
        return Err(Error::InvalidInput("ppl_lcdd output is empty".into()));
    };
    if header != "H-representation" {
        return Err(Error::InvalidInput(format!(
            "ppl_lcdd returned {header:?}, expected H-representation"
        )));
    }

    for line in lines.by_ref() {
        if line == "begin" {
            break;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.first() == Some(&"linearity") {
            if fields.len() < 2 {
                return Err(Error::InvalidInput(
                    "ppl_lcdd linearity header is malformed".into(),
                ));
            }
            let count = fields[1].parse::<usize>().map_err(|err| {
                Error::InvalidInput(format!("ppl_lcdd linearity count parse failed: {err}"))
            })?;
            if fields.len() != count + 2 {
                return Err(Error::InvalidInput(
                    "ppl_lcdd linearity row count mismatch".into(),
                ));
            }
            for field in &fields[2..] {
                let row = field.parse::<usize>().map_err(|err| {
                    Error::InvalidInput(format!("ppl_lcdd linearity row parse failed: {err}"))
                })?;
                if row == 0 {
                    return Err(Error::InvalidInput(
                        "ppl_lcdd linearity rows are one-based".into(),
                    ));
                }
                linearity_rows.insert(row);
            }
        }
    }

    let Some(shape_line) = lines.next() else {
        return Err(Error::InvalidInput(
            "ppl_lcdd output is missing shape line".into(),
        ));
    };
    let shape = shape_line.split_whitespace().collect::<Vec<_>>();
    if shape.len() != 3 || shape[2] != "integer" {
        return Err(Error::InvalidInput(
            "ppl_lcdd H-representation shape line is malformed".into(),
        ));
    }
    let row_count = shape[0]
        .parse::<usize>()
        .map_err(|err| Error::InvalidInput(format!("ppl_lcdd row count parse failed: {err}")))?;
    let col_count = shape[1]
        .parse::<usize>()
        .map_err(|err| Error::InvalidInput(format!("ppl_lcdd column count parse failed: {err}")))?;
    if col_count != ambient_dim + 1 {
        return Err(Error::InvalidInput(format!(
            "ppl_lcdd H-representation has {col_count} columns, expected {}",
            ambient_dim + 1
        )));
    }

    let mut hyperplanes = Vec::new();
    for row_idx in 1..=row_count {
        let Some(row_line) = lines.next() else {
            return Err(Error::InvalidInput(
                "ppl_lcdd output ended before all rows were read".into(),
            ));
        };
        let values = row_line
            .split_whitespace()
            .map(|field| {
                field.parse::<i128>().map_err(|err| {
                    Error::InvalidInput(format!("ppl_lcdd row entry parse failed: {err}"))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        if values.len() != col_count {
            return Err(Error::InvalidInput(
                "ppl_lcdd row length does not match shape line".into(),
            ));
        }
        if values[0] != 0 {
            return Err(Error::InvalidInput(
                "ppl_lcdd returned affine H-representation for a cone".into(),
            ));
        }
        let row = values[1..].to_vec();
        if linearity_rows.contains(&row_idx) {
            hyperplanes.push(row.iter().map(|&value| -value).collect::<Vec<_>>());
        }
        hyperplanes.push(row);
    }

    match lines.next() {
        Some("end") => Ok(normalize_vectors(hyperplanes)),
        Some(other) => Err(Error::InvalidInput(format!(
            "ppl_lcdd output has unexpected row after H-representation: {other}"
        ))),
        None => Err(Error::InvalidInput("ppl_lcdd output is missing end".into())),
    }
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

        assert_eq!(normalized, vec![vec![1, 2, 3]]);
    }

    #[test]
    fn cone_construction_deduplicates_normalized_rows() {
        let mut cone = Cone::from_rays(vec![vec![2, 0], vec![1, 0], vec![0, 3], vec![0, 1]]);
        assert_eq!(cone.rays(), &[vec![1, 0], vec![0, 1]]);

        let mut cone = Cone::from_hyperplanes(vec![vec![2, 0], vec![1, 0], vec![0, 3], vec![0, 1]]);
        assert_eq!(cone.hyperplanes(), &[vec![1, 0], vec![0, 1]]);
    }

    #[test]
    fn cdd_v_representation_writes_cone_rays() {
        let cdd = cdd_v_representation(&[vec![1, 0], vec![0, 1]], 2);

        assert_eq!(
            cdd,
            "V-representation\nbegin\n2 3 integer\n0 1 0\n0 0 1\nend\n"
        );
    }

    #[test]
    fn ppl_lcdd_h_representation_parser_reads_cone_hyperplanes() {
        let output = "H-representation\nbegin\n2 3 integer\n0 1 0\n0 0 1\nend\n";

        let hyperplanes = parse_ppl_lcdd_h_representation(output, 2).unwrap();

        assert_eq!(hyperplanes, vec![vec![1, 0], vec![0, 1]]);
    }

    #[test]
    fn ppl_lcdd_h_representation_parser_expands_lineality_rows() {
        let output = "H-representation\nlinearity 1 1\nbegin\n2 3 integer\n0 1 -1\n0 0 1\nend\n";

        let hyperplanes = parse_ppl_lcdd_h_representation(output, 2).unwrap();

        assert_eq!(hyperplanes, vec![vec![-1, 1], vec![1, -1], vec![0, 1]]);
    }

    #[test]
    fn ppl_lcdd_h_representation_parser_rejects_affine_rows() {
        let output = "H-representation\nbegin\n1 3 integer\n1 1 0\nend\n";

        let err = parse_ppl_lcdd_h_representation(output, 2)
            .expect_err("affine rows are not valid cone hyperplanes");

        assert!(err.to_string().contains("affine H-representation"));
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

    #[test]
    fn test_native_lattice_points_first_quadrant_min_points() {
        let req = LatticeRequest {
            hyperplanes: vec![vec![1, 0], vec![0, 1]],
            grading_vector: vec![1, 1],
            min_points: Some(4),
            max_deg: None,
            max_coord: 10,
            deg_window: 0,
            max_deg_limit: Some(3),
            max_time_sec: None,
            max_solutions: None,
            num_search_workers: None,
            strict: true,
            c: LatticeC::Scalar(0.0),
            check_grading: true,
        };

        let points = find_lattice_points_native(&req)
            .expect("native enumeration should succeed")
            .expect("request is native-supported");

        assert_eq!(
            sort_lattice_points(points, &req.grading_vector),
            vec![
                vec![0, 0],
                vec![0, 1],
                vec![1, 0],
                vec![0, 2],
                vec![1, 1],
                vec![2, 0],
            ]
        );
    }

    #[test]
    fn test_native_lattice_points_accepts_negative_coordinates() {
        let req = LatticeRequest {
            // Cone generated by (-1, 1) and (1, 0): y >= 0, x + y >= 0.
            hyperplanes: vec![vec![0, 1], vec![1, 1]],
            grading_vector: vec![1, 2],
            min_points: None,
            max_deg: Some(2),
            max_coord: 10,
            deg_window: 0,
            max_deg_limit: None,
            max_time_sec: None,
            max_solutions: None,
            num_search_workers: None,
            strict: true,
            c: LatticeC::Scalar(0.0),
            check_grading: true,
        };

        let points = find_lattice_points_native(&req)
            .expect("native enumeration should succeed")
            .expect("request is native-supported");

        assert_eq!(
            sort_lattice_points(points, &req.grading_vector),
            vec![
                vec![0, 0],
                vec![-1, 1],
                vec![1, 0],
                vec![-2, 2],
                vec![0, 1],
                vec![2, 0],
            ]
        );
    }
}
