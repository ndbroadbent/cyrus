//! Gopakumar-Vafa invariant computation utilities.
//!
//! Implements the CYTools mori_cone_cap algorithm and grading vector selection,
//! then wires up the `cygv` crate for actual GV computation.
//!
//! Reference: CYTools `calabiyau.py` (mori_cone_cap) and `cone.py` (grading vector).

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::io::BufWriter;
use std::hash::{Hash, Hasher};
use std::path::PathBuf;

use malachite::Integer;
use malachite::num::arithmetic::traits::Abs;
use nalgebra::{DMatrix, DVector, RowDVector};

use crate::cone::Cone;
use crate::error::{Error, Result};
use crate::integer_math::{gcd_integer, integer_kernel};
use crate::intersection::Intersection;
use crate::lattice::Point;
use crate::polytope::Polytope;
use crate::triangulation::Triangulation;

/// Compute the Mori cone cap generators (rays) using the CYTools algorithm.
///
/// Returns a matrix where each row is a generator (ray) expressed in the
/// divisor/curve coordinates corresponding to `points`.
///
/// This follows `CalabiYau.mori_cone_cap` in CYTools.
pub fn compute_mori_cone_cap_rays(
    tri: &Triangulation,
    points: &[Point],
    polytope: &Polytope,
    in_basis: bool,
    exclude_origin: bool,
    basis: Option<&[usize]>,
) -> Result<Vec<Vec<i64>>> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "mori_cone_cap is only implemented for 4D polytopes".into(),
        ));
    }

    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .ok_or_else(|| Error::InvalidInput("Origin not found in points".into()))?;

    // pts_ext: append 1 to each point (for affine embedding)
    let pts_ext: Vec<Vec<i64>> = points
        .iter()
        .map(|p| {
            let mut v = p.coords().to_vec();
            v.push(1);
            v
        })
        .collect();

    let (facets, twofaces) = compute_faces_4d(points, polytope)?;

    // Collect 2D simplices and circuits
    let mut mori_cap_rays: HashSet<Vec<(usize, i64)>> = HashSet::new();
    let mut simp_2d_all: HashSet<Vec<usize>> = HashSet::new();

    for f in &twofaces {
        if f.len() < 4 {
            let mut f_sorted = f.clone();
            f_sorted.sort_unstable();
            simp_2d_all.insert(f_sorted);
            continue;
        }

        let face_pts: HashSet<usize> = f.iter().copied().collect();
        let mut simp_2d: HashSet<Vec<usize>> = HashSet::new();

        for simplex in tri.simplices() {
            let inter: Vec<usize> = simplex
                .iter()
                .filter(|idx| face_pts.contains(idx))
                .copied()
                .collect();
            if inter.len() == 3 {
                let mut inter_sorted = inter;
                inter_sorted.sort_unstable();
                simp_2d.insert(inter_sorted.clone());
                simp_2d_all.insert(inter_sorted);
            }
        }

        let simps: Vec<Vec<usize>> = simp_2d.into_iter().collect();
        for i in 0..simps.len() {
            for j in i..simps.len() {
                let s1: HashSet<usize> = simps[i].iter().copied().collect();
                let s2: HashSet<usize> = simps[j].iter().copied().collect();
                let comm: Vec<usize> = s1.intersection(&s2).copied().collect();
                if comm.len() == 2 {
                    let diff: Vec<usize> = s1.symmetric_difference(&s2).copied().collect();
                    if diff.len() != 2 {
                        continue;
                    }

                    if let Some(v) = nullspace_vector(&pts_ext, &diff, &comm, false)? {
                        let full_v = build_full_v(&diff, &comm, &v);
                        mori_cap_rays.insert(full_v);
                    }
                }
            }
        }
    }

    // Origin circuits
    for s2d in &simp_2d_all {
        let s2d_set: HashSet<usize> = s2d.iter().copied().collect();
        let mut f1: Option<Vec<usize>> = None;
        let mut f2: Option<Vec<usize>> = None;
        for facet in &facets {
            let facet_set: HashSet<usize> = facet.iter().copied().collect();
            if s2d_set.is_subset(&facet_set) {
                if f1.is_none() {
                    f1 = Some(facet.clone());
                } else {
                    f2 = Some(facet.clone());
                    break;
                }
            }
        }
        let (Some(f1), Some(f2)) = (f1, f2) else {
            continue;
        };

        let f1_set: HashSet<usize> = f1.iter().copied().collect();
        let f2_set: HashSet<usize> = f2.iter().copied().collect();
        let pts_f1: Vec<usize> = f1_set.difference(&f2_set).copied().collect();
        let pts_f2: Vec<usize> = f2_set.difference(&f1_set).copied().collect();

        for p1 in &pts_f1 {
            for p2 in &pts_f2 {
                let diff = vec![*p1, *p2];
                let mut comm = s2d.clone();
                comm.push(origin_idx);
                let Some(v) = nullspace_vector(&pts_ext, &diff, &comm, true)? else {
                    continue;
                };
                let full_v = build_full_v(&diff, &comm, &v);

                // Skip if origin coefficient is zero or positive (CYTools behavior).
                let origin_coeff = full_v
                    .iter()
                    .find(|(idx, _)| *idx == origin_idx)
                    .map(|(_, coeff)| *coeff)
                    .unwrap_or(0);
                if origin_coeff == 0 {
                    continue;
                }
                if origin_coeff > 0 {
                    continue;
                }
                mori_cap_rays.insert(full_v);
            }
        }
    }

    let mut rays: Vec<Vec<i64>> = Vec::new();
    let n_pts = pts_ext.len();
    for r in mori_cap_rays {
        let mut row = vec![0i64; n_pts];
        for (idx, coeff) in r {
            row[idx] = coeff;
        }
        rays.push(row);
    }

    if exclude_origin && !in_basis {
        rays = rays.into_iter().map(|r| r[1..].to_vec()).collect();
    } else if in_basis {
        let basis = basis.ok_or_else(|| {
            Error::InvalidInput("basis indices required for in_basis=true".into())
        })?;
        rays = rays
            .into_iter()
            .map(|r| basis.iter().map(|&idx| r[idx]).collect())
            .collect();
    }

    // Normalize by GCD and drop zero rows (CYTools Cone normalization).
    let mut normalized: Vec<Vec<i64>> = Vec::with_capacity(rays.len());
    for row in rays {
        let mut g = 0i64;
        for &x in &row {
            g = gcd_i64(g, x.abs());
        }
        if g == 0 {
            continue;
        }
        normalized.push(row.into_iter().map(|x| x / g).collect());
    }

    // Deduplicate rows after normalization.
    let mut uniq: HashSet<Vec<i64>> = HashSet::new();
    let mut deduped = Vec::with_capacity(normalized.len());
    for row in normalized {
        if uniq.insert(row.clone()) {
            deduped.push(row);
        }
    }

    deduped.sort();
    Ok(deduped)
}

/// Compute a grading vector for the Mori cone cap.
///
/// The grading vector must be strictly interior to the dual cone.
pub fn compute_grading_vector(rays: &[Vec<i64>]) -> Option<Vec<i64>> {
    if rays.is_empty() {
        return None;
    }

    let (zero_rays, opposite_pairs) = analyze_rays(rays);

    // Grading vector is strictly interior to the dual cone.
    // Avoid DDM by working directly with the dual hyperplanes (the rays).
    let rays_i128: Vec<Vec<i128>> = rays
        .iter()
        .map(|row| row.iter().map(|&x| x as i128).collect())
        .collect();
    let mut dual = Cone::from_hyperplanes(rays_i128);
    let interior = dual.find_interior_point().or_else(|| {
        eprintln!(
            "[WARN] grading vector search failed: rays={}, zero_rays={}, opposite_pairs={}",
            rays.len(),
            zero_rays,
            opposite_pairs
        );
        None
    })?;

    // Scale and round until strictly interior in the dual
    let scales = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0];
    for scale in scales {
        let candidate: Vec<i64> = interior
            .iter()
            .map(|x| (x * scale).round() as i64)
            .collect();
        if candidate.iter().all(|&x| x == 0) {
            continue;
        }
        if is_strictly_dual(rays, &candidate) {
            return Some(candidate);
        }
    }

    None
}

fn analyze_rays(rays: &[Vec<i64>]) -> (usize, usize) {
    let mut zero_rays = 0usize;
    let mut opposite_pairs = 0usize;
    let mut seen: HashSet<Vec<i64>> = HashSet::new();

    for r in rays {
        if r.iter().all(|&x| x == 0) {
            zero_rays += 1;
            continue;
        }
        let norm = normalize_ray(r);
        let neg: Vec<i64> = norm.iter().map(|&x| -x).collect();
        if seen.contains(&neg) {
            opposite_pairs += 1;
        }
        seen.insert(norm);
    }

    (zero_rays, opposite_pairs)
}

fn normalize_ray(r: &[i64]) -> Vec<i64> {
    let mut g = 0i64;
    for &x in r {
        g = gcd_i64(g, x.abs());
    }
    let out: Vec<i64> = if g == 0 { r.to_vec() } else { r.iter().map(|&x| x / g).collect() };
    for &x in &out {
        if x > 0 {
            return out;
        } else if x < 0 {
            return out.iter().map(|&v| -v).collect();
        }
    }
    out
}

fn gcd_i64(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd_i64(b, a % b) }
}

/// Compute GV invariants using cygv.
///
/// This requires:
/// - Mori cone cap generators (rays)
/// - Grading vector (interior to dual cone)
/// - Curve basis matrix (q)
/// - Intersection numbers in basis (intnums)
pub fn compute_gv_invariants(
    rays: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let t0 = std::time::Instant::now();
    eprintln!(
        "[DEBUG] gv start: rays={}, dim={}, max_deg={:?}, min_points={:?}",
        rays.len(),
        rays.first().map(|r| r.len()).unwrap_or(0),
        max_deg,
        min_points
    );
    if min_points.is_none() && max_deg.is_none() {
        return Err(Error::InvalidInput(
            "Either min_points or max_deg must be specified".into(),
        ));
    }

    let dim = rays[0].len();

    // Report generator degree range for diagnostics.
    if let Some(d) = max_deg {
        let mut min_deg: Option<i128> = None;
        let mut max_deg_seen: Option<i128> = None;
        for r in rays {
            let deg: i128 = r
                .iter()
                .zip(grading_vector.iter())
                .map(|(&x, &g)| i128::from(x) * i128::from(g))
                .sum();
            min_deg = Some(min_deg.map_or(deg, |m| m.min(deg)));
            max_deg_seen = Some(max_deg_seen.map_or(deg, |m| m.max(deg)));
        }
        eprintln!(
            "[DEBUG] gv generators: total={}, degree_range={:?}-{:?}, max_deg={}",
            rays.len(),
            min_deg,
            max_deg_seen,
            d
        );
    }

    let n_rays = rays.len();
    eprintln!("[DEBUG] gv generators used: {}", n_rays);

    let grading_vec_i32: Vec<i32> = grading_vector
        .iter()
        .map(|&v| i32::try_from(v))
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| Error::InvalidInput("grading vector does not fit in i32".into()))?;
    let grading = RowDVector::from_row_slice(&grading_vec_i32);
    let mut min_g = i32::MAX;
    let mut max_g = i32::MIN;
    let mut neg_g = 0usize;
    for &v in &grading_vec_i32 {
        min_g = min_g.min(v);
        max_g = max_g.max(v);
        if v < 0 {
            neg_g += 1;
        }
    }
    eprintln!(
        "[DEBUG] gv grading_vec abs max={}, min={}, max={}, neg_count={}, elapsed={:.2?}",
        grading_vec_i32
            .iter()
            .map(|v| v.abs())
            .max()
            .unwrap_or(0),
        min_g,
        max_g,
        neg_g,
        t0.elapsed()
    );

    // Augment generators with lattice points, matching CYTools:
    // lattice_pts = mori.find_lattice_points(min_points=100*h11)
    let factor = env::var("CYRUS_LATTICE_MIN_POINTS_FACTOR")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(100);
    let gen_min_points = factor * dim;
    eprintln!("[DEBUG] gv generator min_points: {}", gen_min_points);

    let lattice_pts = {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        let mut key_rays = rays.to_vec();
        key_rays.sort();
        key_rays.hash(&mut hasher);
        grading_vector.hash(&mut hasher);
        gen_min_points.hash(&mut hasher);
        1000i64.hash(&mut hasher);
        0i64.hash(&mut hasher);
        let key = hasher.finish();

        let cache_dir = env::var("CYRUS_CACHE_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|_| PathBuf::from("target/cyrus-cache"));
        let cache_path = cache_dir.join(format!("lattice_points_{key:x}.json"));

        eprintln!("[DEBUG] lattice cache path: {}", cache_path.display());
        if cache_path.exists() {
            let data = fs::read_to_string(&cache_path).map_err(|e| {
                Error::InvalidInput(format!(
                    "Failed to read lattice point cache {}: {e}",
                    cache_path.display()
                ))
            })?;
            let pts: Vec<Vec<i64>> = serde_json::from_str(&data).map_err(|e| {
                Error::InvalidInput(format!(
                    "Failed to parse lattice point cache {}: {e}",
                    cache_path.display()
                ))
            })?;
            eprintln!(
                "[DEBUG] gv lattice points: {} (cache hit)",
                pts.len()
            );
            pts
        } else {
            let rays_i128: Vec<Vec<i128>> = rays
                .iter()
                .map(|r| r.iter().map(|&x| x as i128).collect())
                .collect();
            let mut cone = Cone::from_rays(rays_i128);
            let pts = cone.find_lattice_points_ortools(
                Some(gen_min_points),
                None,
                grading_vector,
                1000,
                0,
            )?;
            eprintln!(
                "[DEBUG] gv lattice points: {} (cache miss)",
                pts.len()
            );
            if let Err(e) = fs::create_dir_all(&cache_dir) {
                eprintln!(
                    "[WARN] failed to create lattice cache dir {}: {}",
                    cache_dir.display(),
                    e
                );
            } else {
                match fs::File::create(&cache_path) {
                    Ok(file) => {
                        let writer = BufWriter::new(file);
                        if let Err(e) = serde_json::to_writer(writer, &pts) {
                            eprintln!(
                                "[WARN] failed to serialize lattice cache {}: {}",
                                cache_path.display(),
                                e
                            );
                        }
                    }
                    Err(e) => {
                        eprintln!(
                            "[WARN] failed to create lattice cache {}: {}",
                            cache_path.display(),
                            e
                        );
                    }
                }
            }
            pts
        }
    };

    let mut all_generators: Vec<Vec<i64>> = Vec::new();
    for r in rays {
        all_generators.push(r.clone());
    }
    for p in lattice_pts {
        all_generators.push(p);
    }

    let mut uniq: HashSet<Vec<i64>> = HashSet::new();
    let mut uniq_vecs: Vec<Vec<i64>> = Vec::new();
    for v in all_generators {
        if uniq.insert(v.clone()) {
            uniq_vecs.push(v);
        }
    }

    let mut filtered_vecs = uniq_vecs;
    if let Some(d) = max_deg {
        let d_i128 = i128::from(d);
        let before = filtered_vecs.len();
        filtered_vecs = filtered_vecs
            .into_iter()
            .filter(|v| {
                let deg: i128 = v
                    .iter()
                    .zip(grading_vector.iter())
                    .map(|(&x, &g)| i128::from(x) * i128::from(g))
                    .sum();
                deg <= d_i128
            })
            .collect();
        eprintln!(
            "[DEBUG] gv generators filtered by max_deg: {} -> {}",
            before,
            filtered_vecs.len()
        );
        if filtered_vecs.is_empty() {
            return Err(Error::InvalidInput(
                "No generators remain after max_deg filtering".into(),
            ));
        }
    }

    let n_gen = filtered_vecs.len();
    let mut gen_data: Vec<i32> = Vec::with_capacity(dim * n_gen);
    for col in 0..n_gen {
        for row in 0..dim {
            let val = i32::try_from(filtered_vecs[col][row]).map_err(|_| {
                Error::InvalidInput("Mori cone generator does not fit in i32".into())
            })?;
            gen_data.push(val);
        }
    }
    let generators = DMatrix::from_column_slice(dim, n_gen, &gen_data);

    let generator_hash = {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        let mut key_vecs = filtered_vecs.clone();
        key_vecs.sort();
        key_vecs.hash(&mut hasher);
        hasher.finish()
    };

    // q matrix (curve basis), rows = h11, cols = n_pts
    let q_rows = q_matrix.len();
    if q_rows == 0 {
        return Err(Error::InvalidInput("q_matrix is empty".into()));
    }
    if q_rows != dim {
        return Err(Error::InvalidInput(
            "q_matrix row count must match generator dimension".into(),
        ));
    }
    let q_cols = q_matrix[0].len();
    let mut q_data: Vec<i32> = Vec::with_capacity(q_rows * q_cols);
    for row in q_matrix {
        if row.len() != q_cols {
            return Err(Error::InvalidInput("q_matrix rows have inconsistent length".into()));
        }
        for &v in row {
            let v_i32 = i32::try_from(v)
                .map_err(|_| Error::InvalidInput("q_matrix entry does not fit in i32".into()))?;
            q_data.push(v_i32);
        }
    }
    let q = DMatrix::from_row_slice(q_rows, q_cols, &q_data);
    // cygv expects q with shape (n_divisors, h11), i.e. transpose of curve basis.
    let q = q.transpose();
    eprintln!(
        "[DEBUG] gv q shape: {}x{}, elapsed={:.2?}",
        q.nrows(),
        q.ncols(),
        t0.elapsed()
    );

    // Intersection numbers (dok format)
    let mut intnums_map: HashMap<(usize, usize, usize), i32> = HashMap::new();
    for (&(i, j, k), val) in intnums.iter() {
        let (num, den) = val.get().clone().into_numerator_and_denominator();
        if den != 1u32 {
            return Err(Error::InvalidInput(
                "intersection number is not integral".into(),
            ));
        }
        let v_i64: i64 = i64::try_from(&num).map_err(|_| {
            Error::InvalidInput("intersection number does not fit in i64".into())
        })?;
        let v_i32: i32 = v_i64
            .try_into()
            .map_err(|_| Error::InvalidInput("intersection number does not fit in i32".into()))?;
        intnums_map.insert((i, j, k), v_i32);
    }

    let cache_key = {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        grading_vec_i32.hash(&mut hasher);
        max_deg.hash(&mut hasher);
        min_points.hash(&mut hasher);
        generator_hash.hash(&mut hasher);
        q.iter().for_each(|v| v.hash(&mut hasher));
        for (k, v) in &intnums_map {
            k.hash(&mut hasher);
            v.hash(&mut hasher);
        }
        hasher.finish()
    };
    let cache_dir = env::var("CYRUS_CACHE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("target/cyrus-cache"));
    let gv_cache_path = cache_dir.join(format!("gv_invariants_{cache_key:x}.json"));
    eprintln!("[DEBUG] gv cache path: {}", gv_cache_path.display());
    if gv_cache_path.exists() {
        let data = fs::read_to_string(&gv_cache_path).map_err(|e| {
            Error::InvalidInput(format!(
                "Failed to read GV cache {}: {e}",
                gv_cache_path.display()
            ))
        })?;
        #[derive(serde::Deserialize)]
        struct CachedGv {
            charge: Vec<i32>,
            value: String,
        }
        let items: Vec<CachedGv> = serde_json::from_str(&data).map_err(|e| {
            Error::InvalidInput(format!(
                "Failed to parse GV cache {}: {e}",
                gv_cache_path.display()
            ))
        })?;
        let mut out = Vec::with_capacity(items.len());
        for item in items {
            let gv_int = item
                .value
                .parse::<Integer>()
                .map_err(|_| Error::InvalidInput("GV cache integer parse failed".into()))?;
            out.push((item.charge, gv_int));
        }
        eprintln!(
            "[DEBUG] gv invariants: cache hit ({} entries)",
            out.len()
        );
        return Ok(out);
    }

    // No nef-partition for hypersurfaces
    let nefpart: Vec<DVector<i32>> = Vec::new();

    let n_threads = env::var("CYRUS_GV_THREADS")
        .ok()
        .and_then(|v| v.parse::<u32>().ok());
    let pool_size = env::var("CYRUS_GV_POOL_SIZE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1000);

    let invariants = cygv::compute_gv_rat_threefold(
        generators,
        grading,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums_map,
        n_threads,
        pool_size,
    );

    let mut out = Vec::with_capacity(invariants.len());
    for (v, gv) in invariants {
        let gv_str = gv.to_string();
        let gv_int = gv_str
            .parse::<Integer>()
            .map_err(|_| Error::InvalidInput("GV integer conversion failed".into()))?;
        out.push((v.as_slice().to_vec(), gv_int));
    }

    if let Err(e) = fs::create_dir_all(&cache_dir) {
        eprintln!(
            "[WARN] failed to create GV cache dir {}: {}",
            cache_dir.display(),
            e
        );
    } else {
        #[derive(serde::Serialize)]
        struct CachedGv<'a> {
            charge: &'a [i32],
            value: String,
        }
        let payload: Vec<CachedGv<'_>> = out
            .iter()
            .map(|(c, v)| CachedGv {
                charge: c,
                value: v.to_string(),
            })
            .collect();
        match fs::File::create(&gv_cache_path) {
            Ok(file) => {
                let writer = BufWriter::new(file);
                if let Err(e) = serde_json::to_writer(writer, &payload) {
                    eprintln!(
                        "[WARN] failed to serialize GV cache {}: {}",
                        gv_cache_path.display(),
                        e
                    );
                }
            }
            Err(e) => {
                eprintln!(
                    "[WARN] failed to create GV cache {}: {}",
                    gv_cache_path.display(),
                    e
                );
            }
        }
    }

    Ok(out)
}

fn compute_faces_4d(
    points: &[Point],
    polytope: &Polytope,
) -> Result<(Vec<Vec<usize>>, Vec<Vec<usize>>)> {
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "faces_4d only defined for 4D polytopes".into(),
        ));
    }

    let dual_vertices = polytope.dual_vertices()?;
    if dual_vertices.is_empty() {
        return Err(Error::InvalidInput("no dual vertices found".into()));
    }

    // Build facet vertex sets (polytope vertices only) to define 2-faces.
    let poly_vertices = polytope.vertices();
    let mut facet_vertex_sets: Vec<HashSet<usize>> = Vec::with_capacity(dual_vertices.len());
    let mut facets: Vec<Vec<usize>> = Vec::with_capacity(dual_vertices.len());
    for dv in &dual_vertices {
        let mut vert_set: HashSet<usize> = HashSet::new();
        for (idx, vtx) in poly_vertices.iter().enumerate() {
            let dot: i64 = vtx
                .coords()
                .iter()
                .zip(dv.coords().iter())
                .map(|(&a, &b)| a * b)
                .sum();
            if dot == -1 {
                vert_set.insert(idx);
            }
        }
        facet_vertex_sets.push(vert_set);

        // Facet points (triangulation points only, matching CYTools boundary_points()).
        let mut facet_pts: Vec<usize> = Vec::new();
        for (idx, pt) in points.iter().enumerate() {
            let dot: i64 = pt
                .coords()
                .iter()
                .zip(dv.coords().iter())
                .map(|(&a, &b)| a * b)
                .sum();
            if dot == -1 {
                facet_pts.push(idx);
            }
        }
        facet_pts.sort_unstable();
        facets.push(facet_pts);
    }

    // Build 2-faces via intersections of facet vertex sets (CYTools _faces4d).
    let mut twofaces: Vec<Vec<usize>> = Vec::new();
    for i in 0..facet_vertex_sets.len() {
        for j in (i + 1)..facet_vertex_sets.len() {
            let inter_vertices = facet_vertex_sets[i]
                .intersection(&facet_vertex_sets[j])
                .count();
            if inter_vertices >= 3 {
                // Collect triangulation points lying on both facets.
                let mut face_pts: Vec<usize> = Vec::new();
                for (idx, pt) in points.iter().enumerate() {
                    let dot_i: i64 = pt
                        .coords()
                        .iter()
                        .zip(dual_vertices[i].coords().iter())
                        .map(|(&a, &b)| a * b)
                        .sum();
                    if dot_i != -1 {
                        continue;
                    }
                    let dot_j: i64 = pt
                        .coords()
                        .iter()
                        .zip(dual_vertices[j].coords().iter())
                        .map(|(&a, &b)| a * b)
                        .sum();
                    if dot_j == -1 {
                        face_pts.push(idx);
                    }
                }
                face_pts.sort_unstable();
                twofaces.push(face_pts);
            }
        }
    }

    Ok((facets, twofaces))
}

fn nullspace_vector(
    pts_ext: &[Vec<i64>],
    diff_pts: &[usize],
    comm_pts: &[usize],
    require_unique: bool,
) -> Result<Option<Vec<Integer>>> {
    let rows = diff_pts.len() + comm_pts.len();
    let cols = pts_ext[0].len();
    let mut m: Vec<Vec<Integer>> = vec![vec![Integer::from(0); cols]; rows];

    for (r, &idx) in diff_pts.iter().enumerate() {
        for c in 0..cols {
            m[r][c] = Integer::from(pts_ext[idx][c]);
        }
    }
    for (r, &idx) in comm_pts.iter().enumerate() {
        for c in 0..cols {
            m[r + diff_pts.len()][c] = Integer::from(pts_ext[idx][c]);
        }
    }

    // Compute kernel of m^T
    let m_t = transpose(&m);
    let kernel = integer_kernel(&m_t);
    if kernel.is_empty() {
        return Ok(None);
    }
    if require_unique && kernel.len() != 1 {
        return Ok(None);
    }

    let mut v = kernel[0].clone();
    if v[0] < 0 {
        for val in &mut v {
            *val = -val.clone();
        }
    }

    let g = gcd_list(&v);
    if g != 0 {
        for val in &mut v {
            *val /= &g;
        }
    }

    Ok(Some(v))
}

fn build_full_v(diff_pts: &[usize], comm_pts: &[usize], v: &[Integer]) -> Vec<(usize, i64)> {
    let mut full_v: Vec<(usize, i64)> = Vec::new();
    for k in 0..diff_pts.len() {
        full_v.push((
            diff_pts[k],
            i64::try_from(&v[k]).expect("mori coeff fits in i64"),
        ));
    }
    for k in 0..comm_pts.len() {
        if v[k + diff_pts.len()] != 0 {
            full_v.push((
                comm_pts[k],
                i64::try_from(&v[k + diff_pts.len()]).expect("mori coeff fits in i64"),
            ));
        }
    }
    full_v.sort_by_key(|(idx, _)| *idx);
    full_v
}

fn transpose(m: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    if m.is_empty() {
        return Vec::new();
    }
    let rows = m.len();
    let cols = m[0].len();
    let mut t = vec![vec![Integer::from(0); rows]; cols];
    for r in 0..rows {
        for c in 0..cols {
            t[c][r] = m[r][c].clone();
        }
    }
    t
}

fn gcd_list(vals: &[Integer]) -> Integer {
    let mut g = Integer::from(0);
    for v in vals {
        let abs = v.clone().abs();
        if g == 0 {
            g = abs;
        } else if abs != 0 {
            g = gcd_integer(&g, &abs);
        }
    }
    g
}

fn is_strictly_dual(rays: &[Vec<i64>], v: &[i64]) -> bool {
    for r in rays {
        let dot: i128 = r
            .iter()
            .zip(v.iter())
            .map(|(&a, &b)| i128::from(a) * i128::from(b))
            .sum();
        if dot <= 0 {
            return false;
        }
    }
    true
}
