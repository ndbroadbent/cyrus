//! Multi-polytope landscape search.
//!
//! The outer loop over geometries that the GA_v2 research notes call for:
//! a pool of Kreuzer-Skarke polytopes (JSONL), a UCB-style bandit that
//! allocates inner-GA rounds between exploring fresh geometries and
//! exploiting promising ones, and per-polytope persistent state so the
//! whole landscape run stops and resumes like a single-geometry run.
//! Geometry preparation failures mark a polytope dead (loudly) and the
//! search moves on - in a landscape scan, a skipped geometry is lost
//! coverage, never a wrong answer.

use serde::{Deserialize, Serialize};

use crate::geometry::GaGeometry;

/// One polytope record from the pool file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolytopeRecord {
    /// Identifier used for run-dir subdirectories and reporting.
    pub name: String,
    /// Hodge number h11 of the CY hypersurface (informational).
    #[serde(default)]
    pub h11: Option<usize>,
    /// Hodge number h21 (the flux vector length on the mirror side).
    #[serde(default)]
    pub h21: Option<usize>,
    /// Whether the polytope is favorable (informational).
    #[serde(default)]
    pub favorable: Option<bool>,
    /// Full CYTools-ordered lattice point list.
    pub points: Vec<Vec<i64>>,
    /// Mirror (dual) chamber index: 0 = default dual FRST, k > 0 = k wall
    /// flips away. Different chambers have different PFV solution sets.
    #[serde(default)]
    pub chamber: usize,
}

/// Bandit bookkeeping for one polytope (persisted in the global summary).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolytopeStats {
    /// Pool identifier.
    pub name: String,
    /// Rounds of inner GA this polytope has received.
    pub rounds: u64,
    /// Total fitness evaluations spent here.
    pub evaluations: u64,
    /// Best fitness found on this polytope.
    #[serde(with = "crate::population::serde_finite_or_neg_inf")]
    pub best_fitness: f64,
    /// Count of valid vacua ever seen here.
    pub valid_seen: u64,
    /// Geometry preparation failed; never scheduled again.
    pub dead: bool,
    /// Failure description when dead.
    pub dead_reason: Option<String>,
}

impl PolytopeStats {
    /// Fresh stats for a pool entry.
    #[must_use]
    pub const fn new(name: String) -> Self {
        Self {
            name,
            rounds: 0,
            evaluations: 0,
            best_fitness: f64::NEG_INFINITY,
            valid_seen: 0,
            dead: false,
            dead_reason: None,
        }
    }
}

/// Deterministic uniform shuffle of the bandit's visit order.
///
/// A sorted (h21, h11) pool traps the first-visit sweep in the
/// systematically-barren high-h11 tail for hours (the cone-interior
/// condition q.p > 0 over every GV curve is hard to satisfy at large
/// h11). A uniform shuffle is the least-biased order - we do not know
/// which region holds candidates - and interleaves fertile geometries
/// with barren culls so dashboard progress is steady.
pub fn shuffle_visit_order(pool: &mut [PolytopeRecord]) {
    use rand::SeedableRng as _;
    use rand::seq::SliceRandom as _;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(987_654_321);
    pool.shuffle(&mut rng);
}

/// Read a JSONL polytope pool file.
///
/// # Errors
/// Returns a description if the file is unreadable or malformed.
pub fn load_pool(path: &std::path::Path) -> Result<Vec<PolytopeRecord>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    let mut records = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let record: PolytopeRecord = serde_json::from_str(line)
            .map_err(|e| format!("{}:{}: {e}", path.display(), line_no + 1))?;
        records.push(record);
    }
    if records.is_empty() {
        return Err(format!("{} contains no polytopes", path.display()));
    }
    Ok(records)
}

/// UCB-style polytope selection.
///
/// Balances the best fitness seen on a polytope against an exploration
/// bonus for under-visited ones. Fitness tiers span roughly [-2000, +300],
/// so the exploration constant is on that scale.
#[must_use]
pub fn select_next(stats: &[PolytopeStats], total_rounds: u64) -> Option<usize> {
    // First visits go in pool order (pools sort by ascending h21, where the
    // isotropic-flux lattice is densest), but every third round is reserved
    // for UCB exploitation so a fertile polytope found early gets deepened
    // even while the breadth sweep is still running.
    let any_visited = stats.iter().any(|s| !s.dead && s.rounds > 0);
    if !(any_visited && total_rounds % 3 == 2)
        && let Some(idx) = stats.iter().position(|s| !s.dead && s.rounds == 0)
    {
        return Some(idx);
    }
    // UCB over VISITED polytopes only: unvisited ones are the first-visit
    // branch's job (treating them as infinity here made exploit rounds
    // burn budget preparing far-pool geometries instead of deepening known
    // fertile ones). Falls back to the first unvisited when nothing has
    // been visited yet.
    let exploration = 400.0;
    stats
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.dead && s.rounds > 0)
        .map(|(idx, s)| {
            let bonus = exploration * ((total_rounds.max(2) as f64).ln() / s.rounds as f64).sqrt();
            (idx, s.best_fitness + bonus)
        })
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(idx, _)| idx)
        .or_else(|| stats.iter().position(|s| !s.dead && s.rounds == 0))
}

/// Prepare a polytope's geometry, translating failure into a loud
/// dead-marking decision for the scheduler.
pub fn prepare_or_mark_dead(
    record: &PolytopeRecord,
    pool_path: &std::path::Path,
    gv_min_points: u32,
    timeout: std::time::Duration,
    stats: &mut PolytopeStats,
) -> Option<GaGeometry> {
    // Geometry preparation can stall on pathological mirrors (observed: a
    // 5-vertex near-simplex whose dual cone dualization runs 14k
    // hyperplanes). A detached thread cannot be killed and leaks a busy
    // core per timeout, so the probe runs in a KILLABLE subprocess; on
    // success the in-process re-preparation is warm (GV and lattice-point
    // caches persist to disk).
    let exe = std::env::current_exe().expect("current executable path");
    let mut child = std::process::Command::new(exe)
        .arg("--prep-probe")
        .arg(&record.name)
        .arg("--polytope-file")
        .arg(pool_path)
        .arg("--chamber")
        .arg(record.chamber.to_string())
        .arg("--gv-min-points")
        .arg(gv_min_points.to_string())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("spawn prep probe");
    let deadline = std::time::Instant::now() + timeout;
    let outcome = loop {
        match child.try_wait().expect("poll prep probe") {
            Some(status) if status.success() => break Ok(()),
            Some(status) => break Err(format!("geometry preparation failed ({status})")),
            None if std::time::Instant::now() >= deadline => {
                let _ = child.kill();
                let _ = child.wait();
                break Err(format!("geometry preparation exceeded {timeout:?}"));
            }
            None => std::thread::sleep(std::time::Duration::from_millis(100)),
        }
    };
    let outcome = outcome.and_then(|()| {
        // Probe succeeded: re-prepare in-process against warm disk caches.
        GaGeometry::prepare_from_points_in_chamber(
            &record.points,
            gv_min_points,
            record.chamber,
            true,
        )
    });
    match outcome {
        Ok(geom) => Some(geom),
        Err(reason) => {
            eprintln!(
                "[WARN] polytope {} geometry preparation failed; marking dead: {reason}",
                record.name
            );
            stats.dead = true;
            stats.dead_reason = Some(reason);
            None
        }
    }
}

/// Boundary validation of pool data conventions.
///
/// Computes the first record's Hodge data from its geometry and compares
/// against the declared fields. Catches lattice convention mistakes (M vs
/// N polytope: h11/h21 swap) and incomplete point data BEFORE any compute
/// is spent - both happened in production.
pub fn validate_pool_conventions(pool: &[PolytopeRecord]) {
    let Some(record) = pool.first() else { return };
    let (Some(decl_h11), Some(decl_h21)) = (record.h11, record.h21) else {
        eprintln!("[WARN] pool records carry no declared Hodge numbers; skipping validation");
        return;
    };
    use cyrus_core::{Point, Polytope};
    let points: Vec<Point> = record
        .points
        .iter()
        .map(|p| Point::new(p.clone()))
        .collect();
    let computed = (|| -> Result<(usize, usize), String> {
        let primal = Polytope::from_vertices(points)
            .map_err(|e| e.to_string())?
            .compute_dual()
            .map_err(|e| e.to_string())?
            .compute_dual()
            .map_err(|e| e.to_string())?;
        let tri_points = primal
            .points_not_interior_to_facets()
            .map_err(|e| e.to_string())?;
        let (_, _, basis) =
            cyrus_core::compute_glsm_and_linrels(&tri_points).map_err(|e| e.to_string())?;
        let extra = cyrus_core::divisor::batyrev_h11_extra_classes(&primal, &tri_points)
            .map_err(|e| e.to_string())?;
        let h11 = basis.len() + extra;
        let dual = primal.compute_dual().map_err(|e| e.to_string())?;
        let dual_tri = dual
            .points_not_interior_to_facets()
            .map_err(|e| e.to_string())?;
        let (_, _, dual_basis) =
            cyrus_core::compute_glsm_and_linrels(&dual_tri).map_err(|e| e.to_string())?;
        Ok((h11, dual_basis.len()))
    })();
    match computed {
        Ok((h11, h21)) if h11 == decl_h11 && h21 == decl_h21 => {
            eprintln!(
                "[INFO] pool validation OK: {} computes to h11={h11} h21={h21} as declared",
                record.name
            );
        }
        Ok((h11, h21)) => {
            eprintln!(
                "[ERROR] pool convention mismatch: {} declares h11={decl_h11} h21={decl_h21} but computes h11={h11} h21={h21}. Wrong lattice (M vs N) or incomplete points?",
                record.name
            );
            std::process::exit(2);
        }
        Err(e) => {
            eprintln!("[ERROR] pool validation failed on {}: {e}", record.name);
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn select_prefers_unvisited_then_best_with_exploration() {
        let mut stats = vec![
            PolytopeStats::new("a".into()),
            PolytopeStats::new("b".into()),
        ];
        // Both unvisited: any unvisited polytope is acceptable.
        let first = select_next(&stats, 0).expect("selects something");
        assert_eq!(stats[first].rounds, 0);
        stats[0].rounds = 1;
        stats[0].best_fitness = 100.0;
        // b unvisited -> selected next.
        assert_eq!(select_next(&stats, 1), Some(1));
        stats[1].rounds = 1;
        stats[1].best_fitness = -100.0;
        // Both visited once: a's fitness dominates the small bonus gap.
        assert_eq!(select_next(&stats, 2), Some(0));
        // Heavy exploitation of a shrinks its bonus while b's grows with
        // total rounds: exploration eventually wins.
        stats[0].rounds = 10_000;
        assert_eq!(select_next(&stats, 10_001), Some(1));
    }

    #[test]
    fn dead_polytopes_are_never_selected() {
        let mut stats = vec![PolytopeStats::new("a".into())];
        stats[0].dead = true;
        assert_eq!(select_next(&stats, 5), None);
    }

    #[test]
    fn unvisited_stats_roundtrip_through_json() {
        let stats = vec![PolytopeStats::new("a".into())];
        let json = serde_json::to_string(&stats).unwrap();
        let resumed: Vec<PolytopeStats> = serde_json::from_str(&json).unwrap();
        assert_eq!(resumed[0].best_fitness, f64::NEG_INFINITY);
    }

    #[test]
    fn pool_file_roundtrip() {
        let dir = std::env::temp_dir().join("cyrus_ga_pool_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("pool.jsonl");
        std::fs::write(
            &path,
            "{\"name\":\"p0\",\"points\":[[0,0,0,0],[1,0,0,0]]}\n",
        )
        .unwrap();
        let pool = load_pool(&path).unwrap();
        assert_eq!(pool.len(), 1);
        assert_eq!(pool[0].name, "p0");
        assert_eq!(pool[0].points.len(), 2);
    }
}
