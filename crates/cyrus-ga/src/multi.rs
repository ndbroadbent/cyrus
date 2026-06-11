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
    let exploration = 400.0;
    stats
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.dead)
        .map(|(idx, s)| {
            let score = if s.rounds == 0 {
                f64::INFINITY // every live polytope gets at least one round
            } else {
                let bonus =
                    exploration * ((total_rounds.max(2) as f64).ln() / s.rounds as f64).sqrt();
                s.best_fitness + bonus
            };
            (idx, score)
        })
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(idx, _)| idx)
}

/// Prepare a polytope's geometry, translating failure into a loud
/// dead-marking decision for the scheduler.
pub fn prepare_or_mark_dead(
    record: &PolytopeRecord,
    gv_min_points: u32,
    timeout: std::time::Duration,
    stats: &mut PolytopeStats,
) -> Option<GaGeometry> {
    // Geometry preparation can stall on pathological mirrors (observed: a
    // cygv closure hanging a landscape run for minutes). Run it on a worker
    // thread with a deadline; on timeout the thread is detached and the
    // polytope is marked dead - lost coverage, never a wrong answer.
    let points = record.points.clone();
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let result = GaGeometry::prepare_from_points(&points, gv_min_points);
        let _ = tx.send(result);
    });
    let outcome = match rx.recv_timeout(timeout) {
        Ok(result) => result,
        Err(_) => Err(format!("geometry preparation exceeded {timeout:?}")),
    };
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
