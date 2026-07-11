//! Landscape-run resume state: loading, revival, and persistence of the
//! per-arm bandit summary.

use cyrus_ga::multi::{PolytopeRecord, PolytopeStats};

/// Minimum interval between summary writes. The summary is the bandit's
/// resume state: write it compactly (the pretty form was 14 MB over 74k
/// arms), atomically (a kill mid-write must not truncate it - resume loses
/// ALL dead/best bookkeeping), and throttled (it was rewritten every round,
/// ~280 MB/s of serialization and I/O that starved the actual search).
pub const SUMMARY_EVERY: std::time::Duration = std::time::Duration::from_secs(30);

/// Serialize the bandit summary compactly and replace the file atomically.
pub fn write_summary_atomic(path: &std::path::Path, stats: &[PolytopeStats]) {
    let tmp = path.with_extension("json.tmp");
    std::fs::write(
        &tmp,
        serde_json::to_string(stats).expect("serialize summary"),
    )
    .expect("write summary");
    std::fs::rename(&tmp, path).expect("replace summary");
}

/// Load per-arm stats for the pool, merging any existing summary by NAME so
/// a grown pool (new chambers, new slices) resumes cleanly: existing arms
/// keep their history, new arms start fresh.
pub fn load_stats(summary_path: &std::path::Path, pool: &[PolytopeRecord]) -> Vec<PolytopeStats> {
    if !summary_path.exists() {
        return pool
            .iter()
            .map(|p| PolytopeStats::new(p.name.clone()))
            .collect();
    }
    let text = std::fs::read_to_string(summary_path).expect("read summary");
    let loaded: Vec<PolytopeStats> = serde_json::from_str(&text).expect("parse summary");
    let mut by_name: std::collections::HashMap<String, PolytopeStats> =
        loaded.into_iter().map(|s| (s.name.clone(), s)).collect();
    let merged: Vec<PolytopeStats> = pool
        .iter()
        .map(|p| {
            by_name
                .remove(&p.name)
                .unwrap_or_else(|| PolytopeStats::new(p.name.clone()))
        })
        .collect();
    eprintln!(
        "[INFO] resumed landscape run: {} arms ({} with history), {} rounds so far",
        merged.len(),
        merged.iter().filter(|s| s.rounds > 0 || s.dead).count(),
        merged.iter().map(|s| s.rounds).sum::<u64>()
    );
    merged
}

/// Retry transient prep deaths on every start: a timeout or a bare
/// subprocess exit is an environment verdict, not a physics one.
/// Deterministic physics rejections (PFV-barren, no viable orientifold)
/// stay dead. Revived never-searched arms have rounds == 0, so the
/// first-visit sweep re-probes them; genuinely obstructed geometries
/// re-die with a specific captured reason and are not revived again.
pub fn revive_transient_prep_deaths(stats: &mut [PolytopeStats]) {
    let mut revived = 0usize;
    for s in stats.iter_mut() {
        if !s.dead {
            continue;
        }
        let reason = s.dead_reason.as_deref().unwrap_or("");
        let permanent = reason.contains("PFV-barren") || reason.contains("no viable orientifold");
        if !permanent {
            s.dead = false;
            s.dead_reason = None;
            revived += 1;
        }
    }
    if revived > 0 {
        eprintln!("[INFO] revived {revived} arms with transient prep-failure deaths for retry");
    }
}
