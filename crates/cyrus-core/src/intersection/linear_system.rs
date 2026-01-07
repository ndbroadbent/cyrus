//! Linear system building for intersection number computation.
//!
//! Functions for constructing the linear systems from GLSM constraints.

use malachite::Integer;
use std::collections::{HashMap, HashSet};

/// Build 4-form system using CYTools-style linear_relations matrix.
///
/// The linear_relations matrix is in reduced form:
/// - Columns correspond to non-origin points (1-indexed in the fixture, 0-indexed here after offset)
/// - Has -1 in pivot columns (non-basis divisors) and coefficients in basis columns
///
/// For each linear_relations row and each 3-face probe (j,k,l):
///   Σ_m Q_m κ^V_{mjkl} = 0
pub fn build_4form_system_with_linear_relations(
    probe_3faces: &HashSet<Vec<usize>>,
    known_4form: &HashMap<Vec<usize>, f64>,
    var_map: &HashMap<Vec<usize>, usize>,
    n_points: usize,
    origin_idx: Option<usize>,
    linear_relations: &[Vec<Integer>],
) -> (Vec<(usize, usize, f64)>, Vec<f64>) {
    let mut triplets = Vec::new();
    let mut c_vec = Vec::new();
    let mut row_idx = 0;

    let _n_lr_cols = linear_relations.first().map_or(0, std::vec::Vec::len);

    // linear_relations columns map to non-origin points
    // Column i in LR → point (i+1) if origin is at 0
    let col_to_point = |col: usize| -> usize {
        if origin_idx == Some(0) {
            col + 1 // LR col 0 → point 1, etc.
        } else {
            col // No origin adjustment
        }
    };

    let debug_probe = vec![3_usize, 4, 5];

    for probe in probe_3faces {
        let j = probe[0];
        let k = probe[1];
        let l = probe[2];

        let is_debug_probe = probe == &debug_probe;

        for (lr_row_idx, lr_row) in linear_relations.iter().enumerate() {
            let mut row_terms: Vec<(usize, f64)> = Vec::new();
            let mut rhs = 0.0;

            if is_debug_probe && lr_row_idx == 0 {
                eprintln!("\n[DEBUG] === Probe (3,4,5), LR row 0 ===");
            }

            for (col, q_val) in lr_row.iter().enumerate() {
                let q_m: f64 = q_val.to_string().parse().unwrap_or(0.0);
                if q_m.abs() < 1e-12 {
                    continue;
                }

                let m = col_to_point(col);
                if m >= n_points || Some(m) == origin_idx {
                    continue;
                }

                let mut key4 = vec![m, j, k, l];
                key4.sort_unstable();

                if let Some(&known_val) = known_4form.get(&key4) {
                    rhs -= q_m * known_val;
                    if is_debug_probe && lr_row_idx == 0 {
                        eprintln!(
                            "[DEBUG]   m={}: Q={} × κ^V{:?}={:.4} → RHS -= {:.4}",
                            m,
                            q_m,
                            key4,
                            known_val,
                            q_m * known_val
                        );
                    }
                } else if let Some(&var_idx) = var_map.get(&key4) {
                    row_terms.push((var_idx, q_m));
                    if is_debug_probe && lr_row_idx == 0 {
                        eprintln!("[DEBUG]   m={m}: Q={q_m} × κ^V{key4:?}[var#{var_idx}]");
                    }
                } else if is_debug_probe && lr_row_idx == 0 {
                    eprintln!("[DEBUG]   m={m}: Q={q_m} × κ^V{key4:?} = 0 (not in simplex)");
                }
            }

            if is_debug_probe && lr_row_idx == 0 {
                eprintln!("[DEBUG]   => RHS = {rhs:.6}");
                eprintln!("[DEBUG]   => {} variable terms", row_terms.len());
            }

            // Add equation (even if no variables - let solver handle)
            if !row_terms.is_empty() {
                for (var_idx, coeff) in row_terms {
                    triplets.push((row_idx, var_idx, coeff));
                }
                c_vec.push(rhs);
                row_idx += 1;
            }
        }
    }

    eprintln!("[DEBUG] Total equations from linear_relations: {row_idx}");

    (triplets, c_vec)
}

/// Build GLSM system for 4-form variables using 3-face probes.
///
/// **Origin Exclusion Algorithm**: The origin is completely excluded from the system.
/// - GLSM matrix is used WITHOUT the origin column
/// - Equations sum over non-origin indices only: Σ_{m≠0} Q_m κ^V_{mjkl} = 0
/// - No "effective coefficients" (Q_m - Q_0) are used
///
/// For each GLSM row a and each 3-face probe (j,k,l) where j,k,l ≠ 0:
///   Σ_{m≠0} Q^a_m κ^V_{mjkl} = 0
#[allow(clippy::too_many_arguments)]
pub fn build_4form_glsm_system_no_origin(
    probe_3faces: &HashSet<Vec<usize>>,
    known_4form: &HashMap<Vec<usize>, f64>,
    var_map: &HashMap<Vec<usize>, usize>,
    n_points: usize,
    origin_idx: Option<usize>,
    glsm: &[Vec<Integer>],
    glsm_point_offset: Option<i32>,
) -> (Vec<(usize, usize, f64)>, Vec<f64>) {
    let mut triplets = Vec::new();
    let mut c_vec = Vec::new();
    let mut row_idx = 0;

    // Determine GLSM column to point index mapping
    let n_glsm_cols = glsm.first().map_or(0, std::vec::Vec::len);
    let glsm_offset: i32 = glsm_point_offset.unwrap_or(if n_glsm_cols > n_points { -1 } else { 0 });

    // Find the origin column in GLSM (to skip it)
    let origin_glsm_col: Option<usize> = origin_idx.and_then(|oi| {
        let glsm_col = (oi as i32) - glsm_offset;
        if glsm_col >= 0 && (glsm_col as usize) < n_glsm_cols {
            Some(glsm_col as usize)
        } else {
            None
        }
    });

    eprintln!("[DEBUG] n_points={n_points}, n_glsm_cols={n_glsm_cols}, glsm_offset={glsm_offset}");
    eprintln!("[DEBUG] Origin GLSM column (excluded): {origin_glsm_col:?}");

    // Debug: check probe (3,4,5) specifically
    let debug_probe = vec![3_usize, 4, 5];

    for probe in probe_3faces {
        let j = probe[0];
        let k = probe[1];
        let l = probe[2];

        let is_debug_probe = probe == &debug_probe;

        for (glsm_row_idx, q_row) in glsm.iter().enumerate() {
            let mut row_terms: Vec<(usize, f64)> = Vec::new();
            let mut rhs = 0.0;

            if is_debug_probe && glsm_row_idx == 0 {
                eprintln!("\n[DEBUG] === Probe (3,4,5), GLSM row 0 (no origin) ===");
            }

            for glsm_col in 0..q_row.len() {
                // Skip origin column entirely
                if Some(glsm_col) == origin_glsm_col {
                    continue;
                }

                #[allow(clippy::cast_possible_wrap)]
                let m_signed = glsm_col as i32 + glsm_offset;
                if m_signed < 0 {
                    continue;
                }

                #[allow(clippy::cast_sign_loss)]
                let m = m_signed as usize;
                if m >= n_points || Some(m) == origin_idx {
                    continue;
                }

                // Use raw Q_m (NOT effective coefficient Q_m - Q_0)
                let q_m: f64 = q_row[glsm_col].to_string().parse().unwrap_or(0.0);

                if q_m.abs() < 1e-12 {
                    continue;
                }

                // Build sorted 4-tuple (m, j, k, l)
                let mut key4 = vec![m, j, k, l];
                key4.sort_unstable();

                // Check if this is known or a variable
                if let Some(&known_val) = known_4form.get(&key4) {
                    // Known value -> move to RHS
                    rhs -= q_m * known_val;
                    if is_debug_probe && glsm_row_idx == 0 {
                        eprintln!(
                            "[DEBUG]   m={}: Q_m={} * κ^V{:?}={:.4} -> RHS -= {:.4}",
                            m,
                            q_m,
                            key4,
                            known_val,
                            q_m * known_val
                        );
                    }
                } else if let Some(&var_idx) = var_map.get(&key4) {
                    // Variable
                    row_terms.push((var_idx, q_m));
                    if is_debug_probe && glsm_row_idx == 0 {
                        eprintln!("[DEBUG]   m={m}: Q_m={q_m} * κ^V{key4:?}[var#{var_idx}]");
                    }
                } else if is_debug_probe && glsm_row_idx == 0 {
                    eprintln!("[DEBUG]   m={m}: Q_m={q_m} * κ^V{key4:?} = 0 (not in simplex)");
                }
            }

            if is_debug_probe && glsm_row_idx == 0 {
                eprintln!("[DEBUG]   => RHS = {rhs:.6}");
                eprintln!("[DEBUG]   => {} variable terms", row_terms.len());
            }

            // Add equation if it has variables
            if !row_terms.is_empty() {
                for (var_idx, coeff) in row_terms {
                    triplets.push((row_idx, var_idx, coeff));
                }
                c_vec.push(rhs);
                row_idx += 1;
            }
        }
    }

    eprintln!("[DEBUG] Total 4-form equations: {row_idx}");

    (triplets, c_vec)
}
