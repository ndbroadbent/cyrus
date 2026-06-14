//! Per-involution orientifold diagnostic for one pool candidate.
//!
//! Reconstructs the geometry exactly as emit_verification_dir does, then
//! enumerates all involutions and reports, for each, why it is or is not a
//! viable McAllister-class orientifold (O7 count, O7 disjointness, and how
//! many kklt_basis divisors are rigid).
//!
//!   cargo run --release -p cyrus-ga --example orientifold_diag -- <pool.jsonl> <name> <h21>

use cyrus_core::orientifold::{enumerate_involutions, select_orientifold};
use cyrus_core::{Point, Polytope};
use cyrus_ga::multi::load_pool;

fn main() {
    let mut a = std::env::args().skip(1);
    let pool_path = a.next().expect("pool path");
    let name = a.next().expect("name");
    let h21: usize = a.next().expect("h21").parse().expect("h21 int");

    let pool = load_pool(std::path::Path::new(&pool_path)).expect("pool");
    let record = pool.iter().find(|p| p.name == name).expect("record");

    let point_objs: Vec<Point> = record
        .points
        .iter()
        .map(|p| Point::new(p.clone()))
        .collect();
    let polytope = Polytope::from_vertices(point_objs)
        .expect("polytope")
        .compute_dual()
        .expect("dual")
        .compute_dual()
        .expect("double dual");
    let tri_points = polytope
        .points_not_interior_to_facets()
        .expect("tri points");
    let origin_idx = tri_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("origin");
    let (_heights, triangulation) =
        cyrus_core::compute_frst_heights(&tri_points, origin_idx).expect("FRST");
    let (_, _, basis) = cyrus_core::compute_glsm_and_linrels(&tri_points).expect("glsm");
    let points_i64: Vec<Vec<i64>> = tri_points.iter().map(|p| p.coords().to_vec()).collect();
    let linrels = cyrus_core::compute_linear_relations_no_origin(&points_i64);
    let linrels_i64: Vec<Vec<i64>> = linrels
        .iter()
        .map(|r| r.iter().map(|x| i64::try_from(x).expect("i64")).collect())
        .collect();
    let kappa = cyrus_core::compute_intersection_cytools(&triangulation, &tri_points, &linrels_i64)
        .expect("intersection");
    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&polytope, &tri_points).expect("batyrev");
    let h11 = basis.len() + h11_extra;

    println!(
        "{name}: h11={h11} (basis {} + extra {h11_extra}), h21={h21}",
        basis.len()
    );
    println!("basis (divisor indices): {basis:?}");

    let models =
        enumerate_involutions(&polytope, &tri_points, &kappa, &basis, h11, h21).expect("enumerate");
    println!("\n{} involutions enumerated:\n", models.len());
    println!(
        "{:<14}{:>5}{:>10}{:>14}{:>8}",
        "sigma", "#O7", "disjoint", "basis_rigid", "viable"
    );
    for m in &models {
        let rigid = m
            .kklt_basis
            .iter()
            .filter(|&&idx| {
                m.divisor_classes
                    .iter()
                    .find(|(i, _)| *i == idx)
                    .is_some_and(|(_, c)| c.is_rigid())
            })
            .count();
        println!(
            "{:<14}{:>5}{:>10}{:>10}/{:<3}{:>8}",
            format!("{:?}", m.sigma),
            m.o7_points.len(),
            m.o7_pairwise_disjoint,
            rigid,
            m.kklt_basis.len(),
            m.viable
        );
    }

    match select_orientifold(&models) {
        Some(m) => println!(
            "\nSELECTED: sigma={:?}, {} O7-planes, Q_D3={}",
            m.sigma,
            m.o7_points.len(),
            m.q_d3
        ),
        None => println!(
            "\nNO VIABLE ORIENTIFOLD: every involution fails (no O7s, intersecting O7s, or a non-rigid basis divisor)."
        ),
    }

    // For the closest near-misses, show which basis divisors are NOT rigid.
    println!("\nNon-rigid basis divisors per involution that has O7s:");
    for m in models.iter().filter(|m| !m.o7_points.is_empty()) {
        let bad: Vec<(usize, String)> = m
            .kklt_basis
            .iter()
            .filter_map(|&idx| {
                m.divisor_classes
                    .iter()
                    .find(|(i, _)| *i == idx)
                    .filter(|(_, c)| !c.is_rigid())
                    .map(|(_, c)| (idx, format!("{c:?}")))
            })
            .collect();
        if !bad.is_empty() {
            println!("  sigma={:?}: non-rigid {bad:?}", m.sigma);
        }
    }
}
