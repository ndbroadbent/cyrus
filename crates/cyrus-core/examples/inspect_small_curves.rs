#![allow(missing_docs)]

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::time::Instant;

use cyrus_core::{
    F64, Finite, Point, Polytope, Pos, compute_mori_cone_cap_rays, compute_regular_triangulation,
    find_pair_decomposition, remove_pair_decomposable_curve_candidates,
    subcutoff_toric_curve_candidates,
};

fn read_points(path: &Path) -> Vec<Vec<i64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<i64>().expect("integer entry"))
                .collect()
        })
        .collect()
}

fn read_csv_usize(path: &Path) -> Vec<usize> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
        .trim()
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<usize>().expect("usize entry"))
        .collect()
}

fn read_csv_heights(path: &Path) -> Vec<f64> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
        .trim()
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("f64 entry"))
        .collect()
}

fn read_csv_finite(path: &Path) -> Vec<F64<Finite>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
        .trim()
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .map(|s| {
            let value = s.trim().parse::<f64>().expect("f64 entry");
            F64::<Finite>::new(value)
                .unwrap_or_else(|| panic!("{} contains non-finite value {value}", path.display()))
        })
        .collect()
}

fn read_csv_pos(path: &Path) -> Vec<F64<Pos>> {
    read_csv_finite(path)
        .into_iter()
        .map(|value| {
            value.try_to_pos().unwrap_or_else(|| {
                panic!(
                    "{} contains non-positive value {}",
                    path.display(),
                    value.get()
                )
            })
        })
        .collect()
}

fn main() {
    let data_dir = std::env::var("CYRUS_MCALLISTER_DATA_DIR").map_or_else(
        |_| {
            PathBuf::from(
                "/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647",
            )
        },
        PathBuf::from,
    );
    let t0 = Instant::now();

    let points_raw = read_points(&data_dir.join("points.dat"));
    let heights = read_csv_heights(&data_dir.join("heights.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let kahler = read_csv_finite(&data_dir.join("kahler_param.dat"));
    let expected_small = read_points(&data_dir.join("small_curves.dat"));
    let cutoff = read_csv_pos(&data_dir.join("small_curves_cutoff.dat"))[0];

    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("triangulation points");
    let triangulation =
        compute_regular_triangulation(&triangulation_points, &heights).expect("triangulation");
    eprintln!("[TIME] triangulation: {:.2?}", t0.elapsed());

    let rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        false,
        false,
        None,
    )
    .expect("ambient mori cap rays");
    eprintln!(
        "[TIME] ambient mori cap: {:.2?} rays={} dim={}",
        t0.elapsed(),
        rays.len(),
        rays.first().map_or(0, Vec::len)
    );

    let selected =
        subcutoff_toric_curve_candidates(&rays, &basis, &kahler, cutoff).expect("subcutoff curves");
    let filtered = remove_pair_decomposable_curve_candidates(&selected).expect("Hilbert filter");

    let selected_set: HashSet<Vec<i64>> =
        selected.iter().map(|curve| curve.class.clone()).collect();
    let expected_set: HashSet<Vec<i64>> = expected_small.iter().cloned().collect();
    let overlap = selected_set.intersection(&expected_set).count();
    let missing = expected_set.difference(&selected_set).count();
    let extra = selected_set.difference(&expected_set).count();
    let filtered_set: HashSet<Vec<i64>> =
        filtered.iter().map(|curve| curve.class.clone()).collect();

    let min_volume = selected
        .iter()
        .map(|curve| curve.volume.get())
        .fold(f64::INFINITY, f64::min);
    let max_volume = selected
        .iter()
        .map(|curve| curve.volume.get())
        .fold(f64::NEG_INFINITY, f64::max);

    println!("cutoff={}", cutoff.get());
    println!(
        "selected={} expected={} overlap={} missing={} extra={}",
        selected.len(),
        expected_small.len(),
        overlap,
        missing,
        extra
    );
    println!("selected_volume_range=[{min_volume}, {max_volume}]");
    println!(
        "filtered={} filtered_matches_expected={}",
        filtered.len(),
        filtered_set == expected_set
    );

    let pair_decomposable = selected
        .iter()
        .filter(|curve| {
            find_pair_decomposition(curve, &selected)
                .expect("decomposition")
                .is_some()
        })
        .count();
    let pair_decomposable_expected = selected
        .iter()
        .filter(|curve| {
            expected_set.contains(&curve.class)
                && find_pair_decomposition(curve, &selected)
                    .expect("decomposition")
                    .is_some()
        })
        .count();
    let pair_decomposable_extra = selected
        .iter()
        .filter(|curve| {
            !expected_set.contains(&curve.class)
                && find_pair_decomposition(curve, &selected)
                    .expect("decomposition")
                    .is_some()
        })
        .count();
    println!(
        "pair_decomposable={pair_decomposable} expected_decomposable={pair_decomposable_expected} extra_decomposable={pair_decomposable_extra}"
    );

    let volume_by_class: HashMap<Vec<i64>, _> = selected
        .iter()
        .map(|curve| (curve.class.clone(), curve.volume))
        .collect();
    let mut extra_with_volume: Vec<_> = selected
        .iter()
        .filter(|curve| !expected_set.contains(&curve.class))
        .map(|curve| (curve.class.clone(), curve.volume.get()))
        .collect();
    extra_with_volume.sort_by(|a, b| a.1.total_cmp(&b.1));
    for (idx, (ray, volume)) in extra_with_volume.iter().take(10).enumerate() {
        let nonzero: Vec<(usize, i64)> = ray
            .iter()
            .copied()
            .enumerate()
            .filter(|(_, value)| *value != 0)
            .collect();
        println!("extra[{idx}] volume={volume} nonzero={nonzero:?}");
        let target = selected
            .iter()
            .find(|curve| curve.class == *ray)
            .expect("extra curve exists in selected list");
        if let Some((a, b)) = find_pair_decomposition(target, &selected).expect("decomposition") {
            let a_nonzero: Vec<(usize, i64)> = a
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, value)| *value != 0)
                .collect();
            let b_nonzero: Vec<(usize, i64)> = b
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, value)| *value != 0)
                .collect();
            let a_volume = volume_by_class
                .get(&a)
                .expect("decomposition summand is selected")
                .get();
            let b_volume = volume_by_class
                .get(&b)
                .expect("decomposition summand is selected")
                .get();
            println!(
                "  decomposition a={a_nonzero:?} volume={a_volume} b={b_nonzero:?} volume={b_volume}"
            );
        }
    }
}
