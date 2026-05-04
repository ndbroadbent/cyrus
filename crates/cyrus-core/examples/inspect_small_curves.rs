#![allow(missing_docs)]

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

use cyrus_core::{
    F64, Finite, I64, Point, Polytope, Pos, compute_mori_cone_cap_rays,
    compute_regular_triangulation,
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

fn curve_volume_in_basis(curve: &[i64], basis: &[usize], t: &[F64<Finite>]) -> F64<Finite> {
    basis
        .iter()
        .zip(t.iter())
        .fold(F64::<Finite>::ZERO, |acc, (&idx, &ti)| {
            acc + I64::<Finite>::new(curve[idx]).to_f64() * ti
        })
}

fn main() {
    let data_dir = std::env::var("CYRUS_MCALLISTER_DATA_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            PathBuf::from(
                "/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647",
            )
        });
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

    let selected: Vec<(Vec<i64>, F64<Pos>)> = rays
        .iter()
        .filter_map(|ray| {
            let volume = curve_volume_in_basis(ray, &basis, &kahler).try_to_pos()?;
            (volume < cutoff).then_some((ray.clone(), volume))
        })
        .collect();

    let selected_set: HashSet<Vec<i64>> = selected.iter().map(|(ray, _)| ray.clone()).collect();
    let expected_set: HashSet<Vec<i64>> = expected_small.iter().cloned().collect();
    let overlap = selected_set.intersection(&expected_set).count();
    let missing = expected_set.difference(&selected_set).count();
    let extra = selected_set.difference(&expected_set).count();

    let min_volume = selected
        .iter()
        .map(|(_, volume)| volume.get())
        .fold(f64::INFINITY, f64::min);
    let max_volume = selected
        .iter()
        .map(|(_, volume)| volume.get())
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

    let mut extra_with_volume: Vec<(Vec<i64>, f64)> = selected
        .iter()
        .filter(|(ray, _)| !expected_set.contains(ray))
        .map(|(ray, volume)| (ray.clone(), volume.get()))
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
    }
}
