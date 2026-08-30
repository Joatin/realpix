//! Invariants that hold for every cell, independent of the reference implementation.

mod common;

use common::*;
use realpix::{Direction, MAX_DEPTH, from_uniq, nested, ring, to_uniq};

const EXHAUSTIVE_DEPTH: u8 = 6;

#[test]
fn center_round_trips_exhaustively() {
    for depth in 0..=EXHAUSTIVE_DEPTH {
        let n = nested::get(depth);
        let r = ring::get(depth);
        for cell in n.iter() {
            let (lon, lat) = n.center(cell);
            assert_eq!(n.hash(lon, lat), cell, "nested round trip, depth {depth}");
            assert_eq!(
                n.hash_vec(n.center_vec(cell)),
                cell,
                "nested vector round trip, depth {depth}"
            );
            let (lon, lat) = r.center(cell);
            assert_eq!(r.hash(lon, lat), cell, "ring round trip, depth {depth}");
            assert_eq!(
                r.hash_vec(r.center_vec(cell)),
                cell,
                "ring vector round trip, depth {depth}"
            );
        }
    }
}

#[test]
fn center_round_trips_at_every_depth() {
    let mut rng = Rng::new(0xC0FFEE);
    for depth in 0..=MAX_DEPTH {
        let n = nested::get(depth);
        let r = ring::get(depth);
        for _ in 0..2000 {
            let cell = rng.next_u64() % n.n_hash();
            let (lon, lat) = n.center(cell);
            assert_eq!(
                n.hash(lon, lat),
                cell,
                "nested round trip, depth {depth}, cell {cell}"
            );
            assert_eq!(
                n.hash_vec(n.center_vec(cell)),
                cell,
                "nested vec round trip, depth {depth}, cell {cell}"
            );
            let (lon, lat) = r.center(cell);
            assert_eq!(
                r.hash(lon, lat),
                cell,
                "ring round trip, depth {depth}, cell {cell}"
            );
            assert_eq!(
                r.hash_vec(r.center_vec(cell)),
                cell,
                "ring vec round trip, depth {depth}, cell {cell}"
            );
        }
    }
}

#[test]
fn position_round_trips_through_its_cell() {
    // Hashing a position, taking that cell's centre and hashing again must be idempotent.
    let mut rng = Rng::new(7);
    for depth in [0, 1, 5, 12, 20, 29] {
        let n = nested::get(depth);
        for _ in 0..20_000 {
            let (lon, lat) = rng.next_lonlat();
            let cell = n.hash(lon, lat);
            assert!(n.contains(cell));
            let (clon, clat) = n.center(cell);
            assert_eq!(
                n.hash(clon, clat),
                cell,
                "depth {depth}, lon {lon}, lat {lat}"
            );
        }
    }
}

#[test]
fn schemes_are_a_bijection() {
    for depth in 0..=EXHAUSTIVE_DEPTH {
        let n = nested::get(depth);
        let r = ring::get(depth);
        let mut seen = vec![false; n.n_hash() as usize];
        for cell in n.iter() {
            let ring_cell = n.to_ring(cell);
            assert!(
                r.contains(ring_cell),
                "depth {depth}: ring index out of range"
            );
            assert!(
                !seen[ring_cell as usize],
                "depth {depth}: ring index {ring_cell} hit twice"
            );
            seen[ring_cell as usize] = true;
            assert_eq!(r.to_nested(ring_cell), cell, "depth {depth}");
            // Both schemes must describe the same cell on the sphere.
            let d = ang_dist(n.center_vec(cell), r.center_vec(ring_cell));
            assert!(
                d < 1e-15,
                "depth {depth}: schemes disagree on centre by {d} rad"
            );
        }
        assert!(seen.into_iter().all(|s| s));
    }
}

#[test]
fn ring_indices_are_ordered_by_latitude() {
    for depth in 0..=5 {
        let r = ring::get(depth);
        let mut previous = f64::INFINITY;
        for cell in r.iter() {
            let z = r.center_vec(cell)[2];
            assert!(
                z <= previous + 1e-12,
                "depth {depth}, cell {cell}: latitude increased along the ring order"
            );
            previous = previous.min(z);
        }
    }
}

#[test]
fn rings_have_the_expected_lengths() {
    for depth in 1..=5 {
        let r = ring::get(depth);
        let nside = r.nside() as u64;
        // Group cells by latitude: each group must be one ring of the expected length.
        let mut lengths: Vec<u64> = Vec::new();
        let mut current_z: Option<f64> = None;
        let mut count = 0u64;
        for cell in r.iter() {
            let z = r.center_vec(cell)[2];
            if current_z.is_none_or(|c| (z - c).abs() > 1e-12) {
                if count > 0 {
                    lengths.push(count);
                }
                current_z = Some(z);
                count = 0;
            }
            count += 1;
        }
        lengths.push(count);

        assert_eq!(
            lengths.len() as u64,
            4 * nside - 1,
            "depth {depth}: wrong ring count"
        );
        for (i, len) in lengths.iter().enumerate() {
            let ring_index = i as u64 + 1;
            let expected = if ring_index < nside {
                4 * ring_index
            } else if ring_index > 3 * nside {
                4 * (4 * nside - ring_index)
            } else {
                4 * nside
            };
            assert_eq!(*len, expected, "depth {depth}, ring {ring_index}");
        }
    }
}

#[test]
fn cells_have_equal_area() {
    // Uniformly distributed directions must land in every cell equally often. This is the
    // property a merely-invertible but non-equal-area projection would fail.
    let depth = 3;
    let layer = nested::get(depth);
    let n = layer.n_hash() as usize;
    let samples = 400 * n;
    let mut counts = vec![0u32; n];
    let mut rng = Rng::new(0xABCDEF);
    for _ in 0..samples {
        counts[layer.hash_vec(rng.next_vec()) as usize] += 1;
    }

    let expected = samples as f64 / n as f64;
    // Poisson: sigma = sqrt(expected) = 20 for expected = 400, so 6 sigma is 30%.
    let sigma = expected.sqrt();
    let (mut worst, mut chi2) = (0.0f64, 0.0f64);
    for c in &counts {
        let dev = (*c as f64 - expected) / sigma;
        worst = worst.max(dev.abs());
        chi2 += dev * dev;
    }
    assert!(
        worst < 6.0,
        "cell counts deviate by {worst} sigma: cells are not equal area"
    );
    let reduced = chi2 / (n as f64 - 1.0);
    assert!(
        (0.7..1.3).contains(&reduced),
        "reduced chi-squared {reduced} is not consistent with equal-area cells"
    );
}

#[test]
fn neighbours_are_symmetric() {
    for depth in 0..=4 {
        let layer = nested::get(depth);
        for cell in layer.iter() {
            for (k, n) in layer.neighbours(cell).iter().enumerate() {
                let Some(n) = *n else { continue };
                assert!(layer.contains(n), "depth {depth}: neighbour out of range");
                assert_ne!(n, cell, "depth {depth}: cell {cell} is its own neighbour");
                assert!(
                    layer.neighbours(n).contains(&Some(cell)),
                    "depth {depth}: {cell} lists {n} as neighbour {k} but not the other way round"
                );
            }
        }
    }
}

#[test]
fn exactly_twenty_four_cells_have_seven_neighbours() {
    for depth in 0..=4 {
        let layer = nested::get(depth);
        let missing: usize = layer
            .iter()
            .map(|c| layer.neighbours(c).iter().filter(|n| n.is_none()).count())
            .sum();
        assert_eq!(missing, 24, "depth {depth}: wrong number of corner cells");
    }
}

#[test]
fn neighbours_are_adjacent_on_the_sphere() {
    let mut rng = Rng::new(11);
    for depth in [1, 4, 8, 12] {
        let layer = nested::get(depth);
        let bound = 2.5 * nested::MAX_CENTER_TO_VERTEX[depth as usize];
        for _ in 0..2000 {
            let cell = rng.next_u64() % layer.n_hash();
            let c = layer.center_vec(cell);
            for n in layer.neighbours(cell).into_iter().flatten() {
                let d = ang_dist(c, layer.center_vec(n));
                assert!(
                    d < bound,
                    "depth {depth}: neighbour {n} of {cell} is {d} rad away"
                );
            }
        }
        // The per-direction accessor must agree with the batch one.
        let cell = rng.next_u64() % layer.n_hash();
        let all = layer.neighbours(cell);
        for d in Direction::ALL {
            assert_eq!(layer.neighbour(cell, d), all[d.index()]);
        }
    }
}

#[test]
fn vertices_are_shared_and_enclose_the_centre() {
    for depth in 1..=4 {
        let layer = nested::get(depth);
        for cell in layer.iter() {
            let center = layer.center_vec(cell);
            let vertices = layer.vertices(cell);
            for v in vertices {
                // A vertex nudged towards the centre must fall back into this cell.
                let p = [
                    v[0] + (center[0] - v[0]) * 1e-6,
                    v[1] + (center[1] - v[1]) * 1e-6,
                    v[2] + (center[2] - v[2]) * 1e-6,
                ];
                assert_eq!(
                    layer.hash_vec(p),
                    cell,
                    "depth {depth}, cell {cell}: vertex outside its own cell"
                );
                // ... and must be shared with at least one neighbour.
                let shared = layer
                    .neighbours(cell)
                    .into_iter()
                    .flatten()
                    .any(|n| layer.vertices(n).iter().any(|w| ang_dist(*w, v) < 1e-12));
                assert!(
                    shared,
                    "depth {depth}, cell {cell}: vertex not shared with any neighbour"
                );
            }
            // The N/W/S/E ordering: north vertex highest, south lowest.
            assert!(
                vertices[0][2] >= vertices[2][2],
                "depth {depth}, cell {cell}: N below S"
            );
        }
    }
}

#[test]
fn hierarchy_is_consistent() {
    let mut rng = Rng::new(4242);
    for depth in [1, 5, 10, 20, 29] {
        let layer = nested::get(depth);
        for _ in 0..5000 {
            let cell = rng.next_u64() % layer.n_hash();
            let center = layer.center_vec(cell);
            for parent_depth in 0..=depth {
                let parent = layer.parent(cell, parent_depth);
                let parent_layer = nested::get(parent_depth);
                assert!(parent_layer.contains(parent));
                // The parent is the cell that contains this cell's centre.
                assert_eq!(
                    parent_layer.hash_vec(center),
                    parent,
                    "depth {depth} -> {parent_depth}"
                );
                // ... and this cell is inside the parent's descendant range.
                let range = parent_layer.children_range(parent, depth);
                assert!(range.contains(&cell), "depth {depth} -> {parent_depth}");
                assert_eq!(range.end - range.start, 1 << (2 * (depth - parent_depth)));
            }
            if depth < MAX_DEPTH {
                let children = layer.children(cell);
                let child_layer = nested::get(depth + 1);
                for child in children {
                    assert_eq!(child_layer.parent(child, depth), cell);
                }
                assert_eq!(children[0], layer.children_range(cell, depth + 1).start);
            }
        }
    }
}

#[test]
fn uniq_round_trips() {
    let mut rng = Rng::new(9);
    for depth in 0..=MAX_DEPTH {
        let n_hash = nested::get(depth).n_hash();
        for cell in [0, 1, n_hash - 1, rng.next_u64() % n_hash] {
            let uniq = to_uniq(depth, cell);
            assert_eq!(
                from_uniq(uniq),
                (depth, cell),
                "uniq for depth {depth}, cell {cell}"
            );
        }
    }
    // Uniq values never collide across depths.
    assert_ne!(to_uniq(0, 0), to_uniq(1, 0));
    assert_eq!(to_uniq(0, 0), 4);
    assert_eq!(to_uniq(0, 11), 15);
    assert_eq!(to_uniq(1, 0), 16);
}

#[test]
fn invalid_input_is_reported_not_returned() {
    let layer = nested::get(4);
    assert!(layer.checked_center(layer.n_hash()).is_err());
    assert!(layer.checked_center(layer.n_hash() - 1).is_ok());
    assert!(layer.checked_hash(f64::NAN, 0.0).is_err());
    assert!(layer.checked_hash(0.0, f64::INFINITY).is_err());
    assert!(layer.checked_hash(0.0, 2.0).is_err());
    assert!(layer.checked_hash(100.0, 1.0).is_ok());
    assert!(nested::checked_get(MAX_DEPTH + 1).is_err());
    assert!(ring::checked_get(MAX_DEPTH).is_ok());
    assert!(realpix::depth_from_nside(0).is_err());
    assert!(realpix::depth_from_nside(3).is_err());
    assert_eq!(realpix::depth_from_nside(1024), Ok(10));
}

#[test]
#[should_panic(expected = "out of range")]
fn out_of_range_cell_panics() {
    let layer = nested::get(4);
    let _ = layer.center(layer.n_hash());
}

#[test]
fn extreme_inputs_do_not_panic() {
    let layer = nested::get(12);
    use std::f64::consts::{FRAC_PI_2, TAU};
    for lat in [0.0, FRAC_PI_2, -FRAC_PI_2, FRAC_PI_2 + f64::EPSILON] {
        for lon in [0.0, -1e9, 1e9, TAU, -0.0] {
            let cell = layer.hash(lon, lat);
            assert!(layer.contains(cell), "lon {lon}, lat {lat} produced {cell}");
        }
    }
    for v in [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1e-300, 0.0, 1e-300],
    ] {
        assert!(layer.contains(layer.hash_vec(v)));
    }
}
