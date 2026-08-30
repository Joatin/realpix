//! The external edge: the ring of cells surrounding a cell at a deeper resolution.

mod common;

use common::*;
use realpix::{Direction, MAX_DEPTH, nested};
use std::collections::BTreeSet;

/// The external edge, by its definition: every cell at `edge_depth` that neighbours a cell
/// inside `cell` without being inside `cell` itself.
fn brute_force(depth: u8, cell: u64, edge_depth: u8) -> BTreeSet<u64> {
    let layer = nested::get(depth);
    let deep = nested::get(edge_depth);
    let inside = layer.children_range(cell, edge_depth);
    let mut out = BTreeSet::new();
    for descendant in inside.clone() {
        for neighbour in deep.neighbours(descendant).into_iter().flatten() {
            if !inside.contains(&neighbour) {
                out.insert(neighbour);
            }
        }
    }
    out
}

/// The load-bearing test: the walk must produce exactly the set the definition describes,
/// for every cell of the shallow layers and at every resolution up to four levels down.
#[test]
fn external_edge_matches_its_definition() {
    for depth in 0..=2u8 {
        let layer = nested::get(depth);
        for cell in layer.iter() {
            for edge_depth in depth..=(depth + 4) {
                let got = layer.external_edge_cells(cell, edge_depth);
                let expected = brute_force(depth, cell, edge_depth);
                assert_eq!(
                    got.iter().copied().collect::<BTreeSet<_>>(),
                    expected,
                    "depth {depth}, cell {cell}, edge depth {edge_depth}"
                );
            }
        }
    }
}

/// The same, sampled at depths too large to enumerate.
#[test]
fn external_edge_matches_its_definition_deeper() {
    let mut rng = Rng::new(0xED6E);
    for depth in [5u8, 9, 14, 20, 27] {
        let layer = nested::get(depth);
        for _ in 0..40 {
            let cell = layer.hash_vec(rng.next_vec());
            for edge_depth in [depth, depth + 1, depth + 2] {
                let got = layer.external_edge_cells(cell, edge_depth);
                assert_eq!(
                    got.iter().copied().collect::<BTreeSet<_>>(),
                    brute_force(depth, cell, edge_depth),
                    "depth {depth}, cell {cell}, edge depth {edge_depth}"
                );
            }
        }
    }
}

/// Each cell must be emitted exactly once, and the result sorted when collected.
#[test]
fn the_walk_emits_each_cell_once_and_the_collection_is_sorted() {
    let mut rng = Rng::new(11);
    for depth in [0u8, 3, 7, 12] {
        let layer = nested::get(depth);
        for _ in 0..30 {
            let cell = layer.hash_vec(rng.next_vec());
            for edge_depth in depth..=(depth + 3) {
                let mut walked = Vec::new();
                layer.external_edge(cell, edge_depth, |c| walked.push(c));

                let unique: BTreeSet<u64> = walked.iter().copied().collect();
                assert_eq!(
                    unique.len(),
                    walked.len(),
                    "depth {depth}, cell {cell}, edge depth {edge_depth}: duplicate emitted"
                );

                let collected = layer.external_edge_cells(cell, edge_depth);
                assert!(collected.windows(2).all(|w| w[0] < w[1]), "not sorted");
                assert_eq!(collected, unique.into_iter().collect::<Vec<_>>());
            }
        }
    }
}

/// The ring holds `4n + 4` cells, less one corner for each neighbour the cell itself
/// lacks — those are the points where only three base cells meet, so no cell is there.
#[test]
fn the_ring_loses_a_corner_wherever_the_cell_lacks_a_neighbour() {
    for depth in 0..=3u8 {
        let layer = nested::get(depth);
        let mut without_all_eight = 0;
        for cell in layer.iter() {
            let missing = layer
                .neighbours(cell)
                .iter()
                .filter(|n| n.is_none())
                .count() as u64;
            for edge_depth in depth..=(depth + 3) {
                let n = 1u64 << (edge_depth - depth);
                assert_eq!(
                    layer.external_edge_cells(cell, edge_depth).len() as u64,
                    4 * n + 4 - missing,
                    "depth {depth}, cell {cell}, edge depth {edge_depth}"
                );
            }
            if missing > 0 {
                without_all_eight += 1;
                // The sides of the ring are never short; only its corners can be.
                assert!(
                    missing <= 2,
                    "depth {depth}, cell {cell}: {missing} missing"
                );
            }
        }
        // Every base cell touches two of those points; deeper, the same 24 cells per layer
        // that `neighbours` reports seven neighbours for touch one each.
        assert_eq!(
            without_all_eight,
            if depth == 0 { 12 } else { 24 },
            "depth {depth}"
        );
    }
}

/// At its own depth the external edge is just the eight neighbours.
#[test]
fn at_the_same_depth_it_is_the_neighbours() {
    for depth in [0u8, 2, 6, 15, MAX_DEPTH] {
        let layer = nested::get(depth);
        let mut rng = Rng::new(depth as u64 + 1);
        for _ in 0..50 {
            let cell = layer.hash_vec(rng.next_vec());
            let neighbours: BTreeSet<u64> = layer.neighbours(cell).into_iter().flatten().collect();
            assert_eq!(
                layer
                    .external_edge_cells(cell, depth)
                    .into_iter()
                    .collect::<BTreeSet<_>>(),
                neighbours,
                "depth {depth}, cell {cell}"
            );
        }
    }
}

/// Nothing on the ring may be inside the cell, and everything on it must touch the cell:
/// its centre lies within one cell's reach of the cell's boundary.
#[test]
fn the_ring_surrounds_the_cell_without_overlapping_it() {
    let mut rng = Rng::new(505);
    for depth in [2u8, 5, 9] {
        let layer = nested::get(depth);
        for _ in 0..20 {
            let cell = layer.hash_vec(rng.next_vec());
            let center = layer.center_vec(cell);
            let inside = layer.children_range(cell, depth + 3);
            for edge in layer.external_edge_cells(cell, depth + 3) {
                assert!(
                    !inside.contains(&edge),
                    "ring cell {edge} is inside the cell"
                );
                // Just outside means within the cell's own reach plus one small cell.
                let reach = nested::MAX_CENTER_TO_VERTEX[depth as usize]
                    + 2.0 * nested::MAX_CENTER_TO_VERTEX[(depth + 3) as usize];
                let d = ang_dist(nested::get(depth + 3).center_vec(edge), center);
                assert!(d <= reach, "ring cell {edge} is {d} away, beyond {reach}");
            }
        }
    }
}

/// The ring of a cell is the union of the rings of its children, minus the children.
#[test]
fn the_ring_agrees_across_the_hierarchy() {
    let mut rng = Rng::new(9090);
    for depth in [1u8, 4, 8] {
        let layer = nested::get(depth);
        let below = nested::get(depth + 1);
        for _ in 0..25 {
            let cell = layer.hash_vec(rng.next_vec());
            let edge_depth = depth + 3;
            let inside = layer.children_range(cell, edge_depth);

            let mut from_children = BTreeSet::new();
            for child in layer.children(cell) {
                for c in below.external_edge_cells(child, edge_depth) {
                    if !inside.contains(&c) {
                        from_children.insert(c);
                    }
                }
            }
            assert_eq!(
                layer
                    .external_edge_cells(cell, edge_depth)
                    .into_iter()
                    .collect::<BTreeSet<_>>(),
                from_children,
                "depth {depth}, cell {cell}"
            );
        }
    }
}

/// The walk starts at the southern corner and goes once around, so consecutive cells are
/// neighbours of each other.
#[test]
fn the_walk_goes_once_around_the_ring() {
    let layer = nested::get(5);
    let mut rng = Rng::new(3131);
    for _ in 0..30 {
        let cell = layer.hash_vec(rng.next_vec());
        // Pick a cell away from the base-cell corners, where the ring is unbroken.
        if layer.neighbours(cell).iter().any(Option::is_none) {
            continue;
        }
        let edge_depth = 7;
        let deep = nested::get(edge_depth);
        let mut ring = Vec::new();
        layer.external_edge(cell, edge_depth, |c| ring.push(c));

        for pair in ring.windows(2) {
            let neighbours: BTreeSet<u64> =
                deep.neighbours(pair[0]).into_iter().flatten().collect();
            assert!(
                neighbours.contains(&pair[1]),
                "cell {} does not follow {} around the ring",
                pair[1],
                pair[0]
            );
        }
        // And it closes.
        let first_neighbours: BTreeSet<u64> =
            deep.neighbours(ring[0]).into_iter().flatten().collect();
        assert!(
            first_neighbours.contains(ring.last().unwrap()),
            "the ring does not close"
        );
    }
}

/// A worked case that is easy to check by eye: the neighbours of base cell 4, which sits
/// on the equator and is surrounded by eight others.
#[test]
fn a_base_cell_ring_is_what_it_should_be() {
    let base = nested::get(0);
    assert_eq!(base.external_edge_cells(4, 0), vec![0, 3, 5, 7, 8, 11]);
    assert_eq!(base.neighbours(4).into_iter().flatten().count(), 6);
    // Split base cell 4 into four, and its ring grows to 4*2 + 4 = 12, less the two
    // corners that do not exist.
    assert_eq!(base.external_edge_cells(4, 1).len(), 10);
    // Every cell of the ring is outside base cell 4.
    for c in base.external_edge_cells(4, 3) {
        assert_ne!(nested::get(3).parent(c, 0), 4);
    }
}

#[test]
#[should_panic(expected = "edge depth must be between")]
fn rejects_an_edge_depth_above_the_cell() {
    nested::get(5).external_edge_cells(0, 4);
}

#[test]
#[should_panic(expected = "edge depth must be between")]
fn rejects_an_edge_depth_past_the_deepest_layer() {
    nested::get(5).external_edge_cells(0, MAX_DEPTH + 1);
}

#[test]
#[should_panic(expected = "cell index out of range")]
fn rejects_a_cell_out_of_range() {
    nested::get(1).external_edge_cells(48, 2);
}

/// `Direction` still round-trips through the single-direction accessor.
#[test]
fn the_single_direction_accessor_agrees_with_the_full_set() {
    let mut rng = Rng::new(626);
    for depth in [0u8, 4, 11] {
        let layer = nested::get(depth);
        let rings = realpix::ring::get(depth);
        for _ in 0..50 {
            let cell = layer.hash_vec(rng.next_vec());
            let all = layer.neighbours(cell);
            for d in Direction::ALL {
                assert_eq!(layer.neighbour(cell, d), all[d.index()]);
            }
            let ring_cell = layer.to_ring(cell);
            let ring_all = rings.neighbours(ring_cell);
            for d in Direction::ALL {
                assert_eq!(rings.neighbour(ring_cell, d), ring_all[d.index()]);
            }
        }
    }
}
