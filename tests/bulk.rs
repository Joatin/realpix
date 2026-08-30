//! The bulk hashing entry points: they must agree with the scalar ones exactly.

mod common;

use common::*;
use realpix::{MAX_DEPTH, nested, ring};

/// Bulk and scalar must agree bit for bit — the whole point is that it is the same
/// operation, so any divergence is a bug rather than a rounding difference.
#[test]
fn bulk_agrees_with_scalar_at_every_depth() {
    let mut rng = Rng::new(0xB0_1C);
    let positions: Vec<(f64, f64)> = (0..500).map(|_| rng.next_lonlat()).collect();
    let vectors: Vec<[f64; 3]> = positions
        .iter()
        .map(|(lon, lat)| lonlat_to_vec(*lon, *lat))
        .collect();

    let mut out = vec![0u64; positions.len()];
    for depth in 0..=MAX_DEPTH {
        let n = nested::get(depth);
        n.hash_many(&positions, &mut out);
        for (cell, (lon, lat)) in out.iter().zip(&positions) {
            assert_eq!(*cell, n.hash(*lon, *lat), "nested depth {depth}");
        }
        n.hash_many_vec(&vectors, &mut out);
        for (cell, v) in out.iter().zip(&vectors) {
            assert_eq!(*cell, n.hash_vec(*v), "nested depth {depth}, vector");
        }

        let r = ring::get(depth);
        r.hash_many(&positions, &mut out);
        for (cell, (lon, lat)) in out.iter().zip(&positions) {
            assert_eq!(*cell, r.hash(*lon, *lat), "ring depth {depth}");
        }
        r.hash_many_vec(&vectors, &mut out);
        for (cell, v) in out.iter().zip(&vectors) {
            assert_eq!(*cell, r.hash_vec(*v), "ring depth {depth}, vector");
        }
    }
}

/// Awkward inputs must behave the same in bulk as one at a time, including the ones that
/// are not really positions at all.
#[test]
fn bulk_handles_the_awkward_positions() {
    let poles_and_seams = [
        (0.0, std::f64::consts::FRAC_PI_2),
        (0.0, -std::f64::consts::FRAC_PI_2),
        (std::f64::consts::TAU, 0.0),
        (-1.0, 0.0),
        (1e9, 0.7),
        (0.0, 0.0),
    ];
    let vectors = [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        // Not unit length: the scalar entry point normalises, so bulk must too.
        [0.0, 0.0, 7.5],
        [3.0, 4.0, 0.0],
    ];

    for depth in [0u8, 1, 7, 20, MAX_DEPTH] {
        let layer = nested::get(depth);
        let mut out = vec![0u64; poles_and_seams.len()];
        layer.hash_many(&poles_and_seams, &mut out);
        for (cell, (lon, lat)) in out.iter().zip(&poles_and_seams) {
            assert_eq!(*cell, layer.hash(*lon, *lat), "depth {depth}");
            assert!(layer.contains(*cell));
        }

        let mut out = vec![0u64; vectors.len()];
        layer.hash_many_vec(&vectors, &mut out);
        for (cell, v) in out.iter().zip(&vectors) {
            assert_eq!(*cell, layer.hash_vec(*v), "depth {depth}, vector");
            assert!(layer.contains(*cell));
        }
    }
}

#[test]
fn an_empty_batch_is_fine() {
    let mut out: Vec<u64> = Vec::new();
    nested::get(9).hash_many(&[], &mut out);
    nested::get(9).hash_many_vec(&[], &mut out);
    ring::get(9).hash_many(&[], &mut out);
    ring::get(9).hash_many_vec(&[], &mut out);
    assert!(out.is_empty());
}

/// The output slice must not silently be the wrong size — a short one would leave stale
/// cells behind, which is the kind of bug that shows up as a mis-solve much later.
#[test]
#[should_panic(expected = "output slice must match")]
fn rejects_a_mismatched_output_slice() {
    let mut out = [0u64; 2];
    nested::get(5).hash_many(&[(0.0, 0.0), (1.0, 0.5), (2.0, -0.3)], &mut out);
}

#[test]
#[should_panic(expected = "output slice must match")]
fn rejects_a_mismatched_output_slice_for_vectors() {
    let mut out = [0u64; 4];
    ring::get(5).hash_many_vec(&[[1.0, 0.0, 0.0]], &mut out);
}
