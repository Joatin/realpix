//! Verifies the geometric bound the cone search relies on.
//!
//! `nested::MAX_CENTER_TO_VERTEX[d]` must be an upper bound on the distance from any
//! cell's centre to any point on that cell's boundary. Under-estimating it would make the
//! cone search silently drop cells, so this is checked exhaustively where that is
//! affordable and by sampling everywhere else.

mod common;

use common::*;
use realpix::{MAX_DEPTH, nested};

const EXHAUSTIVE_DEPTH: u8 = 6;

fn max_center_to_vertex(layer: &nested::Layer, cell: u64) -> f64 {
    let c = layer.center_vec(cell);
    layer
        .vertices(cell)
        .into_iter()
        .map(|v| ang_dist(c, v))
        .fold(0.0f64, f64::max)
}

#[test]
fn center_to_vertex_bound_holds_exhaustively() {
    for depth in 0..=EXHAUSTIVE_DEPTH {
        let layer = nested::get(depth);
        let bound = nested::MAX_CENTER_TO_VERTEX[depth as usize];
        let mut worst = 0.0f64;
        for cell in layer.iter() {
            worst = worst.max(max_center_to_vertex(&layer, cell));
        }
        assert!(
            worst <= bound,
            "depth {depth}: measured {worst} exceeds the bound {bound}"
        );
        // The bound should also not be wildly loose, or cone searches get sloppy.
        assert!(
            worst > 0.9 * bound,
            "depth {depth}: bound {bound} is much larger than the measured {worst}"
        );
    }
}

#[test]
fn center_to_vertex_bound_holds_at_every_depth() {
    let mut rng = Rng::new(0x5EED);
    for depth in (EXHAUSTIVE_DEPTH + 1)..=MAX_DEPTH {
        let layer = nested::get(depth);
        let bound = nested::MAX_CENTER_TO_VERTEX[depth as usize];
        let mut worst = 0.0f64;
        for _ in 0..20_000 {
            worst = worst.max(max_center_to_vertex(
                &layer,
                rng.next_u64() % layer.n_hash(),
            ));
        }
        // Also probe the cells that are extremal at shallow depths: the ones touching a
        // pole and the ones on the belt/cap transition.
        for cell in [0, 1, 2, 3, layer.n_hash() - 1, layer.n_hash() / 2] {
            worst = worst.max(max_center_to_vertex(&layer, cell));
        }
        assert!(
            worst <= bound,
            "depth {depth}: measured {worst} exceeds the bound {bound}"
        );
    }
}

#[test]
fn cell_boundary_never_leaves_the_bound() {
    // The bound is derived from cell corners; confirm that no point along an edge is
    // farther from the centre than the corners are.
    for depth in 0..=4 {
        let layer = nested::get(depth);
        let bound = nested::MAX_CENTER_TO_VERTEX[depth as usize];
        for cell in layer.iter() {
            let c = layer.center_vec(cell);
            let vertices = layer.vertices(cell);
            for k in 0..4 {
                let (a, b) = (vertices[k], vertices[(k + 1) % 4]);
                for s in 1..32 {
                    let t = s as f64 / 32.0;
                    let p = [
                        a[0] + (b[0] - a[0]) * t,
                        a[1] + (b[1] - a[1]) * t,
                        a[2] + (b[2] - a[2]) * t,
                    ];
                    let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    let d = ang_dist(c, [p[0] / n, p[1] / n, p[2] / n]);
                    assert!(
                        d <= bound,
                        "depth {depth}, cell {cell}: edge point at {d} rad"
                    );
                }
            }
        }
    }
}

#[test]
fn cell_area_matches_the_sphere() {
    for depth in [0, 5, 12, 29] {
        let layer = nested::get(depth);
        let total = layer.cell_area() * layer.n_hash() as f64;
        assert!(
            (total - 4.0 * std::f64::consts::PI).abs() < 1e-9,
            "depth {depth}: cell areas sum to {total}"
        );
    }
}
