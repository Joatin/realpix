//! Multi-order coverage maps: the set algebra against a brute-force reference, the
//! multi-order decomposition, and the ways a coverage is projected back to one depth.

mod common;

use common::*;
use realpix::moc::Moc;
use realpix::{MAX_DEPTH, n_hash, nested};
use std::collections::BTreeSet;

/// The depth the brute-force reference works at. Small enough to enumerate every cell.
const D: u8 = 4;

/// Every cell of the coverage at depth `D`, for comparison against a `BTreeSet`.
fn cells_at(moc: &Moc, depth: u8) -> BTreeSet<u64> {
    let mut out = BTreeSet::new();
    for (d, cell) in moc.cells() {
        assert!(
            d <= depth,
            "coverage holds a cell at depth {d}, deeper than {depth}"
        );
        for c in nested::get(d).children_range(cell, depth) {
            out.insert(c);
        }
    }
    out
}

fn random_set(rng: &mut Rng, depth: u8, count: usize) -> BTreeSet<u64> {
    let limit = n_hash(depth);
    (0..count).map(|_| rng.next_u64() % limit).collect()
}

/// Union, intersection, difference, symmetric difference and complement must each agree
/// with the same operation done cell by cell on a plain set.
#[test]
fn set_algebra_matches_a_brute_force_reference() {
    let mut rng = Rng::new(0xA11CE);
    for _ in 0..60 {
        let (n, m) = (rng.next_u64() % 300, rng.next_u64() % 300);
        let left = random_set(&mut rng, D, 1 + n as usize);
        let right = random_set(&mut rng, D, 1 + m as usize);
        let a = Moc::from_cells(D, left.iter().copied());
        let b = Moc::from_cells(D, right.iter().copied());

        assert_eq!(cells_at(&a, D), left, "round trip through the coverage");
        assert_eq!(
            cells_at(&(&a | &b), D),
            left.union(&right).copied().collect(),
            "union"
        );
        assert_eq!(
            cells_at(&(&a & &b), D),
            left.intersection(&right).copied().collect(),
            "intersection"
        );
        assert_eq!(
            cells_at(&(&a - &b), D),
            left.difference(&right).copied().collect(),
            "difference"
        );
        assert_eq!(
            cells_at(&(&a ^ &b), D),
            left.symmetric_difference(&right).copied().collect(),
            "symmetric difference"
        );
        assert_eq!(
            cells_at(&!&a, D),
            (0..n_hash(D)).filter(|c| !left.contains(c)).collect(),
            "complement"
        );
    }
}

/// The set laws the representation is supposed to make exact.
#[test]
fn set_laws_hold() {
    let mut rng = Rng::new(77);
    let empty = Moc::new();
    let all = Moc::all_sky();
    for _ in 0..30 {
        let a = Moc::from_cone(7, rng.next_vec(), rng.next_f64() * 0.3 + 0.01);
        let b = Moc::from_cone(9, rng.next_vec(), rng.next_f64() * 0.3 + 0.01);
        let c = Moc::from_cells(5, random_set(&mut rng, 5, 40));

        assert_eq!(&a | &a, a, "union is idempotent");
        assert_eq!(&a & &a, a, "intersection is idempotent");
        assert_eq!(&a | &b, &b | &a, "union commutes");
        assert_eq!(&a & &b, &b & &a, "intersection commutes");
        assert_eq!(&(&a | &b) | &c, &a | &(&b | &c), "union associates");
        assert_eq!(&(&a & &b) & &c, &a & &(&b & &c), "intersection associates");
        assert_eq!(&a | &empty, a, "empty is the union identity");
        assert_eq!(&a & &empty, empty, "empty annihilates intersection");
        assert_eq!(&a & &all, a, "all sky is the intersection identity");
        assert_eq!(&a | &all, all, "all sky absorbs union");
        assert_eq!(!&!&a, a, "complement is an involution");
        assert_eq!(
            &a - &b,
            &a & &!&b,
            "difference is intersection with complement"
        );
        assert_eq!(&a ^ &b, &(&a - &b) | &(&b - &a), "symmetric difference");
        // De Morgan.
        assert_eq!(!&(&a | &b), &!&a & &!&b);
        assert_eq!(!&(&a & &b), &!&a | &!&b);
    }
}

/// Two coverages of the same region must be equal whatever depths they were built from —
/// that is the point of normalising to one depth internally.
#[test]
fn equality_ignores_the_depths_a_coverage_was_built_from() {
    // One base cell, described four ways.
    let by_base = Moc::from_cells(0, [3]);
    let by_children = Moc::from_cells(1, 12..16);
    let by_grandchildren = Moc::from_cells(2, 48..64);
    let by_range = Moc::from_ranges(1, std::iter::once(12..16));
    assert_eq!(by_base, by_children);
    assert_eq!(by_base, by_grandchildren);
    assert_eq!(by_base, by_range);
    assert_eq!(by_base, Moc::from_uniq_cells([realpix::to_uniq(0, 3)]));

    // And it decomposes back to the single shallow cell.
    assert_eq!(by_grandchildren.cells().collect::<Vec<_>>(), [(0, 3)]);

    // A mixed-depth build is the union of its parts.
    let mixed = Moc::from_uniq_cells([realpix::to_uniq(0, 3), realpix::to_uniq(2, 100)]);
    assert_eq!(mixed, &by_base | &Moc::from_cells(2, [100]));

    assert_eq!(Moc::all_sky(), Moc::from_cells(0, 0..12));
    assert!(Moc::new().is_empty());
    assert_eq!(Moc::new(), Moc::default());
}

/// `cells` must tile the coverage exactly, and with the largest cells that fit: no cell it
/// emits may have a parent that is itself wholly covered.
#[test]
fn the_multi_order_decomposition_is_exact_and_maximal() {
    let mut rng = Rng::new(2024);
    for _ in 0..40 {
        let moc = Moc::from_cone(8, rng.next_vec(), rng.next_f64() * 0.2 + 0.005);

        let mut covered = 0u64;
        let mut previous: Option<(u8, u64)> = None;
        for (depth, cell) in moc.cells() {
            assert!(depth <= MAX_DEPTH);
            assert!(cell < n_hash(depth), "cell out of range for its depth");
            assert!(
                moc.contains(depth, cell),
                "emitted a cell not in the coverage"
            );
            if depth > 0 {
                assert!(
                    !moc.contains(depth - 1, cell >> 2),
                    "cell ({depth}, {cell}) should have been emitted as its parent"
                );
            }
            // Increasing in position, so the tiles come out in sky order.
            let deep = cell << (2 * (MAX_DEPTH - depth) as u32);
            if let Some((pd, pc)) = previous {
                let previous_deep = pc << (2 * (MAX_DEPTH - pd) as u32);
                assert!(deep > previous_deep, "cells must be emitted in order");
            }
            previous = Some((depth, cell));
            covered += 1u64 << (2 * (MAX_DEPTH - depth) as u32);
        }

        let expected: u64 = moc.deep_ranges().iter().map(|r| r.end - r.start).sum();
        assert_eq!(
            covered, expected,
            "the tiling must cover exactly the coverage"
        );

        // NUNIQ round trip.
        assert_eq!(Moc::from_uniq_cells(moc.uniq_cells()), moc);
    }
}

/// `contains` asks whether a whole cell is inside, and must agree with the point tests.
#[test]
fn containment_agrees_across_the_entry_points() {
    let mut rng = Rng::new(4);
    let layer = nested::get(6);
    for _ in 0..20 {
        let center = rng.next_vec();
        let radius = rng.next_f64() * 0.3 + 0.02;
        let moc = Moc::from_cone(6, center, radius);
        let inside: BTreeSet<u64> = moc.ranges_at(6).into_iter().flatten().collect();

        for cell in layer.iter() {
            assert_eq!(
                moc.contains(6, cell),
                inside.contains(&cell),
                "cell {cell} at the depth the coverage was built at"
            );
            // A cell in the coverage contains its own centre; one outside does not.
            let (lon, lat) = layer.center(cell);
            assert_eq!(moc.contains_lonlat(lon, lat), inside.contains(&cell));
            assert_eq!(
                moc.contains_vec(layer.center_vec(cell)),
                inside.contains(&cell)
            );
        }

        // A cell only partly covered is not contained, but its covered children are.
        let shallow = nested::get(4);
        for cell in shallow.iter() {
            if moc.contains(4, cell) {
                for child in shallow.children_range(cell, 6) {
                    assert!(
                        moc.contains(6, child),
                        "a contained cell must contain its children"
                    );
                }
            }
        }
    }
}

/// `ranges_at` is exact at the depth the coverage was built at, and a superset above it.
#[test]
fn ranges_at_rounds_outward() {
    let mut rng = Rng::new(31);
    for _ in 0..30 {
        let built_at = 8;
        let moc = Moc::from_cone(built_at, rng.next_vec(), rng.next_f64() * 0.15 + 0.01);
        let exact: BTreeSet<u64> = moc.ranges_at(built_at).into_iter().flatten().collect();

        // Exact at the build depth: every cell is fully in the coverage.
        for cell in &exact {
            assert!(moc.contains(built_at, *cell));
        }

        for depth in 0..=built_at {
            let ranges = moc.ranges_at(depth);
            let mut previous_end = 0u64;
            for (i, r) in ranges.iter().enumerate() {
                assert!(r.start < r.end, "empty range");
                assert!(r.end <= n_hash(depth), "range past the end of the layer");
                if i > 0 {
                    assert!(r.start > previous_end, "sorted, disjoint, non-adjacent");
                }
                previous_end = r.end;
            }
            // Superset: every covered cell at the build depth has its ancestor included.
            let shallow: BTreeSet<u64> = ranges.into_iter().flatten().collect();
            for cell in &exact {
                let ancestor = nested::get(built_at).parent(*cell, depth);
                assert!(
                    shallow.contains(&ancestor),
                    "depth {depth} dropped the ancestor of cell {cell}"
                );
            }
        }

        // Below the build depth it stays exact, just expressed in smaller cells.
        let deeper: BTreeSet<u64> = moc.ranges_at(10).into_iter().flatten().collect();
        for cell in &exact {
            for child in nested::get(built_at).children_range(*cell, 10) {
                assert!(deeper.contains(&child));
            }
        }
    }
}

/// A cone coverage must hold the cone, and measure the area the cone covers.
#[test]
fn a_cone_coverage_holds_the_cone() {
    let mut rng = Rng::new(606);
    for depth in [5u8, 8] {
        for _ in 0..15 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.2 + 0.01;
            let moc = Moc::from_cone(depth, center, radius);

            // The convenience constructor must agree with the long way round.
            let layer = nested::get(depth);
            assert_eq!(
                moc,
                Moc::from_ranges(depth, layer.cone_coverage_ranges(center, radius))
            );
            assert_eq!(
                moc,
                Moc::from_cells(depth, layer.cone_coverage_cells(center, radius))
            );

            // Every point of the cone is in the coverage.
            for _ in 0..300 {
                let p = random_point_in_cone(&mut rng, center, radius);
                assert!(
                    moc.contains_vec(p),
                    "a point inside the cone is not covered"
                );
            }

            // The coverage is inclusive, so its area is at least the cone's.
            let cone_area = 2.0 * std::f64::consts::PI * (1.0 - radius.cos());
            assert!(
                moc.area() >= cone_area,
                "area {} is below the cone's own {cone_area}",
                moc.area()
            );
            // And `area` must agree exactly with the cells the coverage actually holds.
            let held = layer.cone_coverage_cells(center, radius).len() as f64;
            assert!(
                (moc.area() - held * layer.cell_area()).abs() < 1e-12,
                "area {} against {held} cells of {}",
                moc.area(),
                layer.cell_area()
            );
        }
    }
}

#[test]
fn area_and_sky_fraction_are_consistent() {
    let all = Moc::all_sky();
    assert!((all.sky_fraction() - 1.0).abs() < 1e-15);
    assert!((all.area() - 4.0 * std::f64::consts::PI).abs() < 1e-12);
    assert_eq!(Moc::new().area(), 0.0);
    assert_eq!(Moc::new().sky_fraction(), 0.0);
    assert!(!Moc::new().contains_lonlat(1.0, 0.0));
    assert!(all.contains_lonlat(1.0, 0.0));

    // One base cell is a twelfth of the sky, however it is expressed.
    let base = Moc::from_cells(0, [7]);
    assert!((base.sky_fraction() - 1.0 / 12.0).abs() < 1e-15);
    assert!((Moc::from_cells(3, 448..452).sky_fraction() - 4.0 / n_hash(3) as f64).abs() < 1e-15);
}

/// The deepest layer is where the representation could overflow, so exercise it directly.
#[test]
fn the_deepest_depth_round_trips() {
    let last = n_hash(MAX_DEPTH) - 1;
    let moc = Moc::from_cells(MAX_DEPTH, [0, 1, last]);
    assert_eq!(
        moc.cells().collect::<Vec<_>>(),
        [(MAX_DEPTH, 0), (MAX_DEPTH, 1), (MAX_DEPTH, last)]
    );
    assert!(moc.contains(MAX_DEPTH, 0));
    assert!(!moc.contains(MAX_DEPTH, 2));
    assert!(
        !moc.contains(0, 0),
        "a base cell is not covered by two of its cells"
    );
    assert_eq!(moc.ranges_at(MAX_DEPTH), vec![0..2, last..last + 1]);
    assert_eq!(Moc::from_uniq_cells(moc.uniq_cells()), moc);
    assert_eq!(&Moc::all_sky() - &moc, moc.complement());
}

#[test]
fn empty_and_degenerate_inputs_behave() {
    assert!(Moc::from_cells(5, []).is_empty());
    assert!(Moc::from_ranges(5, []).is_empty());
    assert!(
        Moc::from_ranges(5, std::iter::once(10..10)).is_empty(),
        "an empty range covers nothing"
    );
    assert!(Moc::from_uniq_cells([]).is_empty());
    assert!(Moc::from_cone(5, lonlat_to_vec(1.0, 0.0), -1.0).is_empty());
    assert_eq!(Moc::new().complement(), Moc::all_sky());
    assert_eq!(Moc::all_sky().complement(), Moc::new());
    assert_eq!(Moc::new().cells().count(), 0);
    assert!(Moc::new().ranges_at(9).is_empty());
    // Duplicates and overlaps normalise away.
    assert_eq!(Moc::from_cells(3, [5, 5, 5]), Moc::from_cells(3, [5]));
    assert_eq!(
        Moc::from_ranges(3, [0..10, 5..20, 20..25]),
        Moc::from_ranges(3, std::iter::once(0..25))
    );
}

#[test]
#[should_panic(expected = "out of range for depth")]
fn rejects_a_cell_beyond_its_depth() {
    Moc::from_cells(1, [48]);
}

#[test]
#[should_panic(expected = "out of range for depth")]
fn rejects_a_range_beyond_its_depth() {
    Moc::from_ranges(1, std::iter::once(40..49));
}

fn random_point_in_cone(rng: &mut Rng, center: [f64; 3], radius: f64) -> [f64; 3] {
    let cos_r = radius.cos();
    let cos_t = 1.0 - rng.next_f64() * (1.0 - cos_r);
    let sin_t = (1.0 - cos_t * cos_t).sqrt();
    let az = rng.next_f64() * std::f64::consts::TAU;
    let a = if center[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let u = norm(cross(a, center));
    let v = cross(center, u);
    norm([
        center[0] * cos_t + (u[0] * az.cos() + v[0] * az.sin()) * sin_t,
        center[1] * cos_t + (u[1] * az.cos() + v[1] * az.sin()) * sin_t,
        center[2] * cos_t + (u[2] * az.cos() + v[2] * az.sin()) * sin_t,
    ])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}
