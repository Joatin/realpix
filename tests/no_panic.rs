//! The functions documented as *not* panicking must hold that line under hostile input.
//!
//! The `checked_*` family exists so a caller can hand the crate a value from outside the
//! program — a file header, a network message, a command line — without risking a panic.
//! The cone searches make the same promise about coordinates. Nothing here asserts a
//! particular answer; the assertion is that these calls return at all.

mod common;

use common::*;
use realpix::{MAX_DEPTH, moc::Moc, nested, ring};

/// Values chosen to break arithmetic: non-finite, out of range, denormal, extreme.
const NASTY: [f64; 14] = [
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    f64::MAX,
    f64::MIN,
    f64::MIN_POSITIVE,
    -f64::MIN_POSITIVE,
    0.0,
    -0.0,
    1e300,
    -1e300,
    1e-300,
    std::f64::consts::PI,
    -std::f64::consts::PI,
];

#[test]
fn checked_constructors_return_errors_rather_than_panicking() {
    for depth in 0..=255u8 {
        assert_eq!(
            nested::checked_get(depth).is_ok(),
            depth <= MAX_DEPTH,
            "nested depth {depth}"
        );
        assert_eq!(
            ring::checked_get(depth).is_ok(),
            depth <= MAX_DEPTH,
            "ring depth {depth}"
        );
    }
}

#[test]
fn checked_hash_and_center_never_panic() {
    for depth in [0u8, 1, 12, MAX_DEPTH] {
        let n = nested::get(depth);
        let r = ring::get(depth);
        for lon in NASTY {
            for lat in NASTY {
                // Only the assertion that these return; the values are checked elsewhere.
                let _ = n.checked_hash(lon, lat);
                let _ = r.checked_hash(lon, lat);
            }
        }
        for cell in [0u64, 1, u64::MAX, u64::MAX / 2, n.n_hash(), n.n_hash() - 1] {
            let ok = n.checked_center(cell).is_ok();
            assert_eq!(ok, n.contains(cell), "nested depth {depth}, cell {cell}");
            let ok = r.checked_center(cell).is_ok();
            assert_eq!(ok, r.contains(cell), "ring depth {depth}, cell {cell}");
        }
    }
}

/// `hash` itself is total: the projection maps any finite input into the grid, so a
/// position that is merely absurd rather than invalid still lands somewhere.
#[test]
fn hash_never_panics_on_a_finite_position() {
    for depth in [0u8, 7, MAX_DEPTH] {
        let n = nested::get(depth);
        let r = ring::get(depth);
        for lon in NASTY.iter().filter(|v| v.is_finite()) {
            for lat in [-1.5f64, -0.5, 0.0, 0.5, 1.5] {
                assert!(n.contains(n.hash(*lon, lat)), "nested {lon}, {lat}");
                assert!(r.contains(r.hash(*lon, lat)), "ring {lon}, {lat}");
            }
        }
        for v in [[0.0, 0.0, 0.0], [1e300, 1e300, 1e300], [1e-300, 0.0, 0.0]] {
            assert!(n.contains(n.hash_vec(v)), "nested {v:?}");
            assert!(r.contains(r.hash_vec(v)), "ring {v:?}");
        }
    }
}

#[test]
fn cone_searches_never_panic() {
    for depth in [0u8, 6, MAX_DEPTH] {
        let n = nested::get(depth);
        let r = ring::get(depth);
        for radius in NASTY {
            for center in [
                [0.0, 0.0, 0.0],
                [f64::NAN, 0.0, 0.0],
                [f64::INFINITY, 1.0, 1.0],
                [1e300, 1e300, 1e300],
                [0.0, 0.0, 1.0],
            ] {
                let mut count = 0u64;
                n.cone_coverage(center, radius, |x| count += x.end - x.start);
                assert!(count <= n.n_hash());
                let mut count = 0u64;
                r.cone_coverage(center, radius, |x| count += x.end - x.start);
                assert!(count <= r.n_hash());
            }
        }
    }
}

/// A coverage built from real data must survive every query without panicking, including
/// the ones that take a depth or a cell.
#[test]
fn moc_queries_never_panic() {
    let mut rng = Rng::new(4242);
    let moc = Moc::from_cone(8, rng.next_vec(), 0.1);
    for depth in 0..=MAX_DEPTH {
        // The deepest valid cell at each depth, and the shallowest.
        for cell in [0u64, realpix::n_hash(depth) - 1] {
            let _ = moc.contains(depth, cell);
        }
        let _ = moc.ranges_at(depth);
    }
    for lon in NASTY {
        for lat in NASTY {
            let _ = moc.contains_lonlat(lon, lat);
        }
    }
    for v in [[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0], [1e300, 0.0, 0.0]] {
        let _ = moc.contains_vec(v);
    }
    // The multi-order export of any coverage is always re-importable.
    assert_eq!(Moc::from_uniq_cells(moc.uniq_cells()), moc);
}

/// `angular_distance` and the coordinate conversions are handed vectors from real data,
/// but must not fall over on the degenerate ones either.
#[test]
fn tangent_helpers_never_panic() {
    for v in [
        [0.0, 0.0, 0.0],
        [f64::NAN, 0.0, 0.0],
        [1e300, 1e300, 1e300],
        [0.0, 0.0, 1.0],
    ] {
        let _ = realpix::vec_to_lonlat(v);
        let _ = realpix::angular_distance(v, [1.0, 0.0, 0.0]);
    }
    for lon in NASTY {
        for lat in NASTY {
            let _ = realpix::lonlat_to_vec(lon, lat);
            let _ = realpix::gnomonic_project(lon, lat, lat, lon);
        }
    }
}
