//! Comparison against golden vectors produced by `healpy`, i.e. by the reference
//! HEALPix C++ implementation. See `tools/gen_golden.py`.

mod common;

use common::*;
use realpix::{nested, ring};

/// Cell indices must match the reference exactly, for both schemes and both entry points.
#[test]
fn ang2pix_matches_reference() {
    let csv = Csv::load("ang2pix.csv");
    assert_eq!(csv.header[0], "depth");
    let (mut checked, mut stable_checked) = (0, 0);
    for row in &csv.rows {
        let depth = u8_at(row, 0);
        let (lon, lat, theta) = (f64_at(row, 1), f64_at(row, 2), f64_at(row, 3));
        let (nest, rng) = (u64_at(row, 4), u64_at(row, 5));
        let v = [f64_at(row, 6), f64_at(row, 7), f64_at(row, 8)];
        let (vec_nest, vec_ring) = (u64_at(row, 9), u64_at(row, 10));
        let stable = u8_at(row, 11) == 1;

        let n = nested::get(depth);
        let r = ring::get(depth);

        // The theta/phi and unit-vector entry points take exactly the inputs the reference
        // was given, so they must agree on every sample, boundaries included.
        assert_eq!(
            n.hash_theta_phi(theta, lon),
            nest,
            "nested hash_theta_phi at depth {depth}, theta {theta}, phi {lon}"
        );
        assert_eq!(
            r.hash_theta_phi(theta, lon),
            rng,
            "ring hash_theta_phi at depth {depth}, theta {theta}, phi {lon}"
        );
        assert_eq!(
            n.hash_vec(v),
            vec_nest,
            "nested hash_vec at depth {depth}, v {v:?}"
        );
        assert_eq!(
            r.hash_vec(v),
            vec_ring,
            "ring hash_vec at depth {depth}, v {v:?}"
        );
        checked += 1;

        // `hash` derives z from sin(lat) rather than cos(theta); the two differ by an ulp,
        // so samples sitting on a cell boundary are excluded (see tools/gen_golden.py).
        if stable {
            assert_eq!(
                n.hash(lon, lat),
                nest,
                "nested hash at depth {depth}, lon {lon}, lat {lat}"
            );
            assert_eq!(
                r.hash(lon, lat),
                rng,
                "ring hash at depth {depth}, lon {lon}, lat {lat}"
            );
            stable_checked += 1;
        }
    }
    assert!(
        checked > 4000,
        "golden file looks truncated: {checked} rows"
    );
    assert!(
        stable_checked > 3000,
        "too few stable samples: {stable_checked}"
    );
}

/// Cell centres must match the reference to well below floating point noise.
#[test]
fn pix2ang_matches_reference() {
    let csv = Csv::load("pix2ang.csv");
    for row in &csv.rows {
        let depth = u8_at(row, 0);
        let cell = u64_at(row, 1);

        let expected = lonlat_to_vec(f64_at(row, 2), f64_at(row, 3));
        let (lon, lat) = nested::get(depth).center(cell);
        let d = ang_dist(lonlat_to_vec(lon, lat), expected);
        assert!(
            d < 1e-12,
            "nested centre at depth {depth}, cell {cell}: off by {d} rad"
        );
        // The unit-vector path must agree with the (lon, lat) path.
        let dv = ang_dist(nested::get(depth).center_vec(cell), expected);
        assert!(
            dv < 1e-12,
            "nested centre_vec at depth {depth}, cell {cell}: off by {dv} rad"
        );

        let expected = lonlat_to_vec(f64_at(row, 4), f64_at(row, 5));
        let (lon, lat) = ring::get(depth).center(cell);
        let d = ang_dist(lonlat_to_vec(lon, lat), expected);
        assert!(
            d < 1e-12,
            "ring centre at depth {depth}, cell {cell}: off by {d} rad"
        );
        let dv = ang_dist(ring::get(depth).center_vec(cell), expected);
        assert!(
            dv < 1e-12,
            "ring centre_vec at depth {depth}, cell {cell}: off by {dv} rad"
        );
    }
}

#[test]
fn scheme_conversion_matches_reference() {
    let csv = Csv::load("nest2ring.csv");
    for row in &csv.rows {
        let depth = u8_at(row, 0);
        let (nest, rng) = (u64_at(row, 1), u64_at(row, 2));
        assert_eq!(
            nested::get(depth).to_ring(nest),
            rng,
            "nest2ring at depth {depth}, cell {nest}"
        );
        assert_eq!(
            ring::get(depth).to_nested(rng),
            nest,
            "ring2nest at depth {depth}, cell {rng}"
        );
    }
}

#[test]
fn neighbours_match_reference() {
    let csv = Csv::load("neighbours.csv");
    for row in &csv.rows {
        let depth = u8_at(row, 0);
        let cell = u64_at(row, 1);
        let got = nested::get(depth).neighbours(cell);
        for (k, got) in got.iter().enumerate() {
            let expected = i64_at(row, 2 + k);
            let expected = if expected < 0 {
                None
            } else {
                Some(expected as u64)
            };
            assert_eq!(
                *got, expected,
                "neighbour {k} of cell {cell} at depth {depth}"
            );
        }
    }
}

#[test]
fn boundaries_match_reference() {
    let csv = Csv::load("boundaries.csv");
    for row in &csv.rows {
        let depth = u8_at(row, 0);
        let cell = u64_at(row, 1);
        let got = nested::get(depth).vertices(cell);
        for (k, got) in got.iter().enumerate() {
            let expected = [
                f64_at(row, 2 + 3 * k),
                f64_at(row, 3 + 3 * k),
                f64_at(row, 4 + 3 * k),
            ];
            let d = ang_dist(*got, expected);
            assert!(
                d < 1e-12,
                "vertex {k} of cell {cell} at depth {depth}: off by {d} rad"
            );
        }
    }
}

/// Our cone search is inclusive, so it must be a superset of healpy's exact disc.
#[test]
fn cone_coverage_contains_reference_disc() {
    let csv = Csv::load("query_disc.csv");
    for row in &csv.rows {
        let depth = u8_at(row, 0);
        let (lon, lat, radius) = (f64_at(row, 1), f64_at(row, 2), f64_at(row, 3));
        let expected: Vec<u64> = row[4]
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();

        let layer = nested::get(depth);
        let got = layer.cone_coverage_cells(lonlat_to_vec(lon, lat), radius);
        for cell in &expected {
            assert!(
                got.binary_search(cell).is_ok(),
                "cone at depth {depth}, ({lon}, {lat}) r={radius} is missing cell {cell}"
            );
        }
        // Sanity: we should not be wildly over-inclusive either.
        assert!(
            got.len() <= 8 * expected.len().max(4),
            "cone at depth {depth} returned {} cells for an exact {} — too loose",
            got.len(),
            expected.len()
        );
    }
}
