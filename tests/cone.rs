//! Cone search: the inclusiveness guarantee, the shape of the returned ranges, and the
//! degenerate cases.

mod common;

use common::*;
use realpix::nested;

/// The guarantee that matters to a sky solver: a source inside the cone is never missed.
#[test]
fn every_point_inside_the_cone_is_covered() {
    let mut rng = Rng::new(0xDEC0DE);
    for depth in [0, 2, 5, 8, 12] {
        let layer = nested::get(depth);
        for _ in 0..40 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.3 + 1e-4;
            let cells = layer.cone_coverage_cells(center, radius);

            for _ in 0..500 {
                let p = rng.next_in_cone(center, radius);
                assert!(ang_dist(p, center) <= radius + 1e-12);

                let cell = layer.hash_vec(p);
                assert!(
                    cells.binary_search(&cell).is_ok(),
                    "depth {depth}: point inside the cone landed in uncovered cell {cell}"
                );
            }
        }
    }
}

/// Every cell whose centre is inside the cone must be covered, at every depth.
#[test]
fn every_cell_centred_in_the_cone_is_covered() {
    let mut rng = Rng::new(1234);
    for depth in [1, 3, 5] {
        let layer = nested::get(depth);
        for _ in 0..30 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.8 + 0.01;
            let cells = layer.cone_coverage_cells(center, radius);
            for cell in layer.iter() {
                if ang_dist(layer.center_vec(cell), center) <= radius {
                    assert!(
                        cells.binary_search(&cell).is_ok(),
                        "depth {depth}: cell {cell} has its centre in the cone but is not covered"
                    );
                }
            }
        }
    }
}

#[test]
fn ranges_are_sorted_disjoint_and_in_bounds() {
    let mut rng = Rng::new(55);
    for depth in [0, 4, 9, 16] {
        let layer = nested::get(depth);
        for _ in 0..50 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 1.5;
            let ranges = layer.cone_coverage_ranges(center, radius);
            let mut previous_end = 0u64;
            for (i, r) in ranges.iter().enumerate() {
                assert!(r.start < r.end, "empty range");
                assert!(r.end <= layer.n_hash(), "range past the end of the layer");
                if i > 0 {
                    assert!(
                        r.start > previous_end,
                        "ranges must be sorted and non-adjacent: {previous_end} then {r:?}"
                    );
                }
                previous_end = r.end;
            }
        }
    }
}

#[test]
fn coverage_is_not_wastefully_large() {
    // Compare against the cells that provably intersect the cone (approximated by a fine
    // sub-sampling of every candidate cell).
    let mut rng = Rng::new(808);
    let mut worst = 0.0f64;
    for depth in [4, 6, 8] {
        let layer = nested::get(depth);
        for _ in 0..20 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.2 + 0.01;
            let covered = layer.cone_coverage_cells(center, radius);
            // A cell counts as touching if any of its corners, its centre, or any of its
            // 64 descendants three depths down falls inside the cone.
            let deep = nested::get(depth + 3);
            let touching = covered
                .iter()
                .filter(|c| {
                    ang_dist(layer.center_vec(**c), center) <= radius
                        || layer
                            .vertices(**c)
                            .iter()
                            .any(|v| ang_dist(*v, center) <= radius)
                        || layer
                            .children_range(**c, depth + 3)
                            .any(|d| ang_dist(deep.center_vec(d), center) <= radius)
                })
                .count();
            // Small discs are dominated by the boundary, so also allow a small absolute
            // slack; the ratio is only meaningful once the disc spans many cells.
            assert!(
                covered.len() as f64 <= touching as f64 * 1.1 + 8.0,
                "depth {depth}, radius {radius}: covered {} cells against {touching} touching",
                covered.len()
            );
            if touching >= 20 {
                worst = worst.max(covered.len() as f64 / touching as f64);
            }
        }
    }
    // Measured worst case at the time of writing: 1.21, on the smallest discs where the
    // boundary band dominates. Larger discs come in under 1.01.
    println!("worst over-inclusion ratio: {worst}");
    assert!(
        worst < 1.3,
        "cone coverage returns up to {worst}x the cells that touch it"
    );
}

#[test]
fn degenerate_radii_behave() {
    let layer = nested::get(6);
    let center = lonlat_to_vec(1.0, 0.4);

    // A zero radius still returns the cell containing the centre.
    let cells = layer.cone_coverage_cells(center, 0.0);
    assert!(
        cells.contains(&layer.hash_vec(center)),
        "zero radius lost its own cell"
    );

    // A radius covering the sphere returns everything, as a single range.
    let ranges = layer.cone_coverage_ranges(center, std::f64::consts::PI);
    assert_eq!(ranges, vec![0..layer.n_hash()]);
    let ranges = layer.cone_coverage_ranges(center, 4.0);
    assert_eq!(ranges, vec![0..layer.n_hash()]);

    // Nonsense inputs yield nothing rather than panicking.
    assert!(layer.cone_coverage_ranges(center, -1.0).is_empty());
    assert!(layer.cone_coverage_ranges(center, f64::NAN).is_empty());
    assert!(
        layer
            .cone_coverage_ranges([f64::NAN, 0.0, 0.0], 0.1)
            .is_empty()
    );
    // A zero axis has finite components but no direction. Normalising it would make every
    // pruning comparison NaN — which answers "no" to both reject and accept — and the
    // descent would then walk the whole layer instead of returning.
    assert!(layer.cone_coverage_ranges([0.0, 0.0, 0.0], 0.1).is_empty());
    assert!(
        nested::get(realpix::MAX_DEPTH)
            .cone_coverage_ranges([0.0, 0.0, 0.0], f64::MIN_POSITIVE)
            .is_empty()
    );
    // An axis so long its own dot product overflows is equally directionless.
    assert!(
        layer
            .cone_coverage_ranges([1e300, 1e300, 1e300], 0.1)
            .is_empty()
    );

    // The lon/lat entry point agrees with the vector one.
    let mut from_lonlat = Vec::new();
    layer.cone_coverage_lonlat(1.0, 0.4, 0.1, |r| from_lonlat.push(r));
    assert_eq!(from_lonlat, layer.cone_coverage_ranges(center, 0.1));
}

#[test]
fn cone_over_a_pole_is_covered() {
    // Cones spanning a pole are where naive longitude-range implementations break.
    let layer = nested::get(7);
    for pole in [1.0, -1.0] {
        let center = [0.0, 0.0, pole];
        let cells = layer.cone_coverage_cells(center, 0.2);
        let mut rng = Rng::new(3);
        for _ in 0..5000 {
            let lat = pole * (std::f64::consts::FRAC_PI_2 - rng.next_f64() * 0.2);
            let lon = rng.next_f64() * std::f64::consts::TAU;
            let p = lonlat_to_vec(lon, lat);
            if ang_dist(p, center) <= 0.2 {
                assert!(
                    cells.binary_search(&layer.hash_vec(p)).is_ok(),
                    "pole {pole} cone missed a point"
                );
            }
        }
        // All four base cells around the pole must be represented.
        let faces: std::collections::BTreeSet<u64> = cells
            .iter()
            .map(|c| nested::get(0).hash_vec(layer.center_vec(*c)))
            .collect();
        assert_eq!(
            faces.len(),
            4,
            "pole {pole}: expected all four polar base cells"
        );
    }
}

// ---------------------------------------------------------------------------------------
// RING scheme
// ---------------------------------------------------------------------------------------

use realpix::ring;

/// The same inclusiveness guarantee as the NESTED search: a source inside the cone is
/// never missed.
#[test]
fn ring_covers_every_point_inside_the_cone() {
    let mut rng = Rng::new(0x2149);
    for depth in [0, 2, 5, 8, 12] {
        let layer = ring::get(depth);
        for _ in 0..40 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.3 + 1e-4;
            let cells = layer.cone_coverage_cells(center, radius);

            for _ in 0..500 {
                let p = rng.next_in_cone(center, radius);
                assert!(
                    cells.binary_search(&layer.hash_vec(p)).is_ok(),
                    "depth {depth}: point inside the cone landed in an uncovered cell"
                );
            }
        }
    }
}

/// Every cell whose centre is inside the cone must be covered, at every depth.
#[test]
fn ring_covers_every_cell_centred_in_the_cone() {
    let mut rng = Rng::new(99);
    for depth in [1, 3, 5] {
        let layer = ring::get(depth);
        for _ in 0..30 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.8 + 0.01;
            let cells = layer.cone_coverage_cells(center, radius);
            for cell in layer.iter() {
                if ang_dist(layer.center_vec(cell), center) <= radius {
                    assert!(
                        cells.binary_search(&cell).is_ok(),
                        "depth {depth}: cell {cell} has its centre in the cone but is not covered"
                    );
                }
            }
        }
    }
}

/// The two schemes describe the same cells, so their coverages must denote the same set —
/// up to each search's own over-inclusion, which is what the containment direction tests.
#[test]
fn ring_and_nested_coverages_agree_on_the_cells_that_matter() {
    let mut rng = Rng::new(7788);
    for depth in [2, 4, 6] {
        let nest = nested::get(depth);
        let rings = ring::get(depth);
        for _ in 0..25 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.25 + 0.005;
            let covered = rings.cone_coverage_cells(center, radius);

            // Everything the NESTED search returns whose centre is genuinely in the cone
            // must be in the RING result too, translated across.
            for cell in nest.cone_coverage_cells(center, radius) {
                if ang_dist(nest.center_vec(cell), center) <= radius {
                    let as_ring = nest.to_ring(cell);
                    assert!(
                        covered.binary_search(&as_ring).is_ok(),
                        "depth {depth}: NESTED cell {cell} (RING {as_ring}) is missing"
                    );
                }
            }
        }
    }
}

#[test]
fn ring_ranges_are_sorted_disjoint_and_in_bounds() {
    let mut rng = Rng::new(4242);
    for depth in [0, 4, 9, 16] {
        let layer = ring::get(depth);
        for _ in 0..50 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 1.5;
            let ranges = layer.cone_coverage_ranges(center, radius);
            let mut previous_end = 0u64;
            for (i, r) in ranges.iter().enumerate() {
                assert!(r.start < r.end, "empty range");
                assert!(r.end <= layer.n_hash(), "range past the end of the layer");
                if i > 0 {
                    assert!(
                        r.start > previous_end,
                        "ranges must be sorted and non-adjacent: {previous_end} then {r:?}"
                    );
                }
                previous_end = r.end;
            }
        }
    }
}

/// A disc covers an arc of each ring it reaches, so the result should be about as many
/// cells as the cells that actually touch it — the search must not fall back to whole
/// rings or whole caps.
#[test]
fn ring_coverage_is_not_wastefully_large() {
    let mut rng = Rng::new(31337);
    let mut worst = 0.0f64;
    for depth in [4, 6, 8] {
        let layer = ring::get(depth);
        let nest = nested::get(depth);
        let deep = nested::get(depth + 3);
        for _ in 0..20 {
            let center = rng.next_vec();
            let radius = rng.next_f64() * 0.2 + 0.01;
            let covered = layer.cone_coverage_cells(center, radius);
            let touching = covered
                .iter()
                .filter(|c| {
                    let n = layer.to_nested(**c);
                    ang_dist(layer.center_vec(**c), center) <= radius
                        || layer
                            .vertices(**c)
                            .iter()
                            .any(|v| ang_dist(*v, center) <= radius)
                        || nest
                            .children_range(n, depth + 3)
                            .any(|d| ang_dist(deep.center_vec(d), center) <= radius)
                })
                .count();
            assert!(
                covered.len() as f64 <= touching as f64 * 1.5 + 8.0,
                "depth {depth}, radius {radius}: covered {} cells against {touching} touching",
                covered.len()
            );
            if touching >= 20 {
                worst = worst.max(covered.len() as f64 / touching as f64);
            }
        }
    }
    // The RING search does no boundary refinement — it takes every cell whose centre is
    // within `radius + rho` — so it over-includes more than the NESTED one (measured
    // worst case 1.33 against that search's 1.21).
    println!("worst RING over-inclusion ratio: {worst}");
    assert!(
        worst < 1.45,
        "ring cone coverage returns up to {worst}x the cells that touch it"
    );
}

#[test]
fn ring_degenerate_radii_behave() {
    let layer = ring::get(6);
    let center = lonlat_to_vec(1.0, 0.4);

    let cells = layer.cone_coverage_cells(center, 0.0);
    assert!(
        cells.contains(&layer.hash_vec(center)),
        "zero radius lost its own cell"
    );

    let ranges = layer.cone_coverage_ranges(center, std::f64::consts::PI);
    assert_eq!(ranges, vec![0..layer.n_hash()]);
    let ranges = layer.cone_coverage_ranges(center, 4.0);
    assert_eq!(ranges, vec![0..layer.n_hash()]);

    assert!(layer.cone_coverage_ranges(center, -1.0).is_empty());
    assert!(layer.cone_coverage_ranges(center, f64::NAN).is_empty());
    assert!(
        layer
            .cone_coverage_ranges([f64::NAN, 0.0, 0.0], 0.1)
            .is_empty()
    );
    // A zero-length axis has no direction to search around.
    assert!(layer.cone_coverage_ranges([0.0, 0.0, 0.0], 0.1).is_empty());

    let mut from_lonlat = Vec::new();
    layer.cone_coverage_lonlat(1.0, 0.4, 0.1, |r| from_lonlat.push(r));
    assert_eq!(from_lonlat, layer.cone_coverage_ranges(center, 0.1));
}

/// A cone over a pole wraps every ring it touches in longitude, and a cone straddling
/// longitude zero wraps some of them — the two cases the arc arithmetic has to get right.
#[test]
fn ring_cone_wraps_in_longitude() {
    let layer = ring::get(7);
    let mut rng = Rng::new(3);

    for pole in [1.0, -1.0] {
        let center = [0.0, 0.0, pole];
        let cells = layer.cone_coverage_cells(center, 0.2);
        for _ in 0..5000 {
            let lat = pole * (std::f64::consts::FRAC_PI_2 - rng.next_f64() * 0.2);
            let lon = rng.next_f64() * std::f64::consts::TAU;
            let p = lonlat_to_vec(lon, lat);
            if ang_dist(p, center) <= 0.2 {
                assert!(
                    cells.binary_search(&layer.hash_vec(p)).is_ok(),
                    "pole {pole} cone missed a point"
                );
            }
        }
    }

    // Centred on longitude zero, so roughly half of every arc is on either side of it.
    for lat in [0.0, 0.5, -1.2] {
        let center = lonlat_to_vec(0.0, lat);
        let cells = layer.cone_coverage_cells(center, 0.15);
        for _ in 0..3000 {
            let p = rng.next_in_cone(center, 0.15);
            assert!(
                cells.binary_search(&layer.hash_vec(p)).is_ok(),
                "cone across longitude zero at lat {lat} missed a point"
            );
        }
    }
}
