//! The optional `latlong` integration.
#![cfg(feature = "latlong")]

use latlong::{Declination, RaDec, RightAscension};
use realpix::nested;

fn ra_dec(ra: f64, dec: f64) -> RaDec<f64> {
    RaDec {
        ra: RightAscension::from_radians(ra),
        dec: Declination::from_radians(dec),
    }
}

#[test]
fn ra_dec_agrees_with_the_plain_api() {
    let layer = nested::get(10);
    for (ra, dec) in [(0.0, 0.0), (1.5, 0.3), (4.9, -1.2), (0.1, 1.5), (6.2, -1.5)] {
        let position = ra_dec(ra, dec);
        assert_eq!(layer.hash_ra_dec(&position), layer.hash(ra, dec));

        let cell = layer.hash_ra_dec(&position);
        let center: RaDec<f64> = layer.center_ra_dec(cell);
        let (lon, lat) = layer.center(cell);
        assert!((center.ra.radians() - lon).abs() < 1e-15);
        assert!((center.dec.radians() - lat).abs() < 1e-15);
    }
}

#[test]
fn negative_right_ascension_wraps() {
    let layer = nested::get(8);
    let wrapped = ra_dec(-0.5, 0.2);
    let plain = ra_dec(-0.5 + std::f64::consts::TAU, 0.2);
    assert_eq!(layer.hash_ra_dec(&wrapped), layer.hash_ra_dec(&plain));
}

#[test]
fn projection_is_relative_to_the_cell_center() {
    let layer = nested::get(8);
    let position = ra_dec(1.5, 0.3);
    let cell = layer.hash_ra_dec(&position);

    // A position inside the cell projects to a small offset from its centre.
    let projected = layer
        .project_ra_dec(cell, &position)
        .expect("inside the cell");
    let cell_size = layer.cell_area().sqrt();
    assert!(projected.x.abs() < cell_size && projected.y.abs() < cell_size);

    // The cell centre itself projects to the origin.
    let center: RaDec<f64> = layer.center_ra_dec(cell);
    let origin = layer.project_ra_dec(cell, &center).unwrap();
    assert!(origin.x.abs() < 1e-12 && origin.y.abs() < 1e-12);

    // A position in a different cell is rejected.
    let elsewhere = ra_dec(3.0, -0.5);
    assert!(layer.project_ra_dec(cell, &elsewhere).is_none());
}

#[test]
fn cone_coverage_accepts_ra_dec() {
    let layer = nested::get(8);
    let position = ra_dec(1.5, 0.3);
    let mut from_ra_dec = Vec::new();
    layer.cone_coverage_ra_dec(&position, 0.05, |r| from_ra_dec.push(r));
    let from_vec = layer.cone_coverage_ranges(realpix::lonlat_to_vec(1.5, 0.3), 0.05);
    assert_eq!(from_ra_dec, from_vec);
}
