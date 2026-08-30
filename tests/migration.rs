//! Compiles every "0.2" entry in the CHANGELOG's 0.1 → 0.2 migration table.
//!
//! Migration advice that does not compile is worse than none, and a table in a markdown
//! file cannot be checked by the compiler on its own. This test is that check: if a symbol
//! named in the table is renamed or removed, this stops building.

use realpix::{depth_from_nside, nested, ring};

#[test]
fn the_changelog_migration_table_compiles() {
    // ConstHealpix::<32>::new()  ->  nested::get(5)
    const LAYER: nested::Layer = nested::get(5);
    assert_eq!(LAYER.nside(), 32);

    // DynamicHealpix::new(nside)?  ->  checked_get(depth_from_nside(nside)?)?
    let dynamic = nested::checked_get(depth_from_nside(32).unwrap()).unwrap();
    assert_eq!(dynamic, LAYER);

    // Pixel<N> -> u64, and the two marker types -> the two modules.
    let cell: u64 = LAYER.hash(1.0, 0.5);
    let _: u64 = ring::get(5).hash(1.0, 0.5);

    // .face_resolution() -> .nside(); .pixels_per_face() -> .n_hash() / 12;
    // .total_pixels() -> .n_hash()
    assert_eq!(LAYER.nside(), 32);
    assert_eq!(LAYER.n_hash() / 12, 32 * 32);
    assert_eq!(LAYER.n_hash(), 12 * 32 * 32);

    // .angle_to_pixel(theta, phi) -> .hash_theta_phi(theta, phi) or .hash(lon, lat)
    let theta = std::f64::consts::FRAC_PI_2 - 0.5;
    assert_eq!(LAYER.hash_theta_phi(theta, 1.0), LAYER.hash(1.0, 0.5));

    // .pixel_to_angle(pixel) -> .center(cell), with the documented convention change:
    // lon is phi, but lat is pi/2 - theta, and the tuple order is (lon, lat).
    let (lon, lat) = LAYER.center(cell);
    assert!((0.0..std::f64::consts::TAU).contains(&lon));
    assert!(lat.abs() <= std::f64::consts::FRAC_PI_2);
    assert_eq!(LAYER.hash(lon, lat), cell);

    // .iter_pixels() -> .iter()
    assert_eq!(LAYER.iter().count() as u64, LAYER.n_hash());

    // gnomonic_project now takes radians rather than RaDec.
    let projected = realpix::gnomonic_project(lon, lat, lon + 1e-4, lat);
    assert!(projected.is_some());

    // The depth limit moved: 0.1 capped out near depth 14 on a u32 pixel count.
    assert!(nested::get(realpix::MAX_DEPTH).n_hash() > u64::from(u32::MAX));
}

#[cfg(feature = "latlong")]
#[test]
fn the_latlong_migration_entries_compile() {
    use latlong::{Declination, RaDec, RightAscension};

    let layer = nested::get(10);
    let rd = RaDec {
        ra: RightAscension::from_radians(1.0_f64),
        dec: Declination::from_radians(0.5_f64),
    };

    // .ra_dec_to_pixel(&rd) -> .hash_ra_dec(&rd)
    let cell = layer.hash_ra_dec(&rd);
    // .pixel_to_ra_dec(pixel)? -> .center_ra_dec(cell)
    let back: RaDec<f64> = layer.center_ra_dec(cell);
    assert_eq!(layer.hash_ra_dec(&back), cell);
    // .project_ra_dec(pixel, &rd) -> .project_ra_dec(cell, &rd)
    assert!(layer.project_ra_dec::<f64>(cell, &rd).is_some());
}

/// The correctness claim the CHANGELOG leads with: 0.1.x indices cannot be converted, and
/// this release is pinned to `healpy`. Guard the property that made the old numbering
/// unusable — that spatial neighbours are adjacent in the numbering's own hierarchy.
#[test]
fn nested_indices_are_hierarchical_as_the_reference_requires() {
    let deep = nested::get(10);
    let shallow = nested::get(6);
    for cell in shallow.iter().step_by(97) {
        // Every descendant of a cell shares its prefix; this is what 0.1.x did not give.
        let range = shallow.children_range(cell, 10);
        for descendant in range.clone().step_by(37) {
            assert_eq!(deep.parent(descendant, 6), cell);
        }
    }
}
