//! Keeps the README's quick-start example compiling and correct.

use realpix::nested;

const LAYER: nested::Layer = nested::get(12);

#[test]
fn readme_quick_start() {
    let cell = LAYER.hash(1.549_729, 0.129_277);
    let (lon, lat) = LAYER.center(cell);
    let corners = LAYER.vertices(cell);
    let neighbours = LAYER.neighbours(cell);

    let r = LAYER.to_ring(cell);
    assert_eq!(cell, realpix::ring::get(12).to_nested(r));

    let mut covered = 0u64;
    LAYER.cone_coverage(
        realpix::lonlat_to_vec(1.55, 0.13),
        5f64.to_radians(),
        |range| covered += range.end - range.start,
    );

    assert_eq!(LAYER.hash(lon, lat), cell);
    assert_eq!(corners.len(), 4);
    assert_eq!(neighbours.iter().flatten().count(), 8);
    assert!(covered > 1000);
}

/// The README's coverage-map example.
#[test]
fn readme_coverage_maps() {
    use realpix::moc::Moc;

    let first = Moc::from_cone(10, realpix::lonlat_to_vec(1.00, 0.5), 0.02);
    let second = Moc::from_cone(10, realpix::lonlat_to_vec(1.01, 0.5), 0.02);

    let observed = &first | &second;
    let overlap = &first & &second;
    let missed = &second - &first;
    let elsewhere = !&observed;

    assert!(observed.area() > 0.0);
    assert!(observed.contains_lonlat(1.0, 0.5));

    assert!(!overlap.is_empty(), "the two fields do overlap");
    assert_eq!(&overlap | &missed, second);
    assert_eq!(&observed | &elsewhere, Moc::all_sky());
    assert!(observed.ranges_at(12).len() > 1);
}

/// The README's external-edge example.
#[test]
fn readme_expanding_a_search() {
    let layer = realpix::nested::get(10);
    let cell = layer.hash(1.549_729, 0.129_277);

    let adjacent = layer.external_edge_cells(cell, 10);
    let finer = layer.external_edge_cells(cell, 12);
    assert_eq!(adjacent.len(), 8);
    assert_eq!(finer.len(), 4 * 4 + 4);
}

/// The README's bulk-hashing example.
#[test]
fn readme_hashing_in_bulk() {
    let layer = realpix::nested::get(12);
    let sources = [(1.549_729, 0.129_277), (1.372_198, -0.143_146)];
    let mut cells = [0u64; 2];
    layer.hash_many(&sources, &mut cells);
    assert_eq!(cells[0], layer.hash(sources[0].0, sources[0].1));
}
