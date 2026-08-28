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
