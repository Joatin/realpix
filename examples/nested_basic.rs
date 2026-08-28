//! Basic NESTED usage: positions to cells, cells back to positions, and the hierarchy.

use realpix::nested;

// Binding the layer to a `const` lets the compiler fold every depth-derived constant.
const LAYER: nested::Layer = nested::get(12);

fn main() {
    // Betelgeuse, in radians.
    let (ra, dec) = (1.549_729_f64, 0.129_277_f64);

    let cell = LAYER.hash(ra, dec);
    let (lon, lat) = LAYER.center(cell);
    println!("depth {}  nside {}", LAYER.depth(), LAYER.nside());
    println!("cell {cell} of {}", LAYER.n_hash());
    println!("centre  ra {lon:.6}  dec {lat:.6}");
    println!(
        "cell size {:.3} arcsec",
        LAYER.cell_area().sqrt().to_degrees() * 3600.0
    );

    // The same cell in the RING scheme.
    println!("ring index {}", LAYER.to_ring(cell));

    // Neighbours, for a local search.
    let neighbours: Vec<u64> = LAYER.neighbours(cell).into_iter().flatten().collect();
    println!("neighbours {neighbours:?}");

    // Where this cell sits in a coarser index, and what it covers in a finer one.
    println!("ancestor at depth 5: {}", LAYER.parent(cell, 5));
    println!(
        "descendants at depth 16: {:?}",
        LAYER.children_range(cell, 16)
    );
}
