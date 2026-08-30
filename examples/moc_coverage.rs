//! Multi-order coverage maps: accumulating the sky a survey has reached, and asking what
//! is left.
//!
//! A plate solver works frame by frame. Each solved frame is a disc of sky; the union of
//! them is where the survey has been, and the parts of a target region that union misses
//! are what still needs pointing at.

use realpix::moc::Moc;

const DEPTH: u8 = 10;

fn main() {
    // Six solved frames, each a 1.5 degree field.
    let radius = 1.5f64.to_radians();
    let frames = [
        (1.55, 0.13),
        (1.58, 0.15),
        (1.61, 0.13),
        (1.55, 0.19),
        (1.61, 0.19),
        (2.90, -0.40),
    ];

    // Collecting unions them in one pass. Folding `|` over the sequence would copy the
    // whole accumulated coverage on every frame, which a real survey would feel.
    let observed: Moc = frames
        .iter()
        .map(|(lon, lat)| Moc::from_cone(DEPTH, realpix::lonlat_to_vec(*lon, *lat), radius))
        .collect();

    let frame_area: f64 = frames
        .iter()
        .map(|(lon, lat)| Moc::from_cone(DEPTH, realpix::lonlat_to_vec(*lon, *lat), radius).area())
        .sum();

    println!("{} frames", frames.len());
    println!(
        "  observed:  {:.4} sr in {} cells at mixed depths ({:.1}% of the sky)",
        observed.area(),
        observed.cells().count(),
        observed.sky_fraction() * 100.0
    );
    println!(
        "  overlap:   {:.4} sr lost to frames covering the same sky",
        frame_area - observed.area()
    );

    // The region we actually care about: 5 degrees around Betelgeuse.
    let target = Moc::from_cone(
        DEPTH,
        realpix::lonlat_to_vec(1.549_729, 0.129_277),
        5.0f64.to_radians(),
    );

    let done = &target & &observed;
    let todo = &target - &observed;
    println!("\ntarget region: {:.4} sr", target.area());
    println!(
        "  covered:   {:.4} sr ({:.1}%)",
        done.area(),
        done.area() / target.area() * 100.0
    );
    println!("  remaining: {:.4} sr", todo.area());

    // What is left is an arbitrary shape, but it still slices a sorted catalogue directly.
    let ranges = todo.ranges_at(DEPTH);
    println!(
        "  the remainder is {} index ranges at depth {DEPTH}, or {} multi-order cells",
        ranges.len(),
        todo.cells().count()
    );

    // Point tests come straight off the coverage.
    for (name, lon, lat) in [
        ("Betelgeuse", 1.549_729, 0.129_277),
        ("Bellatrix", 1.418_772, 0.110_823),
    ] {
        let where_it_is = if observed.contains_lonlat(lon, lat) {
            "already observed"
        } else if target.contains_lonlat(lon, lat) {
            "in the target region, not yet observed"
        } else {
            "outside the target region"
        };
        println!("  {name}: {where_it_is}");
    }
}
