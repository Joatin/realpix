//! Cone search: turning a field of view into index ranges over a sorted catalogue.

use realpix::nested;

const LAYER: nested::Layer = nested::get(10);

/// Stand-in for a star catalogue sorted by NESTED cell index at `LAYER`'s depth.
struct Catalogue {
    entries: Vec<(u64, &'static str)>,
}

impl Catalogue {
    fn in_range(&self, range: std::ops::Range<u64>) -> &[(u64, &'static str)] {
        let start = self.entries.partition_point(|(c, _)| *c < range.start);
        let end = self.entries.partition_point(|(c, _)| *c < range.end);
        &self.entries[start..end]
    }
}

fn main() {
    let stars = [
        ("Betelgeuse", 1.549_729, 0.129_277),
        ("Rigel", 1.372_198, -0.143_146),
        ("Bellatrix", 1.418_772, 0.110_823),
        ("Sirius", 1.767_793, -0.291_752),
        ("Vega", 4.873_563, 0.676_902),
    ];
    let mut entries: Vec<(u64, &'static str)> = stars
        .iter()
        .map(|(name, ra, dec)| (LAYER.hash(*ra, *dec), *name))
        .collect();
    entries.sort();
    let catalogue = Catalogue { entries };

    // Everything within 5 degrees of Betelgeuse.
    let radius = 5.0f64.to_radians();
    let center = realpix::lonlat_to_vec(1.549_729, 0.129_277);

    let mut ranges = 0;
    let mut found = Vec::new();
    LAYER.cone_coverage(center, radius, |range| {
        ranges += 1;
        for (cell, name) in catalogue.in_range(range) {
            // The coverage is inclusive, so confirm the exact distance.
            let d = angular_distance(center, LAYER.center_vec(*cell));
            if d <= radius + 0.01 {
                found.push(*name);
            }
        }
    });

    println!("{ranges} index ranges scanned");
    println!("within 5 deg of Betelgeuse: {found:?}");
}

fn angular_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2])
        .sqrt()
        .atan2(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
}
