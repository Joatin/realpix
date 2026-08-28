#![allow(dead_code)]
//! Shared helpers for the integration tests: golden-file loading and small vector maths.

use std::fs;
use std::path::PathBuf;

/// A parsed golden CSV: the header names plus the rows, split on commas.
pub struct Csv {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl Csv {
    pub fn load(name: &str) -> Csv {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("tests/data");
        path.push(name);
        let text = fs::read_to_string(&path).unwrap_or_else(|e| {
            panic!(
                "cannot read {}: {e}\nrun tools/gen_golden.py",
                path.display()
            )
        });
        let mut lines = text.lines();
        let header = lines
            .next()
            .expect("empty golden file")
            .split(',')
            .map(str::to_owned)
            .collect();
        let rows = lines
            .filter(|l| !l.is_empty())
            .map(|l| l.split(',').map(str::to_owned).collect())
            .collect();
        Csv { header, rows }
    }
}

pub fn f64_at(row: &[String], i: usize) -> f64 {
    row[i].parse().unwrap()
}

pub fn u64_at(row: &[String], i: usize) -> u64 {
    row[i].parse().unwrap()
}

pub fn i64_at(row: &[String], i: usize) -> i64 {
    row[i].parse().unwrap()
}

pub fn u8_at(row: &[String], i: usize) -> u8 {
    row[i].parse().unwrap()
}

/// Angular distance in radians between two unit vectors, well conditioned at small angles.
pub fn ang_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2])
        .sqrt()
        .atan2(dot)
}

pub fn lonlat_to_vec(lon: f64, lat: f64) -> [f64; 3] {
    [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()]
}

/// A tiny deterministic PRNG (splitmix64), so tests never depend on a rand version.
pub struct Rng(pub u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// A direction uniformly distributed over the sphere, as `(lon, lat)`.
    pub fn next_lonlat(&mut self) -> (f64, f64) {
        let lon = self.next_f64() * std::f64::consts::TAU;
        let lat = (2.0 * self.next_f64() - 1.0).asin();
        (lon, lat)
    }

    pub fn next_vec(&mut self) -> [f64; 3] {
        let (lon, lat) = self.next_lonlat();
        lonlat_to_vec(lon, lat)
    }
}
