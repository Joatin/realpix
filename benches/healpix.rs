//! Hot-path benchmarks: the operations a sky solver runs per source, per frame.

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use realpix::{nested, ring};

/// Deterministic pseudo-random directions, so runs are comparable.
fn directions(n: usize) -> Vec<(f64, f64, [f64; 3])> {
    let mut state = 0x1234_5678_9ABC_DEF0u64;
    let mut next = move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    };
    (0..n)
        .map(|_| {
            let lon = next() * std::f64::consts::TAU;
            let lat = (2.0 * next() - 1.0).asin();
            (lon, lat, realpix::lonlat_to_vec(lon, lat))
        })
        .collect()
}

fn conversions(c: &mut Criterion) {
    let points = directions(1024);
    let mut group = c.benchmark_group("position_to_cell");
    for depth in [8u8, 12, 16, 29] {
        let n = nested::get(depth);
        let r = ring::get(depth);
        group.bench_with_input(BenchmarkId::new("nested_hash", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let (lon, lat, _) = points[i % points.len()];
                i += 1;
                black_box(n.hash(lon, lat))
            })
        });
        group.bench_with_input(
            BenchmarkId::new("nested_hash_vec", depth),
            &depth,
            |b, _| {
                let mut i = 0;
                b.iter(|| {
                    let (_, _, v) = points[i % points.len()];
                    i += 1;
                    black_box(n.hash_vec(v))
                })
            },
        );
        group.bench_with_input(BenchmarkId::new("ring_hash", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let (lon, lat, _) = points[i % points.len()];
                i += 1;
                black_box(r.hash(lon, lat))
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("cell_to_position");
    for depth in [8u8, 12, 16, 29] {
        let n = nested::get(depth);
        let cells: Vec<u64> = points
            .iter()
            .map(|(lon, lat, _)| n.hash(*lon, *lat))
            .collect();
        group.bench_with_input(BenchmarkId::new("center", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let c = cells[i % cells.len()];
                i += 1;
                black_box(n.center(c))
            })
        });
        group.bench_with_input(BenchmarkId::new("center_vec", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let c = cells[i % cells.len()];
                i += 1;
                black_box(n.center_vec(c))
            })
        });
        group.bench_with_input(BenchmarkId::new("to_ring", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let c = cells[i % cells.len()];
                i += 1;
                black_box(n.to_ring(c))
            })
        });
        group.bench_with_input(BenchmarkId::new("neighbours", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let c = cells[i % cells.len()];
                i += 1;
                black_box(n.neighbours(c))
            })
        });
    }
    group.finish();
}

fn cone(c: &mut Criterion) {
    let points = directions(64);
    let mut group = c.benchmark_group("cone_coverage");
    for (depth, radius) in [(8u8, 0.05), (12, 0.01), (12, 0.05), (16, 0.005)] {
        let layer = nested::get(depth);
        group.bench_with_input(
            BenchmarkId::new(format!("depth{depth}"), radius),
            &radius,
            |b, radius| {
                let mut i = 0;
                b.iter(|| {
                    let (_, _, v) = points[i % points.len()];
                    i += 1;
                    let mut ranges = 0u64;
                    layer.cone_coverage(v, *radius, |r| ranges += r.end - r.start);
                    black_box(ranges)
                })
            },
        );
    }
    group.finish();
}

/// A solver-shaped workload: one cone per frame, then a cell lookup per detected source.
fn solver_workload(c: &mut Criterion) {
    let layer = nested::get(12);
    let sources = directions(200);
    c.bench_function("solver_frame", |b| {
        b.iter(|| {
            let mut cells = 0u64;
            layer.cone_coverage(sources[0].2, 0.02, |r| cells += r.end - r.start);
            for (_, _, v) in &sources {
                cells ^= layer.hash_vec(*v);
            }
            black_box(cells)
        })
    });
}

criterion_group!(benches, conversions, cone, solver_workload);
criterion_main!(benches);
