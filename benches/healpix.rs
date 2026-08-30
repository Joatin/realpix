//! Benchmarks for the whole public API.
//!
//! Every public function of the crate is benchmarked here, under an id of the form
//! `module::function`; `tests/bench_coverage.rs` checks that none has been missed. The
//! groups follow the shape of the API, and `solver` at the end is the curated hot-path
//! workload to watch for regressions.
//!
//! The `radec::*` benchmarks need the `latlong` feature:
//!
//! ```text
//! cargo bench --all-features
//! ```

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use realpix::moc::Moc;
use realpix::{Direction, MAX_DEPTH, nested, ring};

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

/// The depth at which the per-cell benchmarks run. Deep enough to be realistic for a
/// solver, shallow enough that a whole layer is still enumerable in the tests.
const D: u8 = 12;

// ------------------------------------------------------------------------ free functions

fn depth_helpers(c: &mut Criterion) {
    let mut g = c.benchmark_group("depth");
    g.bench_function("depth::nside", |b| {
        b.iter(|| black_box(realpix::nside(black_box(D))))
    });
    g.bench_function("depth::n_hash", |b| {
        b.iter(|| black_box(realpix::n_hash(black_box(D))))
    });
    g.bench_function("depth::depth_from_nside", |b| {
        b.iter(|| black_box(realpix::depth_from_nside(black_box(4096))))
    });
    g.finish();
}

fn uniq(c: &mut Criterion) {
    let mut g = c.benchmark_group("uniq");
    g.bench_function("uniq::to_uniq", |b| {
        b.iter(|| black_box(realpix::to_uniq(black_box(D), black_box(123_456))))
    });
    let u = realpix::to_uniq(D, 123_456);
    g.bench_function("uniq::from_uniq", |b| {
        b.iter(|| black_box(realpix::from_uniq(black_box(u))))
    });
    g.finish();
}

fn tangent(c: &mut Criterion) {
    let points = directions(1024);
    let mut g = c.benchmark_group("tangent");
    g.bench_function("tangent::lonlat_to_vec", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            black_box(realpix::lonlat_to_vec(lon, lat))
        })
    });
    g.bench_function("tangent::vec_to_lonlat", |b| {
        let mut i = 0;
        b.iter(|| {
            let (_, _, v) = points[i % points.len()];
            i += 1;
            black_box(realpix::vec_to_lonlat(v))
        })
    });
    g.bench_function("tangent::angular_distance", |b| {
        let mut i = 0;
        b.iter(|| {
            let a = points[i % points.len()].2;
            let b2 = points[(i + 1) % points.len()].2;
            i += 1;
            black_box(realpix::angular_distance(a, b2))
        })
    });
    g.bench_function("tangent::gnomonic_project", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            black_box(realpix::gnomonic_project(lon, lat, lon + 1e-4, lat + 1e-4))
        })
    });
    g.finish();
}

fn direction(c: &mut Criterion) {
    let mut g = c.benchmark_group("direction");
    // A `repr(u8)` cast: this measures the call overhead, not the operation.
    g.bench_function("xyf::index", |b| {
        b.iter(|| black_box(black_box(Direction::NE).index()))
    });
    g.finish();
}

// ------------------------------------------------------------------------------- nested

fn nested_layer(c: &mut Criterion) {
    let layer = nested::get(D);
    let mut g = c.benchmark_group("nested/layer");
    // These fold to a constant when the depth is known at compile time; the numbers here
    // are the cost when it is not.
    g.bench_function("nested::get", |b| {
        b.iter(|| black_box(nested::get(black_box(D))))
    });
    g.bench_function("nested::checked_get", |b| {
        b.iter(|| black_box(nested::checked_get(black_box(D))))
    });
    g.bench_function("nested::depth", |b| {
        b.iter(|| black_box(black_box(&layer).depth()))
    });
    g.bench_function("nested::nside", |b| {
        b.iter(|| black_box(black_box(&layer).nside()))
    });
    g.bench_function("nested::n_hash", |b| {
        b.iter(|| black_box(black_box(&layer).n_hash()))
    });
    g.bench_function("nested::cell_area", |b| {
        b.iter(|| black_box(black_box(&layer).cell_area()))
    });
    g.bench_function("nested::contains", |b| {
        b.iter(|| black_box(black_box(&layer).contains(black_box(12_345))))
    });
    g.bench_function("nested::iter", |b| {
        b.iter(|| black_box(black_box(&layer).iter()))
    });
    g.finish();
}

fn nested_hash(c: &mut Criterion) {
    let points = directions(1024);
    let pairs: Vec<(f64, f64)> = points.iter().map(|(a, b, _)| (*a, *b)).collect();
    let vecs: Vec<[f64; 3]> = points.iter().map(|(_, _, v)| *v).collect();
    let mut out = vec![0u64; points.len()];

    let mut g = c.benchmark_group("nested/hash");
    for depth in [8u8, D, 16, MAX_DEPTH] {
        let layer = nested::get(depth);
        g.bench_with_input(BenchmarkId::new("nested::hash", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let (lon, lat, _) = points[i % points.len()];
                i += 1;
                black_box(layer.hash(lon, lat))
            })
        });
        g.bench_with_input(
            BenchmarkId::new("nested::hash_vec", depth),
            &depth,
            |b, _| {
                let mut i = 0;
                b.iter(|| {
                    let (_, _, v) = points[i % points.len()];
                    i += 1;
                    black_box(layer.hash_vec(v))
                })
            },
        );
    }
    let layer = nested::get(D);
    g.bench_function("nested::hash_theta_phi", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            black_box(layer.hash_theta_phi(std::f64::consts::FRAC_PI_2 - lat, lon))
        })
    });
    g.bench_function("nested::checked_hash", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            black_box(layer.checked_hash(lon, lat))
        })
    });
    g.throughput(criterion::Throughput::Elements(points.len() as u64));
    g.bench_function("nested::hash_many", |b| {
        b.iter(|| {
            layer.hash_many(&pairs, &mut out);
            black_box(&out);
        })
    });
    g.bench_function("nested::hash_many_vec", |b| {
        b.iter(|| {
            layer.hash_many_vec(&vecs, &mut out);
            black_box(&out);
        })
    });
    // The loop the bulk calls are meant to match, for comparison.
    g.bench_function("nested::hash_many (as a loop)", |b| {
        b.iter(|| {
            for (o, (lon, lat)) in out.iter_mut().zip(pairs.iter()) {
                *o = layer.hash(*lon, *lat);
            }
            black_box(&out);
        })
    });
    g.finish();
}

fn nested_geometry(c: &mut Criterion) {
    let points = directions(1024);
    let mut g = c.benchmark_group("nested/geometry");
    for depth in [8u8, D, 16, MAX_DEPTH] {
        let layer = nested::get(depth);
        let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();
        g.bench_with_input(BenchmarkId::new("nested::center", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                i += 1;
                black_box(layer.center(cells[i % cells.len()]))
            })
        });
        g.bench_with_input(
            BenchmarkId::new("nested::center_vec", depth),
            &depth,
            |b, _| {
                let mut i = 0;
                b.iter(|| {
                    i += 1;
                    black_box(layer.center_vec(cells[i % cells.len()]))
                })
            },
        );
        g.bench_with_input(
            BenchmarkId::new("nested::vertices", depth),
            &depth,
            |b, _| {
                let mut i = 0;
                b.iter(|| {
                    i += 1;
                    black_box(layer.vertices(cells[i % cells.len()]))
                })
            },
        );
    }
    let layer = nested::get(D);
    let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();
    g.bench_function("nested::vertices_lonlat", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.vertices_lonlat(cells[i % cells.len()]))
        })
    });
    g.bench_function("nested::checked_center", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.checked_center(cells[i % cells.len()]))
        })
    });
    g.finish();
}

fn nested_neighbours(c: &mut Criterion) {
    let points = directions(1024);
    let layer = nested::get(D);
    let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();

    let mut g = c.benchmark_group("nested/neighbours");
    g.bench_function("nested::neighbours", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.neighbours(cells[i % cells.len()]))
        })
    });
    g.bench_function("nested::neighbour", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.neighbour(cells[i % cells.len()], Direction::N))
        })
    });
    for delta in [0u8, 2, 5] {
        g.bench_with_input(
            BenchmarkId::new("nested::external_edge", delta),
            &delta,
            |b, delta| {
                let mut i = 0;
                b.iter(|| {
                    i += 1;
                    let mut n = 0u64;
                    layer.external_edge(cells[i % cells.len()], D + delta, |_| n += 1);
                    black_box(n)
                })
            },
        );
        g.bench_with_input(
            BenchmarkId::new("nested::external_edge_cells", delta),
            &delta,
            |b, delta| {
                let mut i = 0;
                b.iter(|| {
                    i += 1;
                    black_box(layer.external_edge_cells(cells[i % cells.len()], D + delta))
                })
            },
        );
    }
    g.finish();
}

fn nested_hierarchy(c: &mut Criterion) {
    let points = directions(1024);
    let layer = nested::get(D);
    let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();

    let mut g = c.benchmark_group("nested/hierarchy");
    g.bench_function("nested::parent", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.parent(cells[i % cells.len()], 4))
        })
    });
    g.bench_function("nested::children", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.children(cells[i % cells.len()]))
        })
    });
    g.bench_function("nested::children_range", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.children_range(cells[i % cells.len()], D + 4))
        })
    });
    g.bench_function("nested::to_ring", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.to_ring(cells[i % cells.len()]))
        })
    });
    g.finish();
}

fn nested_cone(c: &mut Criterion) {
    let points = directions(64);
    let mut g = c.benchmark_group("nested/cone");
    for (depth, radius) in [(8u8, 0.05), (D, 0.01), (D, 0.05), (16, 0.005)] {
        let layer = nested::get(depth);
        g.bench_with_input(
            BenchmarkId::new(format!("nested::cone_coverage/depth{depth}"), radius),
            &radius,
            |b, radius| {
                let mut i = 0;
                b.iter(|| {
                    let (_, _, v) = points[i % points.len()];
                    i += 1;
                    let mut cells = 0u64;
                    layer.cone_coverage(v, *radius, |r| cells += r.end - r.start);
                    black_box(cells)
                })
            },
        );
    }
    let layer = nested::get(D);
    g.bench_function("nested::cone_coverage_lonlat", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            let mut cells = 0u64;
            layer.cone_coverage_lonlat(lon, lat, 0.02, |r| cells += r.end - r.start);
            black_box(cells)
        })
    });
    g.bench_function("nested::cone_coverage_ranges", |b| {
        let mut i = 0;
        b.iter(|| {
            let (_, _, v) = points[i % points.len()];
            i += 1;
            black_box(layer.cone_coverage_ranges(v, 0.02))
        })
    });
    g.bench_function("nested::cone_coverage_cells", |b| {
        let mut i = 0;
        b.iter(|| {
            let (_, _, v) = points[i % points.len()];
            i += 1;
            black_box(layer.cone_coverage_cells(v, 0.02))
        })
    });
    g.finish();
}

// --------------------------------------------------------------------------------- ring

fn ring_layer(c: &mut Criterion) {
    let layer = ring::get(D);
    let mut g = c.benchmark_group("ring/layer");
    g.bench_function("ring::get", |b| {
        b.iter(|| black_box(ring::get(black_box(D))))
    });
    g.bench_function("ring::checked_get", |b| {
        b.iter(|| black_box(ring::checked_get(black_box(D))))
    });
    g.bench_function("ring::depth", |b| {
        b.iter(|| black_box(black_box(&layer).depth()))
    });
    g.bench_function("ring::nside", |b| {
        b.iter(|| black_box(black_box(&layer).nside()))
    });
    g.bench_function("ring::n_hash", |b| {
        b.iter(|| black_box(black_box(&layer).n_hash()))
    });
    g.bench_function("ring::n_cap", |b| {
        b.iter(|| black_box(black_box(&layer).n_cap()))
    });
    g.bench_function("ring::contains", |b| {
        b.iter(|| black_box(black_box(&layer).contains(black_box(12_345))))
    });
    g.bench_function("ring::iter", |b| {
        b.iter(|| black_box(black_box(&layer).iter()))
    });
    g.finish();
}

fn ring_hash(c: &mut Criterion) {
    let points = directions(1024);
    let pairs: Vec<(f64, f64)> = points.iter().map(|(a, b, _)| (*a, *b)).collect();
    let vecs: Vec<[f64; 3]> = points.iter().map(|(_, _, v)| *v).collect();
    let mut out = vec![0u64; points.len()];

    let mut g = c.benchmark_group("ring/hash");
    for depth in [8u8, D, 16, MAX_DEPTH] {
        let layer = ring::get(depth);
        g.bench_with_input(BenchmarkId::new("ring::hash", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let (lon, lat, _) = points[i % points.len()];
                i += 1;
                black_box(layer.hash(lon, lat))
            })
        });
        g.bench_with_input(BenchmarkId::new("ring::hash_vec", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                let (_, _, v) = points[i % points.len()];
                i += 1;
                black_box(layer.hash_vec(v))
            })
        });
    }
    let layer = ring::get(D);
    g.bench_function("ring::hash_theta_phi", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            black_box(layer.hash_theta_phi(std::f64::consts::FRAC_PI_2 - lat, lon))
        })
    });
    g.bench_function("ring::checked_hash", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            black_box(layer.checked_hash(lon, lat))
        })
    });
    g.throughput(criterion::Throughput::Elements(points.len() as u64));
    g.bench_function("ring::hash_many", |b| {
        b.iter(|| {
            layer.hash_many(&pairs, &mut out);
            black_box(&out);
        })
    });
    g.bench_function("ring::hash_many_vec", |b| {
        b.iter(|| {
            layer.hash_many_vec(&vecs, &mut out);
            black_box(&out);
        })
    });
    g.finish();
}

fn ring_geometry(c: &mut Criterion) {
    let points = directions(1024);
    let mut g = c.benchmark_group("ring/geometry");
    for depth in [8u8, D, 16, MAX_DEPTH] {
        let layer = ring::get(depth);
        let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();
        g.bench_with_input(BenchmarkId::new("ring::center", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                i += 1;
                black_box(layer.center(cells[i % cells.len()]))
            })
        });
        g.bench_with_input(
            BenchmarkId::new("ring::center_vec", depth),
            &depth,
            |b, _| {
                let mut i = 0;
                b.iter(|| {
                    i += 1;
                    black_box(layer.center_vec(cells[i % cells.len()]))
                })
            },
        );
    }
    let layer = ring::get(D);
    let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();
    g.bench_function("ring::vertices", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.vertices(cells[i % cells.len()]))
        })
    });
    g.bench_function("ring::checked_center", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.checked_center(cells[i % cells.len()]))
        })
    });
    g.bench_function("ring::to_nested", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.to_nested(cells[i % cells.len()]))
        })
    });
    g.finish();
}

fn ring_neighbours(c: &mut Criterion) {
    let points = directions(1024);
    let layer = ring::get(D);
    let cells: Vec<u64> = points.iter().map(|(_, _, v)| layer.hash_vec(*v)).collect();

    let mut g = c.benchmark_group("ring/neighbours");
    g.bench_function("ring::neighbours", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.neighbours(cells[i % cells.len()]))
        })
    });
    g.bench_function("ring::neighbour", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.neighbour(cells[i % cells.len()], Direction::N))
        })
    });
    g.finish();
}

fn ring_cone(c: &mut Criterion) {
    let points = directions(64);
    let mut g = c.benchmark_group("ring/cone");
    for (depth, radius) in [(8u8, 0.05), (D, 0.01), (D, 0.05), (16, 0.005)] {
        let layer = ring::get(depth);
        g.bench_with_input(
            BenchmarkId::new(format!("ring::cone_coverage/depth{depth}"), radius),
            &radius,
            |b, radius| {
                let mut i = 0;
                b.iter(|| {
                    let (_, _, v) = points[i % points.len()];
                    i += 1;
                    let mut cells = 0u64;
                    layer.cone_coverage(v, *radius, |r| cells += r.end - r.start);
                    black_box(cells)
                })
            },
        );
    }
    let layer = ring::get(D);
    g.bench_function("ring::cone_coverage_lonlat", |b| {
        let mut i = 0;
        b.iter(|| {
            let (lon, lat, _) = points[i % points.len()];
            i += 1;
            let mut cells = 0u64;
            layer.cone_coverage_lonlat(lon, lat, 0.02, |r| cells += r.end - r.start);
            black_box(cells)
        })
    });
    g.bench_function("ring::cone_coverage_ranges", |b| {
        let mut i = 0;
        b.iter(|| {
            let (_, _, v) = points[i % points.len()];
            i += 1;
            black_box(layer.cone_coverage_ranges(v, 0.02))
        })
    });
    g.bench_function("ring::cone_coverage_cells", |b| {
        let mut i = 0;
        b.iter(|| {
            let (_, _, v) = points[i % points.len()];
            i += 1;
            black_box(layer.cone_coverage_cells(v, 0.02))
        })
    });
    g.finish();
}

// ---------------------------------------------------------------------------------- moc

fn moc_build(c: &mut Criterion) {
    let points = directions(64);
    let layer = nested::get(D);
    let cone_cells = layer.cone_coverage_cells(points[0].2, 0.05);
    let cone_ranges = layer.cone_coverage_ranges(points[0].2, 0.05);
    let uniqs: Vec<u64> = Moc::from_cone(D, points[0].2, 0.05).uniq_cells().collect();

    let mut g = c.benchmark_group("moc/build");
    g.bench_function("moc::new", |b| b.iter(|| black_box(Moc::new())));
    g.bench_function("moc::all_sky", |b| b.iter(|| black_box(Moc::all_sky())));
    g.bench_function("moc::from_cells", |b| {
        b.iter(|| black_box(Moc::from_cells(D, cone_cells.iter().copied())))
    });
    g.bench_function("moc::from_ranges", |b| {
        b.iter(|| black_box(Moc::from_ranges(D, cone_ranges.iter().cloned())))
    });
    g.bench_function("moc::from_uniq_cells", |b| {
        b.iter(|| black_box(Moc::from_uniq_cells(uniqs.iter().copied())))
    });
    for depth in [8u8, D] {
        g.bench_with_input(BenchmarkId::new("moc::from_cone", depth), &depth, |b, _| {
            let mut i = 0;
            b.iter(|| {
                i += 1;
                black_box(Moc::from_cone(depth, points[i % points.len()].2, 0.05))
            })
        });
    }
    g.finish();
}

fn moc_query(c: &mut Criterion) {
    let points = directions(64);
    let moc = Moc::from_cone(D, points[0].2, 0.05);
    let mut g = c.benchmark_group("moc/query");
    g.bench_function("moc::is_empty", |b| {
        b.iter(|| black_box(black_box(&moc).is_empty()))
    });
    g.bench_function("moc::area", |b| {
        b.iter(|| black_box(black_box(&moc).area()))
    });
    g.bench_function("moc::sky_fraction", |b| {
        b.iter(|| black_box(black_box(&moc).sky_fraction()))
    });
    g.bench_function("moc::contains", |b| {
        b.iter(|| black_box(moc.contains(D, black_box(12_345_678))))
    });
    g.bench_function("moc::contains_vec", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(moc.contains_vec(points[i % points.len()].2))
        })
    });
    g.bench_function("moc::contains_lonlat", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            let (lon, lat, _) = points[i % points.len()];
            black_box(moc.contains_lonlat(lon, lat))
        })
    });
    g.finish();
}

fn moc_setops(c: &mut Criterion) {
    let points = directions(64);
    let mut g = c.benchmark_group("moc/setops");
    for depth in [8u8, D] {
        // Two overlapping fields, as a solver would build them.
        let a = Moc::from_cone(depth, points[0].2, 0.05);
        let b2 = Moc::from_cone(depth, points[1].2, 0.05);
        g.bench_with_input(BenchmarkId::new("moc::union", depth), &depth, |b, _| {
            b.iter(|| black_box(&a | &b2))
        });
        g.bench_with_input(
            BenchmarkId::new("moc::intersection", depth),
            &depth,
            |b, _| b.iter(|| black_box(&a & &b2)),
        );
        g.bench_with_input(
            BenchmarkId::new("moc::difference", depth),
            &depth,
            |b, _| b.iter(|| black_box(&a - &b2)),
        );
        g.bench_with_input(
            BenchmarkId::new("moc::symmetric_difference", depth),
            &depth,
            |b, _| b.iter(|| black_box(&a ^ &b2)),
        );
        g.bench_with_input(
            BenchmarkId::new("moc::complement", depth),
            &depth,
            |b, _| b.iter(|| black_box(!&a)),
        );
    }

    // Accumulating a survey. Both spellings are here because the point of `union_all` is
    // that folding `|` over the same sequence is quadratic.
    for count in [50usize, 400] {
        let frames: Vec<Moc> = (0..count)
            .map(|i| Moc::from_cone(9, points[i % points.len()].2, 0.02))
            .collect();
        g.bench_with_input(BenchmarkId::new("moc::union_all", count), &count, |b, _| {
            b.iter(|| black_box(Moc::union_all(frames.iter()).deep_ranges().len()))
        });
        g.bench_with_input(
            BenchmarkId::new("moc::union_all (as a fold)", count),
            &count,
            |b, _| {
                b.iter(|| {
                    let m = frames.iter().fold(Moc::new(), |acc, f| &acc | f);
                    black_box(m.deep_ranges().len())
                })
            },
        );
    }
    g.finish();
}

fn moc_export(c: &mut Criterion) {
    let points = directions(64);
    let moc = Moc::from_cone(D, points[0].2, 0.05);
    let mut g = c.benchmark_group("moc/export");
    g.bench_function("moc::ranges_at", |b| {
        b.iter(|| black_box(moc.ranges_at(black_box(D))))
    });
    g.bench_function("moc::cells", |b| b.iter(|| black_box(moc.cells().count())));
    g.bench_function("moc::uniq_cells", |b| {
        b.iter(|| black_box(moc.uniq_cells().count()))
    });
    g.bench_function("moc::deep_ranges", |b| {
        b.iter(|| black_box(black_box(&moc).deep_ranges().len()))
    });
    g.finish();
}

// -------------------------------------------------------------------------------- radec

#[cfg(feature = "latlong")]
fn radec(c: &mut Criterion) {
    use latlong::{Declination, RaDec, RightAscension};

    let points = directions(1024);
    let coords: Vec<RaDec<f64>> = points
        .iter()
        .map(|(lon, lat, _)| RaDec {
            ra: RightAscension::from_radians(*lon),
            dec: Declination::from_radians(*lat),
        })
        .collect();
    let layer = nested::get(D);
    let cells: Vec<u64> = coords.iter().map(|c| layer.hash_ra_dec(c)).collect();

    let mut g = c.benchmark_group("radec");
    g.bench_function("radec::hash_ra_dec", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.hash_ra_dec(&coords[i % coords.len()]))
        })
    });
    g.bench_function("radec::center_ra_dec", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(layer.center_ra_dec::<f64>(cells[i % cells.len()]))
        })
    });
    g.bench_function("radec::cone_coverage_ra_dec", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            let mut n = 0u64;
            layer.cone_coverage_ra_dec(&coords[i % coords.len()], 0.02, |r| n += r.end - r.start);
            black_box(n)
        })
    });
    g.bench_function("radec::project_ra_dec", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            let j = i % coords.len();
            black_box(layer.project_ra_dec::<f64>(cells[j], &coords[j]))
        })
    });
    g.bench_function("radec::dec_to_theta", |b| {
        let mut i = 0;
        b.iter(|| {
            i += 1;
            black_box(realpix::radec::dec_to_theta(coords[i % coords.len()].dec))
        })
    });
    g.finish();
}

#[cfg(not(feature = "latlong"))]
fn radec(_: &mut Criterion) {}

// ------------------------------------------------------------------------------- solver

/// A solver-shaped workload: one cone per frame, then a cell lookup per detected source.
fn solver_workload(c: &mut Criterion) {
    let layer = nested::get(D);
    let sources = directions(200);
    c.bench_function("solver/frame", |b| {
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

criterion_group!(
    benches,
    depth_helpers,
    uniq,
    tangent,
    direction,
    nested_layer,
    nested_hash,
    nested_geometry,
    nested_neighbours,
    nested_hierarchy,
    nested_cone,
    ring_layer,
    ring_hash,
    ring_geometry,
    ring_neighbours,
    ring_cone,
    moc_build,
    moc_query,
    moc_setops,
    moc_export,
    radec,
    solver_workload
);
criterion_main!(benches);
