[![Cargo Build & Test](https://github.com/Joatin/realpix/actions/workflows/ci.yaml/badge.svg)](https://github.com/Joatin/realpix/actions/workflows/ci.yaml)

# REALPix

**`REALPix`** is a pure-Rust implementation of **HEALPix** — the Hierarchical Equal Area
isoLatitude Pixelisation of the sphere — in both the **NESTED** and **RING** numbering
schemes.

Cell indices are **bit-identical to the reference C++ implementation** (and therefore to
`healpy`), which is enforced by golden test vectors committed to the repository. The crate
is built for the inner loop of astrometry and plate solving: no allocation, no `unsafe`,
and `no_std`-capable.

---

## Features

* ✅ Reference-exact NESTED and RING indexing, verified against `healpy`
* 🌌 `(lon, lat)`, colatitude/longitude, and unit-vector entry points
* 🧭 Cell centres, corners, and the eight neighbours of a cell
* 🔭 Cone search returning **sorted index ranges**, ready to slice a sorted catalogue
* 🌲 Hierarchy: parents, children, and descendant ranges at any depth
* ⚡ ~10 ns per position-to-cell conversion, allocation-free throughout
* 🦀 100% safe Rust (`#![forbid(unsafe_code)]`)
* 🚫 Real `no_std` support via `libm`

---

## Quick start

```rust
use realpix::nested;

// `get` is a `const fn`: binding the layer to a `const` folds every derived constant
// into the call site.
const LAYER: nested::Layer = nested::get(12);

let cell = LAYER.hash(1.549_729, 0.129_277);   // (lon, lat) in radians
let (lon, lat) = LAYER.center(cell);
let corners = LAYER.vertices(cell);            // [N, W, S, E] as unit vectors
let neighbours = LAYER.neighbours(cell);       // [Option<u64>; 8]

// Cross to the RING scheme and back.
let r = LAYER.to_ring(cell);
assert_eq!(cell, realpix::ring::get(12).to_nested(r));

// Which ranges of a catalogue sorted by NESTED index cover this field of view?
LAYER.cone_coverage(realpix::lonlat_to_vec(1.55, 0.13), 5f64.to_radians(), |range| {
    // slice the catalogue by `range`
    let _ = range;
});
```

---

## What is HEALPix?

HEALPix divides the sphere into 12 base cells, each subdivided into an `nside × nside`
grid, for `12 × nside²` cells of **exactly equal area**. `nside` is always a power of two,
so a layer is identified by its **depth** (also called order), with `nside = 2^depth`.

| depth | nside | cells | cell size |
| ----: | ----: | ----: | --------: |
| 0 | 1 | 12 | 58.6° |
| 5 | 32 | 12 288 | 1.8° |
| 10 | 1 024 | 12.6 M | 3.4′ |
| 16 | 65 536 | 51.5 G | 3.2″ |
| 29 | 536 870 912 | 3.5 × 10¹⁸ | 0.4 mas |

Depth 29 is the deepest layer whose NESTED indices fit in a `u64`, and is the maximum this
crate supports.

---

## RING vs NESTED

Both schemes describe the same cells; they differ only in how those cells are numbered.

**NESTED** numbers cells hierarchically, so the `4^k` descendants of a cell form one
contiguous range at every deeper depth. Spatially close cells have close indices, which is
what makes range queries and multi-resolution indexing work. **Use this for catalogues,
quad matching and cone searches.**

**RING** numbers cells ring by ring from the north pole, so indices are ordered by
decreasing latitude. This is the layout spherical-harmonic transforms expect.

---

## Cone search

`cone_coverage` walks the quad-tree from the 12 base cells, emitting a whole subtree as a
single range as soon as it is provably inside the cone, and refining only the cells that
straddle its boundary. The result is:

* **sorted, disjoint and non-adjacent** — ready to binary-search a catalogue with,
* **inclusive** — a superset of the exact cover. A few cells that come close to the
  boundary without touching it can be included, but a cell that intersects the cone is
  never missed. Measured over-inclusion is under 1% for discs spanning many cells and at
  most ~20% for the smallest ones.
* **allocation-free** — results are handed to a closure. `cone_coverage_ranges` and
  `cone_coverage_cells` collect into a `Vec` when the `alloc` feature is on.

---

## Coordinate conventions

Angles are radians throughout.

| | |
| --- | --- |
| `lon`, right ascension | `[0, 2π)`, wrapped on input |
| `lat`, declination | `[-π/2, π/2]` |
| `theta` (colatitude) | `[0, π]`, `0` at the north pole |
| unit vector | `[x, y, z]`, `z = sin(lat)` |

`hash(lon, lat)` and `hash_theta_phi(theta, phi)` differ by one ulp in how `z` is derived,
which only matters for a position landing exactly on a cell boundary. `hash_vec` accepts a
vector of any length and is the most accurate entry point near the poles.

---

## Performance

Measured on an Apple M-series laptop (`cargo bench`), independent of depth:

| operation | time |
| --- | ---: |
| `hash(lon, lat)` → cell | 10.6 ns |
| `hash_vec(v)` → cell | 17.0 ns |
| `center(cell)` → `(lon, lat)` | 15.8 ns |
| `center_vec(cell)` → vector | 17.5 ns |
| `to_ring(cell)` | 5.3 ns |
| `neighbours(cell)` | 22.7 ns |
| `cone_coverage`, depth 12, r = 0.05 rad (127 k cells) | 320 µs |
| `cone_coverage`, depth 8, r = 0.05 rad (563 cells) | 32 µs |

The hot paths perform no allocation, no integer division and no modulo: `nside` is a power
of two, so every derived quantity is a shift or a mask.

---

## Features

| feature | default | effect |
| --- | :---: | --- |
| `std` | ✅ | floating point from `std`; implies `alloc` |
| `alloc` | ✅ | `Vec`-returning convenience methods |
| `libm` | | floating point from `libm`, for `no_std` builds |
| `latlong` | | `RaDec` entry points via the [`latlong`](https://crates.io/crates/latlong) crate (implies `std`) |

For a bare-metal target:

```toml
realpix = { version = "0.2", default-features = false, features = ["libm"] }
```

---

## Correctness

The test suite pins the implementation to the official HEALPix numbering rather than to
itself:

* **Golden vectors** generated from `healpy` (`tools/gen_golden.py`) and committed under
  `tests/data/`: ~4 700 position-to-cell samples across depths 0–29 including poles, the
  belt/cap transition, and base-cell boundaries; cell centres; scheme conversions;
  neighbours; cell corners; and reference discs. Indices must match exactly.
* **Exhaustive round-trips** to depth 6 and sampled at every depth to 29, for both
  schemes and both the angle and vector entry points.
* **Structural invariants**: the schemes are a bijection, RING indices are ordered by
  latitude, ring lengths are `4r`/`4·nside`/`4r`, neighbours are symmetric, exactly 24
  cells have seven neighbours, and cell corners are shared with their neighbours.
* **Equal area**, by Monte Carlo: uniformly random directions must populate every cell
  within Poisson noise.
* **Cone inclusiveness**: for random discs, every point inside the disc lands in a covered
  cell, and the coverage contains `healpy`'s exact disc.
* The geometric bound the cone search prunes with is re-derived and re-checked by the test
  suite at every depth.

Regenerate the golden files with:

```
uv run --python 3.12 --with healpy --with numpy python tools/gen_golden.py         # rewrite
uv run --python 3.12 --with healpy --with numpy python tools/gen_golden.py --check # verify
```

---

## Status

* ✔ NESTED and RING indexing, reference-exact
* ✔ Cell centres, corners, neighbours, hierarchy, `uniq` (NUNIQ) encoding
* ✔ Cone search as index ranges
* ✔ Gnomonic (tangent plane) projection
* 🚧 Polygon and elliptical-cone queries
* 🚧 Bilinear interpolation weights

---

## License

<sup>
Licensed under either of <a href="LICENSE-APACHE">Apache License, Version
2.0</a> or <a href="LICENSE-MIT">MIT license</a> at your option.
</sup>

<br>

<sub>
Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this crate by you, as defined in the Apache-2.0 license, shall
be dual licensed as above, without any additional terms or conditions.
</sub>

---

## Inspiration

* The HEALPix reference implementation (Górski et al. 2005)
* healpy
* astrometry.net
