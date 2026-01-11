[![Cargo Build & Test](https://github.com/Joatin/realpix/actions/workflows/ci.yaml/badge.svg)](https://github.com/Joatin/realpix/actions/workflows/ci.yaml)

# REALPix

**`REALPix`** is a Rust crate providing a **HEALPix-style spherical pixelization**, supporting both **RING** and **NESTED** indexing schemes.

It is designed for **astronomy, astrometry, and fast spatial indexing on the sphere**, with a strong focus on **correctness, performance, and portability**.

---

## Features

* ✅ HEALPix-compatible pixelization
* 🔢 RING and NESTED indexing schemes
* 🌌 RA/Dec and θ/φ support
* ⚡ Fast angle ↔ pixel conversion
* 🧠 Spatial locality with NESTED ordering
* 🦀 Written in safe Rust
* 📦 **`std` enabled by default**
* 🚫 Optional **`no_std`** support

---

## What is HEALPix?

HEALPix (**H**ierarchical **E**qual **A**rea **L**atitude **Pix**elization) divides the sphere into:

* **12 base faces**
* Each face subdivided into an **N × N grid**
* Total pixel count:

```
12 × N²
```

All pixels cover **equal area** on the sphere.

---

## RING vs NESTED ordering

Both schemes describe the **same pixelization**, but differ in how pixels are **numbered**.

### RING ordering

Pixels are numbered in **latitude rings**, starting at the north pole and moving south.

```
North pole
   [ 0  1  2 ]
  [ 3  4  5  6 ]
 [ 7  8  9 10 11 ]
      ...
South pole
```

**Characteristics:**

* Latitude-ordered
* Easy full-sky iteration
* Poor spatial locality
* Commonly used for spherical harmonics

---

### NESTED ordering

Pixels are numbered **hierarchically**, using a quad-tree structure on each face.

```
Base face
┌───────┐
│   0   │
└───────┘

Level 1
┌───┬───┐
│ 0 │ 1 │
├───┼───┤
│ 2 │ 3 │
└───┴───┘

Level 2
┌───┬───┬───┬───┐
│00 │01 │10 │11 │
├───┼───┼───┼───┤
│02 │03 │12 │13 │
├───┼───┼───┼───┤
│20 │21 │30 │31 │
├───┼───┼───┼───┤
│22 │23 │32 │33 │
└───┴───┴───┴───┘
```

**Characteristics:**

* Strong spatial locality
* Hierarchical (multi-resolution)
* Efficient neighbor and range queries
* Ideal for star catalogs and plate solving

> For astrometry, quad matching, and fast spatial indexing, **NESTED ordering is strongly recommended**.

---

## Coordinate conventions

`realpix` supports both:

* **Spherical angles**

    * θ (colatitude): `[0, π]`
    * φ (longitude): `[0, 2π)`
* **Astronomical coordinates**

    * Right Ascension (RA)
    * Declination (Dec)

Standard conversions are used:

```
θ = π/2 − Dec
φ = RA
```

---

## Resolution parameter

HEALPix resolution is controlled by a single parameter (`nside`):

* `nside` is the number of subdivisions **per edge of each base face**
* Must be a **power of two**

```
pixels = 12 × nside²
```

Approximate pixel angular size:

```
pixel size ≈ 2 / nside   radians
```

Examples:

| nside | Pixel size |
| ----: | ---------: |
|    32 |      ~3.6° |
|    64 |      ~1.8° |
|   128 |      ~0.9° |
|   256 |     ~0.45° |

---

## `std` and `no_std`

* **`std` is enabled by default**
* `realpix` can be built in **`no_std` environments**
* No heap allocation is required
* Suitable for:

    * Embedded systems
    * WASM
    * Freestanding / constrained environments

---

## Design goals

* Correct handling of poles and boundaries
* Deterministic, explicit math
* No hidden allocations
* Clear mapping between theory and implementation
* Robust behavior across resolutions

---

## Status

* ✔ RING indexing
* ✔ NESTED indexing
* ✔ RA/Dec ↔ θ/φ conversions
* ✔ Unit-tested across edge cases
* 🚧 Neighbor queries (planned)
* 🚧 Cone / radius searches (planned)

---

## License

MIT OR Apache-2.0

---

## Inspiration

* HEALPix reference implementation
* healpy
* astrometry.net
