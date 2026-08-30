# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) — while the major version is
`0`, a bump of the minor version may break the API.

## [Unreleased]

## [0.2.0]

A complete rewrite of the core. **Read the first section before upgrading**: this release
is not only a new API, it is a correctness fix, and cell indices computed by 0.1.x are not
valid.

### Fixed

- **Cell indices were wrong.** The projection in 0.1.x did not implement HEALPix. Compared
  against this release — which is pinned bit-for-bit to the reference C++ implementation
  and to `healpy` by golden vectors committed under `tests/data/` — the two agree on only
  **2.6%** of random sky positions at `nside = 32`. Indices produced by 0.1.x therefore did
  not interoperate with `healpy`, `astropy`, or any other HEALPix tool, and neighbouring
  cells were not generally adjacent on the sky.

  **If you have stored cell indices from 0.1.x, recompute them from the original
  coordinates.** They cannot be converted; there is no mapping from the old numbering to
  the correct one.

- The `Pixel<N>` wrapper did not provide the type safety it appeared to. Its numbering
  scheme was a `PhantomData` parameter that no constructor checked, so `Pixel::from_u64`
  would mint a pixel of any scheme from any integer, and a NESTED index could be read as a
  RING one without a compile error or a runtime one.

### Changed

The `Healpix` trait, the `NumberingScheme` trait and the `Pixel<N>` wrapper are gone. A
layer is now a plain `Copy` struct obtained from the scheme's module, and a cell is a bare
`u64`. Porting is mechanical:

| 0.1.x | 0.2 |
| --- | --- |
| `ConstHealpix::<32>::new()` | `nested::get(5)` — depth, not `nside`; `const` too |
| `DynamicHealpix::new(nside)?` | `nested::checked_get(depth_from_nside(nside)?)?` |
| `Nested` / `Ring` marker types | the `nested` and `ring` modules |
| `Pixel<N>` | `u64` |
| `pixel.as_u64()` / `Pixel::from_u64(v)` | the value itself |
| `.face_resolution()` | `.nside()` |
| `.pixels_per_face()` | `.n_hash() / 12` |
| `.total_pixels()` → `u32` | `.n_hash()` → `u64` |
| `.angle_to_pixel(theta, phi)` | `.hash_theta_phi(theta, phi)`, or `.hash(lon, lat)` |
| `.pixel_to_angle(pixel)?` → `(theta, phi)` | `.center(cell)` → `(lon, lat)`, see below |
| `.ra_dec_to_pixel(&rd)` | `.hash_ra_dec(&rd)` (feature `latlong`) |
| `.pixel_to_ra_dec(pixel)?` | `.center_ra_dec(cell)` (feature `latlong`) |
| `.iter_pixels()` | `.iter()` |
| `.project_ra_dec(pixel, &rd)` | `.project_ra_dec(cell, &rd)` (feature `latlong`) |
| `gnomonic_project(center, &rd)` | `gnomonic_project(lon, lat, lon, lat)` — radians, not `RaDec` |

Three traps worth calling out:

- **`pixel_to_angle` returned `(theta, phi)`; `center` returns `(lon, lat)`.** Both the
  order and the convention changed: `lon` is `phi`, but `lat` is `π/2 - theta`. Use
  `.hash_theta_phi` / the `theta` accessors if you would rather keep colatitude.
- **Layers are addressed by depth, not by `nside`.** `nside = 2^depth`, so
  `ConstHealpix<32>` becomes `nested::get(5)`. `depth_from_nside` converts at runtime.
- **The depth limit moved.** `total_pixels()` returned a `u32`, capping 0.1.x at about
  depth 14. Cell counts are now `u64` and depths run to 29.

Other changes:

- `Error` gained `InvalidDepth`, `InvalidNside`, `InvalidCell` and `InvalidCoordinate`, and
  is now `#[non_exhaustive]` with a `Display` impl and `core::error::Error`.
  `InvalidFaceResolution` and `InvalidPixel` are gone.
- The crate is now `no_std` by default-off, with floating point from `std` or `libm`.
- Minimum supported Rust version is **1.88**, declared in `Cargo.toml`.

### Added

- **RING as a first-class scheme**, with the same API surface as NESTED and exact
  conversion both ways (`to_ring`, `to_nested`).
- **Cone search** in both schemes, returning sorted, disjoint, non-adjacent index ranges
  ready to slice a catalogue. Allocation-free into a closure, with `Vec`-returning
  variants behind `alloc`.
- **Multi-order coverage maps** (`moc::Moc`): arbitrary sky regions with exact union,
  intersection, difference, symmetric difference and complement, plus NUNIQ import and
  export. Behind the `alloc` feature.
- **Cell geometry**: `center`, `center_vec`, `vertices`, `vertices_lonlat`, `cell_area`.
- **Topology**: `neighbours`, `neighbour`, and `external_edge` — the ring of cells
  surrounding a cell at any deeper resolution.
- **Hierarchy**: `parent`, `children`, `children_range`.
- **Bulk hashing**: `hash_many` and `hash_many_vec` fill a caller-owned buffer.
- `to_uniq` / `from_uniq` for the NUNIQ multi-resolution encoding.
- Vector entry points throughout (`hash_vec`, `center_vec`), and `lonlat_to_vec`,
  `vec_to_lonlat`, `angular_distance`.
- Checked, non-panicking variants: `checked_get`, `checked_hash`, `checked_center`.
- `max_center_to_vertex(depth)` — the radius of the smallest disc centred on a cell that
  contains the whole cell, which is what a search margin has to allow for. This is the
  bound both cone searches prune with, exposed as a function rather than as the lookup
  table behind it, so the table stays an implementation detail.

### Removed

- `Healpix`, `NumberingScheme`, `Pixel`, `ConstHealpix`, `DynamicHealpix`, `Nested`,
  `Ring` (the marker types — the modules of those names replace them).

## [0.1.3] - 2026-01-18

Earlier releases. See the correctness note under 0.2.0 before using them.

[Unreleased]: https://github.com/Joatin/realpix/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Joatin/realpix/compare/v0.1.3...v0.2.0
[0.1.3]: https://crates.io/crates/realpix/0.1.3
