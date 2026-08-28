//! `realpix` — a pure-Rust implementation of the [HEALPix](https://healpix.sourceforge.io/)
//! (Hierarchical Equal Area isoLatitude Pixelisation) tessellation of the sphere.
//!
//! The sphere is split into 12 equal-area base cells, each recursively subdivided into
//! `4^depth` cells, for a total of `12 * 4^depth` cells of exactly equal area. Two cell
//! numbering schemes are provided, both bit-compatible with the reference C++ implementation
//! (and therefore with `healpy`):
//!
//! * [`nested`] — hierarchical (quad-tree) numbering. Spatially local: the cells of a parent
//!   form a contiguous index range at every deeper depth. This is the scheme to use for
//!   catalogue indexing and cone searches.
//! * [`ring`] — iso-latitude numbering, ordered by ring from the north pole. Used by
//!   spherical-harmonic transforms.
//!
//! # Conventions
//!
//! * Angles are **radians**. Positions are given as `(lon, lat)` with `lon` in `[0, 2π)`
//!   (right ascension) and `lat` in `[-π/2, π/2]` (declination), or as unit vectors
//!   `[x, y, z]` with `z = sin(lat)`.
//! * Cell indices are `u64`. `depth` (also called *order*) is a `u8` in `0..=`[`MAX_DEPTH`],
//!   and `nside = 2^depth`.
//!
//! # Example
//!
//! ```
//! use realpix::nested;
//!
//! // `get` is a `const fn`, so the layer constants fold away when the depth is fixed.
//! const LAYER: nested::Layer = nested::get(12);
//!
//! let cell = LAYER.hash(1.234, -0.567);
//! let (lon, lat) = LAYER.center(cell);
//! assert_eq!(cell, LAYER.hash(lon, lat));
//!
//! // The same cell in the RING scheme, and back.
//! let r = LAYER.to_ring(cell);
//! assert_eq!(cell, realpix::ring::get(12).to_nested(r));
//! ```
#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

#[cfg(feature = "alloc")]
extern crate alloc;
#[cfg(feature = "std")]
extern crate std;

pub mod nested;
pub mod ring;

mod base;
mod depth;
mod error;
mod math;
mod proj;
mod tangent;
mod uniq;
mod xyf;

#[cfg(feature = "latlong")]
pub mod radec;

pub use self::depth::{MAX_DEPTH, depth_from_nside, n_hash, nside};
pub use self::error::{Error, Result};
pub use self::tangent::{angular_distance, gnomonic_project, lonlat_to_vec, vec_to_lonlat};
pub use self::uniq::{from_uniq, to_uniq};
pub use self::xyf::Direction;

/// A unit vector on the sphere, `[x, y, z]`, with `z = sin(lat)`.
pub type Vec3 = [f64; 3];
