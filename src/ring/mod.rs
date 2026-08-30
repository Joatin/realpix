//! The RING numbering scheme.
//!
//! Cells are numbered ring by ring from the north pole, so cell indices are ordered by
//! decreasing latitude. This is the layout expected by spherical-harmonic transforms; for
//! spatial indexing prefer [`crate::nested`], whose indices keep spatially close cells
//! close.
//!
//! [`Layer::cone_coverage`] searches a disc in this scheme too, covering the disc about as
//! tightly and in about as many ranges as the NESTED search does, and rather faster — it
//! visits only the rings the disc reaches, where the NESTED descent starts from the whole
//! sphere. Pick the scheme your catalogue is already sorted in.

mod cone;

use crate::Vec3;
use crate::base::Base;
use crate::depth::check_depth;
use crate::error::{Error, Result};
use crate::proj::Loc;
use crate::xyf::Direction;
use core::ops::Range;

/// A RING HEALPix layer at a fixed depth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layer {
    pub(crate) base: Base,
}

/// Returns the RING layer at `depth`.
///
/// # Panics
/// Panics if `depth > `[`crate::MAX_DEPTH`].
///
/// ```
/// const LAYER: realpix::ring::Layer = realpix::ring::get(10);
/// assert_eq!(LAYER.nside(), 1024);
/// ```
#[inline(always)]
pub const fn get(depth: u8) -> Layer {
    Layer {
        base: Base::new(depth),
    }
}

/// Returns the RING layer at `depth`, or [`Error::InvalidDepth`] if it is out of range.
///
/// The non-panicking form of [`get`], for a depth that came from outside the program — a
/// command line, a file header, a network message.
///
/// # Errors
/// Returns [`Error::InvalidDepth`] if `depth > `[`crate::MAX_DEPTH`].
///
/// ```
/// assert_eq!(realpix::ring::checked_get(12).unwrap().depth(), 12);
/// assert!(realpix::ring::checked_get(30).is_err());
/// ```
#[inline]
pub const fn checked_get(depth: u8) -> Result<Layer> {
    match check_depth(depth) {
        Ok(()) => Ok(get(depth)),
        Err(e) => Err(e),
    }
}

impl Layer {
    /// The depth (order) of this layer.
    ///
    /// ```
    /// assert_eq!(realpix::ring::get(12).depth(), 12);
    /// ```
    #[inline(always)]
    pub const fn depth(&self) -> u8 {
        self.base.depth
    }

    /// `nside = 2^depth`, the number of cells along one edge of a base cell.
    ///
    /// ```
    /// assert_eq!(realpix::ring::get(12).nside(), 4096);
    /// ```
    #[inline(always)]
    pub const fn nside(&self) -> u32 {
        self.base.nside
    }

    /// The number of cells in this layer, `12 * 4^depth`.
    ///
    /// Every valid cell index of this layer is below it, so it is also the exclusive end of
    /// [`iter`](Self::iter).
    ///
    /// ```
    /// assert_eq!(realpix::ring::get(1).n_hash(), 48);
    /// ```
    #[inline(always)]
    pub const fn n_hash(&self) -> u64 {
        self.base.n_hash
    }

    /// The number of cells in the north polar cap, `2 * nside * (nside - 1)`.
    ///
    /// RING indices below this are in the north cap, indices from here to
    /// `n_hash() - n_cap()` are in the equatorial belt, and the rest are the south cap.
    /// The three regions have different ring lengths, which is why the boundary is worth
    /// knowing.
    ///
    /// ```
    /// let layer = realpix::ring::get(4);
    /// assert_eq!(layer.n_cap(), 2 * 16 * 15);
    /// // The first cell of the layer is the one nearest the north pole.
    /// assert!(layer.center(0).1 > layer.center(layer.n_cap()).1);
    /// ```
    #[inline(always)]
    pub const fn n_cap(&self) -> u64 {
        self.base.ncap
    }

    /// Returns whether `cell` is a valid index in this layer.
    ///
    /// Use it to guard the methods that panic on an out-of-range cell.
    ///
    /// ```
    /// let layer = realpix::ring::get(1);
    /// assert!(layer.contains(47));
    /// assert!(!layer.contains(48));
    /// ```
    #[inline(always)]
    pub const fn contains(&self, cell: u64) -> bool {
        self.base.contains(cell)
    }

    #[inline(always)]
    fn assert_cell(&self, cell: u64) {
        assert!(
            self.contains(cell),
            "cell index out of range for this depth"
        );
    }

    /// The cell containing the given position, `lon` and `lat` in radians.
    ///
    /// ```
    /// let layer = realpix::ring::get(10);
    /// let cell = layer.hash(1.0, 0.5);
    /// // The same cell the NESTED scheme finds, under RING's numbering.
    /// assert_eq!(cell, realpix::nested::get(10).to_ring(realpix::nested::get(10).hash(1.0, 0.5)));
    /// ```
    #[inline]
    pub fn hash(&self, lon: f64, lat: f64) -> u64 {
        self.base
            .xyf2ring(self.base.loc2xyf(Loc::from_lonlat(lon, lat)))
    }

    /// The cell containing the given position, in the HEALPix-native spherical convention:
    /// `theta` is the colatitude in `[0, π]` (`0` at the north pole) and `phi` the
    /// longitude in radians.
    ///
    /// Equivalent to `hash(phi, π/2 - theta)` up to one ulp; use whichever convention the
    /// caller already holds.
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let layer = realpix::ring::get(10);
    /// // Colatitude 0 is the north pole; the same point as latitude +π/2.
    /// assert_eq!(layer.hash_theta_phi(0.0, 1.0), layer.hash(1.0, FRAC_PI_2));
    /// ```
    #[inline]
    pub fn hash_theta_phi(&self, theta: f64, phi: f64) -> u64 {
        self.base
            .xyf2ring(self.base.loc2xyf(Loc::from_theta_phi(theta, phi)))
    }

    /// The cell containing the direction `v`.
    ///
    /// ```
    /// let layer = realpix::ring::get(10);
    /// // The vector need not be normalised.
    /// assert_eq!(layer.hash_vec([2.0, 0.0, 0.0]), layer.hash_vec([1.0, 0.0, 0.0]));
    /// assert_eq!(layer.hash_vec(realpix::lonlat_to_vec(1.0, 0.5)), layer.hash(1.0, 0.5));
    /// ```
    #[inline]
    pub fn hash_vec(&self, v: Vec3) -> u64 {
        self.base.xyf2ring(self.base.loc2xyf(Loc::from_vec(&v)))
    }

    /// [`hash`](Self::hash) with the arguments validated.
    ///
    /// # Errors
    /// Returns [`Error::InvalidCoordinate`] if either coordinate is not finite or `lat` is
    /// outside `[-π/2, π/2]`.
    ///
    /// ```
    /// let layer = realpix::ring::get(10);
    /// assert_eq!(layer.checked_hash(1.0, 0.5), Ok(layer.hash(1.0, 0.5)));
    /// assert!(layer.checked_hash(1.0, 2.0).is_err());
    /// assert!(layer.checked_hash(f64::NAN, 0.0).is_err());
    /// ```
    #[inline]
    pub fn checked_hash(&self, lon: f64, lat: f64) -> Result<u64> {
        if !lon.is_finite() || !lat.is_finite() || lat.abs() > crate::math::FRAC_PI_2 {
            return Err(Error::InvalidCoordinate);
        }
        Ok(self.hash(lon, lat))
    }

    // ------------------------------------------------------------------------- bulk

    /// Hashes every position in `positions` into `out`, in order.
    ///
    /// Equivalent to calling [`hash`](Self::hash) in a loop, and about as fast: the cost of
    /// a hash is dominated by the `sin` it does per position, which no amount of batching
    /// removes. What this saves is the bookkeeping at the call site — it fills a buffer you
    /// already own, so nothing is allocated and the depth constants are loaded once.
    ///
    /// # Panics
    /// Panics if `out` is not the same length as `positions`.
    ///
    /// ```
    /// let layer = realpix::ring::get(12);
    /// let sources = [(1.549_729, 0.129_277), (1.372_198, -0.143_146)];
    /// let mut cells = [0u64; 2];
    /// layer.hash_many(&sources, &mut cells);
    /// assert_eq!(cells[1], layer.hash(sources[1].0, sources[1].1));
    /// ```
    pub fn hash_many(&self, positions: &[(f64, f64)], out: &mut [u64]) {
        assert_eq!(
            positions.len(),
            out.len(),
            "output slice must match the input length"
        );
        for (cell, (lon, lat)) in out.iter_mut().zip(positions) {
            *cell = self.hash(*lon, *lat);
        }
    }

    /// Hashes every direction in `positions` into `out`, in order.
    ///
    /// Equivalent to calling [`hash_vec`](Self::hash_vec) in a loop; see
    /// [`hash_many`](Self::hash_many) for what batching does and does not buy.
    ///
    /// # Panics
    /// Panics if `out` is not the same length as `positions`.
    ///
    /// ```
    /// let layer = realpix::ring::get(12);
    /// let sources = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    /// let mut cells = [0u64; 2];
    /// layer.hash_many_vec(&sources, &mut cells);
    /// assert_eq!(cells[0], layer.hash_vec(sources[0]));
    /// ```
    pub fn hash_many_vec(&self, positions: &[Vec3], out: &mut [u64]) {
        assert_eq!(
            positions.len(),
            out.len(),
            "output slice must match the input length"
        );
        for (cell, v) in out.iter_mut().zip(positions) {
            *cell = self.hash_vec(*v);
        }
    }

    /// The centre of `cell`, as `(lon, lat)` in radians with `lon` in `[0, 2π)`.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    ///
    /// ```
    /// let layer = realpix::ring::get(10);
    /// let cell = layer.hash(1.0, 0.5);
    /// let (lon, lat) = layer.center(cell);
    /// // The centre of a cell is inside that cell.
    /// assert_eq!(layer.hash(lon, lat), cell);
    /// ```
    #[inline]
    pub fn center(&self, cell: u64) -> (f64, f64) {
        self.center_loc(cell).to_lonlat()
    }

    /// The centre of `cell`, as a unit vector.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    ///
    /// More accurate than `lonlat_to_vec(center(cell))` near the poles, and cheaper: it
    /// never goes through an inverse trigonometric function.
    ///
    /// ```
    /// let layer = realpix::ring::get(10);
    /// let v = layer.center_vec(layer.hash(1.0, 0.5));
    /// assert!((v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - 1.0).abs() < 1e-15);
    /// ```
    #[inline]
    pub fn center_vec(&self, cell: u64) -> Vec3 {
        self.center_loc(cell).to_vec()
    }

    /// [`center`](Self::center) with the cell index validated.
    ///
    /// # Errors
    /// Returns [`Error::InvalidCell`] if `cell` is out of range for this depth.
    ///
    /// ```
    /// let layer = realpix::ring::get(1);
    /// assert!(layer.checked_center(47).is_ok());
    /// assert!(layer.checked_center(48).is_err());
    /// ```
    #[inline]
    pub fn checked_center(&self, cell: u64) -> Result<(f64, f64)> {
        if !self.contains(cell) {
            return Err(Error::InvalidCell {
                cell,
                depth: self.depth(),
            });
        }
        Ok(self.center(cell))
    }

    #[inline]
    fn center_loc(&self, cell: u64) -> Loc {
        self.assert_cell(cell);
        let xyf = self.base.ring2xyf(cell);
        self.base
            .xyf2loc(xyf.face, xyf.x as f64 + 0.5, xyf.y as f64 + 0.5)
    }

    /// The four corners of `cell` as unit vectors, in the order north, west, south, east.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    /// ```
    /// let layer = realpix::ring::get(8);
    /// let cell = layer.hash(1.0, 0.5);
    /// let [north, west, south, east] = layer.vertices(cell);
    /// // The corners bracket the centre in latitude.
    /// assert!(north[2] > layer.center_vec(cell)[2]);
    /// assert!(south[2] < layer.center_vec(cell)[2]);
    /// let _ = (west, east);
    /// ```
    #[inline]
    pub fn vertices(&self, cell: u64) -> [Vec3; 4] {
        self.assert_cell(cell);
        let xyf = self.base.ring2xyf(cell);
        let (x, y) = (xyf.x as f64, xyf.y as f64);
        [
            self.base.xyf2loc(xyf.face, x + 1.0, y + 1.0).to_vec(),
            self.base.xyf2loc(xyf.face, x, y + 1.0).to_vec(),
            self.base.xyf2loc(xyf.face, x, y).to_vec(),
            self.base.xyf2loc(xyf.face, x + 1.0, y).to_vec(),
        ]
    }

    /// The eight neighbours of `cell`, indexed by [`Direction`].
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    /// ```
    /// let layer = realpix::ring::get(6);
    /// let cell = layer.hash(1.0, 0.5);
    /// // Away from the base-cell corners every cell has all eight.
    /// assert_eq!(layer.neighbours(cell).iter().flatten().count(), 8);
    /// // Neighbourhood is symmetric.
    /// for n in layer.neighbours(cell).into_iter().flatten() {
    ///     assert!(layer.neighbours(n).contains(&Some(cell)));
    /// }
    /// ```
    #[inline]
    pub fn neighbours(&self, cell: u64) -> [Option<u64>; 8] {
        self.assert_cell(cell);
        let n = self.base.neighbours_xyf(self.base.ring2xyf(cell));
        let mut out = [None; 8];
        for i in 0..8 {
            out[i] = n[i].map(|x| self.base.xyf2ring(x));
        }
        out
    }

    /// The neighbour of `cell` in one direction, if it exists.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    ///
    /// ```
    /// use realpix::Direction;
    ///
    /// let layer = realpix::ring::get(6);
    /// let cell = layer.hash(1.0, 0.5);
    /// let north = layer.neighbour(cell, Direction::N).unwrap();
    /// // RING numbers cells by decreasing latitude, so north means a lower index.
    /// assert!(north < cell);
    /// ```
    #[inline]
    pub fn neighbour(&self, cell: u64, direction: Direction) -> Option<u64> {
        self.assert_cell(cell);
        let xyf = self.base.ring2xyf(cell);
        self.base
            .neighbour_xyf(xyf, direction.index())
            .map(|x| self.base.xyf2ring(x))
    }

    /// The NESTED index of the cell that `cell` denotes.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    ///
    /// ```
    /// let layer = realpix::ring::get(8);
    /// let cell = layer.hash(1.0, 0.5);
    /// let nested = layer.to_nested(cell);
    /// // The same cell, so the same centre, and it converts back.
    /// assert_eq!(realpix::nested::get(8).center(nested), layer.center(cell));
    /// assert_eq!(realpix::nested::get(8).to_ring(nested), cell);
    /// ```
    #[inline]
    pub fn to_nested(&self, cell: u64) -> u64 {
        self.assert_cell(cell);
        self.base.xyf2nested(self.base.ring2xyf(cell))
    }

    /// Iterates over every cell of this layer, in index order — which for RING is order of
    /// decreasing latitude.
    ///
    /// ```
    /// let layer = realpix::ring::get(1);
    /// assert_eq!(layer.iter(), 0..48);
    /// assert_eq!(layer.iter().count() as u64, layer.n_hash());
    /// ```
    #[inline]
    pub fn iter(&self) -> Range<u64> {
        0..self.base.n_hash
    }
}
