//! The RING numbering scheme.
//!
//! Cells are numbered ring by ring from the north pole, so cell indices are ordered by
//! decreasing latitude. This is the layout expected by spherical-harmonic transforms; for
//! spatial indexing prefer [`crate::nested`].

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
#[inline(always)]
pub const fn get(depth: u8) -> Layer {
    Layer {
        base: Base::new(depth),
    }
}

/// Returns the RING layer at `depth`, or [`Error::InvalidDepth`] if it is out of range.
#[inline]
pub const fn checked_get(depth: u8) -> Result<Layer> {
    match check_depth(depth) {
        Ok(()) => Ok(get(depth)),
        Err(e) => Err(e),
    }
}

impl Layer {
    /// The depth (order) of this layer.
    #[inline(always)]
    pub const fn depth(&self) -> u8 {
        self.base.depth
    }

    /// `nside = 2^depth`.
    #[inline(always)]
    pub const fn nside(&self) -> u32 {
        self.base.nside
    }

    /// The number of cells in this layer, `12 * 4^depth`.
    #[inline(always)]
    pub const fn n_hash(&self) -> u64 {
        self.base.n_hash
    }

    /// The number of cells in the north polar cap, `2 * nside * (nside - 1)`.
    #[inline(always)]
    pub const fn n_cap(&self) -> u64 {
        self.base.ncap
    }

    /// Returns whether `cell` is a valid index in this layer.
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
    #[inline]
    pub fn hash_theta_phi(&self, theta: f64, phi: f64) -> u64 {
        self.base
            .xyf2ring(self.base.loc2xyf(Loc::from_theta_phi(theta, phi)))
    }

    /// The cell containing the direction `v`.
    #[inline]
    pub fn hash_vec(&self, v: Vec3) -> u64 {
        self.base.xyf2ring(self.base.loc2xyf(Loc::from_vec(&v)))
    }

    /// [`hash`](Self::hash) with the arguments validated.
    ///
    /// # Errors
    /// Returns [`Error::InvalidCoordinate`] if either coordinate is not finite or `lat` is
    /// outside `[-π/2, π/2]`.
    #[inline]
    pub fn checked_hash(&self, lon: f64, lat: f64) -> Result<u64> {
        if !lon.is_finite() || !lat.is_finite() || lat.abs() > crate::math::FRAC_PI_2 {
            return Err(Error::InvalidCoordinate);
        }
        Ok(self.hash(lon, lat))
    }

    /// The centre of `cell`, as `(lon, lat)` in radians with `lon` in `[0, 2π)`.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn center(&self, cell: u64) -> (f64, f64) {
        self.center_loc(cell).to_lonlat()
    }

    /// The centre of `cell`, as a unit vector.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn center_vec(&self, cell: u64) -> Vec3 {
        self.center_loc(cell).to_vec()
    }

    /// [`center`](Self::center) with the cell index validated.
    ///
    /// # Errors
    /// Returns [`Error::InvalidCell`] if `cell` is out of range for this depth.
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
    #[inline]
    pub fn neighbour(&self, cell: u64, direction: Direction) -> Option<u64> {
        self.neighbours(cell)[direction.index()]
    }

    /// The NESTED index of the cell that `cell` denotes.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn to_nested(&self, cell: u64) -> u64 {
        self.assert_cell(cell);
        self.base.xyf2nested(self.base.ring2xyf(cell))
    }

    /// Iterates over every cell of this layer, in index order.
    #[inline]
    pub fn iter(&self) -> Range<u64> {
        0..self.base.n_hash
    }
}
