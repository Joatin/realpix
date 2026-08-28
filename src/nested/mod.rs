//! The NESTED numbering scheme.
//!
//! Cells are numbered so that the `4^k` descendants of a cell form a contiguous index
//! range at every deeper depth. That makes NESTED the scheme of choice for spatial
//! indexing: a catalogue sorted by NESTED index at depth `d` can answer a cone search by
//! slicing a handful of index ranges.

mod cone;

pub use cone::MAX_CENTER_TO_VERTEX;

use crate::base::Base;
use crate::depth::check_depth;
use crate::error::{Error, Result};
use crate::proj::Loc;
use crate::xyf::{Direction, spread_bits};
use crate::{MAX_DEPTH, Vec3};
use core::ops::Range;

/// A NESTED HEALPix layer at a fixed depth.
///
/// Obtained from [`get`] (or [`checked_get`]). The struct is `Copy` and holds only
/// precomputed constants, so it is cheap to pass around and free to construct.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layer {
    pub(crate) base: Base,
}

/// Returns the NESTED layer at `depth`.
///
/// This is a `const fn`: binding the layer to a `const` lets the optimiser fold every
/// depth-derived constant into the call site.
///
/// # Panics
/// Panics if `depth > `[`MAX_DEPTH`].
///
/// ```
/// const LAYER: realpix::nested::Layer = realpix::nested::get(10);
/// assert_eq!(LAYER.nside(), 1024);
/// ```
#[inline(always)]
pub const fn get(depth: u8) -> Layer {
    Layer {
        base: Base::new(depth),
    }
}

/// Returns the NESTED layer at `depth`, or [`Error::InvalidDepth`] if it is out of range.
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

    /// The area of one cell, in steradians.
    #[inline]
    pub fn cell_area(&self) -> f64 {
        (4.0 * crate::math::PI) / self.base.n_hash as f64
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

    // ---------------------------------------------------------------- position -> cell

    /// The cell containing the given position.
    ///
    /// `lon` and `lat` are in radians; `lon` may be any finite value and is wrapped into
    /// `[0, 2π)`, `lat` must be in `[-π/2, π/2]`.
    ///
    /// ```
    /// let layer = realpix::nested::get(0);
    /// assert_eq!(layer.hash(0.0, 0.0), 4); // on the equator, at the corner of base cell 4
    /// assert_eq!(layer.hash(0.0, std::f64::consts::FRAC_PI_2), 0); // north pole
    /// ```
    #[inline]
    pub fn hash(&self, lon: f64, lat: f64) -> u64 {
        self.base
            .xyf2nested(self.base.loc2xyf(Loc::from_lonlat(lon, lat)))
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
            .xyf2nested(self.base.loc2xyf(Loc::from_theta_phi(theta, phi)))
    }

    /// The cell containing the direction `v`.
    ///
    /// The vector does not have to be normalised. This is the entry point to prefer when
    /// the caller already holds a direction: it is more accurate than
    /// [`hash`](Self::hash) within 8 degrees of a pole, where deriving the colatitude
    /// from `z` alone loses precision.
    #[inline]
    pub fn hash_vec(&self, v: Vec3) -> u64 {
        self.base.xyf2nested(self.base.loc2xyf(Loc::from_vec(&v)))
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

    // ---------------------------------------------------------------- cell -> position

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
    /// Cheaper and more accurate than [`center`](Self::center) near the poles: no inverse
    /// trigonometric function is involved.
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
    pub(crate) fn center_loc(&self, cell: u64) -> Loc {
        self.assert_cell(cell);
        let xyf = self.base.nested2xyf(cell);
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
        let xyf = self.base.nested2xyf(cell);
        let (x, y) = (xyf.x as f64, xyf.y as f64);
        [
            self.base.xyf2loc(xyf.face, x + 1.0, y + 1.0).to_vec(),
            self.base.xyf2loc(xyf.face, x, y + 1.0).to_vec(),
            self.base.xyf2loc(xyf.face, x, y).to_vec(),
            self.base.xyf2loc(xyf.face, x + 1.0, y).to_vec(),
        ]
    }

    /// The four corners of `cell` as `(lon, lat)` pairs, in the order north, west, south, east.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn vertices_lonlat(&self, cell: u64) -> [(f64, f64); 4] {
        self.assert_cell(cell);
        let xyf = self.base.nested2xyf(cell);
        let (x, y) = (xyf.x as f64, xyf.y as f64);
        [
            self.base.xyf2loc(xyf.face, x + 1.0, y + 1.0).to_lonlat(),
            self.base.xyf2loc(xyf.face, x, y + 1.0).to_lonlat(),
            self.base.xyf2loc(xyf.face, x, y).to_lonlat(),
            self.base.xyf2loc(xyf.face, x + 1.0, y).to_lonlat(),
        ]
    }

    // ------------------------------------------------------------------- neighbours

    /// The eight neighbours of `cell`, indexed by [`Direction`].
    ///
    /// Exactly 24 cells in the layer — the ones sitting on a base-cell corner where only
    /// three base cells meet — have one `None` entry; every other cell has eight
    /// neighbours.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn neighbours(&self, cell: u64) -> [Option<u64>; 8] {
        self.assert_cell(cell);
        let xyf = self.base.nested2xyf(cell);
        let nsm1 = self.base.nside_minus_1;

        if xyf.x > 0 && xyf.x < nsm1 && xyf.y > 0 && xyf.y < nsm1 {
            // Interior of a base cell: the neighbours differ from `cell` only in the
            // interleaved bits, so build them directly without a table lookup.
            let fpix = (xyf.face as u64) << self.base.twice_depth;
            let px0 = spread_bits(xyf.x);
            let py0 = spread_bits(xyf.y) << 1;
            let pxp = spread_bits(xyf.x + 1);
            let pyp = spread_bits(xyf.y + 1) << 1;
            let pxm = spread_bits(xyf.x - 1);
            let pym = spread_bits(xyf.y - 1) << 1;
            [
                Some(fpix + pxm + py0),
                Some(fpix + pxm + pyp),
                Some(fpix + px0 + pyp),
                Some(fpix + pxp + pyp),
                Some(fpix + pxp + py0),
                Some(fpix + pxp + pym),
                Some(fpix + px0 + pym),
                Some(fpix + pxm + pym),
            ]
        } else {
            let n = self.base.neighbours_xyf(xyf);
            let mut out = [None; 8];
            for i in 0..8 {
                out[i] = n[i].map(|x| self.base.xyf2nested(x));
            }
            out
        }
    }

    /// The neighbour of `cell` in one direction, if it exists.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn neighbour(&self, cell: u64, direction: Direction) -> Option<u64> {
        self.assert_cell(cell);
        let xyf = self.base.nested2xyf(cell);
        self.base.neighbours_xyf(xyf)[direction.index()].map(|x| self.base.xyf2nested(x))
    }

    // -------------------------------------------------------------------- hierarchy

    /// The ancestor of `cell` at `parent_depth`.
    ///
    /// # Panics
    /// Panics if `parent_depth > self.depth()` or if `cell` is out of range.
    #[inline]
    pub fn parent(&self, cell: u64, parent_depth: u8) -> u64 {
        self.assert_cell(cell);
        assert!(
            parent_depth <= self.depth(),
            "parent depth must not exceed this layer's depth"
        );
        cell >> ((self.depth() - parent_depth) << 1)
    }

    /// The four children of `cell`, one depth down.
    ///
    /// # Panics
    /// Panics if `cell` is out of range, or if this layer is already at [`MAX_DEPTH`].
    #[inline]
    pub fn children(&self, cell: u64) -> [u64; 4] {
        self.assert_cell(cell);
        assert!(self.depth() < MAX_DEPTH, "no layer below MAX_DEPTH");
        let first = cell << 2;
        [first, first + 1, first + 2, first + 3]
    }

    /// The contiguous range of descendants of `cell` at `child_depth`.
    ///
    /// This is the range a catalogue sorted by NESTED index at `child_depth` should be
    /// sliced by to obtain everything inside `cell`.
    ///
    /// # Panics
    /// Panics if `child_depth < self.depth()`, `child_depth > `[`MAX_DEPTH`], or `cell` is
    /// out of range.
    #[inline]
    pub fn children_range(&self, cell: u64, child_depth: u8) -> Range<u64> {
        self.assert_cell(cell);
        assert!(
            child_depth >= self.depth() && child_depth <= MAX_DEPTH,
            "child depth must be between this layer's depth and MAX_DEPTH"
        );
        let shift = (child_depth - self.depth()) << 1;
        (cell << shift)..((cell + 1) << shift)
    }

    // ---------------------------------------------------------------- scheme change

    /// The RING index of the cell that `cell` denotes.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn to_ring(&self, cell: u64) -> u64 {
        self.assert_cell(cell);
        self.base.xyf2ring(self.base.nested2xyf(cell))
    }

    /// Iterates over every cell of this layer, in index order.
    #[inline]
    pub fn iter(&self) -> Range<u64> {
        0..self.base.n_hash
    }
}
