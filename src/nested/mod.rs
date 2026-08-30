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
use crate::xyf::{Direction, Xyf, spread_bits};
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
#[must_use]
#[inline(always)]
pub const fn get(depth: u8) -> Layer {
    Layer {
        base: Base::new(depth),
    }
}

/// Returns the NESTED layer at `depth`, or [`Error::InvalidDepth`] if it is out of range.
///
/// The non-panicking form of [`get`], for a depth that came from outside the program — a
/// command line, a file header, a network message.
///
/// # Errors
/// Returns [`Error::InvalidDepth`] if `depth > `[`MAX_DEPTH`].
///
/// ```
/// assert_eq!(realpix::nested::checked_get(12).unwrap().depth(), 12);
/// assert!(realpix::nested::checked_get(30).is_err());
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
    /// assert_eq!(realpix::nested::get(12).depth(), 12);
    /// ```
    #[must_use]
    #[inline(always)]
    pub const fn depth(&self) -> u8 {
        self.base.depth
    }

    /// `nside = 2^depth`, the number of cells along one edge of a base cell.
    ///
    /// ```
    /// assert_eq!(realpix::nested::get(12).nside(), 4096);
    /// ```
    #[must_use]
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
    /// assert_eq!(realpix::nested::get(1).n_hash(), 48);
    /// ```
    #[must_use]
    #[inline(always)]
    pub const fn n_hash(&self) -> u64 {
        self.base.n_hash
    }

    /// The area of one cell, in steradians.
    ///
    /// Every cell of a layer has exactly this area — that is the "equal area" in HEALPix.
    ///
    /// ```
    /// let layer = realpix::nested::get(4);
    /// let sphere = layer.cell_area() * layer.n_hash() as f64;
    /// assert!((sphere - 4.0 * std::f64::consts::PI).abs() < 1e-12);
    /// ```
    #[must_use]
    #[inline]
    pub fn cell_area(&self) -> f64 {
        (4.0 * crate::math::PI) / self.base.n_hash as f64
    }

    /// Returns whether `cell` is a valid index in this layer.
    ///
    /// Use it to guard the methods that panic on an out-of-range cell.
    ///
    /// ```
    /// let layer = realpix::nested::get(1);
    /// assert!(layer.contains(47));
    /// assert!(!layer.contains(48));
    /// ```
    #[must_use]
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
    #[must_use]
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
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let layer = realpix::nested::get(10);
    /// // Colatitude 0 is the north pole; the same point as latitude +π/2.
    /// assert_eq!(layer.hash_theta_phi(0.0, 1.0), layer.hash(1.0, FRAC_PI_2));
    /// ```
    #[must_use]
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
    /// ```
    /// let layer = realpix::nested::get(10);
    /// // The vector need not be normalised.
    /// assert_eq!(layer.hash_vec([2.0, 0.0, 0.0]), layer.hash_vec([1.0, 0.0, 0.0]));
    /// assert_eq!(layer.hash_vec(realpix::lonlat_to_vec(1.0, 0.5)), layer.hash(1.0, 0.5));
    /// ```
    #[must_use]
    #[inline]
    pub fn hash_vec(&self, v: Vec3) -> u64 {
        self.base.xyf2nested(self.base.loc2xyf(Loc::from_vec(&v)))
    }

    /// [`hash`](Self::hash) with the arguments validated.
    ///
    /// # Errors
    /// Returns [`Error::InvalidCoordinate`] if either coordinate is not finite or `lat` is
    /// outside `[-π/2, π/2]`.
    ///
    /// ```
    /// let layer = realpix::nested::get(10);
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

    // ---------------------------------------------------------------- cell -> position

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
    /// let layer = realpix::nested::get(12);
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
    /// let layer = realpix::nested::get(12);
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
    /// let layer = realpix::nested::get(10);
    /// let cell = layer.hash(1.0, 0.5);
    /// let (lon, lat) = layer.center(cell);
    /// // The centre of a cell is inside that cell.
    /// assert_eq!(layer.hash(lon, lat), cell);
    /// ```
    #[must_use]
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
    ///
    /// More accurate than `lonlat_to_vec(center(cell))` near the poles, and cheaper: it
    /// never goes through an inverse trigonometric function.
    ///
    /// ```
    /// let layer = realpix::nested::get(10);
    /// let v = layer.center_vec(layer.hash(1.0, 0.5));
    /// assert!((v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - 1.0).abs() < 1e-15);
    /// ```
    #[must_use]
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
    /// let layer = realpix::nested::get(1);
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
    pub(crate) fn center_loc(&self, cell: u64) -> Loc {
        self.assert_cell(cell);
        let xyf = self.base.nested2xyf(cell);
        self.base
            .xyf2loc(xyf.face, xyf.x as f64 + 0.5, xyf.y as f64 + 0.5)
    }

    /// The four corners of `cell` as unit vectors, in the order north, west, south, east.
    ///
    /// Consecutive corners bound one edge of the cell, and each edge is shared with one
    /// neighbour: north-west for `[0]..[1]`, south-west for `[1]..[2]`, south-east for
    /// `[2]..[3]`, and north-east for `[3]..[0]`. That correspondence is what lets you
    /// stroke a cell boundary once rather than twice, or find the cells along the border
    /// of a region.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    /// ```
    /// let layer = realpix::nested::get(8);
    /// let cell = layer.hash(1.0, 0.5);
    /// let [north, west, south, east] = layer.vertices(cell);
    /// // The corners bracket the centre in latitude.
    /// assert!(north[2] > layer.center_vec(cell)[2]);
    /// assert!(south[2] < layer.center_vec(cell)[2]);
    /// let _ = (west, east);
    /// ```
    #[must_use]
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
    /// ```
    /// let layer = realpix::nested::get(8);
    /// let corners = layer.vertices_lonlat(layer.hash(1.0, 0.5));
    /// assert_eq!(corners.len(), 4);
    /// assert!(corners.iter().all(|(lon, _)| (0.0..std::f64::consts::TAU).contains(lon)));
    /// ```
    #[must_use]
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
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let cell = layer.hash(1.0, 0.5);
    /// // Away from the base-cell corners every cell has all eight.
    /// assert_eq!(layer.neighbours(cell).iter().flatten().count(), 8);
    /// // Neighbourhood is symmetric.
    /// for n in layer.neighbours(cell).into_iter().flatten() {
    ///     assert!(layer.neighbours(n).contains(&Some(cell)));
    /// }
    /// ```
    #[must_use]
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
    ///
    /// ```
    /// use realpix::Direction;
    ///
    /// let layer = realpix::nested::get(6);
    /// let cell = layer.hash(1.0, 0.5);
    /// let north = layer.neighbour(cell, Direction::N).unwrap();
    /// assert!(layer.center(north).1 > layer.center(cell).1);
    /// ```
    #[inline]
    pub fn neighbour(&self, cell: u64, direction: Direction) -> Option<u64> {
        self.assert_cell(cell);
        let xyf = self.base.nested2xyf(cell);
        self.base
            .neighbour_xyf(xyf, direction.index())
            .map(|x| self.base.xyf2nested(x))
    }

    /// Calls `sink` with every cell at `edge_depth` that touches `cell` from outside it —
    /// the ring of cells immediately surrounding `cell`, at whatever resolution you ask
    /// for.
    ///
    /// This is the operation to expand a search outwards: the cells of `cell` itself hold
    /// the candidates you already have, and the external edge holds the ones just past the
    /// boundary, which a match near the edge of a cell would otherwise miss.
    ///
    /// `edge_depth` is absolute, like [`children_range`](Self::children_range), and must be
    /// at least this layer's depth. With `n = 2^(edge_depth - depth)` cells along the side
    /// of `cell`, the ring holds `4 * n + 4` cells, less one for each neighbour `cell`
    /// itself lacks: the four sides are always there, but a corner of the ring is missing
    /// wherever `cell` sits on a point where only three base cells meet. That costs the 24
    /// such cells of any layer one corner each, and every base cell two, since at depth 0
    /// each one touches two of those points. At `edge_depth == self.depth()` the ring is
    /// exactly [`neighbours`](Self::neighbours).
    ///
    /// Cells arrive once each, walking the ring, so nothing is allocated and there is
    /// nothing to deduplicate — but they do not arrive in index order. Use
    /// [`external_edge_cells`](Self::external_edge_cells) for that.
    ///
    /// # Panics
    /// Panics if `cell` is out of range, or if `edge_depth` is below this layer's depth or
    /// above [`MAX_DEPTH`].
    ///
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let mut ring = Vec::new();
    /// layer.external_edge(42, 8, |c| ring.push(c));
    /// // Four cells along each side of the 4x4 grid `42` splits into, plus four corners.
    /// assert_eq!(ring.len(), 4 * 4 + 4);
    /// ```
    pub fn external_edge<F: FnMut(u64)>(&self, cell: u64, edge_depth: u8, mut sink: F) {
        self.assert_cell(cell);
        assert!(
            edge_depth >= self.depth() && edge_depth <= MAX_DEPTH,
            "edge depth must be between this layer's depth and MAX_DEPTH"
        );
        let delta = edge_depth - self.depth();
        let n = 1u32 << delta;
        let deep = Base::new(edge_depth);
        let parent = self.base.nested2xyf(cell);
        let (x0, y0, face) = (parent.x << delta, parent.y << delta, parent.face);

        // Step outwards from the cell on the boundary of `cell` that faces each position
        // of the surrounding ring. Every position is reached from exactly one of them, so
        // the walk emits each neighbour once. `neighbour_xyf` handles the base-cell
        // crossings, including the corners where the eighth neighbour does not exist.
        let mut emit = |i: u32, j: u32, direction: Direction| {
            let inside = Xyf {
                x: x0 + i,
                y: y0 + j,
                face,
            };
            if let Some(outside) = deep.neighbour_xyf(inside, direction.index()) {
                sink(deep.xyf2nested(outside));
            }
        };

        // Once around the ring: corner, side, corner, side, ...
        emit(0, 0, Direction::S);
        for j in 0..n {
            emit(0, j, Direction::SW);
        }
        emit(0, n - 1, Direction::W);
        for i in 0..n {
            emit(i, n - 1, Direction::NW);
        }
        emit(n - 1, n - 1, Direction::N);
        for j in (0..n).rev() {
            emit(n - 1, j, Direction::NE);
        }
        emit(n - 1, 0, Direction::E);
        for i in (0..n).rev() {
            emit(i, 0, Direction::SE);
        }
    }

    /// [`external_edge`](Self::external_edge) collected into a `Vec`, sorted by index.
    ///
    /// Sorted and duplicate-free, so it can be binary-searched or merged against a
    /// catalogue sorted by NESTED index at `edge_depth`.
    ///
    /// # Panics
    /// Panics if `cell` is out of range, or if `edge_depth` is below this layer's depth or
    /// above [`MAX_DEPTH`].
    ///
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let cell = layer.hash(1.0, 0.5);
    /// // At its own depth the ring is the eight neighbours.
    /// assert_eq!(layer.external_edge_cells(cell, 6).len(), 8);
    /// // Nothing on the ring is inside the cell.
    /// let inside = layer.children_range(cell, 8);
    /// assert!(layer.external_edge_cells(cell, 8).iter().all(|c| !inside.contains(c)));
    /// ```
    #[must_use]
    #[cfg(feature = "alloc")]
    pub fn external_edge_cells(&self, cell: u64, edge_depth: u8) -> alloc::vec::Vec<u64> {
        let mut out = alloc::vec::Vec::new();
        self.external_edge(cell, edge_depth, |c| out.push(c));
        out.sort_unstable();
        out
    }

    // -------------------------------------------------------------------- hierarchy

    /// The ancestor of `cell` at `parent_depth`.
    ///
    /// # Panics
    /// Panics if `parent_depth > self.depth()` or if `cell` is out of range.
    ///
    /// ```
    /// let layer = realpix::nested::get(10);
    /// let cell = layer.hash(1.0, 0.5);
    /// // The base cell a position falls in, whatever depth you started from.
    /// assert_eq!(layer.parent(cell, 0), realpix::nested::get(0).hash(1.0, 0.5));
    /// assert_eq!(layer.parent(cell, 10), cell);
    /// ```
    #[must_use]
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
    ///
    /// ```
    /// let layer = realpix::nested::get(0);
    /// assert_eq!(layer.children(3), [12, 13, 14, 15]);
    /// // Every child has the cell as its parent.
    /// let below = realpix::nested::get(1);
    /// assert!(layer.children(3).iter().all(|c| below.parent(*c, 0) == 3));
    /// ```
    #[must_use]
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
    ///
    /// ```
    /// let layer = realpix::nested::get(0);
    /// // Two depths down, a base cell holds 16 contiguous cells.
    /// assert_eq!(layer.children_range(3, 2), 48..64);
    /// assert_eq!(layer.children_range(3, 0), 3..4);
    /// ```
    #[must_use]
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
    ///
    /// ```
    /// let layer = realpix::nested::get(8);
    /// let cell = layer.hash(1.0, 0.5);
    /// let ring = layer.to_ring(cell);
    /// // The same cell, so the same centre, and it converts back.
    /// assert_eq!(realpix::ring::get(8).center(ring), layer.center(cell));
    /// assert_eq!(realpix::ring::get(8).to_nested(ring), cell);
    /// ```
    #[must_use]
    #[inline]
    pub fn to_ring(&self, cell: u64) -> u64 {
        self.assert_cell(cell);
        self.base.xyf2ring(self.base.nested2xyf(cell))
    }

    /// Iterates over every cell of this layer, in index order.
    ///
    /// ```
    /// let layer = realpix::nested::get(1);
    /// assert_eq!(layer.iter(), 0..48);
    /// assert_eq!(layer.iter().count() as u64, layer.n_hash());
    /// ```
    #[must_use]
    #[inline]
    pub fn iter(&self) -> Range<u64> {
        0..self.base.n_hash
    }
}
