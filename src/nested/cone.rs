//! Cone (disc) search over the NESTED scheme.
//!
//! The search descends the quad-tree from the 12 base cells, testing each cell's centre
//! against the cone with two precomputed cosines per depth: one that proves the cell is
//! wholly outside, one that proves it is wholly inside. A cell that is wholly inside is
//! emitted as a single index range instead of being descended into, which is what makes
//! the result cheap to intersect with a catalogue sorted by NESTED index.

use super::Layer;
use crate::Vec3;
use crate::base::Base;
use crate::depth::MAX_DEPTH;
use crate::geometry::{N_DEPTHS, TABLE as MAX_CENTER_TO_VERTEX};
use crate::math::{PI, cos, dot};
use crate::merge::Merger;
use core::ops::Range;

/// How many levels below the requested depth the boundary test may be refined.
///
/// Cells straddling the cone boundary are subdivided this far before being accepted, which
/// shrinks the uncertain band by `2^REFINE` and trims most of the cells that come near the
/// cone without touching it.
///
/// Refinement costs up to `4^REFINE` extra centre evaluations per straddling cell, so it
/// is only worth it while the disc is small: once the disc spans many cells the boundary
/// band is already a thin fraction of the result, and the extra work would dominate the
/// query. [`refine_levels`] picks the amount from that ratio.
const REFINE: u8 = 2;

/// Levels of boundary refinement to use for a disc of `radius` at `depth`.
///
/// `radius / rho` is how many cells the disc spans; the thresholds keep the refinement
/// cost bounded while the relative over-inclusion is largest.
#[inline]
fn refine_levels(radius: f64, rho: f64) -> u8 {
    let cells_across = radius / rho;
    if cells_across < 4.0 {
        // A handful of cells in total: refine hard, it costs almost nothing.
        REFINE
    } else if cells_across < 64.0 {
        1
    } else {
        // The boundary band is already a thin shell around a large result.
        0
    }
}

/// Everything the descent needs, built once per query.
struct Cone {
    /// `sin(lat)` of the cone axis, i.e. `cos` of its colatitude.
    center_z: f64,
    /// `sin` of the cone axis's colatitude.
    center_sth: f64,
    /// Longitude of the cone axis, in radians.
    center_phi: f64,
    /// `cos(radius + rho[k])`: a cell at depth `k` whose centre has a smaller dot product
    /// than this cannot touch the cone.
    cos_outer: [f64; N_DEPTHS],
    /// `cos(radius - rho[k])`: a cell at depth `k` whose centre has a larger dot product
    /// than this is entirely inside the cone.
    cos_inner: [f64; N_DEPTHS],
    bases: [Base; N_DEPTHS],
    depth: u8,
    refine_depth: u8,
}

impl Cone {
    /// The dot product of a cell centre with the cone axis, or `None` when the cell was
    /// rejected on latitude alone.
    ///
    /// Working from the cell's own `(z, phi)` rather than building its unit vector saves
    /// the `sin` half of a `sin_cos`, since
    /// `dot = z*z0 + sin(theta)*sin(theta0)*cos(phi - phi0)`. It also puts a free and
    /// exact latitude test in the way: `cos(phi - phi0)` is at most 1, so the first term
    /// plus the whole second term is the largest the dot product could reach at *any*
    /// longitude — and that bound is `cos(theta - theta0)`. A cell failing it is outside
    /// the cone's latitude band altogether and never pays for the `cos`.
    #[inline]
    fn cell_dot(&self, depth: u8, cell: u64) -> Option<f64> {
        let base = &self.bases[depth as usize];
        let xyf = base.nested2xyf(cell);
        let (x, y) = (xyf.x as f64 + 0.5, xyf.y as f64 + 0.5);
        let loc = base.xyf2loc(xyf.face, x, y);

        let across = loc.z * self.center_z;
        let along = loc.sth * self.center_sth;
        if across + along < self.cos_outer[depth as usize] {
            return None;
        }
        Some(across + along * cos(loc.phi - self.center_phi))
    }

    /// Whether any part of `cell` might touch the cone, refined `REFINE` levels down.
    fn may_touch(&self, depth: u8, cell: u64) -> bool {
        let Some(d) = self.cell_dot(depth, cell) else {
            return false;
        };
        if d < self.cos_outer[depth as usize] {
            return false;
        }
        if d >= self.cos_inner[depth as usize] || depth == self.refine_depth {
            return true;
        }
        let first = cell << 2;
        (first..first + 4).any(|child| self.may_touch(depth + 1, child))
    }

    fn descend<F: FnMut(Range<u64>)>(&self, depth: u8, cell: u64, out: &mut Merger<F>) {
        let Some(d) = self.cell_dot(depth, cell) else {
            return; // Too far north or south to touch the cone.
        };
        if d < self.cos_outer[depth as usize] {
            return; // Disjoint from the cone.
        }

        let shift = (self.depth - depth) << 1;
        if d >= self.cos_inner[depth as usize] {
            // Entirely inside: emit the whole subtree as one range.
            out.push((cell << shift)..((cell + 1) << shift));
            return;
        }
        if depth == self.depth {
            if self.may_touch(depth, cell) {
                out.push(cell..cell + 1);
            }
            return;
        }

        let first = cell << 2;
        for child in first..first + 4 {
            self.descend(depth + 1, child, out);
        }
    }
}

impl Layer {
    /// Calls `sink` with the NESTED index ranges covering the cone of angular `radius`
    /// (radians) around the direction `center`.
    ///
    /// The ranges are sorted, disjoint, non-adjacent and half-open, and together they
    /// cover *at least* every cell that intersects the cone: the result is inclusive, so
    /// a few cells that merely come close to the boundary may be included, but no cell
    /// that touches the cone is ever missed. Nothing is allocated.
    ///
    /// A catalogue sorted by NESTED index at this depth can be searched by slicing it with
    /// these ranges.
    ///
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let center = realpix::lonlat_to_vec(1.0, 0.5);
    /// let mut cells = 0u64;
    /// layer.cone_coverage(center, 0.05, |r| cells += r.end - r.start);
    /// assert!(cells > 0);
    /// ```
    pub fn cone_coverage<F: FnMut(Range<u64>)>(&self, center: Vec3, radius: f64, sink: F) {
        let mut merger = Merger::new(sink);

        // The axis must be normalisable. A zero vector has finite components but no
        // direction, and normalising it would put a NaN into every comparison below —
        // where NaN answers "no" to both the reject and the accept test, so nothing would
        // prune and the descent would walk the entire layer.
        let len2 = dot(&center, &center);
        if !radius.is_finite() || radius < 0.0 || !len2.is_finite() || len2 <= 0.0 {
            merger.flush();
            return;
        }
        if radius >= PI {
            merger.push(0..self.n_hash());
            merger.flush();
            return;
        }

        let norm = 1.0 / crate::math::sqrt(len2);
        let center_z = (center[2] * norm).clamp(-1.0, 1.0);
        let mut cone = Cone {
            center_z,
            center_sth: crate::math::sqrt((1.0 - center_z) * (1.0 + center_z)),
            center_phi: crate::math::safe_atan2(center[1], center[0]),
            cos_outer: [-2.0; N_DEPTHS],
            cos_inner: [2.0; N_DEPTHS],
            bases: [Base::new(0); N_DEPTHS],
            depth: self.depth(),
            refine_depth: (self.depth()
                + refine_levels(radius, MAX_CENTER_TO_VERTEX[self.depth() as usize]))
            .min(MAX_DEPTH),
        };
        for (k, &rho) in MAX_CENTER_TO_VERTEX
            .iter()
            .enumerate()
            .take(cone.refine_depth as usize + 1)
        {
            cone.bases[k] = Base::new(k as u8);
            if radius + rho < PI {
                cone.cos_outer[k] = cos(radius + rho);
            }
            if radius > rho {
                cone.cos_inner[k] = cos(radius - rho);
            }
        }

        for base_cell in 0..12u64 {
            cone.descend(0, base_cell, &mut merger);
        }
        merger.flush();
    }

    /// [`cone_coverage`](Self::cone_coverage) taking the cone centre as `(lon, lat)` in radians.
    ///
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let mut from_lonlat = Vec::new();
    /// layer.cone_coverage_lonlat(1.0, 0.5, 0.05, |r| from_lonlat.push(r));
    /// assert_eq!(from_lonlat, layer.cone_coverage_ranges(realpix::lonlat_to_vec(1.0, 0.5), 0.05));
    /// ```
    pub fn cone_coverage_lonlat<F: FnMut(Range<u64>)>(
        &self,
        lon: f64,
        lat: f64,
        radius: f64,
        sink: F,
    ) {
        self.cone_coverage(crate::tangent::lonlat_to_vec(lon, lat), radius, sink)
    }

    /// [`cone_coverage`](Self::cone_coverage), collected into a `Vec`.
    ///
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let ranges = layer.cone_coverage_ranges(realpix::lonlat_to_vec(1.0, 0.5), 0.05);
    /// // Sorted, disjoint and non-adjacent.
    /// assert!(ranges.windows(2).all(|w| w[0].end < w[1].start));
    /// ```
    #[must_use]
    #[cfg(feature = "alloc")]
    pub fn cone_coverage_ranges(&self, center: Vec3, radius: f64) -> alloc::vec::Vec<Range<u64>> {
        let mut out = alloc::vec::Vec::new();
        self.cone_coverage(center, radius, |r| out.push(r));
        out
    }

    /// Every cell covered by [`cone_coverage`](Self::cone_coverage), in index order.
    ///
    /// Sorted, so it can be binary-searched. Prefer
    /// [`cone_coverage_ranges`](Self::cone_coverage_ranges) for a large disc, where the
    /// ranges are far more compact than the cells they expand to.
    ///
    /// ```
    /// let layer = realpix::nested::get(6);
    /// let center = realpix::lonlat_to_vec(1.0, 0.5);
    /// let cells = layer.cone_coverage_cells(center, 0.05);
    /// // The cone's own cell is always covered.
    /// assert!(cells.binary_search(&layer.hash_vec(center)).is_ok());
    /// ```
    #[must_use]
    #[cfg(feature = "alloc")]
    pub fn cone_coverage_cells(&self, center: Vec3, radius: f64) -> alloc::vec::Vec<u64> {
        let mut out = alloc::vec::Vec::new();
        self.cone_coverage(center, radius, |r| out.extend(r));
        out
    }
}
