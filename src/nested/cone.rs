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
use crate::math::{PI, cos, dot};
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

const N_DEPTHS: usize = MAX_DEPTH as usize + 1;

/// Upper bound, in radians, on the angular distance from a cell centre to any point on
/// that cell's boundary, indexed by depth.
///
/// Entries up to depth 9 are the exhaustively measured maxima over every cell of the
/// layer, plus 1%; deeper entries use the limit of `nside * max`, which the measured
/// values approach from below (`1.0690`), rounded up to `1.08`. `tests/geometry.rs`
/// re-derives and checks the whole table.
///
/// The bound is only ever used conservatively: over-estimating it can make a cone search
/// return extra cells, never miss one.
pub const MAX_CENTER_TO_VERTEX: [f64; N_DEPTHS] = center_to_vertex_table();

const fn center_to_vertex_table() -> [f64; N_DEPTHS] {
    let mut t = [0.0f64; N_DEPTHS];
    t[0] = 0.849_479_357_273_609_8;
    t[1] = 0.486_275_076_190_827_7;
    t[2] = 0.256_876_738_598_488_5;
    t[3] = 0.131_729_689_889_850_1;
    t[4] = 0.066_674_909_046_839_1;
    t[5] = 0.033_538_733_479_886_7;
    t[6] = 0.016_819_555_463_525_1;
    t[7] = 0.008_422_310_052_121_9;
    t[8] = 0.004_214_286_344_673_0;
    t[9] = 0.002_107_925_787_852_9;
    let mut d = 10;
    while d < N_DEPTHS {
        t[d] = 1.08 / ((1u64 << d) as f64);
        d += 1;
    }
    t
}

/// Merges the ranges emitted by the descent before handing them to the caller.
///
/// The descent visits cells in increasing index order, so a new range is always either
/// contiguous with (or contained in) the pending one, or entirely beyond it.
struct Merger<F: FnMut(Range<u64>)> {
    pending: Option<Range<u64>>,
    sink: F,
}

impl<F: FnMut(Range<u64>)> Merger<F> {
    #[inline]
    fn push(&mut self, range: Range<u64>) {
        if let Some(pending) = self.pending.as_mut()
            && range.start <= pending.end
        {
            pending.end = pending.end.max(range.end);
            return;
        }
        if let Some(previous) = self.pending.replace(range) {
            (self.sink)(previous);
        }
    }

    #[inline]
    fn flush(mut self) {
        if let Some(pending) = self.pending.take() {
            (self.sink)(pending);
        }
    }
}

/// Everything the descent needs, built once per query.
struct Cone {
    center: Vec3,
    /// `sin(lat)` of the cone axis.
    center_z: f64,
    /// `radius + rho[k]`, as a bound on `|z_cell - z_axis|`. Since
    /// `|sin a - sin b| <= |a - b| <= angular distance`, a cell whose centre differs from
    /// the axis by more than this in `z` cannot touch the cone — and that test costs no
    /// trigonometry at all, so it runs first.
    max_dz: [f64; N_DEPTHS],
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
    /// rejected by the cheap latitude test.
    #[inline]
    fn cell_dot(&self, depth: u8, cell: u64) -> Option<f64> {
        let base = &self.bases[depth as usize];
        let xyf = base.nested2xyf(cell);
        let (x, y) = (xyf.x as f64 + 0.5, xyf.y as f64 + 0.5);
        let z = base.xyf2z(xyf.face, x, y);
        if (z - self.center_z).abs() > self.max_dz[depth as usize] {
            return None;
        }
        Some(dot(&base.xyf2loc(xyf.face, x, y).to_vec(), &self.center))
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
        let mut merger = Merger {
            pending: None,
            sink,
        };

        if !radius.is_finite() || radius < 0.0 || !center.iter().all(|c| c.is_finite()) {
            merger.flush();
            return;
        }
        if radius >= PI {
            merger.push(0..self.n_hash());
            merger.flush();
            return;
        }

        let norm = 1.0 / crate::math::sqrt(dot(&center, &center));
        let mut cone = Cone {
            center: [center[0] * norm, center[1] * norm, center[2] * norm],
            center_z: center[2] * norm,
            max_dz: [2.0; N_DEPTHS],
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
            cone.max_dz[k] = radius + rho;
        }

        for base_cell in 0..12u64 {
            cone.descend(0, base_cell, &mut merger);
        }
        merger.flush();
    }

    /// [`cone_coverage`](Self::cone_coverage) taking the cone centre as `(lon, lat)` in radians.
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
    #[cfg(feature = "alloc")]
    pub fn cone_coverage_ranges(&self, center: Vec3, radius: f64) -> alloc::vec::Vec<Range<u64>> {
        let mut out = alloc::vec::Vec::new();
        self.cone_coverage(center, radius, |r| out.push(r));
        out
    }

    /// Every cell covered by [`cone_coverage`](Self::cone_coverage), in index order.
    #[cfg(feature = "alloc")]
    pub fn cone_coverage_cells(&self, center: Vec3, radius: f64) -> alloc::vec::Vec<u64> {
        let mut out = alloc::vec::Vec::new();
        self.cone_coverage(center, radius, |r| out.extend(r));
        out
    }
}
