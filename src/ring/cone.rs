//! Cone (disc) search over the RING scheme.
//!
//! RING indices are ordered by decreasing latitude, so a disc does not occupy one
//! contiguous run the way it can in NESTED. What it does occupy is one contiguous arc of
//! each iso-latitude ring it reaches — two, when the arc wraps past longitude zero — and
//! that is exactly what this search walks: the rings whose latitude band the disc reaches,
//! and for each of them the half-width in longitude that the disc subtends.
//!
//! The cost is therefore proportional to the number of rings the disc spans, with a
//! constant amount of trigonometry per ring and no descent and no allocation.

use super::Layer;
use crate::Vec3;
use crate::base::Base;
use crate::geometry::TABLE as MAX_CENTER_TO_VERTEX;
use crate::math::{PI, TAU, acos, cos, dot, fmodulo, safe_atan2, sqrt};
use crate::merge::Merger;
use core::ops::Range;

/// The geometry of one iso-latitude ring, `jr` counted from the north pole and 1-based.
struct Ring {
    /// `sin(lat)` of every cell centre on the ring.
    z: f64,
    /// Cells per quadrant; the ring holds `4 * nr` cells.
    nr: i64,
    /// RING index of the ring's first cell.
    first: u64,
    /// Longitude of the `j`-th cell of the ring is `(j + shift) * (π/2) / nr`.
    shift: f64,
}

impl Ring {
    /// The geometry of ring `jr`, which must be in `1..4 * nside`.
    #[inline]
    fn new(base: &Base, jr: i64) -> Self {
        let nside = base.nside as i64;
        if jr < nside {
            // North polar cap: ring `jr` holds `4 * jr` cells, offset by half a cell.
            Ring {
                z: 1.0 - (jr * jr) as f64 * base.fact2,
                nr: jr,
                first: (2 * jr * (jr - 1)) as u64,
                shift: 0.5,
            }
        } else if jr <= 3 * nside {
            // Equatorial belt: every ring holds `4 * nside` cells, every other one offset
            // by half a cell. `kshift` matches the reference implementation's `fodd`.
            let kshift = (jr - nside) & 1;
            Ring {
                z: (2 * nside - jr) as f64 * base.fact1,
                nr: nside,
                first: base.ncap + ((jr - nside) * (nside << 2)) as u64,
                shift: if kshift == 1 { 0.0 } else { 0.5 },
            }
        } else {
            // South polar cap: the north cap mirrored.
            let nr = (nside << 2) - jr;
            Ring {
                z: (nr * nr) as f64 * base.fact2 - 1.0,
                nr,
                first: base.n_hash - (2 * nr * (nr + 1)) as u64,
                shift: 0.5,
            }
        }
    }
}

/// The index of the ring at `z`, to within a ring or so.
///
/// Only ever used to bracket the rings a disc reaches, and the bracket is padded either
/// side before use, so the rounding here does not have to be exact — a ring wrongly
/// included is culled by its own longitude test, at the cost of one `acos`.
#[inline]
fn ring_at_z(base: &Base, z: f64) -> i64 {
    let nside = base.nside as i64;
    if z > crate::math::TRANSITION_Z {
        sqrt((1.0 - z) / base.fact2) as i64
    } else if z >= -crate::math::TRANSITION_Z {
        2 * nside - (z / base.fact1) as i64
    } else {
        (nside << 2) - sqrt((1.0 + z) / base.fact2) as i64
    }
}

impl Layer {
    /// Calls `sink` with the RING index ranges covering the cone of angular `radius`
    /// (radians) around the direction `center`.
    ///
    /// The ranges are sorted, disjoint, non-adjacent and half-open, and together they
    /// cover *at least* every cell that intersects the cone: the result is inclusive, so
    /// a few cells that merely come close to the boundary may be included, but no cell
    /// that touches the cone is ever missed. Nothing is allocated.
    ///
    /// A catalogue sorted by RING index at this depth can be searched by slicing it with
    /// these ranges.
    ///
    /// This covers a disc about as tightly, and in about as many ranges, as the
    /// [`nested`](crate::nested) search does — it over-includes somewhat more on the
    /// smallest discs, having no boundary refinement to fall back on, and converges to the
    /// same cells as the disc grows. It is also markedly cheaper, since it visits only the
    /// rings the disc reaches rather than descending from the whole sphere. Choose the
    /// scheme your catalogue is already sorted in.
    ///
    /// ```
    /// let layer = realpix::ring::get(6);
    /// let center = realpix::lonlat_to_vec(1.0, 0.5);
    /// let mut cells = 0u64;
    /// layer.cone_coverage(center, 0.05, |r| cells += r.end - r.start);
    /// assert!(cells > 0);
    /// ```
    pub fn cone_coverage<F: FnMut(Range<u64>)>(&self, center: Vec3, radius: f64, sink: F) {
        let mut merger = Merger::new(sink);

        let len2 = dot(&center, &center);
        if !radius.is_finite() || radius < 0.0 || !len2.is_finite() || len2 <= 0.0 {
            merger.flush();
            return;
        }

        // Selecting every cell whose *centre* lies within `radius + rho` of the axis is
        // what makes the result inclusive: a cell that touches the cone has some point at
        // most `radius` from the axis and at most `rho` from its own centre.
        let rho = MAX_CENTER_TO_VERTEX[self.depth() as usize];
        let radius = radius + rho;
        if radius >= PI {
            merger.push(0..self.n_hash());
            merger.flush();
            return;
        }

        let norm = 1.0 / sqrt(len2);
        let (cx, cy, cz) = (
            center[0] * norm,
            center[1] * norm,
            (center[2] * norm).clamp(-1.0, 1.0),
        );
        let lon = fmodulo(safe_atan2(cy, cx), TAU);
        let cos_radius = cos(radius);
        // `sin` of the axis colatitude; zero exactly at a pole, where every cell on a
        // ring is equidistant from the axis and the longitude test degenerates.
        let sin_axis = sqrt((1.0 - cz) * (1.0 + cz));

        let base = &self.base;
        let last_ring = ((base.nside as i64) << 2) - 1;
        // `|sin a - sin b| <= |a - b|`, so no ring within the disc differs from the axis
        // by more than `radius` in `z`. The two-ring padding absorbs the rounding in
        // `ring_at_z`; rings that turn out to be too far north or south fall out of the
        // longitude test below.
        let first_ring = (ring_at_z(base, (cz + radius).min(1.0)) - 2).max(1);
        let final_ring = (ring_at_z(base, (cz - radius).max(-1.0)) + 2).min(last_ring);

        for jr in first_ring..=final_ring {
            let ring = Ring::new(base, jr);
            let sin_ring = sqrt((1.0 - ring.z) * (1.0 + ring.z));
            let denominator = sin_ring * sin_axis;

            // Half-width in longitude of the arc this ring shares with the disc, from
            // `cos(radius) = z·z₀ + sin·sin₀·cos(Δlon)`.
            let d_lon = if denominator <= 0.0 {
                // The axis is at a pole: the whole ring is in or out together.
                if ring.z * cz >= cos_radius {
                    PI
                } else {
                    continue;
                }
            } else {
                let cos_d_lon = (cos_radius - ring.z * cz) / denominator;
                if cos_d_lon > 1.0 {
                    continue; // The disc does not reach this ring.
                } else if cos_d_lon <= -1.0 {
                    PI
                } else {
                    acos(cos_d_lon)
                }
            };

            // Longitude to in-ring index: `lon = (j + shift) * (π/2) / nr`.
            let per_cell = (ring.nr as f64) * crate::math::INV_HALF_PI;
            let low = floor_i64((lon - d_lon) * per_cell - ring.shift) + 1;
            let high = floor_i64((lon + d_lon) * per_cell - ring.shift);
            if high < low {
                continue;
            }

            let in_ring = ring.nr << 2;
            let count = high - low + 1;
            if count >= in_ring {
                merger.push(ring.first..ring.first + in_ring as u64);
                continue;
            }

            let start = low.rem_euclid(in_ring);
            let end = start + count;
            if end > in_ring {
                // The arc wraps past longitude zero: emit the low part first, so that the
                // ranges still reach the merger in increasing order.
                merger.push(ring.first..ring.first + (end - in_ring) as u64);
                merger.push(ring.first + start as u64..ring.first + in_ring as u64);
            } else {
                merger.push(ring.first + start as u64..ring.first + end as u64);
            }
        }

        merger.flush();
    }

    /// [`cone_coverage`](Self::cone_coverage) taking the cone centre as `(lon, lat)` in radians.
    ///
    /// ```
    /// let layer = realpix::ring::get(6);
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
    /// let layer = realpix::ring::get(6);
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
    /// let layer = realpix::ring::get(6);
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

/// `floor` as an `i64`, without pulling a rounding routine in from the float backend. The
/// argument is bounded by a few times `nside`, so it always fits.
#[inline]
fn floor_i64(x: f64) -> i64 {
    let truncated = x as i64;
    if x < truncated as f64 {
        truncated - 1
    } else {
        truncated
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The whole search rests on `Ring` describing the same cells that
    /// [`Layer::center`](super::Layer::center) does. Check it cell by cell.
    #[test]
    fn ring_geometry_matches_the_cell_centres() {
        for depth in 0..=5u8 {
            let layer = super::super::get(depth);
            let base = &layer.base;
            let nside = base.nside as i64;
            let mut expected_first = 0u64;

            for jr in 1..(nside << 2) {
                let ring = Ring::new(base, jr);
                assert_eq!(
                    ring.first, expected_first,
                    "depth {depth}, ring {jr}: first"
                );
                expected_first += (ring.nr << 2) as u64;

                for j in 0..(ring.nr << 2) {
                    let cell = ring.first + j as u64;
                    let (lon, lat) = layer.center(cell);
                    let expected_lon =
                        (j as f64 + ring.shift) * crate::math::FRAC_PI_2 / ring.nr as f64;
                    assert!(
                        (crate::math::sin(lat) - ring.z).abs() < 1e-12,
                        "depth {depth}, ring {jr}, cell {cell}: z {} against {}",
                        ring.z,
                        crate::math::sin(lat)
                    );
                    assert!(
                        (lon - expected_lon).abs() < 1e-12,
                        "depth {depth}, ring {jr}, cell {cell}: lon {lon} against {expected_lon}"
                    );
                }
            }
            // The rings together must account for every cell of the layer.
            assert_eq!(expected_first, base.n_hash, "depth {depth}: cell count");
        }
    }

    /// `ring_at_z` is only ever used padded by two rings, so that is the tolerance it has
    /// to meet — at every ring of every region, cap and belt alike.
    #[test]
    fn ring_at_z_lands_within_its_padding() {
        for depth in [0, 1, 5, 10, 20, 29] {
            let base = crate::base::Base::new(depth);
            let last = ((base.nside as i64) << 2) - 1;
            // Every ring at the low depths, a sample of them at the high ones.
            let step = (last / 500).max(1);
            for jr in (1..=last).step_by(step as usize) {
                let got = ring_at_z(&base, Ring::new(&base, jr).z);
                assert!(
                    (got - jr).abs() <= 2,
                    "depth {depth}, ring {jr}: ring_at_z gave {got}"
                );
            }
        }
    }
}
