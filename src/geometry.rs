//! Geometric facts about HEALPix cells, independent of how they are numbered.

use crate::depth::MAX_DEPTH;

pub(crate) const N_DEPTHS: usize = MAX_DEPTH as usize + 1;

/// Upper bound, in radians, on the angular distance from a cell centre to any point on
/// that cell's boundary, at `depth`.
///
/// A cell is not a disc, so it has no single radius; this is the radius of the smallest
/// disc centred on a cell that is guaranteed to contain the whole cell, at any position on
/// the sphere. It is what you need to widen a search by so that a query cannot miss a
/// source lying anywhere in a cell it touches, and it is the bound both cone searches
/// prune with.
///
/// The value is an over-estimate and only ever safe to use as one: too large a bound makes
/// a search return extra cells, never miss one.
///
/// # Panics
/// Panics if `depth > `[`MAX_DEPTH`].
///
/// ```
/// let layer = realpix::nested::get(6);
/// let bound = realpix::max_center_to_vertex(6);
/// let cell = layer.hash(1.0, 0.5);
///
/// // Every corner of a cell is within the bound of its centre.
/// let centre = layer.center_vec(cell);
/// for corner in layer.vertices(cell) {
///     assert!(realpix::angular_distance(centre, corner) <= bound);
/// }
///
/// // It halves with each level, as cells do.
/// assert!(realpix::max_center_to_vertex(7) < bound);
/// ```
#[must_use]
#[inline]
pub const fn max_center_to_vertex(depth: u8) -> f64 {
    assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH");
    TABLE[depth as usize]
}

/// The bound at every depth, for code that walks depths in a loop and would otherwise pay
/// the range check each time.
///
/// Deliberately not public: it is a lookup table because that happens to be the cheapest
/// way to serve [`max_center_to_vertex`], not because the shape is part of the contract.
pub(crate) const TABLE: [f64; N_DEPTHS] = build();

/// Entries up to depth 9 are the exhaustively measured maxima over every cell of the
/// layer, plus 1%; deeper entries use the limit of `nside * max`, which the measured values
/// approach from below (`1.0690`), rounded up to `1.08`. `tests/geometry.rs` re-derives and
/// re-checks the whole table, and `examples/measure_cell_radius.rs` is what measured it.
const fn build() -> [f64; N_DEPTHS] {
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
