//! Multi-order coverage maps: sets of cells drawn from any mix of depths.
//!
//! A [`Moc`] is an arbitrary region of the sphere described as a set of HEALPix cells. It
//! is what you get when you want to keep a coverage around — the sky a survey reached, the
//! union of the fields solved so far, the part of a catalogue worth searching — and then
//! combine it with another one.
//!
//! # Representation
//!
//! Every cell is stored as the range of depth-[`MAX_DEPTH`] cells it contains, so a
//! coverage is one sorted list of disjoint, non-adjacent ranges regardless of the depths it
//! was built from. That makes the set operations plain linear merges over two sorted lists,
//! exact at every depth, and it makes equality meaningful: two `Moc`s covering the same
//! region are equal whatever mix of depths each was built from.
//!
//! The multi-order view is still there when you need it — [`Moc::cells`] decomposes the
//! coverage back into the largest cells that tile it, and [`Moc::uniq_cells`] gives those
//! as [NUNIQ](crate::to_uniq) values, which is how MOCs are normally serialised.
//!
//! Cells are always NESTED: the hierarchy the representation relies on is the one NESTED
//! numbering encodes. Map a RING index across with
//! [`ring::Layer::to_nested`](crate::ring::Layer::to_nested) first.
//!
//! ```
//! use realpix::moc::Moc;
//!
//! // The sky two exposures reached, and the overlap between them.
//! let first = Moc::from_cone(8, realpix::lonlat_to_vec(1.0, 0.5), 0.02);
//! let second = Moc::from_cone(8, realpix::lonlat_to_vec(1.01, 0.5), 0.02);
//!
//! let both = &first | &second;
//! let overlap = &first & &second;
//! assert!(overlap.area() > 0.0);
//! assert!(both.area() > first.area());
//!
//! // Slice a catalogue sorted by NESTED index at depth 10 with the union.
//! for range in both.ranges_at(10) {
//!     let _ = range;
//! }
//! ```

use alloc::vec::Vec;
use core::ops::{BitAnd, BitOr, BitXor, Not, Range, Sub};

use crate::depth::{MAX_DEPTH, n_hash};

/// Number of depth-`MAX_DEPTH` cells inside one cell at `depth`, as a bit shift.
#[inline]
const fn shift_at(depth: u8) -> u32 {
    assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH");
    ((MAX_DEPTH - depth) as u32) << 1
}

/// The number of depth-`MAX_DEPTH` cells on the whole sphere.
const FULL_SKY: u64 = n_hash(MAX_DEPTH);

/// A set of HEALPix cells, of any mix of depths, over which set algebra is exact.
///
/// See the [module documentation](self) for what it is for and how it is stored.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Moc {
    /// Sorted, disjoint, non-adjacent ranges of depth-`MAX_DEPTH` cell indices.
    ranges: Vec<Range<u64>>,
}

impl Moc {
    /// An empty coverage, covering nothing.
    ///
    /// The identity for [`union`](Self::union), so it is the value to fold a sequence of
    /// coverages into.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// assert!(Moc::new().is_empty());
    /// assert_eq!(Moc::new().area(), 0.0);
    /// assert_eq!(Moc::new(), Moc::default());
    /// ```
    #[inline]
    pub const fn new() -> Moc {
        Moc { ranges: Vec::new() }
    }

    /// The whole sphere.
    ///
    /// The identity for [`intersection`](Self::intersection), and the complement of
    /// [`new`](Self::new).
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// assert_eq!(Moc::all_sky().sky_fraction(), 1.0);
    /// assert_eq!(Moc::all_sky(), Moc::from_cells(0, 0..12));
    /// assert_eq!(Moc::all_sky().complement(), Moc::new());
    /// ```
    #[inline]
    pub fn all_sky() -> Moc {
        Moc {
            ranges: alloc::vec![0..FULL_SKY],
        }
    }

    /// A coverage of the given NESTED cells, all at `depth`.
    ///
    /// The cells may arrive in any order and may repeat.
    ///
    /// # Panics
    /// Panics if `depth > `[`MAX_DEPTH`] or a cell is out of range for `depth`.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// // Duplicates and order do not matter; the four children of one base cell are that
    /// // base cell.
    /// assert_eq!(Moc::from_cells(3, [5, 5, 5]), Moc::from_cells(3, [5]));
    /// assert_eq!(Moc::from_cells(1, [15, 12, 14, 13]), Moc::from_cells(0, [3]));
    /// ```
    pub fn from_cells<I: IntoIterator<Item = u64>>(depth: u8, cells: I) -> Moc {
        let shift = shift_at(depth);
        let limit = n_hash(depth);
        Moc::from_deep_ranges(cells.into_iter().map(|cell| {
            assert!(cell < limit, "cell {cell} out of range for depth {depth}");
            (cell << shift)..((cell + 1) << shift)
        }))
    }

    /// A coverage of the given ranges of NESTED cell indices, all at `depth`.
    ///
    /// This is what [`cone_coverage`](crate::nested::Layer::cone_coverage) hands you, so a
    /// cone goes straight in. Ranges may arrive in any order and may overlap.
    ///
    /// # Panics
    /// Panics if `depth > `[`MAX_DEPTH`] or a range reaches past the end of the layer.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let layer = realpix::nested::get(8);
    /// let center = realpix::lonlat_to_vec(1.0, 0.5);
    /// let cone = Moc::from_ranges(8, layer.cone_coverage_ranges(center, 0.05));
    /// assert!(cone.contains_vec(center));
    ///
    /// // Overlapping ranges normalise away.
    /// assert_eq!(
    ///     Moc::from_ranges(3, [0..10, 5..20]),
    ///     Moc::from_ranges(3, [0..20]),
    /// );
    /// ```
    pub fn from_ranges<I: IntoIterator<Item = Range<u64>>>(depth: u8, ranges: I) -> Moc {
        let shift = shift_at(depth);
        let limit = n_hash(depth);
        Moc::from_deep_ranges(ranges.into_iter().filter_map(move |r| {
            assert!(r.end <= limit, "range {r:?} out of range for depth {depth}");
            (r.start < r.end).then(|| (r.start << shift)..(r.end << shift))
        }))
    }

    /// A coverage of the given [NUNIQ](crate::to_uniq) cells, which may be at any mix of
    /// depths.
    ///
    /// # Panics
    /// Panics if a value is not a valid NUNIQ cell.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// // A mix of depths, as a serialised MOC holds it.
    /// let moc = Moc::from_uniq_cells([realpix::to_uniq(0, 3), realpix::to_uniq(2, 100)]);
    /// assert!(moc.contains(0, 3));
    /// assert!(moc.contains(2, 100));
    /// // And it round-trips.
    /// assert_eq!(Moc::from_uniq_cells(moc.uniq_cells()), moc);
    /// ```
    pub fn from_uniq_cells<I: IntoIterator<Item = u64>>(uniqs: I) -> Moc {
        Moc::from_deep_ranges(uniqs.into_iter().map(|uniq| {
            let (depth, cell) = crate::from_uniq(uniq);
            let shift = shift_at(depth);
            (cell << shift)..((cell + 1) << shift)
        }))
    }

    /// The cone of angular `radius` (radians) around `center`, at `depth`.
    ///
    /// A convenience for `Moc::from_ranges(depth, ...)` over
    /// [`cone_coverage`](crate::nested::Layer::cone_coverage). The cone search is
    /// inclusive, so the coverage is too.
    ///
    /// # Panics
    /// Panics if `depth > `[`MAX_DEPTH`].
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let center = realpix::lonlat_to_vec(1.0, 0.5);
    /// let cone = Moc::from_cone(10, center, 0.02);
    /// assert!(cone.contains_vec(center));
    /// // Inclusive, so it is at least as large as the disc itself.
    /// assert!(cone.area() >= 2.0 * std::f64::consts::PI * (1.0 - 0.02f64.cos()));
    /// ```
    pub fn from_cone(depth: u8, center: crate::Vec3, radius: f64) -> Moc {
        // The cone search already emits sorted, disjoint, non-adjacent ranges, so they can
        // be scaled to the deepest layer as they arrive and stored without a second pass.
        let shift = shift_at(depth);
        let mut ranges = Vec::new();
        crate::nested::get(depth).cone_coverage(center, radius, |r| {
            ranges.push((r.start << shift)..(r.end << shift))
        });
        Moc { ranges }
    }

    /// Normalises an arbitrary sequence of depth-`MAX_DEPTH` ranges: sort, then coalesce
    /// everything that touches or overlaps.
    fn from_deep_ranges<I: Iterator<Item = Range<u64>>>(ranges: I) -> Moc {
        let mut ranges: Vec<Range<u64>> = ranges.filter(|r| r.start < r.end).collect();
        ranges.sort_unstable_by_key(|r| r.start);
        let mut out: Vec<Range<u64>> = Vec::with_capacity(ranges.len());
        for r in ranges {
            push_coalesced(&mut out, r);
        }
        Moc { ranges: out }
    }

    /// Whether the coverage is empty.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// assert!(Moc::new().is_empty());
    /// assert!(!Moc::all_sky().is_empty());
    /// // Disjoint fields have an empty overlap.
    /// let a = Moc::from_cone(8, realpix::lonlat_to_vec(0.0, 0.0), 0.01);
    /// let b = Moc::from_cone(8, realpix::lonlat_to_vec(3.0, 0.0), 0.01);
    /// assert!((&a & &b).is_empty());
    /// ```
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    /// The area of the coverage in steradians.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let sphere = 4.0 * std::f64::consts::PI;
    /// assert!((Moc::all_sky().area() - sphere).abs() < 1e-12);
    /// // One base cell is a twelfth of the sky.
    /// assert!((Moc::from_cells(0, [7]).area() - sphere / 12.0).abs() < 1e-12);
    /// ```
    pub fn area(&self) -> f64 {
        let cells: u64 = self.ranges.iter().map(|r| r.end - r.start).sum();
        cells as f64 * (4.0 * core::f64::consts::PI / FULL_SKY as f64)
    }

    /// The fraction of the sphere the coverage takes up, in `0.0..=1.0`.
    ///
    /// [`area`](Self::area) divided by `4π`, for when a proportion reads better than a
    /// solid angle.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// assert_eq!(Moc::all_sky().sky_fraction(), 1.0);
    /// assert_eq!(Moc::new().sky_fraction(), 0.0);
    /// assert!((Moc::from_cells(0, [7]).sky_fraction() - 1.0 / 12.0).abs() < 1e-15);
    /// ```
    pub fn sky_fraction(&self) -> f64 {
        let cells: u64 = self.ranges.iter().map(|r| r.end - r.start).sum();
        cells as f64 / FULL_SKY as f64
    }

    /// Whether the NESTED `cell` at `depth` is entirely inside the coverage.
    ///
    /// # Panics
    /// Panics if `depth > `[`MAX_DEPTH`] or `cell` is out of range for `depth`.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let moc = Moc::from_cells(1, [12, 13]);
    /// assert!(moc.contains(1, 12));
    /// // Only *part* of base cell 3 is covered, so the whole cell is not.
    /// assert!(!moc.contains(0, 3));
    /// // But everything inside a covered cell is covered.
    /// assert!(moc.contains(2, 48));
    /// ```
    pub fn contains(&self, depth: u8, cell: u64) -> bool {
        assert!(
            cell < n_hash(depth),
            "cell {cell} out of range for depth {depth}"
        );
        let shift = shift_at(depth);
        let (start, end) = (cell << shift, (cell + 1) << shift);
        match self.range_covering(start) {
            Some(r) => r.end >= end,
            None => false,
        }
    }

    /// Whether the direction `v` falls inside the coverage.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let center = realpix::lonlat_to_vec(1.0, 0.5);
    /// let cone = Moc::from_cone(8, center, 0.02);
    /// assert!(cone.contains_vec(center));
    /// assert!(!cone.contains_vec(realpix::lonlat_to_vec(3.0, -0.5)));
    /// ```
    #[inline]
    pub fn contains_vec(&self, v: crate::Vec3) -> bool {
        self.range_covering(crate::nested::get(MAX_DEPTH).hash_vec(v))
            .is_some()
    }

    /// Whether the position falls inside the coverage, `lon` and `lat` in radians.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let cone = Moc::from_cone(8, realpix::lonlat_to_vec(1.0, 0.5), 0.02);
    /// assert!(cone.contains_lonlat(1.0, 0.5));
    /// assert!(!cone.contains_lonlat(3.0, -0.5));
    /// ```
    #[inline]
    pub fn contains_lonlat(&self, lon: f64, lat: f64) -> bool {
        self.range_covering(crate::nested::get(MAX_DEPTH).hash(lon, lat))
            .is_some()
    }

    /// The stored range containing the depth-`MAX_DEPTH` cell `deep`, if any.
    #[inline]
    fn range_covering(&self, deep: u64) -> Option<&Range<u64>> {
        // The ranges are sorted and disjoint, so at most one can contain `deep`: the last
        // one that starts at or before it.
        let i = self.ranges.partition_point(|r| r.start <= deep);
        self.ranges.get(i.checked_sub(1)?).filter(|r| deep < r.end)
    }

    /// The union of two coverages: everything in either one.
    ///
    /// Also available as `&a | &b`.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let a = Moc::from_cells(1, [12, 13]);
    /// let b = Moc::from_cells(1, [14, 15]);
    /// assert_eq!(a.union(&b), Moc::from_cells(0, [3]));
    /// assert_eq!(&a | &b, a.union(&b));
    /// ```
    pub fn union(&self, other: &Moc) -> Moc {
        let (a, b) = (&self.ranges, &other.ranges);
        let mut out = Vec::with_capacity(a.len() + b.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() || j < b.len() {
            // Take whichever list starts next, so `out` is built in order.
            let take_a = j >= b.len() || (i < a.len() && a[i].start <= b[j].start);
            let next = if take_a {
                i += 1;
                a[i - 1].clone()
            } else {
                j += 1;
                b[j - 1].clone()
            };
            push_coalesced(&mut out, next);
        }
        Moc { ranges: out }
    }

    /// The intersection of two coverages: everything in both.
    ///
    /// Also available as `&a & &b`.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let a = Moc::from_cells(1, [12, 13, 14]);
    /// let b = Moc::from_cells(1, [14, 15]);
    /// assert_eq!(a.intersection(&b), Moc::from_cells(1, [14]));
    /// assert_eq!(&a & &b, a.intersection(&b));
    /// ```
    pub fn intersection(&self, other: &Moc) -> Moc {
        let (a, b) = (&self.ranges, &other.ranges);
        // No capacity hint: an intersection is usually a small fraction of either input,
        // and reserving for the worst case measured slower than growing into it.
        let mut out = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            let start = a[i].start.max(b[j].start);
            let end = a[i].end.min(b[j].end);
            if start < end {
                push_coalesced(&mut out, start..end);
            }
            // Retire whichever range ends first; the other may still meet the next one.
            if a[i].end < b[j].end {
                i += 1;
            } else {
                j += 1;
            }
        }
        Moc { ranges: out }
    }

    /// The difference: everything in `self` that is not in `other`.
    ///
    /// Also available as `&a - &b`. This is the "what have we not covered yet" operation:
    /// a target region minus the sky already observed.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let target = Moc::from_cells(0, [3]);
    /// let observed = Moc::from_cells(1, [12, 13]);
    /// assert_eq!(&target - &observed, Moc::from_cells(1, [14, 15]));
    /// ```
    pub fn difference(&self, other: &Moc) -> Moc {
        let (a, b) = (&self.ranges, &other.ranges);
        // A difference keeps most of `a`, so this is close to the size actually reached;
        // growing into it from empty measured 44% slower.
        let mut out = Vec::with_capacity(a.len() + b.len());
        let mut j = 0;
        for r in a {
            let mut start = r.start;
            // Skip the subtrahends that end before this range begins. They are behind us
            // for every later range too, since `a` is sorted.
            while j < b.len() && b[j].end <= start {
                j += 1;
            }
            let mut k = j;
            while k < b.len() && b[k].start < r.end {
                if b[k].start > start {
                    push_coalesced(&mut out, start..b[k].start);
                }
                start = start.max(b[k].end);
                if start >= r.end {
                    break;
                }
                k += 1;
            }
            if start < r.end {
                push_coalesced(&mut out, start..r.end);
            }
        }
        Moc { ranges: out }
    }

    /// The symmetric difference: everything in exactly one of the two coverages.
    ///
    /// Also available as `&a ^ &b`. The union less the overlap, so it is what two
    /// exposures reached that the other one did not.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let a = Moc::from_cells(1, [12, 13]);
    /// let b = Moc::from_cells(1, [13, 14]);
    /// assert_eq!(a.symmetric_difference(&b), Moc::from_cells(1, [12, 14]));
    /// assert_eq!(&a ^ &b, a.symmetric_difference(&b));
    /// ```
    pub fn symmetric_difference(&self, other: &Moc) -> Moc {
        self.union(other).difference(&self.intersection(other))
    }

    /// The complement: the rest of the sphere.
    ///
    /// Also available as `!&a`.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let a = Moc::from_cells(0, [3]);
    /// assert_eq!(a.complement().area(), Moc::all_sky().area() - a.area());
    /// assert_eq!(&a | &a.complement(), Moc::all_sky());
    /// assert_eq!(!&!&a, a);
    /// ```
    pub fn complement(&self) -> Moc {
        let mut out = Vec::with_capacity(self.ranges.len() + 1);
        let mut start = 0u64;
        for r in &self.ranges {
            if r.start > start {
                out.push(start..r.start);
            }
            start = r.end;
        }
        if start < FULL_SKY {
            out.push(start..FULL_SKY);
        }
        Moc { ranges: out }
    }

    /// Ranges of NESTED cell indices at `depth` covering this coverage.
    ///
    /// Sorted, disjoint, non-adjacent and half-open, ready to slice a catalogue sorted by
    /// NESTED index at `depth`. A cell only partly inside the coverage is **included**, so
    /// at a depth shallower than the coverage was built at the result is a superset —
    /// matching the inclusive convention of
    /// [`cone_coverage`](crate::nested::Layer::cone_coverage). At a depth at least as deep
    /// as everything in the coverage it is exact.
    ///
    /// # Panics
    /// Panics if `depth > `[`MAX_DEPTH`].
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let moc = Moc::from_cells(1, [12, 13, 14, 15]);
    /// // Exact at or below the depth it was built at.
    /// assert_eq!(moc.ranges_at(1), vec![12..16]);
    /// assert_eq!(moc.ranges_at(2), vec![48..64]);
    /// // Rounded outward above it: the partly covered parent is included.
    /// assert_eq!(Moc::from_cells(1, [12]).ranges_at(0), vec![3..4]);
    /// ```
    pub fn ranges_at(&self, depth: u8) -> Vec<Range<u64>> {
        let shift = shift_at(depth);
        let round_up = (1u64 << shift) - 1;
        let mut out: Vec<Range<u64>> = Vec::with_capacity(self.ranges.len());
        for r in &self.ranges {
            // Round outward: down to the cell holding the start, up past the one holding
            // the last cell inside.
            push_coalesced(&mut out, (r.start >> shift)..((r.end + round_up) >> shift));
        }
        out
    }

    /// The largest NESTED cells that exactly tile this coverage, as `(depth, cell)` pairs
    /// in increasing order.
    ///
    /// This is the multi-order view: a coverage built from a cone at depth 12 comes back as
    /// a few big shallow cells in the middle and progressively smaller ones towards the
    /// boundary, never more cells than it takes to tile the region exactly.
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// // Four sibling cells are emitted as the single parent that they tile.
    /// let moc = Moc::from_cells(1, [12, 13, 14, 15]);
    /// assert_eq!(moc.cells().collect::<Vec<_>>(), [(0, 3)]);
    ///
    /// // Three of them cannot be, so they come back at their own depth.
    /// let moc = Moc::from_cells(1, [12, 13, 14]);
    /// assert_eq!(moc.cells().collect::<Vec<_>>(), [(1, 12), (1, 13), (1, 14)]);
    /// ```
    #[inline]
    pub fn cells(&self) -> Cells<'_> {
        Cells {
            ranges: &self.ranges,
            index: 0,
            next: 0,
        }
    }

    /// The cells of [`cells`](Self::cells) as [NUNIQ](crate::to_uniq) values, which is how
    /// a MOC is normally serialised.
    ///
    /// These come out in sky order, matching [`cells`](Self::cells) — not sorted by NUNIQ
    /// value, which would group them by depth. Sort them if your serialiser wants that
    /// order.
    ///
    /// ```
    /// use realpix::moc::Moc;
    ///
    /// let moc = Moc::from_cells(1, [12, 13, 14, 15]);
    /// assert_eq!(moc.uniq_cells().collect::<Vec<_>>(), [realpix::to_uniq(0, 3)]);
    /// // The pairing with `from_uniq_cells` round-trips any coverage.
    /// assert_eq!(Moc::from_uniq_cells(moc.uniq_cells()), moc);
    /// ```
    #[inline]
    pub fn uniq_cells(&self) -> impl Iterator<Item = u64> + '_ {
        self.cells()
            .map(|(depth, cell)| crate::to_uniq(depth, cell))
    }

    /// The coverage as ranges of depth-[`MAX_DEPTH`] cell indices — the representation it
    /// is stored in.
    ///
    /// Sorted, disjoint and non-adjacent. Pair with [`Moc::from_ranges`] at
    /// [`MAX_DEPTH`] to rebuild the coverage exactly.
    ///
    /// ```
    /// use realpix::{MAX_DEPTH, moc::Moc};
    ///
    /// let moc = Moc::from_cone(8, realpix::lonlat_to_vec(1.0, 0.5), 0.02);
    /// let ranges = moc.deep_ranges().to_vec();
    /// assert!(ranges.windows(2).all(|w| w[0].end < w[1].start));
    /// assert_eq!(Moc::from_ranges(MAX_DEPTH, ranges), moc);
    /// ```
    #[inline]
    pub fn deep_ranges(&self) -> &[Range<u64>] {
        &self.ranges
    }
}

/// Appends `range` to `out`, merging it into the last entry if they touch or overlap.
///
/// `range` must start at or after the last entry already in `out`.
#[inline]
fn push_coalesced(out: &mut Vec<Range<u64>>, range: Range<u64>) {
    match out.last_mut() {
        Some(last) if range.start <= last.end => last.end = last.end.max(range.end),
        _ => out.push(range),
    }
}

/// Iterator over the largest NESTED cells tiling a [`Moc`], returned by [`Moc::cells`].
#[derive(Debug, Clone)]
pub struct Cells<'a> {
    ranges: &'a [Range<u64>],
    index: usize,
    /// The next depth-`MAX_DEPTH` cell to emit from `ranges[index]`.
    next: u64,
}

impl Iterator for Cells<'_> {
    type Item = (u8, u64);

    fn next(&mut self) -> Option<(u8, u64)> {
        let range = loop {
            let range = self.ranges.get(self.index)?;
            if self.next < range.start {
                self.next = range.start;
            }
            if self.next < range.end {
                break range;
            }
            self.index += 1;
        };

        // The largest cell starting at `next` is bounded both by how far `next` is aligned
        // and by how much of the range is left: a cell at depth `MAX_DEPTH - k` spans
        // `4^k` and must start on a multiple of `4^k`.
        let alignment = if self.next == 0 {
            MAX_DEPTH as u32
        } else {
            (self.next.trailing_zeros() >> 1).min(MAX_DEPTH as u32)
        };
        let remaining = range.end - self.next;
        let mut k = alignment;
        while (1u64 << (k << 1)) > remaining {
            k -= 1;
        }

        let cell = self.next >> (k << 1);
        self.next += 1 << (k << 1);
        Some((MAX_DEPTH - k as u8, cell))
    }
}

impl core::iter::FusedIterator for Cells<'_> {}

impl BitOr for &Moc {
    type Output = Moc;
    /// [`Moc::union`].
    #[inline]
    fn bitor(self, rhs: &Moc) -> Moc {
        self.union(rhs)
    }
}

impl BitAnd for &Moc {
    type Output = Moc;
    /// [`Moc::intersection`].
    #[inline]
    fn bitand(self, rhs: &Moc) -> Moc {
        self.intersection(rhs)
    }
}

impl Sub for &Moc {
    type Output = Moc;
    /// [`Moc::difference`].
    #[inline]
    fn sub(self, rhs: &Moc) -> Moc {
        self.difference(rhs)
    }
}

impl BitXor for &Moc {
    type Output = Moc;
    /// [`Moc::symmetric_difference`].
    #[inline]
    fn bitxor(self, rhs: &Moc) -> Moc {
        self.symmetric_difference(rhs)
    }
}

impl Not for &Moc {
    type Output = Moc;
    /// [`Moc::complement`].
    #[inline]
    fn not(self) -> Moc {
        self.complement()
    }
}
