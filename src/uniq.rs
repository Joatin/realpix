//! The `uniq` (NUNIQ) multi-resolution cell encoding.
//!
//! A NESTED cell index only makes sense together with its depth. NUNIQ folds both into a
//! single `u64` as `4 * 4^depth + cell`, which is the encoding used by MOCs. Values are
//! unique across depths, and ordering by NUNIQ groups cells by depth.

use crate::depth::MAX_DEPTH;

/// Encodes a NESTED `(depth, cell)` pair as a single NUNIQ value.
///
/// # Panics
/// Panics if `depth > `[`MAX_DEPTH`] or `cell` is out of range for `depth`.
#[inline]
pub const fn to_uniq(depth: u8, cell: u64) -> u64 {
    assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH");
    assert!(
        cell < (12u64 << (depth << 1)),
        "cell out of range for depth"
    );
    (4u64 << (depth << 1)) + cell
}

/// Decodes a NUNIQ value back into a NESTED `(depth, cell)` pair.
///
/// # Panics
/// Panics if `uniq` is smaller than 4 (the smallest valid NUNIQ value) or encodes a depth
/// beyond [`MAX_DEPTH`].
#[inline]
pub const fn from_uniq(uniq: u64) -> (u8, u64) {
    assert!(uniq >= 4, "not a valid uniq value");
    // The leading set bit is at position 2 * depth + 2.
    let depth = ((63 - uniq.leading_zeros() as u8) - 2) >> 1;
    assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH");
    (depth, uniq - (4u64 << (depth << 1)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_values() {
        assert_eq!(to_uniq(0, 0), 4);
        assert_eq!(to_uniq(0, 11), 15);
        assert_eq!(to_uniq(1, 0), 16);
        assert_eq!(to_uniq(1, 47), 63);
        assert_eq!(to_uniq(2, 0), 64);
        assert_eq!(from_uniq(4), (0, 0));
        assert_eq!(from_uniq(15), (0, 11));
        assert_eq!(from_uniq(16), (1, 0));
        assert_eq!(from_uniq(63), (1, 47));
    }

    #[test]
    fn round_trips_at_the_deepest_depth() {
        let cell = (12u64 << (MAX_DEPTH << 1)) - 1;
        assert_eq!(from_uniq(to_uniq(MAX_DEPTH, cell)), (MAX_DEPTH, cell));
    }

    #[test]
    #[should_panic(expected = "cell out of range")]
    fn rejects_an_out_of_range_cell() {
        to_uniq(1, 48);
    }
}
