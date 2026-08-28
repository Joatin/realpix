//! Depth (order) and `nside` helpers.

use crate::error::{Error, Result};

/// The deepest supported depth.
///
/// At depth 29 there are `12 * 4^29 = 3_458_764_513_820_540_928` cells, which is the largest
/// HEALPix layer whose NESTED indices fit in a `u64` (62 bits used).
pub const MAX_DEPTH: u8 = 29;

/// `nside = 2^depth`, the number of cells along one edge of a base cell.
///
/// # Panics
/// Panics if `depth > `[`MAX_DEPTH`].
#[inline(always)]
pub const fn nside(depth: u8) -> u32 {
    assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH");
    1u32 << depth
}

/// The total number of cells at `depth`, `12 * 4^depth`.
///
/// # Panics
/// Panics if `depth > `[`MAX_DEPTH`].
#[inline(always)]
pub const fn n_hash(depth: u8) -> u64 {
    assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH");
    12u64 << (depth << 1)
}

/// The depth matching `nside`, i.e. `log2(nside)`.
///
/// # Errors
/// Returns [`Error::InvalidNside`] if `nside` is not a power of two in `1..=2^29`.
#[inline]
pub const fn depth_from_nside(nside: u32) -> Result<u8> {
    if nside == 0 || !nside.is_power_of_two() || nside > (1u32 << MAX_DEPTH) {
        return Err(Error::InvalidNside { nside });
    }
    Ok(nside.trailing_zeros() as u8)
}

/// Checks that `depth` is supported.
#[inline]
pub(crate) const fn check_depth(depth: u8) -> Result<()> {
    if depth > MAX_DEPTH {
        Err(Error::InvalidDepth { depth })
    } else {
        Ok(())
    }
}
