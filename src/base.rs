//! Depth-derived constants shared by both numbering schemes.

use crate::depth::MAX_DEPTH;

/// Constants derived from a depth, computed once and reused by every conversion.
///
/// Every field is derivable from `depth` by shifts alone (plus two divisions for the
/// projection factors), so `Base::new` is a `const fn` and folds away entirely when the
/// depth is a compile-time constant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Base {
    pub(crate) depth: u8,
    pub(crate) twice_depth: u8,
    pub(crate) nside: u32,
    pub(crate) nside_minus_1: u32,
    pub(crate) n_hash: u64,
    /// Number of cells in the north polar cap, `2 * nside * (nside - 1)`.
    pub(crate) ncap: u64,
    pub(crate) nside_f64: f64,
    /// `2 / (3 * nside)`, i.e. `2 * nside * fact2`.
    pub(crate) fact1: f64,
    /// `1 / (3 * nside^2)`, i.e. `4 / n_hash`.
    pub(crate) fact2: f64,
}

impl Base {
    #[inline]
    pub(crate) const fn new(depth: u8) -> Self {
        assert!(depth <= MAX_DEPTH, "depth must be <= MAX_DEPTH (29)");
        let nside = 1u32 << depth;
        let nside64 = nside as u64;
        let n_hash = 12u64 << (depth << 1);
        let fact2 = 4.0 / (n_hash as f64);
        Self {
            depth,
            twice_depth: depth << 1,
            nside,
            nside_minus_1: nside - 1,
            n_hash,
            ncap: (nside64 * (nside64 - 1)) << 1,
            nside_f64: nside as f64,
            fact1: ((nside64 << 1) as f64) * fact2,
            fact2,
        }
    }

    /// Cell count of the layer one level up (`self.depth - 1`), used by hierarchy helpers.
    #[inline(always)]
    pub(crate) const fn contains(&self, cell: u64) -> bool {
        cell < self.n_hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_match_their_definitions() {
        for depth in 0..=MAX_DEPTH {
            let b = Base::new(depth);
            let nside = b.nside as f64;
            assert_eq!(b.nside as u64, 1u64 << depth);
            assert_eq!(b.n_hash, 12 * (b.nside as u64) * (b.nside as u64));
            assert_eq!(b.ncap, 2 * (b.nside as u64) * (b.nside as u64 - 1));
            assert_eq!(b.nside_minus_1, b.nside - 1);
            assert_eq!(b.twice_depth, 2 * depth);
            assert!((b.fact2 - 1.0 / (3.0 * nside * nside)).abs() < 1e-18);
            assert!((b.fact1 - 2.0 / (3.0 * nside)).abs() < 1e-18);
        }
    }

    #[test]
    fn deepest_layer_fits_in_u64() {
        let b = Base::new(MAX_DEPTH);
        assert_eq!(b.n_hash, 3_458_764_513_820_540_928);
        assert!(b.n_hash < u64::MAX / 4);
    }
}
