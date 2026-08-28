//! The `(face, x, y)` intermediate representation and the index encodings built on it.
//!
//! Every HEALPix cell is a column `x` and a row `y` inside one of the 12 base cells. Both
//! numbering schemes are a pure re-encoding of that triple: NESTED interleaves the bits of
//! `x` and `y` after the face index, RING walks iso-latitude rings from the north pole.

use crate::base::Base;

/// Ring index of the north corner of each base cell, in units of `nside`.
pub(crate) const JRLL: [i64; 12] = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
/// Longitude index of the north corner of each base cell, in units of `π/4`.
pub(crate) const JPLL: [i64; 12] = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7];

/// A cell as a base-cell index plus its column and row inside that base cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Xyf {
    pub(crate) x: u32,
    pub(crate) y: u32,
    pub(crate) face: u8,
}

/// One of the eight neighbours of a cell, in the order returned by
/// [`nested::Layer::neighbours`](crate::nested::Layer::neighbours).
///
/// The order matches the reference implementation (and `healpy`'s `get_all_neighbours`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Direction {
    /// South-west.
    SW = 0,
    /// West.
    W = 1,
    /// North-west.
    NW = 2,
    /// North.
    N = 3,
    /// North-east.
    NE = 4,
    /// East.
    E = 5,
    /// South-east.
    SE = 6,
    /// South.
    S = 7,
}

impl Direction {
    /// All eight directions, in index order.
    pub const ALL: [Direction; 8] = [
        Direction::SW,
        Direction::W,
        Direction::NW,
        Direction::N,
        Direction::NE,
        Direction::E,
        Direction::SE,
        Direction::S,
    ];

    /// Index of this direction into a neighbour array.
    #[inline(always)]
    pub const fn index(self) -> usize {
        self as usize
    }
}

/// Spreads the low 32 bits of `x` into the even bit positions of the result.
#[inline(always)]
pub(crate) const fn spread_bits(x: u32) -> u64 {
    let mut v = x as u64;
    v = (v | (v << 16)) & 0x0000_FFFF_0000_FFFF;
    v = (v | (v << 8)) & 0x00FF_00FF_00FF_00FF;
    v = (v | (v << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    v = (v | (v << 2)) & 0x3333_3333_3333_3333;
    v = (v | (v << 1)) & 0x5555_5555_5555_5555;
    v
}

/// Gathers the even bits of `v` into the low 32 bits of the result.
#[inline(always)]
pub(crate) const fn compact_bits(v: u64) -> u32 {
    let mut v = v & 0x5555_5555_5555_5555;
    v = (v | (v >> 1)) & 0x3333_3333_3333_3333;
    v = (v | (v >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    v = (v | (v >> 4)) & 0x00FF_00FF_00FF_00FF;
    v = (v | (v >> 8)) & 0x0000_FFFF_0000_FFFF;
    v = (v | (v >> 16)) & 0x0000_0000_FFFF_FFFF;
    v as u32
}

// Neighbour lookup tables, from the reference implementation. Indexed by
// `[direction bucket][face]`, where the bucket encodes which face edge was crossed.
const XOFFSET: [i32; 8] = [-1, -1, 0, 1, 1, 1, 0, -1];
const YOFFSET: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];

#[rustfmt::skip]
const FACEARRAY: [[i8; 12]; 9] = [
    [ 8,  9, 10, 11, -1, -1, -1, -1, 10, 11,  8,  9], // S
    [ 5,  6,  7,  4,  8,  9, 10, 11,  9, 10, 11,  8], // SE
    [-1, -1, -1, -1,  5,  6,  7,  4, -1, -1, -1, -1], // E
    [ 4,  5,  6,  7, 11,  8,  9, 10, 11,  8,  9, 10], // SW
    [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11], // centre
    [ 1,  2,  3,  0,  0,  1,  2,  3,  5,  6,  7,  4], // NE
    [-1, -1, -1, -1,  7,  4,  5,  6, -1, -1, -1, -1], // W
    [ 3,  0,  1,  2,  3,  0,  1,  2,  4,  5,  6,  7], // NW
    [ 2,  3,  0,  1, -1, -1, -1, -1,  0,  1,  2,  3], // N
];

#[rustfmt::skip]
const SWAPARRAY: [[u8; 12]; 9] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3], // S
    [0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6], // SE
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], // E
    [0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5], // SW
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], // centre
    [5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0], // NE
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], // W
    [6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0], // NW
    [3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0], // N
];

impl Base {
    /// `(face, x, y)` -> NESTED index.
    #[inline(always)]
    pub(crate) const fn xyf2nested(&self, xyf: Xyf) -> u64 {
        ((xyf.face as u64) << self.twice_depth) | spread_bits(xyf.x) | (spread_bits(xyf.y) << 1)
    }

    /// NESTED index -> `(face, x, y)`.
    #[inline(always)]
    pub(crate) const fn nested2xyf(&self, cell: u64) -> Xyf {
        let face = (cell >> self.twice_depth) as u8;
        let in_face = cell & ((1u64 << self.twice_depth) - 1);
        Xyf {
            x: compact_bits(in_face),
            y: compact_bits(in_face >> 1),
            face,
        }
    }

    /// `(face, x, y)` -> RING index.
    #[inline]
    pub(crate) fn xyf2ring(&self, xyf: Xyf) -> u64 {
        let nside = self.nside as i64;
        let nl4 = nside << 2;
        let (ix, iy, f) = (xyf.x as i64, xyf.y as i64, xyf.face as usize);

        // Ring index counted from the north pole, 1-based.
        let jr = JRLL[f] * nside - ix - iy - 1;

        let (n_before, nr, kshift) = if jr < nside {
            // North polar cap.
            (2 * jr * (jr - 1), jr, 0)
        } else if jr > 3 * nside {
            // South polar cap.
            let nr = nl4 - jr;
            (self.n_hash as i64 - 2 * nr * (nr + 1), nr, 0)
        } else {
            // Equatorial belt: every other ring is offset by half a cell.
            (
                self.ncap as i64 + (jr - nside) * nl4,
                nside,
                (jr - nside) & 1,
            )
        };

        // The numerator is always even (the half-cell offsets cancel), so the arithmetic
        // shift below is exactly the reference implementation's `/ 2`.
        debug_assert_eq!((JPLL[f] * nr + ix - iy + 1 + kshift) & 1, 0);
        let mut jp = (JPLL[f] * nr + ix - iy + 1 + kshift) >> 1;
        debug_assert!(jp <= 4 * nr);
        if jp < 1 {
            jp += nl4;
        }
        (n_before + jp - 1) as u64
    }

    /// RING index -> `(face, x, y)`.
    #[inline]
    pub(crate) fn ring2xyf(&self, cell: u64) -> Xyf {
        let nside = self.nside as i64;
        let nl2 = nside << 1;
        let pix = cell as i64;

        let (iring, iphi, kshift, nr, face) = if cell < self.ncap {
            // North polar cap: invert the triangular numbering of the rings.
            let iring = (1 + (1 + 2 * pix).isqrt()) >> 1; // counted from the north pole
            let iphi = (pix + 1) - 2 * iring * (iring - 1);
            (iring, iphi, 0, iring, ((iphi - 1) / iring) as usize)
        } else if cell < self.n_hash - self.ncap {
            // Equatorial belt.
            let ip = pix - self.ncap as i64;
            let tmp = ip >> (self.depth + 2);
            let iring = tmp + nside;
            let iphi = ip - tmp * (nside << 2) + 1;
            let kshift = (iring + nside) & 1;
            let ire = iring - nside + 1;
            let irm = nl2 + 2 - ire;
            let ifm = (iphi - (ire >> 1) + nside - 1) >> self.depth;
            let ifp = (iphi - (irm >> 1) + nside - 1) >> self.depth;
            let face = if ifp == ifm {
                ifp | 4
            } else if ifp < ifm {
                ifp
            } else {
                ifm + 8
            };
            (iring, iphi, kshift, nside, face as usize)
        } else {
            // South polar cap: same triangular numbering, counted from the south pole,
            // then flipped back to a north-counted ring index.
            let ip = self.n_hash as i64 - pix;
            let nr = (1 + (2 * ip - 1).isqrt()) >> 1; // counted from the south pole
            let iphi = 4 * nr + 1 - (ip - 2 * nr * (nr - 1));
            (2 * nl2 - nr, iphi, 0, nr, 8 + ((iphi - 1) / nr) as usize)
        };

        let irt = iring - JRLL[face] * nside + 1;
        let mut ipt = 2 * iphi - JPLL[face] * nr - kshift - 1;
        if ipt >= nl2 {
            ipt -= 8 * nside;
        }

        Xyf {
            x: ((ipt - irt) >> 1) as u32,
            y: ((-(ipt + irt)) >> 1) as u32,
            face: face as u8,
        }
    }

    /// The eight neighbours of `xyf`, in [`Direction`] order.
    ///
    /// `None` marks the missing neighbour of the 24 cells that sit on a base-cell corner
    /// where only three base cells meet.
    #[inline]
    pub(crate) fn neighbours_xyf(&self, xyf: Xyf) -> [Option<Xyf>; 8] {
        let nside = self.nside as i32;
        let mut out = [None; 8];
        for i in 0..8 {
            let mut x = xyf.x as i32 + XOFFSET[i];
            let mut y = xyf.y as i32 + YOFFSET[i];
            // Bucket 4 is "same face"; each edge crossing shifts the bucket.
            let mut bucket = 4i32;
            if x < 0 {
                x += nside;
                bucket -= 1;
            } else if x >= nside {
                x -= nside;
                bucket += 1;
            }
            if y < 0 {
                y += nside;
                bucket -= 3;
            } else if y >= nside {
                y -= nside;
                bucket += 3;
            }

            let face = FACEARRAY[bucket as usize][xyf.face as usize];
            if face >= 0 {
                let bits = SWAPARRAY[bucket as usize][xyf.face as usize];
                if bits & 1 != 0 {
                    x = nside - x - 1;
                }
                if bits & 2 != 0 {
                    y = nside - y - 1;
                }
                if bits & 4 != 0 {
                    core::mem::swap(&mut x, &mut y);
                }
                out[i] = Some(Xyf {
                    x: x as u32,
                    y: y as u32,
                    face: face as u8,
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spread_and_compact_are_inverse() {
        for x in [0u32, 1, 2, 3, 0x5555_5555, 0xFFFF_FFFF, 0x1234_5678] {
            assert_eq!(compact_bits(spread_bits(x)), x, "round trip for {x:#x}");
        }
        // Spread must only ever touch even bit positions.
        assert_eq!(spread_bits(0xFFFF_FFFF) & 0xAAAA_AAAA_AAAA_AAAA, 0);
        assert_eq!(spread_bits(1), 1);
        assert_eq!(spread_bits(2), 4);
        assert_eq!(spread_bits(3), 5);
        assert_eq!(spread_bits(0xFFFF_FFFF), 0x5555_5555_5555_5555);
    }

    #[test]
    fn interleaving_is_a_bijection_over_a_face() {
        let base = Base::new(4);
        let nside = base.nside;
        for x in 0..nside {
            for y in 0..nside {
                let xyf = Xyf { x, y, face: 7 };
                let cell = base.xyf2nested(xyf);
                assert!(base.contains(cell));
                assert_eq!(base.nested2xyf(cell), xyf);
            }
        }
    }

    #[test]
    fn face_index_survives_the_round_trip() {
        let base = Base::new(3);
        for face in 0..12u8 {
            let xyf = Xyf { x: 5, y: 2, face };
            assert_eq!(base.nested2xyf(base.xyf2nested(xyf)), xyf);
            assert_eq!(base.ring2xyf(base.xyf2ring(xyf)), xyf);
        }
    }

    #[test]
    fn direction_indices_are_dense_and_ordered() {
        for (i, d) in Direction::ALL.iter().enumerate() {
            assert_eq!(d.index(), i);
        }
    }
}
