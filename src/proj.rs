//! The HEALPix projection: sphere <-> `(face, x, y)`.
//!
//! Both directions follow the reference C++ implementation (`healpix_base`) exactly, so
//! indices are bit-identical to it and to `healpy`.

use crate::base::Base;
use crate::math::{FRAC_PI_2, FRAC_PI_4, INV_HALF_PI, TRANSITION_Z, abs, fmodulo, sqrt};
use crate::xyf::{JPLL, JRLL, Xyf};

/// A position on the sphere as used by the projection: `z = sin(lat)`, `phi = lon`, and
/// optionally `sth = cos(lat)` when it is known to more precision than `sqrt(1 - z^2)`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Loc {
    pub(crate) z: f64,
    pub(crate) phi: f64,
    pub(crate) sth: f64,
    pub(crate) have_sth: bool,
}

impl Base {
    /// Projects a position onto the `(face, x, y)` grid.
    ///
    /// This is the reference `loc2pix` split from the index encoding. It is total: any
    /// finite input yields `face < 12`, `x < nside` and `y < nside`.
    #[inline]
    pub(crate) fn loc2xyf(&self, loc: Loc) -> Xyf {
        let za = abs(loc.z);
        // Longitude in units of 90 degrees, in [0, 4).
        let tt = fmodulo(loc.phi * INV_HALF_PI, 4.0);

        if za <= TRANSITION_Z {
            // Equatorial belt: the projection is linear, so the cell falls out of the
            // intersection of an ascending and a descending diagonal line index.
            let temp1 = self.nside_f64 * (0.5 + tt);
            let temp2 = self.nside_f64 * (loc.z * 0.75);
            let jp = (temp1 - temp2) as i64; // ascending edge line index
            let jm = (temp1 + temp2) as i64; // descending edge line index
            let ifp = jp >> self.depth;
            let ifm = jm >> self.depth;
            let face = if ifp == ifm {
                ifp | 4
            } else if ifp < ifm {
                ifp
            } else {
                ifm + 8
            };
            debug_assert!((0..12).contains(&face));
            let mask = self.nside_minus_1 as i64;
            Xyf {
                x: (jm & mask) as u32,
                y: (mask - (jp & mask)) as u32,
                face: face as u8,
            }
        } else {
            // Polar caps: the cell grid is compressed towards the pole by sqrt(3(1-|z|)).
            let ntt = if tt < 3.0 { tt as i64 } else { 3 };
            let tp = tt - ntt as f64;
            let tmp = if za < 0.99 || !loc.have_sth {
                self.nside_f64 * sqrt(3.0 * (1.0 - za))
            } else {
                // Near the pole `1 - za` cancels catastrophically; use sin(lat) instead.
                self.nside_f64 * loc.sth / sqrt((1.0 + za) / 3.0)
            };
            let max = self.nside_minus_1 as i64;
            let jp = ((tp * tmp) as i64).min(max);
            let jm = (((1.0 - tp) * tmp) as i64).min(max);
            debug_assert!(jp >= 0 && jm >= 0);
            if loc.z >= 0.0 {
                Xyf {
                    x: (max - jm) as u32,
                    y: (max - jp) as u32,
                    face: ntt as u8,
                }
            } else {
                Xyf {
                    x: jp as u32,
                    y: jm as u32,
                    face: (ntt + 8) as u8,
                }
            }
        }
    }

    /// The `z = sin(lat)` of a point on the `(face, x, y)` grid, without any trigonometry.
    ///
    /// This is the cheap half of [`xyf2loc`](Self::xyf2loc): the cone search uses it to
    /// reject cells by latitude alone before paying for a longitude.
    #[inline]
    pub(crate) fn xyf2z(&self, face: u8, x: f64, y: f64) -> f64 {
        let nside = self.nside_f64;
        let jr = (JRLL[face as usize] as f64) * nside - x - y;
        if jr < nside {
            1.0 - jr * jr * self.fact2
        } else if jr > 3.0 * nside {
            let nr = 4.0 * nside - jr;
            nr * nr * self.fact2 - 1.0
        } else {
            (2.0 * nside - jr) * self.fact1
        }
    }

    /// Deprojects a *continuous* position on the `(face, x, y)` grid back to the sphere.
    ///
    /// `x` and `y` range over `[0, nside]`; a cell centre is `(ix + 0.5, iy + 0.5)` and its
    /// corners are the four integer points around it. One routine therefore serves both
    /// cell centres and cell boundaries.
    #[inline]
    pub(crate) fn xyf2loc(&self, face: u8, x: f64, y: f64) -> Loc {
        debug_assert!(face < 12);
        let nside = self.nside_f64;
        // Ring coordinate counted from the north pole, in cells.
        let jr = (JRLL[face as usize] as f64) * nside - x - y;

        let (nr, z, sth) = if jr < nside {
            // North polar cap.
            let tmp = jr * jr * self.fact2;
            (jr, 1.0 - tmp, sqrt(tmp * (2.0 - tmp)))
        } else if jr > 3.0 * nside {
            // South polar cap.
            let nr = 4.0 * nside - jr;
            let tmp = nr * nr * self.fact2;
            (nr, tmp - 1.0, sqrt(tmp * (2.0 - tmp)))
        } else {
            // Equatorial belt.
            let z = (2.0 * nside - jr) * self.fact1;
            (nside, z, sqrt((1.0 - z) * (1.0 + z)))
        };

        let phi = if nr <= 0.0 {
            // Exactly at a pole: longitude is degenerate, take the limit along the face.
            FRAC_PI_4 * JPLL[face as usize] as f64
        } else {
            let mut t = JPLL[face as usize] as f64 * nr + x - y;
            let eight_nr = 8.0 * nr;
            if t < 0.0 {
                t += eight_nr;
            } else if t >= eight_nr {
                t -= eight_nr;
            }
            if nr == nside {
                // Same value as `0.5 * FRAC_PI_2 * t / nr`, without the division.
                0.75 * FRAC_PI_2 * t * self.fact1
            } else {
                0.5 * FRAC_PI_2 * t / nr
            }
        };

        Loc {
            z,
            phi,
            sth,
            have_sth: true,
        }
    }
}

impl Loc {
    /// Builds a `Loc` from a longitude/latitude pair.
    ///
    /// `have_sth` is left false so that the projection matches the reference `ang2pix`
    /// bit for bit.
    #[inline(always)]
    pub(crate) fn from_lonlat(lon: f64, lat: f64) -> Self {
        Loc {
            z: crate::math::sin(lat),
            phi: lon,
            sth: 0.0,
            have_sth: false,
        }
    }

    /// Builds a `Loc` from a colatitude/longitude pair, mirroring the reference `ang2pix`.
    #[inline(always)]
    pub(crate) fn from_theta_phi(theta: f64, phi: f64) -> Self {
        Loc {
            z: crate::math::cos(theta),
            phi,
            sth: 0.0,
            have_sth: false,
        }
    }

    /// Builds a `Loc` from a vector, mirroring the reference `vec2pix`.
    ///
    /// The vector is normalised, exactly as the reference does, so callers holding a
    /// direction of arbitrary length get the same answer as with a unit vector.
    #[inline(always)]
    pub(crate) fn from_vec(v: &[f64; 3]) -> Self {
        let inv_len = 1.0 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        let phi = crate::math::safe_atan2(v[1], v[0]);
        let z = v[2] * inv_len;
        if abs(z) > 0.99 {
            // Near a pole, sin(theta) from the vector itself is far better conditioned
            // than anything derived from z.
            Loc {
                z,
                phi,
                sth: sqrt(v[0] * v[0] + v[1] * v[1]) * inv_len,
                have_sth: true,
            }
        } else {
            Loc {
                z,
                phi,
                sth: 0.0,
                have_sth: false,
            }
        }
    }

    /// Longitude and latitude, in radians, with `lon` in `[0, 2π)`.
    #[inline(always)]
    pub(crate) fn to_lonlat(self) -> (f64, f64) {
        let lat = if self.have_sth && abs(self.z) > 0.99 {
            crate::math::atan2(self.z, self.sth)
        } else {
            crate::math::asin(self.z)
        };
        (fmodulo(self.phi, crate::math::TAU), lat)
    }

    /// Unit vector, without ever going through an inverse trigonometric function.
    #[inline(always)]
    pub(crate) fn to_vec(self) -> [f64; 3] {
        debug_assert!(self.have_sth);
        let (sin_phi, cos_phi) = crate::math::sin_cos(self.phi);
        [self.sth * cos_phi, self.sth * sin_phi, self.z]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::{FRAC_PI_2, PI, TAU};

    #[test]
    fn projection_stays_inside_the_grid() {
        // Whatever it is handed, `loc2xyf` must produce indices that are in range: the
        // table lookups downstream depend on it.
        let base = Base::new(6);
        let nside = base.nside;
        for z in [
            -1.0,
            -0.99,
            -2.0 / 3.0,
            -0.1,
            0.0,
            0.1,
            2.0 / 3.0,
            0.99,
            1.0,
        ] {
            for phi in [-100.0, -1e-9, 0.0, 0.1, FRAC_PI_2, PI, TAU, 1e9] {
                for have_sth in [false, true] {
                    let sth = sqrt((1.0 - z) * (1.0 + z));
                    let xyf = base.loc2xyf(Loc {
                        z,
                        phi,
                        sth,
                        have_sth,
                    });
                    assert!(xyf.face < 12, "face {} for z {z}, phi {phi}", xyf.face);
                    assert!(
                        xyf.x < nside && xyf.y < nside,
                        "({}, {}) for z {z}, phi {phi}",
                        xyf.x,
                        xyf.y
                    );
                }
            }
        }
    }

    #[test]
    fn deprojection_returns_a_unit_vector() {
        let base = Base::new(5);
        for face in 0..12u8 {
            for x in [0.0, 0.5, 3.0, 31.5, 32.0] {
                for y in [0.0, 0.5, 3.0, 31.5, 32.0] {
                    let v = base.xyf2loc(face, x, y).to_vec();
                    let norm = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                    assert!(
                        (norm - 1.0).abs() < 1e-12,
                        "face {face}, ({x}, {y}) -> {v:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_poles_deproject_exactly() {
        let base = Base::new(4);
        // The north corner of every northern base cell is the pole itself.
        for face in 0..4u8 {
            let loc = base.xyf2loc(face, 16.0, 16.0);
            assert_eq!(loc.z, 1.0);
            assert_eq!(loc.sth, 0.0);
        }
        for face in 8..12u8 {
            let loc = base.xyf2loc(face, 0.0, 0.0);
            assert_eq!(loc.z, -1.0);
            assert_eq!(loc.sth, 0.0);
        }
    }
}
