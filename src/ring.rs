use crate::hpd::Hpd;
use crate::numbering_scheme::NumberingScheme;
use crate::pixel::Pixel;
use core::f64::consts::{PI, TAU};

/// Precomputed constants from the original HEALPix C code.
static JRLL: [i32; 12] = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
static JPLL: [i32; 12] = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ring;

impl Ring {
    #[inline]
    fn north_polar(h: &Hpd, nside: f64) -> (f64, f64) {
        let x = (h.x as f64 + 0.5) / nside;
        let y = (h.y as f64 + 0.5) / nside;
        let z = 1.0 - (x + y).powi(2) / 3.0;
        let phi = (PI / 2.0) * (x - y) / (x + y);
        (z, phi.rem_euclid(TAU))
    }

    #[inline]
    fn equatorial(h: &Hpd, nside: f64) -> (f64, f64) {
        let x = (h.x as f64 + 0.5) / nside;
        let y = (h.y as f64 + 0.5) / nside;
        let z = (2.0 / 3.0) * (2.0 - (x + y));
        let phi = (PI / 2.0) * (x - y) + (h.f as f64 - 4.0) * (PI / 2.0);
        (z, phi.rem_euclid(TAU))
    }

    #[inline]
    fn south_polar(h: &Hpd, nside: f64) -> (f64, f64) {
        let x = (h.x as f64 + 0.5) / nside;
        let y = (h.y as f64 + 0.5) / nside;
        let z = -1.0 + (x + y).powi(2) / 3.0;
        let phi = (PI / 2.0) * (x - y) / (x + y);
        (z, phi.rem_euclid(TAU))
    }

    /// Convert a discrete face-coordinate (`Hpd`) to a global ring pixel index.
    pub fn hpd2ring(nside: i64, h: Hpd) -> i64 {
        let nl4 = 4 * nside;
        let jr = JRLL[h.f as usize] as i64 * nside - h.x - h.y - 1;

        if jr < nside {
            // North polar cap
            let mut jp = (JPLL[h.f as usize] as i64 * jr + h.x - h.y + 1) / 2;
            if jp >= nl4 {
                jp -= nl4;
            } else if jp < 1 {
                jp += nl4;
            }
            2 * jr * (jr - 1) + jp - 1
        } else if jr > 3 * nside {
            // South polar cap
            let jr = nl4 - jr;
            let mut jp = (JPLL[h.f as usize] as i64 * jr + h.x - h.y + 1) / 2;
            if jp >= nl4 {
                jp -= nl4;
            } else if jp < 1 {
                jp += nl4;
            }
            12 * nside * nside - 2 * (jr + 1) * jr + jp - 1
        } else {
            // Equatorial region
            let mut jp =
                (JPLL[h.f as usize] as i64 * nside + h.x - h.y + 1 + ((jr - nside) & 1)) / 2;
            if jp >= nl4 {
                jp -= nl4;
            } else if jp < 1 {
                jp += nl4;
            }
            2 * nside * (nside - 1) + (jr - nside) * nl4 + jp - 1
        }
    }

    /// Convert a global ring pixel index to discrete face coordinates (`Hpd`).
    fn ring2hpd(nside: i64, pix: i64) -> Hpd {
        let ncap = 2 * nside * (nside - 1);
        let npix = 12 * nside * nside;

        if pix < ncap {
            // North polar cap
            let iring = (1 + i64::isqrt(1 + 2 * pix)) >> 1;
            let iphi = pix + 1 - 2 * iring * (iring - 1);
            let face = (iphi - 1) / iring;
            let irt = iring - (JRLL[face as usize] as i64 * nside) + 1;
            let mut ipt = 2 * iphi - JPLL[face as usize] as i64 * iring - 1;
            if ipt >= 2 * nside {
                ipt -= 8 * nside;
            }
            Hpd {
                x: (ipt - irt) >> 1,
                y: (-(ipt + irt)) >> 1,
                f: face as i32,
            }
        } else if pix < (npix - ncap) {
            // Equatorial region
            let ip = pix - ncap;
            let iring = (ip / (4 * nside)) + nside;
            let iphi = (ip % (4 * nside)) + 1;
            let kshift = (iring + nside) & 1;
            let ire = iring - nside + 1;
            let irm = 2 * nside + 2 - ire;
            let ifm = (iphi - irm / 2 + nside - 1) / nside;
            let ifp = (iphi - ire / 2 + nside - 1) / nside;
            let face = if ifp == ifm {
                ifp | 4
            } else if ifp < ifm {
                ifp
            } else {
                ifm + 8
            };
            let irt = iring - (JRLL[face as usize] as i64 * nside) + 1;
            let mut ipt = 2 * iphi - JPLL[face as usize] as i64 * nside - kshift - 1;
            if ipt >= 2 * nside {
                ipt -= 8 * nside;
            }
            Hpd {
                x: (ipt - irt) >> 1,
                y: (-(ipt + irt)) >> 1,
                f: face as i32,
            }
        } else {
            // South polar cap
            let ip = npix - pix;
            let iring = (1 + i64::isqrt(2 * ip - 1)) >> 1;
            let iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
            let face = 8 + (iphi - 1) / iring;
            let irt = 4 * nside - iring - (JRLL[face as usize] as i64 * nside) + 1;
            let mut ipt = 2 * iphi - JPLL[face as usize] as i64 * iring - 1;
            if ipt >= 2 * nside {
                ipt -= 8 * nside;
            }
            Hpd {
                x: (ipt - irt) >> 1,
                y: (-(ipt + irt)) >> 1,
                f: face as i32,
            }
        }
    }
}

impl NumberingScheme for Ring {
    fn angle_to_pixel<N: NumberingScheme>(face_resolution: u32, theta: f64, phi: f64) -> Pixel<N> {
        let z = theta.cos();
        let za = z.abs();
        let nside = face_resolution as f64;
        let phi = phi.rem_euclid(TAU);

        let hpd = if za <= 2.0 / 3.0 {
            // Equatorial region
            let tt = phi / (PI / 2.0);
            let temp1 = 0.5 + tt;
            let temp2 = z * 0.75;
            let jp = temp1 - temp2;
            let jm = temp1 + temp2;
            let ifp = jp.floor() as u32;
            let ifm = jm.floor() as u32;
            let f = if ifp == ifm {
                ifp | 4
            } else if ifp < ifm {
                ifp
            } else {
                ifm + 8
            };
            let ix = ((jm - ifm as f64) * nside)
                .floor()
                .max(0.0)
                .min(nside - 1.0) as i64;
            let iy = ((1.0 + ifp as f64 - jp) * nside)
                .floor()
                .max(0.0)
                .min(nside - 1.0) as i64;
            Hpd {
                x: ix,
                y: iy,
                f: f as i32,
            }
        } else {
            // Polar regions
            let mut tt = phi / TAU * 4.0;
            if tt >= 4.0 {
                tt = 3.999999;
            }
            let tmp = (3.0 * (1.0 - za)).sqrt();
            let mut jp = tt * tmp;
            let mut jm = (1.0 - tt) * tmp;
            if jp > 1.0 {
                jp = 1.0;
            }
            if jm > 1.0 {
                jm = 1.0;
            }
            let (jp, jm) = if z >= 0.0 {
                (1.0 - jm, 1.0 - jp)
            } else {
                (jp, jm)
            };
            let f = if z >= 0.0 {
                tt.floor() as i32
            } else {
                (tt.floor() as i32) + 8
            };
            Hpd {
                x: (jp * nside).floor().max(0.0).min(nside - 1.0) as i64,
                y: (jm * nside).floor().max(0.0).min(nside - 1.0) as i64,
                f,
            }
        };

        Pixel::from_u64(Self::hpd2ring(face_resolution as i64, hpd) as u64)
    }

    fn pixel_to_angle<N: NumberingScheme>(
        face_resolution: u32,
        pixel: Pixel<N>,
    ) -> crate::Result<(f64, f64)> {
        let nside = face_resolution as i64;
        let hpd = Self::ring2hpd(nside, pixel.as_u64() as i64);
        let (z, phi) = if hpd.f < 4 {
            Ring::north_polar(&hpd, nside as f64)
        } else if hpd.f < 8 {
            Ring::equatorial(&hpd, nside as f64)
        } else {
            Ring::south_polar(&hpd, nside as f64)
        };
        Ok((z.acos(), phi))
    }
}
