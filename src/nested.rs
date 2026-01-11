use crate::numbering_scheme::NumberingScheme;
use crate::pixel::Pixel;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Nested;

impl NumberingScheme for Nested {
    fn angle_to_pixel<N: NumberingScheme>(face_resolution: u32, theta: f64, phi: f64) -> Pixel<N> {
        let z = theta.cos();
        let za = z.abs();
        let nside = face_resolution as f64;
        let phi = phi.rem_euclid(core::f64::consts::TAU);

        let (face, ix, iy) = if za <= 2.0 / 3.0 {
            // Equatorial region: 4 faces
            let tt = phi / (core::f64::consts::PI / 2.0); // [0,4)
            let face = 4 + tt.floor() as u32;
            let f = tt - tt.floor();

            // Map to face coordinates
            let x = nside * f;
            let y = nside * (0.75 * (1.0 - z) - f);

            // Clamp to valid range
            let ix = x.floor().max(0.0).min(nside - 1.0) as u32;
            let iy = y.floor().max(0.0).min(nside - 1.0) as u32;

            (face, ix, iy)
        } else {
            // Polar caps: 0-3 = north, 8-11 = south
            let tmp = nside * (3.0 * (1.0 - za)).sqrt();
            let f = (phi / (core::f64::consts::PI / 2.0)).fract();
            let ix = (f * tmp).floor().max(0.0).min(nside - 1.0) as u32;
            let iy = (tmp - ix as f64 - 1.0).floor().max(0.0).min(nside - 1.0) as u32;

            let face = if z > 0.0 {
                // North polar cap
                (f * 4.0).floor() as u32
            } else {
                // South polar cap
                8 + (f * 4.0).floor() as u32
            };

            (face, ix, iy)
        };

        Pixel::from_u64(
            (face as u64) * ((face_resolution * face_resolution) as u64) + interleave(ix, iy),
        )
    }

    fn pixel_to_angle<N: NumberingScheme>(
        face_resolution: u32,
        pixel: Pixel<N>,
    ) -> crate::Result<(f64, f64)> {
        let total_pixels = (12 * face_resolution * face_resolution) as u64;
        let pixel_per_face = (face_resolution * face_resolution) as u64;
        let pixel = pixel.as_u64();
        if pixel >= total_pixels {
            return Err(crate::Error::InvalidPixel);
        }

        let face = (pixel / pixel_per_face) as u32;
        let ipf = (pixel % pixel_per_face) as u32;
        let (ix, iy) = deinterleave(ipf);

        let nside = face_resolution as f64;
        let x = (ix as f64 + 0.5) / nside;
        let y = (iy as f64 + 0.5) / nside;

        let (z, phi) = if face < 4 {
            // North polar cap
            let z = 1.0 - (x + y).powi(2) / 3.0;
            let phi = (core::f64::consts::PI / 2.0) * (x - y) / (x + y);
            (z, phi)
        } else if face < 8 {
            // Equatorial region
            let z = (2.0 / 3.0) * (2.0 - (x + y));
            let phi = (core::f64::consts::PI / 2.0) * (x - y)
                + (face as f64 - 4.0) * (core::f64::consts::PI / 2.0);
            (z, phi)
        } else {
            // South polar cap
            let z = -1.0 + (x + y).powi(2) / 3.0;
            let phi = (core::f64::consts::PI / 2.0) * (x - y) / (x + y);
            (z, phi)
        };

        Ok((z.acos(), phi.rem_euclid(core::f64::consts::TAU)))
    }
}

#[inline]
fn interleave(x: u32, y: u32) -> u64 {
    spread_bits(x) | (spread_bits(y) << 1)
}

#[inline]
fn deinterleave(v: u32) -> (u32, u32) {
    (compact_bits(v), compact_bits(v >> 1))
}

#[inline]
fn spread_bits(x: u32) -> u64 {
    let mut r = 0u64;
    for i in 0..16 {
        r |= ((x as u64 >> i) & 1) << (2 * i);
    }
    r
}

#[inline]
fn compact_bits(x: u32) -> u32 {
    let mut r = 0u32;
    for i in 0..16 {
        r |= ((x >> (2 * i)) & 1) << i;
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::{PI, TAU};

    type N = Nested;

    #[test]
    fn pixel_count_is_correct() {
        let nside = 32;
        let total = 12 * nside * nside;

        let max_pixel = Pixel::<N>::from_u64(total as u64 - 1).as_u64();

        assert_eq!(max_pixel, total as u64 - 1);
    }

    #[test]
    fn angle_to_pixel_is_in_range() {
        let nside = 64;
        let total = (12 * nside * nside) as u64;

        for i in 0..1000 {
            let theta = PI * (i as f64 / 1000.0);
            let phi = TAU * ((i * 37) as f64 / 1000.0);

            let pix = Nested::angle_to_pixel::<N>(nside, theta, phi).as_u64();
            assert!(pix < total, "pixel out of range: {}", pix);
        }
    }

    #[test]
    fn poles_map_to_polar_faces() {
        let nside = 32;

        let pix_north = Nested::angle_to_pixel::<N>(nside, 0.0, 0.0).as_u64();
        let face_north = pix_north / (nside * nside) as u64;

        assert!(face_north < 4, "north pole not in north cap");

        let pix_south = Nested::angle_to_pixel::<N>(nside, PI, 0.0).as_u64();
        let face_south = pix_south / (nside * nside) as u64;

        assert!(face_south >= 8, "south pole not in south cap");
    }

    #[test]
    fn equator_maps_to_equatorial_faces() {
        let nside = 32;

        for i in 0..4 {
            let phi = (i as f64 + 0.5) * (PI / 2.0);
            let pix = Nested::angle_to_pixel::<N>(nside, PI / 2.0, phi).as_u64();
            let face = pix / (nside * nside) as u64;

            assert!(
                (4..8).contains(&face),
                "equator pixel in wrong face: {}",
                face
            );
        }
    }

    #[test]
    fn phi_wraparound_consistency() {
        let nside = 64;
        let theta = PI / 3.0;

        let pix1 = Nested::angle_to_pixel::<N>(nside, theta, 0.01).as_u64();
        let pix2 = Nested::angle_to_pixel::<N>(nside, theta, TAU + 0.01).as_u64();

        assert_eq!(pix1, pix2);
    }

    #[test]
    fn locality_small_perturbation() {
        let nside = 128;

        let theta = 0.6 * PI;
        let phi = 1.0;

        let pix1 = Nested::angle_to_pixel::<N>(nside, theta, phi).as_u64();
        let pix2 = Nested::angle_to_pixel::<N>(nside, theta + 1e-5, phi + 1e-5).as_u64();

        assert!(
            (pix1 as i64 - pix2 as i64).abs() < 20,
            "nested locality violated"
        );
    }

    #[test]
    fn invalid_pixel_is_rejected() {
        let nside = 16;
        let total = (12 * nside * nside) as u64;

        let bad_pixel = Pixel::<N>::from_u64(total);
        assert!(Nested::pixel_to_angle::<N>(nside, bad_pixel).is_err());
    }
}
