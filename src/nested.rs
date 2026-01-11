use crate::numbering_scheme::NumberingScheme;
use crate::pixel::Pixel;

pub struct Nested;

impl NumberingScheme for Nested {
    fn angle_to_pixel<N: NumberingScheme>(face_resolution: u32, theta: f64, phi: f64) -> Pixel<N> {
        let z = theta.cos();
        let za = z.abs();
        let nside = face_resolution as f64;
        let phi = phi.rem_euclid(std::f64::consts::TAU);

        let (face, ix, iy) = if za <= 2.0 / 3.0 {
            // Equatorial region: 4 faces
            let tt = phi / (std::f64::consts::PI / 2.0); // [0,4)
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
            let f = (phi / (std::f64::consts::PI / 2.0)).fract();
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
            let phi = (std::f64::consts::PI / 2.0) * (x - y) / (x + y);
            (z, phi)
        } else if face < 8 {
            // Equatorial region
            let z = (2.0 / 3.0) * (2.0 - (x + y));
            let phi = (std::f64::consts::PI / 2.0) * (x - y)
                + (face as f64 - 4.0) * (std::f64::consts::PI / 2.0);
            (z, phi)
        } else {
            // South polar cap
            let z = -1.0 + (x + y).powi(2) / 3.0;
            let phi = (std::f64::consts::PI / 2.0) * (x - y) / (x + y);
            (z, phi)
        };

        Ok((z.acos(), phi.rem_euclid(std::f64::consts::TAU)))
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
fn spread_bits(mut x: u32) -> u64 {
    let mut r = 0u64;
    for i in 0..16 {
        r |= ((x as u64 >> i) & 1) << (2 * i);
    }
    r
}

#[inline]
fn compact_bits(mut x: u32) -> u32 {
    let mut r = 0u32;
    for i in 0..16 {
        r |= ((x >> (2 * i)) & 1) << i;
    }
    r
}
