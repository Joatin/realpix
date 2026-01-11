use crate::numbering_scheme::NumberingScheme;
use crate::pixel::Pixel;

pub struct Ring;

impl NumberingScheme for Ring {
    fn angle_to_pixel<N: NumberingScheme>(face_resolution: u32, theta: f64, phi: f64) -> Pixel<N> {
        todo!()
    }

    fn pixel_to_angle<N: NumberingScheme>(
        face_resolution: u32,
        pixel: Pixel<N>,
    ) -> crate::Result<(f64, f64)> {
        todo!()
    }
}
