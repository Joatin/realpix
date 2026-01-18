use crate::pixel::Pixel;

/// A trait representing a numbering scheme for mapping between spherical angles
/// (given as theta and phi) and pixel indices or coordinates on a face
/// of a subdivided spherical surface.
///
/// This trait provides two primary methods for converting:
/// 1. From spherical angles to a pixel index (`angle_to_pixel`).
/// 2. From a pixel index back to spherical angles (`pixel_to_angle`).
///
/// The trait is generic over the implementing type `N`, which should also implement
/// the `NumberingScheme` trait.
pub trait NumberingScheme: PartialEq {
    fn angle_to_pixel<N: NumberingScheme>(pixels_per_face: u32, theta: f64, phi: f64) -> Pixel<N>;
    fn pixel_to_angle<N: NumberingScheme>(
        pixels_per_face: u32,
        pixel: Pixel<N>,
    ) -> crate::Result<(f64, f64)>;
}
