use crate::numbering_scheme::NumberingScheme;
use crate::pixel::Pixel;

/// A trait that defines the basic operations for a HEALPix (Hierarchical Equal Area isoLatitude Pixelization)
/// grid structure. It provides methods for determining resolution, pixel counts, and coordinate transformations
/// between angular positions and pixel indices based on a numbering scheme.
///
/// This trait is used to represent and manipulate data on a HEALPix grid, which is a popular choice for
/// astronomical data analysis and geospatial applications due to its equal area and hierarchical properties.
pub trait Healpix {
    /// Returns the resolution of the HEALPix grid in number of pixels per face.
    fn face_resolution(&self) -> u32;

    /// Returns the total number of pixels in the HEALPix grid per face.
    fn pixels_per_face(&self) -> u32;

    /// Returns the total number of pixels in the HEALPix grid.
    fn total_pixels(&self) -> u32;

    /// Converts angular coordinates to pixel indices based on the HEALPix grid numbering scheme.
    fn angle_to_pixel<N: NumberingScheme>(&self, theta: f64, phi: f64) -> Pixel<N> {
        N::angle_to_pixel(self.face_resolution(), theta, phi)
    }

    /// Converts pixel indices to angular coordinates based on the HEALPix grid numbering scheme.
    fn pixel_to_angle<N: NumberingScheme>(&self, pixel: Pixel<N>) -> crate::Result<(f64, f64)> {
        N::pixel_to_angle(self.face_resolution(), pixel)
    }
}
