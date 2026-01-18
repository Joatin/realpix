use crate::gnomonic_project::gnomonic_project;
use crate::numbering_scheme::NumberingScheme;
use crate::pixel::Pixel;
use latlong::{Declination, Float, RaDec, RightAscension, TangentPosition};

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

    fn ra_dec_to_pixel<N: NumberingScheme, T: Float>(&self, ra_dec: &RaDec<T>) -> Pixel<N> {
        let theta = core::f64::consts::FRAC_PI_2 - ra_dec.dec.radians().to_f64();
        let phi = ra_dec
            .ra
            .radians()
            .to_f64()
            .rem_euclid(core::f64::consts::TAU);
        N::angle_to_pixel(self.face_resolution(), theta, phi)
    }

    fn pixel_to_ra_dec<N: NumberingScheme, T: Float>(
        &self,
        pixel: Pixel<N>,
    ) -> crate::Result<RaDec<T>> {
        let (theta, phi) = N::pixel_to_angle(self.face_resolution(), pixel)?;
        let dec = core::f64::consts::FRAC_PI_2 - theta;
        let ra = phi.rem_euclid(core::f64::consts::TAU);
        Ok(RaDec {
            ra: RightAscension::from_radians(T::from(ra)),
            dec: Declination::from_radians(T::from(dec)),
        })
    }

    fn iter_pixels<N: NumberingScheme>(&self) -> impl Iterator<Item = Pixel<N>> + '_ {
        (0..self.total_pixels()).map(|index| Pixel::from_u64(index as u64))
    }

    fn project_ra_dec<N: NumberingScheme, T: Float>(
        &self,
        pixel: Pixel<N>,
        ra_dec: &RaDec<T>,
    ) -> Option<TangentPosition<T>> {
        let pixel_2 = self.ra_dec_to_pixel::<N, T>(ra_dec);
        if pixel != pixel_2 {
            None
        } else {
            let center_ra_dec = self.pixel_to_ra_dec::<N, T>(pixel).unwrap();
            gnomonic_project::<T>(center_ra_dec, ra_dec)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{ConstHealpix, Healpix, Nested};

    const HEALPIX: ConstHealpix<32> = ConstHealpix::new();

    #[test]
    fn iter_pixels_should_iterate_over_all_pixels() {
        let count = HEALPIX.iter_pixels::<Nested>().count();
        assert_eq!(count, HEALPIX.total_pixels() as usize);
    }
}
