use crate::Healpix;

/// A struct representing a dynamic HEALPix (Hierarchical Equal Area isoLatitude Pixelization) grid.
///
/// The `DynamicHealpix` struct is used to model a hierarchical, equal-area pixelization of a sphere,
/// as commonly used in astrophysics for storing and analyzing spherical data.
///
/// # Fields
///
/// * `face_resolution` - The resolution of an individual face. This determines the level of detail
///   for each face of the HEALPix grid.
/// * `pixels_per_face` - The total number of pixels contained on a single face of the HEALPix grid.
///   This value is typically a function of the `face_resolution`.
/// * `total_pixels` - The total number of pixels in the entire HEALPix grid, which is the product
///   of `pixels_per_face` and the number of faces in the HEALPix system (commonly 12 for the standard
///   HEALPix grid).
///
/// # Example
///
/// ```rust
/// use realpix::DynamicHealpix;
///
/// let healpix = DynamicHealpix::new(16);
/// ```
///
/// This structure allows for flexibility in defining HEALPix grids, enabling both standardized
/// systems and custom configurations for specific applications.
pub struct DynamicHealpix {
    face_resolution: u32,
}

impl DynamicHealpix {
    pub fn new(face_resolution: u32) -> crate::Result<Self> {
        if !face_resolution.is_power_of_two() {
            return Err(crate::Error::InvalidFaceResolution);
        }
        Ok(Self { face_resolution })
    }
}

impl Healpix for DynamicHealpix {
    fn face_resolution(&self) -> u32 {
        self.face_resolution
    }

    fn pixels_per_face(&self) -> u32 {
        self.face_resolution * self.face_resolution
    }

    fn total_pixels(&self) -> u32 {
        12 * (self.face_resolution * self.face_resolution)
    }
}
