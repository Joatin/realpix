use crate::healpix::Healpix;

/// A constant generic struct representing a HEALPix (Hierarchical Equal Area isoLatitude
/// Pixelation) structure with a set face resolution.
///
/// HEALPix is a spatial indexing method commonly used for representing spherical data, such as in
/// astrophysics or geographic information systems (GIS). The constant generic parameter
/// `FACE_RESOLUTION` determines the resolution of the HEALPix faces, which directly impacts
/// the granularity of the spatial division.
///
/// # Type Parameters
/// * `FACE_RESOLUTION` - A `u32` constant specifying the resolution of the HEALPix face.
///                       Higher values indicate finer levels of detail.
///
/// # Example
/// ```
/// use realpix::ConstHealpix;
///
/// let healpix = ConstHealpix::<16>::new();
/// ```
pub struct ConstHealpix<const FACE_RESOLUTION: u32> {
    _private: (),
}

impl<const FACE_RESOLUTION: u32> ConstHealpix<FACE_RESOLUTION> {
    /// Creates a new instance of the struct.
    ///
    /// This constant function ensures that `FACE_RESOLUTION` is a power of two
    /// before constructing the instance. If `FACE_RESOLUTION` is not a power of two,
    /// the function will panic at compile time with an appropriate error message.
    ///
    /// # Panics
    /// Panics if `FACE_RESOLUTION` is not a power of two.
    ///
    /// # Returns
    /// A new instance of the struct.
    ///
    /// # Example
    /// ```
    /// use realpix::ConstHealpix;
    ///
    /// const HEALPIX: ConstHealpix<16> = ConstHealpix::new();
    /// ```
    pub const fn new() -> Self {
        assert!(
            FACE_RESOLUTION.is_power_of_two(),
            "FACE_RESOLUTION must be a power of two"
        );
        Self { _private: () }
    }
}

impl<const FACE_RESOLUTION: u32> Healpix for ConstHealpix<FACE_RESOLUTION> {
    fn face_resolution(&self) -> u32 {
        FACE_RESOLUTION
    }

    fn pixels_per_face(&self) -> u32 {
        FACE_RESOLUTION * FACE_RESOLUTION
    }

    fn total_pixels(&self) -> u32 {
        12 * (FACE_RESOLUTION * FACE_RESOLUTION)
    }
}
