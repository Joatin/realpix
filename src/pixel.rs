use crate::numbering_scheme::NumberingScheme;
use core::marker::PhantomData;

/// A generic structure representing a pixel with an associated numbering scheme.
///
/// The `Pixel` struct is a generic type that encapsulates a value representing
/// a healpix pixel
///
/// # Type Parameters
/// - `N`: The numbering scheme used to associate the pixel value with a specific
///   coordinate or logical representation. This must implement the `NumberingScheme` trait.
///
/// # Example
/// ```rust
/// use realpix::{Pixel, Nested};
///
/// let pixel = Pixel::<Nested>::from_u64(42);
/// ```
///
/// # Notes
/// This struct does not contain runtime information associated with the numbering scheme,
/// as its purpose is to leverage the type system to enforce correctness.
pub struct Pixel<N: NumberingScheme>(u64, PhantomData<N>);

impl<N: NumberingScheme> Pixel<N> {
    pub fn from_u64(v: u64) -> Self {
        Self(v, PhantomData)
    }

    pub fn as_u64(&self) -> u64 {
        self.0
    }
}
