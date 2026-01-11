#![no_std]

#[cfg(feature = "std")]
extern crate std;
#[cfg(feature = "std")]
extern crate std as core;

mod const_healpix;
mod dynamic_healpix;
mod error;
mod healpix;
mod hpd;
mod nested;
mod numbering_scheme;
mod pixel;
mod result;
mod ring;

pub use self::const_healpix::ConstHealpix;
pub use self::dynamic_healpix::DynamicHealpix;
pub use self::error::Error;
pub use self::healpix::Healpix;
pub use self::nested::Nested;
pub use self::numbering_scheme::NumberingScheme;
pub use self::pixel::Pixel;
pub use self::result::Result;
pub use self::ring::Ring;
