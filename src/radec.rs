//! Integration with the [`latlong`] crate's astronomical coordinate types.
//!
//! Enabled by the `latlong` feature. Right ascension maps to HEALPix longitude and
//! declination to latitude, both in radians.

use latlong::{Declination, Float, RaDec, RightAscension, TangentPosition};

use crate::math::{FRAC_PI_2, TAU, fmodulo};
use crate::{nested, ring};

#[inline]
fn to_lonlat<T: Float>(ra_dec: &RaDec<T>) -> (f64, f64) {
    (
        fmodulo(ra_dec.ra.radians().to_f64(), TAU),
        ra_dec.dec.radians().to_f64(),
    )
}

#[inline]
fn from_lonlat<T: Float>(lon: f64, lat: f64) -> RaDec<T> {
    RaDec {
        ra: RightAscension::from_radians(T::from(lon)),
        dec: Declination::from_radians(T::from(lat)),
    }
}

macro_rules! ra_dec_api {
    ($layer:ty) => {
        impl $layer {
            /// The cell containing the given right ascension / declination.
            #[inline]
            pub fn hash_ra_dec<T: Float>(&self, ra_dec: &RaDec<T>) -> u64 {
                let (lon, lat) = to_lonlat(ra_dec);
                self.hash(lon, lat)
            }

            /// The centre of `cell` as a right ascension / declination pair.
            ///
            /// # Panics
            /// Panics if `cell` is out of range for this depth.
            #[inline]
            pub fn center_ra_dec<T: Float>(&self, cell: u64) -> RaDec<T> {
                let (lon, lat) = self.center(cell);
                from_lonlat(lon, lat)
            }
        }
    };
}

ra_dec_api!(nested::Layer);
ra_dec_api!(ring::Layer);

impl nested::Layer {
    /// Projects `ra_dec` onto the tangent plane touching the sphere at the centre of `cell`.
    ///
    /// Returns `None` when the position does not fall inside `cell`, or when it lies on or
    /// behind the horizon of the tangent point.
    ///
    /// # Panics
    /// Panics if `cell` is out of range for this depth.
    #[inline]
    pub fn project_ra_dec<T: Float>(
        &self,
        cell: u64,
        ra_dec: &RaDec<T>,
    ) -> Option<TangentPosition<T>> {
        let (lon, lat) = to_lonlat(ra_dec);
        if self.hash(lon, lat) != cell {
            return None;
        }
        let (center_lon, center_lat) = self.center(cell);
        crate::tangent::gnomonic_project(center_lon, center_lat, lon, lat).map(|(x, y)| {
            TangentPosition {
                x: T::from(x),
                y: T::from(y),
            }
        })
    }

    /// [`cone_coverage`](nested::Layer::cone_coverage) centred on a right ascension /
    /// declination, with `radius` in radians.
    pub fn cone_coverage_ra_dec<T: Float, F: FnMut(core::ops::Range<u64>)>(
        &self,
        ra_dec: &RaDec<T>,
        radius: f64,
        sink: F,
    ) {
        let (lon, lat) = to_lonlat(ra_dec);
        self.cone_coverage_lonlat(lon, lat, radius, sink)
    }
}

/// The colatitude `theta` matching a declination, in radians.
#[inline]
pub fn dec_to_theta<T: Float>(dec: Declination<T>) -> f64 {
    FRAC_PI_2 - dec.radians().to_f64()
}
