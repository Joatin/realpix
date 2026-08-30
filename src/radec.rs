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
            ///
            /// Right ascension maps to longitude and declination to latitude, so this is
            /// [`hash`](Self::hash) with the units unwrapped for you.
            ///
            /// ```
            /// use latlong::{Declination, RaDec, RightAscension};
            ///
            /// let layer = realpix::nested::get(10);
            /// let betelgeuse = RaDec {
            ///     ra: RightAscension::from_radians(1.549_729_f64),
            ///     dec: Declination::from_radians(0.129_277_f64),
            /// };
            /// assert_eq!(layer.hash_ra_dec(&betelgeuse), layer.hash(1.549_729, 0.129_277));
            /// ```
            #[inline]
            pub fn hash_ra_dec<T: Float>(&self, ra_dec: &RaDec<T>) -> u64 {
                let (lon, lat) = to_lonlat(ra_dec);
                self.hash(lon, lat)
            }

            /// The centre of `cell` as a right ascension / declination pair.
            ///
            /// # Panics
            /// Panics if `cell` is out of range for this depth.
            ///
            /// ```
            /// use latlong::RaDec;
            ///
            /// let layer = realpix::nested::get(10);
            /// let cell = layer.hash(1.0, 0.5);
            /// let center: RaDec<f64> = layer.center_ra_dec(cell);
            /// // The centre of a cell is inside that cell.
            /// assert_eq!(layer.hash_ra_dec(&center), cell);
            /// ```
            #[inline]
            pub fn center_ra_dec<T: Float>(&self, cell: u64) -> RaDec<T> {
                let (lon, lat) = self.center(cell);
                from_lonlat(lon, lat)
            }

            /// `cone_coverage` centred on a right ascension / declination, with `radius`
            /// in radians.
            ///
            /// ```
            /// use latlong::{Declination, RaDec, RightAscension};
            ///
            /// let layer = realpix::nested::get(8);
            /// let field = RaDec {
            ///     ra: RightAscension::from_radians(1.55_f64),
            ///     dec: Declination::from_radians(0.13_f64),
            /// };
            /// let mut cells = 0u64;
            /// layer.cone_coverage_ra_dec(&field, 5f64.to_radians(), |r| cells += r.end - r.start);
            /// assert!(cells > 0);
            /// ```
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
    ///
    /// ```
    /// use latlong::{Declination, RaDec, RightAscension};
    ///
    /// let layer = realpix::nested::get(10);
    /// let star = RaDec {
    ///     ra: RightAscension::from_radians(1.549_729_f64),
    ///     dec: Declination::from_radians(0.129_277_f64),
    /// };
    /// let cell = layer.hash_ra_dec(&star);
    /// // Inside its own cell, the star projects close to the tangent point.
    /// let p = layer.project_ra_dec::<f64>(cell, &star).unwrap();
    /// assert!(p.x.abs() < 0.01 && p.y.abs() < 0.01);
    /// // A different cell rejects it.
    /// let other = layer.neighbour(cell, realpix::Direction::N).unwrap();
    /// assert!(layer.project_ra_dec::<f64>(other, &star).is_none());
    /// ```
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
}

/// The colatitude `theta` matching a declination, in radians.
///
/// Colatitude is `0` at the north pole and `π` at the south, the convention
/// [`hash_theta_phi`](crate::nested::Layer::hash_theta_phi) takes.
///
/// ```
/// use latlong::Declination;
///
/// let equator = Declination::from_radians(0.0_f64);
/// assert!((realpix::radec::dec_to_theta(equator) - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
/// ```
#[inline]
pub fn dec_to_theta<T: Float>(dec: Declination<T>) -> f64 {
    FRAC_PI_2 - dec.radians().to_f64()
}
