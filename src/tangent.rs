//! Gnomonic (tangent plane) projection.

use crate::Vec3;
use crate::math::{cos, sin, sin_cos};

/// Projects `(lon, lat)` onto the plane tangent to the sphere at `(center_lon, center_lat)`.
///
/// The returned `(x, y)` are in units of the sphere's radius, with `x` towards increasing
/// longitude and `y` towards the north pole. Returns `None` when the position is on or
/// behind the horizon of the tangent point, where the projection is undefined.
///
/// ```
/// use realpix::gnomonic_project;
/// let (x, y) = gnomonic_project(0.0, 0.0, 0.01, 0.0).unwrap();
/// assert!((x - 0.01f64.tan()).abs() < 1e-12 && y.abs() < 1e-15);
/// ```
#[inline]
pub fn gnomonic_project(
    center_lon: f64,
    center_lat: f64,
    lon: f64,
    lat: f64,
) -> Option<(f64, f64)> {
    let d_lon = lon - center_lon;
    let (sin_d, cos_d) = sin_cos(d_lon);
    let (sin_lat, cos_lat) = sin_cos(lat);
    let (sin_c, cos_c) = sin_cos(center_lat);

    let denom = sin_lat * sin_c + cos_lat * cos_c * cos_d;
    if denom <= 0.0 {
        return None;
    }
    Some((
        cos_lat * sin_d / denom,
        (cos_c * sin_lat - sin_c * cos_lat * cos_d) / denom,
    ))
}

/// The angular distance, in radians, between two directions.
///
/// Uses the cross/dot form, which stays accurate for small angles where `acos` of a dot
/// product would lose most of its significant digits.
///
/// ```
/// let a = realpix::lonlat_to_vec(0.0, 0.0);
/// let b = realpix::lonlat_to_vec(1e-9, 0.0);
/// assert!((realpix::angular_distance(a, b) - 1e-9).abs() < 1e-24);
/// ```
#[inline]
pub fn angular_distance(a: Vec3, b: Vec3) -> f64 {
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    crate::math::atan2(
        crate::math::sqrt(crate::math::dot(&c, &c)),
        crate::math::dot(&a, &b),
    )
}

/// Converts `(lon, lat)` in radians to a unit vector.
#[inline]
pub fn lonlat_to_vec(lon: f64, lat: f64) -> Vec3 {
    let (sin_lon, cos_lon) = sin_cos(lon);
    let cos_lat = cos(lat);
    [cos_lat * cos_lon, cos_lat * sin_lon, sin(lat)]
}

/// Converts a unit vector to `(lon, lat)` in radians, with `lon` in `[0, 2π)`.
#[inline]
pub fn vec_to_lonlat(v: Vec3) -> (f64, f64) {
    let lon = crate::math::fmodulo(crate::math::safe_atan2(v[1], v[0]), crate::math::TAU);
    let lat = crate::math::atan2(v[2], crate::math::sqrt(v[0] * v[0] + v[1] * v[1]));
    (lon, lat)
}
