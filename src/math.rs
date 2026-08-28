//! Floating point routines, backed by `std` or by `libm` in `no_std` builds.

#[cfg(all(not(feature = "std"), not(feature = "libm")))]
compile_error!(
    "realpix requires floating point support: enable either the `std` (default) or the `libm` feature"
);

pub(crate) use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

/// `1 / (π/2)`: longitude in radians times this gives longitude in units of 90 degrees.
pub(crate) const INV_HALF_PI: f64 = 2.0 / PI;
/// Boundary between the equatorial belt and the polar caps, in `z = sin(lat)`.
pub(crate) const TRANSITION_Z: f64 = 2.0 / 3.0;

#[cfg(feature = "std")]
mod imp {
    #[inline(always)]
    pub(crate) fn abs(x: f64) -> f64 {
        x.abs()
    }
    #[inline(always)]
    pub(crate) fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }
    /// Truncated remainder, i.e. C's `fmod`.
    #[inline(always)]
    pub(crate) fn fmod(x: f64, y: f64) -> f64 {
        x % y
    }
    #[inline(always)]
    pub(crate) fn sin(x: f64) -> f64 {
        x.sin()
    }
    #[inline(always)]
    pub(crate) fn cos(x: f64) -> f64 {
        x.cos()
    }
    #[inline(always)]
    pub(crate) fn sin_cos(x: f64) -> (f64, f64) {
        x.sin_cos()
    }
    #[inline(always)]
    pub(crate) fn asin(x: f64) -> f64 {
        x.asin()
    }
    #[inline(always)]
    pub(crate) fn atan2(y: f64, x: f64) -> f64 {
        y.atan2(x)
    }
}

#[cfg(not(feature = "std"))]
mod imp {
    #[inline(always)]
    pub(crate) fn abs(x: f64) -> f64 {
        libm::fabs(x)
    }
    #[inline(always)]
    pub(crate) fn sqrt(x: f64) -> f64 {
        libm::sqrt(x)
    }
    /// Truncated remainder, i.e. C's `fmod`.
    #[inline(always)]
    pub(crate) fn fmod(x: f64, y: f64) -> f64 {
        libm::fmod(x, y)
    }
    #[inline(always)]
    pub(crate) fn sin(x: f64) -> f64 {
        libm::sin(x)
    }
    #[inline(always)]
    pub(crate) fn cos(x: f64) -> f64 {
        libm::cos(x)
    }
    #[inline(always)]
    pub(crate) fn sin_cos(x: f64) -> (f64, f64) {
        libm::sincos(x)
    }
    #[inline(always)]
    pub(crate) fn asin(x: f64) -> f64 {
        libm::asin(x)
    }
    #[inline(always)]
    pub(crate) fn atan2(y: f64, x: f64) -> f64 {
        libm::atan2(y, x)
    }
}

pub(crate) use imp::*;

/// `x mod y`, always in `[0, y)` for `y > 0`. Matches the reference `fmodulo` exactly,
/// including its use of an exact `fmod` (so no rounding error is introduced).
#[inline(always)]
pub(crate) fn fmodulo(x: f64, y: f64) -> f64 {
    if x >= 0.0 {
        if x < y { x } else { fmod(x, y) }
    } else {
        let t = fmod(x, y) + y;
        if t == y { 0.0 } else { t }
    }
}

/// `atan2` that returns `0` instead of `NaN` at the origin.
#[inline(always)]
pub(crate) fn safe_atan2(y: f64, x: f64) -> f64 {
    if x == 0.0 && y == 0.0 {
        0.0
    } else {
        atan2(y, x)
    }
}

/// Dot product of two 3-vectors.
#[inline(always)]
pub(crate) fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmodulo_stays_in_range() {
        for (x, y, expected) in [
            (0.0, 4.0, 0.0),
            (1.5, 4.0, 1.5),
            (4.0, 4.0, 0.0),
            (5.5, 4.0, 1.5),
            (-0.5, 4.0, 3.5),
            (-4.0, 4.0, 0.0),
            (-8.5, 4.0, 3.5),
        ] {
            let got = fmodulo(x, y);
            assert!((got - expected).abs() < 1e-15, "fmodulo({x}, {y}) = {got}");
            assert!((0.0..y).contains(&got));
        }
    }

    #[test]
    fn safe_atan2_handles_the_origin() {
        assert_eq!(safe_atan2(0.0, 0.0), 0.0);
        assert!((safe_atan2(1.0, 0.0) - FRAC_PI_2).abs() < 1e-15);
    }
}
