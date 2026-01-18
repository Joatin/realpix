use latlong::{Float, RaDec, TangentPosition};

/// Gnomonic projection: RA/Dec → tangent plane (x, y)
pub fn gnomonic_project<T: Float>(
    center_ra_dec: RaDec<T>,
    ra_dec: &RaDec<T>,
) -> Option<TangentPosition<T>> {
    let delta_ra = ra_dec.ra - center_ra_dec.ra;
    let sin_dec = ra_dec.dec.sin();
    let cos_dec = ra_dec.dec.cos();
    let sin_center = center_ra_dec.dec.sin();
    let cos_center = center_ra_dec.dec.cos();

    let denom = sin_dec * sin_center + cos_dec * cos_center * delta_ra.cos();
    if denom <= T::from(0.0) {
        return None; // behind tangent plane
    }

    let x = cos_dec * delta_ra.sin() / denom;
    let y = (cos_center * sin_dec - sin_center * cos_dec * delta_ra.cos()) / denom;

    Some(TangentPosition { x, y })
}
