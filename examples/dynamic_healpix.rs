use realpix::{DynamicHealpix, Healpix, Nested};

pub fn main() {
    let healpix = DynamicHealpix::new(64).unwrap();

    let pixel = healpix.angle_to_pixel::<Nested>(0.0, 0.0);
    let (theta, phi) = healpix.pixel_to_angle(pixel).unwrap();
}
