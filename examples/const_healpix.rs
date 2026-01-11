use realpix::{ConstHealpix, Healpix, Nested};

const HEALPIX: ConstHealpix<64> = ConstHealpix::new();

pub fn main() {
    let pixel = HEALPIX.angle_to_pixel::<Nested>(0.0, 0.0);
    let (theta, phi) = HEALPIX.pixel_to_angle(pixel).unwrap();
}
