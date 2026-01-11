use thiserror::Error;

/// Enum `Error` represents possible error types that can occur in a specific context.
///
/// Variants:
/// - `InvalidFaceResolution`: This variant indicates that an invalid or unsupported face resolution value has been provided.
/// This could occur in scenarios where the input does not meet the required constraints or parameters.
#[derive(Error, Debug)]
pub enum Error {
    #[error("face resolution must be a power of two")]
    InvalidFaceResolution,

    #[error("The pixel is out of bounds")]
    InvalidPixel,
}
