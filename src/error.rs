//! Error type.

use core::fmt;

/// Alias for the `Result` type returned by this crate.
pub type Result<T> = core::result::Result<T, Error>;

/// Errors returned by the checked entry points of this crate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// A depth outside `0..=`[`crate::MAX_DEPTH`] was given.
    InvalidDepth {
        /// The offending depth.
        depth: u8,
    },
    /// An `nside` that is not a power of two in `1..=2^29` was given.
    InvalidNside {
        /// The offending nside.
        nside: u32,
    },
    /// A cell index outside `0..12 * 4^depth` was given.
    InvalidCell {
        /// The offending cell index.
        cell: u64,
        /// The depth the index was interpreted at.
        depth: u8,
    },
    /// A coordinate was not finite, or was outside its valid range.
    InvalidCoordinate,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidDepth { depth } => {
                write!(
                    f,
                    "depth {depth} is out of range, expected 0..={}",
                    crate::MAX_DEPTH
                )
            }
            Error::InvalidNside { nside } => {
                write!(f, "nside {nside} is not a power of two in 1..=536870912")
            }
            Error::InvalidCell { cell, depth } => {
                write!(
                    f,
                    "cell {cell} is out of range at depth {depth}, expected 0..{}",
                    crate::n_hash(*depth)
                )
            }
            Error::InvalidCoordinate => write!(f, "coordinate is not finite or out of range"),
        }
    }
}

impl core::error::Error for Error {}
