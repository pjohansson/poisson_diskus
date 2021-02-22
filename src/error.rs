use std::fmt;

#[derive(Clone, Debug)]
/// Errors encountered when sampling coordinates.
pub enum Error<const D: usize> {
    /// Invalid input box size to sample coordinates in.
    ///
    /// All box size lengths must be positive, real values.
    InvalidBoxSize { value: f64, box_size: Vec<f64> },
    /// Invalid input minimum distance between coordinates.
    ///
    /// Minimum distance must be a positive, real value. This will also be yielded
    /// for infinite and not-a-number values.
    InvalidRmin(f64),
    /// Generated a coordinate that was outside the box.
    ///
    /// This should not happen. Please file an issue if you encounter this.
    GenCoordOutOfBounds([f64; D]),
    /// The active list is inconsistent.
    ///
    /// This should not happen. Please file an issue if you encounter this.
    InvalidActiveList,
}

impl<const D: usize> fmt::Display for Error<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidBoxSize { value, box_size } => write!(
                f,
                "invalid value '{}' in input box size '{:?}': must be a positive number",
                value, box_size
            ),
            Error::InvalidRmin(rmin) => write!(
                f,
                "invalid value for input `rmin` ({}): must be a positive number",
                rmin
            ),
            Error::GenCoordOutOfBounds(coord) => write!(
                f,
                "generated and used an out-of-box coordinate ({:?}): this should not happen, please file an issue",
                coord
            ),
            Error::InvalidActiveList => write!(
                f,
                "active list contains an index which does not lead to a sample: this should not happen, please file an issue"
            ),
        }
    }
}

impl<const D: usize> std::error::Error for Error<D> {}

#[cfg(test)]
impl<const D: usize> Error<D> {
    pub(crate) fn is_invalid_box_size(&self) -> bool {
        match self {
            Error::InvalidBoxSize { .. } => true,
            _ => false,
        }
    }

    pub(crate) fn is_invalid_rmin(&self) -> bool {
        match self {
            Error::InvalidRmin(_) => true,
            _ => false,
        }
    }
}
