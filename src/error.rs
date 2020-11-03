use std::fmt;

#[derive(Clone, Debug)]
/// Errors encountered when sampling coordinates.
pub enum Error<const DIM: usize> {
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
    GenCoordOutOfBounds([f64; DIM]),
    /// The active list is inconsistent.
    ///
    /// This should not happen. Please file an issue if you encounter this.
    InvalidActiveList,
    /// Tried to create a grid with mismatching dimensions between the box size
    /// and shape arrays.
    ///
    /// This should not happen. Please file an issue if you encounter this.
    UnmatchedDims,
}

impl<const DIM: usize> fmt::Display for Error<DIM> {
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
            Error::UnmatchedDims => write!(
                f,
                "created a shape of dimensions which do not match the given box size: this should not happen, please file an issue"
            ),
        }
    }
}

impl<const DIM: usize> std::error::Error for Error<DIM> {}

#[cfg(test)]
impl<const DIM: usize> Error<DIM> {
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
