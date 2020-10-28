use std::fmt;

#[derive(Clone, Debug)]
/// Errors encountered when sampling coordinates.
///
/// None of these should be due to invalid input from the user but due to inconsistencies
/// within the library itself. Please open an issue if you ever encounter these.
pub enum Error {
    GenCoordOutOfBounds(Vec<f64>),
    InvalidActiveList,
    UnmatchedDims,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::GenCoordOutOfBounds(coord) => write!(
                f,
                "generated and used an out-of-box coordinate ({:?})",
                coord
            ),
            Error::InvalidActiveList => write!(
                f,
                "active list contains an index which does not lead to a sample"
            ),
            Error::UnmatchedDims => write!(
                f,
                "created a shape of dimensions which do not match the given box size"
            ),
        }
    }
}

impl std::error::Error for Error {}
