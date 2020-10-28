use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    UnmatchedShapeAndSize { shape: Vec<usize>, size: Vec<f64> },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::UnmatchedShapeAndSize { size, shape } => write!(
                f,
                "could not create grid: shape ({:?}) and size ({:?}) vectors do not have equal dimension",
                shape, size
            ),
        }
    }
}

impl std::error::Error for Error {}

#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    /// Conversion multipliers to convert to and from bin positions and the 1d index in `data`.
    dim_multipliers: Vec<usize>,
    /// Number of bins around each bin which are within range of the minimum distance.
    pub(crate) num_adjacent: usize,
    /// Shape of grid.
    pub(crate) shape: Vec<usize>,
    /// Spacing of bins in grid.
    pub(crate) spacing: Vec<f64>,
    /// Data of grid as a one-dimensional vector.
    pub data: Vec<Option<usize>>,
}

impl Grid {
    pub fn new(shape: &[usize], size: &[f64]) -> Result<Self, Error> {
        if shape.len() != size.len() {
            return Err(Error::UnmatchedShapeAndSize {
                shape: shape.to_vec(),
                size: size.to_vec(),
            });
        }

        let spacing = size
            .iter()
            .zip(shape.iter())
            .map(|(s, &n)| s / (n as f64))
            .collect();

        Ok(Grid {
            dim_multipliers: create_dimension_multipliers(shape),
            num_adjacent: (shape.len() as f64).sqrt().ceil() as usize,
            shape: shape.to_vec(),
            spacing,
            data: vec![None; shape.iter().product()],
        })
    }

    pub fn get_index(&self, position: &[usize]) -> Option<usize> {
        if position.iter().zip(self.shape.iter()).any(|(p, s)| p >= s) {
            return None;
        }

        Some(
            position
                .iter()
                .zip(self.dim_multipliers.iter())
                .map(|(p, m)| p * m)
                .sum(),
        )
    }

    pub fn get_position(&self, index: usize) -> Option<Vec<usize>> {
        if index < self.data.len() {
            Some(
                self.dim_multipliers
                    .iter()
                    .scan(index, |i, d| {
                        let n = *i / d;
                        *i = *i % d;

                        Some(n)
                    })
                    .collect(),
            )
        } else {
            None
        }
    }

    pub fn get_position_from_coord(&self, coord: &[f64]) -> Option<Vec<usize>> {
        coord
            .iter()
            .zip(self.spacing.iter())
            .map(|(c, dx)| (c / dx).floor() as isize)
            .zip(self.shape.iter().map(|n| *n as isize))
            .map(|(i, n)| {
                if i >= 0 && i < n {
                    Some(i as usize)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_index_from_coord(&self, coord: &[f64]) -> Option<usize> {
        self.get_position_from_coord(coord)
            .and_then(|position| self.get_index(&position))
    }
}

fn create_dimension_multipliers(shape: &[usize]) -> Vec<usize> {
    let mut dim_multipliers = vec![1; shape.len()];

    for (i, s) in shape.iter().rev().enumerate() {
        for d in dim_multipliers.iter_mut().rev().skip(i + 1) {
            *d *= s;
        }
    }

    dim_multipliers
}

#[cfg(test)]
mod tests {
    use super::*;

    impl Grid {
        fn with_shape(shape: &[usize]) -> Self {
            let size: Vec<f64> = shape.iter().map(|v| *v as f64).collect();

            Grid::new(shape, &size).unwrap()
        }
    }

    #[test]
    fn grid_with_shape_is_consistent_with_new() {
        let shape = [1, 2, 3, 4];

        let grid_with_shape = Grid::with_shape(&shape);
        let grid_new = Grid::new(&shape, &[1.0, 2.0, 3.0, 4.0]).unwrap();

        assert_eq!(grid_with_shape, grid_new);
    }

    #[test]
    fn grid_new_requires_matching_shape_and_size_vector_lengths() {
        // Mathing lengths
        assert!(Grid::new(&[], &[]).is_ok());
        assert!(Grid::new(&[3], &[5.0]).is_ok());
        assert!(Grid::new(&[3, 5], &[5.0, 7.0]).is_ok());

        // Non-matching lengths
        assert!(Grid::new(&[3], &[5.0, 7.0]).is_err());
        assert!(Grid::new(&[3, 5], &[5.0]).is_err());
    }

    #[test]
    fn grid_is_initialized_with_correct_shape() {
        let shape = [1, 2, 3, 4];
        let grid = Grid::with_shape(&shape);

        assert_eq!(grid.shape.as_ref(), shape);
    }

    #[test]
    fn grid_is_initialized_with_correct_number_of_bins_for_each_dimension() {
        let shape = [3, 5, 7, 11];

        let grid = Grid::with_shape(&shape);
        assert_eq!(grid.data.len(), 3 * 5 * 7 * 11);
    }

    #[test]
    fn grid_is_initialized_with_none_in_bins() {
        let shape = [1, 2, 3, 4];
        let grid = Grid::with_shape(&shape);

        assert!(grid.data.iter().all(|d| d.is_none()));
    }

    #[test]
    fn spacing_is_set_from_size_and_shape() {
        let shape = [1, 2, 3, 4];
        let size = [4.0, 4.0, 4.0, 4.0];

        let grid = Grid::new(&shape, &size).unwrap();

        assert_eq!(&grid.spacing, &[4.0, 2.0, 4.0 / 3.0, 1.0]);
    }

    #[test]
    fn grid_indexing_is_last_to_first_order() {
        let shape = [3, 5, 7];
        let grid = Grid::with_shape(&shape);

        assert_eq!(grid.get_index(&[0, 0, 0]).unwrap(), 0);
        assert_eq!(grid.get_index(&[0, 0, 1]).unwrap(), 1);
        assert_eq!(grid.get_index(&[0, 0, 6]).unwrap(), 6);

        assert_eq!(grid.get_index(&[0, 1, 0]).unwrap(), 7);
        assert_eq!(grid.get_index(&[0, 1, 1]).unwrap(), 8);
        assert_eq!(grid.get_index(&[0, 1, 6]).unwrap(), 13);
        assert_eq!(grid.get_index(&[0, 4, 6]).unwrap(), 34);

        assert_eq!(grid.get_index(&[1, 0, 0]).unwrap(), 35);
        assert_eq!(grid.get_index(&[2, 4, 6]).unwrap(), 3 * 5 * 7 - 1);
    }

    #[test]
    fn indexing_outside_of_shape_yields_none() {
        let shape = [3, 5, 7];
        let grid = Grid::with_shape(&shape);

        assert!(grid.get_index(&[0, 0, 7]).is_none());
        assert!(grid.get_index(&[0, 5, 0]).is_none());
        assert!(grid.get_index(&[3, 0, 0]).is_none());
    }

    #[test]
    fn position_to_index_is_consistent_with_index_to_position() {
        let shape = [3, 5, 7];
        let grid = Grid::with_shape(&shape);

        // Check that the indices yielded by the  positions were tested in the `get_index` test
        // produce the original positions with the `get_position` method
        for position in &[
            [0, 0, 0],
            [0, 0, 1],
            [0, 0, 6],
            [0, 1, 0],
            [0, 1, 1],
            [0, 1, 6],
            [0, 4, 6],
            [1, 0, 0],
            [2, 4, 6],
        ] {
            let index = grid.get_index(position).unwrap();
            assert_eq!(grid.get_position(index).unwrap(), position);
        }
    }

    #[test]
    fn position_from_out_of_bounds_index_yields_none() {
        let shape = [3, 5, 7];
        let grid = Grid::with_shape(&shape);

        assert!(grid.get_position(grid.data.len()).is_none());
    }

    #[test]
    fn dimension_multipliers_correspond_to_last_to_first_order() {
        let shape = [3, 5, 7];
        let [_, ny, nz] = shape;

        let dim_multipliers = create_dimension_multipliers(&shape);

        assert_eq!(dim_multipliers, [ny * nz, nz, 1]);
    }

    #[test]
    fn getting_position_from_coordinate_works_for_in_box_coords() {
        let shape = [2, 4];
        let size = [10.0, 10.0];

        let grid = Grid::new(&shape, &size).unwrap();

        assert_eq!(grid.get_position_from_coord(&[0.0, 0.0]).unwrap(), &[0, 0]);
        assert_eq!(grid.get_position_from_coord(&[4.9, 2.4]).unwrap(), &[0, 0]);
        assert_eq!(grid.get_position_from_coord(&[5.1, 2.4]).unwrap(), &[1, 0]);
        assert_eq!(grid.get_position_from_coord(&[5.1, 2.6]).unwrap(), &[1, 1]);
        assert_eq!(grid.get_position_from_coord(&[7.5, 5.1]).unwrap(), &[1, 2]);
        assert_eq!(grid.get_position_from_coord(&[7.5, 7.6]).unwrap(), &[1, 3]);
    }

    #[test]
    fn getting_position_from_out_of_box_coordinate_yields_none() {
        let shape = [2, 4];
        let size = [10.0, 10.0];

        let grid = Grid::new(&shape, &size).unwrap();

        assert!(grid.get_position_from_coord(&[-0.5, 0.0]).is_none());
        assert!(grid.get_position_from_coord(&[0.0, -0.5]).is_none());
        assert!(grid.get_position_from_coord(&[10.1, 5.0]).is_none());
        assert!(grid.get_position_from_coord(&[5.0, 10.1]).is_none());
    }

    #[test]
    fn getting_index_from_coordinate_is_consistent_with_from_position() {
        let shape = [2, 4];
        let size = [10.0, 10.0];

        let grid = Grid::new(&shape, &size).unwrap();

        for coord in &[
            [0.0, 0.0],
            [4.9, 2.4],
            [5.1, 2.4],
            [5.1, 2.6],
            [7.5, 5.1],
            [7.5, 7.6],
        ] {
            let position = grid.get_position_from_coord(coord).unwrap();
            assert_eq!(grid.get_index(&position), grid.get_index_from_coord(coord));
        }

        for coord in &[[-5.0, 0.0], [0.0, -5.0], [15.0, 0.0], [0.0, 15.0]] {
            assert!(grid.get_index_from_coord(coord).is_none());
        }
    }

    #[test]
    fn num_adjacent_bins_is_square_root_of_dims() {
        assert_eq!(Grid::with_shape(&vec![1; 1]).num_adjacent, 1);
        assert_eq!(Grid::with_shape(&vec![1; 2]).num_adjacent, 2);
        assert_eq!(Grid::with_shape(&vec![1; 3]).num_adjacent, 2);
        assert_eq!(Grid::with_shape(&vec![1; 4]).num_adjacent, 2);
        assert_eq!(Grid::with_shape(&vec![1; 5]).num_adjacent, 3);
        assert_eq!(Grid::with_shape(&vec![1; 6]).num_adjacent, 3);
    }
}
