#[derive(Clone, Debug, PartialEq)]
pub struct PbcInfo<const D: usize> {
    pub box_size: [f64; D],
    pub inv_box_size: [f64; D],
}

impl<const D: usize> PbcInfo<D> {
    pub fn new(box_size: [f64; D]) -> Self {
        let mut inv_box_size = box_size.clone();
        inv_box_size.iter_mut().for_each(|v| *v = 1.0 / *v);

        PbcInfo {
            box_size,
            inv_box_size,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Grid<const D: usize> {
    /// Conversion multipliers to convert to and from bin positions and the 1d index in `data`.
    dim_multipliers: [usize; D],
    /// Number of bins around each bin which are within range of the minimum distance.
    pub(crate) num_adjacent: usize,
    /// Shape of grid.
    pub(crate) shape: [usize; D],
    /// Spacing of bins in grid.
    pub(crate) spacing: [f64; D],
    /// Information about the size of the grid, used to calculate PBC information.
    pub(crate) pbc: PbcInfo<D>,
    /// Data of grid as a one-dimensional vector.
    pub data: Vec<Option<usize>>,
}

impl<const D: usize> Grid<D> {
    pub fn new(shape: &[usize; D], box_size: &[f64; D]) -> Self {
        let mut spacing = [0.0; D];
        spacing
            .iter_mut()
            .zip(shape.iter().zip(box_size.iter()))
            .for_each(|(sp, (&n, size))| *sp = size / n as f64);

        Grid {
            dim_multipliers: create_dimension_multipliers(shape),
            num_adjacent: (D as f64).sqrt().ceil() as usize,
            shape: shape.clone(),
            spacing,
            pbc: PbcInfo::new(box_size.clone()),
            data: vec![None; shape.iter().product()],
        }
    }

    /// Return the 1d data index of the bin with the input position.
    ///
    /// Returns `None` if no valid bin is available.
    pub fn get_index(&self, position: &[isize], use_pbc: bool) -> Option<usize> {
        if use_pbc {
            Some(
                position
                    .iter()
                    .zip(self.shape.iter())
                    .map(|(&i, &n)| i.rem_euclid(n as isize))
                    .zip(self.dim_multipliers.iter())
                    .map(|(p, m)| p as usize * m)
                    .sum(),
            )
        } else {
            if position
                .iter()
                .zip(self.shape.iter())
                .any(|(&p, &s)| p < 0 || p >= s as isize)
            {
                None
            } else {
                Some(
                    position
                        .iter()
                        .zip(self.dim_multipliers.iter())
                        .map(|(&p, m)| p as usize * m)
                        .sum(),
                )
            }
        }
    }

    pub fn get_position(&self, index: usize) -> Option<Vec<isize>> {
        if index < self.data.len() {
            Some(
                self.dim_multipliers
                    .iter()
                    .scan(index, |i, d| {
                        let n = *i / d;
                        *i = *i % d;

                        Some(n as isize)
                    })
                    .collect(),
            )
        } else {
            None
        }
    }

    /// Return the grid position of a coordinate.
    pub fn get_position_from_coord(&self, coord: &[f64]) -> Option<Vec<isize>> {
        coord
            .iter()
            .zip(self.spacing.iter())
            .map(|(c, dx)| (c / dx).floor() as isize)
            .zip(self.shape.iter().map(|n| *n as isize))
            .map(|(i, n)| if i >= 0 && i < n { Some(i) } else { None })
            .collect()
    }

    /// Return the 1d data index of a coordinate on the grid.
    pub fn get_index_from_coord(&self, coord: &[f64], use_pbc: bool) -> Option<usize> {
        self.get_position_from_coord(coord)
            .and_then(|position| self.get_index(&position, use_pbc))
    }
}

fn create_dimension_multipliers<const D: usize>(shape: &[usize; D]) -> [usize; D] {
    let mut dim_multipliers = [1; D];

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

    impl<const D: usize> Grid<D> {
        fn with_shape(shape: &[usize; D]) -> Self {
            use std::convert::TryInto;

            let size: [f64; D] = shape
                .iter()
                .map(|v| *v as f64)
                .collect::<Vec<f64>>()
                .try_into()
                .expect("error in test: did not create `size` array consistent with `shape` array");

            Grid::new(shape, &size)
        }
    }

    #[test]
    fn grid_with_shape_is_consistent_with_new() {
        let shape = [1, 2, 3, 4];

        let grid_with_shape = Grid::with_shape(&shape);
        let grid_new = Grid::new(&shape, &[1.0, 2.0, 3.0, 4.0]);

        assert_eq!(grid_with_shape, grid_new);
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

        let grid = Grid::new(&shape, &size);

        assert_eq!(&grid.spacing, &[4.0, 2.0, 4.0 / 3.0, 1.0]);
    }

    #[test]
    fn grid_indexing_is_last_to_first_order() {
        let shape = [3, 5, 7];
        let grid = Grid::with_shape(&shape);

        assert_eq!(grid.get_index(&[0, 0, 0], false).unwrap(), 0);
        assert_eq!(grid.get_index(&[0, 0, 1], false).unwrap(), 1);
        assert_eq!(grid.get_index(&[0, 0, 6], false).unwrap(), 6);

        assert_eq!(grid.get_index(&[0, 1, 0], false).unwrap(), 7);
        assert_eq!(grid.get_index(&[0, 1, 1], false).unwrap(), 8);
        assert_eq!(grid.get_index(&[0, 1, 6], false).unwrap(), 13);
        assert_eq!(grid.get_index(&[0, 4, 6], false).unwrap(), 34);

        assert_eq!(grid.get_index(&[1, 0, 0], false).unwrap(), 35);
        assert_eq!(grid.get_index(&[2, 4, 6], false).unwrap(), 3 * 5 * 7 - 1);
    }

    #[test]
    fn indexing_out_of_shape_uses_pbc_to_get_bin_if_use_pbc_is_set() {
        let shape = [3, 5];
        let grid = Grid::with_shape(&shape);

        assert_eq!(
            grid.get_index(&[3, 0], true),
            grid.get_index(&[0, 0], false)
        );
        assert_eq!(
            grid.get_index(&[4, 0], true),
            grid.get_index(&[1, 0], false)
        );
        assert_eq!(
            grid.get_index(&[5, 0], true),
            grid.get_index(&[2, 0], false)
        );
        assert_eq!(
            grid.get_index(&[6, 0], true),
            grid.get_index(&[0, 0], false)
        );
        assert_eq!(
            grid.get_index(&[-1, 0], true),
            grid.get_index(&[2, 0], false)
        );
        assert_eq!(
            grid.get_index(&[-2, 0], true),
            grid.get_index(&[1, 0], false)
        );
        assert_eq!(
            grid.get_index(&[-3, 0], true),
            grid.get_index(&[0, 0], false)
        );
        assert_eq!(
            grid.get_index(&[-4, 0], true),
            grid.get_index(&[2, 0], false)
        );

        assert_eq!(
            grid.get_index(&[2, 5], true),
            grid.get_index(&[2, 0], false)
        );
        assert_eq!(
            grid.get_index(&[2, 6], true),
            grid.get_index(&[2, 1], false)
        );
        assert_eq!(
            grid.get_index(&[2, 9], true),
            grid.get_index(&[2, 4], false)
        );
        assert_eq!(
            grid.get_index(&[2, 10], true),
            grid.get_index(&[2, 0], false)
        );
        assert_eq!(
            grid.get_index(&[2, -1], true),
            grid.get_index(&[2, 4], false)
        );
        assert_eq!(
            grid.get_index(&[2, -2], true),
            grid.get_index(&[2, 3], false)
        );
        assert_eq!(
            grid.get_index(&[2, -5], true),
            grid.get_index(&[2, 0], false)
        );
        assert_eq!(
            grid.get_index(&[2, -6], true),
            grid.get_index(&[2, 4], false)
        );
    }

    #[test]
    fn indexing_outside_of_shape_yields_none_if_use_pbc_not_set() {
        let shape = [3, 5, 7];
        let grid = Grid::with_shape(&shape);

        assert!(grid.get_index(&[-1, 0, 0], false).is_none());
        assert!(grid.get_index(&[0, -1, 0], false).is_none());
        assert!(grid.get_index(&[0, 0, -1], false).is_none());
        assert!(grid.get_index(&[0, 0, 7], false).is_none());
        assert!(grid.get_index(&[0, 5, 0], false).is_none());
        assert!(grid.get_index(&[3, 0, 0], false).is_none());
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
            // Try with and without pbc
            let index = grid.get_index(position, true).unwrap();
            assert_eq!(grid.get_position(index).unwrap(), position);

            let index = grid.get_index(position, false).unwrap();
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

        let grid = Grid::new(&shape, &size);

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

        let grid = Grid::new(&shape, &size);

        assert!(grid.get_position_from_coord(&[-0.5, 0.0]).is_none());
        assert!(grid.get_position_from_coord(&[0.0, -0.5]).is_none());
        assert!(grid.get_position_from_coord(&[10.1, 5.0]).is_none());
        assert!(grid.get_position_from_coord(&[5.0, 10.1]).is_none());
    }

    #[test]
    fn getting_index_from_coordinate_is_consistent_with_from_position() {
        let shape = [2, 4];
        let size = [10.0, 10.0];

        let grid = Grid::new(&shape, &size);

        for coord in &[
            [0.0, 0.0],
            [4.9, 2.4],
            [5.1, 2.4],
            [5.1, 2.6],
            [7.5, 5.1],
            [7.5, 7.6],
        ] {
            let position = grid.get_position_from_coord(coord).unwrap();
            assert_eq!(
                grid.get_index(&position, true),
                grid.get_index_from_coord(coord, true)
            );
            assert_eq!(
                grid.get_index(&position, false),
                grid.get_index_from_coord(coord, false)
            );
        }

        for coord in &[[-5.0, 0.0], [0.0, -5.0], [15.0, 0.0], [0.0, 15.0]] {
            assert!(grid.get_index_from_coord(coord, true).is_none());
            assert!(grid.get_index_from_coord(coord, false).is_none());
        }
    }

    #[test]
    fn num_adjacent_bins_is_square_root_of_dims() {
        assert_eq!(Grid::with_shape(&[1; 1]).num_adjacent, 1);
        assert_eq!(Grid::with_shape(&[1; 2]).num_adjacent, 2);
        assert_eq!(Grid::with_shape(&[1; 3]).num_adjacent, 2);
        assert_eq!(Grid::with_shape(&[1; 4]).num_adjacent, 2);
        assert_eq!(Grid::with_shape(&[1; 5]).num_adjacent, 3);
        assert_eq!(Grid::with_shape(&[1; 6]).num_adjacent, 3);
    }
}
