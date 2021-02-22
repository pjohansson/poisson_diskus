use rand::Rng;
use std::collections::HashSet;

use crate::{
    coord::calc_distance,
    error::Error,
    grid::Grid,
    sample::{gen_init_coord, get_active_index, NBallGen},
};

type Coord = Vec<f64>;

/// Generate samples from a Poisson disc distribution within the given box.
///
/// The `box_size` array may have any non-zero length. Samples are generated for the given number
/// of dimensions and are separated from each other by a minimum distance `rmin`. A set number of
/// attempts `num_attempts` is made for each sample candidate (30 is suggested as a good value by
/// Bridson, but this can be increased to produce tighter samples).
///
/// This function uses `rand::thread_rng()` as a random number generator. To use another generator,
/// use the [`bridson_rng`] function.
pub fn bridson(box_size: &[f64], rmin: f64, num_attempts: usize) -> Result<Vec<Coord>, Error> {
    let mut rng = rand::thread_rng();

    bridson_rng(box_size, rmin, num_attempts, &mut rng)
}

/// Generate samples from a Poisson disc distribution using a specific random number generator.
///
/// See [`bridson`] for more information.
pub fn bridson_rng<R: Rng>(
    box_size: &[f64],
    rmin: f64,
    num_attempts: usize,
    rng: &mut R,
) -> Result<Vec<Coord>, Error> {
    // Validate input numbers as positive and bounded
    validate_rmin(rmin)?;
    validate_box_size(box_size)?;

    if box_size.is_empty() {
        return Ok(vec![]);
    }

    let shape = get_grid_shape(rmin, box_size);
    let mut grid = Grid::new(&shape, box_size).map_err(|_| Error::UnmatchedDims)?;

    let mut sphere_gen = NBallGen::new(rmin, box_size.len());

    let x0 = gen_init_coord(box_size, rng);
    let grid_index = grid
        .get_index_from_coord(&x0)
        .ok_or(Error::GenCoordOutOfBounds(x0.clone()))?;

    let mut active_inds = HashSet::new();
    let mut samples = Vec::new();

    add_sample_to_list_and_grid(x0, grid_index, &mut samples, &mut active_inds, &mut grid);

    while let Some(grid_index) = get_active_index(&active_inds, rng) {
        let x0 =
            get_sample_from_grid(grid_index, &samples, &grid).ok_or(Error::InvalidActiveList)?;

        match get_sample_around(
            x0,
            &samples,
            &grid,
            num_attempts,
            rmin,
            &mut sphere_gen,
            rng,
        ) {
            Some(coord) => {
                let sample_grid_index = grid
                    .get_index_from_coord(&coord)
                    .ok_or(Error::GenCoordOutOfBounds(coord.clone()))?;

                add_sample_to_list_and_grid(
                    coord,
                    sample_grid_index,
                    &mut samples,
                    &mut active_inds,
                    &mut grid,
                );
            }
            None => {
                active_inds.remove(&grid_index);
            }
        }
    }

    Ok(samples)
}

fn add_sample_to_list_and_grid(
    coord: Coord,
    grid_index: usize,
    samples: &mut Vec<Vec<f64>>,
    active_inds: &mut HashSet<usize>,
    grid: &mut Grid,
) {
    let sample_index = samples.len();

    samples.push(coord);
    active_inds.insert(grid_index);
    grid.data[grid_index] = Some(sample_index);
}

/// Get the coordinate sample if it exists in the grid.
fn get_sample_from_grid<'a>(
    grid_index: usize,
    samples: &'a [Coord],
    grid: &Grid,
) -> Option<&'a Coord> {
    grid.data
        .get(grid_index)
        .cloned()
        .flatten()
        .and_then(|sample_index| samples.get(sample_index))
}

fn get_sample_around<R: Rng>(
    x0: &Coord,
    samples: &[Coord],
    grid: &Grid,
    num_attempts: usize,
    rmin: f64,
    sphere_gen: &mut NBallGen,
    rng: &mut R,
) -> Option<Coord> {
    for _ in 0..num_attempts {
        let x1 = sphere_gen.gen_around(x0, rng);

        if check_if_coord_is_valid(&x1, samples, grid, rmin) {
            return Some(x1);
        }
    }

    None
}

fn check_if_coord_is_valid(coord: &Coord, samples: &[Coord], grid: &Grid, rmin: f64) -> bool {
    match grid.get_position_from_coord(coord) {
        Some(position) => {
            let index_ranges: Vec<(usize, usize)> = position
                .iter()
                .map(|&i| (i.saturating_sub(grid.num_adjacent), i + grid.num_adjacent))
                .collect();

            let mut position_buf = Vec::with_capacity(position.len());

            recurse_and_check(
                &mut position_buf,
                &index_ranges,
                coord,
                samples,
                grid,
                &position,
                rmin,
            )
        }
        None => false,
    }
}

/// Recurse through the position range of all dimensions and verify the grid position.
fn recurse_and_check(
    position: &mut Vec<usize>,
    index_ranges: &[(usize, usize)],
    coord: &Coord,
    samples: &[Coord],
    grid: &Grid,
    original_position: &[usize],
    rmin: f64,
) -> bool {
    match index_ranges.split_first() {
        Some((&(imin, imax), tail)) => {
            for i in imin..=imax {
                position.push(i);

                let result = if tail.is_empty() {
                    check_coord_at_position(coord, position.as_ref(), samples, grid, rmin)
                } else {
                    recurse_and_check(
                        position,
                        tail,
                        coord,
                        samples,
                        grid,
                        original_position,
                        rmin,
                    )
                };

                // We are searching for any sample that is too close to the given coordinate,
                // thus we break the loop and return false at the first sight of one.
                if !result {
                    return false;
                }

                position.pop();
            }
        }
        None => (),
    }

    true
}

/// If the grid has a sample at the position, check if it is too close.
///
/// This function returns true if the grid does not have a sample at the position
/// or if the position is further away from the given coordinate than the minimum
/// distance. Only if there is a sample and it is closer to the coordinate than
/// the minimum distance is false returned, since we are excluding such points
/// from the output.
fn check_coord_at_position(
    coord: &Coord,
    grid_position: &[usize],
    samples: &[Coord],
    grid: &Grid,
    rmin: f64,
) -> bool {
    match grid
        .get_index(grid_position)
        .and_then(|grid_index| get_sample_from_grid(grid_index, samples, grid))
    {
        Some(existing_coord) => calc_distance(coord, existing_coord) >= rmin,
        None => true,
    }
}

fn get_grid_shape(rmin: f64, box_size: &[f64]) -> Vec<usize> {
    get_max_bin_size(rmin, box_size)
        .iter()
        .map(|v| v.ceil() as usize)
        .collect()
}

fn get_max_bin_size(rmin: f64, box_size: &[f64]) -> Vec<f64> {
    let ndim = box_size.len();
    let max_size = rmin / (ndim as f64).sqrt();

    box_size.iter().map(|v| v / max_size).collect()
}

/*************************
 * User input validation *
 *************************/

fn validate_number(value: f64) -> bool {
    value > 0.0 && value.is_finite()
}

fn validate_rmin(rmin: f64) -> Result<(), Error> {
    if validate_number(rmin) {
        Ok(())
    } else {
        Err(Error::InvalidRmin(rmin))
    }
}

fn validate_box_size(box_size: &[f64]) -> Result<(), Error> {
    for &value in box_size {
        if !validate_number(value) {
            return Err(Error::InvalidBoxSize {
                value,
                box_size: box_size.to_vec(),
            });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn add_sample_at_grid_position(position: &[usize], samples: &mut Vec<Coord>, grid: &mut Grid) {
        let coord = position
            .iter()
            .zip(grid.spacing.iter())
            .map(|(&i, dx)| (i as f64 + 0.5) * dx)
            .collect();

        let grid_index = grid.get_index(position).unwrap();
        let mut buf = HashSet::new();

        add_sample_to_list_and_grid(coord, grid_index, samples, &mut buf, grid);
    }

    #[test]
    fn checking_coordinate_in_empty_grid_returns_true() {
        let shape = [4, 4];
        let size = [8.0, 8.0];

        let coord = vec![3.0, 3.0]; // bin: [1, 1]
        let rmin = 2.0 * 2.0_f64.sqrt();

        let grid = Grid::new(&shape, &size).unwrap();

        assert!(check_if_coord_is_valid(&coord, &[], &grid, rmin))
    }

    #[test]
    fn checking_coordinate_in_grid_with_distant_coord_returns_true() {
        let shape = [4, 4];
        let size = [8.0, 8.0];

        let coord = vec![3.0, 3.0]; // bin: [1, 1]
        let rmin = 2.0 * 2.0_f64.sqrt();

        let mut grid = Grid::new(&shape, &size).unwrap();
        let mut samples = Vec::new();

        add_sample_at_grid_position(&[3, 3], &mut samples, &mut grid);

        assert!(check_if_coord_is_valid(&coord, &samples, &grid, rmin))
    }

    #[test]
    fn checking_coordinate_uses_num_adjacent_for_its_search_space() {
        // Same test as above, but increase the num_adjacent value to look in the distant bin
        let shape = [4, 4];
        let size = [8.0, 8.0];

        let coord = vec![1.0, 1.0]; // bin: [0, 0]
        let rmin = 2.0 * 2.0_f64.sqrt();

        let mut grid = Grid::new(&shape, &size).unwrap();

        let mut samples = Vec::new();
        add_sample_at_grid_position(&[3, 3], &mut samples, &mut grid);

        // Increase num_adjacent to look in [3, 3]
        grid.num_adjacent = 3;

        // Set the coordinate close enough to the candidate (cheating!)
        samples[0] = vec![1.0, 2.0];

        assert!(!check_if_coord_is_valid(&coord, &samples, &grid, rmin))
    }

    #[test]
    fn checking_coordinate_in_grid_with_close_coord_returns_false() {
        let shape = [4, 4];
        let size = [8.0, 8.0];

        let coord = vec![3.0, 3.0]; // bin: [1, 1]
        let rmin = 2.0 * 2.0_f64.sqrt();

        let mut grid = Grid::new(&shape, &size).unwrap();
        let mut samples = Vec::new();

        add_sample_at_grid_position(&[2, 1], &mut samples, &mut grid);

        assert!(!check_if_coord_is_valid(&coord, &samples, &grid, rmin))
    }

    #[test]
    fn checking_coordinate_in_grid_with_coord_in_adjacent_box_but_not_within_rmin_returns_true() {
        let shape = [4, 4];
        let size = [8.0, 8.0];

        let coord = vec![2.01, 2.01]; // bin: [1, 1]
        let rmin = 2.0 * 2.0_f64.sqrt();

        let mut grid = Grid::new(&shape, &size).unwrap();
        let mut samples = Vec::new();

        add_sample_at_grid_position(&[2, 1], &mut samples, &mut grid);
        samples[0] = vec![5.99, 3.99]; // bin: 2, 1

        assert!(check_if_coord_is_valid(&coord, &samples, &grid, rmin))
    }

    /*************************
     * USER INPUT VALIDATION *
     *************************/

    #[test]
    fn non_positive_rmin_yields_error() {
        let box_size = [5.0, 5.0];
        let num_attempts = 10;

        assert!(bridson(&box_size, 1.0, 10).is_ok());

        let rmin_zero = 0.0;
        let rmin_neg = -1.0;
        let rmin_nan = f64::NAN;
        let rmin_inf = f64::INFINITY;
        let rmin_neg_inf = f64::NEG_INFINITY;

        assert!(bridson(&box_size, rmin_zero, num_attempts)
            .unwrap_err()
            .is_invalid_rmin());

        assert!(bridson(&box_size, rmin_neg, num_attempts)
            .unwrap_err()
            .is_invalid_rmin());

        assert!(bridson(&box_size, rmin_nan, num_attempts)
            .unwrap_err()
            .is_invalid_rmin());

        assert!(bridson(&box_size, rmin_inf, num_attempts)
            .unwrap_err()
            .is_invalid_rmin());

        assert!(bridson(&box_size, rmin_neg_inf, num_attempts)
            .unwrap_err()
            .is_invalid_rmin());
    }

    #[test]
    fn empty_box_size_yields_no_coords() {
        assert!(bridson(&[], 1.0, 10).unwrap().is_empty());
    }

    #[test]
    fn non_positive_box_size_yields_error() {
        let rmin = 1.0;
        let num_attempts = 10;

        assert!(bridson(&[5.0, 5.0], rmin, num_attempts).is_ok());

        // Use invalid box size values at the front and back of the given
        // array, to ensure that all are tested.
        assert!(bridson(&[-5.0, 5.0], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[5.0, -5.0], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[-5.0, -5.0], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[f64::NAN, 5.0], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[5.0, f64::NAN], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[f64::INFINITY, 5.0], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[5.0, f64::INFINITY], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[5.0, f64::NEG_INFINITY], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());

        assert!(bridson(&[f64::NEG_INFINITY, 5.0], rmin, num_attempts)
            .unwrap_err()
            .is_invalid_box_size());
    }

    #[test]
    fn non_positive_num_attempts_works() {
        let box_size = [5.0, 5.0];
        let rmin = 1.0;

        assert!(bridson(&box_size, rmin, 0).is_ok());
        assert!(bridson(&box_size, rmin, 1).is_ok());
        assert!(bridson(&box_size, rmin, 10).is_ok());
    }
}
