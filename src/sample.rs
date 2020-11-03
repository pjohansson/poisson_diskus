use rand::{
    distributions::{Distribution, Uniform},
    seq::IteratorRandom,
    Rng,
};
use rand_distr::StandardNormal;
use std::collections::HashSet;

use crate::coord::{add_coords, Coord};

/// Generator of points inside the annulus of an n-ball.
pub struct NBallGen<const DIM: usize> {
    distance: Uniform<f64>,
}

impl<const DIM: usize> NBallGen<DIM> {
    /// Create a generator which samples in the annulus between [rmin, 2.0 * rmin).
    pub fn new(rmin: f64) -> Self {
        NBallGen {
            distance: Uniform::new(rmin, 2.0 * rmin),
        }
    }

    /// Sample a point in the annulus.
    fn sample<R: Rng>(&mut self, rng: &mut R) -> Coord<DIM> {
        let at_distance = self.distance.sample(rng);

        let mut coords: [f64; DIM] = sample_const_num_values(&StandardNormal, rng);

        let current_distance = coords.iter().map(|v| v.powi(2)).sum::<f64>().sqrt();

        coords
            .iter_mut()
            .for_each(|v| *v = *v * (at_distance / current_distance));

        coords
    }

    /// Generate a coordinate in the annulus around a given point.
    pub fn gen_around<R: Rng>(&mut self, x0: &Coord<DIM>, rng: &mut R) -> Coord<DIM> {
        add_coords(x0, &self.sample(rng))
    }
}

/// Return a random index from the set, or `None` if it is empty.
pub fn get_active_index<R: Rng>(inds: &HashSet<usize>, rng: &mut R) -> Option<usize> {
    inds.iter().choose(rng).cloned()
}

/// Generate a random point inside the box as an initial value for the Bridson algorithm.
pub fn gen_init_coord<R: Rng, const DIM: usize>(box_size: &Coord<DIM>, rng: &mut R) -> Coord<DIM> {
    let mut xinit = [0.0; DIM];

    xinit
        .iter_mut()
        .zip(box_size.iter())
        .for_each(|(x, length)| *x = length * rng.gen::<f64>());

    xinit
}

/// Sample N values from a distribution and return as an array.
fn sample_const_num_values<T: Copy + Default, D: Distribution<T>, R: Rng, const N: usize>(
    distr: &D,
    rng: &mut R,
) -> [T; N] {
    let mut coord = [T::default(); N];

    for (c, v) in coord.iter_mut().zip(distr.sample_iter(rng)) {
        *c = v;
    }

    coord
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::coord::calc_distance;

    const NUM_ROUNDS: usize = 1000;

    #[test]
    fn n_ball_samples_points_from_annulus_in_1_dim() {
        type CoordN = Coord<1>;

        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        let mut sphere_gen = NBallGen::new(rmin);

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);
        let x0: CoordN = sample_const_num_values(&x0_range, &mut rng);
        assert_eq!(x0.len(), 1); // Sanity check

        for _ in 0..NUM_ROUNDS {
            let x1 = sphere_gen.gen_around(&x0, &mut rng);
            let distance = calc_distance(&x0, &x1);

            assert!(distance >= rmin && distance < rmax);
        }
    }

    #[test]
    fn n_ball_samples_points_from_annulus_in_2_dims() {
        type CoordN = Coord<2>;

        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        let mut sphere_gen = NBallGen::new(rmin);

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);
        let x0: CoordN = sample_const_num_values(&x0_range, &mut rng);
        assert_eq!(x0.len(), 2); // Sanity check

        for _ in 0..NUM_ROUNDS {
            let x1 = sphere_gen.gen_around(&x0, &mut rng);
            let distance = calc_distance(&x0, &x1);

            assert!(distance >= rmin && distance < rmax);
        }
    }

    #[test]
    fn n_ball_samples_points_from_annulus_in_3_dims() {
        type CoordN = Coord<3>;

        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        let mut sphere_gen = NBallGen::new(rmin);

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);
        let x0: CoordN = sample_const_num_values(&x0_range, &mut rng);
        assert_eq!(x0.len(), 3); // Sanity check

        for _ in 0..NUM_ROUNDS {
            let x1 = sphere_gen.gen_around(&x0, &mut rng);
            let distance = calc_distance(&x0, &x1);

            assert!(distance >= rmin && distance < rmax);
        }
    }

    #[test]
    fn n_ball_samples_points_from_annulus_in_4_dims() {
        type CoordN = Coord<4>;

        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        let mut sphere_gen = NBallGen::new(rmin);

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);
        let x0: CoordN = sample_const_num_values(&x0_range, &mut rng);
        assert_eq!(x0.len(), 4); // Sanity check

        for _ in 0..NUM_ROUNDS {
            let x1 = sphere_gen.gen_around(&x0, &mut rng);
            let distance = calc_distance(&x0, &x1);

            assert!(distance >= rmin && distance < rmax);
        }
    }

    #[test]
    fn n_ball_samples_points_from_annulus_in_5_dims() {
        type CoordN = Coord<5>;

        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        let mut sphere_gen = NBallGen::new(rmin);

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);
        let x0: CoordN = sample_const_num_values(&x0_range, &mut rng);
        assert_eq!(x0.len(), 5); // Sanity check

        for _ in 0..NUM_ROUNDS {
            let x1 = sphere_gen.gen_around(&x0, &mut rng);
            let distance = calc_distance(&x0, &x1);

            assert!(distance >= rmin && distance < rmax);
        }
    }

    #[test]
    fn n_ball_samples_points_from_annulus_in_6_dims() {
        type CoordN = Coord<6>;

        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        let mut sphere_gen = NBallGen::new(rmin);

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);
        let x0: CoordN = sample_const_num_values(&x0_range, &mut rng);
        assert_eq!(x0.len(), 6); // Sanity check

        for _ in 0..NUM_ROUNDS {
            let x1 = sphere_gen.gen_around(&x0, &mut rng);
            let distance = calc_distance(&x0, &x1);

            assert!(distance >= rmin && distance < rmax);
        }
    }

    #[test]
    fn sample_const_values_works_for_0_dim() {
        let mut rng = rand::thread_rng();
        let distr = Uniform::new(512.0_f64, 1024.0_f64);

        for _ in 0..NUM_ROUNDS {
            assert_eq!(sample_const_num_values(&distr, &mut rng), []);
        }
    }

    #[test]
    fn sample_const_values_works_for_32_dim() {
        let mut rng = rand::thread_rng();

        let min = 512.0;
        let max = 1024.0;

        let distr = Uniform::new(min, max);

        for _ in 0..NUM_ROUNDS {
            let sample: [f64; 32] = sample_const_num_values(&distr, &mut rng);

            assert!(sample.iter().all(|&v| v >= min && v < max));
        }
    }

    #[test]
    fn generating_initial_point_uses_box_size() {
        let xmax = 10.0;
        let ymax = 20.0;
        let box_size = [xmax, ymax];

        let mut rng = rand::thread_rng();

        let samples = (0..100)
            .map(|_| gen_init_coord(&box_size, &mut rng))
            .collect::<Vec<_>>();

        // Assert that all samples are inside box
        assert!(samples
            .iter()
            .all(|&[x, y]| x >= 0.0 && x < xmax && y >= 0.0 && y < ymax));

        // The probability of not generating a point inside the any specific half-corner
        // over 100 rounds is (3/4)^100 = 3.2e-13. Thus we use this criteria to assert
        // that we are generating over the entire box.
        let xmid = xmax / 2.0;
        let ymid = ymax / 2.0;

        assert!(samples.iter().any(|&[x, y]| x < xmid && y < ymid)); // lower-left
        assert!(samples.iter().any(|&[x, y]| x >= xmid && y < ymid)); // lower-right
        assert!(samples.iter().any(|&[x, y]| x < xmid && y >= ymid)); // upper-left
        assert!(samples.iter().any(|&[x, y]| x >= xmid && y >= ymid)); // upper-right
    }
}
