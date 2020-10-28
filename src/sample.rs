use rand::{
    distributions::{Distribution, Uniform},
    seq::IteratorRandom,
    Rng,
};
use rand_distr::StandardNormal;
use std::collections::HashSet;

use crate::coord::add_coords;

/// Generator of points inside the annulus of an n-ball.
pub struct NBallGen {
    distance: Uniform<f64>,
    ndim: usize,
}

impl NBallGen {
    /// Create a generator which samples in the annulus between [rmin, 2.0 * rmin).
    pub fn new(rmin: f64, ndim: usize) -> Self {
        NBallGen {
            distance: Uniform::new(rmin, 2.0 * rmin),
            ndim,
        }
    }

    /// Sample a point in the annulus.
    fn sample<R: Rng>(&mut self, rng: &mut R) -> Vec<f64> {
        let at_distance = self.distance.sample(rng);

        let coords_unscaled = rng
            .sample_iter(StandardNormal)
            .take(self.ndim)
            .collect::<Vec<f64>>();

        let current_distance = coords_unscaled
            .iter()
            .map(|v| v.powi(2))
            .sum::<f64>()
            .sqrt();

        coords_unscaled
            .into_iter()
            .map(|v| v * (at_distance / current_distance))
            .collect()
    }

    /// Generate a coordinate in the annulus around a given point.
    pub fn gen_around<R: Rng>(&mut self, x0: &[f64], rng: &mut R) -> Vec<f64> {
        add_coords(x0, &self.sample(rng))
    }
}

pub fn get_active_index<R: Rng>(inds: &HashSet<usize>, rng: &mut R) -> Option<usize> {
    inds.iter().choose(rng).cloned()
}

pub fn gen_init_coord<R: Rng>(box_size: &[f64], rng: &mut R) -> Vec<f64> {
    box_size
        .iter()
        .map(|length| length * rng.gen::<f64>())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::coord::calc_distance;

    const NUM_ROUNDS: usize = 1000;

    #[test]
    fn n_ball_samples_points_from_annulus() {
        let rmin = 3.0;
        let rmax = 2.0 * rmin;

        let mut rng = rand::thread_rng();

        // Sample around coordinates which are not at origin
        let x0_range = Uniform::new(512.0, 1024.0);

        for ndim in 1..6 {
            let mut sphere_gen = NBallGen::new(rmin, ndim);
            let x0 = x0_range
                .sample_iter(&mut rng)
                .take(ndim)
                .collect::<Vec<_>>();

            for _ in 0..NUM_ROUNDS {
                let x1 = sphere_gen.gen_around(&x0, &mut rng);
                let distance = calc_distance(&x0, &x1);

                assert!(distance >= rmin && distance < rmax);
            }
        }
    }
}
