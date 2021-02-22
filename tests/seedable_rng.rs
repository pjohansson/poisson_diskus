use poisson_diskus::*;
use rand::{rngs::StdRng, SeedableRng};

#[test]
fn seedable_rng_can_be_used() {
    let mut rng = StdRng::seed_from_u64(10);

    let box_size = [2.5, 3.5, 4.5];
    let rmin = 1.0;
    let num_attempts = 30;

    assert!(bridson_rng(&box_size, rmin, num_attempts, true, &mut rng).is_ok());
}
