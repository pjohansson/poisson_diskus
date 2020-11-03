use poisson_diskus::*;

const NUM_ROUNDS: usize = 30;
const NUM_ROUNDS_LARGEDIM: usize = 4;

fn calc_distance(x0: &[f64], x1: &[f64]) -> f64 {
    x0.iter()
        .zip(x1.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn check_coords<const D: usize>(coords: &[[f64; D]], rmin: f64) {
    eprintln!("num_coords: {}", coords.len());
    assert!(coords.len() > 2);

    for i in 0..coords.len() {
        let x0 = coords.get(i).unwrap();

        for j in (i + 1)..coords.len() {
            let x1 = coords.get(j).unwrap();

            assert!(calc_distance(x0, x1) >= rmin);
        }
    }
}

#[test]
fn one_dimension() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS {
        let coords = bridson(&[100.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn two_dimensions() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS {
        let coords = bridson(&[10.0, 10.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn two_dimensions_skewed() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS {
        let coords = bridson(&[3.0, 37.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn three_dimensions() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS {
        let coords = bridson(&[5.0, 5.0, 5.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn three_dimensions_skewed() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS {
        let coords = bridson(&[3.0, 11.0, 7.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn four_dimensions() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS {
        let coords = bridson(&[3.0, 3.0, 3.0, 3.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn five_dimensions() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS_LARGEDIM {
        let coords = bridson(&[2.5, 2.5, 2.5, 2.5, 2.5], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}

#[test]
fn six_dimensions() {
    let rmin = 1.0;

    for _ in 0..NUM_ROUNDS_LARGEDIM {
        let coords = bridson(&[3.0, 2.0, 2.0, 2.0, 2.0, 2.0], rmin, 30, true).unwrap();
        check_coords(&coords, rmin);
    }
}
