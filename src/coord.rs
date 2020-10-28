pub fn add_coords(x0: &[f64], x1: &[f64]) -> Vec<f64> {
    x0.iter().zip(x1.iter()).map(|(a, b)| a + b).collect()
}

pub fn calc_distance(x0: &[f64], x1: &[f64]) -> f64 {
    x0.iter()
        .zip(x1.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_coords_works_as_expected() {
        assert_eq!(&add_coords(&[3.0, 5.0], &[7.0, 11.0]), &[10.0, 16.0]);
    }

    #[test]
    fn calc_distance_works_as_expected() {
        let x0 = [1.0, 2.0];
        let x1 = [2.0, 3.0];

        assert_eq!(calc_distance(&x0, &x1), 2.0_f64.sqrt());
    }
}
