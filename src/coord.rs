use crate::grid::PbcInfo;

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

/// Return the distance between two coordinates, with periodic boundary conditions taken into account.
pub fn calc_distance_pbc(x0: &[f64], x1: &[f64], pbc: &PbcInfo) -> f64 {
    // x0.iter()
    //     .zip(x1.iter())
    //     .zip(pbc.box_size.iter().zip(pbc.inv_box_size.iter()))
    //     .map(|((a, b), (size, inv_size))| {
    //         let dx = b - a;
    //         dx - size * (dx * inv_size).round()
    //     })
    calc_distance_vec_pbc(x0, x1, pbc)
        .iter()
        .map(|v| v.powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Return x1 - x0 with periodic boundary conditions taken into account.
pub fn calc_distance_vec_pbc(x0: &[f64], x1: &[f64], pbc: &PbcInfo) -> Vec<f64> {
    x0.iter()
        .zip(x1.iter())
        .zip(pbc.box_size.iter().zip(pbc.inv_box_size.iter()))
        .map(|((a, b), (size, inv_size))| {
            let dx = b - a;
            dx - size * (dx * inv_size).round()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_approx_eq {
        ($a:expr, $b:expr) => {
            if $a.len() != $b.len() {
                panic!("not equal length: {:?}, {:?}", $a, $b);
            }

            for (a0, b0) in $a.iter().zip($b.iter()) {
                if b0 - a0 > 1e-6 || a0 - b0 > 1e-6 {
                    panic!("not close enough: {:?}, {:?}", $a, $b);
                }
            }
        };
    }

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

    #[test]
    fn calc_dist_vec_with_pbc_one_dimension_same_image() {
        let box_size = [2.0];
        let pbc = PbcInfo::new(&box_size);

        assert_eq!(&calc_distance_vec_pbc(&[0.5], &[0.5], &pbc), &[0.0]);
        assert_eq!(&calc_distance_vec_pbc(&[0.5], &[1.0], &pbc), &[0.5]);
        assert_eq!(&calc_distance_vec_pbc(&[1.0], &[0.5], &pbc), &[-0.5]);
        assert_approx_eq!(&calc_distance_vec_pbc(&[0.5], &[1.4], &pbc), &[0.9]);
        assert_approx_eq!(&calc_distance_vec_pbc(&[1.4], &[0.5], &pbc), &[-0.9]);
        assert_approx_eq!(&calc_distance_vec_pbc(&[1.6], &[0.5], &pbc), &[0.9]);
        assert_approx_eq!(&calc_distance_vec_pbc(&[0.5], &[1.6], &pbc), &[-0.9]);
    }

    #[test]
    fn calc_dist_vec_with_pbc_one_dimension_next_image() {
        let box_size = [2.0];
        let pbc = PbcInfo::new(&box_size);

        // 3.5 shifts to 1.5 with pbc
        assert_eq!(&calc_distance_vec_pbc(&[1.0], &[3.5], &pbc), &[0.5]);
        assert_eq!(&calc_distance_vec_pbc(&[3.5], &[1.0], &pbc), &[-0.5]);

        // -0.5 shifts to 1.5 with pbc
        assert_eq!(&calc_distance_vec_pbc(&[1.0], &[-0.5], &pbc), &[0.5]);
        assert_eq!(&calc_distance_vec_pbc(&[-0.5], &[1.0], &pbc), &[-0.5]);
    }

    #[test]
    fn calc_dist_vec_with_pbc_one_dimension_second_image() {
        let box_size = [2.0];
        let pbc = PbcInfo::new(&box_size);

        // -3.5 shifts to 0.5 with pbc
        assert_eq!(&calc_distance_vec_pbc(&[-3.5], &[1.0], &pbc), &[0.5]);
        assert_eq!(&calc_distance_vec_pbc(&[1.0], &[-3.5], &pbc), &[-0.5]);

        // 4.5 shifts to 0.5 with pbc
        assert_eq!(&calc_distance_vec_pbc(&[4.5], &[1.0], &pbc), &[0.5]);
        assert_eq!(&calc_distance_vec_pbc(&[1.0], &[4.5], &pbc), &[-0.5]);
    }

    #[test]
    fn calc_dist_vec_with_pbc_two_dimensions() {
        let box_size = [2.0, 4.0];
        let pbc = PbcInfo::new(&box_size);

        assert_eq!(
            &calc_distance_vec_pbc(&[1.0, 2.0], &[3.5, 2.0], &pbc),
            &[0.5, 0.0]
        );
        assert_eq!(
            &calc_distance_vec_pbc(&[1.0, 2.0], &[3.5, 3.0], &pbc),
            &[0.5, 1.0]
        );
        assert_eq!(
            &calc_distance_vec_pbc(&[1.0, 2.0], &[3.5, 1.0], &pbc),
            &[0.5, -1.0]
        );
        assert_eq!(
            &calc_distance_vec_pbc(&[1.0, 2.0], &[3.5, 4.0], &pbc),
            &[0.5, -2.0]
        );
        assert_eq!(
            &calc_distance_vec_pbc(&[1.0, 2.0], &[3.5, 5.0], &pbc),
            &[0.5, -1.0]
        );
        assert_eq!(
            &calc_distance_vec_pbc(&[1.0, 7.0], &[3.5, 2.0], &pbc),
            &[0.5, -1.0]
        );
    }
}
