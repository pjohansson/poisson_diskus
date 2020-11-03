use crate::grid::PbcInfo;

pub type Coord<const D: usize> = [f64; D];

/// Add two points of identical dimension.
pub fn add_coords<const D: usize>(x0: &Coord<D>, x1: &Coord<D>) -> Coord<D> {
    let mut buf: Coord<D> = x0.clone();

    buf.iter_mut().zip(x1.iter()).for_each(|(a, b)| *a += b);

    buf
}

/// Calculate the Euclidian distance between two points of identical dimension.
pub fn calc_distance<const D: usize>(x0: &Coord<D>, x1: &Coord<D>) -> f64 {
    x0.iter()
        .zip(x1.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Return the distance between two coordinates, with periodic boundary conditions taken into account.
pub fn calc_distance_pbc(x0: &[f64], x1: &[f64], pbc: &PbcInfo) -> f64 {
    x0.iter()
        .zip(x1.iter())
        .zip(pbc.box_size.iter().zip(pbc.inv_box_size.iter()))
        .map(|((a, b), (size, inv_size))| {
            let dx = b - a;
            dx - size * (dx * inv_size).round()
        })
        .map(|v| v.powi(2))
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_coords_works_as_expected_0_dims() {
        assert_eq!(&add_coords(&[], &[]), &[]);
    }

    #[test]
    fn add_coords_works_as_expected_2_dims() {
        assert_eq!(&add_coords(&[3.0, 5.0], &[7.0, 11.0]), &[10.0, 16.0]);
    }

    #[test]
    fn add_coords_works_as_expected_6_dims() {
        assert_eq!(
            &add_coords(
                &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                &[6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
            ),
            &[6.0, 8.0, 10.0, 12.0, 14.0, 16.0]
        );
    }

    #[test]
    fn calc_distance_works_as_expected_0_dims() {
        assert_eq!(calc_distance(&[], &[]), 0.0);
    }

    #[test]
    fn calc_distance_works_as_expected_2_dims() {
        let x0 = [1.0, 2.0];
        let x1 = [2.0, 3.0];

        assert_eq!(calc_distance(&x0, &x1), 2.0_f64.sqrt());
    }

    #[test]
    fn calc_distance_with_pbc_works_as_expected() {
        let dx = 2.0;
        let dy = 4.0;
        let box_size = [dx, dy];
        let pbc = PbcInfo::new(&box_size);

        let x0 = 1.0;
        let y0 = 2.0;

        assert_eq!(
            calc_distance_pbc(&[x0, y0], &[2.0, 3.0], &pbc),
            2.0_f64.sqrt()
        );
        assert_eq!(
            calc_distance_pbc(&[x0 + dx, y0], &[2.0, 3.0], &pbc),
            2.0_f64.sqrt()
        );
        assert_eq!(
            calc_distance_pbc(&[x0 + 5.0 * dx, y0], &[2.0, 3.0], &pbc),
            2.0_f64.sqrt()
        );
        assert_eq!(
            calc_distance_pbc(&[x0 + 5.0 * dx, y0 + 3.0 * dy], &[2.0, 3.0], &pbc),
            2.0_f64.sqrt()
        );
        assert_eq!(
            calc_distance_pbc(&[x0 - 5.0 * dx, y0 - 3.0 * dy], &[2.0, 3.0], &pbc),
            2.0_f64.sqrt()
        );
    }

    #[test]
    fn calc_distance_works_as_expected_6_dims() {
        let x0 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let x1 = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

        assert_eq!(calc_distance(&x0, &x1), 6.0_f64.sqrt());
    }
}
