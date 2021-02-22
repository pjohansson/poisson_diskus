use poisson_diskus::bridson;
use std::time::Instant;

fn calc_distance(x0: &[f64], x1: &[f64], box_size: &[f64], inv_box: &[f64]) -> f64 {
    x0.iter()
        .zip(x1.iter())
        .zip(box_size.iter().zip(inv_box.iter()))
        .map(|((a, b), (size, inv_size))| {
            let dx = b - a;
            dx - size * (dx * inv_size).round()
        })
        .map(|v| v.powi(2))
        .sum::<f64>()
        .sqrt()
}

fn main() {
    // let start = Instant::now();

    let box_size = [100.0, 200.0];
    let rmin = 0.5;
    let num_attempts = 30;
    let use_pbc = true;

    let coords = bridson(&box_size, rmin, num_attempts, use_pbc).unwrap();

    // let duration = start.elapsed();
    // eprintln!(
    //     "generated {} points in {} seconds",
    //     coords.len(),
    //     duration.as_secs_f64()
    // );

    // let max_dist = 5.0 * rmin;
    // let inv_box = box_size.iter().map(|v| 1.0 / v).collect::<Vec<_>>();

    // let num_bins = 400;
    // let dr = max_dist / num_bins as f64;

    // let rs = (0..num_bins).map(|i| i as f64 * dr).collect::<Vec<_>>();
    // let mut rdf = vec![0.0f64; num_bins];

    // let start = Instant::now();

    // 'outer: for i in 0..coords.len() {
    //     let x0 = &coords[i];

    //     for (&x, &size) in x0.iter().zip(box_size.iter()) {
    //         if x < rmin || x + rmin >= size {
    //             continue 'outer;
    //         }
    //     }

    //     for j in (i + 1)..coords.len() {
    //         let x1 = &coords[j];

    //         let distance = calc_distance(x0, x1, &box_size, &inv_box);
    //         let i = (distance / dr).round() as usize;

    //         if i < num_bins {
    //             rdf[i] += 1.0;
    //         }
    //     }
    // }

    // let duration = start.elapsed();
    // eprintln!("calculated rdf in {} seconds", duration.as_secs_f64());

    // for (r, v) in rs.iter().zip(rdf.iter()) {
    //     let area = std::f64::consts::PI * ((r + dr).powi(2) - r.powi(2));
    //     println!("{:12.9} {:12.9}", r + 0.5 * dr, v / area);
    // }
}
