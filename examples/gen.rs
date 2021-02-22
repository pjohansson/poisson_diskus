use poisson_diskus::bridson;
use std::time::Instant;

fn main() {
    let start = Instant::now();

    let coords = bridson(&[20.0, 10.0], 0.5, 30, true).unwrap();

    let duration = start.elapsed();
    eprintln!(
        "generated {} points in {} seconds",
        coords.len(),
        duration.as_secs_f64()
    );

    for coord in coords {
        for v in coord {
            print!("{} ", v);
        }

        println!("");
    }
}
