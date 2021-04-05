use criterion::{criterion_group, criterion_main, Criterion};
use poisson_diskus::bridson;

fn criterion_benchmark(c: &mut Criterion) {
    let box_size_2d = [10.0, 10.0];
    let box_size_3d = [10.0, 10.0, 10.0];
    let box_size_4d = [10.0, 10.0, 10.0, 10.0];

    let rmin_2d = 0.5;
    let rmin_3d = 2.0;
    let rmin_4d = 3.0;

    let num_attempts = 30;

    c.bench_function("bridson 2d no pbc", |b| b.iter(|| bridson(&box_size_2d, rmin_2d, num_attempts, false)));
    c.bench_function("bridson 3d no pbc", |b| b.iter(|| bridson(&box_size_3d, rmin_3d, num_attempts, false)));
    c.bench_function("bridson 4d no pbc", |b| b.iter(|| bridson(&box_size_4d, rmin_4d, num_attempts, false)));
    c.bench_function("bridson 2d with pbc", |b| b.iter(|| bridson(&box_size_2d, rmin_2d, num_attempts, true)));
    c.bench_function("bridson 3d with pbc", |b| b.iter(|| bridson(&box_size_3d, rmin_3d, num_attempts, true)));
    c.bench_function("bridson 4d with pbc", |b| b.iter(|| bridson(&box_size_4d, rmin_4d, num_attempts, true)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);