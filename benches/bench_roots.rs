use criterion::{Criterion, black_box};
use criterion::{criterion_group, criterion_main};

use roots::roots;

fn criterion_bench(c: &mut Criterion) {
    c.bench_function("roots 1", |b| b.iter(|| roots(black_box(&[3.2, 2.0, 1.0]))));
}

criterion_group!(benches, criterion_bench);
criterion_main!(benches);
