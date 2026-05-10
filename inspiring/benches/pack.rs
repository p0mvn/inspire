//! Phase 10 benchmarks — paper Table 5 reproduction.
//!
//! Phase 4 status: scaffold only. The benchmark groups for the two
//! parameter sets (`(log d, log q, log p, ℓ, z) = (10, 28, 6, 8, 2^4)` and
//! `(11, 56, 15, 3, 2^19)`) will be filled in alongside the
//! [`bench/REPORT.md`](../bench/REPORT.md) document in Phase 10.

use criterion::{criterion_group, criterion_main, Criterion};

fn pack_param_set_1_skeleton(c: &mut Criterion) {
    c.bench_function("pack/skeleton/param_set_1", |b| {
        b.iter(|| {
            // Phase 10 fills in the actual benchmark body. For Phase 4 we
            // emit a no-op so `cargo bench --no-run` succeeds and the
            // criterion harness wiring is exercised.
            std::hint::black_box(0u64)
        });
    });
}

fn pack_param_set_2_skeleton(c: &mut Criterion) {
    c.bench_function("pack/skeleton/param_set_2", |b| {
        b.iter(|| std::hint::black_box(0u64));
    });
}

criterion_group!(
    benches,
    pack_param_set_1_skeleton,
    pack_param_set_2_skeleton
);
criterion_main!(benches);
