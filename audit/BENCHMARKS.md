# Baseline benchmarks

Measured on CPython 3.12.14 (Linux-6.18.44-x86_64-with-glibc2.41) against source `5257291431bb45db0274dc48edf24694ecfe2e2d`. Values are microseconds per operation, including allocation and result-object construction where the public API does so. Raw samples and calibrated iteration counts are in [benchmark-results.json](benchmark-results.json).

| Operation | Median µs | Min–max µs |
| --- | ---: | ---: |
| `vector3_add` | 0.467 | 0.461–0.496 |
| `vector3_dot` | 0.272 | 0.269–0.291 |
| `vector3_cross` | 0.481 | 0.468–0.497 |
| `vector3_normalize` | 0.732 | 0.719–0.830 |
| `matrix4_raw_multiply` | 5.271 | 5.237–5.348 |
| `matrix4_object_multiply` | 7.740 | 7.619–8.100 |
| `matrix4_vector4_multiply` | 1.723 | 1.656–1.770 |
| `matrix4_inverse` | 6.309 | 6.176–6.578 |
| `quaternion_multiply` | 0.635 | 0.630–0.676 |
| `quaternion_rotate_vector3` | 1.781 | 1.765–1.888 |
| `quaternion_slerp` | 2.119 | 2.094–2.253 |
| `quadratic_bezier_scalar` | 0.191 | 0.185–0.197 |
| `legendre_l2_m0` | 0.616 | 0.564–0.634 |
| `spherical_harmonic_l2_m0` | 0.982 | 0.978–1.010 |

The harness performs a warmup, timeit autorange calibration, and five repetitions per operation. GC is disabled during timing by timeit and restored afterward. Fixed inputs isolate representative small operations; these measurements are not application frame rates, worst-case latency, allocation profiles, or cross-machine performance guarantees. The shared cloud machine was not pinned to an isolated CPU, and frequency/load can change results. Compare the samples, not only a single median.

Broken functions (inverse2, refraction, plane constructors, project, SQUAD, SH sampling/object generation, and cubic Bezier) are intentionally excluded from performance claims. Tested working low-order Legendre and spherical harmonics are benchmarked; higher-order failures are not timed as successful math.

Raw versus public-object 4×4 multiply is useful for identifying possible wrapper/ctypes overhead, but their difference is not a causal profile. No optimization or refactoring was performed. Before comparing a later change, rerun the same harness without coverage on similar hardware/load and confirm the operation remains mathematically correct.
