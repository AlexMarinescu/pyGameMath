# Core performance baseline and profiling

Base: master `8cb07d4d6ac96980c25bd698869734005d45d2ca` (merged PR #32).
Branch: `perf/phase3a-baseline`. This phase changes benchmark tooling and audit
artifacts only; mathematical code, APIs, packaging and test markers are unchanged.

## Reproduction and measurement scope

```
/workspace/.venvs/pyGameMath/bin/python benchmarks/core_baseline.py --output audit/phase3a-benchmarks.json --trials 7 --target-seconds .02 --historical --profiles
/workspace/.venvs/pyGameMath/bin/python benchmarks/core_baseline.py --output /tmp/gem-repeat.json --trials 7 --target-seconds .02
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider --junitxml=/tmp/phase3a-tests.xml
```

The [harness guide](../benchmarks/README.md) specifies all deterministic input
ranges, setup boundaries, timing and profiler limitations. No RNG/assets or
benchmark dependencies are added. Historical comparison requires a local git
object; all current measurements also work in source archives without git.
The seven-trial primary run includes 108 current-core cases. Each case warms
once, doubles its calibrated operation count until at least .02 seconds, then
records seven wall-time samples. Expensive single calls exceed this target.
JSON records operation counts, samples, extrema, median absolute deviation
(MAD), setup/import timing and source/harness SHA-256 fingerprints.

GC is disabled only within timeit. Public result allocation and ctypes export
are included; raw-kernel and constructor cases are measured separately. Input
generation is outside timed calls. Ray transforms include a fresh duplicate,
whose cost is separately measured. Fresh-interpreter import timing excludes
process launch and does not flush filesystem caches. No-op timing is reported
without subtraction. cProfile and tracemalloc run separately after timings;
their instrumented runtimes are not latency measurements, and cumulative times
overlap. Traced peak bytes are not total process RSS.

Environment: **CPython 3.12.14**, Linux-6.18.44-x86_64-with-glibc2.41,
x86_64, CPU **INTEL(R) XEON(R) PLATINUM 8573C**, 3 visible logical CPUs,
six 1.17.0. No CPU affinity or frequency/load isolation was imposed.
Setup took 0.0861 seconds; median internal cold
import time was 18.96 ms. These
are not universal throughput or startup guarantees.

## Current-core results

All values below are microseconds per public operation or named raw kernel,
including returned results. MAD and min/max show within-run variability.
Construction overhead cannot be inferred by simply subtracting independent
medians; profiles provide the more direct evidence.

| Operation | Median µs | MAD µs | Min–max µs |
| --- | ---: | ---: | ---: |
| `vector2_allocate` | 0.228 | 0.007 | 0.221–0.264 |
| `vector2_add` | 0.540 | 0.008 | 0.503–0.592 |
| `vector2_subtract` | 0.574 | 0.017 | 0.536–0.631 |
| `vector2_scalar_multiply` | 0.599 | 0.041 | 0.558–0.734 |
| `vector2_dot` | 0.260 | 0.004 | 0.245–0.265 |
| `vector2_magnitude` | 1.047 | 0.043 | 1.004–1.120 |
| `vector2_normalize` | 2.632 | 0.056 | 2.456–2.983 |
| `matrix2_allocate` | 1.885 | 0.153 | 1.596–2.737 |
| `matrix2_raw_multiply` | 1.170 | 0.032 | 1.091–1.211 |
| `matrix2_multiply` | 2.695 | 0.038 | 2.657–2.780 |
| `matrix2_vector` | 0.863 | 0.018 | 0.818–0.918 |
| `matrix2_determinant` | 0.127 | 0.005 | 0.115–0.131 |
| `matrix2_inverse` | 1.933 | 0.016 | 1.857–2.173 |
| `matrix2_raw_inverse` | 0.758 | 0.028 | 0.715–0.855 |
| `matrix2_transpose` | 2.177 | 0.087 | 2.005–2.334 |
| `vector2_transform` | 0.633 | 0.006 | 0.574–0.645 |
| `vector3_allocate` | 0.227 | 0.006 | 0.216–0.239 |
| `vector3_add` | 0.570 | 0.008 | 0.560–0.627 |
| `vector3_subtract` | 0.597 | 0.018 | 0.571–0.664 |
| `vector3_scalar_multiply` | 0.596 | 0.020 | 0.560–0.628 |
| `vector3_dot` | 0.272 | 0.006 | 0.266–0.300 |
| `vector3_magnitude` | 1.288 | 0.040 | 1.210–1.569 |
| `vector3_normalize` | 2.837 | 0.099 | 2.576–3.026 |
| `matrix3_allocate` | 2.377 | 0.081 | 2.175–2.458 |
| `matrix3_raw_multiply` | 2.405 | 0.020 | 2.367–2.678 |
| `matrix3_multiply` | 4.889 | 0.088 | 4.696–5.033 |
| `matrix3_vector` | 1.172 | 0.036 | 1.121–1.284 |
| `matrix3_determinant` | 0.289 | 0.002 | 0.278–0.326 |
| `matrix3_inverse` | 14.967 | 0.433 | 14.420–15.883 |
| `matrix3_raw_inverse` | 14.019 | 0.865 | 12.732–15.152 |
| `matrix3_transpose` | 2.963 | 0.121 | 2.821–3.124 |
| `vector3_transform` | 0.934 | 0.022 | 0.890–0.956 |
| `matrix3_raw_inverse_scale_1e-200` | 19.178 | 2.789 | 14.905–26.574 |
| `matrix3_raw_inverse_scale_1e+200` | 19.258 | 0.497 | 18.384–19.882 |
| `vector4_allocate` | 0.218 | 0.005 | 0.212–0.231 |
| `vector4_add` | 0.550 | 0.015 | 0.531–0.583 |
| `vector4_subtract` | 0.553 | 0.019 | 0.534–0.667 |
| `vector4_scalar_multiply` | 0.588 | 0.012 | 0.576–0.653 |
| `vector4_dot` | 0.299 | 0.011 | 0.288–0.361 |
| `vector4_magnitude` | 1.565 | 0.046 | 1.463–2.881 |
| `vector4_normalize` | 3.380 | 0.110 | 3.261–3.762 |
| `matrix4_allocate` | 3.401 | 0.073 | 3.225–4.038 |
| `matrix4_raw_multiply` | 5.043 | 0.224 | 4.719–5.397 |
| `matrix4_multiply` | 8.146 | 0.151 | 7.758–9.124 |
| `matrix4_vector` | 1.631 | 0.020 | 1.506–1.685 |
| `matrix4_determinant` | 3.380 | 0.079 | 3.300–3.891 |
| `matrix4_inverse` | 33.318 | 0.625 | 32.372–37.112 |
| `matrix4_raw_inverse` | 30.495 | 0.607 | 29.439–35.474 |
| `matrix4_transpose` | 4.185 | 0.110 | 4.074–4.853 |
| `vector4_transform` | 1.389 | 0.040 | 1.317–1.430 |
| `matrix4_raw_inverse_scale_1e-200` | 34.626 | 1.384 | 33.168–42.054 |
| `matrix4_raw_inverse_scale_1e+200` | 106.100 | 4.230 | 100.270–145.692 |
| `vector3_cross` | 0.601 | 0.007 | 0.556–0.643 |
| `vector3_magnitude_large` | 1.523 | 0.098 | 1.367–1.679 |
| `vector3_normalize_large` | 3.821 | 0.564 | 3.257–6.284 |
| `vector3_magnitude_tiny` | 1.343 | 0.066 | 1.248–1.426 |
| `vector3_normalize_tiny` | 3.235 | 0.109 | 2.899–4.195 |
| `vector3_magnitude_mixed` | 1.383 | 0.034 | 1.349–1.588 |
| `vector3_normalize_mixed` | 3.357 | 0.094 | 3.232–3.465 |
| `quaternion_normalize_large` | 4.064 | 0.148 | 3.753–4.236 |
| `quaternion_normalize_tiny` | 3.906 | 0.030 | 3.848–4.387 |
| `quaternion_allocate` | 0.171 | 0.003 | 0.163–0.175 |
| `quaternion_multiply` | 0.622 | 0.032 | 0.576–0.675 |
| `quaternion_normalize` | 3.865 | 0.111 | 3.734–4.580 |
| `quaternion_inverse` | 0.381 | 0.012 | 0.369–0.503 |
| `quaternion_rotate_vector` | 2.211 | 0.286 | 1.840–2.614 |
| `quaternion_slerp` | 3.725 | 0.252 | 3.014–4.138 |
| `quaternion_slerp_near` | 3.614 | 0.198 | 3.336–3.950 |
| `quaternion_squad` | 6.815 | 0.017 | 6.798–7.553 |
| `quaternion_squad4` | 11.259 | 1.802 | 9.457–82.022 |
| `quaternion_to_matrix` | 6.772 | 0.191 | 6.522–7.558 |
| `quaternion_from_matrix` | 0.800 | 0.036 | 0.765–1.082 |
| `matrix4_translate` | 11.852 | 0.682 | 11.067–16.360 |
| `matrix3_rotate` | 11.807 | 0.654 | 11.142–15.142 |
| `vector3_affine_transform` | 1.067 | 0.027 | 1.013–1.115 |
| `bezier_quadratic` | 3.196 | 0.147 | 3.025–4.333 |
| `bezier_cubic` | 4.586 | 0.138 | 4.424–5.211 |
| `bezier_cubic_scalar` | 0.275 | 0.011 | 0.258–0.297 |
| `bezier_subdivide_0.1` | 135.841 | 20.799 | 105.461–177.964 |
| `bezier_subdivide_0.01` | 420.809 | 17.183 | 389.386–459.478 |
| `bezier_subdivide_0.0001` | 3456.741 | 257.643 | 3073.190–3714.384 |
| `bezier_subdivide_1e-08` | 333532.249 | 2507.242 | 327283.135–338661.475 |
| `bezier_subdivide_1e-15` | 1404000.322 | 4849.743 | 1391144.414–1425363.517 |
| `bezier_path_8_segments` | 8278.609 | 295.351 | 7935.586–9922.150 |
| `sh_basis_0_0` | 0.916 | 0.104 | 0.719–1.020 |
| `sh_basis_2_1` | 1.497 | 0.376 | 1.019–17.631 |
| `sh_basis_8_3` | 2.751 | 0.280 | 2.296–29.198 |
| `sh_basis_12_6` | 3.381 | 0.148 | 3.219–4.152 |
| `sh_project_n64_b1` | 68.452 | 4.127 | 64.174–77.728 |
| `sh_project_n64_b3` | 254.986 | 14.618 | 217.994–281.452 |
| `sh_project_n64_b5` | 593.668 | 15.585 | 564.574–830.450 |
| `sh_project_n256_b1` | 254.961 | 6.431 | 240.212–271.932 |
| `sh_project_n256_b3` | 937.027 | 55.350 | 845.568–1040.109 |
| `sh_project_n256_b5` | 2314.475 | 103.021 | 2186.958–2791.172 |
| `sh_project_n1024_b1` | 947.435 | 32.211 | 904.823–1095.607 |
| `sh_project_n1024_b3` | 3769.962 | 118.927 | 3611.861–3930.260 |
| `sh_project_n1024_b5` | 9011.714 | 325.056 | 8490.738–10220.706 |
| `sh_convolve_l2` | 6.467 | 0.229 | 6.235–6.814 |
| `sh_rotate_l2_scalar` | 18.971 | 1.039 | 17.797–23.394 |
| `sh_rotate_l2_rgb` | 45.781 | 2.471 | 42.563–55.157 |
| `sh_reconstruct_l2` | 21.615 | 0.644 | 18.224–23.259 |
| `plane_normalize` | 2.850 | 0.147 | 2.703–3.289 |
| `plane_dot` | 0.140 | 0.003 | 0.134–0.152 |
| `plane_from_points` | 7.682 | 0.297 | 7.064–8.691 |
| `ray_duplicate` | 1.126 | 0.124 | 0.955–1.290 |
| `ray_translate` | 6.663 | 0.377 | 6.116–8.902 |
| `ray_matrix_rotate` | 7.637 | 0.378 | 7.145–8.054 |
| `ray_quaternion_rotate` | 10.028 | 0.379 | 9.131–10.501 |

Ray intersections are **unavailable**: the supported Ray class has no intersection
method, only stored intersection-placeholder state. No fabricated intersection
algorithm or unfinished E07 transport was timed. Plane construction, evaluation
and normalization, and Ray copying/matrix/quaternion/translation operations are
covered. This is a representative core baseline, not an exhaustive application
workload or unsupported-domain benchmark.

## Historical inverse comparison

Historical Matrix code is loaded from Phase 2F-1 commit
`2cbd899a47546b8c40f8244f57799969c86b0d87` into an isolated namespace, then timed
in this same process with the same triangular matrices and seven-trial policy.
The old code imports current helper modules; this comparison isolates the matrix
inverse implementation rather than emulating a complete historical runtime.
The triangular dataset differs from the dense current-core table above; ratios
below use matched inputs only.

| Size | Old raw µs | Current raw µs | Raw ratio | Old wrapper µs | Current wrapper µs | Wrapper ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 1.411 | 9.754 | 6.91× | 3.355 | 11.762 | 3.51× |
| 4 | 5.906 | 21.095 | 3.57× | 8.682 | 22.967 | 2.65× |

This confirms the stabilization cost reported in Phase 2F-2; earlier snapshots
are context, not controlled cross-run speed comparisons. Independent triangular
inverse formulas and both multiplication orders pass at scales 1e-300, 1e-200,
1, 1e200 and 1e300. Historical arithmetic fails or produces inaccurate/nonfinite
answers at extreme scales; statuses and current residuals are in JSON.
The speed difference buys power-of-two scaling and exact represented-coefficient
singularity handling without epsilon thresholds. Removing those protections to
recover old speed would violate the established contract. Ill-conditioning and
unrepresentable outputs remain separate from this well-conditioned scale study.

## Profiles and allocation pressure

The JSON stores the 20 highest cumulative frames for nine selected expensive or
representative paths, with primitive/total calls and self/cumulative times.
Profile operation counts are bounded by a .1-second uninstrumented estimate,
up to 2000 repetitions. Deep Bezier uses one call; its millions of inner calls
still give clear structural evidence. Instrumentation magnifies generator/loop
costs, so percentages are diagnostic rather than unprofiled latency predictions.

| Profiled operation | Unprofiled median µs | Single-call traced peak bytes |
| --- | ---: | ---: |
| `matrix3_raw_inverse` | 14.019 | 1505 |
| `matrix4_raw_inverse` | 30.495 | 3769 |
| `matrix4_multiply` | 8.146 | 856 |
| `quaternion_squad4` | 11.259 | 968 |
| `bezier_subdivide_1e-15` | 1404000.322 | 21047976 |
| `bezier_path_8_segments` | 8278.609 | 107760 |
| `sh_project_n1024_b5` | 9011.714 | 110776 |
| `sh_rotate_l2_rgb` | 45.781 | 4296 |
| `ray_quaternion_rotate` | 10.028 | 1928 |

- **Inverse3/4:** `_scaled_cofactor_inverse` dominates. Finite checks traverse all
  entries; each row constructs ratios/integer rows, exact determinant state,
  exponent/scaled lists and rescaled outputs. `max`, `any`, generators and
  conversions appear repeatedly. The cofactor kernel is only part of the cost.
  Optimize this plumbing only with preserved exact-zero and extreme-scale tests.
- **Matrix4 multiplication:** Python triple loops dominate the raw product;
  `Matrix.__init__` and `conv_list_2d` add result construction/ctypes conversion.
  Public exports cannot simply be removed or left stale.
- **Quaternion SQUAD4:** three SLERPs call multiply/add, dot and constructors
  repeatedly. Temporary quaternion allocation and Python function dispatch
  matter more than a new interpolation formula. Do not regress shortest-path,
  norm accuracy or existing sign branches.
- **Bezier:** each interior control distance recomputes chord/unit/projection
  state; `_chord_distance`, generator tuples and `max` consume most profiled
  time. `_split` creates de Casteljau levels for every internal node; final
  tuple samples are converted to Vectors. Allocation grows with the emitted
  path, so local arithmetic improvement cannot remove the exponential cap cost.
- **SH projection:** three compensated sums per coefficient traverse sample
  generators. Validation also scans RGB values, weights and all precomputed
  basis values. Cost is O(N*B²*3), with O(N*B²) validation; precomputation/setup
  is excluded from projection timing, and basis evaluation is separately timed.
  Compensated `fsum` accuracy is part of the correctness tradeoff.
- **L2 RGB rotation:** small analytical tensors incur repeated nested `sum`
  generators and validation checks. No directional resampling/reintegration is
  appropriate; any optimization must preserve analytical band energy and
  active-rotation contracts.
- **Ray quaternion rotation:** two quaternion sandwiches plus normalization and
  duplication account for most work. Removing normalization would change
  established behavior for nonunit inputs; this profile is not permission to do so.

## Adaptive point growth

The same 3D cubic is sampled with decreasing distance tolerance; the public
parameter stores its square. Input bounds and curve are unchanged.

| Distance tolerance | Points | Median ms | Maximum possible points |
| --- | ---: | ---: | ---: |
| 0.1 | 7 | 0.136 | 65537 |
| 0.01 | 19 | 0.421 | 65537 |
| 0.0001 | 161 | 3.457 | 65537 |
| 1e-08 | 16395 | 333.532 | 65537 |
| 1e-15 | 65537 | 1404.000 | 65537 |

Depth 16 permits 2^16 leaf chords and 65537 endpoint samples per cubic.
A full tree visits 131071 nodes; two interior control-distance calls per node
explain the measured 262142 `_chord_distance` calls in the deepest profile.
The tiny tolerance saturates this cap; accuracy is best effort at exhaustion,
not an absolute error guarantee. Multiple segments multiply this bound, with
shared endpoints removed. A caller's tolerance and segment count are therefore
important performance controls already exposed by the API.

## SH scaling

| Samples | Bands / coefficients | Median ms | ns / sample / coefficient |
| --- | --- | ---: | ---: |
| 64 | 1 / 1 | 0.068 | 1069.6 |
| 64 | 3 / 9 | 0.255 | 442.7 |
| 64 | 5 / 25 | 0.594 | 371.0 |
| 256 | 1 / 1 | 0.255 | 995.9 |
| 256 | 3 / 9 | 0.937 | 406.7 |
| 256 | 5 / 25 | 2.314 | 361.6 |
| 1024 | 1 / 1 | 0.947 | 925.2 |
| 1024 | 3 / 9 | 3.770 | 409.1 |
| 1024 | 5 / 25 | 9.012 | 352.0 |

Per-sample/per-coefficient values include all three RGB channels. Fixed validation
and output overhead distort small-N ratios, but larger workloads follow the
expected linear sample and quadratic band growth. Basis evaluation adds a
separate degree-dependent recurrence cost when datasets are prepared; none is
silently attributed to the precomputed projection API.

## Phase 3B–3E priorities and acceptance criteria

These are candidate targets for controlled follow-up measurements, not promised
speedups or authorization to change contracts. Each phase must preserve the full
regression suite and unrelated expected failures, use paired same-environment
runs, report distributions/allocations and investigate regressions beyond noise.
Use at least seven trials and two independent executions; require the improvement
to exceed median variability before calling it a gain. Stop for decisions where
an optimization would change ownership, validation, numerical accuracy or APIs.

| Priority / phase | Candidate work | Proposed acceptance |
| --- | --- | --- |
| 1 — 3B numerical kernels | Reduce inverse scaling/check allocations and repeated traversals; then profile fixed-size matrix loops | Aim for >=20% median inverse improvement on matched typical inputs, with no material extreme-scale regression. Preserve exact singular ZeroDivisionError, no epsilon rejection, signed-infinity rescaling, both inverse identities, row conventions and ctypes state. Keep independent 1e±300 references. |
| 2 — 3C Bezier sampling | Reuse chord state per subdivision node, reduce temporary levels/tuple churn without changing flatness or depth | Aim for >=20% runtime or >=25% traced-peak reduction on difficult paths. Preserve endpoints/order, finite-chord/coincident cases, depth 16, squared tolerance and best-effort cap semantics. Independently measure approximation error, ownership, repeated builders and connected deduplication. |
| 3 — 3D SH workloads | Reduce projection generator/validation overhead, then small analytical L2 tensor overhead | Aim for >=20% median improvement at N=1024/B=3 and 5 or L2 RGB rotation. Preserve compensated accumulation accuracy using independent known answers (1e-12 absolute around zero, relative where meaningful), basis signs, channel independence, energy/composition/inverse checks, unit norm boundary and input preservation. Retain O(NB²) measured scaling. |
| 4 — 3E object overhead | Profile and reduce vector/quaternion/matrix temporary objects and conversion overhead where contracts permit | Aim for >=10% on a clearly named frequent-operation basket without a >5% reproducible unrelated regression. Preserve fresh-return/mutation semantics, angle units, SLERP 1e-12 unit accuracy, nonunit behavior, ctypes synchronization and public signatures. Ray pivot/intersection-state behavior remains fixed. |

A real application's call-frequency profile may reorder these priorities:
Bezier is the largest individual measured cost, but inversion and SH kernels
may dominate frequent workloads. This machine provides hypotheses, not universal
application speed rankings. No optimization implementation starts in Phase 3A.

## Regression verification and clean-source repeat

Verification used pytest 9.1.1. Full suite: **1778 passed, 4 xfailed, 0 unexpected failures, 0 skips**. All four
remaining expected-failure identities match Phase 2F-7. No tests or markers change.
The harness was rerun sequentially from a clean staged-tree source export under
`/tmp/phase3a-clean`, outside the development checkout, without .git/build/cache
files. It runs all current cases without historical objects and produces matching
core/harness fingerprints and deterministic Bezier point counts. The repeat uses
the same interpreter/dependency versions; comparison statistics are recorded in
`phase3a-reproducibility.json`. This is a same-machine clean-source repeat, not
independent hardware validation or a clean-OS claim.


The full clean-source repeat had **6.18% median absolute timing difference**, a
**22.37% 95th-percentile difference**, and a **112.48% maximum**; 37 of 108
cases differed by more than 10%. Representative point counts and inputs are
reproducible, but uniform timing repeatability is not established on this shared
host. The largest outliers are explicitly retained in the machine-readable
comparison. Do not interpret small improvements as signal against this noise.
A follow-up interleaved recheck uses three rounds of seven trials with .05-second
calibration, concentrating on inversion, projection and the anomalous geometry
paths:

```
/workspace/.venvs/pyGameMath/bin/python benchmarks/recheck_variability.py --output audit/phase3a-variability-recheck.json
```

For optimization acceptance, increase trial duration or use controlled CPU/load
conditions until variability is below the proposed improvement target; a nominal
5–10% difference here alone is insufficient evidence.

Interleaved recheck (same source, no optimization):

| Case | Round medians µs | Max/min ratio |
| --- | --- | ---: |
| `vector2_dot` | 0.305, 0.288, 0.276 | 1.105 |
| `matrix3_raw_inverse` | 16.397, 13.761, 15.138 | 1.191 |
| `matrix4_raw_inverse` | 33.753, 36.119, 33.820 | 1.070 |
| `sh_project_n1024_b3` | 3803.606, 3940.139, 3996.111 | 1.051 |
| `sh_rotate_l2_scalar` | 21.019, 18.932, 21.652 | 1.144 |
| `ray_duplicate` | 1.144, 0.950, 1.143 | 1.204 |
| `ray_translate` | 7.733, 6.903, 5.933 | 1.304 |
| `ray_matrix_rotate` | 7.931, 8.500, 7.940 | 1.072 |
| `ray_quaternion_rotate` | 12.882, 11.064, 9.604 | 1.341 |

These rechecks distinguish workload stability from the full-run outliers; they
do not erase the original measurements or control the shared machine. Raw
seven-trial samples for every round remain available in the recheck JSON.

Both CLI scripts executed from clean source exports. The variability recheck
also ran from `/tmp/phase3a-final-clean`, with all nine cases completing all
three rounds. No installed gem implementation, core source or test changed.
