# Matrix inversion overhead optimization

Base: master `89cd4b97784d32625de65b8f465cfcbf6e102943`, merged PR #34.
Branch: `perf/phase3b-matrix-inverse`.

## Changes and equivalence

The finite inverse paths still perform an exact represented-binary64 singularity
check, power-of-two row scaling, the existing cofactor inverse, and inverse-column
exponent rescaling. No ordinary-scale shortcut, epsilon or near-singular threshold
is introduced. No Gaussian elimination replaces the implementation.

- One `isfinite` check replaces `isnan` followed by `isinf` for each input entry.
  Nonfinite values still use the legacy unscaled cofactor path. Established invalid
  examples retain TypeError, IndexError and OverflowError.
- For a float ratio n/2^k, row denominator clearing now shifts n by K-k,
  where K is the largest denominator exponent in that row. This is exactly the
  former n*(2^K//2^k) integer arithmetic, including negative numerators and zeros.
  Integer denominators from as_integer_ratio are powers of two by construction.
- The private 4x4 integer determinant expands only its first row: six bottom-row
  2x2 minors produce the four needed cofactors. The previous public det4 computed
  all cofactor rows even though the zero check only used the first. Integer
  arithmetic makes the reduction exact. Public determinant functions are unchanged.
- The 3x3 floating kernel binds nine entries locally, reuses the same first minor
  and constructs result rows directly. Floating multiply/subtract/add/divide order
  is retained; signed-zero and nonfinite outputs are independently checked.
- The 4x4 kernel combines adjugate transpose and elementwise division into one
  fresh nested result, avoiding two general-purpose helper traversals/allocations.
  Cofactor expressions and determinant arithmetic are unchanged.
- Row maxima use map(abs) rather than Python generator expressions. Exponents,
  ldexp scaling/rescaling and signed-infinity overflow handling are unchanged.

Numerator/denominator conversion and scaled lists remain deliberately present.
The exact zero check is not removed. Matrix wrapper construction and float32
ctypes synchronization remain on the public path; no lazy exports or ownership
changes are introduced. Matrix2 inversion and unrelated algorithms are unchanged.

## Correctness and compatibility

178 independent cases add deterministic diagonally dominant nonsymmetric
matrices at 1e-300, 1e-150, 1, 1e150 and 1e300; exact Fraction inverses; both
multiplication orders; returning/in-place ownership and ctypes exports; exact
integer determinants with exponents up to 1099; exact linear-dependence checks;
binary near-singular matrices; mixed row exponents; short/nonnumeric/huge-integer
inputs; and NaN/positive/negative infinity behavior. Existing tolerances are
retained, including 3e-14 inverse/identity checks. Existing minimum-subnormal and
signed-infinity rescaling regressions remain passing.

The test oracle uses Fraction Gaussian elimination and recursive determinants;
these are independent references, not production algorithms or relaxed tolerances.
The clean wheel/sdist regression adds isolated installed-wheel returning/in-place
inverse checks at ordinary/extreme scales with expected bidiagonal coefficients
and public c_matrix checks. No dependency is added.

Compatibility: no public signatures, return types, row-vector conventions,
storage layout, mutation rules, singular exceptions or numerical policies change.
Severely ill-conditioned/unrepresentable cases remain outside the existing accuracy
guarantee; float32 exports still have their independent range/precision limits.

## Commands and measurement boundaries

```
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider --junitxml=/tmp/phase3b.xml
/workspace/.venvs/pyGameMath/bin/python benchmarks/matrix_inverse.py --output audit/phase3b-benchmarks.json
```

The matched comparison loads the post-2G baseline source in the same interpreter;
it does not use historical Phase 2F-1 latency or an unrelated host as its baseline.
Data preparation is excluded. Three rounds alternate before/after ordering, each
with seven calibrated trials targeting .04 seconds. A no-coverage, shared-host
wall timer and timeit's GC policy match Phase 3A; no affinity/frequency isolation
is imposed. Results retain samples, MAD/extrema, medians, counts and source hashes.

Four ordinary dense matrices per size use a fixed seed (307+size), plus uniform
1e-300/1e300 scaling. One mixed-row bidiagonal input per size has binary row
exponents [-300,300,-200,200]. Timings are normalized by inverse count in the
batch; identical Python driver/list overhead is included. Returning wrappers
are prepared outside timing. In-place measurements include fresh receiver
construction equally for both versions to prevent inverse-of-inverse drift.

Profiles run after timings: 100 batches per case; instrumented durations are
not latency measurements. Single-call tracemalloc peak is per batch, not RSS.
Cumulative times overlap and must not be summed. All raw profiles and timing
samples accompany this report. No Phase 3C work is included.

The comparison was executed twice sequentially without a concurrent full test
run. Six paired-round blocks per case provide repeated observations. The
reproducible descriptive summary is generated with:

```
/workspace/.venvs/pyGameMath/bin/python benchmarks/matrix_inverse.py --output audit/phase3b-repeat.json
/workspace/.venvs/pyGameMath/bin/python benchmarks/summarize_matrix_inverse.py audit/phase3b-performance-summary.json audit/phase3b-benchmarks.json audit/phase3b-repeat.json
/workspace/.venvs/pyGameMath/bin/python benchmarks/verify_inverse_baseline.py audit/phase3b-equivalence.json
```

The summary uses 10000 seeded paired-block resamples for a descriptive
95% percentile interval of median speedup. Shared-host load and independence
of blocks are not guaranteed, so this is conditional evidence from these
measurements, not a universal population confidence interval or latency promise.

## Matched performance results

Environment: CPython 3.12.14, Linux-6.18.44-x86_64-with-glibc2.41,
INTEL(R) XEON(R) PLATINUM 8573C, six 1.17.0. Source hashes identify both implementations.

Medians below combine the six block medians. Speedup is the median of paired
before/after ratios, which need not equal the ratio of separate latency medians.
Values are per inverse, including identical batch driver costs.

| Case | Before µs | After µs | Paired speedup | Six-block range | Descriptive 95% interval |
| --- | ---: | ---: | ---: | ---: | ---: |
| `matrix3_ordinary_raw` | 13.967 | 11.398 | 1.206× | 1.131–1.478 | 1.149–1.394 |
| `matrix3_ordinary_returning` | 16.210 | 14.226 | 1.159× | 0.996–1.430 | 1.065–1.330 |
| `matrix3_ordinary_in_place` | 17.284 | 15.089 | 1.145× | 1.033–1.179 | 1.048–1.168 |
| `matrix3_scale_1e-300_raw` | 17.199 | 12.290 | 1.390× | 1.225–1.929 | 1.295–1.701 |
| `matrix3_scale_1e-300_returning` | 18.999 | 14.209 | 1.299× | 1.075–1.488 | 1.166–1.429 |
| `matrix3_scale_1e-300_in_place` | 20.333 | 16.270 | 1.244× | 1.191–1.339 | 1.205–1.316 |
| `matrix3_scale_1e+300_raw` | 28.032 | 25.825 | 1.077× | 1.035–1.186 | 1.040–1.161 |
| `matrix3_scale_1e+300_returning` | 31.335 | 28.543 | 1.089× | 0.248–1.190 | 0.649–1.168 |
| `matrix3_scale_1e+300_in_place` | 31.459 | 29.568 | 1.048× | 0.903–1.219 | 0.964–1.151 |
| `matrix3_mixed_rows_raw` | 11.476 | 9.195 | 1.266× | 1.172–1.302 | 1.187–1.301 |
| `matrix3_mixed_rows_returning` | 14.036 | 11.532 | 1.205× | 1.152–1.603 | 1.170–1.431 |
| `matrix3_mixed_rows_in_place` | 15.060 | 13.025 | 1.139× | 1.030–1.202 | 1.036–1.197 |
| `matrix4_ordinary_raw` | 33.425 | 21.403 | 1.539× | 1.476–1.637 | 1.477–1.617 |
| `matrix4_ordinary_returning` | 35.269 | 24.579 | 1.444× | 1.339–1.474 | 1.367–1.471 |
| `matrix4_ordinary_in_place` | 37.965 | 26.879 | 1.394× | 1.271–1.881 | 1.323–1.667 |
| `matrix4_scale_1e-300_raw` | 34.615 | 22.043 | 1.572× | 1.354–1.761 | 1.416–1.673 |
| `matrix4_scale_1e-300_returning` | 38.999 | 25.795 | 1.524× | 1.424–1.577 | 1.466–1.556 |
| `matrix4_scale_1e-300_in_place` | 41.769 | 28.255 | 1.457× | 1.141–1.509 | 1.245–1.500 |
| `matrix4_scale_1e+300_raw` | 199.296 | 75.544 | 2.693× | 2.455–3.524 | 2.475–3.119 |
| `matrix4_scale_1e+300_returning` | 194.596 | 77.028 | 2.534× | 2.449–3.810 | 2.458–3.199 |
| `matrix4_scale_1e+300_in_place` | 205.267 | 80.752 | 2.539× | 2.246–2.949 | 2.366–2.754 |
| `matrix4_mixed_rows_raw` | 25.379 | 17.439 | 1.493× | 1.394–1.589 | 1.407–1.558 |
| `matrix4_mixed_rows_returning` | 28.444 | 20.807 | 1.363× | 1.311–1.448 | 1.320–1.441 |
| `matrix4_mixed_rows_in_place` | 30.879 | 23.033 | 1.364× | 0.995–1.573 | 1.167–1.482 |

Ordinary raw Matrix3/Matrix4 gains are about **1.21x/1.54x** across two executions;
all six blocks improved, and their descriptive intervals exclude 1. This is
repeatable same-environment evidence above the within-trial noise in these cases,
not a comparison with the historical Phase 3A host snapshot. Returning Matrix4
also improves in every block, about 1.44x overall. Returning Matrix3 is more
modest (about 1.16x), with one essentially tied block. At 1e300, Matrix4 raw
improves about 2.69x; large exact integers make the removed determinant work
particularly expensive.

Matrix3 large-scale in-place results do **not** establish a speedup: the descriptive
interval includes 1 and block ratios include regressions. Its median roughly
1.05x is noise-sensitive. No blanket claim applies to every wrapper or regime.
All raw samples, MAD, trial counts and per-round values are retained. Phase 3A's
large environmental variation remains a reason to require matched repeats;
small differences alone cannot justify claims. No material extreme-scale
correctness regression was observed in independent references or bitwise checks.

## Profiling and allocations

| Case (batch) | Before traced peak bytes | After traced peak bytes |
| --- | ---: | ---: |
| `matrix3_ordinary_raw` | 2864 | 2248 |
| `matrix4_ordinary_raw` | 6204 | 4320 |
| `matrix4_ordinary_returning` | 7372 | 5728 |
| `matrix4_scale_1e+300_raw` | 20936 | 9280 |
| `matrix3_mixed_rows_raw` | 1264 | 1560 |

- Finite classification drops from two math calls per finite entry to one;
  the guard/generator remains because nonfinite behavior must be preserved.
- The exact-integer row setup still converts every input, materializes ratios,
  clears denominators and calculates an exact zero check. Shift operations
  replace denominator quotient/multiplication; bit_length calls remain visible.
- Matrix4's 1e300 profile attributes much of the old cost to public det4's
  unnecessary integer cofactors. The reduced exact helper removes that work
  without changing singularity classification. No scaled path is bypassed.
- The floating 3x3 kernel avoids repeated nested indexing and zero-fill output;
  Matrix4 avoids separate transpose and division outputs. Ordinary Matrix4
  raw traced peak falls by about 30%, and its returning-wrapper peak by about
  22%. At 1e300 the raw batch peak falls by about 56%.
- Allocation does not improve uniformly: the mixed-row Matrix3 peak rises
  from 1264 to 1560 bytes in the primary trace. Local scalar lifetimes and
  integer/shift temporaries still contribute. Peak traces are limited observations,
  not an RSS guarantee or a count of every native allocation.
- Wrapper and c_matrix construction remain unchanged and visible in profiles;
  these costs reduce the relative public-method improvement. The in-place
  benchmark intentionally includes receiver construction and cannot isolate
  only synchronization cost by subtracting medians.

cProfile cumulative frames overlap, and instrumentation changes costs. The
unprofiled calibrated timing data establish latency; profile rows identify
where work moved. The second run's full profiles are retained as well.

## Verification and stop boundary

Full pytest: **1999 passed, 0 xfailed, 0 unexpected failures, 0 skips**.
Baseline: 1821 passed; 178 added cases explain the difference. No existing test
or tolerance is weakened or removed. The complete run includes clean wheel/sdist
and isolated installed-wheel inverse/export checks. pytest 9.1.1 runs under
CPython 3.12.14 with six 1.17.0. Existing JUnit metadata warnings are not failures.

The separate preservation script verifies **500 bit-identical finite inverses**
at ordinary/extreme scales and AST-equivalent definitions for all 33 unaffected
public/unrelated functions and the Matrix class. This differential evidence
supplements, rather than replaces, independent Fraction/recursive determinant
regressions. Source and result hashes are saved with the artifacts.

Only gem/matrix.py changes in production. New tests, installed-wheel checks,
benchmark/proof tools and audit documentation support that change. No Phase 3C
implementation, merge or package publication is included.
