# Vector and Quaternion overhead optimization

Base: master `86dcd405da75c10da976ac67c0b02d774f066920`, merged PR #35.
Branch: `perf/phase3c-vector-quaternion`.

## Changes and numerical equivalence

Finite normalization still validates inputs, scales by the largest absolute
component and uses chained `math.hypot`. The scaled nonzero finite array has
maximum absolute component exactly 1: a largest component divides by its own
absolute value to give +/-1, and every other quotient has magnitude at most 1.
Calling magnitude on that array previously repeated finite checks, computed
the same maximum, divided every component by 1 and multiplied the norm by 1.
The optimized path skips only that redundant work, retaining hypot iteration
and final division order. Zero storage is allocated only for the zero fallback.
Nonfinite arithmetic, empty/generic dimensions and exact-zero policies remain.

Quaternion conjugation constructs its fresh result directly, retaining double
negation of the scalar component, including signed zero and numeric coercion.
Native Quaternion/Vector rotation uses the same two Hamilton kernels and
conjugate sandwich without three temporary Quaternion wrappers. Nonunit input
rotation still scales by squared norm; no normalization is added. Subclasses
retain the existing operator/conjugate path.

A private blend helper creates one native Quaternion directly instead of two
scalar-product wrappers plus their addition. Each component keeps the same
multiply/multiply/add order. Float acceptance and subclass dispatch use the
existing operators when the native shortcut does not apply. Unused identity
allocations in LERP and no-invert SLERP are removed. Interpolation weights,
shortest-path sign handling, nearby-angle accuracy, half-turn ties, no-invert
linear threshold and sign-sensitive legacy SQUAD behavior are unchanged.
Standard squad4 continues composing three accurate SLERPs.

No magnitude, arithmetic, dot/cross, equality, clamp, transform, inverse, power,
logarithm, matrix-conversion or public class implementation is changed. There
are no API, ownership, degree/radian, dependency or ctypes contract changes.

## Independent correctness and compatibility

101 new regressions use 800-digit Decimal normalization references, exact
Fraction Hamilton products/inverses, known-axis interpolation angles and an
independent piecewise reference for the legacy SQUAD approximation. They cover
generic dimensions, zeros, subnormals, 1e-300/1e300 and overrange norms, mixed
scales, nonunit rotation/inversion, near and sign-equivalent SLERP, both SQUADs,
fresh/in-place storage, signed-zero conjugation, subclass dispatch, invalid and
nonfinite arithmetic, and existing transform/clamp behavior. Existing antipodal,
matrix/ctypes, equality and numerical regressions retain their original bounds.

The preservation script additionally verifies **2,018 bit-identical outputs**
against the two independently loaded baseline modules. These include 606 raw
normalizations, 200 each of rotation and all five interpolation entry points,
and conjugate/nonfinite examples. Invalid examples retain exception classes.
AST comparison checks all 51 unaffected top-level definitions, including both
public classes. This differential evidence supplements the independent tests.
No rounding difference is intended or observed for these deterministic inputs.

## Reproduction

```sh
python benchmarks/vector_quaternion.py --output audit/phase3c-benchmarks.json
python benchmarks/vector_quaternion.py --output audit/phase3c-repeat.json
python benchmarks/vector_quaternion.py --output audit/phase3c-followup.json --rounds 7 --target-seconds .04 --case vector2_subtract --case vector3_dot --case quaternion_squad4
python benchmarks/summarize_vector_quaternion.py audit/phase3c-performance-summary.json audit/phase3c-benchmarks.json audit/phase3c-repeat.json --follow-up audit/phase3c-followup.json
python benchmarks/verify_vector_quaternion.py audit/phase3c-equivalence.json
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase3c-tests.xml
```

Commands were executed with `/workspace/.venvs/pyGameMath/bin/python`. The paired
harness requires the baseline git object. Benchmark tools remain outside gem.
Core dependencies and packaging remain unchanged.

## Measurement scope and environment

Two consecutive executions each use three alternating before/after rounds,
seven trials per case and a calibrated minimum .02-second timed batch. Each
call warms before calibration; doubling determines counts saved in JSON. Timeit
disables GC during timing. No concurrent tests/benchmark workloads were run.
Inputs/setup are excluded and result allocation is included. In-place cases
include the same fresh receiver/list setup on both sides to avoid input drift.

The 57 cases reuse Phase 3A deterministic Vector2/3/4 and Quaternion data,
plus raw normalization/conjugation/multiplication, raw function SLERP, clamp,
affine/in-place transforms and zero fallbacks. Ordinary components are roughly
0.25-3.75; extreme vectors use 1e300, 1e-300 and mixed scales, quaternions up to
4e300. Axis rotations use 73/-41 degrees and a 73.001-degree nearby control.
Interpolation uses t=.37; added clamp/in-place inputs are [-2,.5,10], with
[0,5] bounds. No RNG is used for timings. Raw SLERP still accepts
Quaternion objects and returns a Quaternion; it is a function-entry measurement,
not a storage-only arithmetic kernel. Allocation and import costs are not inferred
by subtracting unrelated medians. Full per-trial counts/samples and source hashes
are saved in both run files; no historical cross-host timings are compared.

Environment: CPython 3.12.14, Linux-6.18.44-x86_64-with-glibc2.41,
CPU INTEL(R) XEON(R) PLATINUM 8573C, six 1.17.0; pytest 9.1.1.
CPU affinity/frequency and other host tenants are uncontrolled. Descriptive
intervals use 10,000 seeded paired-round bootstrap samples over six blocks.
They are conditional on these measurements; block independence and stationary
load are not assured. They do not establish universal speedups.

## Paired results

Microseconds per named operation. Medians aggregate round medians; speedup is
the median of paired before/after ratios, so it need not equal the ratio of
the displayed aggregate times. MAD is the median within-round relative MAD.

| Operation | Before us | After us | Paired speedup | Descriptive 95% interval | MAD before/after % |
| --- | ---: | ---: | ---: | --- | ---: |
| `vector2_allocate` | 0.210 | 0.212 | 0.979x | 0.951-1.024 | 2.4/2.4 |
| `vector2_add` | 0.512 | 0.504 | 1.012x | 0.978-1.044 | 2.4/2.6 |
| `vector2_subtract` | 0.505 | 0.518 | 0.975x | 0.967-0.981 | 2.1/2.0 |
| `vector2_scalar_multiply` | 0.536 | 0.529 | 1.017x | 0.989-1.040 | 2.1/2.9 |
| `vector2_dot` | 0.233 | 0.236 | 0.982x | 0.974-1.029 | 2.7/3.2 |
| `vector2_magnitude` | 1.022 | 1.035 | 1.004x | 0.977-1.038 | 2.6/2.6 |
| `vector2_normalize` | 2.467 | 1.473 | 1.664x | 1.527-1.724 | 2.1/3.4 |
| `vector2_transform` | 0.583 | 0.577 | 0.997x | 0.957-1.053 | 3.5/2.0 |
| `vector3_allocate` | 0.216 | 0.216 | 0.986x | 0.953-1.030 | 2.6/1.2 |
| `vector3_add` | 0.530 | 0.551 | 0.959x | 0.931-1.012 | 3.5/3.9 |
| `vector3_subtract` | 0.539 | 0.549 | 0.989x | 0.958-1.012 | 2.2/4.0 |
| `vector3_scalar_multiply` | 0.578 | 0.556 | 1.018x | 0.987-1.117 | 3.0/2.5 |
| `vector3_dot` | 0.267 | 0.278 | 0.960x | 0.943-0.976 | 3.5/3.7 |
| `vector3_magnitude` | 1.258 | 1.273 | 1.001x | 0.977-1.028 | 4.7/2.9 |
| `vector3_normalize` | 2.858 | 1.793 | 1.661x | 1.594-1.717 | 2.7/2.8 |
| `vector3_transform` | 0.934 | 0.933 | 1.013x | 0.993-1.045 | 3.0/2.5 |
| `vector4_allocate` | 0.221 | 0.225 | 0.994x | 0.965-1.019 | 1.9/3.4 |
| `vector4_add` | 0.549 | 0.560 | 0.986x | 0.970-1.045 | 1.6/1.9 |
| `vector4_subtract` | 0.571 | 0.554 | 1.015x | 0.979-1.040 | 1.4/1.6 |
| `vector4_scalar_multiply` | 0.578 | 0.574 | 0.995x | 0.974-1.033 | 2.0/3.5 |
| `vector4_dot` | 0.289 | 0.282 | 1.024x | 0.995-1.058 | 2.0/2.2 |
| `vector4_magnitude` | 1.428 | 1.421 | 1.000x | 0.974-1.053 | 2.7/2.4 |
| `vector4_normalize` | 3.296 | 2.059 | 1.621x | 1.583-1.696 | 3.3/2.5 |
| `vector4_transform` | 1.325 | 1.282 | 1.033x | 1.006-1.090 | 2.9/1.8 |
| `vector3_cross` | 0.480 | 0.477 | 0.999x | 0.986-1.050 | 2.4/2.4 |
| `vector3_magnitude_large` | 1.229 | 1.210 | 1.031x | 0.995-1.041 | 3.1/1.8 |
| `vector3_normalize_large` | 2.804 | 1.739 | 1.631x | 1.608-1.687 | 2.6/1.7 |
| `vector3_magnitude_tiny` | 1.214 | 1.214 | 1.000x | 0.988-1.031 | 2.3/2.6 |
| `vector3_normalize_tiny` | 2.830 | 1.705 | 1.644x | 1.567-1.705 | 2.8/2.6 |
| `vector3_magnitude_mixed` | 1.279 | 1.288 | 0.980x | 0.929-1.016 | 2.5/3.8 |
| `vector3_normalize_mixed` | 2.821 | 1.790 | 1.585x | 1.535-1.681 | 2.4/2.4 |
| `quaternion_normalize_large` | 3.863 | 2.839 | 1.455x | 1.318-1.579 | 2.8/1.9 |
| `quaternion_normalize_tiny` | 3.856 | 2.679 | 1.450x | 1.316-1.537 | 2.5/2.4 |
| `quaternion_allocate` | 0.172 | 0.172 | 0.996x | 0.927-1.052 | 2.0/1.8 |
| `quaternion_multiply` | 0.601 | 0.616 | 0.986x | 0.937-1.109 | 1.8/3.0 |
| `quaternion_normalize` | 3.944 | 2.929 | 1.389x | 1.125-1.527 | 2.2/1.7 |
| `quaternion_inverse` | 0.386 | 0.378 | 1.013x | 0.985-4.814 | 2.0/3.2 |
| `quaternion_rotate_vector` | 2.020 | 1.208 | 1.644x | 1.206-1.713 | 3.5/2.5 |
| `quaternion_slerp` | 3.182 | 2.499 | 1.281x | 1.236-1.318 | 2.5/2.8 |
| `quaternion_slerp_near` | 3.327 | 2.528 | 1.324x | 1.268-1.446 | 4.2/2.5 |
| `quaternion_squad` | 5.920 | 3.292 | 1.893x | 1.698-1.927 | 3.0/2.2 |
| `quaternion_squad4` | 9.795 | 8.356 | 1.205x | 0.951-1.308 | 3.7/9.5 |
| `quaternion_to_matrix` | 6.481 | 6.379 | 0.993x | 0.796-1.198 | 2.5/3.8 |
| `quaternion_from_matrix` | 0.746 | 0.738 | 1.014x | 0.919-1.027 | 2.4/3.7 |
| `vector2_raw_normalize` | 2.172 | 1.238 | 1.813x | 1.598-1.976 | 2.7/3.1 |
| `vector3_raw_normalize` | 2.589 | 1.547 | 1.634x | 1.166-1.784 | 2.9/3.0 |
| `vector4_raw_normalize` | 3.022 | 1.828 | 1.700x | 1.647-9.055 | 3.5/2.9 |
| `quaternion_raw_multiply` | 0.415 | 0.410 | 1.023x | 0.997-1.274 | 3.4/2.8 |
| `quaternion_raw_conjugate` | 0.301 | 0.139 | 2.102x | 1.959-2.185 | 4.1/3.3 |
| `quaternion_conjugate` | 0.418 | 0.231 | 1.725x | 1.597-1.924 | 4.9/2.5 |
| `quaternion_raw_slerp` | 3.219 | 2.523 | 1.275x | 1.049-1.331 | 3.0/3.0 |
| `vector3_clamp` | 0.579 | 0.576 | 1.006x | 0.911-1.344 | 2.2/2.2 |
| `vector3_in_place_normalize` | 2.858 | 1.817 | 1.617x | 1.494-3.378 | 2.8/2.6 |
| `vector3_affine_transform` | 1.144 | 1.121 | 0.999x | 0.891-1.039 | 4.0/2.9 |
| `vector3_in_place_transform` | 1.434 | 1.387 | 1.025x | 0.985-1.207 | 3.1/4.5 |
| `vector3_normalize_zero` | 1.069 | 1.022 | 1.019x | 0.982-1.255 | 3.2/2.7 |
| `quaternion_normalize_zero` | 0.524 | 0.521 | 1.030x | 0.989-1.140 | 3.4/2.0 |

Normalization, conjugation, rotation, shortest-path SLERP and legacy SQUAD
show gains exceeding their descriptive intervals in these runs. Standard squad4
has a 1.205x paired median but its initial interval includes 1. Longer follow-up
trials below provide stronger evidence for this workload. Zero normalization has no demonstrated gain. Most unchanged
controls sit near parity. Outlier blocks (for example raw Vector4 normalization)
show substantial host variation; the raw trials/ranges remain available.

## Longer-trial follow-up

The six primary blocks suggest unchanged Vector2 subtraction and Vector3 dot
are 2.5% and 4.0% slower, respectively, with conditional intervals excluding 1.
Their implementations are AST-identical, so these shifts do not establish a
code regression. A separate, recorded seven-round run doubles the batch target
and checks these controls and the inconclusive squad4 result:

```sh
python benchmarks/vector_quaternion.py --output audit/phase3c-followup.json --rounds 7 --target-seconds .04 --case vector2_subtract --case vector3_dot --case quaternion_squad4
```

| Operation | Before us | After us | Paired speedup | Descriptive 95% interval |
| --- | ---: | ---: | ---: | --- |
| `vector2_subtract` | 0.545 | 0.534 | 1.002x | 0.957-1.069 |
| `vector3_dot` | 0.271 | 0.277 | 0.995x | 0.981-1.006 |
| `quaternion_squad4` | 9.708 | 7.611 | 1.249x | 1.209-1.331 |

The unchanged controls return to parity in the longer follow-up. No material
regression is demonstrated. All seven squad4 follow-up blocks are faster,
with a 1.249x paired median. This targeted follow-up is exploratory and does
not replace the original six-block result; both sets remain available. The
same conditional-bootstrap limitations apply.

## Profiling and allocation

cProfile runs separately over 1,000 calls and tracemalloc measures one call.
Instrumented times are not unprofiled latency, cumulative frames overlap, and
traced peak is neither RSS nor a count of every native allocation.

| Workload | Peak before bytes | Peak after bytes | Quaternion constructors before/after per call |
| --- | ---: | ---: | ---: |
| `vector3_normalize` | 632 | 472 | 0/0 |
| `quaternion_rotate_vector` | 712 | 400 | 3/0 |
| `quaternion_slerp` | 976 | 544 | 3/1 |
| `quaternion_squad` | 1256 | 544 | 14/3 |
| `quaternion_squad4` | 1400 | 752 | 9/3 |

Normalization removes its nested magnitude call, one repeated finite/max scan
and redundant division/rescaling. Its returning Vector construction remains. Native
rotation still allocates raw Hamilton product/conjugate lists and one Vector,
but dispatch and Quaternion wrapper allocation disappear. Native interpolation
retains angle/weight calls and reduces wrapper construction. Legacy SQUAD also
loses discarded identity wrappers from its linear/no-invert branches. Matrix
conversion remains dominated by Matrix allocation and two ctypes snapshots;
it is left intact. Unchanged cross/matrix allocation traces vary with allocator
state, illustrating why single peak observations are not allocation guarantees.

## Verification

Full pytest: **2,100 passed, 0 expected failures, 0 unexpected failures, 0 skips**.
The baseline has 1,999 passes; 101 new cases explain the increase. No existing
tests or tolerances are removed/weakened. Clean offline wheel and sdist builds
copy source into a temporary directory; wheel imports execute with `-I` outside
the repository. Additional installed-wheel checks exercise stable/zero
normalization, known quaternion rotation, SLERP/squad4 norms and preservation.
Machine-readable results are in `phase3c-test-results.json`. A fresh local clone of
the baseline with the staged patch also executes all 57 before/after benchmark
workloads and reproduces the preservation JSON exactly.

Production changes are confined to gem/vector.py and gem/quaternion.py. Tests,
matched benchmark/proof tools and audit documentation support them. This phase
does not include Phase 3D implementation, merging or package publication.
