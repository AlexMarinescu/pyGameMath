# Bezier and Legendre numerical repairs

Base: master `a52c8559ae9bcac6f86d84a9a6ef98b1b2274ac4` (merged PR #52).
Branch: `fix/phase4g2r-curves-legendre`.

The six Bezier underflow and two Legendre overflow regressions from
[Phase 4G-2](PHASE4G2-CURVES-LEGENDRE.md) now pass as ordinary tests.
The complete suite passes **2,984 tests**, with no expected failures,
unexpected failures, errors or skips. Repairs preserve ordinary arithmetic,
public signatures, input ownership and established mathematical conventions.
Only `gem/bezier.py` and `gem/legendre.py` change at runtime.

## Bezier range handling

[Evaluation](../gem/bezier.py#L51) retains the existing Bernstein implementation
for ordinary parameters, endpoints, extrapolation and generic arithmetic.
A narrow fallback handles nonzero float parameters below 2^-511 for quadratics
and 2^-340 for cubics, with finite builtin numeric controls or native matching
Vectors. These bounds describe binary64's normal/subnormal weight range, not an
epsilon, parameter clamp or geometric tolerance. The cubic bound conservatively
covers the point where t³ becomes subnormal. At either bound, 1−t rounds to exactly 1.

The [fallback](../gem/bezier.py#L27) applies small factors to controls before
forming a potentially lost weight. It calculates p2·t·t for the quadratic tail,
p2·(3t)·t for the cubic's quadratic term and p3·t·t·t for the cubic tail.
All multiplicative factors have magnitude less than one in this branch; finite
controls therefore cannot overflow during these products. Keeping 3 with the
first factor avoids a final amplification of an already underflowed intermediate.
The constant and linear terms retain their rounded weights and summation order.

| Reproduction | Before | After / independent Fraction answer |
|---|---:|---:|
| Quadratic (t=1e-200, controls 0,0,1e300) | 0 | 1e-100 |
| Cubic (t=1e-150, controls 0,0,0,1e300) | 0 | 1e-150 |

Scalar outputs remain numeric, and Vector outputs preserve dimensions and own
fresh component storage. This includes empty, Vector1/2/3/4 and larger matching
dimensions. Custom Vector subclasses, mismatched dimensions, Fraction arithmetic,
nonfinite controls and oversized-integer conversions retain their previous
dispatch/exception behavior. Signed tiny parameters are not clamped.

### Alternatives

The benchmark includes weighted de Casteljau, `(1−t)a+tb`, for scalar and Vector
evaluation. It rescues both minimal reproductions but changes 737 of the 1,280
ordinary comparison results at the bit level. It also adds intermediate arrays,
interpolation stages and different custom-arithmetic calls. The reference candidate
does not include the production dispatch/finite guards, so its timings are an
algorithm comparison rather than a drop-in compatibility comparison.

The difference form, `a+t(b−a)`, can overflow b−a even inside [0,1]. For controls
[-1e308,1e308,-1e308] at t=.5, the exact quadratic and weighted de Casteljau return
zero; the difference form returns NaN. A universal replacement is therefore not
selected. The existing path sampler's de Casteljau implementation is unchanged.

The repair fixes premature weight loss, not every summation or cancellation
problem. Near cancellation, error must be assessed against the absolute weighted
input scale. Genuinely unrepresentable tiny terms still round to zero. No universal
extreme-range accuracy bound or new generic-arithmetic contract is introduced.

## Legendre overflow handling

[Ordinary extrapolation](../gem/legendre.py#L48) retains the original recurrence
for finite steps. If the result overflows while x and both previous values remain
finite, it decomposes the products with frexp, aligns exponents, subtracts/divides
at the scaled magnitude, then rescales with ldexp. A genuinely unrepresentable
result remains signed infinity rather than acquiring a new OverflowError.

| Independent polynomial answer | Before | After |
|---|---:|---:|
| P2(±1e154) = (3x²−1)/2 | inf | 1.5e308 |
| P3(4e102) = (5x³−3x)/2 | inf | 1.6e308 |
| P4(8e76) = (35x⁴−30x²+3)/8 | inf | 1.792e308 |

Builtin integer arguments such as ±10^154 are also repaired. Fraction polynomial
coefficients provide independent references for these values and finite/infinite
output boundaries. `run()` still preserves scratch fields; `calculatePML()`
retains its explicit seed/PM1/PML updates and None return.

The original loop remains available without per-step checks for associated inputs,
ordinary x in [-1,1], and a conservative ordinary-extrapolation region: integer
degree ≤64 and |x|≤2. For m=0, bounded |x| and state magnitude S, each recurrence
step has absolute magnitude below 8S, including generous rounding headroom.
Starting with seeds at most two, even undivided products remain below
2^(3l+8), at most 2^200. This bound excludes binary64 overflow without guessing
a numerical tolerance. Independent tests straddle degrees 63/64/65 and x=2.

Distributing division before every multiplication is another possible repair,
but the investigated form changes 221 ordinary-polynomial comparison results.
The selected retry preserves existing finite-step rounding instead. Condon–Shortley
phase, unnormalized associated values, associated-domain errors, nonfinite behavior
and invalid-order exceptions are unchanged. High-order overflow/underflow,
severe cancellation and propagation after an unrepresentable preceding state
retain their documented limitations; this is not general extreme-order stabilization.

## Verification and compatibility

The original eight failures were reproduced with `--runxfail` before changes.
Their markers were removed only after all eight passed with xfail disabled.
No numerical assertion or existing tolerance was weakened. The historical audit
report remains intact; its recorded counts describe its original base.

[New regressions](../tests/test_curve_numerical_repairs.py) add 91 cases, including
many deterministic inputs within each parameterized case:

- Exact Fraction Bernstein references for every term, negative parameters,
  minimum subnormals, representable normal/subnormal outputs, true underflow,
  range-guard boundaries and cancellation budgets.
- Scalar and Vector dimension/ownership checks, endpoints, extrapolation,
  Fraction/custom dispatch, mismatched dimensions and historical exceptions.
- Exact ordinary-polynomial references at extreme arguments and both sides of
  the finite-output boundary, integer inputs and repeated scratch-helper use.
- Associated signs, state preservation, conservative-bound boundaries and
  unchanged invalid/nonfinite characterization.

Matched comparisons against the actual master source verify bit-identical outputs
for **1,280 ordinary Bezier calls**, **1,007 Legendre calls** and **1,140 SH basis
values**. SH consumers remain on the unchanged in-domain recurrence; no SH
algorithm is rewritten. Existing addition-theorem, projection, analytical rotation,
HDR-reference and generated-asset tests pass in the full suite.

| Verification | Passed | Xfailed | Failed/errors | Skips |
|---|---:|---:|---:|---:|
| Unchanged master | 2,885 | 8 | 0 | 0 |
| Original reproduction-only run, xfail disabled | 0 | 0 | 8 intentional failures | 0 |
| Corrected original cases | 8 | 0 | 0 | 0 |
| Focused curves/Legendre suite | 637 | 0 | 0 | 0 |
| Complete suite | 2,984 | 0 | 0 | 0 |
| Isolated installed-wheel smoke | 10 | 0 | 0 | 0 |

Wheel and sdist builds used a fresh staging copy, with unchanged packaging
metadata. Both archives contain byte-identical corrected core sources. A separate
venv, containing only the installed gem wheel and established six dependency,
passes all original range cases plus integer extrapolation with `python -I`,
outside the source tree. `pip check` reports no broken requirements.
The [smoke script](../benchmarks/curve_wheel_smoke.py) asserts that imports come
from that environment, not the checkout, and checks compatibility reexports.

Counts, former-failure identities, artifact/source hashes and protected-file checks
are in [phase4g2r-test-results.json](phase4g2r-test-results.json). Public API docs
and COMPATIBILITY record the range correction. No release metadata, dependencies,
other mathematical modules or experimental implementations change.

## Performance measurements

Environment: CPython 3.12.14, Linux x86_64, Intel Xeon Platinum 8370C at advertised
2.80GHz, six 1.17.0, pytest 9.1.1. Other interpreters/platforms were not verified.
[Harness](../benchmarks/curve_repairs.py) loads immutable master source from git;
setup is excluded. Both final runs use nine interleaved trials, 30,000 calls per
operation/trial, 200 warmup calls, process_time_ns, disabled GC and deterministic
workload shuffling (472101+trial). Before/after/candidate order reverses on alternating
trials. Loop overhead is not subtracted. Inputs and operation counts are recorded.

Selected timings below are first-run medians ± median absolute deviation, in
microseconds. Ratios are median **paired after/before** for run one / repeat;
they need not equal the ratio of separate timing medians.

| Workload | Before µs ± MAD | After µs ± MAD | Paired ratio, runs 1 / 2 |
|---|---:|---:|---:|
| Quadratic scalar, ordinary | 0.333 ± 0.057 | 0.383 ± 0.032 | 1.253 / 1.281 |
| Cubic scalar, ordinary | 0.373 ± 0.014 | 0.445 ± 0.015 | 1.162 / 1.282 |
| Quadratic Vector3, ordinary | 1.059 ± 0.166 | 1.102 ± 0.093 | 1.066 / 1.092 |
| Cubic Vector3, ordinary | 1.245 ± 0.107 | 1.298 ± 0.138 | 1.013 / 1.107 |
| P2(.37), constructor+run | 0.710 ± 0.041 | 0.755 ± 0.014 | 1.103 / 1.137 |
| P12(.37), constructor+run | 1.959 ± 0.098 | 1.990 ± 0.074 | 1.027 / 1.035 |
| P12^5(.37), constructor+run | 2.024 ± 0.267 | 2.494 ± 0.735 | 1.000 / 1.014 |
| P12(2), constructor+run | 2.250 ± 0.446 | 2.073 ± 0.164 | 1.067 / 1.068 |
| Quadratic scalar, tiny parameter | 0.298 ± 0.013 | 1.578 ± 0.077 | 5.496 / 5.296 |
| Cubic Vector3, tiny parameter | 1.203 ± 0.046 | 5.657 ± 0.450 | 4.487 / 4.191 |
| P2(1e154), constructor+run | 0.687 ± 0.024 | 1.762 ± 0.105 | 2.528 / 2.552 |
| P2(1e154), bound run | 0.419 ± 0.026 | 1.584 ± 0.169 | 3.697 / 3.563 |

Ordinary scalar Bezier and low-degree ordinary Legendre guards have measurable
cost. Vector differences are smaller and more affected by trial noise; unchanged
associated paths show no consistent material algorithmic slowdown. Tiny-parameter
Bezier validation/temporary lists and overflow-retry exponent work cost several
times the old failing path. These comparisons pay for corrected answers, not
equivalent historical accuracy. No speedup or universal latency guarantee is claimed.
Raw samples, all 29 workloads, candidate timings and variability are retained in
[run one](phase4g2r-performance.json) and [repeat](phase4g2r-performance-repeat.json).

## Reproduction commands

From the checkout using its audit Python:

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase4g2r-full.xml
python -m pytest tests/test_curve_numerical_repairs.py tests/test_curves_legendre_audit.py tests/test_bezier.py tests/test_bezier_sampling.py tests/test_legendre.py -q
python benchmarks/curve_repairs.py --trials 9 --iterations 30000 --output /tmp/phase4g2r-performance.json
python benchmarks/curve_repairs.py --trials 9 --iterations 30000 --output /tmp/phase4g2r-performance-repeat.json
python -m pytest tests/test_core_packaging.py -q
git diff --check
git diff --exit-code a52c8559ae9bcac6f86d84a9a6ef98b1b2274ac4 -- setup.py setup.cfg MANIFEST.in requirements-audit.txt requirements-docs.txt
```

For an explicit wheel smoke, build `python setup.py sdist bdist_wheel` from a
clean staging copy with the existing build tooling, install the wheel and six
into a new venv, then execute:

```sh
cd /tmp
/tmp/phase4g2r-wheel-env/bin/python -I /workspace/pyGameMath/benchmarks/curve_wheel_smoke.py
/tmp/phase4g2r-wheel-env/bin/python -m pip check
```

The verification used `/tmp/phase4g2r-package-final/source` for the clean build;
artifact paths/hashes and the build interpreter are recorded in the JSON.
The feature scope remains frozen through gem 1.0. No new spline, differentiation,
integration or interpolation API is added, and Phase 4G-3 is not started.
