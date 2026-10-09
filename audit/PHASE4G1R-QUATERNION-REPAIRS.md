# Quaternion numerical repairs

Base: master `f540db6cc5d48ec9fdd60d0972f26f636c66ab19` (merged PR #48).
PR base includes documentation-workflow updates through
`1ded38f0d36d4a70f264e20ff0324031bb917fe5`; mathematical baseline is unchanged.
Branch: `fix/phase4g1r-quaternion-repairs`.

The eight A01/A02 regressions now pass. Power/log preserve subnormal imaginary
directions, and large integer powers of exactly represented cyclic unit controls
return their algebraic answers. Only `gem/quaternion.py` changes at runtime.
General enormous-exponent accuracy is not claimed.

## A01: direction independent of the rounded subnormal norm

When a nonzero imaginary norm is below the minimum normal binary64 value
(`2.2250738585072014e-308`), components are divided by their largest magnitude
before calculating the axis direction with chained hypot. Direction division
therefore uses an ordinary-scale norm rather than a rounded subnormal norm.
The principal angle still uses the existing atan2 calculation. This boundary
classifies representational range; it is not an epsilon, axis cutoff or new
zero-input rule. No whole-quaternion normalization occurs.

For `t=2^-1074` and `q=[-1,t,t,0]`:

| Operation | Before | After | Independent answer |
|---|---|---|---|
| `pow(.5)` imaginary XY | [1,1] | [0.7071067811865475,0.7071067811865475] | [1/sqrt(2),1/sqrt(2)] |
| `log()` imaginary XY | [pi,pi] | [2.221441469079183,2.221441469079183] | [pi/sqrt(2),pi/sqrt(2)] |

Four original failures cover two axes and both methods. Additional Decimal
references cover signed/zero axis components, subnormal through ordinary tiny
scales, positive/negative fractional powers, and values adjacent to the
minimum-normal boundary. Positive-identity-adjacent subnormal outputs permit
one minimum-subnormal ulp: rounded angle/results cannot have arbitrary relative
precision on that grid. Negative-identity-adjacent finite axis results use tight
normal-scale tolerances. Ordinary power/log calculations retain their previous
operation order and allocate no extra direction list.

## A02: exact cycles with guarded Hamilton squaring

The integer path recognizes unit basis controls and controls with all four
components exactly +/-0.5. Including +/-identity, they form a closed set of
24 exactly represented quaternions. Their orders are 1, 2, 3, 4 or 6, all
dividing 12. Independent Fraction tests verify every one of the 576 products
in this set, exact cycles and the period divisor. Products use dyadic fractions
whose intermediate values remain exactly representable in binary64.

Integer exponents are converted without rounding, reduced modulo the exact
period 12, then evaluated with Hamilton exponentiation by squaring. Negative
exponents use a conjugate copy. Unlike reduction with a rounded trigonometric
period, integer reduction preserves the algebraic phase. The bounded remainder
also avoids squaring through hundreds of exponent bits. Integral float inputs
use their actual represented integer value, not an inferred decimal intention.

At exponent `10**16`, `[0,1,0,0]` now gives `[1,0,0,0]`, and
`[.5,.5,.5,.5]` gives `[-.5,-.5,-.5,-.5]`. At the negative exponent the latter
gives `[-.5,.5,.5,.5]`. Four original failing cases now pass. Additional exact
references cover all 22 non-real controls, signs, odd exponents, integral
floats, the 2^53 boundary, 201-bit integers and `1e308`.

This is deliberately a guarded path. Neighbouring float values are not snapped
to these controls, and there is no new norm tolerance. General unit rotations,
fractional powers and unsupported nonunit legacy paths retain principal-angle
evaluation. Existing operand conversion happens before dispatch, preserving
TypeError/OverflowError behavior rather than broadening exponent protocols.

### Why not square every floating-point unit input?

The benchmark harness also measures an **unguarded investigation prototype**;
it is not part of the runtime. At exponent 10^16:

| Input | Exact represented input norm minus one (Decimal) | Unguarded output norm |
|---|---:|---:|
| [0,1,0,0] | 0 | 1 |
| [.5,.5,.5,.5] | 0 | 1 |
| [cos(.3),sin(.3),0,0] | -4.5475e-17 | 0.5199591499448871 |
| [-sqrt(.5),sqrt(.5),0,0] | +6.8358e-17 | 3.0350351814145453 |

Input approximation and multiplication roundoff amplify under repeated squaring.
Implicit renormalization would change the established contract; leaving the
result unnormalized would break existing finite/unit-result behavior. General
very large exponents therefore keep the previous range-reduction path, whose
phase can still be inaccurate. These repairs do not establish a universal large
fractional/integer power guarantee. Any broader strategy needs separate accuracy
and compatibility analysis.

## Compatibility and verification

Public signatures, [w,x,y,z], principal signs/angles, q^0/q^1, zero errors,
negative-identity parity/fractional rejection, return types and input storage are
preserved. Exact cyclic integer results can replace ordinary trigonometric
residues as well as huge-exponent errors. Nonfinite/invalid paths retain their
existing behavior; no generalized input validator is added. The core stays pure
Python with the existing six dependency. Python 2.7 remains unverified.

Sixteen deterministic general-rotation fixtures (seed `471000+seed`, 0–15)
compare ordinary powers/log bit-for-bit with the pre-repair formula. Separate
nonunit fixtures and nextafter neighbours retain their legacy outputs. These
compatibility checks supplement independent Decimal/Fraction mathematical
references; they are not the only correctness oracles. The benchmark harness
also verifies 288 bit-for-bit comparisons against the actual merged-master
implementation, including scaled nonunit fixtures, outside the timed loop.

| Run | Passed | Expected failures | Unexpected failures / skips |
|---|---:|---:|---:|
| Original eight cases with `--runxfail`, before repair | 0 | 0 | 8 deliberate reproductions / 0 |
| Audit + original power tests + new repair tests | 486 | 0 | 0 / 0 |
| Final full suite | 2,615 | 0 | 0 / 0 |

The final full run took 10.26s. The new module adds 231 cases; the eight existing
regressions lose only their three parameterized strict-defect decorators. No
test is deleted or weakened. Machine-readable counts/case identities are in
[phase4g1r-test-results.json](phase4g1r-test-results.json).

```sh
python -m pytest tests/test_core_algebra_audit.py -k 'subnormal_quaternion or large_integer_quaternion' --runxfail -q
python -m pytest tests/test_quaternion_numerical_repairs.py tests/test_core_algebra_audit.py tests/test_quaternion_powers.py -q -o junit_family=legacy --junitxml=/tmp/phase4g1r-focused.xml
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase4g1r-full.xml
```

## Focused performance measurements

CPython 3.12.14, Linux x86_64, Intel Xeon Platinum 8370C @ 2.80GHz,
pytest 9.1.1 and six 1.17.0. Two separate passes use nine interleaved trials
each on the same host/interpreter. Setup and data generation are excluded;
each side gets 200 warmup calls. GC is disabled only during timings.
The clock is process CPU time, excluding off-CPU scheduling; loop overhead is
included on both sides. Workload order is shuffled deterministically with
`471001+trial`; before/after order alternates. There is no simultaneous pytest
run during the reported comparisons.

Each ordinary trial performs 30,000 operations per side; the large cyclic
workload performs 2,000. General data uses quaternion half-angle .4 and axis
[1/3,-2/3,2/3]; the remaining inputs/exponents appear literally in the harness
and JSON. Raw and wrapper timings are distinguished.

| Workload | Pass 1 median us, before→after | Pass 2 median us, before→after | Paired after/before median, passes 1 / 2 | Combined paired range |
|---|---:|---:|---:|---:|
| `pow_general_8` | 1.026→1.071 | 0.956→1.074 | 1.092 / 1.162 | 0.943–1.908 |
| `pow_general_half` | 0.887→1.053 | 0.958→1.106 | 1.119 / 1.073 | 0.674–1.475 |
| `pow_general_half_raw` | 0.922→1.009 | 0.868→1.054 | 1.272 / 1.235 | 0.850–1.900 |
| `pow_copy_1` | 0.716→0.635 | 0.659→0.729 | 0.945 / 1.061 | 0.313–3.545 |
| `pow_nonunit_2` | 0.883→0.999 | 0.972→1.117 | 1.116 / 1.154 | 0.907–1.542 |
| `pow_cyclic_8` | 0.935→3.710 | 0.951→3.815 | 4.090 / 3.919 | 3.135–5.042 |
| `pow_cyclic_large` | 1.033→2.887 | 1.121→3.343 | 2.772 / 2.841 | 1.991–3.679 |
| `pow_subnormal_half` | 1.489→2.201 | 1.540→2.264 | 1.497 / 1.453 | 1.189–2.023 |
| `log_general` | 0.408→0.430 | 0.444→0.504 | 1.034 / 1.122 | 0.925–1.960 |
| `log_general_raw` | 0.406→0.428 | 0.462→0.450 | 1.055 / 1.011 | 0.519–1.315 |
| `log_subnormal` | 0.973→1.554 | 1.033→1.654 | 1.627 / 1.613 | 1.002–1.852 |

Corrected cyclic powers and scaled subnormal paths cost more than the former
incorrect calculations. Cyclic slowdown is consistently above timing noise in
both passes. Ordinary guard overhead is small in absolute time, but variability
is substantial; the unchanged copy path also fluctuates. These samples do not
justify a precise universal slowdown percentage or a speedup claim. Subnormal
directions allocate one temporary three-component list; exact integer paths
allocate bounded Hamilton product lists; ordinary paths add no temporary list.

Raw samples, median absolute deviations, operation counts, source hash and drift
outputs: [pass 1](phase4g1r-performance.json) and
[pass 2](phase4g1r-performance-repeat.json). Reproduce separately from tests:

```sh
python benchmarks/quaternion_repairs.py --output audit/phase4g1r-performance.json
python benchmarks/quaternion_repairs.py --output audit/phase4g1r-performance-repeat.json
```

The harness loads the immutable baseline source from git into an isolated module
in the same process; it does not install a historical environment or duplicate
the old algorithm in core. Host variation remains, and timings are not release
gates. Mathematical modules other than quaternion, dependencies and packaging
remain unchanged. No release or Phase 4G-2 work is included.
