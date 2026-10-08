# Quaternion interpolation

Q05 restores the legacy three-control nested blend by replacing the
callable t expression with multiplication. Q10 uses accurate shortest-path
spherical interpolation and stable difference/sum norms for its angle.
`squad4` adds conventional four-control SQUAD without changing the legacy
signature. All results preserve caller objects and component storage.

Base: master `4eafc86731aac609a6bc81113207c0e5b9444308` (PR #20).
Branch: `fix/phase2e4-quaternion-interpolation`.

## Mathematical evidence and compatibility

Historical source at `4253839` and `b7eca0f` contains the same three-control
expression; the archived quaternion wiki is a placeholder. The retained
contract uses q0=start, q2=end, and q1=additional blend control, with
slerp_no_invert throughout. Unit length is not guaranteed in its linear
branches, and independent sign flips can change the path.

Conventional squad4 uses q0/q1 endpoints and s0/s1 SQUAD controls, with
accurate shortest-path SLERP. The controls are not ordinary neighbouring
keyframes. SLERP sign correction and exact half-turn ties are preserved;
inputs are not normalized and parameters are not clamped. Public examples,
formulas and accuracy domain are in [conventions](CONVENTIONS.md#quaternion-interpolation).
See [compatibility](COMPATIBILITY.md#phase-2e-4-quaternion-interpolation).

Independent same-axis measurements sample 1,001 parameters per separation.
At one degree, old norm error is 9.519279e-6; near the old 5.125117-degree
cutoff it is 2.500313e-4, with spatial-angle error 2.869875e-6 radians.
Corrected cutoff norm/component errors are 2.22e-16. The old midpoint norm
is sqrt((1+dot)/2). A median seven-by-50,000-call one-degree benchmark
measured 2.610 microseconds baseline and 3.949 corrected (about 51% slower).
The earlier acos-based prototype estimate did not include the final stable
angle computation. Timing is environment-specific, not a performance promise.
[Measurements](phase2e4-interpolation-measurements.json) can be reproduced:

```sh
python audit/benchmark-quaternion-interpolation.py
```

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Baseline: 1,128 passed,
36 expected failures. Q05 reproduces TypeError under --runxfail. The initial
140 additional cases produced 113 failures and 27 passes before correction.
The final 141 cases cover endpoint/identical behavior, independent
known-axis curves, orthogonal great-circle midpoints, control influence,
sign-equivalence and half-turn ties, dense deterministic unit-norm checks,
local continuity, tiny angles, extrapolation and independent result storage.

Full suite: **1,271 passed, 0 failed, 34 expected failures** (1,305 cases).
Only the Q05 defect and Q10 accuracy-question markers are removed. The exact
remaining set matches baseline minus those two cases. Diagnostic --runxfail:
1,271 passed, exactly those 34 failures (29 defects and five contracts).
[Detailed results](phase2e4-test-results.json) list every remaining case.

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/interp-final.xml
python -m pytest --runxfail -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/interp-runxfail.xml
```

## Changed files

| File | Change |
| --- | --- |
| `gem/quaternion.py` | Accurate SLERP, legacy SQUAD fix, separate squad4 |
| `tests/test_quaternion.py` | Convert Q05/Q10 cases |
| `tests/test_quaternion_interpolation.py` | 141 independent regression cases |
| `audit/CONVENTIONS.md` | Interpolation contracts and public examples |
| `audit/COMPATIBILITY.md` | Numerical changes, additive API and cost |
| `audit/PHASE2-DECISIONS.md` | Q05/Q10 contract scope |
| `audit/PHASE2E4.md` | Evidence, verification and file summary |
| `audit/phase2e4-test-results.json` | Counts and remaining failure identities |
| `audit/benchmark-quaternion-interpolation.py` | Reproducible baseline comparison |
| `audit/phase2e4-interpolation-measurements.json` | Accuracy and timing measurements |
