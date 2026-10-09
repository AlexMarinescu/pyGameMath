# Final Vector and viewport contracts

Base: master `c1e3fe745668b205dbd7e2247350fbb6b4f0962c`, merged PR #33.
Branch: `fix/phase2g-final-contracts`.

## Corrections and compatibility

Vector equality checks declared dimensions before reading components. Empty
Vectors compare equal, empty/nonempty and different dimensions compare unequal
in both orders, and exact component comparison is retained. Inequality has
complementary results using the same direct comparison structure; both methods
preserve NotImplemented for other types. No tolerance or storage validation is
introduced. New boolean outcomes replace None, inconsistent prefix comparison
and possible IndexError for these newly resolved domains.

Returning clamp copies value storage before the historical max-then-min steps.
It preserves the method receiver, input values and bound lists, including shared
storage. In-place clamp still changes and returns the receiver, replacing its
list without altering other owners of the original data. Aliased bounds remain
original inputs; the old incidental bound mutation is removed. Signatures,
supported dimensions, comparison arithmetic and reversed-bound validation policy
are unchanged. Callers relying on returning clamp's mutation must explicitly
assign its result or use i_clamp on the intended receiver.

C02 repairs only Vector component access in getViewPort. Historical source
`4253839:gem/common.py` explicitly normalizes the entire Vector, then computes
`(normalized.x+1)*width/2+original.x` and the corresponding Y expression.
The archived Common Functions wiki is a placeholder, and Phase 1B identifies
the formula as decision-sensitive. The retained contract preserves this source
formula, including original XY offsets and Z/W participation in normalization,
rather than replacing it with an inferred NDC mapping. Direct normalization
uses the existing stable Vector algorithm. The exact-zero guard remains,
raising ZeroDivisionError for zero Vectors. Result is a fresh four-element list;
input data are preserved. Supported Vector2/3/4 cases receive literal known
answers; broader unsupported/nonfinite domains receive no new policy.

This helper does not set OpenGL state or implement project/unproject. Public
[contracts and migration examples](../docs/VECTOR_VIEWPORT_CONTRACTS.md) explain
the unusual normalized-coordinate/offset formula separately. The guide ships
in the source distribution. Core imports, pure Python, quaternion/matrix
conventions, ray state and experimental compatibility shims are unchanged.

## Original failures and regression coverage

Before correction:

```
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_vector_common.py -k 'viewport_vector or clamp_preserves_input or empty_equality or equality_dimensions' --runxfail -q -p no:cacheprovider
```

**4 failed, 45 deselected**: different dimensions compare incorrectly, empty
equality returns None, clamp changes caller values, and getViewPort raises
TypeError by subscripting Vector. The four original test bodies remain intact;
only their strict expected-failure markers are removed after the corrections.

Thirty-nine new cases cover exact Vector2/3/4 comparisons, all pairs of
empty/2/3/4 dimensions, both operand orders, complementary inequality,
unsupported operand protocols, fresh clamp output, receiver/input ownership,
shared lists and aliased bounds, literal viewport coordinates (including negative
offsets and zero/negative viewport extents), repeated-call independence and
zero-vector errors. No tests are deleted, skipped or weakened.

The clean wheel/sdist regression now additionally imports these APIs and checks
empty/dimension comparisons, clamp storage and a literal viewport answer under
python -I outside the source tree. It also preserves the existing core/shim
checks and E07 absence checks. Installation uses --no-index --no-deps with the
existing six dependency supplied explicitly; no runtime dependency is added.

## Verification

CPython 3.12.14, pytest 9.1.1, six 1.17.0.

```
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_final_vector_contracts.py tests/test_vector_common.py -q -p no:cacheprovider
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider --junitxml=/tmp/phase2g.xml
```

Focused contract suite: **88 passed**.
Full suite: **1821 passed, 0 xfailed, 0 unexpected failures, 0 skips**.
Baseline was 1778 passed and four xfailed: four ordinary passes plus 39 added
cases account for the increase. The installed-wheel smoke test is included in
the full run. All four former expected-failure identities are verified present
and passing in JUnit output. Existing per-case record_property/JUnit warnings
are unchanged and do not represent mathematical failures. Machine-readable
accounting is in `phase2g-test-results.json`.

## Performance safeguards

```
/workspace/.venvs/pyGameMath/bin/python benchmarks/final_contracts.py audit/phase2g-benchmarks.json
```

The comparison loads the baseline Vector source from the base git commit into
an isolated namespace in the same interpreter. Fifteen Vector2/3/4 equality,
inequality and clamp cases run three rounds with alternating before/after
order, seven trials and .05-second calibration, reusing the Phase 3A measurement
routine. Public result allocation is included. Both clamp variants receive the
same fresh input-list copy per call to prevent old mutation changing later
inputs; that reset copy is included equally. Empty/mixed-domain results are
correctness changes and are not falsely compared as equivalent old operations.

Final raw samples, median timings and paired ratio ranges are in
`phase2g-benchmarks.json`. No unrelated algorithm is optimized. Inequality uses
direct component checks to avoid adding method-dispatch overhead to the
contract repair. Clamp's ownership guarantee necessarily adds a small list
copy. Shared-host variation documented in Phase 3A limits conclusions about
small changes; ratios below that variation are measurements, not speed claims.

| Case | Before median µs | After median µs | Median paired ratio | Paired range |
| --- | ---: | ---: | ---: | ---: |
| `equality_2_equal` | 0.181 | 0.181 | 1.016 | 0.994–1.044 |
| `inequality_2_equal` | 0.175 | 0.180 | 1.043 | 0.982–1.080 |
| `equality_2_different_last` | 0.171 | 0.179 | 1.053 | 1.037–1.073 |
| `inequality_2_different_last` | 0.169 | 0.174 | 1.031 | 1.025–1.051 |
| `clamp_2` | 0.495 | 0.531 | 1.097 | 1.048–1.214 |
| `equality_3_equal` | 0.196 | 0.203 | 1.023 | 0.467–1.045 |
| `inequality_3_equal` | 0.194 | 0.200 | 1.023 | 1.013–1.034 |
| `equality_3_different_last` | 0.195 | 0.201 | 1.023 | 0.933–1.063 |
| `inequality_3_different_last` | 0.192 | 0.213 | 1.125 | 1.045–1.288 |
| `clamp_3` | 0.656 | 0.659 | 1.006 | 0.920–1.068 |
| `equality_4_equal` | 0.248 | 0.241 | 1.012 | 0.915–1.053 |
| `inequality_4_equal` | 0.228 | 0.280 | 1.224 | 1.093–1.396 |
| `equality_4_different_last` | 0.251 | 0.242 | 1.001 | 0.411–1.071 |
| `inequality_4_different_last` | 0.220 | 0.222 | 1.043 | 0.981–1.051 |
| `clamp_4` | 0.652 | 0.767 | 1.177 | 1.028–1.260 |

These measurements include noise and the small dimension-check/list-copy costs.
No speedup is claimed. Compare ratio ranges and trial variability with the
Phase 3A full-repeat median difference of 6.18% and 22.37% 95th-percentile
variation before attributing a small change to code. The final benchmark run
had no concurrent test execution. No extra inequality-method dispatch is introduced; direct component checks
retain the previous structure. Observed median overheads include +22.4% for
equal Vector4 inequality (paired range +9.3% to +39.6%), +12.5% for Vector3
last-component inequality (+4.5% to +28.8%), and +17.7% for Vector4 clamp
(+2.8% to +26.0%). These are recorded potential regressions, not hidden by
aggregate medians. The wide ranges and Phase 3A noise prevent precise cost
attribution; correctness adds a dimension check and independent clamp storage.
Controlled longer-duration profiling can revisit these costs separately.
