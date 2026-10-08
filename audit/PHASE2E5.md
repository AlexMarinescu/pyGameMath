# Quaternion API contract review

Base: master `5a106997ac7392cd9d4f05a7f2d8c9f9092d7e03` (PR #21).
Branch: `fix/phase2e5-quaternion-contracts`.

## Q11 evidence

The original `4253839` helper is misspelled `quat_roate_from_axis_angle`.
Commit `b7eca0f` corrects the name; its implementation still constructs a
rotation quaternion and applies it to its own normalized axis using a
Hamilton conjugate sandwich. The old docstring instead describes creating
an arbitrary-axis rotation quaternion. The historical wiki is a placeholder.

For nonzero axis n, rotation about n leaves n unchanged. Consequently the
helper returns approximately `[0,n/|n|]`, with ordinary roundoff, independently
of theta. At +Z and theta=90 degrees it returns `[0,0,0,1]`; the conventional
rotation quaternion is `[sqrt(0.5),0,0,sqrt(0.5)]`. Interpreting the legacy
result as an orientation produces a half-turn even when theta=0.

Executing `b7eca0f` and current source on Vector axis `[0,0,2]` with angles
0,90,180 confirms the same Quaternion values. Historical code normalized
the caller in place and list inputs crashed. Phase 2E-2 fixed list inputs
and preserved axis storage using a temporary normalization; those ownership
and representation decisions remain established.

Repository caller search finds only tests for this helper, including
explicit legacy-value/ownership regressions. No production, launcher or
benchmark caller uses it. No downstream compatibility claim follows from
that search. Changing its output would change the public semantic contract
without changing its Python result type.

`quat_from_axis_angle(axis,theta)` already returns the conventional
`[cos(theta/2), normalized_axis*sin(theta/2)]` rotation Quaternion in degrees,
without mutating its Vector/list input. A second constructor is unnecessary
to obtain the conventional result.

## Compatibility strategy

Retain `quat_rotate_from_axis_angle` and its current signature, numerical
sandwich, accepted representations, ownership and boundary behavior. Describe
it explicitly as a legacy axis-rotation result and point rotation callers
to `quat_from_axis_angle`. No runtime deprecation warning or
redundant new constructor is added.

Q11 now verifies the retained legacy-result contract, with separate
independent constructor regressions. The [public API guide](../docs/QUATERNIONS.md)
shows the distinction and correct vector-rotation sandwich. README removes
unsupported quaternion cross-product/exponential claims. Executable
quaternion code is unchanged, verified by an AST comparison excluding
docstrings. No mathematical algorithm is altered.

## Remaining quaternion review

| Area | Established behavior / deferred question |
| --- | --- |
| Components and composition | [w,x,y,z], Hamilton product, row-vector rotation matrices, q/-q rotation equivalence; no redesign |
| Helper return types | Axis-angle constructor returns Quaternion; X/Y/Z angle helpers return lists; vector rotation returns Vector; q*Vector returns Quaternion; log returns list |
| Ownership | Constructor retains supplied component storage; returning arithmetic creates independent storage, i-methods change only receiver; axis helpers preserve lists/Vectors since Phase 2E-2 |
| Units and forward axes | Arbitrary-axis helpers use degrees, X/Y/Z angle helpers use radians; Quaternion forward is +Z while Vector front is -Z; preserve both |
| Rotation domains | Vector rotation/matrix conversion assume unit rotations, with no automatic normalization; nonunit sandwich scales by squared norm; broader validation remains separate |
| Zero and nonunit arithmetic | Inverse supports ordinary nonzero general quaternions; zero inverse raises ZeroDivisionError; zero normalization and extreme norms remain N01/N02 work, outside this phase |
| Powers and logarithms | Established unit-only principal-angle domain, zero rejection, negative-identity branch and result ownership from Phase 2E-3; no changes |
| Interpolation | Accurate shortest-path SLERP/squad4 with 1e-12 unit-input guarantees on [0,1]; sign-preserving half-turn ties; legacy no-invert/SQUAD approximation and antipodal ambiguity retained |
| Numeric protocols | Float-only scalar multiplication/division wrappers; no new reflected operators, malformed-shape or nonfinite policy |
| Public documentation | Consolidate current APIs and examples; README advertises quaternion Cross Product and Exponential although neither public API exists; correct those claims without adding algorithms |

The Q01-Q10 correctness changes remain settled. Remaining N01/N02 numerical
failures, nonfinite/invalid-input policies and no-invert antipodal ambiguity
prevent a claim of universal quaternion robustness. Control generation and
an exponential API remain outside scope.

## Baseline verification

Full suite: **1,271 passed, 0 failed, 34 expected failures** (1,305 cases).
Exposing the original `test_arbitrary_axis_helper_returns_rotation_quaternion`
with `--runxfail` reproduced its assertion failure: `[0,0,0,1]` differs
from the axis-angle constructor. Q11's incompatible constructor expectation
is replaced with a passing legacy-result regression.

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/contracts-baseline.xml
python -m pytest tests/test_api_edges.py::test_arbitrary_axis_helper_returns_rotation_quaternion --runxfail -q -p no:cacheprovider
```


## Final verification

Full suite: **1,282 passed, 0 failed, 33 expected failures** (1,315 cases).
Ten added cases cover literal standard/mixed-axis answers, signed angles,
zero angle, list/Vector storage, fresh results, the documented vector
rotation example, and the established constructor-list aliasing behavior.
Existing constructor and legacy-helper regressions remain intact.

Only Q11 leaves the expected-failure set; all 33 unrelated identities match
the baseline. Diagnostic --runxfail exposes exactly those 33 failures with
1,282 passes. Remaining cases comprise 29 defects and four contract
questions; quaternion zero normalization and extreme magnitude still fall
under the separate N01/N02 work. [Detailed results](phase2e5-test-results.json)
record every remaining identity. No quaternion Q01-Q11 marker remains.

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/contracts-final.xml
python -m pytest --runxfail -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/contracts-runxfail.xml
```

## Changed files

| File | Change |
| --- | --- |
| `gem/quaternion.py` | Clarify legacy helper docstring; executable code unchanged |
| `tests/test_api_edges.py` | Replace incorrect Q11 constructor expectation |
| `tests/test_quaternion_contracts.py` | Ten independent result/ownership cases |
| `README.rst` | Remove unsupported quaternion claims and link API guide |
| `docs/QUATERNIONS.md` | Public APIs, units, domains, ownership and examples |
| `audit/CONVENTIONS.md` | Record Q11 resolution and quaternion scope |
| `audit/COMPATIBILITY.md` | Compatibility implications |
| `audit/PHASE2-DECISIONS.md` | Close Q11 while retaining numerical questions |
| `audit/PHASE2E5.md` | Historical evidence, review and verification |
| `audit/phase2e5-test-results.json` | Counts and remaining failure identities |
