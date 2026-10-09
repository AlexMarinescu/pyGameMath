# Independent core mathematics review

Base: master `ef5a9a36555163e7e6c5b675400f0d0e86087b29` (after PR #47).
Branch: `audit/phase4g1-core-algebra`.

Two confirmed numerical defects remain in supported quaternion powers/logarithms.
They have independent failing references; neither implementation is changed here.
No additional defect was confirmed in ordinary Vector, Matrix2/3/4, conversion,
transformation, camera or viewport operations in the cases examined. This is a
bounded audit, not proof of correctness for every input.

## Review method and coverage

The current `gem/vector.py`, `gem/matrix.py`, `gem/quaternion.py` and
`gem/common.py` were read before consulting previous findings. Review then used
the API reference, historical wiki snapshots, CONVENTIONS, COMPATIBILITY,
PHASE2-DECISIONS, numerical-stability and optimization reports, and existing
regressions to distinguish supported behavior from unresolved contracts.

[Independent tests](../tests/test_core_algebra_audit.py) add 124 cases:
116 passing invariants and eight strict expected failures. References use only
the standard library, with pytest as the existing test dependency:

| Area | Independent checks |
|---|---|
| Vector2/3/4 | Exact dyadic arithmetic/dot; self-aliased in-place addition; independent returning storage; exact equality; cross-product incidence; reflection; non-orthogonal barycentric answer |
| Norms and viewport | 800-digit Decimal norms/directions; smallest subnormal, 1e-300, mixed 1e300/1e-300, overflowing true norm and signed zero; historical whole-vector viewport formula; input preservation |
| Matrix2/3/4 | Leibniz permutation determinants; Fraction Gauss–Jordan inverses; rational products; exact small-integer associativity; row-vector application; self-multiplication; ctypes float32 snapshots |
| Scaled Matrix3/4 inverse | Nonsymmetric matrices at binary exponents -996, -500, 0, 500, 996 (ordinary condition number); independent left/right identity products; exact singularity; failed in-place inversion preserves receiver/export |
| Quaternion algebra | Exact Fraction Hamilton products; nonunit ordinary inverse in both orders; noncommuting composition and reversed row-matrix order; input storage |
| Rotation/conversion | Independent active Rodrigues calculation; random axes and angles; quaternion/matrix vector rotation; length preservation; q/-q equivalence; X/Y/Z exact half-turn branches |
| Interpolation | Analytic same-axis SLERP and squad4 half-angles; negative equivalent endpoint signs; exact half-turn sign ties; one independently derived legacy long-arc SQUAD midpoint |
| Projection | Analytical perspective/orthographic windows with noncommuting modelview; independent windows supplied to unproject; mixed raw/wrapper matrices; non-square offset viewports; negative clip W and unclamped depth |
| Boundaries/helpers | Zero clip W; forbidden general Vector3/Matrix4 promotion; degenerate lookAt; negative-identity log; unsupported operators; angle units and historical keyword names |

Deterministic seeds are `470100+seed` (Vector, seed 0–7), `470200+seed`
(Matrix, 0–7), `470300+seed` (rotation, 0–11) and `470400+seed` (camera, 0–7).
Ordinary random components are small dyadic values, integer matrices, axes in
[-2,2], points in [-4,4], and rotation angles in [-pi,pi]. Camera Z is [-8,-3]
before the explicitly described modelview. Extreme-scale cases are separate.
Fraction inversion reuses the existing independent test helper, not gem inversion.
New tests do not depend exclusively on round trips or before/after equality.

Existing regressions additionally cover refraction/Snell's law, near/far depth,
projective matrices, pivot/shear mappings, ownership, SLERP/SQUAD, conversion
round trips and invalid/nonfinite characterization. They ran unchanged. Bezier,
Legendre, SH, Plane, Ray and advanced algorithms were not newly audited.

## Confirmed findings, ranked for repair

### 4G1-A01 — Subnormal imaginary-axis scaling in power and log

**Severity: medium (P2).** Both functions lose the axis scale for nonzero
subnormal imaginary components near negative identity. Fractional powers can
return nonunit results; logarithms have the wrong imaginary magnitude. Ordinary
rotation inputs are unaffected by this reproduction.

Locations: [quaternion.py:169](../gem/quaternion.py#L169),
[189](../gem/quaternion.py#L189), [201](../gem/quaternion.py#L201) and
[207](../gem/quaternion.py#L207). The original nested hypot and subsequent
component division are the relevant operations, not the principal branch choice.

Minimal reproduction:

```python
import math
from gem.quaternion import Quaternion
t = math.ldexp(1.0, -1074)       # 4.9406564584124654e-324
q = Quaternion([-1.0, t, t, 0.0])
print(q.pow(0.5).data)
print(q.log())
```

| Operation | Actual | Independent expected result (binary64 rounding) |
|---|---|---|
| `pow(0.5)` | `[6.123233995736766e-17, 1.0, 1.0, 0.0]` | approximately `[0, 0.7071067811865476, 0.7071067811865476, 0]` |
| `log()` | `[0, pi, pi, 0]` | `[0, 2.221441469079183, 2.221441469079183, 0]` |

The input norm rounds to 1. Its exact norm differs from unity far below a
binary64 rounding unit, like the existing supported tiny-imaginary tests; this
does not require a new nonunit-domain tolerance. The mathematical imaginary
norm is sqrt(2)*t, but rounding it to the subnormal grid produces t. Dividing
each component by this rounded t produces axis [1,1,0], whose norm is sqrt(2).
The principal angle rounds to pi. The correct axis is [1/sqrt(2),1/sqrt(2),0],
derived independently using 800-digit Decimal calculations. The actual power
has norm 1.4142135623730951 and the log has imaginary norm pi*sqrt(2).
The signed three-component fixture similarly produces axis [1,-1,1] instead
of [1,-1,1]/sqrt(3). Errors are order one, not a last-bit tolerance dispute.

Coverage gap: `test_tiny_imaginary_direction` previously reached 1e-320 but
used a 1:-2:2 axis with a looser relative tolerance. It did not exercise the
minimum-subnormal grid and an irrational axis length. Four new strict failures
cover powers/logarithms and two independent axes. Ownership/type checks pass
before the numerical assertion fails.

Proposed repair: calculate the imaginary **direction** with maximum-component
scaling and hypot before division, independently of the rounded unscaled norm.
The existing stable Vector normalization provides a reference strategy. Retain
principal atan2, exact-zero handling, negative-identity exceptions, signs,
return types and caller storage; no arbitrary axis, epsilon or implicit
whole-quaternion normalization is needed. Verify representable tiny results as
well as negative-identity-adjacent powers/logs.

Compatibility/performance: corrected numbers change only failing range cases;
no API extension is necessary. A small fixed-size scaled direction calculation
adds work to power/log paths. Measure its focused cost during the repair rather
than redesigning norm/interpolation implementations here.

### 4G1-A02 — Large integer powers lose the algebraic rotation phase

**Severity: medium (P2).** Large finite integer exponents produce materially
wrong rotations even for exactly represented unit quaternions with short,
exact power cycles. Unit output length alone does not establish correctness.

Location: [quaternion.py:183–188](../gem/quaternion.py#L183). Multiplying the
rounded principal angle by a large exponent loses its phase; the overflow-only
fallback does not run for this finite product.

Minimal reproduction:

```python
from gem.quaternion import Quaternion
print(Quaternion([0.0, 1.0, 0.0, 0.0]).pow(10**16).data)
```

Actual: `[0.981564736578555, -0.19112997647012858, -0.0, -0.0]`.
Expected: `[1,0,0,0]`, since i^2=-1, i^4=1 and 10^16 is divisible by four.
The exponent itself is exactly representable in binary64. The absolute dot
with the expected quaternion is 0.981564736578555, so this is not q/-q sign
equivalence; the represented rotation is about 22 degrees away from identity.

A second exactly unit control `[.5,.5,.5,.5]` has order six. At 10^16 its
expected result is `[-.5,-.5,-.5,-.5]`; actual is
`[0.014858533001856982,-0.5772865331292115,-0.5772865331292115,-0.5772865331292115]`.
Negative exponents fail too. Exact Fraction Hamilton products over the short
cycles derive all four expected answers without trig or gem powers.

Coverage gap: ordinary integer tests use small exponents. The extreme finite
exponent test checks finiteness and unit length, not an independent orientation.
Four new strict failures cover two exact cycles and both exponent signs.

Proposed repair: investigate an integer-power path using logarithmic Hamilton
squaring, with explicit analysis of roundoff/norm drift for general unit
rotations and preservation of negative powers. Exact cyclic controls are a
necessary oracle. Simply applying fmod after an inaccurate angle product cannot
recover lost phase; reducing with a rounded period is also insufficient as a
general solution. Keep the fractional principal branch and q^0/q^1 rules.

Compatibility/performance: an integer route may change ordinary rounding and
cost from fixed trig work to O(log(abs(exponent))) products. Compare ordinary
results and ownership, and measure only this path in the repair. It must not add
implicit normalization or silently expand the nonunit domain. General huge
fractional-exponent accuracy remains a separate accuracy-policy question; these
exact-cycle failures do not establish a universal precision guarantee.

## Other observations and disposition

These reproductions are **not** added as newly confirmed contract regressions.
The existing API reference explicitly bounds these domains or leaves policy open.

| Classification | Observation and independent mathematical answer | Existing boundary / action |
|---|---|---|
| Known documented limitation | `Quaternion([1e200,0,0,0]).inverse()` returns zeros instead of `[1e-200,0,0,0]` | Quaternion inverse explicitly has no extreme squared-product protection; separate stability review |
| Known documented limitation | `inverse2([[1e200,0],[0,1e200]])` returns zero rows instead of diagonal 1e-200 | Matrix2 inverse explicitly remains unscaled |
| Known documented limitation | Barycentric point `[.25e-100,.25e-100]` in the axis triangle with side 1e-100 raises ZeroDivisionError; weights should be [.5,.25,.25] | Dot/cross/barycentric range stability is explicitly not guaranteed |
| Known documented limitation | Matrix3 with rows `[1,1,1]`, `[1,1+2^-52,1]`, `[1,1,1+2^-52]` passes exact nonsingularity but floating inverse raises ZeroDivisionError | Cofactor cancellation/severe conditioning is explicitly outside accuracy guarantees; determinant exactly 2^-104 |
| Compatibility/accuracy-policy question | Identity `project(Vector(4,[s,s,s,s]),..., [0,0,10,20])` at s=1e-310 returns infinities; direct quotients give [10,20,1] | Nonzero-W reciprocal overflows; near-zero-W/extreme projection policy remains open. Recommend assessing direct component division without inventing a cutoff in a separate review |
| False positive | Matrix2 rotation ignores the supplied pivot | Documented origin-only wrapper behavior; use homogeneous rotate2 for a pivot |
| False positive | Matrix3 `translate3` replaces the final row rather than adding 3D translation | Explicit legacy linear operation; not a general affine 3D translation API |
| False positive | Direct mutation of matrix rows leaves ctypes stale | ctypes is an explicitly documented snapshot; supported in-place operations resynchronize |
| False positive | Opposite endpoint signs at an exact half-turn produce different SLERP midpoints | Documented dot==0 branch tie; new tests preserve the supplied sign rather than impose canonicalization |
| False positive | Vector front is -Z while Quaternion forward is +Z; angle units vary | Established documented conventions, not silently unified |
| Unverified concern | Uniform exception types for malformed shapes, nonfinite values and extra storage | No new validated contract is established; native characterization tests remain unchanged |
| Unverified concern | General enormous fractional powers or ill-conditioned inversion accuracy | No independent exhaustive bound established; no speculative defect claim or new threshold |

Missing wrapper indexing, integer Quaternion scalar operators, unclamped depth,
constructor storage retention and legacy SQUAD nonunit outputs are also
established behavior, not repair targets. Bezier, SH and geometry policies stay
outside this audit. Prior conventions and decisions are not reopened here.

## Verification

Environment: CPython 3.12.14, Linux x86_64, pytest 9.1.1, six 1.17.0.
Python 2.7 and other platforms were not verified. No benchmark campaign ran;
performance recommendations above are algorithmic expectations, not measurements.

| Run | Passed | Expected failures | Unexpected failures | Skips |
|---|---:|---:|---:|---:|
| Unchanged master baseline | 2,260 | 0 | 0 | 0 |
| New audit tests | 116 | 8 | 0 | 0 |
| Complete suite | 2,376 | 8 | 0 | 0 |
| New tests with `--runxfail` | 116 | 0 | 8 deliberately exposed findings | 0 |

The eight expected failures are newly introduced evidence for A01/A02; there
were none on master. No existing test was deleted, skipped, weakened or edited.
`--runxfail` reproduces exactly these eight numerical assertions, with no other
failure. Results and case identities are in
[phase4g1-test-results.json](phase4g1-test-results.json).

Commands from the repository (using the current environment's Python):

```sh
python -m pytest -q
python -m pytest tests/test_core_algebra_audit.py -q -o junit_family=legacy --junitxml=/tmp/phase4g1-focused.xml
python -m pytest tests/test_core_algebra_audit.py --runxfail -q -o junit_family=legacy --junitxml=/tmp/phase4g1-unmarked.xml
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase4g1-full.xml
git diff --exit-code ef5a9a36555163e7e6c5b675400f0d0e86087b29 -- gem setup.py setup.cfg MANIFEST.in pyproject.toml requirements-docs.txt
git diff --check
```

All mathematical source, compatibility shims, dependencies, packaging/release
metadata and existing test files are unchanged against the stated base. Changes
are this report, one new test module and its machine-readable results. No runtime
repair, new public API, release, deployment or Phase 4G-2 work is included.
