# Independent curves and Legendre review

Base: master `db303e99c75807b68f262a69b4c4a00ce5a54420` (merged PR #51).
Branch: `audit/phase4g2-curves-legendre`.

Two numerical defects are confirmed: premature Bezier weight underflow and
intermediate overflow in low-degree ordinary Legendre extrapolation. Both have
finite, representable mathematical answers and independent failing references.
No production implementation is changed. Ordinary cases, subdivision, builders,
ownership and the selected associated-Legendre cases pass the new checks.
This is a bounded review, not an exhaustive accuracy guarantee.

## Implementations and contracts

The current implementation was read before consulting prior findings. Review
then used the API inventory/reference, curve tutorial, frozen historical source
at `5257291431bb45db0274dc48edf24694ecfe2e2d`, wiki manifest, CONVENTIONS,
COMPATIBILITY, PHASE2-DECISIONS and the Phase 2F-3A/3B/4 and 3D reports.
The historical wiki contains no Bezier or Legendre page; current core contracts
come from the validated migration and subsequent documentation.

| Implementation | Established behavior |
|---|---|
| `gem.bezier.quadraticBezierPoint`, `cubicBezierPoint` | Degree 2/3 Bernstein evaluation; unclamped parameter, including extrapolation; scalar or matching Vector controls; fresh Vector/storage without input mutation. Evaluation is not limited to Vector2/3. Compatible generic arithmetic takes the fallback branch. |
| `BezierPath.setControlPoints`, `getControlPoints`, `calculateBezerPoint` | Caller-owned control list retained and exposed; cubic layout 3k+1; stored curve count; historical evaluation spelling and indexing. Direct caller edits may stale the count. |
| `interpolate`, `samplePoints` | Finite scalars or uniform Vector2/3 sources. Historical endpoint/interior tangent construction; append-only interpolation, rebuilding source thinning, None returns, fewer-than-two no-op. Squared thinning thresholds are heuristics, not spacing guarantees. |
| `getDrawingPoints`, `findDrawingPoints`, `findDrawingPointsAdded` | Ordered adaptive cubic samples; independent storage; standalone endpoints retained, connected boundaries emitted once; insertion retains existing elements. Positive finite squared tolerance; depth 16; capped output is best effort. |
| Private `_coordinates`, `_point`, `_split`, `_chord_distance`, `_flatness`, `_subdivide` | Representation/finite validation, fresh output, de Casteljau restriction, finite-segment distance and bounded stack subdivision. Tested as implementation details, not promoted to public APIs. |
| `gem.legendre.Legendre(l,m,x).run()` | Integer 0≤m≤l, associated x∈[-1,1]; ordinary m=0 also evaluates outside that interval. Unnormalized Condon–Shortley phase. Local recurrence preserves fields and scratch state. |
| `mGreaterThan0`, `calculatePM1`, `calculatePML(i)` | Explicit mutable scratch helpers; deterministic seed initialization, None returns; target i≥m for supported calculations. |
| `gem.quaternion.quat_squad`, `Quaternion.squad`, `squad4` | Existing orientation-spline blends. Legacy three-control no-invert branch is sign-sensitive and need not be unit length; four-control SQUAD uses shortest-path SLERP. No generated intermediate controls. Bounded same-axis checks supplement the prior quaternion audit. |
| Transitional experimental imports | Identical core objects, including the private historical BezierPath alias; no duplicate algorithms. |

There is **no public Bezier derivative/tangent evaluator**, general-degree curve
API, B-spline, Hermite, Catmull–Rom or NURBS implementation. The tutorial's local
tangent example is not a gem entry point. Derivative identities below test the
existing evaluated polynomial and generated controls; they do not add an API.
SH projection, normalization algorithms and SH rotation are not newly audited.
Existing SH tests still run unchanged in the complete suite.

## Independent validation

[New tests](../tests/test_curves_legendre_audit.py) add 278 cases: 270 passing
checks and eight strict expected failures. References use the standard library
and the existing pytest dependency; no numerical dependency is added.

| Area | Independent reference and coverage |
|---|---|
| Evaluation | Exact Fraction Bernstein sums; endpoints, extrapolation, reversal, coincident controls, scalar Fraction fallback and Vector1/2/3/4 ownership. Dyadic non-diagonal affine maps verify affine invariance. |
| Derivatives | Five-point differentiation, exact for polynomials through degree 3 in real arithmetic, compared with the Bernstein polynomial of degree-scaled control differences; endpoint and interior parameters. |
| Subdivision | Left/right restricted controls calculated from independent prefix/suffix Bernstein sums; exact dyadic evaluation of the restricted polynomial, including t=0/1. Public interval insertion covers empty intervals and preserves endpoint object identities. |
| Flatness/sampling | Fraction segment projection with Decimal square root, including 1e-250/1e250 scales, coincidence and collinear overshoot. Dense independent curve points verify ordered samples and chord error for ordinary 2D/3D cases at two tolerances. Existing depth-16 regressions remain unchanged. |
| Builders | Independently calculated noncollinear tangent controls, repeated/coincident sources, rebuilt output, input storage and geometric joins. Historical append behavior remains covered by existing tests. |
| Orientation splines | Analytic same-axis half-angle polynomials for conventional SQUAD and legacy spherical branches about X/Y/Z; endpoints, selected interior parameters, unit output where promised and input ownership. |
| Legendre values | Exact differentiated Rodrigues coefficients, Fraction argument arithmetic and 180-digit Decimal associated factors; selected orders at degrees 0,1,2,3,7,12,20,32,64,128 across seven signed/boundary arguments. Existing tests cover every order through degree 12. |
| Legendre structure | Parity, independent recurrence right-hand sides, scratch helpers, repeatability, low-degree extreme arguments, factored associated seeds near both poles and ordinary extrapolation. |
| Normalization/high degree | Independent integral identity ∫P_l^m(x)²dx=2(l+m)!/[(2l+1)(l−m)!], checked with Simpson quadrature; exact finite Taylor expansion about 1 and endpoints at degrees up to 1,000. |

Randomness is local and deterministic: seeds `472000+seed`, seed 0–3, controls
in [-4,4] on a 1/8 grid. Ordinary parameters include -0.5 through 1.5. Extreme
evaluation cases use separate scales and tiny parameters. No round trip or one
gem algorithm is the sole oracle for another. Of 231 selected Rodrigues
comparisons, the largest relative error for a nonzero answer was approximately
2.02e-14 at (l,m,x)=(32,1,0.875).

At exact cancellation zeros, the uniform-scale tests allow roundoff proportional
to the sum of absolute weighted controls: 16·2^-53 times that scale. A relative-only
zero comparison would misclassify normal summation residuals. High-degree
near-endpoint tests use an exploratory O(l² ulp) budget, not a new public accuracy
contract. Existing test tolerances were not changed.

## Confirmed findings, ranked for repair

### 4G2-A01 — Bezier weights underflow before multiplication by controls

**Severity: medium (P2).** Supported finite scalar/Vector evaluations can return
zero for a nonzero, normal binary64 answer. The two degree-specific cases affect
both the optimized native Vector branch and scalar fallback.

Locations: [bezier.py:19–23](../gem/bezier.py#L19) and
[32–35](../gem/bezier.py#L32), where t²/t³ are materialized before multiplication
by a control. Six strict failures cover scalars, Vector2 and Vector3.

Minimal reproductions:

```python
from gem.bezier import quadraticBezierPoint, cubicBezierPoint
print(quadraticBezierPoint(1e-200, 0., 0., 1e300))
print(cubicBezierPoint(1e-150, 0., 0., 0., 1e300))
```

| Evaluation | Actual | Independent expected, rounded to binary64 |
|---|---:|---:|
| Quadratic | 0.0 | 1e-100 |
| Cubic | 0.0 | 1e-150 |

With only the final control nonzero, the exact polynomials are p2·t² and p3·t³.
The weight underflows to zero even though multiplication with the large control
would produce a representable normal value. Exact Fraction arithmetic uses the
actual represented inputs before final conversion. There is no cancellation;
relative parameter sensitivity is only the degree (2 or 3), so this is avoidable
range loss rather than severe conditioning. Vector channels with alternating
signs lose the same information; input ownership checks pass before the failure.

Relevant contract: [Bezier API](../docs/api/bezier.md) and Phase 2F-3A preserve
scalar/Vector polynomial evaluation without parameter clamping. Sampling's
explicit extreme-arithmetic limitation is a separate contract. This finding
does not assert a universal floating-point accuracy bound for evaluation.

Coverage gap: earlier uniform-scale and ordinary-parameter checks do not combine
a tiny parameter with a large control. A uniform small scale alone also cannot
expose a weight that should be rescued by multiplication with a large control.

Proposed repair: investigate range-safe weighted products or a carefully selected
de Casteljau evaluation path. Merely reversing three multiplications is not a
general solution; cancellation, extrapolation, overflowing control differences,
generic arithmetic dispatch and endpoint behavior require independent checks.
Preserve signatures, dimension behavior and independent output storage. Add
near-boundary tests for representable and genuinely unrepresentable results.

Compatibility/performance: fixes change the erroneous numbers without adding an
API. Reordering arithmetic can change ordinary last-bit results and custom
operator dispatch; a universal de Casteljau replacement adds interpolation work.
Measure scalar and native Vector paths and compare ordinary rounding before
selecting a repair. No repair or performance claim is made here.

### 4G2-A02 — Low-degree Legendre extrapolation overflows a recurrence numerator

**Severity: medium (P2).** Ordinary P2 outside [-1,1] returns infinity while its
true result remains finite. Both signs of the argument fail.

Location: [legendre.py:36–37](../gem/legendre.py#L36), the multiplication in the
numerator before division by `index-self.m`.

```python
from gem.legendre import Legendre
print(Legendre(2, 0, 1e154).run())
print(Legendre(2, 0, -1e154).run())
```

Actual: `inf` in both cases. Expected: `1.5e308`, rounded from the exact
represented-input polynomial P2(x)=(3x²−1)/2. Scratch fields remain untouched.

At l=2,m=0 the numerator is approximately 3e308, which overflows before division
by two. The mathematical result is below binary64's maximum. Relative sensitivity
is approximately two; this is not a high-order or ill-conditioned example.
Exact Rodrigues coefficients and Fraction arithmetic establish the expected
answer independently of the recurrence under audit.

Relevant contract: [Legendre API](../docs/api/legendre.md), Phase 2F-4 and
CONVENTIONS explicitly support ordinary-polynomial extrapolation. The documented
extreme-order limitation does not describe this degree-two failure. No associated
domain extension or new invalid-input behavior is needed to reproduce it.

Coverage gap: ordinary extrapolation previously used small arguments, and existing
extreme-order observations concern large l/m. Neither exercises a finite final
answer with an overflowing intermediate numerator at low degree.

Proposed repair: assess scale-aware recurrence arithmetic or distributing the
division before risky products, with range and rounding analysis across ordinary
and associated inputs. Simply clamping x would change a supported domain.
Keep phase, unnormalized values, scratch semantics and historical exceptions.

Compatibility/performance: finite answers replace erroneous infinity. Altered
operation order can change ordinary rounding; additional range checks/scaling
could cost more in low-degree and SH-consumer workloads. Evaluate those tradeoffs
in a dedicated repair, without changing SH basis mathematics in this audit.

## Documented limitations and compatibility observations

These are not additional strict expected failures or proposed new contracts.

| Classification | Evidence | Disposition |
|---|---|---|
| Documented sampling limit | Existing actual depth-16 regressions emit 65,537 samples while exceeding tiny tolerance. | Preserve best-effort cap; no unconditional approximation guarantee. |
| Documented extreme sampling arithmetic | A straight Vector2 cubic with X controls [-1e308,-1e307,1e307,1e308] emits an unnecessary midpoint because chord subtraction overflows. Output is finite and the requested geometric path is unchanged. | Extreme sampling arithmetic was excluded in Phase 2F-3B; future range review, not a new accuracy-contract failure. |
| Documented high-order range limit | `Legendre(200,200,.2).run()` returns infinity; exact diagonal magnitude is outside binary64. | Do not require an unrepresentable finite answer or add implicit normalization. |
| Documented high-degree accuracy limit | P1000(nextafter(1,0)) is 0.9999999999213982 versus exact-polynomial reference 0.9999999999444333 (absolute error 2.3035e-11). | Record endpoint recurrence roundoff; no arbitrary universal high-degree tolerance or new defect threshold. |
| Intentional compatibility | Explicit control list/getter aliases, stale count after direct edits, append-only interpolation and coincident standalone endpoint duplicates. | Preserve documented ownership/layout behavior. |
| Continuity distinction | Noncollinear interpolation controls join at one point with parallel incoming/outgoing derivatives of unequal lengths. | Geometric G1 where nondegenerate, not a promised parameter-space C1 or constant-speed spline. |
| Intentional spline branch | Legacy SQUAD retains no-invert SLERP's linear approximations and sign-sensitive behavior. | Unit-length/shortest-path promises belong to squad4; do not redesign legacy behavior. |
| Unsupported domain | Associated abs(x)>1 raises ValueError; noninteger iterated order raises TypeError; l<m returns stored PML. Negative degree/order and nonfinite inputs have no unified policy. | Characterize native behavior; do not manufacture standardized exceptions or silently broaden support. |
| Unverified concern | Arbitrarily high-order associated stability, malformed Vector storage, custom arithmetic types, and bit-identical results across interpreters/platforms. | No exhaustive bound or compatibility claim established. |

The private sampler validates the whole control list on each segment, making
multi-segment validation quadratic in segment count. This is a source-level cost
observation, not a measured performance defect. No benchmark campaign or
optimization is included.

## Verification and review scope

Environment: CPython 3.12.14, Linux x86_64, pytest 9.1.1, six 1.17.0. Other
interpreters, Python 2.7 and other platforms were not verified.

| Run | Passed | Xfailed | Unexpected failures | Skips |
|---|---:|---:|---:|---:|
| Unchanged master baseline | 2,615 | 0 | 0 | 0 |
| New audit tests | 270 | 8 | 0 | 0 |
| Complete suite | 2,885 | 8 | 0 | 0 |

The reproduction-only `--runxfail` run fails exactly eight numerical assertions
(six A01, two A02); 270 unrelated new cases are deselected. Exit 1 is intentional
for that diagnostic command. There are no collection errors, skips, retired tests
or unexpected passes. Existing mathematical tests, including isolated-wheel
packaging/import checks, run unchanged. Machine-readable counts, case identities,
reproduction messages, environment and 23 protected-file SHA-256 hashes are in
[phase4g2-test-results.json](phase4g2-test-results.json).

Reproduce from the checkout with its existing audit environment:

```sh
python -m pytest -q
python -m pytest tests/test_curves_legendre_audit.py -q -o junit_family=legacy --junitxml=/tmp/phase4g2-focused.xml
python -m pytest tests/test_curves_legendre_audit.py --runxfail -k 'tiny_parameter_large_control or low_degree_extrapolation_avoids' -q --tb=short -o junit_family=legacy --junitxml=/tmp/phase4g2-reproductions.xml
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/phase4g2-full.xml
git diff --exit-code db303e99c75807b68f262a69b4c4a00ce5a54420 -- gem setup.py setup.cfg pyproject.toml MANIFEST.in requirements-audit.txt requirements-docs.txt
git diff --check
```

Changes are confined to this report, one new independent test module and its
verification JSON. Every runtime/compatibility module, dependency/packaging file
and existing test remains unchanged. Release scope stays frozen through gem 1.0;
the repair recommendations do not add planned algorithms or authorize production
changes. No implementation repair, SH audit, release or subsequent phase is included.
