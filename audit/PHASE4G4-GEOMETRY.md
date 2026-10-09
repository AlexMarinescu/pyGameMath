# Phase 4G-4 — Independent ray, plane and geometry audit

Two reproducible numerical defects remain in the implemented geometry: equal-index
grazing refraction loses a nonzero normal component, and Newell polygon normals
can cancel after translation despite nondegenerate represented vertices. Both are
low-priority numerical repairs (P3), not requests for new geometry APIs.

Base: `cdad68f0e7525fdfdc67b0a3aaea3b8112d003c6`, latest master after PR #55.
Branch: `audit/phase4g4-ray-plane-geometry`. Production code is unchanged.

## Inventory and contract evidence

Reviewed source: `gem/ray.py`, `gem/plane.py`, relevant Vector/Matrix kernels and
quaternion vector rotation. Contract evidence: current Ray/Plane/Vector/Matrix API
references, architecture inventory/conventions, geometry and numerical tutorials,
`audit/CONVENTIONS.md`, `COMPATIBILITY.md`, `PHASE2-DECISIONS.md`, Phase 1/1B findings,
Phase 2D-1/2D-2/2D-4/2F-1/2F-2 reports and the Phase 4G-1 core review. Existing
plane/ray, angle/refraction, numerical, transformation and core-algebra regressions
were inspected. Archived Plane/Ray wiki pages contain only “Coming soon”; later
recorded project contracts establish their representation and ownership.

| Implemented supported surface | Established contract | Boundary |
| --- | --- | --- |
| `Ray(startVector,dirVector)` | Retain Vector3 references; normalize caller direction in place; distance is original direction magnitude | Nonzero direction required; replaces direction's component list |
| `Ray.duplicate()` | Copy every stored Vector/list independently, exact distance/state, no constructor rerun | Does not preserve cross-field alias identity |
| `roateUsingMatrix`, `rotateUsingQuaternion` | Matrix3/unit Hamilton quaternion rotation about coordinate origin; replace Vector3 start/dir; normalize dir; None return | Historical spelling retained; no implicit quaternion normalization |
| `Ray.translate(Matrix4)` | Pure translation with local position w=1/direction w=0, Vector3 output, None return | General Matrix*Vector promotion, scale/shear/projective ray semantics absent |
| Ray state/output | `.start`, `.dir`, `.distance`, `.end`; `output()` prints | `.end` is historical hit placeholder/state, not a computed endpoint; transforms leave it untouched |
| `Plane.fromCoeffs`, `fromPoints` | Scalars in a*x+b*y+c*z+d=0; normal=[a,b,c] at coefficient scale; points use unit cross normal and d=-n·a | Direct mutable field edits can desynchronize normal |
| Plane `clone`, `flip`, `normalize` and in-place forms | Fresh returning objects/storage; in-place forms return self; normalize all four coefficients together | Zero normal/degenerate construction raises ZeroDivisionError |
| Plane `dot(Vector4)`, `point_location(plane,point)` | Supplied W participates in dot; sign classification is exact, uses explicit plane argument | Signed distance only for unit normal and w=1; no epsilon |
| `bestFitNormal`, `bestFitD` | Wrapped Newell unit normal; signed D=mean(n·p), construct d=-D; no receiver/input mutation | Nonplanar approximation, not least squares; repeated closing vertex contributes again to D |
| Raw Plane `flip`, `normalize` | Five-element coefficient/normal list → fresh list; coefficient list → four-element tuple | Native malformed/extreme/nonfinite policies not generalized |
| `Vector.reflect`, `refract` | Matching dimensions, unit normal; refraction also unit incident, IOR=n1/n2 and opposing normal | Fresh output; TIR zero-vector sentinel; no implicit normalization/flipping |
| `Vector.barycentric(a,b,c)` | Fresh [u,v,w], negative exterior weights allowed, orthogonal projection for off-plane 3D points | Gram arithmetic; zero denominator raises; extreme stability unpromised |
| `toAngle`, `lperp`, `rperp` | Raw 2D sequences; atan2 radians; fresh ±90° perpendicular Vectors | Wrappers have no indexing protocol |
| Vector direction sign helpers | Positive/negative dot sign, not collinearity | Zero dot is neither direction; unsupported operand NotImplemented |
| Existing geometry Matrix/Vector transforms and `lookAt` | Row-vector products, row-major storage, final-row affine translation; local transform-helper promotion; negative-Z view | Matching dimensions for general multiplication; explicit W retained; no helper perspective divide |

There is **no** Ray evaluation method, ray/plane intersection method, plane
projection/distance-query method, ray-triangle routine, AABB, capsule, collision
system or spatial tree. Mathematical ray evaluation is composed locally as
start+dir*t; t is distance for unit dir. Stored `.distance` does not enforce a
parameter interval. The tutorial-local plane query solves
`t=-(n·origin+d)/(n·dir)`, distinguishing exactly parallel/coplanar lines and negative
parameters. New tests exercise that local equation through existing operations;
they do not introduce or claim an intersection API.

Private `_require_nonzero` and quaternion component kernels are implementation
details, not newly supported public geometry functions. Automated source inventory
of the actual Ray/Plane declarations is included in the verification JSON.

## Independent validation and coverage

New tests use exact Fraction signed triangle areas, translated triangle-fan area
vectors, scalar affine equations and 110-digit Decimal unit directions. Rodrigues
rotation references use standard-library trig and independently normalized axes;
gem cross/normalization/quaternion matrices do not supply expected answers.

Deterministic datasets include seeds 440001–440020 (integer graph planes),
440100–440119 (nondegenerate signed-area triangles), 440200–440214 (off-plane
projection coordinates) and 440300–440319 (Snell frames). Integer points are small,
triangle signed area is at least 8, dyadic weights range from -1 to 2, and embedded
3D/4D coordinates are affine functions of XY. Snell ratios range .5–2, with angles
kept away from the critical boundary in passing randomized tests. Adversarial
finding tests are separate and deterministic.

The 227 new passing cases cover:

- Exact graph-plane incidence and signed distance, homogeneous W=-2/0/1/3,
  positive/negative coefficient scale, point-side sign, clone/flip and normalization.
- Planar/nonplanar polygon area, cyclic wrapping, winding reversal, repeated
  closure weighting and input/receiver preservation.
- Barycentric reconstruction/orientation in Vector2/3/4, exterior weights,
  independent repeated result storage and off-plane 3D projection behavior.
- Fraction reflection answers, involution, lengths, Vector2/3/4 and zero incident;
  independently constructed oblique refraction/Snell frames and TIR.
- Axis and arbitrary Rodrigues ray rotations, negative/zero/positive parameters,
  exact distance/state preservation, caller reference/list ownership and duplication.
- Noncommuting literal rotation/translation, transformed plane incidence and line
  equation, exact parallel/coplanar/behind cases and structured nearly parallel rays.
- Finite normalization scales 1e-300–1e300; three-point cross products at coordinate
  scales 1e-150–1e150 where intermediates remain representable; non-axis lookAt and
  structured near-parallel up vectors through 1e-300.
- Signed-zero exact degeneracy errors, intentional shared-Vector constructor
  aliasing, raw 2D perpendicular/angle helpers and repeatability.

Existing tests independently cover projection/unprojection, pivot/shear mappings,
Matrix ctypes snapshots and no general homogeneous promotion. No existing test,
tolerance, marker or configuration was edited. The 14 new strict cases below
preserve independently confirmed failures for a separate repair phase.

## Severity-ranked findings

| Priority | ID | Defect | Strict cases | Practical consequence |
| --- | --- | --- | ---: | --- |
| P3 | 4G4-A02 | Translated Newell normal cancellation | 8 | Valid polygon rejected or given the wrong orientation when global origin dominates local edges |
| P3 | 4G4-A01 | Equal-index grazing refraction component loss | 6 | Identical media unexpectedly remove a small normal component |

### 4G4-A01 — Refraction loses grazing incidence for identical media

Source: [vector.py:65](../gem/vector.py#L65), through lines 70–71.

```python
from gem.vector import Vector, refract
incident = Vector(3, [1.0, -1e-9, 0.0])
normal = Vector(3, [0.0, 1.0, 0.0])
result = refract(1.0, incident, normal)
# actual [1.0, 0.0, 0.0]; expected [1.0, -1e-9, 0.0]
```

Both supplied lengths round to 1 and the normal opposes incidence. For equal
indices, Snell's law gives the unchanged direction. Algebraically k=d² and
`d+sqrt(k)=0` for d≤0, even without needing to infer sine from the incident norm.
Here `1-d*d` rounds to 1, so k rounds to zero and the returned vector subtracts
its entire nonzero normal component. At 1e-100 the same loss occurs. At 1e-7 the
normal component is -9.996002811937585e-8, already a relative component error of
about 4e-4. The minimal 1e-9 case loses 100% of that component; overall directional
error is only about 1e-9. This is not a large-angle or length failure.

Contract: established n1/n2 refraction, unit/opposing inputs, fresh matching output,
no new normalization/tolerance. Equal indices have an unambiguous identity law;
this is not an antipodal or TIR convention question. Current tests use ordinary
angles/normal incidence and moderate media ratios, leaving grazing equal-index
identity uncovered. Six strict cases cover two small components in Vector2/3/4.
Inputs and storage are checked before the failing component assertion.

Repair recommendation: assess a range-aware discriminant factorization and
transverse/normal decomposition; an exact equal-index identity path is also a
focused option. For example, avoid losing d² through nested unit subtraction.
Retain fresh output, signs, TIR sentinel and exact branch decisions without an
arbitrary clamping epsilon. Independently test critical boundaries and ratios near
one before adopting algebraic reorderings. A new near-critical tolerance or broad
extreme-IOR policy would require a separate contract decision.

Compatibility/performance: the corrected tiny component changes numerical results
where callers currently receive tangent directions. Algebraic reorderings can
change ordinary last bits/critical-side rounding; added operations or an identity
branch have costs to measure in a repair. No runtime performance claim is made here.

### 4G4-A02 — Newell normal depends incorrectly on a large translation

Source: [plane.py:105](../gem/plane.py#L105), through lines 108–109.

```python
from gem.plane import Plane
from gem.vector import Vector
O = float(2**52)
points = [Vector(3, [O,O,O]), Vector(3, [O+1,O,O]), Vector(3, [O,O+1,O])]
normal = Plane().bestFitNormal(points)
# actual ZeroDivisionError; expected Vector3 [0.0, 0.0, 1.0]
```

Every input coordinate and unit edge difference is exactly represented. An
independent exact triangle fan gives area vector [0,0,1], so this is not a
collinear polygon. Newell's Z sum multiplies unit differences by absolute-coordinate
sums near 2^53; these sums lose their unit increment and contributions cancel to
zero. The exact-zero guard consequently reports a geometric degeneracy that is
absent in the represented inputs. No product overflows or underflows.

A non-axis example with vertices `[O,O,O]`, `[O+1,O,O+1]`, `[O,O+1,O+1]` has exact
area [-1,-1,1] and expected unit normal [-1,-1,1]/sqrt(3); actual is [-1,0,0]. Its
edge Gram eigenvalues are 1 and 3, so the local triangle is well-conditioned.
Absolute-coordinate perturbations at this scale can of course alter the triangle;
this finding uses the exact represented vertices and does not promise precision
below their spacing. A broad arbitrary-dynamic-range geometry guarantee is absent.

Contract: wrapped unit Newell normals for ordered polygon vertices, orientation
set by winding. Translation cannot change the plane orientation. Phase 2D-1 tests
cover wrapping/incidence for small/moderate offsets, not this large-origin/unit-edge
cancellation. Eight strict cases cover axis/oblique planes, reversed winding and
implicit/repeated closure. Positive 2^52 reproduces; analogous local triangles at
2^50, 2^51 and negative 2^52 were checked and pass. These observations do not define
a universal magnitude threshold.

Repair recommendation: consider translating vertices to a local reference before
Newell accumulation, optionally compensated summation where justified. Preserve
open/repeated closure, winding, nonplanar area-vector semantics and exact degenerate
errors. Simply normalizing the damaged zero/wrong vector cannot restore its area.
Do not silently replace Newell with a least-squares fitting algorithm.

Compatibility/performance: valid translated polygons would stop raising or return
the correct orientation. Centering/subtraction can change ordinary rounding and
adds per-vertex work/temporary-storage choices; compare ordinary outputs and profile
in a repair. Public signatures, weighting and input ownership need no change.

## Other observations and disposition

The following are not newly marked strict contract failures. Established stability
guarantees are narrow; source/API numerical references explicitly leave dot/cross/
barycentric and wider extreme geometric arithmetic unprotected or policy open.

| Classification | Executed observation / independent answer | Disposition |
| --- | --- | --- |
| Documented limitation | Axis triangle side 1e-100, point [.25e-100,.25e-100], barycentric raises ZeroDivisionError; exact weights [.5,.25,.25] | Previously recorded in 4G-1; unscaled Gram products are not guaranteed stable |
| Documented limitation | Three-point plane edges 1e-200 underflow their cross to zero and raise; exact plane is z=0 | Stable normalization does not stabilize preceding cross products |
| Documented limitation | `normalize([1e308]*4)` returns four zeros; normalized coefficients would each be 1/sqrt(3) | True normal length exceeds binary64; Plane normalization explicitly has no scaled guarantee for every extreme input |
| Documented arithmetic limitation | Reflect [1e308,0,0] against unit +X returns [-inf,nan,nan]; exact [-1e308,0,0] | Unscaled doubled product overflows; finite-norm guarantees do not cover every algorithm; assess separately from ordinary reflection correctness |
| Accuracy-policy question | `bestFitD` on two points with z=1e308 and unit +Z returns inf; exact mean 1e308 | Broader extreme polygon accumulation policy is unresolved; recommend separate scale-aware mean review |
| Accuracy-policy question | IOR=2^54 at exact normal incidence returns zero; exact direction unchanged | Very large-ratio cancellation; no ordinary-media or extreme-IOR accuracy contract established; no new broad domain promise |
| Compatibility question | Zero/nonzero `.end` cannot establish whether a hit is valid | QD04 remains open; transforms intentionally leave the stored state untouched |
| False positive | Ray `.distance` does not clamp t and `.end` is not start+dir*distance | Stored original magnitude and historical placeholder, not interval/endpoint APIs |
| False positive | Ray constructor normalizes a shared start/direction Vector, moving its aliased origin | Explicit reference retention and mutation; new characterization test preserves this behavior |
| False positive | A repeated closing vertex changes nonplanar mean D | Each supplied vertex contributes once; documented weighting |
| False positive | Nonplanar best-fit plane does not contain every vertex | Area-normal/mean approximation, no least-squares/interpolation promise |
| False positive | Same/opposite direction checks accept noncollinear vectors | Documented dot-sign tests, not geometric parallelism |
| Unsupported behavior | Raw lists for Ray construction, nonunit quaternion rotation, Matrix4 scale/shear/projective Ray transforms, wrapper indexing | Existing Vector3/rigid/unit/sequence prerequisites retained |
| Unverified concern | Uniform malformed/nonfinite exceptions, arbitrary ill-conditioned predicates, universal extreme rotation precision | No exhaustive accuracy bound or policy established; no speculative findings/tests |

Nearly parallel divisions can amplify input error. Structured exact examples are
covered, but do not establish universal robust intersection predicates. There is no
collision/visibility API to audit and no basis for promoting retired transport.

## Verification and installed package

Environment: CPython 3.12.14, Linux x86_64 (kernel 6.18.44, glibc 2.41), pytest 9.1.1,
six 1.17.0. Only this interpreter/platform was executed; Python 2.7 support was not
verified. Source and installed reproductions record their actual module paths and
SHA-256 values in [phase4g4-verification.json](phase4g4-verification.json).

| Run | Passed | Failed/errors | Strict xfail | XPASS | Skipped |
| --- | ---: | ---: | ---: | ---: | ---: |
| Unchanged master baseline | 3,313 | 0 | 0 | 0 | 0 |
| New audit cases (from full suite) | 227 | 0 | 14 | 0 | 0 |
| Full suite | 3,540 | 0 | 14 | 0 | 0 |
| Audit with `--runxfail` | 227 | 14 deliberately exposed | 0 | 0 | 0 |

Full pytest elapsed time: 15.71 s; baseline 16.87 s. These are verification runtimes,
not a controlled performance comparison. The 14 xfail identities exactly match
those 14 exposed failures. Every baseline test identity remains present and passes.
JUnit xunit1 retains fixture properties without the xunit2 warning from prior runs;
configuration and assertions are unchanged.

Executed commands from the repository root:

```sh
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -o junit_family=xunit1 --junitxml=/tmp/phase4g4-baseline.xml
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_geometry_audit.py --runxfail -q -o junit_family=xunit1 --junitxml=/tmp/phase4g4-exposed.xml
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -o junit_family=xunit1 --junitxml=/tmp/phase4g4-full.xml
PYTHONPATH=. /workspace/.venvs/pyGameMath/bin/python audit/geometry_reproductions.py --output /tmp/phase4g4-reproductions.json
```

The installed-wheel tests in the full suite pass. Separately, a clean copy of the
unchanged package/build files was staged in `/tmp/phase4g4-package/source`, built
with `setup.py sdist bdist_wheel`, and installed into a fresh venv with only six.
All 17 runtime Python module bytes in wheel and sdist match source; all five audited
installed module hashes match source. Ordinary plane/ray/quaternion/translation
smoke assertions pass outside the checkout, and both findings reproduce identically.

```sh
cd /tmp/phase4g4-package/source
/opt/codex/runtimes/codex-primary-runtime/dependencies/python/bin/python3.12 setup.py sdist bdist_wheel
/workspace/.venvs/pyGameMath/bin/python -m venv /tmp/phase4g4-wheel-env
/tmp/phase4g4-wheel-env/bin/python -m pip install --no-index --no-deps /tmp/phase4d-dependencies/six-1.17.0-py2.py3-none-any.whl /tmp/phase4g4-package/source/dist/gem-0.1.12-py3-none-any.whl
cd /tmp
/tmp/phase4g4-wheel-env/bin/python -I /workspace/pyGameMath/audit/geometry_reproductions.py --installed --output /tmp/phase4g4-wheel.json
```

`-I` prevents source-tree import leakage; the reproduction script does not insert
the checkout into sys.path. The historical package version is unchanged, not a
release assertion. No mandatory dependency was added. Preservation checks verify
all 66 preexisting runtime/test/packaging files in that protected set byte-for-byte,
and the only changes are the four additions below. Public production signatures
and all preexisting tracked files remain unchanged.

## Files added and repair priorities

| New file | Purpose |
| --- | --- |
| `tests/test_geometry_audit.py` | 227 independent passes and 14 confirmed strict failures |
| `audit/geometry_reproductions.py` | Executable minimal findings, ordinary installed smoke and classified numerical observations |
| `audit/phase4g4-verification.json` | Test identities/counts, source/installed findings, inventory, environment and preservation evidence |
| `audit/PHASE4G4-GEOMETRY.md` | Inventory, severity, coverage gaps and repair/compatibility recommendations |

Recommend a separate focused repair for A01/A02, retaining all established geometry
contracts. Extreme mean/IOR policies and ray hit-state design remain separate
reviews. No repair, collision primitive, general robust-predicate engine or new
mathematical API is implemented by this audit. The gem 1.0 feature scope is unchanged.
