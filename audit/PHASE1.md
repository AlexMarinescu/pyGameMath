# Phase 1 correctness and modernization audit

Audited all 14 Python files under `gem/` (including both empty package initializers), plus `launcher.py`, packaging metadata, README, and legacy CI. Source baseline: `5257291431bb45db0274dc48edf24694ecfe2e2d`. Work is isolated on `audit/phase1-correctness`. Library implementation, public API, `gem` namespace, runtime dependencies, and default branch remain unchanged. Phase 2 and Phase 3 have not begun.

## Evidence and priority

The test baseline contains **271 cases: 179 pass and 92 reproduce failures**, grouped under **44 finding IDs**. Normal runs mark these known failures as strict expected failures; `--runxfail` runs their assertions normally and exits nonzero. They were executed against the unchanged source before any correction. A zero exit status with xfails does not mean the library is correct.

Testing used CPython 3.12.14, pytest 9.1.1, coverage 7.16.2, and six 1.17.0. Statement coverage is **82.24%**, branch coverage **62.99%**; combined coverage is 78.34%. Failures prevent reaching much of the experimental code. Coverage is evidence of execution, not proof of completeness. No other Python interpreter, graphics driver, or external integration was tested.

High severity means ordinary supported-looking inputs give wrong math or make an operation unusable. Medium means edge-case robustness, stale secondary state, hidden side effects, or accuracy loss. Experimental failures are prioritized below core failures even where their local severity is high. IDs correspond to `@pytest.mark.defect(...)` and [machine-readable results](test-results.json).

### Core high-priority findings

| ID / severity | Source | Reproduced behavior / cause |
| --- | --- | --- |
| V01 / High | [vector.py:253](../gem/vector.py#L253), :263 | Equality and inequality return inside the first iteration: `[1,2,3] == [1,9,3]` is true. Dimensions are ignored; empty equality returns `None`. |
| M01 / High | [matrix.py:252](../gem/matrix.py#L252) | 2×2 inverse overwrites the first row and never fills the second. For `[[1,2],[3,4]]`, returns `[[-1.5,0.5],[0,0]]` instead of `[[-2,1],[1.5,−0.5]]`; 20 seeded nonsingular cases also fail. |
| Q01 / High | [quaternion.py:77](../gem/quaternion.py#L77) | Inverse divides by squared norm but omits conjugation. For `[0,1,0,0]`, `q*q.inverse()` is `[-1,0,0,0]`, not identity. |
| C01 / High | [common.py:79](../gem/common.py#L79), :83 | Named angle conversions are reversed and use 3.14. `radiansToDegrees(pi)` gives about 0.0548 rather than 180. |
| M02 / High | [matrix.py:39](../gem/matrix.py#L39) | Scalar division transposes off-diagonal entries. Important dependency: inverse4 intentionally constructs cofactors then uses this transposing helper to get the adjugate; fixing division alone would break inverse4. |
| M03 / High | [matrix.py:384](../gem/matrix.py#L384) | Only Python 2 division special methods exist; `Matrix / 2.0` raises TypeError on Python 3. Integer scalar support is also inconsistent. |
| V02 / High | [vector.py:55](../gem/vector.py#L55) | Refraction uses eta³ instead of eta², then uses unsupported scalar-left multiplication. At eta=1.5 and sin(theta)=0.6 it falsely returns a zero/TIR result although `k=0.19>0`. Normal incidence raises TypeError. |
| G01 / High | [plane.py:45](../gem/plane.py#L45) | `fromCoeffs(0,2,0,−4)` calls a cross product on scalar differences and raises AttributeError. Normal should derive from `(a,b,c)`. |
| G02 / High | [plane.py:13](../gem/plane.py#L13), :86, :92 | Computes magnitude after normalizing; offset is not scaled. `[0,2,0,−4]` becomes `[0,1,0,−4]` rather than `[0,1,0,−2]`. Both class normalization methods leave `normal` inconsistent. |
| G03 / High | [plane.py:53](../gem/plane.py#L53) | Stores point Vectors in scalar coefficient fields; positive `d` conflicts with coefficient equation. A point used to construct the plane does not satisfy its own plane-dot incidence check. |
| P01 / High | [matrix.py:670](../gem/matrix.py#L670) | `project` subscripts Matrix despite no indexing API. Both nested-list and Matrix arguments fail. Unreachable return also declares size 3 with 4 values and mismatches window depth mapping. |
| P02 / High | [matrix.py:700](../gem/matrix.py#L700) | `unproject` uses projection*modelview, but row vectors require modelview*projection. Noncommuting translation/projection recovers X=−5 instead of X=1. Identity tests alone miss this. |
| M04 / High | [matrix.py:488](../gem/matrix.py#L488) | In-place 4×4 translation with Vector3 builds a 3×3 helper matrix; wrapper indexes past its rows. Returning translation works on the same input. |
| M05 / High | [matrix.py:60](../gem/matrix.py#L60) | 3×3 XY shear writes index `[0][3]`; IndexError. XY4 puts entries in the homogeneous column, unlike other shear helpers; its intended axis mapping needs review. |
| M06 / High | [matrix.py:130](../gem/matrix.py#L130) | Pivot rotation translation formula is wrong: rotating pivot `(2,3)` around itself by 90° sends it to `(4,6)`. Matrix2's rotation happens to discard translation and passes the origin-rotation probe. |
| V03 / High | [vector.py:150](../gem/vector.py#L150) | Transform mixes column-style multiplication with final-row translation and double-counts a diagonal: identity sends `[2,3,4]` to `[2,3,5]`. |
| Q02 / High | [quaternion.py:442](../gem/quaternion.py#L442) | Matrix conversion chooses a largest-component index without updating its value and uses mutually exclusive comparisons rather than a full maximum. Valid half-turn matrices raise ZeroDivisionError/ValueError; tests cover three basis axes and a mixed axis. |
| Q03 / High | [quaternion.py:159](../gem/quaternion.py#L159) | Power multiplies zero-initialized result vector components, never the original components; `q.pow(1)` is wrong and identity power divides by zero. Unit-only assumptions are undocumented. |
| Q04 / High | [quaternion.py:173](../gem/quaternion.py#L173) | Unit quaternion logarithm has scalar component 1 instead of 0; identity fallback aliases input storage. No general-norm logarithm or exponential implementation exists despite README claims. |
| Q05 / High | [quaternion.py:244](../gem/quaternion.py#L244) | SQUAD evaluates `t(1−t)` as a call; even identical rotations raise TypeError. Three-control signature leaves intended spline semantics unresolved. |
| Q06 / High | [quaternion.py:86](../gem/quaternion.py#L86), :134 | Advertised list-axis branches call `.normalize()` on a list; both fail. Vector-axis branches mutate caller-owned axes. |
| Q08 / High | [quaternion.py:324](../gem/quaternion.py#L324) | In-place Quaternion*Vector reads `other.data` rather than `.vector`; AttributeError. Returning form works and returns a Quaternion product. |
| Q09 / High | [quaternion.py:337](../gem/quaternion.py#L337) | Quaternion division has only Python 2 special methods, so Python 3 `/` fails. Scalar multiplication rejects ints while vector multiplication accepts them. |
| Q11 / High, contract-dependent | [quaternion.py:134](../gem/quaternion.py#L134) | Docstring says it creates a rotation quaternion, but it returns the axis rotated about itself as a pure quaternion. For +Z/90° result is `[0,0,0,1]` instead of the rotation quaternion. Test encodes the docstring; intended historical contract must be approved before correction. |
| R02 / High | [ray.py:25](../gem/ray.py#L25) | Quaternion rotation stores Quaternion products in `start` and `dir`, leaving them without the ray's expected `.vector` field. Its normalization call succeeds on the wrong type; product alone also omits the conjugate sandwich. |
| R03 / High | [ray.py:31](../gem/ray.py#L31) | Translation never moves origin; multiplies direction by a 4×4 matrix without homogeneous promotion and raises IndexError. A translation must leave direction unchanged. |
| G04 / High | [plane.py:99](../gem/plane.py#L99) | Polygon best-fit normal indexes `i+1` past the final vertex instead of wrapping; a closed four-vertex polygon crashes. |
| C02 / High | [common.py:71](../gem/common.py#L71) | Viewport helper needs `.normalize()` but then subscripts the resulting Vector, which has no indexing API. Formula's intended viewport semantics also need clarification. |

### Core medium-priority findings

| ID / severity | Source | Reproduced behavior / cause |
| --- | --- | --- |
| N01 / Medium | [vector.py:111](../gem/vector.py#L111), [quaternion.py:60](../gem/quaternion.py#L60) | `length is not 0` compares identity, not value. Zero floats enter the division branch. Zero-vector/zero-quaternion normalization raises instead of taking the existing intended zero/identity fallback. |
| N02 / Medium | [vector.py:105](../gem/vector.py#L105), [quaternion.py:53](../gem/quaternion.py#L53) | Naive sum of squares overflows at 1e200 and underflows at 1e−200 despite representable norms. Hypot-style scaling would be a later numerical correction. |
| N03 / Medium | [matrix.py:286](../gem/matrix.py#L286) | Cofactor inversion of a well-conditioned scaled identity fails at 1e−100 (determinant underflow) and produces NaN at 1e100 (overflow). This differs from ill-conditioning; ordinary-scale diagonals near 1e−12 pass. |
| Q07 / Medium | [quaternion.py:248](../gem/quaternion.py#L248) | `toMatrix` changes `.matrix` after construction but leaves `.c_matrix` as identity. A 90° rotation's exported ctypes state is wrong. |
| M07 / Medium | [matrix.py:391](../gem/matrix.py#L391) | Legacy in-place division changes `.matrix` without refreshing `.c_matrix`; separate from division's transpose defect. |
| Q10 / Medium, accuracy contract | [quaternion.py:214](../gem/quaternion.py#L214) | Near-angle SLERP returns unnormalized LERP. Unit inputs 1° apart at t=0.5 give norm 0.9999904807207345. Test uses a 1e−12 unit-norm requirement; preserve the approximation only if explicitly documented and acceptable. |
| R01 / Medium | [ray.py:15](../gem/ray.py#L15) | `duplicate` shares origin and direction and loses original distance (5 becomes 1). Independence and distance preservation have separate probes. Constructor mutation is characterized, not silently changed. |
| V04 / Medium, ownership contract | [vector.py:139](../gem/vector.py#L139), :307 | `clamp` mutates the caller's list while returning a new Vector. Test encodes the class docstring's “new clamped vector” expectation; changing aliasing semantics needs compatibility review. |

### Experimental findings (after core correctness)

| ID / local severity | Source | Reproduced behavior / cause |
| --- | --- | --- |
| E01 / High | [bezier.py:6](../gem/experimental/bezier.py#L6) | Cubic final term is `ttt + p3`, not `ttt*p3`. Control points 0,1,2,3 fail all five endpoint/interior cases. |
| E02 / High | [bezier.py:21](../gem/experimental/bezier.py#L21), :6 | Bezier functions put scalar operands on the left; gem.Vector has no reflected multiplication. Scalar quadratic curves pass, vector curves fail. |
| E03 / High | [bezier.py:41](../gem/experimental/bezier.py#L41), :86, :134 | `/3` gives float curveCount on Python 3; `range` cannot consume it. Control-point validity is unchecked. |
| E04 / High | [bezier.py:88](../gem/experimental/bezier.py#L88), :160 | Sampling fails on absent `sqrMagnitude` or earlier scalar-left multiplication. Further inspection identifies midpoint `(t0−t1)/2`, nonexistent `math.abs`, two-argument `list.append`, unsafe `[-2]` on one sample, and no recursion bound. These downstream issues are source-confirmed but masked by earlier failures in the present end-to-end probes. |
| E05 / High | [legendre.py:15](../gem/experimental/legendre.py#L15), :29 | Every recurrence step reinitializes/multiplies P and recomputes PM1, corrupting degree ≥m+3. P3(0.2) is −0.6 instead of −0.28. Running the same nonzero-order object twice changes the result. Addition theorem fails at l=3 and l=4; low orders pass. |
| E06 / High | [sph_sample.py:52](../gem/experimental/sph_sample.py#L52) | `math.PI` is absent. Subsequent assignment targets `.dir.vec` rather than `.dir.vector`, leaving directions zero. Seeded generation probe fails on the first issue; the second is source-confirmed. |
| E07 / High | [sph_object.py:33](../gem/experimental/sph_object.py#L33) | Uses builtin `object[i]`, leaves coefficient arrays None, rescales only the final vertex, and adds rather than multiplies shadowed scale. Shadow/collision implementation is unfinished. Current probe fails at builtin subscripting. |
| E08 / High | [sph_irradiance_map.py:12](../gem/experimental/sph_irradiance_map.py#L12), :25 | Allocates `[height][width]` but indexes `[width][height]`; 2×1 input crashes. Both angular coordinates and integration weight depend only on width, so rectangular support requires an explicit mapping contract. |

## What passed, and what remains uncertain

Twenty seeded integer matrices per size were checked against recursive determinants and exact Fraction Gauss–Jordan inverses. 3×3 and 4×4 inverses pass both left/right identity checks at ordinary scales. Their cofactor arithmetic has no pivoting or conditioning policy. The inverse4 transpose dependency is deliberate in its current algebra and must be handled together with M02.

Known-answer and property tests also cover vector arithmetic, cross-product orthogonality, barycentric reconstruction for a reference triangle, reflection, LERP endpoints, quaternion associativity/conjugation, rotation length preservation, quaternion/matrix agreement, matrix transpose/composition, FOV/depth mappings, lookAt, low-order Legendre functions, SH addition/orthonormality, and nine-coefficient irradiance updates. Raw ctypes snapshots and Python-3 in-place operations are included.

Source review identified additional design hazards, not all labeled as confirmed bugs:

- Constructors accept inconsistent dimensions and alias caller lists without validation (`Vector.__init__`, `Matrix.__init__`, `Quaternion.__init__`). Kernels can truncate, crash, or propagate nonfinite values. Direct mutation of public `.matrix` bypasses the ctypes refresh.
- Scalar type acceptance is inconsistent; reflected numeric operators and sequence indexing are absent. Adding these is an API extension, not a license to rewrite the types.
- Zero axes, coincident lookAt points, parallel up vectors, degenerate triangles, singular matrices, invalid frusta, empty best-fit data, and invalid SH degrees/orders lack a documented domain-error policy. Singular inverse and collinear barycentric currently raise ZeroDivisionError; tests characterize these instead of choosing new exceptions.
- Float identity comparisons emit SyntaxWarnings. Acos inputs lack clamp/domain guidance. General quaternion log/power, nonunit rotations, and the antipodal no-invert SLERP path need mathematical contracts. No quaternion exponential implementation was found.
- “Same direction” checks dot-product sign (same hemisphere), not parallelism; changing it would be a semantic break. `sign(0)` returns +1. `sinc` approximates small inputs by 1, with measured error below 2e−11 at 1e−5; this is an approximation, not automatically a defect.
- Bezier recursive sampling, SH object transport/shadowing, and most full rendering workflows remain unusable; reaching downstream algorithms requires fixes in Phase 2. Exhaustive invalid-input and extreme-conditioning analysis is not claimed.

## Packaging, compatibility, documentation, and efficiency review

A local wheel builds successfully, but `setup.py` lists only `packages=['gem']`: **no experimental package entries appear in the wheel**. This makes editable/source behavior differ from installed behavior. The wheel uses the existing distribution name `gem` and imports `gem`; names were not changed. Runtime remains pure Python plus six and standard-library ctypes, with no NumPy requirement. Coverage's optional development installation may use its C accelerator; that is not a runtime library dependency.

Metadata advertises Python 2.7 and 3.2–3.5 but lacks `python_requires`, modern build metadata, a test extra, and version support policy. Current legacy Travis configuration only runs `launcher.py test`, which ignores the argument and prints examples without assertions; it is not a test suite. Documentation links reference an older repository identity and off-repository wiki. Several README feature claims exceed the implemented API, and angle/storage/ownership/error conventions are undocumented. These are modernization findings; no packaging or CI redesign was performed in Phase 1.

Efficiency review found list(range(...)) allocation in small matrix kernels, repeated small wrappers and zero buffers, repeated ctypes float32 conversion in every Matrix construction, and extra allocations in quaternion interpolation/conjugation. Vector reference-buffer slicing is already a deliberate optimization. Do not optimize before fixing correctness. The measured 4×4 raw multiply versus object multiply (5.27 vs 7.74 microseconds here) suggests wrapper/conversion cost, but does not isolate its cause; removing public ctypes state could break callers.

See [benchmark results](BENCHMARKS.md) and the [compatibility assessment](COMPATIBILITY.md). No changes to supported Python versions, typing, performance, or packaging are approved or implemented here.

## Reproduction and review

Use [README](README.md) for executable commands and expected exit statuses. Tests contain independent Fraction oracles, deterministic random seeds, known answers, and quadrature; no mandatory NumPy or property-testing dependency is introduced. Remove a finding's defect marker only together with a separately reviewed fix; strict XPASS then flags accidental early corrections.

Review changes as two local commits: (1) tests and optional audit dependencies, (2) audit documentation, measured test results, and benchmark harness/baseline. These are Phase 1 review units only. Subsequent fixes should be separate small PRs grouped by subsystem and dependency (especially M02/inverse4), with their own compatibility notes.

Remote pull-request creation is currently blocked: GitHub Git reads work, but a direct read of `api.github.com` fails at the environment proxy with CONNECT 403, and `gh api` returns Forbidden. API egress access is needed before opening draft PRs; this is not evidence that a new token is required. No default-branch update, merge, or PyPI publication was attempted.
