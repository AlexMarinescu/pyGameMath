# Compatibility impact assessment

## Impact of Phase 1

Only audit tests, optional development-tool requirements, benchmark code, and reports are added. No existing library file, namespace, runtime dependency, manifest, packaging rule, CI workflow, or public API is changed. The test suite requires a modern audit interpreter (the pinned pytest 9 requires Python ≥3.10); this does **not** change gem's supported-version policy. Only Python 3.12 was validated. Existing distribution name `gem` and module import paths are preserved.

Known failures are explicit strict xfails, not disabled assertions. `--runxfail` reproduces their nonzero outcomes. Three expectations need policy decisions: Q11 uses the arbitrary-axis helper's stated rotation-quaternion contract; V04 expects the returning clamp operation not to mutate its input; Q10 requires unit-norm SLERP within 1e−12. They remain failures under those stated expectations, rather than silently approved behavioral changes.

## Review boundaries for later corrections

| Candidate correction | Compatibility exposure | Safe review boundary |
| --- | --- | --- |
| Equality/inequality | Collections and branches may have relied on erroneous first-component comparison; zero/dimension behavior changes | Keep exact component comparison, not approximate equality; decide dimension/empty behavior explicitly. |
| Inverse2, quaternion inverse | Correct numerical answers change outputs for ordinary existing inputs | Preserve types, names, component order; known-answer and both-side identity tests. |
| Angle conversion helpers | Callers may have compensated for swapped names | Do not change any other API's angle unit; document the named-function correction and migration examples. |
| Matrix division | Numeric result changes; Python 3 operators become available | Change inverse4's adjugate construction in the same PR; retain legacy methods and refresh ctypes state. |
| Refraction | False TIR and failing ordinary cases become valid outputs | Keep argument order and ratio convention eta; document unit incident/normal prerequisites. |
| Plane construction and normalization | Fields currently contain both scalars and points; offset sign varies by method | Preserve scalar `a,b,c,d` API with an explicit equation; decide `bestFitD` sign separately. Do not silently reinterpret caller data. |
| Matrix/vector transform, projection | Composition order and homogeneous dimensions affect every caller | Preserve row vectors and storage order. Define project input types, return dimension, and window depth before implementation. |
| Ray copying and transforms | Ownership/mutation and `distance` meaning can affect caller state | Keep `roateUsingMatrix` working; any corrected spelling is additive. Define homogeneous promotion internally and preserve distance semantics. |
| Quaternion rotation helpers, pow/log/SQUAD | Current outputs/types and mixed units are inconsistent; SQUAD has only three controls | Correct crashes separately from choosing general/unit domains or spline mathematics; any signature change requires approval. |
| ctypes synchronization | Refreshing stale buffers changes previously wrong graphics data | Retain the public field and float32 export; do not remove ctypes or automatically redesign mutable storage. |
| Magnitudes, inverse scaling, zero/domain handling | New edge-case behavior or exception types may be observable | Preserve ordinary finite behavior; decide zero, NaN/Inf, degeneracy, and unsupported-size contracts first. |
| SLERP accuracy | Normalizing the nearby branch changes values; component LERP is intentionally unnormalized | Preserve LERP; document whether SLERP promises unit rotations and what input normalization it requires. |
| Experimental algorithms | Mostly currently failing; packing them in wheels exposes previously absent modules | Keep `gem.experimental` paths and legacy misspellings; separate usable algorithms from unfinished shadow transport. |
| Packaging, version support, CI, typing | Installation and downstream support can change even without math changes | Propose only after correctness work; retain pure-Python core and optional external accelerators, if any. |

## Changes that are not authorized by this audit

No switch to column vectors, quaternion `[x,y,z,w]`, one universal angle unit, NumPy-backed storage, compiled mandatory acceleration, renamed distribution/import namespace, automatic approximate equality, or altered default coordinate axes. No library changes were made and no later-phase corrections have been started. Publication, default-branch merge, and breaking API changes require explicit approval.

Every later PR should identify which findings it addresses, show the pre-fix failure, remove only the corresponding strict-xfail annotations, and report newly passing tests alongside all remaining known failures. Source-review-only downstream experimental issues require new focused reproductions before their corrections.

## Phase 1B historical validation update

The [full wiki reconciliation](PHASE1B-WIKI.md) supersedes the implication that every failing Phase 1 expectation is an approved correction. The two substantive wiki pages require matched operand dimensions and distinguish new results from receiver mutation; four other pages are placeholders. No storage/order, unit, or public signature change is authorized.

V04, Q10, Q11, and the implicit-promotion portion of R03 are now explicitly unresolved contract groups. V01 mixed/empty inputs and P01 input/depth subcases are also questions. G01's normal expectation is narrowed to coefficient-aligned direction without choosing raw versus unit normal length. These tests still execute under strict contract-question xfails, separately recorded from confirmed-defect xfails. The original Phase 1 commits/results are preserved.

Matrix default identity and i-method self returns remain compatible with established source. Wiki typos and contradictory example outputs are documentation errors, not approval for behavior changes. Pure-Python runtime requirements, existing gem imports, and PR #9 remain unchanged.

## Phase 2A: V01, M01, Q01

Equal-size, nonempty Vector equality and inequality now examine all components using the existing exact numeric comparisons. For example, `[1,2,3] == [1,9,3]` changes from true to false, and inequality changes from false to true. Branches or collections that relied on the first-component bug will behave differently. No tolerance or approximate equality is introduced; supported positive dimensions and non-Vector `NotImplemented` behavior are preserved. Empty comparisons retain the existing `None` result. Mixed dimensions and malformed storage remain outside the documented equal-size domain: no dimension check, exception policy, or support extension is introduced, and traversing more components can expose an IndexError on mismatched storage. Both V01 contract-question tests remain strict xfails pending QD01.

The 2x2 inverse now returns `[[d,-b],[-c,a]]/(a*d-b*c)`. For `[[1,2],[3,4]]` the result changes from `[[-1.5,0.5],[0,0]]` to `[[-2,1],[1.5,-0.5]]`. The helper, returning Matrix method, and in-place Matrix method share this correction. Nested row storage, ordinary matrix multiplication, row-vector application, ctypes exports, return types, and receiver mutation rules are preserved. Singular matrices still raise ZeroDivisionError; larger inverses and the matrix-division/inverse4 dependency are unchanged.

Quaternion inverse now computes conjugation divided by squared magnitude: `[w,-x,-y,-z]/(w*w+x*x+y*y+z*z)`. This corrects finite ordinary nonzero inputs, including nonunit quaternions as explicitly requested for Phase 2A. For `[0,1,0,0]` the inverse changes to `[0,-1,0,0]`; real-only inputs keep their mathematical values. The `[w,x,y,z]` ordering, Hamilton product, list helper result, Quaternion method result, and input ownership are preserved. Zero still raises ZeroDivisionError. No policy for extreme scales, nonfinite values, other quaternion domains, or rotation helpers is added.

Public signatures, the `gem` namespace, pure-Python implementation, and runtime dependencies are unchanged. No NumPy, Cython, or compiled requirement is added. Only these three confirmed defects are corrected; unrelated defect and contract-question markers remain. See [Phase 2A verification and file summary](PHASE2A.md).
