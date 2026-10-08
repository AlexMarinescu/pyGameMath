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
