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

## Phase 2B: M02, M03, M07

Matrix/scalar division now divides each element in its original position. For `[[2,4],[6,8]] / 2.0`, the result is `[[1,2],[3,4]]`; the legacy helper and method previously produced `[[1,3],[2,4]]`. Callers that compensated for the erroneous transpose should remove that compensation. This correction affects non-symmetric matrices; symmetric matrices retain their values.

`Matrix.__truediv__` and `Matrix.__itruediv__` expose the existing `__div__` and `__idiv__` implementations to Python 3, so `/` and `/=` now work with supported float scalars. Both legacy method names remain available with their existing signatures. Scalar acceptance is preserved: the wrapper accepts floats and delegates other operands with NotImplemented, while the raw helper continues to use the supplied numbers' division behavior (including integer divisors). Integer-wrapper support, reflected division, and broader numeric protocols remain deferred under QD12. Python 2 was not executed; its legacy methods are retained without changing interpreter support metadata.

Returning division constructs a fresh Matrix and leaves the receiver and caller-owned rows unchanged. Both legacy and Python 3 in-place division return the receiver, replace its matrix list, and now regenerate its public float32 `c_matrix` snapshot. Consumers that observed stale pre-division exports will see the corrected values. Direct external list mutation still does not refresh snapshots. Zero float division still raises ZeroDivisionError without changing receiver matrix or export state.

The 4x4 inverse still constructs cofactors, but now explicitly transposes them to form the adjugate before elementwise division by determinant. Its correct ordinary-scale results are preserved independently of the old matrix-division bug. Existing 2x2 and 3x3 inverse algorithms, singular exceptions, and unresolved scaling/conditioning policies are unchanged. Both multiplication orders against identity and an independent Fraction oracle verify the coupled change.

Row-major nested storage, row-vector application, ordinary matrix multiplication, other public signatures, `gem` imports, ctypes export type, and pure-Python dependencies are preserved. Only M02/M03/M07 defect markers are removed; every unrelated expected failure and contract question remains. See [Phase 2B verification](PHASE2B.md).

## Phase 2C: C01, V02

The two named angle helpers now perform the conversions their names specify, using math.pi rather than 3.14. `radiansToDegrees(math.pi)` returns 180 instead of about 0.0548; `degreesToRadians(180)` returns math.pi instead of about 10318.47. Callers that used the wrong helper to compensate for swapped behavior must switch to the correctly named helper, and the old 3.14 approximation is removed. Public function names and signatures, including the misleading legacy keyword names `degrees` for radiansToDegrees and `radians` for degreesToRadians, are preserved. All other degree/radian APIs and their numerical algorithms remain unchanged.

Refraction now uses IOR squared rather than cubed and orders scalar multiplication locally as Vector*scalar, avoiding TypeError without adding reflected operators. For IOR=1.5 and unit incident `[0.6,-0.8,0]` with normal `[0,1,0]`, the result changes from false total internal reflection `[0,0,0]` to `[0.9,-sqrt(0.19),0]`; normal incidence now returns a valid transmitted direction. Total internal reflection continues to return a fresh zero Vector of the input dimension. The function name, argument order and keyword names `refract(IOR, incidentVec, normal)` are unchanged, and caller vectors are not mutated.

Historical documentation did not establish ratio direction or unit/orientation requirements. The user explicitly approved IOR=n1/n2, unit inputs, and a normal pointing into the incident medium and opposing incidence before implementation. Air-to-glass (1.0 to 1.5) uses ratio 1.0/1.5 and a normal toward air; glass-to-air uses 1.5/1.0 and a normal toward glass. Passing an absolute glass index when entering glass would use the wrong ratio. The caller must prepare normalized, correctly oriented vectors; no automatic normalization, flipping, validation exceptions, or near-critical rounding/clamping policy is introduced. See [CONVENTIONS.md](CONVENTIONS.md#phase-2c-current-angle-and-refraction-conventions).

Pure-Python `gem` imports, runtime dependencies, public APIs, other algorithms, and interpreter support policy are preserved. Only C01/V02 strict-defect markers are removed; all unrelated expected failures and unresolved contracts remain. See [Phase 2C validation and file summary](PHASE2C.md).

## Plane construction and normalization

Plane coefficients are consistently scalar values for `a*x+b*y+c*z+d=0`, with `.normal=Vector([a,b,c])`. `fromCoeffs(0,2,0,-4)` now succeeds and retains normal `[0,2,0]`; normalize the plane when a unit normal is needed. `fromPoints` replaces its incompatible Vector-valued coefficient fields with unit-normal scalar coefficients and sets the negative dot-product offset. For points on y=2 with +Y orientation, coefficients are `[0,1,0,-2]`. Code that extracted construction points from `.a/.b/.c` must keep those inputs separately.

Normalization now scales d along with the normal and refreshes `.normal` in both returning and in-place methods. `[0,2,0,-4]` becomes `[0,1,0,-2]`, preserving y=2 rather than moving the plane to y=4. Constructors still return None; returning normalization creates an independent Plane, and in-place normalization returns self. Existing flip, clone, homogeneous dot, and side-classification signatures are preserved.

Polygon normal calculation now wraps the final edge and works with either implicit closure or a repeated first vertex. Winding controls orientation. `bestFitD` keeps its signed mean-dot result, including negative values; polygon construction passes `d=-bestFitD(vertices, normal)`. Its existing sample weighting and nonunit-normal scaling are unchanged. There is no new polygon-constructor method or nonplanar least-squares algorithm.

Zero-normal normalization, collinear points, and empty best-fit inputs retain characterized ZeroDivisionError behavior. Degeneracy, malformed inputs, nonfinite values, and extreme-scale guarantees remain separate policy questions. Only `gem/plane.py` changes in production; other algorithms and dependencies are unchanged. Verification and changed files are listed in [PHASE2D1.md](PHASE2D1.md).

## Matrix translation and vector transformations

In-place Matrix4 translation with Vector3 now constructs the correctly sized helper instead of raising IndexError. Returning and in-place translations both postmultiply the receiver by the row-vector translation matrix; nonidentity receivers therefore retain their existing composition order. The public float32 `c_matrix` snapshot refreshes through `__imul__`. Inputs, signatures, return types, Matrix3 translation variants, and Vector4's ignored fourth offset component are preserved.

Vector transformation now uses `position[j]*matrix[j][i]` and counts each component once. Identity changes from the incorrect `[2,3,5]` back to `[2,3,4]`; a +Z quarter-turn maps `[2,3,4]` to `[-3,2,4]`. Callers that transposed inputs or subtracted the last matrix row to compensate for the old implementation should remove those workarounds. Same-size transforms apply the whole row-vector product rather than adding an unconditional extra final row.

The transform helper supports local w=1 promotion for an N-component affine position with an (N+1)×(N+1) matrix. Thus Vector3 transformation with a translated Matrix4's nested rows now applies translation. Explicit homogeneous positions use their supplied w, directions with w=0 receive no affine translation, and projective matrices compute a potentially changed output w. No perspective divide occurs. Implicit projective output is not perspective-correct Cartesian output; use the appropriate projection/unprojection path for that conversion. P01/P02 remain unresolved. General Matrix*Vector multiplication keeps its existing dimension precondition and receives no promotion.

The free helper returns a list; the returning method returns a fresh Vector of the receiver's size; the in-place method replaces the receiver's storage and returns self. Raw position lists and nested matrix lists remain the accepted representations. Neither method changes a separately supplied position list or matrix. Dimensions 1, 2, 3, 4, and 8 are exercised without restricting the existing generic vector kernel. Broader shape/error and numeric-domain policies remain deferred.

Production changes are confined to `gem/matrix.py` and `gem/vector.py`. Runtime dependencies, pure-Python storage, public signatures, and unrelated algorithms remain unchanged. See [Phase 2D-2 verification](PHASE2D2.md) and [transformation conventions](CONVENTIONS.md#vector-transformations).

## Projection and unprojection

`project` now produces usable window coordinates instead of failing while subscripting Matrix wrappers. It accepts Matrix4 or raw 4×4 nested lists for either matrix argument. `unproject` adds the same raw-list and mixed-form support alongside its existing wrapper form. Project's object input remains explicit Vector4; there is no Vector3 promotion or new position representation. Public signatures, row-major storage, and general multiplication and transform helpers are preserved.

Projection returns exactly `[winx,winy,winz]` in a fresh Vector3, omitting the unusable implementation's fourth reciprocal-W entry. NDC depth now maps from [-1,1] to [0,1], matching unprojection's existing `2*winz-1` input mapping. Identity projection of `[0,0,0,1]` through viewport `[0,0,100,100]` returns `[50,50,0.5]`. Viewports have lower-left origin and upward Y; callers using top-left screen coordinates must convert Y explicitly. Coordinates and depth outside the normal interval are not clamped.

Unprojection now inverts modelview*projection rather than projection*modelview. Identity and commuting matrices retain their ordinary results; noncommuting transforms produce corrected coordinates. In the original translation/orthographic regression, X changes from -5 to the correct 1. Callers that reversed operands or precompensated for the old order should remove that compensation.

Zero clip W raises ZeroDivisionError in project. Unproject preserves its fresh zero-Vector3 sentinel for zero homogeneous output W and the existing singular-inverse ZeroDivisionError. Project does not require an invertible combined matrix. Inputs and public ctypes snapshots are preserved. Numerical conditioning, near-zero-W tolerances, malformed shapes, and invalid viewport policies are not expanded.

Only P01/P02 defect markers and the two now-resolved P01 input/depth question cases are converted. All unrelated expected failures and questions remain. [CONVENTIONS.md](CONVENTIONS.md#projection-and-unprojection) defines the contract; [PHASE2D3.md](PHASE2D3.md) records verification and changed files.

## Pivot rotation and shear

`rotate2(point,theta)` now computes the affine offset needed to keep the pivot fixed. Rotating `(2,3)` about itself by 90 degrees changes from the incorrect `(4,6)` to `(2,3)`; `(3,3)` rotates to `(2,4)`. Zero-angle rotation is now identity for nonzero pivots. The raw point argument, degree units, returned 3×3 nested matrix, and row-vector order are preserved. Code that compensated for the old offset should remove that compensation. Wrap the helper in Matrix3 for homogeneous 2D application; Matrix2 remains origin-only, and Matrix3/Matrix4 rotation methods remain axis-angle operations.

XY3 now writes the Z column instead of raising IndexError. XY4 likewise shears Z instead of modifying homogeneous W. With factors `(1,2)`, `[2,3,4,1]` now becomes `[2,3,12,1]` rather than `[2,3,4,9]`. This changes numerical behavior for callers using XY4 as a projective W modification; it now represents the same geometric shear as XY3. No replacement projective API is introduced. Zero factors retain identity, and supplied W is preserved for positions and directions.

YZ/XZ retain their existing mappings and numerical algorithms. All three size-3 forms are XYZ linear transforms, not 2D homogeneous shears. Returning and in-place methods retain postmultiplication, signatures, input ownership, and float32 ctypes refresh behavior. Only XY3, XY4, and rotate2 change mathematical code in `gem/matrix.py`; the remaining shear helpers receive explanatory docstrings. Runtime dependencies, public imports, and unrelated algorithms remain unchanged.

Only M05/M06 markers are removed. [CONVENTIONS.md](CONVENTIONS.md#pivot-rotation-and-shear) records the mappings and rotation dispatch; [PHASE2D4.md](PHASE2D4.md) records verification and changed files.

## Quaternion/matrix conversions

Matrix-to-quaternion conversion now tracks the largest squared-component candidate and uses its value when reconstructing the quaternion. This corrects both the missing candidate update and the early mutually exclusive selection. Legitimate half-turns about basis or mixed axes return rotations instead of ZeroDivisionError/ValueError. Other inputs with an axis component larger than W now recover the correct unit orientation. Quaternion order [w,x,y,z], row-vector signs, Matrix wrapper input, and Quaternion return type are preserved.

Recovered components can differ by an overall sign from an original quaternion. Compare rotations modulo q/-q equivalence, rather than requiring exact components or positive W. No canonical sign, normalization, or input-validation contract is introduced. Only the upper-left 3×3 block is read, preserving existing Matrix3/Matrix4 behavior and ignoring other entries. Input rows and exports are not modified.

Quaternion-to-matrix conversion now refreshes the returned float32 `c_matrix` after populating the rotation rows. A +Z quarter-turn exports its actual rotation instead of identity. The Python matrix formula and numerical values are unchanged, including existing nonunit and zero-quaternion arithmetic. Each call returns a fresh Matrix4 with independent rows and export storage and preserves the quaternion component list. Direct caller edits still leave snapshots stale until an existing refresh operation runs.

Changes are confined to `gem/quaternion.py`, Q02/Q07 tests, and documentation. Forward axes, handedness, other quaternion algorithms, public signatures, gem imports, and runtime dependencies are unchanged. Proper rotations remain the conversion domain; policies for invalid matrices, nonfinite values, normalization, and other quaternion domains stay separate. See [conversion conventions](CONVENTIONS.md#quaternionmatrix-conversions) and [Phase 2E-1 verification](PHASE2E1.md).

## Quaternion axes and arithmetic

The two axis-angle helpers now accept Python-list axes without calling a Vector method on a list. Vector inputs are normalized in temporary storage rather than replacing the caller's component list. For example, axis `[1,2,2]` remains `[1,2,2]` after conversion while quaternion construction uses `[1/3,2/3,2/3]`. Code relying on the old incidental Vector normalization must normalize its axis explicitly. List and Vector inputs yield the same ordinary numerical results; list-based returned quaternion components use normal mutable list storage.

`quat_rotate_from_axis_angle` retains its historical Quaternion result from rotating the normalized axis about itself. Its list branch now passes a temporary Vector into the existing sandwich; no general list multiplication is added. The Q11 return-semantics question and strict xfail remain. Both axis-angle helpers keep degrees, the X/Y/Z angle helpers keep radians, and zero-axis normalization errors and unsupported-axis NotImplemented behavior are preserved. Other axis functions are unchanged.

In-place Quaternion/Vector multiplication now reads `.vector`, avoiding AttributeError and matching the returning Hamilton product. `[1,2,3,4] * Vector([5,6,7])` yields Quaternion `[-56,2,12,4]` in both forms. This is a product with a pure quaternion, not vector rotation; operand/result types and ordering are unchanged. In-place operations retain receiver identity, replace component storage, and preserve separately held source lists and Vector inputs.

Quaternion `__truediv__`/`__itruediv__` alias the retained legacy division methods. Python 3 float division now works; the float-only wrapper acceptance rule, including float subclasses, is unchanged. Integers and other unsupported operands remain NotImplemented/TypeError rather than becoming accepted scalars. Returning division is independent; in-place division returns self. Zero float division raises ZeroDivisionError without changing data or its list identity. Legacy methods remain callable, but Python 2 was not executed and interpreter support metadata is unchanged.

Only Q06/Q08/Q09 defect markers are removed. Runtime dependencies, pure-Python gem imports, public signatures, Hamilton conventions, and unrelated algorithms remain unchanged. [CONVENTIONS.md](CONVENTIONS.md#quaternion-axes-and-arithmetic) records ownership and operand rules; [PHASE2E2.md](PHASE2E2.md) records verification and changed files.

## Phase 2E-3: quaternion powers and logarithms

Q03 now preserves imaginary components and handles identity without division
by zero. Finite real powers use the principal unit-quaternion angle;
negative identity supports integer parity and rejects fractional powers.
Q04 changes the logarithm's scalar component from one to zero and always
returns independent list storage, including for identity. Zero inputs and
the undefined negative-identity logarithm raise ValueError.

These numerical and error changes replace defective results. Public
signatures, Quaternion power results, list logarithm results, [w,x,y,z]
ordering, and caller ownership are preserved. Unit input is a prerequisite;
there is no implicit normalization, new norm tolerance, general nonunit
support, sign canonicalization, or exponential API. Broader accuracy and
invalid-input policies remain separate decisions.

## Phase 2E-4: quaternion interpolation

Q05 fixes the callable-parameter expression while retaining the legacy
three-control blend, signatures, sign-sensitive no-invert branches and
caller ownership. It is not reinterpreted as four-control SQUAD. Its linear
branches still produce nonunit results and retain antipodal degeneracy.

Q10 replaces the dot>0.999 unnormalized approximation with accurate
spherical interpolation. Near-angle output components and norm change;
there is no implicit input normalization. Stable difference/sum angle
calculation preserves tiny rotations; identical endpoints return fresh
storage. Shortest-path sign correction, half-turn tie behavior, unclamped
parameters and ordinary Quaternion results are preserved. Accuracy is
specified for unit inputs over [0,1]; LERP/no-invert semantics are unchanged.

`squad4(q0,q1,s0,s1,t)` adds a separate public function in gem.quaternion.
All arguments are unit quaternions except t; s0/s1 are explicit SQUAD
controls. It uses accurate shortest-path SLERP, returns fresh storage, and
preserves inputs. No method overload, control-generation or exponential API
is added. See [API examples](CONVENTIONS.md#quaternion-interpolation).

At the old near-angle cutoff, measured norm error falls from about 2.50e-4
to 2.22e-16. The final stable-angle implementation costs about 51% more in
a local one-degree microbenchmark (3.949 versus 2.610 microseconds/call);
this replaces the approximate prototype estimate and is not an application
performance guarantee. Full measurements and reproduction are in
[Phase 2E-4](PHASE2E4.md).


## Phase 2E-5: quaternion API contracts

Q11 documents the historical `quat_rotate_from_axis_angle` result rather
than changing its values: approximately `[0,normalized_axis]`, the pure
Quaternion obtained by rotating that axis about itself. Rotation callers
should use the existing `quat_from_axis_angle` constructor. Signatures,
accepted axes, nonmutating inputs and numerical code are unchanged; no
redundant API or runtime deprecation warning is added.

The former constructor-expectation test now protects the retained legacy
result. Independent known answers distinguish the two helpers and verify
caller storage. README quaternion cross-product/exponential claims are
removed because those APIs do not exist; vector cross products are retained.
The [API guide](../docs/QUATERNIONS.md) consolidates established return types,
units, forward axes, ownership and domains. Existing zero normalization,
extreme norms and other unresolved numerical policies remain unchanged.
