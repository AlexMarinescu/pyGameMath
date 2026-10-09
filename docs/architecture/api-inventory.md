# Current public API inventory

Audited reference: master `dc418923c692b609a9fc611c66433e5950e0a321` after PR #38.
Sources, retained reexports, tests and examples were inspected directly. The
[declaration catalog](#source-declaration-catalog) lists every nonprivate declared
core function/class method, constructors/operators and division aliases. It is an
inventory, not a promise that every callable accepts arbitrary shapes/types.

**Status vocabulary:** supported core means the canonical implemented module;
established contract means behavior documented in decisions and tested in its
stated domain; historical means preserved behavior needing narrower documentation
or compatibility review; transitional means a retained import alias; retired means
removed functionality; proposed means unimplemented roadmap work. No blanket 1.0
stability designation follows from passing tests. No runtime deprecation warnings
or removal date are introduced. Most core modules have no `__all__`; incidental
imports such as math/six are not intended mathematical APIs.

## `gem` package

[Source](../../gem/__init__.py): empty initializer. Import classes/functions from
their modules; there are no root-level Vector/Matrix exports, separate Vector2/3/4
or Matrix2/3/4 classes. Numeric suffixes in discussion denote dimensions of the
existing `Vector(size, data=None)` and `Matrix(size, data=None)` wrappers.
Package names remain `gem` and transitional `gem.experimental`. Packaging/import
tests: [test_core_packaging.py](../../tests/test_core_packaging.py). Version/support
metadata is historical; see [compatibility](compatibility.md).

## `gem.common` — scalar and interoperability utilities

[Source](../../gem/common.py). Implemented historical utility surface:

- `convertArr(l,n)` chunks a flat sequence into slices (a final chunk may be short);
  `list_2d_to_1d` flattens rectangular rows; `convertM4to3` copies the leading 3x3
  raw matrix block. Results are new lists.
- `mulV4` returns four componentwise products, not a dot or matrix product.
- `GLfloat` aliases `ctypes.c_float`; `conv_list`/`conv_list_2d` build new ctypes
  arrays using a supplied ctypes type. They do not require OpenGL.
- `sinc(x)` returns 1 for abs(x)<1e-4, otherwise sin(x)/x; this historical small-angle
  approximation is not the angular-probe integration implementation.
- `scalarLerp(a,b,time)` returns unclamped linear interpolation. `sign(x)` returns
  +1.0 for x>=0, including zero, and -1.0 otherwise.
- Corrected angle conversions retain misleading keyword names.
  `getViewPort` takes a Vector and dimensions, returns a fresh rectangle list using
  the unusual whole-vector-normalization/XY-offset formula, and preserves inputs.

Constraints/gaps: raw sequence shape/type/error policies are not unified. Viewport
is not project/unproject or glViewport. Tests:
[test_vector_common.py](../../tests/test_vector_common.py),
[test_angles_refraction.py](../../tests/test_angles_refraction.py),
[test_final_vector_contracts.py](../../tests/test_final_vector_contracts.py),
[test_api_edges.py](../../tests/test_api_edges.py) and
[test_wiki_contracts.py](../../tests/test_wiki_contracts.py).
Examples: [viewport guide](../VECTOR_VIEWPORT_CONTRACTS.md); no standalone utility
tutorial. Raw ctypes/conversion helpers need a fuller API reference.

## `gem.vector` — component vectors and geometric helpers

[Source](../../gem/vector.py). `Vector` exposes `.size` and `.vector`; construction
retains supplied storage, default construction allocates zeros. Operators support
Vector addition/subtraction and int/float scalar arithmetic; reflected scalar
operators are not supplied. Exact comparison includes dimension/empty semantics.
Returning arithmetic/clone/normalize/clamp/transform uses independent storage;
augmented assignment and i-prefixed methods mutate the receiver. `zero`/`one`
also mutate and return self despite lacking the i prefix.

Important groups:

- List kernels `zero_vector`, `one_vector`, `vec_add/sub/mul/div/neg`, `s_vec_add/sub`,
  `dot`, `magnitude`, `normalize`, `maxV/minV/maxS/minS`; sizes are explicit.
  Dot/magnitude/scalar extrema return numbers; the other raw kernels return lists.
- Wrapper extrema, direction-sign predicates, swizzles xy/yz/xz/xw/yw/zw and
  xyw/yzw/xzw/xyz; swizzles require sufficient components and return new Vectors.
- `cross` returns Vector3; `reflect`, `refract`, `lerp` accept Vectors. Refraction
  has the unit-vector n1/n2 contract and zero total-internal-reflection sentinel.
- `toAngle`, `lperp`, `rperp` take raw 2D indexable sequences; toAngle returns
  radians, perpendicular helpers return Vector2.
- `Vector.barycentric(a,b,c)` returns a fresh three-weight list via dot products,
  with ZeroDivisionError for degenerate triangles. It already exists; future
  triangle work extends geometry coverage rather than introducing barycentrics
  as an entirely absent feature.
- `transform` uses raw position/nested matrix lists, local affine promotion and
  no perspective division. `clamp` explicitly takes size/value/bound lists;
  receiver values are not implicit arguments to that method.
- right/left/up/down/front/back return fixed Vector3 directions; front is -Z.

Generic component operations accept well-formed sizes beyond 2/3/4, but specialized
geometry is dimension-specific. Constructors do not validate supplied length;
cross-dimension arithmetic errors remain historical. Stable norms/normalization do
not stabilize every dot/cross/barycentric calculation. Mutable implementation
reference lists `REFRENCE_VECTOR_2/3/4` and `IREFRENCE_VECTOR_2/3/4` remain module-visible
historical constants, not an invitation to change global default buffers.

Tests: [test_vector_common.py](../../tests/test_vector_common.py),
[test_transformations.py](../../tests/test_transformations.py),
[test_angles_refraction.py](../../tests/test_angles_refraction.py),
[test_final_vector_contracts.py](../../tests/test_final_vector_contracts.py),
[test_numerical_robustness.py](../../tests/test_numerical_robustness.py),
[test_vector_quaternion_optimization.py](../../tests/test_vector_quaternion_optimization.py).
Examples: [launcher.py](../../launcher.py), [ownership/viewport guide](../VECTOR_VIEWPORT_CONTRACTS.md),
and the HDR example uses Vector3. Gaps: complete method/reference pages, dimension
and scalar protocol review, degenerate barycentric/invalid-data policies.

## `gem.matrix` — small matrices, transformations and projection

[Source](../../gem/matrix.py). `Matrix` exposes `.size`, caller-retained `.matrix`
and a fresh float32 `.c_matrix` snapshot; default storage is identity. Matrix
products are ordinary products; Matrix*Vector computes row-vector application.
Wrapper multiplication accepts Matrix/Vector, not a general scalar product.
Division accepts floats only (legacy methods plus Python 3 aliases); raw
`matrix_div` uses the operand's arithmetic. Supported in-place operations return
self, replace rows and refresh ctypes; returning variants preserve inputs.

Important groups:

- Raw list construction/arithmetic: `zero_matrix`, `identity`, `scale`,
  `matrix_multiply`, `matrix_vector_multiply`, `matrix_div`, `transpose`.
  Raw matrix/vector multiplication takes nested rows and a Vector, returns Vector.
- `det2/3/4`, `inverse2/3/4` operate on raw nested lists; wrapper `det`, `inverse`,
  `i_inverse` dispatch for dimensions 2/3/4 only. Matrix3/4 scaled inverses preserve
  numerical robustness and singular ZeroDivisionError; public determinants do
  not have equivalent extreme-scale guarantees.
- `translate2/3/4`, `rotate2/3/4`, `rotate_origin2` return raw matrices.
  rotate2 is homogeneous pivot rotation; rotate_origin2 is a radians 3x3 helper.
  Matrix2 wrapper rotation is origin-only, Matrix3/4 axis-angle accepts Vector.
  Legacy translate3 is last-row replacement, not general affine 3D translation.
- scale/rotate/translate/transpose and all three shear planes have returning and
  in-place wrappers. Shear helpers exist in 3x3/4x4 forms; size-4 preserves w.
- `orthographic`, `perspective`, `perspectiveX`, `lookAt` return Matrix4. FOV is
  degrees, aspect=width/height, negative-Z camera and OpenGL NDC depth.
- `project` accepts explicit Vector4 and Matrix4/raw rows, returning Vector3.
  `unproject` takes scalar window XYZ, the same matrix forms and viewport,
  returning Vector3. Viewport is lower-left `[x,y,width,height]`, window depth
  [0,1], without clamping. Their zero-w behavior intentionally differs.

General Matrix*Vector does not implicitly promote or uniformly validate dimension
mismatches. No arbitrary-size inverse, generalized nonfinite policy, invalid-frustum
validator or automatic ctypes tracking of direct row edits is established.
Tests: [test_matrix.py](../../tests/test_matrix.py),
[test_matrix_division.py](../../tests/test_matrix_division.py),
[test_pivot_shear.py](../../tests/test_pivot_shear.py),
[test_projection.py](../../tests/test_projection.py),
[test_inverse_optimization.py](../../tests/test_inverse_optimization.py),
[test_transformations.py](../../tests/test_transformations.py).
Examples: launcher, [current conventions](conventions.md), quaternion/ray guides.
Gaps: full wrapper/kernel reference, legacy translation/shape constraints and
external upload examples verified against a real graphics consumer.

## `gem.quaternion` — algebra, rotations and interpolation

[Source](../../gem/quaternion.py); [public guide](../QUATERNIONS.md).
`Quaternion(data=None)` retains supplied `.data`, with default identity in
[w,x,y,z] order. Addition/subtraction accept Quaternion; multiplication accepts
Quaternion, Vector or float; scalar division accepts float only. Quaternion*Vector
returns a Quaternion product. In-place variants return the receiver and replace
storage. Returning algebra and orientation results are independent.

Important groups:

- Raw-list kernels `quat_identity/add/sub/mul_quat/mul_vect/mul_float/div_float/neg`,
  dot/magnitude/normalize/conjugate/inverse return lists except scalar dot/norm.
  Inverse uses conjugate/norm²; normalization is stable with identity zero fallback,
  but inverse's squared sum is not extreme-scale stabilized.
- `quat_from_axis_angle` returns the recommended rotation Quaternion, with degrees
  and nonmutating Vector3/list axis normalization. `quat_rotate_from_axis_angle`
  is retained legacy pure-axis output. `quat_rotate` takes raw unit axis coordinates
  and a Vector3 point, returning Vector3; X/Y/Z angle helpers take radians and
  return four-component lists. `quat_rotate_vector` applies a unit sandwich.
- `quat_pow`/`Quaternion.pow` return fresh Quaternion on the unit domain;
  `quat_log`/`Quaternion.log` return a fresh four-element list. Nonunit inputs are
  unsupported, not consistently rejected by validation. Zero/negative-identity
  branches and principal angle are documented in conventions.
- LERP, accurate shortest-path SLERP, historical no-invert SLERP and legacy
  three-control SQUAD have wrapper/free entry points. The free `squad4` adds the
  separate conventional four-control API; there is no Quaternion.squad4 method.
- `quat_to_matrix`/`toMatrix` return Matrix4 with synchronized export;
  `quat_from_matrix` accepts a Matrix wrapper's proper rotation block. Neither
  implicitly normalizes/orthogonalizes inputs. Unit orientations permit sign equivalence.
- getForward/getBack/getLeft/getRight/getUp/getDown return rotated Vector3 axes;
  identity forward is +Z, distinct from Vector.front.

There is no quaternion exponential, cross-product or intermediate-control-generation
API. No runtime deprecation is applied to the historical helper. Tests:
[test_quaternion.py](../../tests/test_quaternion.py),
[test_quaternion_contracts.py](../../tests/test_quaternion_contracts.py),
[test_quaternion_operations.py](../../tests/test_quaternion_operations.py),
[test_quaternion_matrix.py](../../tests/test_quaternion_matrix.py),
[test_quaternion_powers.py](../../tests/test_quaternion_powers.py),
[test_quaternion_interpolation.py](../../tests/test_quaternion_interpolation.py).
Examples: public guide, conventions, ray guide and analytical SH/HDR rotation.
Gaps: a complete function reference and explicit unsupported-domain/scalar policies;
preserve established units, legacy SQUAD and sign handling during documentation work.

## `gem.plane` — coefficient planes and polygon approximation

[Source](../../gem/plane.py). `Plane()` starts with zero scalar a/b/c/d and zero
Vector3 normal. `.normal` matches coefficient scale after supported construction/
normalization; arbitrary public field edits do not synchronize other fields.
fromCoeffs/fromPoints mutate and return None. clone/flip/normalize return fresh
Plane; i_flip/i_normalize return self. dot requires Vector4 and uses supplied w.
`point_location(self,plane,point)` keeps its unusual separate plane argument and
indexable XYZ input, returning -1/0/+1 for ordinary finite values.

`bestFitNormal` returns unit Newell Vector3 with vertex wrapping;
`bestFitD` returns signed D=mean(n·p), so construction uses d=-D. Helpers do not
change the receiver. Nonplanar input is an approximation, not a least-squares
guarantee. Raw `flip` takes `[a,b,c,d,normal]`, returns a list; raw `normalize`
takes coefficients and returns a four-coefficient tuple. Zero-normal and degenerate
construction retain ZeroDivisionError. Broader degeneracy/extreme/nonfinite policy
is unresolved. Tests: [test_planes.py](../../tests/test_planes.py),
[test_plane_ray.py](../../tests/test_plane_ray.py), numerical robustness tests.
Example: audit/current conventions. Gap: dedicated practical plane tutorial,
point-location reference and clarified wider validation policy.

## `gem.ray` — stored ray geometry and rigid transforms

[Source](../../gem/ray.py); [public guide](../RAYS.md).
`Ray(startVector,dirVector)` retains caller Vectors, stores original direction
length as `.distance` and normalizes direction in place. `.end` starts as zero
intersection placeholder/state. `duplicate` deep-copies all Vector fields and
preserves exact distance without construction. `.start`, `.dir`, `.end`, `.distance`
are public mutable state.

`roateUsingMatrix` (historical spelling) applies Matrix3 rotation about the origin;
`rotateUsingQuaternion` uses a unit Hamilton sandwich. Both replace start/direction
and normalize direction. `translate` locally promotes Vector3/Matrix4 position and
direction with w=1/0 for pure translation; all transforms preserve distance and
leave `.end` untouched. Methods mutate the receiver and return None; output prints.
Matrix/quaternion inputs and previously referenced transformed Vectors are preserved.

No intersection methods, hit-validity model or general scale/shear/projective ray
contract exists. External hit state must be managed explicitly; zero/nonzero end
values cannot reliably distinguish a placeholder from an actual hit. Tests:
[test_rays.py](../../tests/test_rays.py), [test_plane_ray.py](../../tests/test_plane_ray.py),
numerical robustness tests. Example: ray guide. Gap: explicit intersection/distance
design before adding future primitive queries; no renamed method in this phase.

## `gem.bezier` — evaluation, cubic paths and adaptive sampling

[Source](../../gem/bezier.py). Canonical supported functions
`quadraticBezierPoint`/`cubicBezierPoint` return scalar or new Vector using Bernstein
evaluation, without clamping t or mutating controls. Evaluation is not limited
to the sampling API's Vector2/3 validation.

`BezierPath` preserves the spelling `calculateBezerPoint`. setControlPoints retains
the supplied list and returns None; getControlPoints exposes it directly.
curveCount=(len(controls)-1)//3, with valid cubic layout 3k+1. Public historical
fields include controlPoints, curveCount, minimum_sqr_distance, segments_per_curve
and misspelled divison_threshold; the latter two do not govern corrected adaptive
subdivision and must not be described as active sampling controls.

`findDrawingPoints` returns a fresh ordered list of scalar/Vector samples including
endpoints. `getDrawingPoints` returns nested per-segment lists, omitting duplicate
shared boundaries. `findDrawingPointsAdded` inserts fresh interior subinterval
samples into caller pointList and returns inserted count. Sampling checks finite
scalar or matching Vector2/3 controls, positive finite squared tolerance and
depth-16 best effort with control-to-chord-segment flatness.

`interpolate` intentionally appends generated controls; `samplePoints` rebuilds
using squared-distance thinning heuristics, retains endpoints and preserves inputs.
Both return None and no-op below two source points. Accumulating independent
control sets is not automatically a valid connected 3k+1 path. Malformed sampling
layouts, invalid tolerance/intervals and indices have explicit errors; wider
evaluation overflow behavior is not newly defined.

Tests: [test_bezier.py](../../tests/test_bezier.py),
[test_bezier_sampling.py](../../tests/test_bezier_sampling.py),
[test_bezier_sh_optimization.py](../../tests/test_bezier_sh_optimization.py).
Examples: performance harness and [Phase 2F-3B report](../../audit/PHASE2F3B.md),
not a standalone public curve tutorial. Gaps: tutorial, return-shape examples,
builder-vs-sampling distinction and legacy field guidance in a full reference.

## `gem.legendre` — ordinary and associated Legendre functions

[Source](../../gem/legendre.py). `Legendre(l,m,x).run()` returns a scalar, preserving
scratch fields P/PM1/PML. l/m/x remain caller-settable numeric state.
`mGreaterThan0`, `calculatePM1`, `calculatePML(i)` explicitly update scratch fields
and return None; repeated helper initialization is deterministic.
Associated domain: integer 0<=m<=l, x in [-1,1], unnormalized with Condon–Shortley
phase. Ordinary m=0 evaluation also supports x outside that interval. Invalid,
negative-order and extreme-order behavior is not generalized; in particular
l<m retains the historical scratch PML result rather than new validation.

Tests: [test_legendre.py](../../tests/test_legendre.py) uses independent Rodrigues
references through degree 12, boundaries, parity, recurrence, repeated state and
SH addition theorem. [test_experimental.py](../../tests/test_experimental.py) retains
low-order SH checks. Consumer/example: core SH and its guide/HDR workflow.
Gap: dedicated polynomial reference/example and future high-order numeric policy.

## `gem.spherical_harmonics` — real SH and environment lighting

[Source](../../gem/spherical_harmonics.py); [public guide](../SPHERICAL_HARMONICS.md).
Canonical orthonormal Condon–Shortley basis uses index l(l+1)+m, theta/phi radians.
Factorial/K/SPH are historical scalar helpers with stated integer domains and
no broad invalid/nonfinite validation. Legendre remains module-visible because
compatibility sph imports historically exposed it; canonical polynomial import
is gem.legendre.

- `SPHSample` stores theta/phi, mutable values list and a supplied Vector reference
  as dir (non-Vector input produces zero Vector3). GenerateSamples creates jittered
  equal-solid-angle samples with global RNG; seed externally for repeatability.
- `project_radiance` integrates RGB against precomputed complete basis arrays.
  Default weights mean uniform sphere; explicit nonnegative finite solid angles
  are not renormalized. Output is a fresh coefficient-by-RGB list.
- `project_angular_probe` accepts rectangular row/column/RGB angular disks,
  pixel-center mapping and Jacobian weighting; not latitude-longitude or mirrored
  balls. Image decoding is outside core.
- `reconstruct` accepts canonical RGB complete bands and a unit Vector3/XYZ triple,
  returns RGB without normalizing or convolving. `convolve_diffuse` creates separate
  first-three-band irradiance coefficients. `legacy_to_canonical` converts exactly
  nine legacy RGB entries with sign and rounded-scale correction.
- `rotate_coefficients` analytically rotates canonical scalar/RGB arrays of 1/4/9
  entries, active f'(d)=f(R^-1d), with finite Quaternion norm drift <=1e-12 and
  temporary normalization. Inputs are preserved; no resampling or Matrix orientation
  adapter is provided.
- Historical `SPH_IrradianceMapCoeff(fileU,width,height)` reads native float32
  RGB angular-disk data, storing hdr and nine **legacy radiance** coeffs despite its
  name. load/calculateCoefficients rebuild; updateCoefficients accumulates;
  output prints. Constructor performs I/O; malformed dimensions/short files reject.

New coefficient APIs explicitly validate finite/layout requirements, while unit
reconstruction direction is a prerequisite rather than a full norm-tolerance
check. Extreme basis orders and universal integration bounds are not established.
The basis-layout cache is bounded/private, not a public cache-management API.

Tests: [test_spherical_harmonics.py](../../tests/test_spherical_harmonics.py),
[test_sh_rotation.py](../../tests/test_sh_rotation.py),
[test_hdr_sh_example.py](../../tests/test_hdr_sh_example.py),
[test_experimental.py](../../tests/test_experimental.py), Legendre/optimization tests.
Example: [complete headless HDR/GLSL reference](../../examples/hdr_sh/README.md),
procedural fixture and committed numerical/image outputs. Example latitude-longitude
and RGBE adapters are not core APIs or a renderer. Gaps: full reference, higher-band
accuracy guidance and maintained integration tutorials. Visibility/shadow transport
is retired, not provided by SH projection.

## Transitional imports and retired code

| Historical path | Current exports/replacement | Status/action |
|---|---|---|
| `gem.experimental` | Empty package marker | Retained for imports below; directory not fully eliminated |
| `gem.experimental.bezier` | cubicBezierPoint, quadraticBezierPoint, BezierPath from gem.bezier | Thin reexport; use core |
| `gem.experimental._bezier_legacy` | BezierPath from gem.bezier | Private historical alias retained; use core |
| `gem.experimental.legendre` | Legendre from gem.legendre | Thin reexport; use core |
| `gem.experimental.sph` | Factorial, K, SPH, incidental Legendre | Thin reexport; use SH/polynomial core modules |
| `gem.experimental.sph_sample` | SPHSample, GenerateSamples | Thin reexport; use gem.spherical_harmonics |
| `gem.experimental.sph_irradiance_map` | SPH_IrradianceMapCoeff | Thin reexport; import migration does not convert basis |
| `gem.experimental.sph_object` | No replacement for SPHVertex/SPHObject/GenereateCoeffs | Removed unfinished E07; importing fails |

Shims preserve object identity/signatures with no algorithm duplication. Their
removal boundary/window remains unresolved. There are **no remaining experimental
algorithm implementations** to treat as stable, but the compatibility namespace
continues to ship. Retired tests/source evidence remain archived; they were not
silently skipped. See [migration guidance](../EXPERIMENTAL_MIGRATION.md),
[Phase 2F-6 disposition](../../audit/PHASE2F6.md),
[test_core_packaging.py](../../tests/test_core_packaging.py) and compatibility tests
in Bezier/Legendre/SH suites. External consumers/serialized retired-class references
cannot be inferred from repository consumers alone.

## Examples, tools and gaps across the package

[launcher.py](../../launcher.py) prints Vector/Matrix examples; its ignored `test`
argument does not run assertions. [examples/hdr_sh](../../examples/hdr_sh/README.md)
is a source-tree, CPU/headless environment-lighting reference with GLSL formulas,
not an installed core renderer or GPU validation. [benchmarks](../../benchmarks/README.md)
are separate development tools, not supported mathematical runtime APIs.

The historical [wiki snapshot](../../audit/wiki-snapshot/Home.md) is evidence,
not the current authoritative reference: several pages are placeholders or contain
corrected defects. The wiki/README and packaging are unchanged here. Future
documentation should specify all return shapes, numeric prerequisites and
ownership exceptions before promoting additional APIs to a 1.0 stability promise.
Remaining decisions are centralized on the [development page](../development/roadmap.md#open-decisions).

## Source declaration catalog

The following signatures mirror source, including explicit self on methods.
Underscore helpers are private and omitted, except constructors/operators that
define wrapper behavior. Division aliases are included. Declarations describe
implemented names only; status and constraints are given above.

### gem.bezier declarations

| Kind | Source signature |
|---|---|
| Function | `cubicBezierPoint(t, p0, p1, p2, p3)` |
| Function | `quadraticBezierPoint(t, p0, p1, p2)` |
| Class | `BezierPath` |
| Method | `BezierPath.__init__(self)` |
| Method | `BezierPath.setControlPoints(self, newControlPoints)` |
| Method | `BezierPath.getControlPoints(self)` |
| Method | `BezierPath.calculateBezerPoint(self, curveIndex, t)` |
| Method | `BezierPath.interpolate(self, segmentPoints, scale)` |
| Method | `BezierPath.samplePoints(self, sourcePoints, minSqrDistance, maxSqrDistance, scale)` |
| Method | `BezierPath.getDrawingPoints(self)` |
| Method | `BezierPath.findDrawingPoints(self, curveIndex)` |
| Method | `BezierPath.findDrawingPointsAdded(self, curveIndex, t0, t1, pointList, insertionIndex)` |

### gem.common declarations

| Kind | Source signature |
|---|---|
| Function | `convertArr(l, n)` |
| Function | `mulV4(v1, v2)` |
| Function | `conv_list(listIn, cType)` |
| Function | `conv_list_2d(listIn, cType)` |
| Function | `list_2d_to_1d(inlist)` |
| Function | `convertM4to3(matrix)` |
| Function | `sinc(x)` |
| Function | `scalarLerp(a, b, time)` |
| Function | `getViewPort(coords, width, height)` |
| Function | `radiansToDegrees(degrees)` |
| Function | `degreesToRadians(radians)` |
| Function | `sign(x)` |

### gem.legendre declarations

| Kind | Source signature |
|---|---|
| Class | `Legendre` |
| Method | `Legendre.__init__(self, l, m, x)` |
| Method | `Legendre.mGreaterThan0(self)` |
| Method | `Legendre.calculatePM1(self)` |
| Method | `Legendre.calculatePML(self, i)` |
| Method | `Legendre.run(self)` |

### gem.matrix declarations

| Kind | Source signature |
|---|---|
| Function | `zero_matrix(size)` |
| Function | `identity(size)` |
| Function | `scale(size, value)` |
| Function | `matrix_multiply(matrixA, matrixB)` |
| Function | `matrix_vector_multiply(matrix, vec)` |
| Function | `matrix_div(mat, scalar)` |
| Function | `transpose(mat)` |
| Function | `shearXY3(x, y)` |
| Function | `shearYZ3(y, z)` |
| Function | `shearXZ3(x, z)` |
| Function | `shearXY4(x, y)` |
| Function | `shearYZ4(y, z)` |
| Function | `shearXZ4(x, z)` |
| Function | `translate2(vector)` |
| Function | `translate3(vector)` |
| Function | `translate4(vector)` |
| Function | `rotate2(point, theta)` |
| Function | `rotate3(axis, theta)` |
| Function | `rotate4(axis, theta)` |
| Function | `rotate_origin2(theta)` |
| Function | `det2(mat)` |
| Function | `det3(mat)` |
| Function | `det4(mat)` |
| Function | `inverse2(mat)` |
| Function | `inverse3(mat)` |
| Function | `inverse4(mat)` |
| Class | `Matrix` |
| Method | `Matrix.__init__(self, size, data=None)` |
| Method | `Matrix.__mul__(self, other)` |
| Method | `Matrix.__imul__(self, other)` |
| Method | `Matrix.__div__(self, other)` |
| Method | `Matrix.__idiv__(self, other)` |
| Alias | `Matrix.__truediv__ = Matrix.__div__` |
| Alias | `Matrix.__itruediv__ = Matrix.__idiv__` |
| Method | `Matrix.i_scale(self, value)` |
| Method | `Matrix.scale(self, value)` |
| Method | `Matrix.det(self)` |
| Method | `Matrix.i_inverse(self)` |
| Method | `Matrix.inverse(self)` |
| Method | `Matrix.i_rotate(self, axis, theta)` |
| Method | `Matrix.rotate(self, axis, theta)` |
| Method | `Matrix.i_translate(self, vecA)` |
| Method | `Matrix.translate(self, vecA)` |
| Method | `Matrix.i_transpose(self)` |
| Method | `Matrix.transpose(self)` |
| Method | `Matrix.shearXY(self, x, y)` |
| Method | `Matrix.i_shearXY(self, x, y)` |
| Method | `Matrix.shearYZ(self, y, z)` |
| Method | `Matrix.i_shearYZ(self, y, z)` |
| Method | `Matrix.shearXZ(self, x, z)` |
| Method | `Matrix.i_shearXZ(self, x, z)` |
| Function | `orthographic(left, right, bottom, top, zNear, zFar)` |
| Function | `perspective(fov, aspect, znear, zfar)` |
| Function | `perspectiveX(fov, aspect, znear, zfar)` |
| Function | `lookAt(eye, center, up)` |
| Function | `project(obj, model, proj, viewport)` |
| Function | `unproject(winx, winy, winz, modelview, projection, viewport)` |

### gem.plane declarations

| Kind | Source signature |
|---|---|
| Function | `flip(plane)` |
| Function | `normalize(pdata)` |
| Class | `Plane` |
| Method | `Plane.__init__(self)` |
| Method | `Plane.clone(self)` |
| Method | `Plane.fromCoeffs(self, a, b, c, d)` |
| Method | `Plane.fromPoints(self, a, b, c)` |
| Method | `Plane.i_flip(self)` |
| Method | `Plane.flip(self)` |
| Method | `Plane.dot(self, vec)` |
| Method | `Plane.i_normalize(self)` |
| Method | `Plane.normalize(self)` |
| Method | `Plane.bestFitNormal(self, vecList)` |
| Method | `Plane.bestFitD(self, vecList, bestFitNormal)` |
| Method | `Plane.point_location(self, plane, point)` |

### gem.quaternion declarations

| Kind | Source signature |
|---|---|
| Function | `quat_identity()` |
| Function | `quat_add(quat, quat1)` |
| Function | `quat_sub(quat, quat1)` |
| Function | `quat_mul_quat(quat, quat1)` |
| Function | `quat_mul_vect(quat, vect)` |
| Function | `quat_mul_float(quat, scalar)` |
| Function | `quat_div_float(quat, scalar)` |
| Function | `quat_neg(quat)` |
| Function | `quat_dot(quat1, quat2)` |
| Function | `quat_magnitude(quat)` |
| Function | `quat_normalize(quat)` |
| Function | `quat_conjugate(quat)` |
| Function | `quat_inverse(quat)` |
| Function | `quat_from_axis_angle(axis, theta)` |
| Function | `quat_rotate(origin, axis, theta)` |
| Function | `quat_rotate_x_from_angle(theta)` |
| Function | `quat_rotate_y_from_angle(theta)` |
| Function | `quat_rotate_z_from_angle(theta)` |
| Function | `quat_rotate_from_axis_angle(axis, theta)` |
| Function | `quat_rotate_vector(quat, vec)` |
| Function | `quat_pow(quat, exp)` |
| Function | `quat_log(quat)` |
| Function | `quat_lerp(quat0, quat1, t)` |
| Function | `quat_slerp(quat0, quat1, t)` |
| Function | `quat_slerp_no_invert(quat0, quat1, t)` |
| Function | `quat_squad(quat0, quat1, quat2, t)` |
| Function | `squad4(q0, q1, s0, s1, t)` |
| Function | `quat_to_matrix(quat)` |
| Class | `Quaternion` |
| Method | `Quaternion.__init__(self, data=None)` |
| Method | `Quaternion.__add__(self, other)` |
| Method | `Quaternion.__iadd__(self, other)` |
| Method | `Quaternion.__sub__(self, other)` |
| Method | `Quaternion.__isub__(self, other)` |
| Method | `Quaternion.__mul__(self, other)` |
| Method | `Quaternion.__imul__(self, other)` |
| Method | `Quaternion.__div__(self, other)` |
| Method | `Quaternion.__idiv__(self, other)` |
| Alias | `Quaternion.__truediv__ = Quaternion.__div__` |
| Alias | `Quaternion.__itruediv__ = Quaternion.__idiv__` |
| Method | `Quaternion.i_negate(self)` |
| Method | `Quaternion.negate(self)` |
| Method | `Quaternion.i_identity(self)` |
| Method | `Quaternion.identity(self)` |
| Method | `Quaternion.magnitude(self)` |
| Method | `Quaternion.dot(self, quat2)` |
| Method | `Quaternion.i_normalize(self)` |
| Method | `Quaternion.normalize(self)` |
| Method | `Quaternion.i_conjugate(self)` |
| Method | `Quaternion.conjugate(self)` |
| Method | `Quaternion.inverse(self)` |
| Method | `Quaternion.pow(self, e)` |
| Method | `Quaternion.log(self)` |
| Method | `Quaternion.lerp(self, quat1, time)` |
| Method | `Quaternion.slerp(self, quat1, time)` |
| Method | `Quaternion.slerp_no_invert(self, quat1, time)` |
| Method | `Quaternion.squad(self, quat1, quat2, time)` |
| Method | `Quaternion.toMatrix(self)` |
| Method | `Quaternion.getForward(self)` |
| Method | `Quaternion.getBack(self)` |
| Method | `Quaternion.getLeft(self)` |
| Method | `Quaternion.getRight(self)` |
| Method | `Quaternion.getUp(self)` |
| Method | `Quaternion.getDown(self)` |
| Function | `quat_from_matrix(matrix)` |

### gem.ray declarations

| Kind | Source signature |
|---|---|
| Class | `Ray` |
| Method | `Ray.__init__(self, startVector, dirVector)` |
| Method | `Ray.duplicate(self)` |
| Method | `Ray.roateUsingMatrix(self, matrix)` |
| Method | `Ray.rotateUsingQuaternion(self, quat1)` |
| Method | `Ray.translate(self, matrix)` |
| Method | `Ray.output(self)` |

### gem.spherical_harmonics declarations

| Kind | Source signature |
|---|---|
| Function | `Factorial(n)` |
| Function | `K(l, m)` |
| Function | `SPH(l, m, theta, phi)` |
| Class | `SPHSample` |
| Method | `SPHSample.__init__(self, theta, phi, dirc, sampleNumber)` |
| Function | `GenerateSamples(sqrtNumSamples, numBands)` |
| Function | `project_radiance(samples, radiances, weights=None)` |
| Function | `project_angular_probe(hdr, numBands=3)` |
| Function | `reconstruct(coefficients, direction)` |
| Function | `convolve_diffuse(radiance_coefficients)` |
| Function | `legacy_to_canonical(coefficients)` |
| Class | `SPH_IrradianceMapCoeff` |
| Method | `SPH_IrradianceMapCoeff.__init__(self, fileU, width, height)` |
| Method | `SPH_IrradianceMapCoeff.load(self)` |
| Method | `SPH_IrradianceMapCoeff.calculateCoefficients(self)` |
| Method | `SPH_IrradianceMapCoeff.updateCoefficients(self, hdr, domega, x, y, z)` |
| Method | `SPH_IrradianceMapCoeff.output(self)` |
| Function | `rotate_coefficients(coefficients, orientation)` |

### gem.vector declarations

| Kind | Source signature |
|---|---|
| Function | `zero_vector(size)` |
| Function | `one_vector(size)` |
| Function | `lerp(vecA, vecB, time)` |
| Function | `cross(vecA, vecB)` |
| Function | `reflect(incidentVec, normal)` |
| Function | `refract(IOR, incidentVec, normal)` |
| Function | `toAngle(vector)` |
| Function | `lperp(vector)` |
| Function | `rperp(vector)` |
| Function | `vec_add(size, vecA, vecB)` |
| Function | `s_vec_add(size, vecA, scalar)` |
| Function | `vec_sub(size, vecA, vecB)` |
| Function | `s_vec_sub(size, vecA, scalar)` |
| Function | `vec_mul(size, vecA, scalar)` |
| Function | `vec_div(size, vecA, scalar)` |
| Function | `vec_neg(size, vecA)` |
| Function | `dot(size, vecA, vecB)` |
| Function | `magnitude(size, vecA)` |
| Function | `normalize(size, vecA)` |
| Function | `maxV(size, vecA, vecB)` |
| Function | `minV(size, vecA, vecB)` |
| Function | `maxS(size, vecA)` |
| Function | `minS(size, vecA)` |
| Function | `clamp(size, value, minS, maxS)` |
| Function | `transform(size, position, matrix)` |
| Class | `Vector` |
| Method | `Vector.__init__(self, size, data=None)` |
| Method | `Vector.__repr__(self)` |
| Method | `Vector.__add__(self, other)` |
| Method | `Vector.__iadd__(self, other)` |
| Method | `Vector.__sub__(self, other)` |
| Method | `Vector.__isub__(self, other)` |
| Method | `Vector.__mul__(self, scalar)` |
| Method | `Vector.__imul__(self, scalar)` |
| Method | `Vector.__div__(self, scalar)` |
| Method | `Vector.__truediv__(self, scalar)` |
| Method | `Vector.__idiv__(self, scalar)` |
| Method | `Vector.__itruediv__(self, scalar)` |
| Method | `Vector.__eq__(self, vecB)` |
| Method | `Vector.__ne__(self, vecB)` |
| Method | `Vector.__neg__(self)` |
| Method | `Vector.clone(self)` |
| Method | `Vector.one(self)` |
| Method | `Vector.zero(self)` |
| Method | `Vector.negate(self)` |
| Method | `Vector.maxV(self, vecB)` |
| Method | `Vector.maxS(self)` |
| Method | `Vector.minV(self, vecB)` |
| Method | `Vector.minS(self)` |
| Method | `Vector.magnitude(self)` |
| Method | `Vector.clamp(self, size, value, minS, maxS)` |
| Method | `Vector.i_clamp(self, size, value, minS, maxS)` |
| Method | `Vector.i_normalize(self)` |
| Method | `Vector.normalize(self)` |
| Method | `Vector.dot(self, vecB)` |
| Method | `Vector.isInSameDirection(self, otherVec)` |
| Method | `Vector.isInOppositeDirection(self, otherVec)` |
| Method | `Vector.barycentric(self, a, b, c)` |
| Method | `Vector.transform(self, position, matrix)` |
| Method | `Vector.i_transform(self, position, matrix)` |
| Method | `Vector.xy(self)` |
| Method | `Vector.yz(self)` |
| Method | `Vector.xz(self)` |
| Method | `Vector.xw(self)` |
| Method | `Vector.yw(self)` |
| Method | `Vector.zw(self)` |
| Method | `Vector.xyw(self)` |
| Method | `Vector.yzw(self)` |
| Method | `Vector.xzw(self)` |
| Method | `Vector.xyz(self)` |
| Method | `Vector.right(self)` |
| Method | `Vector.left(self)` |
| Method | `Vector.front(self)` |
| Method | `Vector.back(self)` |
| Method | `Vector.up(self)` |
| Method | `Vector.down(self)` |
