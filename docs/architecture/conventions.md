# Mathematical and ownership conventions

Current reference: master `3e714fe949b7a6b7724d5c0da3395ee92483265f`, including
the final SLERP repair (PR #59). These describe implemented contracts checked against source and
regressions, not a claim that all numeric inputs are validated. Historical
[audit conventions](../../audit/CONVENTIONS.md) contain baseline defects followed
by corrections; read their later decisions when comparing old behavior. This
page summarizes current behavior. API spellings/signatures remain unchanged.

## Vectors, dimensions and exact comparison

`Vector(size, data=None)` stores ordered components in `.vector`, with declared
dimension `.size`; 2/3/4 mean XY/XYZ/XYZW. Generic component arithmetic also works
for other well-formed sizes. Cross products and directional conveniences are 3D;
perpendicular helpers are 2D. `toAngle`, `lperp`, `rperp` take raw indexable
coordinate sequences, not a subscriptable Vector wrapper.

Equal dimensions compare every component exactly, with no tolerance. Different
dimensions compare unequal; two empty Vectors compare equal. `!=` complements
supported Vector equality. Other operand types receive `NotImplemented` from
the special methods. Arithmetic/dot expect matching well-formed storage; mixed
dimension error behavior is not standardized merely because equality is defined.
`isInSameDirection`/`isInOppositeDirection` test the sign of dot, not collinearity.

`lerp` and `scalarLerp` implement `a+t*(b-a)` without clamping. Reflection requires
a unit normal. Refraction uses eta=n1/n2 (incident/transmitted refractive index),
unit incident/normal vectors of matching dimension and a normal toward the
incident medium, opposing incidence. It does not normalize or flip inputs.
Total internal reflection returns a fresh zero Vector, not a reflected ray.

## Matrices and row-vector application

`.matrix[row][column]` is nested row-major storage. Products are the ordinary
matrix product `C[i][j]=sum(A[i][k]*B[k][j])`. Despite syntax `M*v`, the Vector
kernel computes the row-vector product **v M**. Mathematically v(A B) applies A
then B; in wrapper syntax `(A*B)*v == B*(A*v)` for matching dimensions.

```python
from gem.matrix import Matrix
from gem.vector import Vector

m = Matrix(2, [[1, 2], [3, 4]])
assert (m * Vector(2, [5, 6])).vector == [23, 34]
assert Matrix(2).matrix == [[1.0, 0.0], [0.0, 1.0]]
```

General Matrix*Vector requires matching dimensions; it does not promote Vector3
to Vector4. The implementation does not uniformly validate mismatches before
indexing. Matrix*Matrix explicitly rejects different declared sizes. Determinants
and inverses dispatch for sizes 2/3/4, not arbitrary-size linear algebra.

Homogeneous translation occupies the final row. Position `[x,y,z,1]` receives
translation; direction `[x,y,z,0]` does not. Vector4 offsets to Matrix4 translation
use their first three components. Matrix3/Vector2 translation is affine 2D;
legacy `translate3` overwrites the last row of a 3x3 matrix and is **not general
3D translation**. Matrix2 translation is unsupported.

`vector.transform(size, position, matrix)` takes raw position and nested matrix
lists. Same-size inputs use the ordinary row product. With an (N+1)x(N+1)
matrix it locally supplies w=1 and returns N components; this is intended for
affine positions. Explicit homogeneous inputs compute output w normally. There
is no perspective divide: implicit projective results are homogeneous numerators,
not perspective-correct Cartesian coordinates. Wrapper `transform` returns a new
Vector; `i_transform` changes only the receiver.

`rotate2(point, theta)` is a 3x3 homogeneous pivot rotation in degrees, positive
counterclockwise: `p'=pivot+(p-pivot)R`. Matrix2's wrapper rotation uses its linear
2x2 block (origin-only), while Matrix3/4 wrappers dispatch to 3D axis-angle
rotation. Shear XY means `Z'=Z+xX+yY`; YZ means `X'=X+yY+zZ`; XZ means
`Y'=Y+xX+zZ`. Other coordinates remain unchanged. Size-3 forms act on XYZ;
size-4 forms preserve homogeneous w. Wrapper transforms postmultiply the receiver.

## Orientation and angles

Cross products use right-handed X cross Y = Z. Positive +Z rotations move +X
toward +Y. `lookAt` and projection cameras view toward negative Z with a
right-handed frame. These conventions do not enforce a universal application
world frame. Vector `front()` is -Z, but identity Quaternion `getForward()` is
+Z. Neither is silently redefined.

| Existing API | Unit |
|---|---|
| `rotate2`, `rotate3`, `rotate4`, `Matrix.rotate`/`i_rotate` | Degrees |
| `rotate_origin2` | Radians |
| `perspective`, `perspectiveX` FOV | Degrees; vertical/horizontal respectively, aspect=width/height |
| `quat_from_axis_angle`, `quat_rotate`, `quat_rotate_from_axis_angle` | Degrees |
| `quat_rotate_x_from_angle`, `quat_rotate_y_from_angle`, `quat_rotate_z_from_angle` | Radians |
| `vector.toAngle`; SH theta/phi | Radians |

`radiansToDegrees(degrees)` multiplies by 180/pi;
`degreesToRadians(radians)` multiplies by pi/180. Their misleading historical
parameter names are retained; function names describe actual conversion.

## Quaternion algebra, domains and interpolation

Components are `[w,x,y,z]`, identity `[1,0,0,0]`, with Hamilton multiplication.
Unit q and -q represent the same rotation. q1*q2 applies q2 first; equivalent
row-vector rotation matrices compose in reversed order. Vector rotation is
`q*(0,v)*conjugate(q)`, returning Vector3. q*Vector alone returns a Quaternion
Hamilton product, not a rotated Vector. Nonunit sandwich inputs scale the result
by the squared norm; rotation APIs do not implicitly normalize q.

Use `quat_from_axis_angle` with a nonzero Vector3/list axis; it normalizes temporary
storage. Legacy `quat_rotate_from_axis_angle` rotates the normalized axis about
itself and returns an approximately pure Quaternion `[0,axis]`, not an orientation
constructor. `quat_rotate` instead takes raw unit axis coordinates without
normalizing them. X/Y/Z constructors return lists, not Quaternion wrappers.

```python
from gem.quaternion import quat_from_axis_angle, quat_rotate_vector
from gem.vector import Vector

axis = [0, 0, 2]
q = quat_from_axis_angle(axis, 90)
v = quat_rotate_vector(q, Vector(3, [1, 0, 0]))
assert axis == [0, 0, 2]
assert abs(v.vector[0]) < 1e-14 and abs(v.vector[1]-1) < 1e-14
```

Conversions assume unit quaternions/proper rotation matrices. `toMatrix` returns
Matrix4; `quat_from_matrix` reads a Matrix wrapper's upper-left 3x3 block and
returns Quaternion. There is no automatic orthogonalization, sign canonicalization
or general matrix validation. Quaternion inverse is conjugate/norm² for nonzero
ordinary inputs; its direct squared sum is not an extreme-scale stable inverse.

Powers/logarithms support unit inputs without implicit normalization or a new
unit tolerance check; general nonunit behavior is unsupported, not uniformly
rejected. Powers use a principal angle atan2(|imaginary|,w) in [0,pi] with finite
real exponents and fresh Quaternion output. q^0=identity and q^1=q on that domain.
Log returns a fresh list `[0,axis*angle]`. Exact zero inputs raise ValueError.
Negative identity supports integer powers by parity; fractional powers/log reject
its unspecified axis. Near-zero imaginary direction is retained without a cutoff.

LERP is component-linear, not normalized. Accurate `quat_slerp`/`Quaternion.slerp`
uses shortest-path sign correction for unit inputs, preserving input storage.
It neither normalizes inputs nor clamps t. Legacy `slerp_no_invert` remains
sign-sensitive, with linear approximation for dot outside (-.95,.95); unit norm
is not guaranteed, and an antipodal midpoint can be zero.

Legacy three-control SQUAD is nested no-invert interpolation: start=quat0,
end=quat2, control=quat1. Conventional `squad4(q0,q1,s0,s1,t)` separately uses
accurate SLERP between endpoints and between SQUAD controls, then blends with
2t(1-t). Controls are not simply neighbouring keyframes. Unit-output regression
tolerance is 1e-12 for supported unit inputs in [0,1], not a new validation rule.
No quaternion exponential, cross product or automatic control generation exists.

## Normalization, inversion and floating-point limits

Stable hypot-style norms and scaled normalization avoid avoidable intermediate
overflow/underflow for finite supported inputs. Direct zero Vector normalization
returns zero; zero Quaternion normalization returns identity. Returning results
are fresh; i-prefixed variants return the same receiver. These fallbacks do not
extend to geometric degeneracy: zero axes, degenerate lookAt frames, ray directions
and plane normals retain exact-zero ZeroDivisionError guards, without epsilons.

Matrix3/4 inverses retain power-of-two row scaling and an exact represented-binary64
coefficient singularity check before floating cofactor evaluation. Singular matrices
raise ZeroDivisionError; no arbitrary near-singular cutoff exists. Extreme uniform
scales around 1e-300 to 1e300 and selected mixed exponents have independent tests.
This does not guarantee accuracy for severely ill-conditioned matrices, floating
cofactor cancellation or all mixed-scale cases. Matrix2 inversion and public
determinants lack the same stabilization.

Ordinary Python float arithmetic is binary64. A true finite-input norm or inverse
coefficient outside its range may become signed infinity. Dot products, cross
products and all other helpers do not inherit blanket extreme-value guarantees.
No universal NaN/Infinity input policy is established: legacy paths retain their
arithmetic, while newer SH APIs explicitly reject nonfinite data. c_matrix is
binary32 and can overflow/round values that remain representable in Python.

## Ownership and mutation

Vector/Quaternion constructors retain supplied component storage; Matrix retains
supplied rows. Shape/type validation is limited. Returning arithmetic, normalize,
clamp, quaternion conversion and supported transform results allocate independent
output storage. i-prefixed methods and augmented assignment change receivers,
usually replacing component lists rather than updating external references to old
lists. Vector `zero()`/`one()` also mutate despite lacking the i prefix.

Returning clamp preserves value/bound lists. In-place clamp replaces receiver
storage without changing a separate caller list or another wrapper sharing old
storage. Direct field edits remain caller-managed; Matrix ctypes buffers and Plane
normal/coefficient snapshots do not track arbitrary edits automatically.

Plane coefficients satisfy `a*x+b*y+c*z+d=0`, with `.normal=[a,b,c]` at the
same scale. Three-point construction uses a unit normal and d=-dot(n,point).
Normalization scales all four coefficients together. `dot` takes Vector4,
including the supplied w; signed distance applies to a position with w=1 and a
unit normal. Newell polygon normals wrap vertices; bestFitD is signed D=mean(n·p),
converted to coefficient d=-D. Nonplanar polygons yield an approximation.

Ray construction retains start/direction Vectors, records the original direction
length as distance and normalizes the caller direction in place. Duplication deep
copies all Vector fields without rerunning construction. Rigid transforms mutate
the ray, replace start/direction, preserve distance and leave `.end` unchanged.
`.end` is intersection placeholder/state, not an automatically derived geometric
endpoint. It has no validity flag; even a zero value could be a hit at the origin.
Matrix4 translation promotes local positions with w=1 and directions with w=0.

## Projection and legacy viewport

`project` takes explicit Vector4 and Matrix4 wrappers/raw nested lists, including
mixed forms. Compose modelview*projection, divide by clip w, then map NDC to a
lower-left viewport `[x,y,width,height]`, Y upward. OpenGL NDC z in [-1,1]
becomes window z=(z+1)/2 in [0,1] for visible points; no clamping occurs. Near/far
camera depths -near/-far map to 0/1. The result is fresh Vector3.

`unproject` reverses viewport/depth mapping, multiplies by the inverse combined
matrix and divides by output w. Singular inversion raises ZeroDivisionError.
Project zero clip w raises ZeroDivisionError; unproject zero output w retains a
zero Vector3 sentinel. Invalid viewports/frusta and near-zero thresholds remain
separate policies, not a universal input-validation contract.

```python
from gem.matrix import Matrix, project, unproject
from gem.vector import Vector

viewport = [10, 20, 100, 200]
window = project(Vector(4, [0, 0, 0, 1]), Matrix(4), Matrix(4), viewport)
assert window.vector == [60.0, 120.0, 0.5]
assert unproject(*window.vector, Matrix(4), Matrix(4), viewport).vector == [0.0, 0.0, 0.0]
```

`common.getViewPort` is separate: normalize the whole Vector, scale XY by
width/height, then add original XY as offsets. It returns a four-element list,
preserves inputs and rejects zero with ZeroDivisionError. It is neither glViewport
nor conventional project/NDC mapping; Z/W affect its normalization.

## Curves, polynomials and sampling

Bezier evaluation uses Bernstein polynomials for scalars or Vector controls;
t is not clamped. Cubic paths use 3k+1 controls and retain explicit supplied
control lists. Sampling accepts finite scalars or uniform control representations
of Vector2/3 controls. `minimum_sqr_distance` is a positive finite **squared
coordinate-distance** tolerance; sampling uses its square root. Midpoint de
Casteljau subdivision checks maximum interior-control distance to the endpoint
**segment**, with clamped chord projection to retain collinear overshoot. Depth
is capped at 16; best-available output may exceed tolerance. No unconditional
approximation-error guarantee is claimed.

Standalone samples include both endpoints in increasing t order. `getDrawingPoints`
returns per-segment lists with shared boundary points omitted after the first.
`interpolate` appends controls intentionally; `samplePoints` rebuilds generated
controls using ordered source vertices and min/max squared-distance thinning
heuristics, not spacing guarantees. Both return None and no-op below two sources.

Legendre is unnormalized P_l^m with Condon–Shortley phase. Associated inputs use
integer 0<=m<=l and x in [-1,1]; m=0 ordinary polynomials also evaluate outside
that interval. `run()` preserves scratch state; explicit historical helpers update
it. Invalid/extreme-order behavior is not newly standardized.

## Spherical harmonics and radiometry

Real orthonormal SH includes Condon–Shortley phase, cosine for positive m, sine
for negative m. Theta is polar angle from +Z; phi is azimuth +X toward +Y, in
radians. Index=l(l+1)+m, complete coefficient count=bands². Canonical L2 order
has polynomial signs `[1,-Y,Z,-X,XY,-YZ,3Z²-1,-XZ,X²-Y²]` with orthonormal scales.
RGB projection/reconstruction returns fresh RGB lists; analytical rotation also
accepts scalar arrays of length 1/4/9. Scalar reconstruction is not a current API.

GenerateSamples uses global-RNG jittered strata and N=(sqrtNumSamples)² directions;
reset the seed for repeatability. Omitted projection weights mean 4pi/N, valid
for uniform sphere samples. Explicit weights are solid angles in steradians,
finite/nonnegative and not forced to sum to 4pi. Reconstruction requires supplied
unit XYZ/Vector3 directions, without implicit normalization or convolution.

Angular-disk RGB pixels use centers u=2(col+.5)/width-1, v=1-2(row+.5)/height,
r=hypot(u,v), theta=pi*r and phi=atan2(v,u); ignore r>1. Right/up/center correspond
to +X/+Y/+Z. Rectangular images stretch the disk. Weights are
4pi²/(width*height)*sinc(theta), without renormalization. This is approximate
quadrature, not a mirrored-ball photograph mapping. Latitude-longitude/RGBE
adapters exist only in the [HDR example](../../examples/hdr_sh/README.md), with
separate mapping and solid-angle rules.

The legacy probe class stores nine **radiance** coefficients in a rounded,
positive-X/Y basis. Explicit `legacy_to_canonical` changes odd-order signs and
scale factors; an import-path change alone does not convert coefficients.
`calculateCoefficients` rebuilds, while `updateCoefficients` accumulates.

Active analytical rotation means f_rotated(d)=f_original(R^-1 d): +90 degrees
about Z moves a +X light feature toward +Y. No directional resampling occurs.
Orientation must be finite Quaternion with norm deviation <=1e-12; a temporary
copy is normalized, larger deviations/zero/nonfinite inputs raise ValueError.
This local validation does not redefine all Quaternion APIs.

Diffuse convolution returns separate irradiance coefficients with degree factors
pi, 2pi/3, pi/4, supporting at most three bands. Apply it once. Reconstruct these
directly; Lambertian reflected radiance additionally equals albedo*irradiance/pi.
Lighting stays linear until display conversion; low-order SH truncation can ring
and yield negative reconstructions without changing coefficients.

```python
import math
from gem.spherical_harmonics import convolve_diffuse, reconstruct

# Only L0: orthonormal Y00=1/sqrt(4*pi); constant RGB radiance [1,2,3].
radiance = [[math.sqrt(4*math.pi)*c for c in [1, 2, 3]]]
irradiance = reconstruct(convolve_diffuse(radiance), [0, 0, 1])
assert all(abs(value-math.pi*c) < 1e-12 for value, c in zip(irradiance, [1, 2, 3]))
```

## Evidence and unresolved policies

Representative regressions: [vector/viewport](../../tests/test_final_vector_contracts.py),
[transforms](../../tests/test_transformations.py), [pivot/shear](../../tests/test_pivot_shear.py),
[projection](../../tests/test_projection.py), [quaternions](../../tests/test_quaternion_contracts.py),
[numerical limits](../../tests/test_numerical_robustness.py), [planes](../../tests/test_planes.py),
[rays](../../tests/test_rays.py), [Bezier sampling](../../tests/test_bezier_sampling.py),
[Legendre](../../tests/test_legendre.py), [SH rotation](../../tests/test_sh_rotation.py)
and [HDR reference](../../tests/test_hdr_sh_example.py). Their full-suite success
does not settle wider shape/scalar/error policies or Ray hit-state design.
See [open decisions](../development/roadmap.md#open-decisions).
