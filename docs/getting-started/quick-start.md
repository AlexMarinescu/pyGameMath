# First mathematical operations

Install [current repository code](installation.md), then run these independent
Python snippets using that environment. Each example is checked against an
isolated installed wheel, with no source-tree gem imports, renderer or GPU.
Examples use simple cross-version syntax; only the documented modern reference
interpreter was executed. Python 2.7 compatibility does not follow from syntax.

## Vectors and ownership

```python
from gem.vector import Vector, cross

values = [3.0, 4.0, 0.0]
v = Vector(3, values)
unit = v.normalize()
assert v.vector is values            # construction retains supplied storage
assert values == [3.0, 4.0, 0.0]     # returning normalization preserves it
assert unit.vector is not values
print(v.magnitude())                 # 5.0
print(unit.vector)                   # [0.6, 0.8, 0.0]
print(v.dot(Vector(3, [1, 0, 0])))    # 3.0
print(cross(Vector(3, [1, 0, 0]), Vector(3, [0, 1, 0])).vector)
# [0, 0, 1]
```

Use matching dimensions for arithmetic. Vector equality is exact, including
dimension semantics. `normalize()` returns a new Vector; `i_normalize()` changes
the receiver and returns it. Direct zero-Vector normalization yields zero, but
geometric callers can still reject a degenerate direction. Constructor aliasing
and receiver mutation are deliberate distinctions, not a universal copy policy.

## Matrix transformation order

Matrices store nested rows. `M * v` implements the mathematical row-vector product
vM. A product translation*rotation applies translation **then** rotation.
This example starts at +X, translates two units along +X and rotates +90 degrees
about +Z, producing +3Y:

```python
from gem.matrix import Matrix
from gem.vector import Vector

translation = Matrix(4).translate(Vector(3, [2, 0, 0]))
rotation = Matrix(4).rotate(Vector(3, [0, 0, 1]), 90)  # degrees
position = Vector(4, [1, 0, 0, 1])
after = (translation * rotation) * position
assert abs(after.vector[0]) < 1e-14
assert abs(after.vector[1] - 3.0) < 1e-14
assert after.vector[2:] == [0.0, 1.0]
print([round(value, 6) for value in after.vector])    # [0.0, 3.0, 0.0, 1.0]

other_order = (rotation * translation) * position
assert abs(other_order.vector[0] - 2.0) < 1e-14
assert abs(other_order.vector[1] - 1.0) < 1e-14
assert position.vector == [1, 0, 0, 1]
```

Use explicit w=1 positions and w=0 directions with general Matrix4 multiplication;
it does not implicitly promote Vector3. Translation is stored in the final row.
The separate Vector transform helper has local affine promotion and no automatic
perspective divide. Read [projection and transform conventions](../architecture/conventions.md)
before combining camera/projective matrices.

## Matrix inversion

```python
from gem.matrix import Matrix

matrix = Matrix(2, [[4.0, 7.0], [2.0, 6.0]])
inverse = matrix.inverse()
expected = [[0.6, -0.7], [-0.2, 0.4]]
assert all(abs(inverse.matrix[i][j] - expected[i][j]) < 1e-14
           for i in range(2) for j in range(2))
for identity in (matrix * inverse, inverse * matrix):
    assert all(abs(identity.matrix[i][j] - (1 if i == j else 0)) < 1e-14
               for i in range(2) for j in range(2))
assert matrix.matrix == [[4.0, 7.0], [2.0, 6.0]]
```

Inverse/determinant dispatch supports dimensions 2/3/4. Singular inversion raises
ZeroDivisionError. Matrix3/4 inverses have scaled extreme-value handling; that
does not establish identical robustness for public determinants or Matrix2 inverse.
Use `i_inverse()` only when receiver mutation is intended.

## Quaternion orientation

```python
import math
from gem.quaternion import quat_from_axis_angle, quat_rotate_vector
from gem.vector import Vector

axis = [0, 0, 2]
q = quat_from_axis_angle(axis, 90)      # degrees, temporary axis normalization
assert axis == [0, 0, 2]
assert abs(q.data[0] - math.sqrt(0.5)) < 1e-14  # [w,x,y,z]
turned = quat_rotate_vector(q, Vector(3, [1, 0, 0]))
assert abs(turned.vector[0]) < 1e-14 and abs(turned.vector[1] - 1.0) < 1e-14
print([round(value, 6) for value in turned.vector])   # [0.0, 1.0, 0.0]
```

The recommended constructor is `quat_from_axis_angle`, not the similarly named
legacy `quat_rotate_from_axis_angle`. Rotation uses a unit Quaternion sandwich;
q*Vector alone returns a Quaternion product, not rotated Vector. Other helpers
use radians, so units are specified per API. See the [quaternion guide](../QUATERNIONS.md)
for interpolation, return types and the legacy forward-axis distinction.

## Bezier evaluation and sampling

```python
from gem.bezier import BezierPath, cubicBezierPoint, quadraticBezierPoint
from gem.vector import Vector

controls = [Vector(2, [0, 0]), Vector(2, [1, 2]),
            Vector(2, [2, 2]), Vector(2, [3, 0])]
assert cubicBezierPoint(0.5, *controls).vector == [1.5, 1.5]
assert quadraticBezierPoint(0.5, 0.0, 2.0, 4.0) == 2.0
path = BezierPath()
path.setControlPoints(controls)
path.minimum_sqr_distance = 0.0001     # coordinate-distance tolerance 0.01
points = path.findDrawingPoints(0)
assert points[0].vector == [0, 0] and points[-1].vector == [3, 0]
assert controls[1].vector == [1, 2]
```

Parameters are not clamped. Cubic paths use 3k+1 controls and retain supplied
control lists. Adaptive sampling supports finite scalar/Vector2/Vector3 controls,
with squared-distance tolerance and a depth-16 cap; capped output can exceed
tolerance. `getDrawingPoints()` returns nested per-segment lists rather than a
single flat path. Source-point thinning is a different builder operation.

## Constant-radiance SH reference

This exact L0 reference avoids introducing sampling error into the first example.
Canonical Y00=1/sqrt(4pi), so constant RGB radiance [1,2,3] has the coefficient
sqrt(4pi)*[1,2,3]. The pipeline convolves once to get irradiance pi*[1,2,3]:

```python
import math
from gem.spherical_harmonics import convolve_diffuse, reconstruct

radiance = [[math.sqrt(4 * math.pi) * channel for channel in [1, 2, 3]]]
irradiance_coefficients = convolve_diffuse(radiance)
irradiance = reconstruct(irradiance_coefficients, [0, 0, 1])  # unit normal
assert all(abs(value - math.pi * channel) < 1e-12
           for value, channel in zip(irradiance, [1, 2, 3]))
albedo = [0.5, 0.5, 0.5]
reflected = [a * value / math.pi for a, value in zip(albedo, irradiance)]
assert all(abs(value - expected) < 1e-12
           for value, expected in zip(reflected, [0.5, 1.0, 1.5]))
assert radiance != irradiance_coefficients
```

Coefficients are canonical coefficient-by-RGB lists; full bands have bands²
entries, index l(l+1)+m. Radiance, irradiance and reflected Lambertian radiance
are distinct. Reconstruction does not normalize the supplied direction or convolve
again. All calculations above remain linear, without tone mapping or display encoding.
Historical probe coefficients require explicit `legacy_to_canonical`; changing
import paths alone does not change their basis.

See [the SH guide](../SPHERICAL_HARMONICS.md) for projection, sampling, basis signs
and analytical rotation, and the [HDR/SH example](../../examples/hdr_sh/README.md)
for procedural environments, shader-ready exports and CPU visualization.

## ctypes matrix export

```python
import ctypes
from gem.matrix import Matrix
from gem.vector import Vector

matrix = Matrix(4).translate(Vector(3, [2, -3, 4]))
assert list(matrix.c_matrix[3]) == [2.0, -3.0, 4.0, 1.0]
pointer = ctypes.cast(matrix.c_matrix, ctypes.POINTER(ctypes.c_float))
assert [pointer[index] for index in range(12, 16)] == [2.0, -3.0, 4.0, 1.0]
matrix.i_translate(Vector(3, [1, 0, 0]))
assert list(matrix.c_matrix[3]) == [3.0, -3.0, 4.0, 1.0]
```

The export is a row-major float32 snapshot, not a live view of Python rows.
Supported in-place operations refresh it; direct row edits do not. Keep the
owning array alive while a foreign consumer uses its pointer, and reacquire the
pointer after operations replace the snapshot. Match the consumer's layout and
transpose policy explicitly; this example verifies memory values, not an OpenGL upload.

## Next steps

Use [the API inventory](../architecture/api-inventory.md) to find actual signatures
and documented constraints, [conventions](../architecture/conventions.md) for
composition/units and [compatibility](../architecture/compatibility.md) for support
evidence. The [canonical roadmap](../../ROADMAP.md) describes future features;
none is implied by these examples. Return to [getting started](README.md) or
[installation](installation.md) for the beginner journey.
