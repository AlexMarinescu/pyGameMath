# Quaternion algebra, rotation and interpolation

Import `Quaternion` and the functions below from `gem.quaternion`.
[Source](../../gem/quaternion.py), [operations](../../tests/test_quaternion_operations.py),
[powers/logarithms](../../tests/test_quaternion_powers.py),
[interpolation](../../tests/test_quaternion_interpolation.py),
[conversion](../../tests/test_quaternion_matrix.py).

`.data` stores [w,x,y,z]; default identity is [1,0,0,0]. Supplied storage is retained
by reference. Quaternion has no component indexing, iteration, value equality,
unary-minus or power operator; access `.data`, use `negate()` and `pow(e)`.
No public Quaternion constants or direct ctypes export exist; `toMatrix()` provides
a Matrix4 with its own float32 snapshot.

Raw algebra kernels take four-element numeric sequences, with no uniform shape
validation, and return fresh lists except scalar dot/magnitude. Wrapper algebra
returns fresh Quaternion objects; in-place methods/augmented assignment replace
`.data` and return self, preserving separate external references to old lists.
Unsupported operator operands return NotImplemented (normally causing TypeError
if Python cannot dispatch elsewhere). There are no reflected scalar operators.

Hamilton product (w,u)(s,v)=(ws−u·v,wv+su+u×v). q and −q describe the same unit
rotation. q1*q2 applies q2 first; their row-vector matrices compose in reversed
order. Ordinary inverse is conjugate/norm²; it is not stabilized for extreme
squared-component overflow/underflow. Exact zero inverse raises ZeroDivisionError.
Stable norm and normalization are separate: zero normalization returns identity,
and finite scaled hypot avoids avoidable intermediate range errors. Nonfinite
paths retain legacy arithmetic, not a universal validation policy.

## Class and methods

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `Quaternion` | Hamilton component wrapper; use constructor and documented domains below. |
| `Quaternion.__init__(self, data=None)` | `data=None`: fresh identity list; supplied [w,x,y,z] storage retained without shape validation. Initialization returns None. |
| `Quaternion.__add__(self, other)` | `other`: Quaternion only; component addition; fresh Quaternion. Unsupported types return NotImplemented. |
| `Quaternion.__iadd__(self, other)` | `other`: Quaternion only; component addition; replace receiver data, return self. Unsupported types return NotImplemented. |
| `Quaternion.__sub__(self, other)` | `other`: Quaternion only; component subtraction; fresh Quaternion. Unsupported types return NotImplemented. |
| `Quaternion.__isub__(self, other)` | `other`: Quaternion only; component subtraction; replace receiver data, return self. Unsupported types return NotImplemented. |
| `Quaternion.__mul__(self, other)` | `other`: Quaternion → Hamilton product; Vector3 → Quaternion product with [0,v]; float → scalar product (int unsupported). Fresh Quaternion. Vector product alone is not rotated Vector; unsupported operands return NotImplemented. |
| `Quaternion.__imul__(self, other)` | `other`: Quaternion → Hamilton product; Vector3 → Quaternion product with [0,v]; float → scalar product (int unsupported). Replace data, return self. Vector product alone is not rotated Vector; unsupported operands return NotImplemented. |
| `Quaternion.__div__(self, other)` | `other`: float only; component division, fresh Quaternion. Unsupported types return NotImplemented; zero divisor raises ZeroDivisionError. |
| `Quaternion.__idiv__(self, other)` | `other`: float only; component division, replace data and return self. Unsupported types return NotImplemented; zero divisor raises ZeroDivisionError. |
| `Quaternion.__truediv__ = Quaternion.__div__` | Python 3 alias of `__div__`; retains signature, float-only acceptance, ownership and errors. |
| `Quaternion.__itruediv__ = Quaternion.__idiv__` | Python 3 alias of `__idiv__`; retains signature, float-only acceptance, ownership and errors. |
| `Quaternion.i_negate(self)` | Negate components, replace receiver data, return self. |
| `Quaternion.negate(self)` | Negate all four components; fresh Quaternion, preserving receiver. |
| `Quaternion.i_identity(self)` | Replace receiver data with identity; return self. |
| `Quaternion.identity(self)` | Fresh identity Quaternion, independent of receiver values. |
| `Quaternion.magnitude(self)` | Stable numeric norm of data; no mutation. |
| `Quaternion.dot(self, quat2)` | `quat2`: Quaternion; numeric four-component dot; unsupported type returns NotImplemented. |
| `Quaternion.i_normalize(self)` | Stable normalization replaces data, exact zero → identity; return self. |
| `Quaternion.normalize(self)` | Fresh normalized Quaternion; exact zero → identity; stable scaled norm. |
| `Quaternion.i_conjugate(self)` | Conjugate components, replace data, return self. |
| `Quaternion.conjugate(self)` | Fresh [w,−x,−y,−z] Quaternion. |
| `Quaternion.inverse(self)` | Fresh conjugate/norm² Quaternion; exact zero raises ZeroDivisionError, no extreme-range protection. |
| `Quaternion.pow(self, e)` | `e`: finite real exponent on unit inputs; fresh Quaternion; see principal-branch domain below. |
| `Quaternion.log(self)` | Unit input; fresh four-element logarithm list, not Quaternion; see branch restrictions below. |
| `Quaternion.lerp(self, quat1, time)` | `quat1`: Quaternion, `time`: unclamped numeric interpolation parameter; fresh component-linear Quaternion, not normalized. |
| `Quaternion.slerp(self, quat1, time)` | `quat1`: unit Quaternion, `time`: numeric; fresh accurate shortest-path Quaternion; no input normalization/clamping. |
| `Quaternion.slerp_no_invert(self, quat1, time)` | `quat1`: unit Quaternion, `time`: numeric; fresh sign-sensitive interpolation with near-angle linear approximation, not guaranteed unit. |
| `Quaternion.squad(self, quat1, quat2, time)` | `quat1`: intermediate control, `quat2`: endpoint, `time`: numeric; fresh legacy three-control blend (receiver starts). See formula below. |
| `Quaternion.toMatrix(self)` | Unit receiver prerequisite; fresh row-vector Matrix4 with synchronized ctypes; no implicit normalization. |
| `Quaternion.getForward(self)` | Fresh Vector3 from conjugate-sandwich rotation of [0,0,1]; unit receiver prerequisite; `getForward` uses +Z unlike Vector.front (−Z). |
| `Quaternion.getBack(self)` | Fresh Vector3 from conjugate-sandwich rotation of [0,0,-1]; unit receiver prerequisite; `getForward` uses +Z unlike Vector.front (−Z). |
| `Quaternion.getLeft(self)` | Fresh Vector3 from conjugate-sandwich rotation of [-1,0,0]; unit receiver prerequisite; `getForward` uses +Z unlike Vector.front (−Z). |
| `Quaternion.getRight(self)` | Fresh Vector3 from conjugate-sandwich rotation of [1,0,0]; unit receiver prerequisite; `getForward` uses +Z unlike Vector.front (−Z). |
| `Quaternion.getUp(self)` | Fresh Vector3 from conjugate-sandwich rotation of [0,1,0]; unit receiver prerequisite; `getForward` uses +Z unlike Vector.front (−Z). |
| `Quaternion.getDown(self)` | Fresh Vector3 from conjugate-sandwich rotation of [0,-1,0]; unit receiver prerequisite; `getForward` uses +Z unlike Vector.front (−Z). |

## Functions

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `quat_identity()` | Fresh [1.0,0.0,0.0,0.0] list; no operands. |
| `quat_add(quat, quat1)` | Raw `quat,quat1`: fresh component sum list. |
| `quat_sub(quat, quat1)` | Raw `quat,quat1`: fresh component difference list. |
| `quat_mul_quat(quat, quat1)` | Raw `quat,quat1`: fresh Hamilton product list, in operand order. |
| `quat_mul_vect(quat, vect)` | Raw four-component `quat`, XYZ `vect` sequence: fresh quaternion list quat*[0,vect]. |
| `quat_mul_float(quat, scalar)` | Raw `quat`, numeric `scalar`: fresh component product list; unlike wrapper, kernel does not explicitly restrict scalar to float. |
| `quat_div_float(quat, scalar)` | Raw `quat`, numeric `scalar`: fresh component quotient list; zero division raises ZeroDivisionError. |
| `quat_neg(quat)` | Raw `quat`: fresh negated four-component list. |
| `quat_dot(quat1, quat2)` | Raw `quat1,quat2`: scalar sum of four products, inputs preserved. |
| `quat_magnitude(quat)` | Raw `quat`: stable four-component norm scalar, no mutation; true unrepresentable length may be infinity. |
| `quat_normalize(quat)` | Raw `quat`: fresh stable normalized list; exact zero maps to identity. |
| `quat_conjugate(quat)` | Raw `quat`: fresh list [w,−x,−y,−z]. |
| `quat_inverse(quat)` | Raw `quat`: fresh conjugate/squared-norm list; exact zero raises ZeroDivisionError; extreme products may overflow/underflow. |
| `quat_from_axis_angle(axis, theta)` | `axis`: nonzero Vector3 or list; `theta`: degrees. Normalize a temporary axis, return fresh rotation Quaternion [cos(theta/2),axis*sin(theta/2)]. Preserve caller storage; zero axis ZeroDivisionError, unsupported axis type returns NotImplemented. |
| `quat_rotate(origin, axis, theta)` | `origin`: Vector3 to rotate about coordinate origin; `axis`: raw unit XYZ coordinates; `theta`: degrees. Fresh Vector3 from sandwich; axis is not normalized, so nonunit axes do not yield pure rotations. |
| `quat_rotate_x_from_angle(theta)` | `theta`: radians; fresh raw list [cos(theta/2),sin(theta/2),0,0], not a Quaternion wrapper. |
| `quat_rotate_y_from_angle(theta)` | `theta`: radians; fresh raw list [cos(theta/2),0,sin(theta/2),0]. |
| `quat_rotate_z_from_angle(theta)` | `theta`: radians; fresh raw list [cos(theta/2),0,0,sin(theta/2)]. |
| `quat_rotate_from_axis_angle(axis, theta)` | Legacy: nonzero Vector3/list `axis`, `theta` degrees; normalize temporary axis, rotate it about itself, return approximately pure Quaternion [0,axis]. Not an orientation constructor. Preserve inputs; zero axis ZeroDivisionError, unsupported type NotImplemented. |
| `quat_rotate_vector(quat, vec)` | `quat`: unit Quaternion, `vec`: Vector3; fresh Vector3 imaginary part of q*[0,v]*conjugate(q). Preserve both inputs; nonunit q scales by norm², no implicit normalization. |
| `quat_pow(quat, exp)` | `quat`: unit Quaternion, `exp`: finite real; fresh Quaternion [cos(exp*a),axis*sin(exp*a)], a=atan2(imaginary norm,w). See exact-zero and negative-identity exceptions below. |
| `quat_log(quat)` | `quat`: unit Quaternion; fresh list [0,axis*a], a=atan2(imaginary norm,w). Identity returns zeros; exact zero and negative identity raise ValueError; inputs preserved. |
| `quat_lerp(quat0, quat1, t)` | `quat0,quat1`: Quaternions; numeric `t`; fresh (1−t)*q0+t*q1, no normalization or clamping. |
| `quat_slerp(quat0, quat1, t)` | Unit `quat0,quat1`, numeric `t`; fresh accurate shortest-path spherical interpolation; negative dot flips only a temporary endpoint, no normalization or clamping. |
| `quat_slerp_no_invert(quat0, quat1, t)` | Unit `quat0,quat1`, numeric `t`; fresh sign-sensitive spherical blend for −0.95<dot<0.95, otherwise component-linear blend. No unit-output guarantee; antipodal midpoint can be zero. |
| `quat_squad(quat0, quat1, quat2, t)` | Unit `quat0` start, `quat2` end, `quat1` control, numeric `t`; legacy nested no-invert formula below. Fresh Quaternion, inputs preserved. |
| `squad4(q0, q1, s0, s1, t)` | Unit endpoints `q0,q1`, SQUAD controls `s0,s1`, `t` ordinarily [0,1]; conventional nested accurate shortest-path SLERP below. Fresh Quaternion; no auto-generation of controls or validation/clamping of t. |
| `quat_to_matrix(quat)` | Unit `quat`: Quaternion; fresh Matrix4, row-vector rotation and independent synchronized c_matrix. Input preserved; no normalization. |
| `quat_from_matrix(matrix)` | `matrix`: Matrix3/4 wrapper with proper rotation in upper-left 3×3; fresh Quaternion, largest-component branch. Preserve input; q/−q equivalent. No automatic orthogonalization, validation or sign canonicalization. |

## Conversion and boundary examples

For a unit [w,x,y,z], the row-vector rotation block is
[[1−2(y²+z²), 2(xy+wz), 2(xz−wy)],
 [2(xy−wz), 1−2(x²+z²), 2(yz+wx)],
 [2(xz+wy), 2(yz−wx), 1−2(x²+y²)]].
Matrix conversion selects the largest squared quaternion component from the
rotation trace/diagonal candidates, then recovers the other components using
sums/differences. This avoids a zero-trace branch failure for legitimate half-turns.
The result's sign is not unique; compare rotations or absolute unit dot, not
component equality.

```python
from gem.quaternion import Quaternion, quat_from_axis_angle, quat_from_matrix

q = quat_from_axis_angle([0, 1, 0], 180)
m = q.toMatrix()
back = quat_from_matrix(m)
assert abs(abs(q.dot(back)) - 1.0) < 1e-14
assert abs(m.matrix[0][0] + 1.0) < 1e-14
assert abs(m.matrix[2][2] + 1.0) < 1e-14
assert abs(float(m.c_matrix[0][0]) + 1.0) < 1e-6
assert Quaternion([0, 0, 0, 0]).normalize().data == [1.0, 0.0, 0.0, 0.0]
assert Quaternion().__mul__(2) is NotImplemented   # float-only scalar operation
assert Quaternion([-1, 0, 0, 0]).pow(3).data == [-1.0, 0.0, 0.0, 0.0]
try:
    Quaternion([-1, 0, 0, 0]).log()
except ValueError:
    pass
else:
    raise AssertionError('negative identity has no unique log axis')
```

## Powers, logarithms and interpolation domains

Powers/logarithms support unit quaternions without implicitly normalizing or
checking a new norm tolerance. General nonunit behavior is unsupported, not
uniformly rejected. Principal angle atan2(hypot(x,y,z),w) lies in [0,pi]; imaginary
components retain near-zero direction without an epsilon cutoff. Signs are not
canonicalized. q^0=identity and q^1=q on the supported domain. Exact zero inputs
raise ValueError. Negative identity [−1,0,0,0] supports integer powers by parity;
fractional powers and log reject the nonunique axis with ValueError. No quaternion
exponential API exists. Finite exponents are prerequisites, not a general validator.

Accurate SLERP uses the difference/sum norm angle to resolve tiny separations,
then sine weights; identical orientations use the continuous limit. For unit inputs
and t in [0,1], unit norm is tested within 1e-12; this is not implicit normalization.
Legacy no-invert SLERP retains sign-sensitive branches and linear approximation.

Legacy three-control formula, with N=no-invert SLERP:
a=N(q0,q2,t); b=N(q0,q1,t); result=N(a,b,2t(1−t)).
Conventional `squad4`, with S=accurate shortest-path SLERP:
a=S(q0,q1,t); b=S(s0,s1,t); result=S(a,b,2t(1−t)).
The latter controls are SQUAD controls, not simply adjacent keyframes. Neither
function generates Shoemake controls or changes the legacy signature.

```python
import math
from gem.quaternion import Quaternion, quat_from_axis_angle, quat_rotate_vector, quat_rotate_from_axis_angle, squad4
from gem.vector import Vector

axis = [0, 0, 2]
q = quat_from_axis_angle(axis, 90)
assert axis == [0, 0, 2]
v = quat_rotate_vector(q, Vector(3, [1, 0, 0]))
assert all(abs(a - b) < 1e-14 for a, b in zip(v.vector, [0, 1, 0]))
legacy = quat_rotate_from_axis_angle(axis, 90)
assert all(abs(a - b) < 1e-14 for a, b in zip(legacy.data, [0, 0, 0, 1]))
assert abs(q.log()[3] - math.pi / 4) < 1e-14
assert abs(q.pow(2).data[3] - 1) < 1e-14
mid = Quaternion().slerp(q, 0.5)
assert abs(mid.data[0] - math.cos(math.pi / 8)) < 1e-14
s = squad4(Quaternion(), q, Quaternion(), q, 0.5)
assert all(abs(a - b) < 1e-14 for a, b in zip(s.data, mid.data))
assert abs(s.magnitude() - 1) < 1e-12
nonunit = Quaternion([2, 1, 0, 0])
identity = nonunit * nonunit.inverse()
assert all(abs(a - b) < 1e-14 for a, b in zip(identity.data, [1, 0, 0, 0]))
```

See [historical quaternion guide](../QUATERNIONS.md), [matrix](matrix.md),
[decisions](decisions.md) and [API index](index.md).

See the [graphics gallery example](../examples/gallery/quaternions.md) for an executable visualization.
