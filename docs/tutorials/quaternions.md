# Object orientation and quaternion interpolation

**Intermediate.** Prerequisites: [vectors](vectors.md), [transforms](transformations.md)
and sine/cosine. Learn to rotate a direction, compose rotations and interpolate
an object/camera orientation without confusing component blends with unit rotations.

## What a rotation quaternion represents

Quaternion storage is [w,x,y,z]. For a unit axis n and angle theta,
q=[cos(theta/2),n*sin(theta/2)]. Unit q rotates v by the Hamilton sandwich
q*[0,v]*conjugate(q). q and -q represent the same orientation; component equality
is not the correct rotation comparison. Four components with a unit constraint
avoid choosing an Euler-angle sequence, but they do not remove interpolation,
sign or branch choices.

`quat_from_axis_angle` takes degrees and normalizes temporary Vector/list axis
storage. Raw X/Y/Z angle helpers take radians and return lists, not Quaternion
wrappers. The legacy `quat_rotate_from_axis_angle` rotates its normalized axis
about itself and returns a pure Quaternion; it is not the constructor to use here.

Hamilton qZ*qX applies qX first; equivalent row matrices compose RX*RZ. Identity
Quaternion.getForward is +Z, while Vector.front is -Z. Set an application's local
forward axis explicitly rather than assuming those conveniences share a convention.

## Composition and a halfway orientation

```python
# Axis-angle construction uses degrees; interpolation here uses unit inputs.
import math
from gem.quaternion import Quaternion, quat_from_axis_angle, quat_from_matrix, quat_rotate_vector
from gem.vector import Vector

qz = quat_from_axis_angle([0, 0, 1], 90)
qx = quat_from_axis_angle([1, 0, 0], 90)
y = Vector(3, [0, 1, 0])
a = quat_rotate_vector(qz * qx, y)
b = quat_rotate_vector(qx * qz, y)
assert all(abs(x - e) < 1e-14 for x, e in zip(a.vector, [0, 0, 1]))
assert all(abs(x - e) < 1e-14 for x, e in zip(b.vector, [-1, 0, 0]))
mid = Quaternion().slerp(qz, 0.5)
turned = quat_rotate_vector(mid, Vector(3, [1, 0, 0]))
s = math.sqrt(0.5)
assert all(abs(x - e) < 1e-14 for x, e in zip(turned.vector, [s, s, 0]))
assert abs(mid.magnitude() - 1) < 1e-12
matrix_result = qz.toMatrix() * Vector(4, [1, 0, 0, 0])
assert all(abs(x - e) < 1e-14 for x, e in zip(matrix_result.vector, [0, 1, 0, 0]))
mapped = (qx.toMatrix() * qz.toMatrix()) * Vector(4, [0, 1, 0, 0])
assert all(abs(x - e) < 1e-14 for x, e in zip(mapped.vector, [0, 0, 1, 0]))
sign_equivalent = Quaternion().slerp(qz.negate(), 0.5)
assert abs(abs(sign_equivalent.dot(mid)) - 1) < 1e-14
back = quat_from_matrix(qz.toMatrix())
assert abs(abs(back.dot(qz)) - 1) < 1e-14
assert y.vector == [0, 1, 0]
print([round(x, 6) for x in turned.vector])
```

Output: `[0.707107, 0.707107, 0.0]`. For this single-axis interval, halfway
SLERP gives 45° and constant angular speed when t advances uniformly. Accurate
shortest-path SLERP corrects endpoint signs in temporary storage and resolves tiny
angles using difference/sum norms. It does not normalize inputs or clamp t.
Quaternion*Vector alone returns a Quaternion product, not a rotated Vector.

## Conventional SQUAD4 versus retained SQUAD

Let S be accurate shortest-path SLERP. Conventional SQUAD4 computes
S(S(q0,q1,t),S(s0,s1,t),2t(1-t)). Endpoints are q0,q1; s0,s1 are intermediate
SQUAD control quaternions, not arbitrary neighbouring keyframes. gem does not
generate those controls automatically. Here choose controls equal to endpoints
to obtain an independently checkable single-axis case.

Legacy three-control SQUAD instead uses sign-sensitive no-invert interpolation N:
N(N(q0,q2,t),N(q0,q1,t),2t(1-t)), with start q0, control q1, end q2.
It retains near-angle linear approximation and is not guaranteed unit length.

```python
import math
from gem.quaternion import Quaternion, quat_from_axis_angle, squad4

start = Quaternion()
end = quat_from_axis_angle([0, 0, 1], 90)
saved = end.data[:]
standard = squad4(start, end, start, end, 0.5)
legacy = start.squad(start, end, 0.5)
assert abs(standard.data[0] - math.cos(math.pi / 8)) < 1e-14
assert abs(standard.data[3] - math.sin(math.pi / 8)) < 1e-14  # 45-degree orientation
assert abs(legacy.data[3] - math.sin(math.pi / 16)) < 1e-14  # 22.5-degree orientation
assert abs(standard.magnitude() - 1) < 1e-12
assert standard is not start and standard.data is not start.data
assert end.data == saved
```

For the legacy midpoint, the first inner blend is 45°, the second 0°, and the
outer blend weight is 1/2, giving 22.5°. Its spherical branches are active for
these particular angles; that answer does not prove a general norm guarantee.
Neither algorithm alone guarantees constant angular speed on a multi-control path.

## Limits and applications

Use orientation interpolation for camera turns or object poses, keeping position
interpolation separate. Conversions assume proper rotation matrices/unit quaternions;
scale/shear is not automatically removed. Nonunit sandwich rotation scales vectors
by norm squared, and zero normalization's identity fallback does not turn a zero
axis into a valid axis-angle input. Wider invalid/nonfinite behavior is not uniform.

See [Quaternion API](../api/quaternion.md), [legacy distinctions](../QUATERNIONS.md),
[interpolation regressions](../../tests/test_quaternion_interpolation.py) and
[numerical accuracy](numerical.md). Continue with [motion paths](bezier.md),
[SH rotation](lighting.md) or the [index](index.md).

See the [visual example](../examples/gallery/quaternions.md) and its reproducible assets.
