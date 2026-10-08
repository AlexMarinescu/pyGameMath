# Ray geometry

`gem.ray.Ray(startVector,dirVector)` retains caller-owned Vector references.
Construction records the original direction magnitude as `distance` and
normalizes the caller's direction in place. Ordinary geometry uses Vector3
positions and nonzero directions. The existing zero-direction normalization
error and wider malformed/nonfinite policies are unchanged.

The geometric position at distance d is `start + dir*d`. The stored
`.end` field is separate **intersection placeholder/state**, initialized to
a zero Vector3. It is not automatically set to `start + dir*distance`.

## Copying and transforms

`duplicate()` returns a new Ray with independent start, direction and end
Vectors/component lists and the exact stored distance. It does not invoke
the constructor or re-normalize stored direction. Modifying either copy's
Vector fields cannot modify the other's storage.

`rotateUsingQuaternion(quat1)` rotates start and direction about the
coordinate origin using the established Hamilton sandwich
`q*(0,v)*conjugate(q)`. Inputs must be unit rotation quaternions in [w,x,y,z]
order. Start/direction remain Vector3 values, direction is normalized, and
distance is retained. The supplied quaternion and previously referenced
input Vectors are preserved.

`roateUsingMatrix(matrix)` retains its historical spelling and Matrix3
rotation behavior, including rotation about the coordinate origin and
direction normalization. Proper rotation matrices use the established
row-vector convention; distance is retained. No new method name or
Matrix4 rotation overload is added.

`translate(matrix)` supports pure Matrix4 translations on Vector3 geometry
using local homogeneous positions `[x,y,z,1]` and directions `[x,y,z,0]`.
It returns spatial Vector3 fields, changes origin, and preserves direction
and distance. No perspective divide or general Matrix*Vector promotion is
introduced. Existing matching-dimension multiplication remains the fallback;
this does not establish scale/shear/projective ray semantics.

All transform methods mutate the receiver and retain their None returns.
They replace transformed start/direction fields rather than mutate the
previously referenced Vectors. Duplication and transforms preserve matrix
or quaternion inputs. Public field edits remain caller-managed.

```python
from gem.ray import Ray
from gem.vector import Vector
from gem.matrix import Matrix
from gem.quaternion import quat_from_axis_angle

origin = Vector(3, [1, 2, 3])
displacement = Vector(3, [0, 0, 5])
r = Ray(origin, displacement)  # distance=5; displacement becomes [0,0,1]
copy = r.duplicate()
r.translate(Matrix(4).translate(Vector(3, [2, -3, 4])))
# r.start=[3,-1,7], r.dir=[0,0,1], distance=5; origin remains [1,2,3].
r.rotateUsingQuaternion(quat_from_axis_angle([0, 0, 1], 90))
# r.start is approximately [1,3,7]; copy retains its original geometry.
# r.end remains the zero intersection placeholder throughout.
```

## Intersection state remains unresolved

Duplication copies `.end` exactly. Construction and all transforms retain
its historical semantics: transforms leave the existing object, component
storage and values unchanged, whether zero or nonzero. A zero Vector might
mean an unset placeholder or an actual hit at the origin; there is no
validity flag. A nonzero value alone likewise does not establish a valid hit.
No value-based inference or new intersection behavior is introduced.

A future API decision must define validity, hit ownership and whether hit
coordinates transform with geometry. Applications using `.end` as a hit
position must currently manage that state explicitly. Derived geometric
tips can be computed independently and must not be confused with `.end`.
Scale, shear, projective transforms and broader numerical policies remain
outside this contract.
