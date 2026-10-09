# Ray state and rigid transformations

Import `Ray` from `gem.ray`.
[Source](../../gem/ray.py), [geometry/ownership tests](../../tests/test_rays.py),
[ray guide](../RAYS.md).

Ray stores `.start` (Vector3 origin), `.dir` (Vector3 normalized direction),
`.distance` (original supplied direction magnitude) and `.end` (historical zero
Vector3 intersection placeholder/state). `.end` is **not** start+dir*distance.
There is no hit validity flag, intersection routine, shadow/visibility transport,
public module constant or ctypes export. Distance remains stored scalar state,
not a newly defined parameter interval or validated hit distance.

Construction retains caller start/direction Vectors and normalizes the caller's
direction in place (replacing that Vector's component list). Rigid transforms
replace start/dir, preserve distance and leave end unchanged; they return None.
All rotations are about the coordinate origin, not the ray start. Scaling/shear,
projective transforms and general invalid inputs are outside the established
ray contract.

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `Ray` | Stateful ray wrapper; use constructor and established Vector3 prerequisites. |
| `Ray.__init__(self, startVector, dirVector)` | `startVector,dirVector`: Vector3 references retained; distance=original dir magnitude, require nonzero direction (ZeroDivisionError), i_normalize caller direction, create end zero Vector3. Returns None. |
| `Ray.duplicate(self)` | Fresh Ray without rerunning constructor; independently copy start,dir,end Vectors and component lists, preserve exact stored distance and direction. Original unchanged. |
| `Ray.roateUsingMatrix(self, matrix)` | Historical misspelling retained. `matrix`: Matrix3 rotation; replace start/dir by row-vector products, normalize rotated direction, return None. Origin pivot; matrix/end/distance preserved, zero resulting direction raises ZeroDivisionError. |
| `Ray.rotateUsingQuaternion(self, quat1)` | `quat1`: unit Quaternion; Hamilton conjugate sandwich rotates start and dir about origin; retain Vector3 types, normalize resulting direction; return None. Preserve quaternion, end and distance; zero resulting dir ZeroDivisionError. No automatic normalization of quat1. |
| `Ray.translate(self, matrix)` | `matrix`: pure translation Matrix4; locally append w=1 to Vector3 start and w=0 to dir, transform and retain XYZ. Replace start/dir, preserve direction under pure translation, distance and end; return None. Other matching-dimensional cases use general Matrix*Vector without new promotion; no projective divide. |
| `Ray.output(self)` | Print Ray, start and dir to stdout; return None, no mutation. |

```python
from gem.ray import Ray
from gem.vector import Vector
from gem.matrix import Matrix
from gem.quaternion import quat_from_axis_angle

start = Vector(3, [1, 0, 0])
direction = Vector(3, [2, 0, 0])
r = Ray(start, direction)
assert r.start is start and r.dir is direction
assert direction.vector == [1.0, 0.0, 0.0] and r.distance == 2
copy = r.duplicate()
assert copy.start is not start and copy.start.vector is not start.vector
assert copy.end is not r.end and copy.distance == 2
assert r.rotateUsingQuaternion(quat_from_axis_angle([0, 0, 1], 90)) is None
assert abs(r.start.vector[1] - 1) < 1e-14
assert r.translate(Matrix(4).translate(Vector(3, [2, 3, 4]))) is None
assert all(abs(a - b) < 1e-14 for a, b in zip(r.start.vector, [2, 4, 4]))
assert r.distance == 2 and r.end.vector == [0.0, 0.0, 0.0]
assert start.vector == [1, 0, 0]
```

Whether end represents a valid hit, and whether/how it should transform, remains
an [open decision](decisions.md). A nonzero end does not establish validity,
and zero can be a hit at the origin. See [matrix](matrix.md), [quaternion](quaternion.md),
[plane](plane.md) and [index](index.md).
