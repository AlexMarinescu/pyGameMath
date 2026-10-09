# World coordinates, camera space and window coordinates

**Intermediate.** Prerequisites: [transform order](transformations.md), division
and trigonometry. Learn to trace a perspective point through every stage, then
recover world coordinates. These calculations support screen markers and picking,
without requiring a renderer.

## The camera frame and clip coordinates

`lookAt(eye,center,up)` creates a right-handed view matrix looking down camera -Z.
f=normalize(center-eye), s=normalize(f×up), u=s×f form its frame. Here eye=[0,0,5],
center=[0,0,0], up=+Y, so world [1,1,0] becomes camera [1,1,-5]. Inputs are preserved;
eye=center, zero up or parallel up/view direction is degenerate.

For vertical FOV=90°, aspect=2, near=1 and far=10, perspective XY scale is
[1/2,1]. Depth entries are c=-(far+near)/(far-near)=-11/9 and
d=-2*far*near/(far-near)=-20/9. The row product of camera [1,1,-5,1]
produces clip [1/2,1,35/9,5]. After division by W, NDC=[1/10,1/5,7/9].
For viewport [10,20,200,100], window coordinates are [120,80,8/9].

Viewport origin is lower-left, with Y increasing upward. OpenGL-style NDC Z
[-1,1] maps to window Z=(NDC_Z+1)/2, without clamping. This differs from APIs
using NDC depth [0,1] or top-left raster origins; convert explicitly when integrating.

## Trace projection and reverse it

```python
from gem.matrix import Matrix, lookAt, perspective, project, unproject
from gem.vector import Vector

eye = Vector(3, [0, 0, 5])
view = lookAt(eye, Vector(3, [0, 0, 0]), Vector(3, [0, 1, 0]))
projection = perspective(90, 2, 1, 10)  # vertical degrees; aspect=width/height
viewport = [10, 20, 200, 100]
world = Vector(4, [1, 1, 0, 1])
camera = view * world
clip = (view * projection) * world
assert camera.vector == [1.0, 1.0, -5.0, 1.0]
assert all(abs(a - b) < 1e-14 for a, b in zip(clip.vector, [0.5, 1, 35.0 / 9, 5]))
ndc = [value / clip.vector[3] for value in clip.vector[:3]]
assert all(abs(a - b) < 1e-14 for a, b in zip(ndc, [0.1, 0.2, 7.0 / 9]))
window = project(world, view, projection.matrix, viewport)
assert all(abs(a - b) < 1e-13 for a, b in zip(window.vector, [120, 80, 8.0 / 9]))
back = unproject(window.vector[0], window.vector[1], window.vector[2], view.matrix, projection, viewport)
assert all(abs(a - b) < 1e-13 for a, b in zip(back.vector, [1, 1, 0]))
near = project(Vector(4, [0, 0, 4, 1]), view, projection, viewport)
far = project(Vector(4, [0, 0, -5, 1]), view, projection, viewport)
assert abs(near.vector[2]) < 1e-14 and abs(far.vector[2] - 1) < 1e-14
assert world.vector == [1, 1, 0, 1] and eye.vector == [0, 0, 5]
try:
    project(Vector(4, [0, 0, 5, 1]), view, projection, viewport)
except ZeroDivisionError:
    pass
else:
    raise AssertionError('eye-plane projection must reject zero clip W')
print([round(x, 6) for x in window.vector])
```

Output: `[120.0, 80.0, 0.888889]`. The independent clip calculation checks more
than a round trip, since two incorrect conversions could otherwise cancel.
Raw 4×4 lists and Matrix wrappers can be mixed; project requires explicit Vector4.
If an object model matrix is also present, row composition is model*view*projection.
The inverse then recovers object coordinates, not automatically world coordinates.

## Boundaries and integration choices

A point at the eye plane has clip W=0: project raises ZeroDivisionError. Points
behind the camera or outside its frustum are not rejected/clipped here; a sensible
window result alone does not establish visibility. Near-zero W amplifies errors.
Unproject inverts the combined matrix, so singular transforms raise ZeroDivisionError;
zero output W retains the zero Vector3 sentinel. Neither convention is a universal
invalid-input policy. Viewport extent zero and invalid frusta retain native errors.

Perspective depth is nonlinear in camera distance. For top-left input pixels,
first map Y to this lower-left coordinate system; pixel centers versus edges are
an application choice, not an implicit adjustment. `common.getViewPort` uses a
separate historical normalized-coordinate formula and is not used here.

Continue with [screen picking](geometry.md#build-a-picking-ray),
[quaternion cameras](quaternions.md) and [memory interoperability](interop.md).
References: [projection API](../api/projection.md), [conventions](../architecture/conventions.md),
[independent projection tests](../../tests/test_projection.py), [tutorial index](index.md).
