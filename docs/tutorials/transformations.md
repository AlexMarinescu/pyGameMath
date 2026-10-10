# Transform a group of object points

**Intermediate.** Prerequisites: [vectors](vectors.md) and basic matrix products.
Learn to compose scale, rotation and translation, distinguish points from
directions and recover object-local coordinates. These calculations are useful
for object placement, attachments and camera-relative tools.

## gem's row-vector convention

Rows are stored in `.matrix[row][column]`. The matrix product C=A*B is ordinary
C_ij=sum_k A_ik B_kj, but the wrapper expression `M * v` evaluates the mathematical
row product **vM**. Therefore `(S*R*T)*p` applies scale, then rotation, then
translation. Translation is in the final row. Do not import column-vector
composition order from another engine.

A homogeneous position [x,y,z,1] receives translation; a direction [x,y,z,0]
does not. Scale changes direction length; rotation alone preserves it. For XYZ
scale [2,1,1], +90° about Z and translation [10,-2,3],
(x,y,z) → (2x,y,z) → (-y,2x,z) → (10-y,2x-2,z+3).
This coordinate derivation gives independent answers for the point group.

## Worked local-to-world conversion

```python
# With row vectors, the first matrix in the product acts first.
from gem.matrix import Matrix
from gem.vector import Vector

local = [Vector(4, [0, 0, 0, 1]), Vector(4, [1, 0, 0, 1]), Vector(4, [0, 1, 0, 1])]
scale = Matrix(4).scale(Vector(3, [2, 1, 1]))
rotation = Matrix(4).rotate(Vector(3, [0, 0, 1]), 90)  # degrees
translation = Matrix(4).translate(Vector(3, [10, -2, 3]))
model = scale * rotation * translation
world = [model * point for point in local]
expected = [[10, -2, 3, 1], [10, 0, 3, 1], [9, -2, 3, 1]]
assert all(abs(a - b) < 1e-14 for p, e in zip(world, expected) for a, b in zip(p.vector, e))
world_direction = model * Vector(4, [1, 0, 0, 0])
assert all(abs(a - b) < 1e-14 for a, b in zip(world_direction.vector, [0, 2, 0, 0]))
recovered = [model.inverse() * point for point in world]
assert all(abs(a - b) < 1e-13 for p, e in zip(recovered, local) for a, b in zip(p.vector, e.vector))
other_order = (translation * rotation) * Vector(4, [0, 0, 0, 1])
assert all(abs(a - b) < 1e-14 for a, b in zip(other_order.vector, [2, 10, 3, 1]))
assert local[1].vector == [1, 0, 0, 1]
print([[round(x, 6) for x in point.vector[:3]] for point in world])
```

Output: `[[10.0, -2.0, 3.0], [10.0, 0.0, 3.0], [9.0, -2.0, 3.0]]`.
The local X axis becomes twice-length +Y; moving the translation before rotation
instead rotates its offset to [2,10,3]. Thus changing order changes the frame in
which a transformation acts. Inverse(model) converts the world coordinates back
to the original local frame, without mutating the points or model.

## Homogeneous and ownership boundaries

General Matrix*Vector requires matching dimensions. Use explicit Vector4 here;
it does not promote Vector3 automatically. Vector.transform has separate local
affine promotion rules and performs no perspective divide. Those rules are not
an interchangeable camera projection API.

Returning matrix transforms allocate fresh rows; in-place variants postmultiply
the receiver and synchronize ctypes. Nonuniform scale requires a separate
inverse-transpose treatment for surface normals: applying the position/direction
matrix indiscriminately does not preserve their perpendicularity.

Zero scale can make the model singular; inverse then raises ZeroDivisionError.
Very ill-conditioned transforms and float32 exports have numerical limits even
when Python rows are representable. See [matrix API](../api/matrix.md),
[transform conventions](../api/transformations.md) and [accuracy](numerical.md).
Continue with [camera coordinates](camera.md), [interoperability](interop.md) or
[the tutorial index](index.md).

See the [visual example](../examples/gallery/transforms.md) and its reproducible assets.
