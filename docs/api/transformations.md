# Transform builders and vector application

Import matrix builders/methods from `gem.matrix`; the raw vector transform
kernel and Vector methods are in `gem.vector`.
[Source](../../gem/matrix.py), [transform tests](../../tests/test_transformations.py),
[pivot/shear](../../tests/test_pivot_shear.py).

All matrix methods here **postmultiply** the receiver: result=M*builder.
Returning methods allocate fresh Matrix/ctypes storage; `i_` methods replace rows
and ctypes and return self. Inputs remain unchanged. Wrapper scale/rotate/translate
require a Vector argument and otherwise raise TypeError. Raw builders use
indexable coordinate sequences, not Vector wrappers. Angles are degrees unless
explicitly marked radians. Positive +Z rotation maps +X toward +Y; this is not
an enforced application-wide world frame.

Axis-angle builders use Rodrigues rotation: for normalized axis n, the
row-vector matrix is cos(theta)*I+(1−cos(theta))*n*n^T−sin(theta)*[n]_cross,
where [n]_cross*v=n×v in column notation. This is the transpose of the
conventional column-vector active matrix.

## Matrix builders and methods

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `scale(size, value)` | Raw `value`: coordinate scale factors; fresh diagonal size×size list. First up to three diagonal entries use value[x]; entries x>=3 are 1.0, preserving homogeneous W for size 4. |
| `shearXY3(x, y)` | Numeric `x,y`; fresh 3×3 rows implementing Z′=Z+xX+yY; other coordinates unchanged. |
| `shearYZ3(y, z)` | Numeric `y,z`; fresh 3×3 rows implementing X′=X+yY+zZ; other coordinates unchanged. |
| `shearXZ3(x, z)` | Numeric `x,z`; fresh 3×3 rows implementing Y′=Y+xX+zZ; other coordinates unchanged. |
| `shearXY4(x, y)` | Numeric `x,y`; fresh 4×4 rows implementing Z′=Z+xX+yY; other coordinates unchanged, homogeneous W preserved. |
| `shearYZ4(y, z)` | Numeric `y,z`; fresh 4×4 rows implementing X′=X+yY+zZ; other coordinates unchanged, homogeneous W preserved. |
| `shearXZ4(x, z)` | Numeric `x,z`; fresh 4×4 rows implementing Y′=Y+xX+zZ; other coordinates unchanged, homogeneous W preserved. |
| `translate2(vector)` | Raw `vector` XY sequence → fresh homogeneous 3×3 translation; last row [x,y,1]. |
| `translate3(vector)` | Raw `vector` XYZ sequence → fresh 3×3 with last row [x,y,z]. Historical linear-row replacement, not general 3D affine translation. |
| `translate4(vector)` | Raw `vector` XYZ (extra W ignored) → fresh homogeneous 4×4; last row [x,y,z,1]. |
| `rotate2(point, theta)` | Raw 2D `point` pivot; `theta` degrees; fresh 3×3 homogeneous p′=pivot+(p−pivot)R. Pivot stationary; W=0 directions receive origin-only linear rotation. |
| `rotate3(axis, theta)` | Raw 3D `axis`, `theta` degrees; temporary normalized axis; fresh 3×3 Rodrigues row matrix. Exact zero axis raises ZeroDivisionError. |
| `rotate4(axis, theta)` | Raw 3D `axis`, `theta` degrees; temporary normalization; fresh 4×4 Rodrigues row matrix preserving W. Zero axis raises ZeroDivisionError. |
| `rotate_origin2(theta)` | `theta` radians; fresh homogeneous 3×3 counterclockwise XY rotation about origin. |
| `Matrix.i_scale(self, value)` | `value`: Vector with enough coordinate scale factors; postmultiply by the scale diagonal. Replace rows/ctypes and return self. No zero-factor rejection; size 4 preserves W. |
| `Matrix.scale(self, value)` | `value`: Vector with enough coordinate scale factors; postmultiply by the scale diagonal. Fresh Matrix. No zero-factor rejection; size 4 preserves W. |
| `Matrix.i_rotate(self, axis, theta)` | `axis`: Vector; `theta`: degrees. Size 2 uses the linear block of rotate2, ignoring pivot translation (origin-only); size 3/4 uses normalized axis-angle builder. Other sizes NotImplementedError. Replace rows/ctypes, return self. Nonzero 3D axis prerequisite; zero raises ZeroDivisionError. |
| `Matrix.rotate(self, axis, theta)` | `axis`: Vector; `theta`: degrees. Size 2 uses the linear block of rotate2, ignoring pivot translation (origin-only); size 3/4 uses normalized axis-angle builder. Other sizes NotImplementedError. Fresh Matrix. Nonzero 3D axis prerequisite; zero raises ZeroDivisionError. |
| `Matrix.i_translate(self, vecA)` | `vecA`: Vector. Matrix3/Vector2 → affine 2D; Matrix3/Vector3 → legacy last-row replacement; Matrix4/Vector3 or Vector4 → final-row XYZ translation. Matrix2/other receiver sizes NotImplementedError. Replace rows/ctypes and return self. Unsupported argument dimensions have historical, nonuniform errors. |
| `Matrix.translate(self, vecA)` | `vecA`: Vector. Matrix3/Vector2 → affine 2D; Matrix3/Vector3 → legacy last-row replacement; Matrix4/Vector3 or Vector4 → final-row XYZ translation. Matrix2/other receiver sizes NotImplementedError. Fresh Matrix. Unsupported argument dimensions have historical, nonuniform errors. |
| `Matrix.shearXY(self, x, y)` | Numeric `x,y`; postmultiply Z′=Z+xX+yY shear for sizes 3/4 (size 4 preserves W), otherwise NotImplementedError. Fresh Matrix. |
| `Matrix.i_shearXY(self, x, y)` | Numeric `x,y`; postmultiply Z′=Z+xX+yY shear for sizes 3/4 (size 4 preserves W), otherwise NotImplementedError. Replace rows/ctypes, return self. |
| `Matrix.shearYZ(self, y, z)` | Numeric `y,z`; postmultiply X′=X+yY+zZ shear for sizes 3/4 (size 4 preserves W), otherwise NotImplementedError. Fresh Matrix. |
| `Matrix.i_shearYZ(self, y, z)` | Numeric `y,z`; postmultiply X′=X+yY+zZ shear for sizes 3/4 (size 4 preserves W), otherwise NotImplementedError. Replace rows/ctypes, return self. |
| `Matrix.shearXZ(self, x, z)` | Numeric `x,z`; postmultiply Y′=Y+xX+zZ shear for sizes 3/4 (size 4 preserves W), otherwise NotImplementedError. Fresh Matrix. |
| `Matrix.i_shearXZ(self, x, z)` | Numeric `x,z`; postmultiply Y′=Y+xX+zZ shear for sizes 3/4 (size 4 preserves W), otherwise NotImplementedError. Replace rows/ctypes, return self. |

## Vector transform

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `transform(size, position, matrix)` | `size`: output dimension; raw `position` and nested `matrix`; fresh list row product. Same-size computes all components, including explicit W. Matrix of size+1 supplies local W=1 and returns size numerators. No perspective divide or input mutation. |

`Vector.transform(self, position, matrix)` and `Vector.i_transform` use the receiver
size and explicit raw position, not its stored coordinates; see [Vector methods](vector.md).
General Matrix*Vector still requires matching dimensions. Local implicit promotion
is intended for affine position transforms. With projective matrices, truncated
numerators without output-W division are not perspective-correct coordinates;
use [project/unproject](projection.md) or explicit homogeneous processing.

```python
from gem.matrix import Matrix, rotate2, shearXY4
from gem.vector import Vector, transform

pivot = [2, 3]
r = Matrix(3, rotate2(pivot, 90))
assert all(abs(a - b) < 1e-14 for a, b in zip((r * Vector(3, [3, 3, 1])).vector, [2, 4, 1]))
assert all(abs(a - b) < 1e-14 for a, b in zip((r * Vector(3, [2, 3, 1])).vector, [2, 3, 1]))
t = Matrix(4).translate(Vector(3, [2, -3, 4]))
assert transform(3, [1, 2, 3], t.matrix) == [3.0, -1.0, 7.0]
assert (t * Vector(4, [1, 2, 3, 0])).vector == [1.0, 2.0, 3.0, 0.0]
assert (Matrix(4, shearXY4(2, -1)) * Vector(4, [3, 4, 5, 7])).vector == [3.0, 4.0, 7.0, 7.0]
rot = Matrix(4).rotate(Vector(3, [0, 0, 1]), 90)
a = (t * rot) * Vector(4, [0, 0, 0, 1])
b = (rot * t) * Vector(4, [0, 0, 0, 1])
assert all(abs(x - y) < 1e-14 for x, y in zip(a.vector, [3, 2, 4, 1]))
assert b.vector == [2.0, -3.0, 4.0, 1.0]
```

Homogeneous pivot construction is not exposed as an additional Matrix overload.
Unsupported translate dimensions and direct shape edits need a separate policy,
not an invented exception guarantee. See [decisions](decisions.md) and [index](index.md).
