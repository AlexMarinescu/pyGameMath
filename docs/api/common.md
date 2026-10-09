# Scalar, viewport and interoperability helpers

Import these functions and `GLfloat` from `gem.common`.
[Source](../../gem/common.py), [utility tests](../../tests/test_vector_common.py),
[viewport contracts](../../tests/test_final_vector_contracts.py),
[angles](../../tests/test_angles_refraction.py).

These historical raw sequence utilities preserve inputs and return independent
containers unless noted. Shape/type validation is limited; short/malformed
sequences can raise IndexError/TypeError from direct indexing or ctypes conversion.
Do not infer a uniform invalid-input API from ordinary examples.
`GLfloat` is exactly `ctypes.c_float`, a binary32 type alias, not an OpenGL dependency.
It is the only intended public common constant/type alias; imported math, sm and
ct are implementation dependencies, not mathematical exports.

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `convertArr(l, n)` | `l`: flat sliceable sequence, `n`: positive integer chunk width; return new list of slices l[i:i+n], final chunk may be shorter. n=0 raises ValueError from range; negative step yields empty result for ordinary input. Slice ownership follows supplied sequence semantics. |
| `mulV4(v1, v2)` | `v1,v2`: raw four-component numeric sequences; fresh list of four componentwise products, not dot/matrix multiplication. |
| `conv_list(listIn, cType)` | `listIn`: numeric/value list; `cType`: ctypes element type; fresh initialized cType*len array. Native ctypes conversion errors propagate; inputs preserved. |
| `conv_list_2d(listIn, cType)` | `listIn`: nonempty rectangular row lists; `cType`: ctypes type; fresh row-major nested ctypes array with lengths inferred from row count/first row. Empty outer list raises IndexError; mismatched shapes retain index behavior. |
| `list_2d_to_1d(inlist)` | `inlist`: nonempty rectangular rows; fresh row-major flat list of numeric elements. Width from first row; no ragged validation or deep copying of arbitrary objects. |
| `convertM4to3(matrix)` | Raw nested `matrix`: at least leading 3×3 entries, normally 4×4; fresh copy of leading 3×3 rows. Not a Matrix wrapper or affine transform adapter. |
| `sinc(x)` | Numeric `x`: 1.0 when abs(x)<1e−4, otherwise sin(x)/x; historical approximation at small angles, no mutation. Not the exact center-limit implementation of angular-probe quadrature. |
| `scalarLerp(a, b, time)` | Numeric `a,b,time`; scalar a+time*(b−a), no clamping or mutation. |
| `getViewPort(coords, width, height)` | `coords`: finite nonzero Vector2/3/4, numeric width/height; fresh four-element list from whole-vector normalization with original XY offsets. Zero Vector raises ZeroDivisionError; inputs preserved. Formula below. |
| `radiansToDegrees(degrees)` | Numeric `degrees` (misleading retained parameter name) actually represents radians; return input*180/pi. |
| `degreesToRadians(radians)` | Numeric `radians` (misleading retained parameter name) actually represents degrees; return input*pi/180. |
| `sign(x)` | Numeric `x`; return +1.0 if x≥0, including signed zero, otherwise −1.0. NaN follows else branch; no new nonfinite policy. |

## Historical viewport

For L=norm(coords), output is
[(x/L+1)*width/2+x, (y/L+1)*height/2+y, width, height].
Z/W affect the normalization for Vector3/4; original XY are offsets.
This is neither OpenGL glViewport nor conventional NDC-to-window projection.
Use [project/unproject](projection.md) for those coordinate conversions.
Widths/heights and wider malformed/nonfinite domains are not standardized here.

## ctypes ownership

ctypes helpers create independent arrays; converting to binary32 can round/overflow
binary64 values. Keep the array alive while foreign code uses a pointer. Matrix
snapshots are documented in [matrix](matrix.md); no foreign graphics upload is
verified by these memory checks.

```python
import ctypes
import math
from gem.common import GLfloat, conv_list_2d, list_2d_to_1d, convertArr, getViewPort, radiansToDegrees, degreesToRadians
from gem.vector import Vector

rows = [[1, 2], [3, 4]]
array = conv_list_2d(rows, GLfloat)
assert GLfloat is ctypes.c_float and list(array[1]) == [3.0, 4.0]
rows[1][0] = 99
assert array[1][0] == 3.0
assert list_2d_to_1d([[1, 2], [3, 4]]) == [1, 2, 3, 4]
assert convertArr([1, 2, 3], 2) == [[1, 2], [3]]
coords = Vector(2, [3, 4])
assert getViewPort(coords, 100, 200) == [83.0, 184.0, 100, 200]
assert coords.vector == [3, 4]
assert radiansToDegrees(math.pi) == 180.0
assert degreesToRadians(180) == math.pi
```

See [vector](vector.md), [decisions](decisions.md) and [index](index.md).
