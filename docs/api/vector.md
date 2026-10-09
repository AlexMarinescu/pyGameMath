# Vector mathematics

Import `Vector` and the functions below from `gem.vector`.
[Source](../../gem/vector.py), [conventions](../architecture/conventions.md),
[evidence](../../tests/test_vector_common.py), [ownership/extremes](../../tests/test_vector_quaternion_optimization.py).

`Vector(size, data=None)` is the class for Vector2/3/4: these are dimensions,
not distinct class names. `.size` is the declared dimension and `.vector` is
component storage. No wrapper indexing or iteration protocol is implemented:
use `v.vector[index]` and `for component in v.vector`. There is no Vector ctypes
export; use [common conversion helpers](common.md).

Unless stated otherwise, `size` is a nonnegative integer and raw `vecA`/`vecB`
are indexable component sequences with at least `size` entries. Wrapper operands
use matching well-formed dimensions. Basic kernels allocate fresh lists and
preserve inputs; they do not validate every shape/type. Arithmetic mismatches can
index-fail or ignore extra storage rather than uniformly raising ValueError.
These are prerequisites, not an expanded validation contract.

## Class and methods

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `Vector` | Dimensioned component wrapper; instantiate with the constructor below. |
| `Vector.__init__(self, size, data=None)` | `size` declares dimension; omitted `data` creates zero storage, supplied component storage is retained by reference without length validation. Initialization returns None; see fields above. |
| `Vector.__repr__(self)` | Return a string containing the declared size and storage; no mutation. |
| `Vector.__add__(self, other)` | `other`: Vector or int/float (bool inherits int). Component +; return a fresh Vector. Unsupported types return NotImplemented; no reflected scalar operator. |
| `Vector.__iadd__(self, other)` | `other`: Vector or int/float (bool inherits int). Component +; replace receiver storage and return self. Unsupported types return NotImplemented; no reflected scalar operator. |
| `Vector.__sub__(self, other)` | `other`: Vector or int/float (bool inherits int). Component -; return a fresh Vector. Unsupported types return NotImplemented; no reflected scalar operator. |
| `Vector.__isub__(self, other)` | `other`: Vector or int/float (bool inherits int). Component -; replace receiver storage and return self. Unsupported types return NotImplemented; no reflected scalar operator. |
| `Vector.__mul__(self, scalar)` | `scalar`: int/float; component multiplication; fresh Vector. Unsupported types return NotImplemented; scalar-left multiplication is absent. |
| `Vector.__imul__(self, scalar)` | `scalar`: int/float; component multiplication; replace storage, return self. Unsupported types return NotImplemented; scalar-left multiplication is absent. |
| `Vector.__div__(self, scalar)` | `scalar`: int/float; component /; fresh Vector. Unsupported operands return NotImplemented; division by exact zero raises ZeroDivisionError. |
| `Vector.__truediv__(self, scalar)` | `scalar`: int/float; component /; fresh Vector. Unsupported operands return NotImplemented; division by exact zero raises ZeroDivisionError. |
| `Vector.__idiv__(self, scalar)` | `scalar`: int/float; component /; replace storage, return self. Unsupported operands return NotImplemented; division by exact zero raises ZeroDivisionError. |
| `Vector.__itruediv__(self, scalar)` | `scalar`: int/float; component /; replace storage, return self. Unsupported operands return NotImplemented; division by exact zero raises ZeroDivisionError. |
| `Vector.__eq__(self, vecB)` | `vecB`: Vector. Return exact equality bool including declared dimensions; two empty Vectors are equal. Unsupported types return NotImplemented. No tolerance; same-size comparison reads each declared component. |
| `Vector.__ne__(self, vecB)` | `vecB`: Vector. Return exact inequality bool including declared dimensions; two empty Vectors are equal. Unsupported types return NotImplemented. No tolerance; same-size comparison reads each declared component. |
| `Vector.__neg__(self)` | Unary minus; return a fresh component-negated Vector. |
| `Vector.clone(self)` | Copy the component sequence by slicing into a new Vector (ordinary numeric components are independent). |
| `Vector.one(self)` | Replace storage with ones; return self. |
| `Vector.zero(self)` | Replace storage with zeros; return self. |
| `Vector.negate(self)` | Return a fresh component-negated Vector, equivalent to unary minus. |
| `Vector.maxV(self, vecB)` | `vecB`: matching Vector; fresh componentwise maximum Vector. Inputs preserved; comparison-based NaN behavior is historical. |
| `Vector.maxS(self)` | Return the largest component by comparisons; empty storage raises IndexError. |
| `Vector.minV(self, vecB)` | `vecB`: matching Vector; fresh componentwise minimum Vector. Inputs preserved; comparison-based NaN behavior is historical. |
| `Vector.minS(self)` | Return the smallest component by comparisons; empty storage raises IndexError. |
| `Vector.magnitude(self)` | Return the stable numeric length; no mutation. |
| `Vector.clamp(self, size, value, minS, maxS)` | Explicit `size,value,minS,maxS` raw component lists are passed to module clamp; receiver values are not implicit. Return a fresh Vector; receiver and lists preserved. |
| `Vector.i_clamp(self, size, value, minS, maxS)` | Explicit `size,value,minS,maxS` raw component lists are passed to module clamp; receiver values are not implicit. Replace only receiver storage, return self; original caller lists preserved. Receiver .size is not updated if explicit size differs. |
| `Vector.i_normalize(self)` | Replace storage with stable normalized components; exact zero becomes zero; return self. |
| `Vector.normalize(self)` | Return a fresh normalized Vector; exact zero returns a fresh zero Vector. |
| `Vector.dot(self, vecB)` | `vecB`: Vector; numeric sum of products, preserving inputs; unsupported type returns NotImplemented. No stable extreme-product summation guarantee. |
| `Vector.isInSameDirection(self, otherVec)` | `otherVec`: Vector; return dot(otherVec) > 0, a bool sign test, not collinearity. Zero dot yields False; unsupported type returns NotImplemented. |
| `Vector.isInOppositeDirection(self, otherVec)` | `otherVec`: Vector; return dot(otherVec) < 0, a bool sign test, not collinearity. Zero dot yields False; unsupported type returns NotImplemented. |
| `Vector.barycentric(self, a, b, c)` | `a,b,c`: Vector triangle vertices; receiver is point p. Return fresh [u,v,w], sum 1, from the dot-product Gram system. Exact zero denominator raises ZeroDivisionError; outside-triangle weights may be negative. No clamping or mutation. |
| `Vector.transform(self, position, matrix)` | `position`: raw coordinates; `matrix`: nested square rows; use receiver size, not its stored position. Fresh Vector; receiver preserved. See [transform rules](transformations.md#vector-transform). |
| `Vector.i_transform(self, position, matrix)` | `position`: raw coordinates; `matrix`: nested square rows; use receiver size, not its stored position. Replace receiver storage, return self; raw inputs preserved. See [transform rules](transformations.md#vector-transform). |
| `Vector.xy(self)` | Return fresh Vector2 with components X, Y in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.yz(self)` | Return fresh Vector2 with components Y, Z in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.xz(self)` | Return fresh Vector2 with components X, Z in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.xw(self)` | Return fresh Vector2 with components X, W in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.yw(self)` | Return fresh Vector2 with components Y, W in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.zw(self)` | Return fresh Vector2 with components Z, W in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.xyw(self)` | Return fresh Vector3 with components X, Y, W in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.yzw(self)` | Return fresh Vector3 with components Y, Z, W in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.xzw(self)` | Return fresh Vector3 with components X, Z, W in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.xyz(self)` | Return fresh Vector3 with components X, Y, Z in that order; receiver must contain every referenced component, otherwise IndexError. No mutation. |
| `Vector.right(self)` | Return a fresh fixed Vector3 [1,0,0] regardless of receiver dimension. `front` differs from Quaternion.getForward (+Z). |
| `Vector.left(self)` | Return a fresh fixed Vector3 [-1,0,0] regardless of receiver dimension. `front` differs from Quaternion.getForward (+Z). |
| `Vector.front(self)` | Return a fresh fixed Vector3 [0,0,-1] regardless of receiver dimension. `front` differs from Quaternion.getForward (+Z). |
| `Vector.back(self)` | Return a fresh fixed Vector3 [0,0,1] regardless of receiver dimension. `front` differs from Quaternion.getForward (+Z). |
| `Vector.up(self)` | Return a fresh fixed Vector3 [0,1,0] regardless of receiver dimension. `front` differs from Quaternion.getForward (+Z). |
| `Vector.down(self)` | Return a fresh fixed Vector3 [0,-1,0] regardless of receiver dimension. `front` differs from Quaternion.getForward (+Z). |

## Raw kernels and geometry helpers

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `zero_vector(size)` | Fresh zero-filled list; `size` controls length. Sizes 2/3/4 copy the historical reference buffers. |
| `one_vector(size)` | Fresh one-filled list; sizes 2/3/4 copy the historical reference buffers. |
| `lerp(vecA, vecB, time)` | `vecA,vecB`: matching Vectors; numeric `time` is unclamped. Return fresh a+(b-a)*time Vector. |
| `cross(vecA, vecB)` | `vecA,vecB`: Vector3; fresh Vector3 a×b, right-handed X×Y=Z; inputs preserved. |
| `reflect(incidentVec, normal)` | `incidentVec,normal`: matching Vectors, normal unit prerequisite; fresh i−2(i·n)n; no automatic normalization. |
| `refract(IOR, incidentVec, normal)` | `IOR`: n1/n2; matching unit `incidentVec,normal`, normal opposes incidence. k=1−IOR²(1−(n·i)²); fresh IOR*i−(IOR*(n·i)+sqrt(k))*n, or zero Vector when k<0. No input normalization/flipping. |
| `toAngle(vector)` | Raw 2D `vector` sequence; atan2(y,x) radians, a scalar; no mutation. |
| `lperp(vector)` | Raw 2D `vector`; fresh Vector2 [-y,x], a +90° rotation. |
| `rperp(vector)` | Raw 2D `vector`; fresh Vector2 [y,-x], a −90° rotation. |
| `vec_add(size, vecA, vecB)` | Fresh list vecA[i]+vecB[i]. |
| `s_vec_add(size, vecA, scalar)` | Fresh list vecA[i]+scalar. |
| `vec_sub(size, vecA, vecB)` | Fresh list vecA[i]−vecB[i]. |
| `s_vec_sub(size, vecA, scalar)` | Fresh list vecA[i]−scalar. |
| `vec_mul(size, vecA, scalar)` | Fresh list vecA[i]*scalar. |
| `vec_div(size, vecA, scalar)` | Fresh list vecA[i]/scalar; exact zero divisor raises ZeroDivisionError. |
| `vec_neg(size, vecA)` | Fresh list −vecA[i]. |
| `dot(size, vecA, vecB)` | Numeric sum vecA[i]*vecB[i]; no mutation or extreme-scale product guarantee. |
| `magnitude(size, vecA)` | Stable finite-input length using scaled chained hypot. Empty/zero length is 0.0; representational overflow may return infinity; nonfinite paths retain square-sum arithmetic. |
| `normalize(size, vecA)` | Fresh stable normalized list using scaled hypot; exact zero returns zeros. Nonfinite paths retain legacy arithmetic without a new error policy. |
| `maxV(size, vecA, vecB)` | Fresh componentwise maxima of raw vecA and vecB. |
| `minV(size, vecA, vecB)` | Fresh componentwise minima of raw vecA and vecB. |
| `maxS(size, vecA)` | Largest component scalar; seeds from vecA[0], so empty input raises IndexError. |
| `minS(size, vecA)` | Smallest component scalar; seeds from vecA[0], so empty input raises IndexError. |
| `clamp(size, value, minS, maxS)` | Raw `value,minS,maxS` lists and explicit `size`. Copy value, cap above maxS then below minS per component, return Vector(size). Preserve all lists; no validation of ordered bounds. Extra value storage survives the slice. |

## Exposed reference buffers

`gem.vector.REFRENCE_VECTOR_2`, `REFRENCE_VECTOR_3`, `REFRENCE_VECTOR_4` are
mutable lists of respectively 2/3/4 zeros. `IREFRENCE_VECTOR_2`,
`IREFRENCE_VECTOR_3`, `IREFRENCE_VECTOR_4` contain ones. Spellings are historical.
They are implementation reference buffers, not immutable mathematical constants;
mutating them changes subsequent zero/one defaults. Do not edit them in applications.
There are no other intended public vector constants.

## Ownership and stable norms

Returning operations preserve caller buffers. In-place operations replace the
receiver list, so a separate wrapper sharing the old list is not updated.
Finite lengths/normalization avoid unnecessary intermediate overflow/underflow;
this does not stabilize dot, cross or barycentric products. No general NaN/Infinity
policy or mismatched-dimension error policy is established.

```python
from gem.vector import Vector, cross, refract

values = [3.0, 4.0, 0.0]
v = Vector(3, values)
unit = v.normalize()
assert v.vector is values and values == [3.0, 4.0, 0.0]
assert unit.vector is not values and unit.vector == [0.6, 0.8, 0.0]
assert v.i_clamp(3, values, [0, 0, 0], [2, 2, 2]) is v
assert v.vector == [2, 2, 0.0] and values == [3.0, 4.0, 0.0]
assert Vector(0) == Vector(0) and Vector(2) != Vector(3)
assert Vector(2).__eq__(object()) is NotImplemented
assert cross(Vector(3, [1, 0, 0]), Vector(3, [0, 1, 0])).vector == [0, 0, 1]
assert refract(1.0 / 1.5, Vector(3, [0, -1, 0]), Vector(3, [0, 1, 0])).vector == [0.0, -1.0, 0.0]
assert Vector(3).normalize().vector == [0.0, 0.0, 0.0]
```

Barycentric weights describe affine coordinates, not an intersection test:

```python
from gem.vector import Vector

p = Vector(2, [0.5, 0.5])
weights = p.barycentric(Vector(2, [0, 0]), Vector(2, [2, 0]), Vector(2, [0, 2]))
assert weights == [0.5, 0.25, 0.25]
```

See [transforms](transformations.md), [utilities/viewport](common.md),
[decisions](decisions.md) and [API index](index.md).

See the [graphics gallery example](../examples/gallery/vectors.md) for an executable visualization.
