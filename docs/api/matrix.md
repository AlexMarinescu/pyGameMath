# Matrices: storage, algebra and inversion

Import `Matrix` and raw kernels from `gem.matrix`.
[Source](../../gem/matrix.py), [multiplication/inverse tests](../../tests/test_matrix.py),
[scaled inversion](../../tests/test_inverse_optimization.py),
[ctypes/division](../../tests/test_matrix_division.py).

Matrix2/3/4 denote `Matrix(2)`, `Matrix(3)`, `Matrix(4)`, not separate classes.
`.size` declares dimension; `.matrix` is nested row-major storage, retained when
supplied. Index/iterate `.matrix` or its rows; Matrix has no wrapper indexing or
iteration protocol. `.c_matrix` is an independently allocated binary32 ctypes
array snapshot created at construction and refreshed by supported in-place methods.
Direct nested-row edits do not refresh it. No mathematical public constants exist.

Raw kernels expect well-formed square nested numeric lists; no universal shape
validation is promised. `matrix_vector_multiply` takes a Vector wrapper, not a
raw list. Returning kernels allocate new rows, returning wrappers new Matrices;
in-place wrappers replace receiver rows and ctypes storage and return self.

C[i][j]=sum(A[i][k]*B[k][j]) is the ordinary matrix product. Despite `M*v` syntax,
vector application is vM; `(A*B)*v` applies A then B. Matching dimensions are a
prerequisite for Matrix*Vector, with no implicit promotion or uniform mismatch
check. Matrix*Matrix explicitly raises ValueError for differing declared sizes.

## Class and algebra kernels

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `zero_matrix(size)` | `size`: dimension; fresh nested zero rows with independent row storage. |
| `identity(size)` | `size`: dimension; fresh nested identity rows with independent row storage. |
| `matrix_multiply(matrixA, matrixB)` | `matrixA,matrixB`: matching raw square nested lists; ordinary product A*B into fresh rows; uses len(matrixA), no broad validation. |
| `matrix_vector_multiply(matrix, vec)` | Raw square `matrix`, matching `vec`: Vector; fresh Vector of vec.size with row product vec*matrix. No promotion; inputs preserved. |
| `matrix_div(mat, scalar)` | Raw `mat`, numeric `scalar`: fresh rows with mat[i][j]/scalar, no transpose. Kernel arithmetic accepts integers unlike wrapper division; exact zero raises ZeroDivisionError. |
| `transpose(mat)` | Raw square `mat`: fresh nested transpose rows; no mutation. |
| `det2(mat)` | Raw 2×2 `mat`; numeric determinant from direct products/cofactors. No scaling or overflow protection; input preserved. |
| `det3(mat)` | Raw 3×3 `mat`; numeric determinant from direct products/cofactors. No scaling or overflow protection; input preserved. |
| `det4(mat)` | Raw 4×4 `mat`; numeric determinant from direct products/cofactors. No scaling or overflow protection; input preserved. |
| `inverse2(mat)` | Raw 2×2 `mat`; fresh nested inverse. Singular inputs raise ZeroDivisionError. Direct adjugate/determinant arithmetic; no extreme-scale stabilization. Inputs preserved. |
| `inverse3(mat)` | Raw 3×3 `mat`; fresh nested inverse. Singular inputs raise ZeroDivisionError. Power-of-two row scaling plus cofactor inverse, exact represented-coefficient singularity check; limits below. Inputs preserved. |
| `inverse4(mat)` | Raw 4×4 `mat`; fresh nested inverse. Singular inputs raise ZeroDivisionError. Power-of-two row scaling plus cofactor inverse, exact represented-coefficient singularity check; limits below. Inputs preserved. |
| `Matrix` | Small square-matrix wrapper; constructor below covers dimensions and ownership. |
| `Matrix.__init__(self, size, data=None)` | `size`: dimension; `data=None`: identity; supplied nested rows are retained, not copied or comprehensively shape-validated. Create float32 .c_matrix; initialization returns None. |
| `Matrix.__mul__(self, other)` | `other`: Matrix of equal size → fresh Matrix (mismatch ValueError), or matching Vector → fresh Vector computing vM. Other types return NotImplemented; no scalar multiplication/reflected operators. |
| `Matrix.__imul__(self, other)` | `other`: equal-size Matrix only; replace rows with A*B, synchronize ctypes and return self; mismatch ValueError, unsupported type NotImplemented. |
| `Matrix.__div__(self, other)` | `other`: float only, not int. Componentwise division; fresh Matrix. Exact zero divisor raises ZeroDivisionError; unsupported type returns NotImplemented. |
| `Matrix.__idiv__(self, other)` | `other`: float only, not int. Componentwise division; replace rows/ctypes and return self. Exact zero divisor raises ZeroDivisionError; unsupported type returns NotImplemented. |
| `Matrix.__truediv__ = Matrix.__div__` | Python 3 alias of `__div__`; same signature, scalar domain, result/ownership and zero-divisor behavior. |
| `Matrix.__itruediv__ = Matrix.__idiv__` | Python 3 alias of `__idiv__`; same signature, scalar domain, result/ownership and zero-divisor behavior. |
| `Matrix.det(self)` | Return numeric determinant for sizes 2/3/4 without mutation (despite historical docstring). Other sizes raise NotImplementedError; ordinary cofactor arithmetic, no scaled determinant promise. |
| `Matrix.i_inverse(self)` | Sizes 2/3/4; replace rows/ctypes, return self. Singular inputs raise ZeroDivisionError; other sizes NotImplementedError. Numerical limits below apply. |
| `Matrix.inverse(self)` | Sizes 2/3/4; fresh Matrix, input preserved. Singular inputs raise ZeroDivisionError; other sizes NotImplementedError. Numerical limits below apply. |
| `Matrix.i_transpose(self)` | Swap row/column entries in a square matrix; replace rows/ctypes and return self. |
| `Matrix.transpose(self)` | Swap row/column entries in a square matrix; fresh Matrix, receiver preserved. |

## Numerical domain and ctypes snapshots

For finite Matrix3/4 inputs, binary row-denominator clearing detects singularity
exactly for represented binary64 coefficients before row power-of-two scaling
and floating cofactors. No epsilon or condition threshold rejects small valid
matrices. Approximate uniform scales 1e-300 to 1e300 and selected mixed exponents
have regression evidence. Severe conditioning, cofactor cancellation, extreme
within-row dynamic range and unrepresentable results remain outside accuracy
guarantees; exact nonsingularity does not ensure successful floating inversion.
Rescaling an unrepresentable inverse coefficient yields signed infinity.
Nonfinite inputs use the legacy cofactor path with no new validation policy.
Matrix2 inverse and public determinants do not inherit this scaling protection.

`.c_matrix` uses `gem.common.GLfloat` (ctypes.c_float), with row-major nested
ctypes arrays, not a live binary64 view. Float32 may round or overflow values that
are representable in Python. Reacquire pointers after replacing snapshots and
keep owning arrays alive. See the [installed ctypes example](../getting-started/quick-start.md#ctypes-matrix-export).

```python
from gem.matrix import Matrix, inverse3

m = Matrix(2, [[4.0, 7.0], [2.0, 6.0]])
expected = [[0.6, -0.7], [-0.2, 0.4]]
n = m.inverse()
assert all(abs(n.matrix[i][j] - expected[i][j]) < 1e-14
           for i in range(2) for j in range(2))
for identity in (m * n, n * m):
    assert all(abs(identity.matrix[i][j] - (1 if i == j else 0)) < 1e-14
               for i in range(2) for j in range(2))
small = [[1e-300 if i == j else 0.0 for j in range(3)] for i in range(3)]
large = inverse3(small)
assert all(abs(large[i][i] / 1e300 - 1.0) < 1e-14 for i in range(3))
assert m.i_transpose() is m
assert m.matrix == [[4.0, 2.0], [7.0, 6.0]]
assert list(m.c_matrix[0]) == [4.0, 2.0]
try:
    Matrix(3, [[0.0] * 3 for _ in range(3)]).inverse()
except ZeroDivisionError:
    pass
else:
    raise AssertionError('singular inverse must fail')
```

Transform methods are documented on [transformations](transformations.md),
projection on [projection](projection.md). See [decisions](decisions.md) and [index](index.md).

See the [graphics gallery example](../examples/gallery/transforms.md) for an executable visualization.
