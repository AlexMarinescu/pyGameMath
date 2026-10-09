# Numerical accuracy and the limits of a result

**Beginner → advanced.** Prerequisites: [vector operations](vectors.md) and,
for inversion, [matrices](transformations.md). Learn to distinguish rounding,
representability and ill-conditioning; choose a meaningful comparison; and read
an exception/fallback in its specific API context. These habits help maintain
CPU reference calculations and diagnose graphics artifacts.

## Absolute and relative comparisons

Binary64 floats approximate most decimal fractions. Algebraically equivalent
expressions can round differently. A useful finite-value test is
abs(actual-reference) ≤ max(abs_tol,rel_tol*max(abs(actual),abs(reference))).
Absolute tolerance sets the scale around zero; relative tolerance scales with
large magnitudes. Choose tolerances from units and expected numerical error,
not from a desire to make a failing test pass. Vector equality stays exact;
this local scalar helper is not an approximate-equality method added to gem.

```python
import math
from gem.vector import Vector
from gem.quaternion import Quaternion

# Local finite-value test, not a new gem API.
def close(a, b, absolute=1e-12, relative=1e-12):
    return abs(a - b) <= max(absolute, relative * max(abs(a), abs(b)))

assert close(1e-15, 0)                   # absolute comparison near zero
assert close(1e8 + 1e-5, 1e8)           # relative comparison at large scale
assert not close(1.01, 1.0)
assert Vector(2, [0.1 + 0.2, 0]) != Vector(2, [0.3, 0])
large = Vector(2, [1e300, 1e300])
tiny = Vector(2, [1e-300, 1e-300])
assert close(large.magnitude() / 1e300, math.sqrt(2))
assert close(tiny.magnitude() / 1e-300, math.sqrt(2))
assert all(close(x, math.sqrt(0.5)) for x in large.normalize().vector)
assert all(close(x, math.sqrt(0.5)) for x in tiny.normalize().vector)
assert large.vector == [1e300, 1e300]
assert Vector(3).normalize().vector == [0.0, 0.0, 0.0]
assert Quaternion([0, 0, 0, 0]).normalize().data == [1.0, 0.0, 0.0, 0.0]
assert math.isinf(Quaternion([1e308] * 4).magnitude())  # true length exceeds binary64
assert Quaternion([1e308] * 4).normalize().data == [0.5, 0.5, 0.5, 0.5]
```

The naive square sum would overflow for 1e300 and underflow for 1e-300.
Scaled hypot-style length/normalization preserves these directions. A true norm
outside binary64 can still be infinity; normalization scales components first,
so it can succeed without representing that norm. These guarantees apply to
supported finite norms, not automatically to dot/cross products or every algorithm.

Direct zero normalization returns zero Vector/identity Quaternion, but geometric
callers reject zero directions/axes or plane normals with ZeroDivisionError.
There is no epsilon-based zero detection. A zero vector has no geometric angle,
and the fallback does not justify treating it as a usable direction.

## Singularity is different from conditioning

A singular matrix has no inverse. A near-singular matrix may have one but can
amplify small data/rounding errors enormously. Matrix3/4 inversion uses exact
represented-coefficient singularity detection, power-of-two row scaling and
floating cofactors. There is no arbitrary near-singular cutoff. Row scaling
reduces avoidable range errors; it does not make a severely conditioned problem
accurate or guarantee every mixed-scale cofactor computation succeeds.

```python
from gem.matrix import Matrix, inverse3, identity
from gem.vector import Vector

small = [[1e-300 if i == j else 0.0 for j in range(3)] for i in range(3)]
inv = inverse3(small)
assert all(abs(inv[i][i] / 1e300 - 1) < 1e-14 for i in range(3))
# Independent diagonal answer; verify both product orders as additional evidence.
a = Matrix(3, small)
b = Matrix(3, inv)
assert all(abs(product.matrix[i][j] - identity(3)[i][j]) < 1e-14
           for product in (a * b, b * a) for i in range(3) for j in range(3))
near_singular = Matrix(3, [[1, 0, 0], [0, 1, 0], [0, 0, 1e-12]])
assert abs(near_singular.inverse().matrix[2][2] / 1e12 - 1) < 1e-14
# A change of only 1e-12 in the last equation changes its solution by 1.
x0 = near_singular.inverse() * Vector(3, [0, 0, 0])
x1 = near_singular.inverse() * Vector(3, [0, 0, 1e-12])
assert x0.vector == [0.0, 0.0, 0.0]
assert all(abs(a - b) < 1e-14 for a, b in zip(x1.vector, [0, 0, 1]))
try:
    Matrix(3, [[1, 0, 0], [0, 1, 0], [0, 0, 0]]).inverse()
except ZeroDivisionError:
    pass
else:
    raise AssertionError('singular input must raise')
```

The small diagonal inverse is exactly 1e300 mathematically; tolerances allow
binary rounding. The third equation in the second matrix is 1e-12*x3=b3, so a
tiny absolute perturbation in b3 produces a unit change in x3. That is input
conditioning, not a scaling-algorithm defect. Binary32 ctypes snapshots may lose
precision or overflow even when Python rows are finite; do not use a float32
snapshot as the oracle for an extreme binary64 inverse.

## Domain limits and regression evidence

Matrix2 inversion and public determinant functions do not inherit Matrix3/4
scaling protection. Quaternion inverse uses a direct squared norm, so its extreme
range behavior differs from Quaternion normalization. Signed infinity is permitted
when inverse rescaling cannot represent a genuine coefficient. General nonfinite
inputs retain API-specific behavior; newer SH APIs reject checked nonfinite
coefficients, but there is no library-wide NaN/Infinity exception policy. The
finite comparison helper above is not a valid infinity/NaN predicate.

Independent known answers, both inverse product orders, invariants and ownership
tests complement each other. Round trips alone can conceal coupled defects;
plausible images alone cannot validate linear lighting coefficients. Exact equality
is useful for exact contractual values, while tolerance tests belong to documented
numerical cases. Do not weaken existing tolerances to hide a regression.

Further evidence: [stable-norm tests](../../tests/test_numerical_robustness.py),
[inversion tests](../../tests/test_inverse_optimization.py), [SH reference checks](../../tests/test_hdr_sh_example.py).
See [Vector API](../api/vector.md), [Matrix limits](../api/matrix.md),
[Quaternion domains](../api/quaternion.md), [open policies](../api/decisions.md),
[lighting](lighting.md), [interoperability](interop.md) and [index](index.md).
