# Ordinary and associated Legendre functions

Import `Legendre` from `gem.legendre`.
[Source](../../gem/legendre.py), [independent Rodrigues references through degree 12](../../tests/test_legendre.py).

Supported degree/order are integers 0≤m≤l. Associated values require x in [-1,1];
ordinary m=0 polynomials also evaluate outside that interval. Functions are
**unnormalized**, with Condon–Shortley phase:
P_m^m(x)=(−1)^m(2m−1)!!(1−x²)^(m/2),
P_(m+1)^m=x(2m+1)P_m^m, and
P_l^m=((2l−1)xP_(l−1)^m−(l+m−1)P_(l−2)^m)/(l−m).
SH applies its own normalization; do not substitute a different phase convention.
No module constants exist.

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `Legendre` | Mutable degree/order/argument and scratch container; construction does not evaluate. |
| `Legendre.__init__(self, l, m, x)` | Store numeric `l,m,x` directly; initialize P=1.0, PM1=PML=0.0; return None. No generalized domain validation. |
| `Legendre.mGreaterThan0(self)` | Rebuild scratch P=P_m^m from current m/x; even m=0 resets P to 1; return None. Deterministic repeated initialization. |
| `Legendre.calculatePM1(self)` | Reset P seed, set PM1=x(2m+1)P; return None; mutates P and PM1. |
| `Legendre.calculatePML(self, i)` | `i`: target integer degree ≥m; reset P/PM1 and set PML from recurrence; return None; explicit scratch mutation. Wider invalid-degree behavior is historical. |
| `Legendre.run(self)` | Numeric P_l^m result from local recurrence; preserve P/PM1/PML and l/m/x. Historical l<m returns stored PML rather than raising a new validation exception. |

Fields `.l`, `.m`, `.x`, `.P`, `.PM1`, `.PML` are mutable and not copied from
another polynomial. `run()` does not use helper-mutated state for valid degrees.
Noninteger orders can produce TypeError from range; associated x outside [-1,1]
can raise ValueError from sqrt. Negative orders/degrees, nonfinite inputs,
high-order overflow and generalized invalid-input behavior remain unsupported
rather than covered by a new validation policy.

```python
import math
from gem.legendre import Legendre

p = Legendre(2, 0, 0.5)
assert p.run() == -0.125            # (3*x*x-1)/2
assert (p.P, p.PM1, p.PML) == (1.0, 0.0, 0.0)
assert p.run() == p.run()
a = Legendre(2, 1, 0.5)
assert abs(a.run() + 3 * 0.5 * math.sqrt(0.75)) < 1e-14
assert Legendre(12, 0, 1).run() == 1.0
assert a.calculatePML(2) is None
assert abs(a.PML - a.run()) < 1e-14
```

See [SH normalization](spherical-harmonics.md), [shims](legacy.md),
[decisions](decisions.md) and [index](index.md).
