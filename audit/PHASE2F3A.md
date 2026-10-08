# Bezier evaluation and core migration

E01–E03 are corrected in `gem.bezier`, the canonical supported import.
The cubic final term is p3*t³; quadratic and cubic evaluators multiply
controls by scalar weights on the right, supporting gem Vectors without
reflected operators. Evaluation returns fresh Vector storage and preserves
controls. Scalars and same-dimensional Vectors are supported; parameters
remain unclamped, including extrapolation.

`BezierPath` evaluates cubic segments stored as 3*k+1 controls with shared
endpoints. `curveCount` uses integer division, including legacy interpolation.
`setControlPoints` and `getControlPoints` retain caller-owned list storage;
invalid layouts and broader ownership policies are not newly defined.
The historical `calculateBezerPoint` spelling and signatures are retained.

## Migration

```python
from gem.bezier import BezierPath, cubicBezierPoint, quadraticBezierPoint

assert cubicBezierPoint(0.5, 0, 1, 2, 3) == 1.5
path = BezierPath()
path.setControlPoints([0, 1, 2, 3, 4, 5, 6])
assert path.curveCount == 2
assert path.calculateBezerPoint(1, 0.5) == 4.5
```

`gem.experimental.bezier` reexports the same evaluator functions. Its
`BezierPath` is a compatibility subclass that retains legacy interpolation
and adaptive-sampling methods in `_bezier_legacy`; it inherits evaluation
rather than duplicating algorithms. Core `BezierPath` deliberately exposes
only validated evaluation. Legacy adaptive sampling remains affected by
E04 and is not promoted as supported functionality.

Packaging includes `gem.bezier` and the transitional experimental package,
so legacy imports also work in wheels. Other experimental modules become
included in that package but remain unvalidated; no algorithms are moved or
repaired there. No mandatory dependencies or version changes are introduced.

## Verification

Baseline on master 00d002421eb1aed091e143d0126387fb631911df:
**1383 passed, 22 xfailed**. Running experimental defect tests with
`--runxfail` reproduced all seven E01/E02/E03 failures (18 total experimental
failures). E03's original test also called adaptive sampling; its passing
regression now isolates integer count and direct evaluation. A separate
legacy test verifies integer iteration with a sampling stub. Both E04
end-to-end failures remain expected failures.

Final full suite: **1441 passed, 0 failed, 15 xfailed** on Python 3.12.14,
pytest 9.1.1. Exact expected-failure identities match baseline minus the seven
corrected cases: 12 defect cases and three unresolved contract questions.
See `phase2f3a-test-results.json` for identities. Independent de Casteljau
references cover Vector2/3/4, negative and extrapolated parameters, reversal,
degree elevation, ownership, segment boundaries and compatibility imports.

A wheel built with setuptools 84.0.0 / wheel 0.48.0 contains both import
paths. Isolated wheel extraction verifies scalar/vector evaluation and
compatibility function identity. No build artifacts are checked in.

## Remaining migration

Repair E04 separately before promoting sampling. Subsequent phases can
validate and promote Legendre, spherical-harmonics sampling and irradiance
into coherent core modules. Incomplete shadow transport requires a separate
review. Remove `gem.experimental` only in a dedicated final cleanup after
migration and compatibility decisions are resolved.
