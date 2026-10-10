# Bezier curves as motion paths

**Intermediate.** Prerequisites: [vectors](vectors.md), polynomial weights and
basic differentiation. Learn how controls shape a curve, how to derive a tangent
locally and how to sample an approximation without promising constant speed.
Use this for camera rails, trails and procedural paths rather than treating it
as a ready-made animation timing system.

## Controls, points and derivatives

With $u=1-t$, the quadratic and cubic forms are

$$
B_2(t)=u^2p_0+2utp_1+t^2p_2,\qquad
B_3(t)=u^3p_0+3u^2tp_1+3ut^2p_2+t^3p_3.
$$

Endpoints are p0/p3; interior controls pull the
curve and determine endpoint tangents but normally are not points on the curve.
For a cubic, the derivative is

$$
B'_3(t)=3u^2(p_1-p_0)+6ut(p_2-p_1)+3t^2(p_3-p_2).
$$

There is no public tangent helper: the local function below is this derivative
written with existing Vector arithmetic, not a new gem API.

Choose p0=[0,0,0], p1=[1,2,0], p2=[2,2,0], p3=[3,0,0]. Then
B(t)=[3t,6t(1-t),0] and B′(t)=[3,6-12t,0]. The explicit polynomials give
independent point/tangent checks. Increasing t advances X monotonically, yet
speed ||B′(t)|| changes from sqrt(45) at t=0 to 3 at t=1/2.

## Evaluate and sample a motion path

```python
# Control points shape the curve; equal parameter steps need not have equal length.
import math
from gem.bezier import BezierPath, cubicBezierPoint, quadraticBezierPoint
from gem.vector import Vector

controls = [Vector(3, [0, 0, 0]), Vector(3, [1, 2, 0]), Vector(3, [2, 2, 0]), Vector(3, [3, 0, 0])]

def tangent(t, points):
    u = 1 - t
    return ((points[1] - points[0]) * (3 * u * u)
            + (points[2] - points[1]) * (6 * u * t)
            + (points[3] - points[2]) * (3 * t * t))

times = [0, 0.25, 0.5, 0.75, 1]
motion = [cubicBezierPoint(t, *controls) for t in times]
assert all(p.vector == [3 * t, 6 * t * (1 - t), 0] for t, p in zip(times, motion))
assert tangent(0.5, controls).vector == [3.0, 0.0, 0.0]
assert tangent(0.5, controls).normalize().vector == [1.0, 0.0, 0.0]
assert abs(tangent(0, controls).magnitude() - math.sqrt(45)) < 1e-14
assert quadraticBezierPoint(0.5, 0, 2, 4) == 2
path = BezierPath()
path.setControlPoints(controls)
path.minimum_sqr_distance = 0.0025  # squared tolerance: 0.05 coordinate units
samples = path.findDrawingPoints(0)
assert samples[0].vector == [0, 0, 0] and samples[-1].vector == [3, 0, 0]
assert all(a.vector[0] < b.vector[0] for a, b in zip(samples, samples[1:]))
polyline_length = sum((b - a).magnitude() for a, b in zip(samples, samples[1:]))
assert polyline_length > 3.0
assert path.getControlPoints() is controls
assert samples[0].vector is not controls[0].vector
assert controls[1].vector == [1, 2, 0]
print([point.vector for point in motion])
```

Output: `[[0, 0, 0], [0.75, 1.125, 0.0], [1.5, 1.5, 0.0], [2.25, 1.125, 0.0], [3, 0, 0]]`
(numeric zero formatting may include `.0`). The midpoint tangent is +X, but the
endpoint tangent tilts upward. An orientation built from a tangent additionally
needs an up/frame policy; this guide does not invent one. A zero derivative has
no unique direction even though direct Vector zero normalization returns zero.

## Parameter time, geometric distance and approximation

Equal t steps are parameter-space samples, not equal traveled distances. Arc
length is integral_0^t ||B′(s)|| ds. For approximate distance-based traversal, one
can build cumulative chord lengths from sampled points and invert that table;
this would be a polyline approximation with its own timing/error policy, not
constant-speed Bezier traversal supplied by gem. No such traversal is claimed here.

Adaptive midpoint de Casteljau sampling tests interior-control distance to the
endpoint **segment**, supporting finite scalars and uniform Vector2/3 controls.
Squared tolerance controls geometric flatness, not a time step. Recursion-equivalent
depth is capped at 16; difficult intervals emit best available endpoints and may
exceed tolerance. Coincident endpoints and collinear overshoot are handled; lower
tolerance often yields more samples but is not an unconditional numerical guarantee.

Cubic paths use 3k+1 controls. `getDrawingPoints` returns nested per-segment lists,
omitting duplicate shared endpoints after the first. Controls are retained by
reference; sampled Vectors are fresh. `interpolate` intentionally appends generated
controls, whereas `samplePoints` replaces generated controls on each build using
separate min/max squared-distance thinning heuristics. These are not strict spacing
limits and do not replace adaptive flatness tolerance. `segments_per_curve` and
`divison_threshold` are historical fields, not active subdivision controls.

See [Bezier API](../api/bezier.md), [sampling regressions](../../tests/test_bezier_sampling.py),
[ownership](../architecture/conventions.md#curves-polynomials-and-sampling) and
[numerical accuracy](numerical.md). Continue with [quaternion orientation](quaternions.md),
[geometry](geometry.md) or the [index](index.md).

See the [visual example](../examples/gallery/bezier.md) and its reproducible assets.
