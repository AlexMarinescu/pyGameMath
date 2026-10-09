# Displacement, facing and movement

**Beginner.** Prerequisites: coordinates, addition and square roots;
[Vector quick start](../getting-started/quick-start.md#vectors-and-ownership).
Learn to find a target direction, measure its angle from an object's forward axis,
and move without overshooting. This is useful for steering targets, aim cones and
simple position updates; it does not implement a physics or navigation system.

## From two positions to a direction

For position p and target a, displacement v=a−p contains both direction and length.
Distance is ||v|| and unit direction is u=v/||v|| for nonzero v. Moving with speed
s over timestep dt gives p′=p+u*min(s*dt,||v||). The minimum prevents passing a
stationary target in one step; it is a locally derived movement rule, not a gem
steering API. Speed and dt must be nonnegative in this example.

For a unit forward direction f, f·u=cos(theta). Positive dot means the target is
in the forward hemisphere; it does not mean the vectors are collinear. To compute
an angle, clamp the floating-point dot into [-1,1] before acos. This clamp protects
a local trigonometric calculation from rounding, not a change to Vector equality.

Cross products encode orientation: X×Y=Z. Their magnitude is ||a||||b||sin(theta),
so parallel inputs give zero rather than a usable orientation axis. Coordinate
labels do not choose a universal world frame: Vector.front is -Z, while identity
Quaternion.getForward is +Z. Here we explicitly choose forward=+X.

## Worked target update

```python
import math
from gem.vector import Vector, cross

position = Vector(3, [1.0, 2.0, 0.0])
target = Vector(3, [4.0, 6.0, 0.0])
forward = Vector(3, [1.0, 0.0, 0.0])
displacement = target - position
remaining = displacement.magnitude()
assert displacement.vector == [3.0, 4.0, 0.0] and remaining == 5.0
direction = displacement.normalize()
cosine = forward.dot(direction)
angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
in_front = cosine > 0.0
speed, dt = 2.0, 0.5
moved = position + direction * min(speed * dt, remaining)
assert in_front and abs(angle - 53.13010235415598) < 1e-12
assert all(abs(a - b) < 1e-14 for a, b in zip(moved.vector, [1.6, 2.8, 0.0]))
assert position.vector == [1.0, 2.0, 0.0] and target.vector == [4.0, 6.0, 0.0]
assert cross(forward, Vector(3, [0, 1, 0])).vector == [0, 0, 1]
arrived = target + Vector(3).normalize() * 1.0
assert arrived == target
print(in_front, round(angle, 6), [round(x, 6) for x in moved.vector])
```

Output: `True 53.130102 [1.6, 2.8, 0.0]`. The 3–4–5 displacement gives a unit
[0.6,0.8,0] direction; one unit of movement adds [0.6,0.8,0]. If already at the
target, zero normalization returns zero and no motion occurs; an angle to that
zero displacement is undefined, so handle arrival before asking for a facing angle.
Returning arithmetic/normalize preserves source storage. Constructors can retain
caller lists; see [Vector ownership](../api/vector.md).

## Limits and next steps

Stable magnitude/normalization handle large/tiny finite components without naive
square-sum overflow. Dot/cross products do not inherit a blanket extreme-value
guarantee. Floating results are approximate; exact Vector equality is appropriate
for deliberately exact values, not for comparing independently rounded directions.
No uniform malformed-dimension or NaN/Infinity policy is established. See the
[numerical guide](numerical.md) before widening inputs.

An aim cone uses dot≥cos(half-angle) for unit directions, avoiding acos when only
a yes/no answer is needed. Real steering additionally needs acceleration, timestep
and obstacle policies. Continue with [object transforms](transformations.md) or
[quaternion orientation](quaternions.md); consult [Vector API](../api/vector.md),
[common angle utilities](../api/common.md) and [tutorial index](index.md).
