# Quaternion API

Import helpers from `gem.quaternion` and vectors from `gem.vector`.
Components use `[w,x,y,z]`, with identity `[1,0,0,0]` and Hamilton
multiplication. Quaternion/matrix conversions preserve the library's
row-vector convention: positive rotation about +Z sends +X toward +Y.
Unit q and -q represent the same rotation. `q1*q2` applies q2 first; the
corresponding row matrices compose in reverse order.

## Constructing and applying rotations

Use `quat_from_axis_angle(axis,theta)` to construct a rotation Quaternion.
Theta is in **degrees**. A finite nonzero three-component axis may be a
Python list or Vector3; temporary normalization preserves its original
values and storage.

`quat_rotate_from_axis_angle(axis,theta)` is a **legacy helper**. It rotates
its normalized axis about itself and returns the resulting pure Quaternion,
approximately `[0,normalized_axis]`. A rotation about its own axis leaves
that axis unchanged, so the result does not encode theta as an orientation.
Its original numerical sandwich and public signature are retained, including
floating-point roundoff. It is not an axis-angle rotation constructor.

```python
from gem.quaternion import (
    quat_from_axis_angle, quat_rotate_from_axis_angle, quat_rotate_vector,
)
from gem.vector import Vector

axis = [0, 0, 2]
legacy = quat_rotate_from_axis_angle(axis, 90)  # approximately [0,0,0,1]
rotation = quat_from_axis_angle(axis, 90)      # [sqrt(.5),0,0,sqrt(.5)]
point = Vector(3, [1, 0, 0])
rotated = quat_rotate_vector(rotation, point) # approximately [0,1,0]
# axis and point retain their values and storage.
```

Applying `legacy` as an orientation would produce a half-turn about Z,
even for theta=0. Construct rotations with `quat_from_axis_angle` instead.
No runtime deprecation warning or additional constructor is needed.

## Helper representations and units

| API | Inputs / units | Result |
| --- | --- | --- |
| `Quaternion(data)` | Four component values [w,x,y,z] | Quaternion retaining supplied storage |
| `quat_from_axis_angle(axis,theta)` | Vector3/list axis; degrees | Rotation Quaternion |
| `quat_rotate_from_axis_angle(axis,theta)` | Vector3/list axis; degrees | Legacy pure Quaternion axis result |
| `quat_rotate_x/y/z_from_angle(theta)` | Radians | Four-element component list |
| `quat_rotate(origin,axis,theta)` | Vector3 origin; raw unit axis coordinates; degrees | Rotated Vector3 |
| `quat_rotate_vector(quat,vec)` | Unit Quaternion, Vector3 | Rotated Vector3 |
| `Quaternion.toMatrix()` | Unit Quaternion | Matrix4 with synchronized ctypes export |
| `quat_from_matrix(matrix)` | Matrix3/Matrix4 proper rotation block | Quaternion, with sign equivalence |
| `Quaternion.log()` / `quat_log(quat)` | Unit Quaternion | Fresh four-element list |

The X/Y/Z names in the table abbreviate three separate existing helpers;
they are not a new combined function. `quat_rotate` retains raw axis
handling and does not normalize that axis. `q*Vector` is the Hamilton
product with the pure quaternion `(0,v)`, returning a Quaternion; it is
not vector rotation. Rotation uses `q*(0,v)*conjugate(q)`.

Quaternion `getForward()` uses +Z, `getBack()` uses -Z, and right/up use
+X/+Y. Vector `front()` uses -Z; this separate established convention is
preserved.

## Ownership and domains

`Quaternion(data)` retains the supplied mutable component list; make an
explicit list copy when isolation is needed. Direct `.data` edits are
observable through that shared storage. Returning arithmetic, interpolation,
conjugation, inversion, powers and logarithms produce independent storage.
In-place operations retain receiver identity, replace its component list,
and return self. Axis-angle helpers preserve caller axes. Vector rotation
preserves its input Vector. Scalar multiplication/division wrappers accept
floats, including float subclasses; other scalar protocols are not added.

Rotation and proper-rotation matrix conversion require unit inputs; they do
not normalize implicitly. A nonunit quaternion sandwich scales a Vector by
the quaternion's squared magnitude. Nonunit/zero matrix conversion retains
historical arithmetic, which is not a general rotation guarantee.

Inverse supports ordinary nonzero general quaternions as conjugate divided
by squared norm. Zero inverse raises ZeroDivisionError. Direct zero-quaternion normalization returns a fresh identity Quaternion;
in-place normalization preserves receiver identity. Finite norms and
normalization use stable scaled hypot calculations. A finite input norm
beyond binary64 range may be infinity, while scaled normalization still
produces a finite unit direction. Geometric zero-axis callers retain
ZeroDivisionError through exact-zero guards. NaN/Inf input arithmetic
retains its historical path, without a new policy.
Unsupported axis representations in the two axis-angle helpers return
NotImplemented; zero-axis normalization raises ZeroDivisionError. Shape,
nonfinite, near-degenerate and extreme-scale validation policies remain
separate questions.

## Powers, logarithms and interpolation

Powers and logarithms support unit quaternions without implicit normalization
or a norm tolerance. Their principal angle is `atan2(|imaginary|,w)` in
[0,pi]. Powers support finite real exponents and return fresh Quaternions;
log returns `[0,axis*angle]`. Zero inputs raise ValueError. Identity log is
four zeros. Negative identity supports integer powers by parity, but its
fractional powers and logarithm raise ValueError because there is no unique
axis. Quaternion signs are not canonicalized.

SLERP is accurate shortest-path spherical interpolation for unit inputs,
with norm and known-axis component errors checked within 1e-12 on [0,1].
At exactly zero endpoint dot, supplied signs select between equally short
half-turn paths. Parameters are not clamped; the accuracy guarantee does
not extend to arbitrary extrapolation. LERP remains unnormalized.

Legacy `squad(q1,q2,t)` uses the receiver as start, q2 as end, and q1 as
additional blend control. It uses sign-sensitive `slerp_no_invert`, whose
linear branches may produce nonunit or antipodally degenerate results.
`squad4(q0,q1,s0,s1,t)` is separate conventional SQUAD using accurate SLERP;
s0/s1 are unit SQUAD controls, not simply neighbouring keyframes. See
[formulas and examples](../audit/CONVENTIONS.md#quaternion-interpolation).

There is no quaternion cross-product or exponential API. Exponential and
SQUAD intermediate-control generation are not provided. The remaining
numerical policies do not reopen the established rotation, ownership,
power/logarithm or interpolation conventions.
