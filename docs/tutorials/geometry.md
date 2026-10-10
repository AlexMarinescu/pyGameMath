# Planes, rays and a local picking query

**Intermediate.** Prerequisites: [dot products](vectors.md) and, for screen picking,
[camera unprojection](camera.md). Derive ray-plane intersection and interpret hits
without treating Ray's historical state as a collision API.

## Solve the plane equation along a ray

Plane coefficients satisfy n·p+d=0, where n=[a,b,c]. A ray's mathematical half-line
is p(t)=o+t*u, t≥0. Substitution gives t=-(n·o+d)/(n·u) when denominator is nonzero.
For unit u, t measures geometric distance. Multiplying all plane coefficients by
the same nonzero scalar leaves t unchanged. Only a unit plane normal makes its
position evaluation a signed distance; normalization must scale d as well.

Ray construction retains start/direction Vector references and normalizes the
caller direction in place. Its `.distance` stores the original direction length;
it is not automatically a hit distance or an interval bound. `.end` is an
intersection placeholder/state with no validity flag. The local query below
returns a separate status, parameter and hit point and does not mutate that state.

An exactly zero denominator means parallel; a zero numerator too means the whole
line lies in the plane, not a unique hit. t<0 lies behind the ray origin. This
example uses exact-zero comparisons and an unbounded forward half-line. Choosing
a near-parallel tolerance or finite range is an application-specific numerical
policy, not an implicit rule added to gem.

## Local ray-plane calculation

```python
# Use the plane equation n dot p + d = 0 for the independent query.
from gem.plane import Plane
from gem.ray import Ray
from gem.vector import Vector

# Tutorial-local calculation, not a gem intersection method.
def query_plane(ray, plane):
    value = plane.dot(Vector(4, ray.start.vector + [1.0]))
    denominator = plane.normal.dot(ray.dir)
    if denominator == 0.0:
        return ('coplanar' if value == 0.0 else 'parallel'), None, None
    t = -value / denominator
    if t < 0.0:
        return 'behind', t, None
    return 'hit', t, ray.start + ray.dir * t

plane = Plane()
plane.fromCoeffs(0, 0, 2, -2)  # z=1, deliberately nonunit normal
start = Vector(3, [0, 0, 3])
direction = Vector(3, [0, 0, -2])
ray = Ray(start, direction)
status, t, hit = query_plane(ray, plane)
assert status == 'hit' and t == 2.0 and hit.vector == [0.0, 0.0, 1.0]
assert plane.dot(Vector(4, hit.vector + [1])) == 0
assert ray.distance == 2 and ray.end.vector == [0.0, 0.0, 0.0]
assert direction.vector == [0.0, 0.0, -1.0] and ray.start is start
assert query_plane(ray, plane.normalize())[1] == t
parallel = Ray(Vector(3, [0, 0, 3]), Vector(3, [1, 0, 0]))
coplanar = Ray(Vector(3, [0, 0, 1]), Vector(3, [1, 0, 0]))
behind = Ray(Vector(3, [0, 0, 0]), Vector(3, [0, 0, -1]))
assert query_plane(parallel, plane)[0] == 'parallel'
assert query_plane(coplanar, plane)[0] == 'coplanar'
assert query_plane(behind, plane)[0] == 'behind'
copy = ray.duplicate()
assert copy.start.vector is not ray.start.vector and copy.end is not ray.end
print(status, t, hit.vector)
```

Output: `hit 2.0 [0.0, 0.0, 1.0]`. Plane evaluation starts at 4 and decreases by
2 per unit ray distance, giving t=2. The constructor changed the supplied
Vector direction, while duplication copies every Vector and its storage without
rerunning normalization. Coefficients/normal must remain synchronized: arbitrary
field edits can make this local query inconsistent.

## Build a picking ray

Unproject the same window XY at depth 0 and 1 to obtain near/far positions. Their
difference gives a direction through that screen location; this example chooses
the near-plane position as ray origin. It is a local construction, not an exposed
picking helper or general scene intersection system.

```python
from gem.matrix import lookAt, perspective, unproject
from gem.ray import Ray
from gem.vector import Vector

view = lookAt(Vector(3, [0, 0, 5]), Vector(3, [0, 0, 0]), Vector(3, [0, 1, 0]))
projection = perspective(90, 2, 1, 10)
viewport = [10, 20, 200, 100]
near = unproject(110, 70, 0, view, projection, viewport)
far = unproject(110, 70, 1, view, projection, viewport)
assert all(abs(x - e) < 1e-13 for x, e in zip(near.vector, [0, 0, 4]))
assert all(abs(x - e) < 1e-13 for x, e in zip(far.vector, [0, 0, -5]))
picking = Ray(near, far - near)
assert all(abs(x - e) < 1e-14 for x, e in zip(picking.dir.vector, [0, 0, -1]))
assert abs(picking.distance - 9) < 1e-13
point = picking.start + picking.dir * 4.0
assert all(abs(x) < 1e-13 for x in point.vector)
assert picking.end.vector == [0.0, 0.0, 0.0]
```

The center pixel line travels along -Z and reaches the world origin four units
after the near plane. Its stored distance happens to be the near/far separation;
an application must explicitly enforce any finite query interval. Unprojection's
zero-W sentinel cannot itself identify a valid picking ray; degenerate direction
construction rejects zero with ZeroDivisionError.

## Limits and uses

Ray/Plane support placement, construction and rigid transforms, not AABB/triangle
intersection, collision acceleration or scene visibility. Ray transforms rotate
about the coordinate origin, preserve distance and leave end untouched. Do not
set end automatically just because the query returned a point; choose application
hit ownership/validity explicitly. Near-parallel division, coefficient scale and
ill-conditioned projection inverses can amplify errors; no broad epsilon policy
is established. Normalizing a zero-normal plane or constructing a plane from
collinear points raises ZeroDivisionError; fromCoeffs itself does not enforce a
nonzero normal. The local query assumes a valid plane with synchronized fields.
Nonfinite and malformed inputs are outside this example's domain.

See [Plane API](../api/plane.md), [Ray API](../api/ray.md),
[hit-state decisions](../api/decisions.md), [camera](camera.md),
[numerical accuracy](numerical.md) and [tutorial index](index.md).
