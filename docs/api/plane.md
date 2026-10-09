# Plane representation and operations

Import `Plane`, `flip`, `normalize` from `gem.plane`.
[Source](../../gem/plane.py), [incidence/scaling/polygon tests](../../tests/test_planes.py).

Coefficients satisfy a*x+b*y+c*z+d=0. `.a,.b,.c,.d` are scalar fields;
`.normal` is Vector3 [a,b,c] at the same coefficient scale. The default plane is
all zeros, an unnormalized placeholder, not a valid geometric plane. Direct
field edits can desynchronize normal; construction/normalization synchronize it.
There is no public module constant or plane ctypes export.

Three-point construction uses normalize((b−a)×(c−a)) and d=−n·a; winding sets
orientation. Newell polygon normals wrap last vertex to first; bestFitD is signed
D=mean(n·p), not the equation's d. To construct from a polygon, explicitly use
fromCoeffs(n.x,n.y,n.z,−D). Nonplanar polygons produce an approximation; repeating
the first vertex contributes again to the mean offset.

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `flip(plane)` | `plane`: raw [a,b,c,d,normal] with normal a Vector; return fresh five-element list [−a,−b,−c,−d,−normal], preserving input. |
| `normalize(pdata)` | `pdata`: raw [a,b,c,d]; divide all four by stable magnitude of first three; return four-element tuple. Exact zero normal raises ZeroDivisionError. Inputs preserved; not a scaled coefficient-normalization guarantee for every extreme value. |
| `Plane` | Mutable scalar-coefficient/normal wrapper; see constructor and representation above. |
| `Plane.__init__(self)` | Initialize coefficients to zero and a fresh zero Vector3 normal; return None. |
| `Plane.clone(self)` | Fresh Plane with copied coefficients, independently cloned normal and component list. |
| `Plane.fromCoeffs(self, a, b, c, d)` | Numeric `a,b,c,d`; store scalar coefficients unchanged, replace normal with fresh Vector3 [a,b,c]; mutate receiver, return None. |
| `Plane.fromPoints(self, a, b, c)` | `a,b,c`: noncollinear Vector3 positions; replace coefficients with unit-normal plane, d=−n·a; return None. Inputs preserved; exact degenerate cross product raises ZeroDivisionError. |
| `Plane.i_flip(self)` | Negate all coefficients and replace normal with its negative; return self. |
| `Plane.flip(self)` | Fresh oppositely oriented Plane and normal, preserving receiver; same geometric locus. |
| `Plane.dot(self, vec)` | `vec`: Vector4, compute a*x+b*y+c*z+d*w, numeric result. Supplied w is used; signed geometric distance requires w=1 and unit normal. No mutation. |
| `Plane.i_normalize(self)` | Divide a,b,c,d by same normal magnitude, rebuild normal, return self. Zero normal raises ZeroDivisionError. |
| `Plane.normalize(self)` | Fresh normalized Plane/normal, input preserved; divide all four coefficients together; zero normal raises ZeroDivisionError. |
| `Plane.bestFitNormal(self, vecList)` | `vecList`: ordered Vector3 polygon vertices; fresh unit Newell normal, no receiver/input mutation; wrapped edges and reversed winding respected. Empty/degenerate normal raises ZeroDivisionError. |
| `Plane.bestFitD(self, vecList, bestFitNormal)` | `vecList`: Vector3 vertices; `bestFitNormal`: Vector3; scalar signed mean dot(n,p), preserving inputs/receiver. Unit normal gives geometric offset D; nonunit result scales with n. Empty list raises ZeroDivisionError. |
| `Plane.point_location(self, plane, point)` | Explicit `plane`: Plane (not implicitly self), `point`: XYZ tuple/indexable list. Return 1/−1/0 by exact sign of a*x+b*y+c*z+d; no epsilon. Unordered NaN result prints diagnostic and returns None; no broader nonfinite policy. |

```python
from gem.plane import Plane
from gem.vector import Vector

p = Plane()
assert p.fromCoeffs(0, 0, 2, -4) is None  # z=2, normal scale=2
n = p.normalize()
assert (n.a, n.b, n.c, n.d) == (0.0, 0.0, 1.0, -2.0)
assert n.dot(Vector(4, [1, 3, 2, 1])) == 0
assert p.c == 2 and p.normal.vector == [0, 0, 2]
vertices = [Vector(3, [0, 0, -3]), Vector(3, [1, 0, -3]), Vector(3, [1, 1, -3]), Vector(3, [0, 1, -3])]
normal = p.bestFitNormal(vertices)
D = p.bestFitD(vertices, normal)
polygon = Plane()
polygon.fromCoeffs(normal.vector[0], normal.vector[1], normal.vector[2], -D)
assert D == -3 and polygon.d == 3
assert all(polygon.dot(Vector(4, v.vector + [1])) == 0 for v in vertices)
```

No automatic polygon constructor or distance-query method is added. Wider
malformed/nonfinite policies need separate decisions; see [decisions](decisions.md),
[Ray](ray.md) and [index](index.md).
