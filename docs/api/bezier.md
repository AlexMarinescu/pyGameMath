# Bezier evaluation and adaptive paths

Import from canonical `gem.bezier`.
[Source](../../gem/bezier.py), [evaluation](../../tests/test_bezier.py),
[sampling/builders](../../tests/test_bezier_sampling.py).

For u=1−t, quadratic evaluation is u²p0+2utp1+t²p2; cubic is
u³p0+3u²tp1+3ut²p2+t³p3. Scalar controls produce numeric results; matching Vector
controls produce fresh Vectors, preserving input storage. Evaluation does not
clamp t, restrict Vector dimension to 2/3, or enforce sampling's finite policy.
Generic compatible arithmetic follows the fallback path; arbitrary types are
not guaranteed. No public module constants exist.

| Exact source declaration | Parameters, result and behavior |
|---|---|
| `cubicBezierPoint(t, p0, p1, p2, p3)` | Numeric `t`; `p0,p1,p2,p3` scalar or matching Vector controls; cubic Bernstein result, fresh Vector or numeric value; no mutation/clamping. |
| `quadraticBezierPoint(t, p0, p1, p2)` | Numeric `t`; `p0,p1,p2` scalar or matching Vector controls; quadratic Bernstein result, fresh Vector or numeric value; no mutation/clamping. |
| `BezierPath` | Cubic segment/control and adaptive sampling container; constructor below. |
| `BezierPath.__init__(self)` | Empty control list, curveCount=0; minimum_sqr_distance=0.01, segments_per_curve=10, historical divison_threshold=−0.99. Returns None. |
| `BezierPath.setControlPoints(self, newControlPoints)` | Retain caller `newControlPoints` list, set curveCount=(len−1)//3, return None; no initial shape validation. Caller changes can stale curveCount. |
| `BezierPath.getControlPoints(self)` | Return the exact controlPoints list, not a copy; subsequent caller mutation changes the path. |
| `BezierPath.calculateBezerPoint(self, curveIndex, t)` | `curveIndex`: zero-based cubic index, numeric `t`; evaluate controls[3*i:3*i+4]. Numeric/new Vector result; preserves inputs. Misspelling retained; indexing errors are historical. |
| `BezierPath.interpolate(self, segmentPoints, scale)` | `segmentPoints`: ordered finite scalars or uniform Vector2/3 sources, finite numeric `scale`: tangent length factor. Append fresh generated cubic controls, update count, return None. Fewer than two sources no-op. No input storage mutation; accumulated independent sets need not form a valid connected path. |
| `BezierPath.samplePoints(self, sourcePoints, minSqrDistance, maxSqrDistance, scale)` | Ordered `sourcePoints`, squared thresholds `minSqrDistance,maxSqrDistance`, finite `scale`; rebuild controls via thinning then interpolation, return None. Fewer than two sources no-op before validation; others require finite 0≤min≤max and max>0 (ValueError). Rules below. |
| `BezierPath.getDrawingPoints(self)` | Fresh nested list per segment; include first endpoint, omit shared first endpoint of later segments. Empty controls → []; sampling validation otherwise applies. Does not mutate path/controls. |
| `BezierPath.findDrawingPoints(self, curveIndex)` | `curveIndex`: integer valid segment; fresh ordered sample list with both endpoints; scalar or fresh Vector2/3 points. ValueError for malformed/nonfinite controls or tolerance, IndexError for invalid index. Rules below. |
| `BezierPath.findDrawingPointsAdded(self, curveIndex, t0, t1, pointList, insertionIndex)` | `curveIndex`, interval `t0,t1` with 0≤t0≤t1≤1, caller `pointList` already holding endpoints, `insertionIndex` in [0,len]. Insert fresh ordered interior samples, preserve existing elements, return inserted count. Invalid interval ValueError; insertion/index IndexError; standard sampling validation applies. |

## Fields, tolerance and bounded subdivision

`controlPoints` is mutable control storage; a valid k-segment path has 3k+1 entries.
`curveCount` is the stored segment count maintained by builders/setter.
`minimum_sqr_distance` is a positive finite **squared coordinate-distance**
tolerance; sampling compares geometric distance against its square root.
`segments_per_curve` and misspelled `divison_threshold` are retained historical
fields but do not govern current subdivision.

Sampling validates finite scalar or uniform Vector2/3 controls. It subdivides by
midpoint de Casteljau and tests maximum interior-control distance to the endpoint
**segment**. Clamping the chord projection handles collinear overshoot;
coincident endpoints use distance to that endpoint. Depth is at most 16 per segment,
implemented with an explicit stack; depth-exhausted intervals emit best available
endpoints, which can exceed tolerance. No unconditional floating-point error bound
is claimed. Output follows increasing t, and sampled Vectors/storage are independent.

`samplePoints` preserves first/last source vertices. An interior source is retained
if its squared distance from the last retained vertex reaches minSqrDistance, or
skipping it makes the next source exceed maxSqrDistance from that vertex. These
are thinning heuristics, not gap limits or approximation bounds; sparse sources
can still violate a desired gap. Rebuilding replaces stale generated controls.
`interpolate` remains intentionally append-only; its scale controls endpoint and
interior tangent offsets, not adaptive tolerance.

```python
from gem.bezier import BezierPath, cubicBezierPoint, quadraticBezierPoint
from gem.vector import Vector

controls = [Vector(2, [0, 0]), Vector(2, [1, 2]), Vector(2, [2, 2]), Vector(2, [3, 0])]
assert cubicBezierPoint(0.5, *controls).vector == [1.5, 1.5]
assert quadraticBezierPoint(0.5, 0, 2, 4) == 2
path = BezierPath()
assert path.setControlPoints(controls) is None
assert path.getControlPoints() is controls
path.minimum_sqr_distance = 0.0001
points = path.findDrawingPoints(0)
assert points[0].vector == [0, 0] and points[-1].vector == [3, 0]
assert points[0] is not controls[0] and points[0].vector is not controls[0].vector
source = [Vector(2, [0, 0]), Vector(2, [1, 0]), Vector(2, [2, 0])]
assert path.samplePoints(source, 0, 1, 0.25) is None
first = [p.vector[:] for p in path.controlPoints]
path.samplePoints(source, 0, 1, 0.25)
assert [p.vector for p in path.controlPoints] == first
assert source[1].vector == [1, 0]
```

Transitional imports reexport the same objects; see [compatibility](legacy.md),
[decisions](decisions.md) and [index](index.md).

See the [graphics gallery example](../examples/gallery/bezier.md) for an executable visualization.
