"""Quadratic and cubic Bezier evaluation for scalars and gem Vectors.

Parameters are not clamped. BezierPath evaluates cubic segments with a
shared endpoint layout of 3*k+1 control points. Adaptive sampling uses
squared coordinate-distance tolerance and a maximum subdivision depth of 16.
"""

import math
from gem.vector import Vector


# Binary64 parameter ranges in which the final Bernstein weight can become
# subnormal. The cubic bound conservatively covers t**3 < 2**-1022.
_QUADRATIC_SMALL_PARAMETER = math.ldexp(1.0, -511)
_CUBIC_SMALL_PARAMETER = math.ldexp(1.0, -340)


def _finite_scalar_controls(values):
    try:
        return all(type(value) in (int, float) and math.isfinite(value)
                   for value in values)
    except OverflowError:
        # Huge integers retain the existing arithmetic/conversion behavior.
        return False


def _small_parameter_component(t, values):
    # Here 1-t rounds to exactly 1. Apply factors smaller than one to the
    # control first, so a lost t**2/t**3 weight cannot erase a finite term.
    if len(values) == 3:
        return values[0] + values[1]*(2*t) + (values[2]*t)*t
    return (values[0] + values[1]*(3*t) + (values[2]*(3*t))*t
            + ((values[3]*t)*t)*t)


def _small_parameter_point(t, controls):
    if all(type(point) is Vector for point in controls):
        dimension = controls[0].size
        if all(point.size == dimension for point in controls):
            rows = [point.vector for point in controls]
            if all(_finite_scalar_controls(row) for row in rows):
                return Vector(dimension, [_small_parameter_component(
                    t, [row[i] for row in rows]) for i in range(dimension)])
    elif _finite_scalar_controls(controls):
        return _small_parameter_component(t, controls)
    # Custom arithmetic, mismatched dimensions and nonfinite inputs retain
    # the original dispatch and operand order.
    return None


def cubicBezierPoint(t, p0, p1, p2, p3):
    """Evaluate the cubic Bernstein polynomial without changing controls."""
    if type(t) is float and 0 < abs(t) < _CUBIC_SMALL_PARAMETER:
        result = _small_parameter_point(t, (p0, p1, p2, p3))
        if result is not None:
            return result
    u = 1 - t
    if (type(p0) is Vector and type(p1) is Vector and type(p2) is Vector
            and type(p3) is Vector and p0.size == p1.size == p2.size == p3.size
            and type(t) in (int, float)):
        a, b, c, d = p0.vector, p1.vector, p2.vector, p3.vector
        w0, w1, w2, w3 = u * u * u, 3 * u * u * t, 3 * u * t * t, t * t * t
        return Vector(p0.size, [(a[i] * w0 + b[i] * w1 + c[i] * w2 + d[i] * w3)
                                for i in range(p0.size)])
    return (p0 * (u * u * u) + p1 * (3 * u * u * t)
            + p2 * (3 * u * t * t) + p3 * (t * t * t))


def quadraticBezierPoint(t, p0, p1, p2):
    """Evaluate the quadratic Bernstein polynomial."""
    if type(t) is float and 0 < abs(t) < _QUADRATIC_SMALL_PARAMETER:
        result = _small_parameter_point(t, (p0, p1, p2))
        if result is not None:
            return result
    u = 1 - t
    if (type(p0) is Vector and type(p1) is Vector and type(p2) is Vector
            and p0.size == p1.size == p2.size and type(t) in (int, float)):
        a, b, c = p0.vector, p1.vector, p2.vector
        w0, w1, w2 = u * u, 2 * u * t, t * t
        return Vector(p0.size, [a[i] * w0 + b[i] * w1 + c[i] * w2
                                for i in range(p0.size)])
    return p0 * (u * u) + p1 * (2 * u * t) + p2 * (t * t)


class BezierPath(object):
    """Evaluate and sample cubic segments; explicit control lists remain caller-owned.

    Use setControlPoints with 3*k+1 points for k complete cubic segments.
    The historical calculateBezerPoint spelling is preserved.
    """
    def __init__(self):
        self.segments_per_curve = 10
        self.minimum_sqr_distance = 0.01
        self.divison_threshold = -0.99
        self.controlPoints = []
        self.curveCount = 0

    def setControlPoints(self, newControlPoints):
        self.controlPoints = newControlPoints
        self.curveCount = (len(self.controlPoints) - 1) // 3

    def getControlPoints(self):
        return self.controlPoints

    def calculateBezerPoint(self, curveIndex, t):
        nodeIndex = curveIndex * 3
        return cubicBezierPoint(t, self.controlPoints[nodeIndex],
                                self.controlPoints[nodeIndex + 1],
                                self.controlPoints[nodeIndex + 2],
                                self.controlPoints[nodeIndex + 3])

    def interpolate(self, segmentPoints, scale):
        """Append cubic controls using the historical tangent construction.

        Existing controls are retained; this is not a path replacement API.
        Generated points are fresh and caller-owned lists are not extended.
        Fewer than two source points leave the path unchanged.
        """
        if len(segmentPoints) < 2:
            return
        coords, dimension = _coordinates(segmentPoints)
        if not math.isfinite(scale):
            raise ValueError("scale must be finite")
        generated = []
        for index, point in enumerate(coords):
            if index == 0:
                tangent = tuple(b-a for a,b in zip(point,coords[1]))
                generated.extend([point, tuple(a+scale*b for a,b in zip(point,tangent))])
            elif index == len(coords)-1:
                tangent = tuple(b-a for a,b in zip(coords[index-1],point))
                generated.extend([tuple(a-scale*b for a,b in zip(point,tangent)),point])
            else:
                previous, following = coords[index-1], coords[index+1]
                tangent = tuple(b-a for a,b in zip(previous,following))
                magnitude = math.hypot(*tangent)
                tangent = tuple(a/magnitude for a in tangent) if magnitude else tuple(0 for a in tangent)
                before = math.hypot(*(a-b for a,b in zip(point,previous)))
                after = math.hypot(*(a-b for a,b in zip(following,point)))
                generated.extend([tuple(a-scale*before*b for a,b in zip(point,tangent)),
                                  point, tuple(a+scale*after*b for a,b in zip(point,tangent))])
        self.controlPoints = list(self.controlPoints) + [_point(p,dimension) for p in generated]
        self.curveCount = (len(self.controlPoints)-1)//3

    def samplePoints(self, sourcePoints, minSqrDistance, maxSqrDistance, scale):
        """Thin ordered source vertices and replace the generated cubic path.

        Thresholds are squared coordinate distances, not spacing or error
        guarantees. Retain endpoints. Fewer than two points are a no-op.
        """
        if len(sourcePoints) < 2:
            return
        coords, _ = _coordinates(sourcePoints)
        if not (math.isfinite(minSqrDistance) and math.isfinite(maxSqrDistance)
                and 0 <= minSqrDistance <= maxSqrDistance and maxSqrDistance > 0):
            raise ValueError("require finite 0 <= minSqrDistance <= maxSqrDistance and max > 0")
        retained = [0]
        for index in range(1,len(coords)-1):
            last = coords[retained[-1]]
            distance_squared = sum((a-b)*(a-b) for a,b in zip(coords[index],last))
            next_squared = sum((a-b)*(a-b) for a,b in zip(coords[index+1],last))
            if distance_squared >= minSqrDistance or next_squared > maxSqrDistance:
                retained.append(index)
        retained.append(len(coords)-1)
        generated = BezierPath()
        generated.interpolate([sourcePoints[i] for i in retained], scale)
        self.controlPoints = generated.controlPoints
        self.curveCount = generated.curveCount

    def getDrawingPoints(self):
        """Return per-curve lists, omitting repeated shared boundary samples."""
        if not self.controlPoints:
            return []
        self._sampling_controls(0)
        result = []
        for index in range(self.curveCount):
            points = self.findDrawingPoints(index)
            result.append(points if index == 0 else points[1:])
        return result

    def findDrawingPoints(self, curveIndex):
        """Sample a cubic in increasing parameter order, including endpoints.

        minimum_sqr_distance is a positive finite squared distance. The
        control-to-chord segment criterion is conservative before the depth
        cap; depth-limited output may exceed tolerance. No absolute numerical
        error guarantee is made.
        """
        controls, dimension = self._sampling_controls(curveIndex)
        samples = _subdivide(controls, self._sampling_tolerance())
        return [_point(coords, dimension) for coords in samples]

    def findDrawingPointsAdded(self, curveIndex, t0, t1, pointList, insertionIndex):
        """Insert ordered interior samples for [t0,t1]; return their count.

        pointList already contains interval endpoints. Existing list elements
        are retained; only freshly allocated interior points are inserted.
        """
        controls, dimension = self._sampling_controls(curveIndex)
        if not (0 <= t0 <= t1 <= 1):
            raise ValueError("sampling interval must satisfy 0 <= t0 <= t1 <= 1")
        if not 0 <= insertionIndex <= len(pointList):
            raise IndexError("sampling insertion index out of range")
        # Restrict the Bernstein polynomial to the requested interval.
        if t1 < 1:
            controls, _ = _split(controls, t1)
        if t0 > 0 and t1 > 0:
            _, controls = _split(controls, t0 / t1)
        samples = _subdivide(controls, self._sampling_tolerance())[1:-1]
        points = [_point(coords, dimension) for coords in samples]
        pointList[insertionIndex:insertionIndex] = points
        return len(points)

    def _sampling_tolerance(self):
        value = self.minimum_sqr_distance
        if not math.isfinite(value) or value <= 0:
            raise ValueError("minimum_sqr_distance must be positive and finite")
        return math.sqrt(value)

    def _sampling_controls(self, curveIndex):
        count = len(self.controlPoints)
        if count < 4 or (count - 1) % 3:
            raise ValueError("cubic paths require 3*k+1 controls")
        if not isinstance(curveIndex, int) or not 0 <= curveIndex < (count-1)//3:
            raise IndexError("curve index out of range")
        coordinates, dimension = _coordinates(self.controlPoints)
        index = curveIndex * 3
        return coordinates[index:index+4], dimension


def _coordinates(points):
    """Validate sampling representations without changing evaluation APIs."""
    dimension = points[0].size if isinstance(points[0], Vector) else None
    if dimension is not None and dimension not in (2, 3):
        raise ValueError("sampling requires Vector2 or Vector3 controls")
    result = []
    for point in points:
        if dimension is None:
            if not isinstance(point, (int, float)):
                raise ValueError("sampling controls must share a representation")
            coords = (point,)
        else:
            if not isinstance(point, Vector) or point.size != dimension:
                raise ValueError("sampling controls must share a dimension")
            coords = tuple(point.vector)
        if not all(math.isfinite(value) for value in coords):
            raise ValueError("sampling controls must be finite")
        result.append(coords)
    return result, dimension


def _point(coords, dimension):
    return coords[0] if dimension is None else Vector(dimension, list(coords))


def _split(controls, t=0.5):
    levels = [controls]
    while len(levels[-1]) > 1:
        level = levels[-1]
        levels.append([tuple((1-t)*a+t*b for a,b in zip(left,right))
                       for left,right in zip(level,level[1:])])
    return ([level[0] for level in levels],
            [level[-1] for level in reversed(levels)])


def _chord_distance(point, start, end):
    chord = tuple(b-a for a,b in zip(start,end))
    length = math.hypot(*chord)
    offset = tuple(p-a for p,a in zip(point,start))
    if length == 0:
        return math.hypot(*offset)
    unit = tuple(value/length for value in chord)
    projection = max(0, min(length, sum(a*b for a,b in zip(offset,unit))))
    return math.hypot(*(a-projection*b for a,b in zip(offset,unit)))


def _flatness(polygon):
    """Reuse a segment's chord without changing control-to-segment distances."""
    start, end = polygon[0], polygon[-1]
    chord = tuple(b-a for a,b in zip(start,end))
    length = math.hypot(*chord)
    if length == 0:
        return max(math.hypot(*(p-a for p,a in zip(point,start)))
                   for point in polygon[1:-1])
    unit = tuple(value/length for value in chord)
    distances = []
    for point in polygon[1:-1]:
        offset = tuple(p-a for p,a in zip(point,start))
        projection = max(0, min(length, sum(a*b for a,b in zip(offset,unit))))
        distances.append(math.hypot(*(a-projection*b for a,b in zip(offset,unit))))
    return max(distances)


def _subdivide(controls, tolerance):
    # Stack order emits left intervals first; bounded without Python recursion.
    stack = [(controls, 0)]
    result = [controls[0]]
    while stack:
        polygon, depth = stack.pop()
        flatness = _flatness(polygon)
        if flatness <= tolerance or depth == 16:
            result.append(polygon[-1])
        else:
            left, right = _split(polygon)
            stack.append((right, depth+1))
            stack.append((left, depth+1))
    return result
