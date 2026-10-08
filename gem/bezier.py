"""Quadratic and cubic Bezier evaluation for scalars and gem Vectors.

Parameters are not clamped. BezierPath evaluates cubic segments with a
shared endpoint layout of 3*k+1 control points; adaptive sampling is not
part of this supported module.
"""


def cubicBezierPoint(t, p0, p1, p2, p3):
    """Evaluate the cubic Bernstein polynomial without changing controls."""
    u = 1 - t
    return (p0 * (u * u * u) + p1 * (3 * u * u * t)
            + p2 * (3 * u * t * t) + p3 * (t * t * t))


def quadraticBezierPoint(t, p0, p1, p2):
    """Evaluate the quadratic Bernstein polynomial."""
    u = 1 - t
    return p0 * (u * u) + p1 * (2 * u * t) + p2 * (t * t)


class BezierPath(object):
    """Evaluate cubic segments; control-point storage remains caller-owned.

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
