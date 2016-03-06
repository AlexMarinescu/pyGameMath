import math
import six.moves as sm


# Cubic Bezier Cruve
def cubicBezierPoint(t, p0, p1, p2, p3):
    u = 1 - t
    tt = t * t
    uu = u * u
    uuu = uu * u
    ttt = tt * t

    p = uuu * p0
    p += 3 * uu * t * p1
    p += 3 * u * tt * p2
    p += ttt + p3

    return p

# Quadratic Bezier Curve
def quadraticBezierPoint(t, p0, p1, p2):
    u = 1 - t
    tt = t * t
    uu = u * u

    p = uu * p0
    p += 2 * u * t * p1
    p += tt * p2

    return p

class BezierPath(object):
    def __init__(self):
        self.segments_per_curve = 10
        self.minimum_sqr_distance = 0.01
        self.divison_threshold = -0.99

        self.controlPoints = []
        self.curveCount = 0

    def setControlPoints(self, newControlPoints):
        self.controlPoints = newControlPoints
        self.curveCount = (len(self.controlPoints) - 1) / 3

    def getControlPoints(self):
        return self.controlPoints

    def interpolate(self, segmentPoints, scale):

        segmentPointsLen = len(segmentPoints)

        if segmentPointsLen < 2:
            return

        for i in sm.range(segmentPointsLen):
            if i == 0:
                p1 = segmentPoints[i]
                p2 = segmentPoints[i + 1]

                tangent = p2 - p1
                q1 = p1 + scale * tangent

                self.controlPoints.append(p1)
                self.controlPoints.append(q1)
            elif i == segmentPointsLen - 1:
                p0 = segmentPoints[i - 1]
                p1 = segmentPoints[i]

                tangent = p1 - p0
                q0 = p1 - scale * tangent

                self.controlPoints.append(q0)
                self.controlPoints.append(p1)
            else:
                p0 = segmentPoints[i - 1]
                p1 = segmentPoints[i]
                p2 = segmentPoints[i + 1]
                tangent = (p2 - p0).normalize()
                q0 = p1 - scale * tangent * (p1 - p0).magnitude()
                q1 = p1 + scale * tangent * (p2 - p1).magnitude()

                self.controlPoints.append(q0)
                self.controlPoints.append(p1)
                self.controlPoints.append(q1)

        self.curveCount = (len(self.controlPoints) - 1) / 3

    def samplePoints(self, sourcePoints, minSqrDistance, maxSqrDistance, scale):
        sourcePointsLen = len(sourcePoints)

        if (sourcePointsLen < 2):
            return

        samplePoints = []

        samplePoints.append(sourcePoints[0])

        potentialSamplePoint = sourcePoints[1]

        for i in sm.range(2, sourcePointsLen, 1):
            # These variables need better names
            pointDiff = potentialSamplePoint - sourcePoints[i]
            pointDIff2 = samplePoints[len(samplePoints) - 2] - sourcePoints[i]

            if (pointDiff.sqrMagnitude() > minSqrDistance and
                    pointDIff2.sqrMagnitude() > maxSqrDistance):
                samplePoints.append(potentialSamplePoint)
            potentialSamplePoint = sourcePoints[i]

        p1 = samplePoints[len(samplePoints) - 1]
        p0 = samplePoints[len(samplePoints) - 2]
        tangent = (p0 - potentialSamplePoint).normalize()
        d2 = (potentialSamplePoint - p1).magnitude()
        d1 = (p1 - p0).magnitude()
        p1 = p1 + tangent * ((d1 - d2) / 2)

        samplePoints.append(p1)
        samplePoints.append(potentialSamplePoint)

        self.interpolate(samplePoints, scale)

    def calculateBezerPoint(self, curveIndex, t):
        nodeIndex = curveIndex * 3

        p0 = self.controlPoints[nodeIndex]
        p1 = self.controlPoints[nodeIndex + 1]
        p2 = self.controlPoints[nodeIndex + 2]
        p3 = self.controlPoints[nodeIndex + 3]

        return cubicBezierPoint(t, p0, p1, p2, p3)

    # Find the drawing points of the curve using recusive divison
    # Less points but same accuracy
    def getDrawingPoints(self):
        drawingPoints = []

        for i in sm.range(self.curveCount):
            bezierCurveDrawingPoints = self.findDrawingPoints(i)

            if i != 0:
                del bezierCurveDrawingPoints[0]

            drawingPoints.append(bezierCurveDrawingPoints)

        return drawingPoints

    def findDrawingPoints(self, curveIndex):
        pointList = []

        left = self.calculateBezerPoint(curveIndex, 0)
        right = self.calculateBezerPoint(curveIndex, 1)

        pointList.append(left)
        pointList.append(right)

        self.findDrawingPointsAdded(curveIndex, 0, 1, pointList, 1)

        return pointList

    def findDrawingPointsAdded(self, curveIndex, t0, t1, pointList, insertionIndex):
        left = self.calculateBezerPoint(curveIndex, t0)
        right = self.calculateBezerPoint(curveIndex, t1)

        if (left - right).sqrMagnitude() < self.minimum_sqr_distance:
            return 0

        tMid = (t0 - t1) / 2
        mid = self.calculateBezerPoint(curveIndex, tMid)

        leftDirection = (left - mid).normalize()
        rightDirection = (right - mid).normalize()

        if (leftDirection.dot(rightDirection) > self.divison_threshold or
                math.abs(tMid - 0.5) < 0.0001):
            pointsAddedCount = 0

            pointsAddedCount += self.findDrawingPointsAdded(curveIndex, t0, tMid, pointList, insertionIndex)
            pointList.append(insertionIndex + pointsAddedCount, mid)
            pointsAddedCount += 1
            pointsAddedCount += self.findDrawingPointsAdded(curveIndex, tMid, t1, pointList, insertionIndex + pointsAddedCount)

            return pointsAddedCount

        return 0
