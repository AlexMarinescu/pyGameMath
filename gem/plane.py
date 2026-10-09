import six.moves as sm
from gem import vector

def flip(plane):
    ''' Flips the plane.'''
    fA = -plane[0]
    fB = -plane[1]
    fC = -plane[2]
    fD = -plane[3]
    fNormal = -plane[4]
    return [fA, fB, fC, fD, fNormal]

def normalize(pdata):
    ''' Divide all four coefficients by the magnitude of (a, b, c). '''
    vec = vector.Vector(3, data=pdata)
    length = vec.magnitude()
    return pdata[0] / length, pdata[1] / length, pdata[2] / length, pdata[3] / length

class Plane(object):
    ''' Plane a*x + b*y + c*z + d = 0, with normal = Vector([a, b, c]). '''
    def __init__(self):
        ''' Plane class constructor. '''
        self.normal = vector.Vector(3, data=[0.0, 0.0, 0.0])
        self.a = 0
        self.b = 0
        self.c = 0
        self.d = 0

    def clone(self):
        '''Create a new Plane with similar propertise.'''
        nPlane = Plane()
        nPlane.normal = self.normal.clone()
        nPlane.a = self.a
        nPlane.b = self.b
        nPlane.c = self.c
        nPlane.d = self.d
        return nPlane

    def fromCoeffs(self, a, b, c, d):
        ''' Set scalar coefficients, preserving their scale and orientation. '''
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.normal = vector.Vector(3, data=[a, b, c])

    def fromPoints(self, a, b, c):
        ''' Set a unit-normal plane through three noncollinear 3D Vectors. '''
        normal = vector.cross(b - a, c - a)
        vector._require_nonzero(normal.size, normal.vector)
        self.normal = normal.normalize()
        self.a, self.b, self.c = self.normal.vector
        self.d = -self.normal.dot(a)

    def i_flip(self):
        ''' Flip the plane in its place. '''
        data = flip([self.a, self.b, self.c, self.d, self.normal])
        self.a = data[0]
        self.b = data[1]
        self.c = data[2]
        self.d = data[3]
        self.normal = data[4]
        return self

    def flip(self):
        ''' Return a flipped plane. '''
        nPlane = Plane()
        data = flip([self.a, self.b, self.c, self.d, self.normal])
        nPlane.a = data[0]
        nPlane.b = data[1]
        nPlane.c = data[2]
        nPlane.d = data[3]
        nPlane.normal = data[4]
        return nPlane

    def dot(self, vec):
        ''' Return the dot product between a plane and 4D vector. '''
        return self.a * vec.vector[0] + self.b * vec.vector[1] + self.c * vec.vector[2] + self.d * vec.vector[3]

    def i_normalize(self):
        ''' Normalize all coefficients and synchronize the normal in place. '''
        pdata = [self.a, self.b, self.c, self.d]
        self.a, self.b, self.c, self.d = normalize(pdata)
        self.normal = vector.Vector(3, data=[self.a, self.b, self.c])
        return self

    def normalize(self):
        ''' Return the normalized plane.'''
        nPlane = Plane()
        pdata = [self.a, self.b, self.c, self.d]
        nPlane.a, nPlane.b, nPlane.c, nPlane.d = normalize(pdata)
        nPlane.normal = vector.Vector(3, data=[nPlane.a, nPlane.b, nPlane.c])
        return nPlane

    def bestFitNormal(self, vecList):
        ''' Return a unit Newell normal for an ordered polygon of 3D Vectors.

        The last vertex wraps to the first; a repeated first vertex is allowed.
        Reversing vertex order reverses the normal.
        '''
        output = vector.Vector(3).zero()
        if len(vecList):
            origin = vecList[0].vector
            current = (0.0, 0.0, 0.0)
        for i in sm.range(len(vecList)):
            point = vecList[(i + 1) % len(vecList)].vector
            # Newell's sum is translation invariant. Use local coordinates so
            # large world offsets cannot obscure the polygon's edge geometry.
            following = (point[0] - origin[0], point[1] - origin[1], point[2] - origin[2])
            output.vector[0] += (current[2] + following[2]) * (current[1] - following[1])
            output.vector[1] += (current[0] + following[0]) * (current[2] - following[2])
            output.vector[2] += (current[1] + following[1]) * (current[0] - following[0])
            current = following
        vector._require_nonzero(output.size, output.vector)
        return output.normalize()

    def bestFitD(self, vecList, bestFitNormal):
        ''' Return signed D = average(normal.dot(point)); use d = -D.

        For a unit normal, D is the geometric offset in normal.dot(point) = D.
        It is not forced nonnegative. Each supplied vertex contributes once.
        '''
        val = 0.0
        for vec in vecList:
            val += vec.dot(bestFitNormal)
        return val / len(vecList)

    def point_location(self, plane, point):
        ''' Returns the location of the point. Point is a tuple. '''
        # If s > 0 then the point is on the same side as the normal. (front)
        # If s < 0 then the point is on the opposide side of the normal. (back)
        # If s = 0 then the point lies on the plane.
        s = plane.a * point[0] + plane.b * point[1] + plane.c * point[2] + plane.d

        if s > 0:
            return 1
        elif s < 0:
            return -1
        elif s == 0:
            return 0
        else:
            print("Not a clue where the point is.")

