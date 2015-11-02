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
    ''' Return the normalized plane.'''
    vec = vector.Vector(3, data=pdata)
    vecN = vec.normalize()

    length = vecN.magnitude()

    if length is not 0:
        return vecN.vector[0], vecN.vector[1], vecN.vector[2], pdata[3] / length
    else:
        print("Plane fail to normalize due to zero division.")
        return 0.0, 0.0, 0.0, 0.0

class Plane(object):
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
        ''' Create the plane from A,B,C,D. '''
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.normal = vector.cross(b - a, c - a).normalize()

    def fromPoints(self, a, b, c):
        '''Calculate the plane from A,B,C.'''
        self.a = a
        self.b = b
        self.c = c
        self.normal = vector.cross(b - a, c - a).normalize()
        self.d = self.normal.dot(self.a)

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
        ''' Normalize the vector in place. '''
        pdata = [self.a, self.b, self.c, self.d]
        self.a, self.b, self.c, self.d = normalize(pdata)
        return self

    def normalize(self):
        ''' Return the normalized plane.'''
        nPlane = Plane().clone()
        pdata = [self.a, self.b, self.c, self.d]
        nPlane.a, nPlane.b, nPlane.c, nPlane.d = normalize(pdata)
        return nPlane

    def bestFitNormal(self, vecList):
        ''' Pass in a list of vectors to find the best fit normal. '''
        output = vector.Vector(3).zero()
        for i in sm.range(len(vecList)):
            output.vector[0] += (vecList[i].vector[2] + vecList[i + 1].vector[2]) * (vecList[i].vector[1] - vecList[i + 1].vector[1])
            output.vector[1] += (vecList[i].vector[0] + vecList[i + 1].vector[0]) * (vecList[i].vector[2] - vecList[i + 1].vector[2])
            output.vector[2] += (vecList[i].vector[1] + vecList[i + 1].vector[1]) * (vecList[i].vector[0] - vecList[i + 1].vector[0])
        return output.normalize()

    def bestFitD(self, vecList, bestFitNormal):
        ''' Returns the best fit D from a list of vectors using the best fit normal. '''
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

