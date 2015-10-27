import six
from src import vector

def flip(plane):
    ''' Flips the plane.'''
    fA = -plane.a
    fB = -plane.b
    fC = -plane.c
    fD = -plane.D
    fNormal = -plane.normal
    return [fA, fB, fC, fD, fNormal]

def normalize(pdata):
    ''' Return the normalized plane.'''
    vec = vector.Vector(3, data=[pdata[0], pdata[1], pdata[2]])
    vecN = vec.normalize()

    length = vecN.magnitude()

    if length is not 0:
        return [vecN.vector[0], vecN.vector[1], vecN.vector[2], pdata[3] / length]
    else:
        print("Plane fail to normalize due to zero division.")
        return [0.0, 0.0, 0.0, 0.0]

class Plane(object):
    def __init__(self, *args):
        ''' Define a plane using 3 vectors, return a 4 item 1D array. '''
        if isinstance(args[0], list):
            self.data = args[0]
        elif isinstance(args[0], vector.Vector):
            # Define using 3 vectors
            v1 = args[1] - args[0]
            v2 = args[2] - args[0]
            n = v1.cross(v2)
            n.i_normalize()
            self.data = [n.vector[0], n.vector[1], n.vector[2], (n * -1.0).dot(args[0])]
        else:
            self.data = [0.0, 0.0, 0.0, 0.0]

        self.normal = vector.Vector(3, data=[0.0, 0.0, 0.0])

    def i_flip(self):
        ''' Flip the plane in its place. '''
        data = flip(self)
        self.a = data[0]
        self.b = data[1]
        self.c = data[2]
        self.d = data[3]
        self.normal = data[4]
        return self

    def flip(self):
        ''' Return a flipped plane. '''
        return Plane(flip(self))

    def dot(self, vec):
        ''' Return the dot product between a plane and 4D vector. '''
        return self.data[0] * vec.vector[0] + self.data[1] * vec.vector[1] + self.data[2] * vec.vector[2] + self.data[3] * vec.vector[3]

    def i_normalize(self):
        ''' Normalize the vector in place. '''
        self.data = normalize(self.data)
        return self

    def normalize(self):
        ''' Return the normalized plane.'''
        return Plane(normalize(self.data))

    def bestFitNormal(self, vecList):
        ''' Pass in a list of vectors to find the best fit normal. '''
        output = vector.Vector(3).zero()
        for i in six.range(len(vecList)):
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
        ''' Returns the location of the point. '''
        # If s > 0 then the point is on the same side as the normal. (front)
        # If s < 0 then the point is on the opposide side of the normal. (back)
        # If s = 0 then the point lies on the plane.
        s = plane.data[0] * point[0] + plane.data[1] * point[1] + plane.data[2] * point[2] + plane.data[3]

        if s > 0:
            return 1
        elif s < 0:
            return -1
        elif s == 0:
            return 0
        else:
            print("Not a clue where the point is.")

