import math
import sys

# Vector N Dimensions Class
class Vector (object):

    def __init__(self, array):
        self.vec = array
        self.len = len(array)

    # Clone the vector.
    def duplicate(self):
        return Vector(self.vec)

    # Overloading +=
    def __iadd__(self, vector):
        for x in xrange(self.len):
            self.vec[x] += vector.vec[x]
        return Vector(self.vec)
    
    # Overloading -=
    def __isub__(self, vector):
        for x in xrange(self.len):
            self.vec[x] -= vector.vec[x]
        return Vector(self.vec)
    
    # Overloading *=
    def __imul__(self, scalar):
        for x in xrange(self.len):
            self.vec[x] *= scalar
        return Vector(self.vec)
        
    # Overloading *
    def __mul__(self, scalar):
        for x in xrange(self.len):
            self.vec[x] *= scalar
        return Vector(self.vec)
    
    # Overloading /=
    def __itruediv__(self, scalar):
        for x in xrange(self.len):
            self.vec[x] /= scalar
        return Vector(self.vec)
    
    # Overloading +
    def __add__(self, vectorB):
        vectorA = clone()
        for x in xrange(self.len):
            vectorA.vec[x] += vectorB.vec[x]
        return vectorA

    # Overloading -
    def __sub__(self, vectorB):
        vectorA = clone()
        for x in xrange(self.len):
            vectorA.vec[x] -= vectorB.vec[x]
        return vectorA

    # Overloading -(a)
    def __neg__(self):
        vectorA = clone()
        for x in xrange(self.len):
            vectorA.vec[x] = -vectorA.vec[x]
        return vectorA

    # Overloading !=
    def __ne__(self, vector):
        true_or_false = False
        for x in xrange(self.len):
            if (self.vec[x] == vector.vec[x]):
                true_or_false = True
            else:
                true_of_false = False
        return true_or_false

    # Overloading ==
    def __eq__(self, vector):
        true_or_false = False
        for x in xrange(self.len):
            if (self.vec[x] == vector.vec[x]):
                true_or_false = True
            else:
                true_of_false = False
        return true_or_false

    def magnitude(self):
        result = 0
        for x in xrange(self.len):
            result += self.vec[x] * self.vec[x]
        return math.sqrt(result)

    def scale(self, value):
        for x in xrange(self.len):
            self.vec[x] *= value
            
    def zero(self):
        for x in xrange(self.len):
            self.vec[x] = 0.0
            
    def invert(self):
        for x in xrange(self.len):
            self.vec[x] = -self.vec[x]

    def normalize(self):
        length = self.magnitude()

        if(length != 0):
            for x in xrange(self.len):
                self.vec[x] /= length

    def output(self):
        print "Vector ", self.len, "D: ", self.vec

# Dot Product (Multivectors not implemented).
def dot(vectorA, vectorB):
    if vectorA.len == vectorB.len:
        result = 0
        for x in xrange(vectorA.len):
            result += vectorA.vec[x] * vectorB.vec[x]
        return result
    else:
        print "Blah multivectors... Don't feel like it."

# Cross Product (only 3D and 7D)
# Though the 7D one won't be implemented yet, maybe later.
def cross(vectorA, vectorB):
    vectorC = Vector([0.0,0.0,0.0])
    if vectorA.len == 3 and vectorB.len == 3:
        vectorC.vec[0] = vectorA.vec[1] * vectorB.vec[2] - vectorA.vec[2] * vectorB.vec[1]
        vectorC.vec[1] = vectorA.vec[2] * vectorB.vec[0] - vectorA.vec[0] * vectorB.vec[2]
        vectorC.vec[2] = vectorA.vec[0] * vectorB.vec[1] - vectorA.vec[1] * vectorB.vec[0]
        return vectorC
    else:
        print "Only 3D vectors."

# The formulas for the following two functions are from the GLSL specification.
# Reflection
def reflect(incidentVector, normal):
    return incidentVector - (2.0 * dot(incidentVector, normal) * normal)
    
# Refraction
def refract(indexofrefraction, incidentVector, normal):
    # This is gonna be used in a few places and its quite long so just store the result once.
    dotNI = dot(normal, incidentVector)
    k = 1.0 - indexofrefraction * indexofrefraction * indexofrefraction * (1.0 - dotNI * dotNI)
    if k < 0.0:
        container = []
        # The length can be the normal or incidentVectors because they will be the same length.
        for x in xrange(normal.len):
            container.append(0.0)
        return Vector(container)
    else:
        return indexofrefraction * incidentVector - (indexofrefraction * dotNI + math.sqrt(k)) * normal
    
# Vector Stack Class
class VectorStack(object):

    def __init__(self, vector_list):
        self.stack = vector_list

    def __getitem__(self, key):
        return self.stack[key]

    def __setitem__(self, key, value):
        self.stack[key] = value

    def push(self, vector):
        self.stack.append(vector)

    def assign(self, location, vector):
        self.stack[location] = vector

    def pop(self,*args):
        self.stack.pop(args[0])

    def output(self):
        print self.stack
