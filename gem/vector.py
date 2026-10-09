import math
import six.moves as sm

REFRENCE_VECTOR_2 = [0.0 for x in sm.range(2)]
REFRENCE_VECTOR_3 = [0.0 for x in sm.range(3)]
REFRENCE_VECTOR_4 = [0.0 for x in sm.range(4)]

def zero_vector(size):
    ''' Return a zero filled vector list of the requested size '''
    # Return a copy of the reference vector, this is faster than making a new one
    if size == 2:
        return REFRENCE_VECTOR_2[:]
    if size == 3:
        return REFRENCE_VECTOR_3[:]
    if size == 4:
        return REFRENCE_VECTOR_4[:]
    else:
        return [0.0 for _ in sm.range(size)]

IREFRENCE_VECTOR_2 = [1.0 for x in sm.range(2)]
IREFRENCE_VECTOR_3 = [1.0 for x in sm.range(3)]
IREFRENCE_VECTOR_4 = [1.0 for x in sm.range(4)]

def one_vector(size):
    ''' Return a one filled vector list of the requested size '''
    # Return a copy of the reference vector, this is faster than making a new one
    if size == 2:
        return IREFRENCE_VECTOR_2[:]
    if size == 3:
        return IREFRENCE_VECTOR_3[:]
    if size == 4:
        return IREFRENCE_VECTOR_4[:]
    else:
        return [1.0 for _ in sm.range(size)]

# Vector Functions
def lerp(vecA, vecB, time):
    '''Linear interpolation between two vectors.'''
    #Note: The commented out part is the precise version of linear interpolation
    #return (vecA * time) + (vecB * (1.0 - time))
    return vecA + ((vecB - vecA) * time)

def cross(vecA, vecB):
    ''' Cross product between two 3D vectors, returns a vector.'''
    vecC = Vector(3)
    vecC.vector[0] = vecA.vector[1] * vecB.vector[2] - vecA.vector[2] * vecB.vector[1]
    vecC.vector[1] = vecA.vector[2] * vecB.vector[0] - vecA.vector[0] * vecB.vector[2]
    vecC.vector[2] = vecA.vector[0] * vecB.vector[1] - vecA.vector[1] * vecB.vector[0]
    return vecC

def reflect(incidentVec, normal):
    '''Reflect a vector'''
    return incidentVec - (normal * (2.0 * incidentVec.dot(normal)))

def refract(IOR, incidentVec, normal):
    ''' Return a refracted direction for IOR = n1 / n2.

    incidentVec and normal must be unit vectors of matching dimension.
    The incident direction points toward the interface; normal points toward
    the incident medium, so normal.dot(incidentVec) <= 0. Inputs are not
    normalized or flipped. Total internal reflection returns a zero Vector
    of the same dimension. Neither input is mutated.
    '''
    dotNI = normal.dot(incidentVec)
    k = 1.0 - IOR * IOR * (1.0 - dotNI * dotNI)

    if k < 0.0:
        return Vector(normal.size)
    else:
        scalar = IOR * dotNI + math.sqrt(k)
        return (incidentVec * IOR) - (normal * scalar)

# 2D - get angle of the vector
def toAngle(vector):
    return math.atan2(vector[1], vector[0])

# 2D -  90 degree rotation
def lperp(vector):
    return Vector(2, data = [-vector[1], vector[0]])

# 2D -  -90 degree rotation
def rperp(vector):
    return Vector(2, data = [vector[1], -vector[0]])

def vec_add(size, vecA, vecB):
    return [(vecA[i] + vecB[i]) for i in sm.range(size)]

def s_vec_add(size, vecA, scalar):
    return [(vecA[i] + scalar) for i in sm.range(size)]

def vec_sub(size, vecA, vecB):
    return [(vecA[i] - vecB[i]) for i in sm.range(size)]

def s_vec_sub(size, vecA, scalar):
    return [(vecA[i] - scalar) for i in sm.range(size)]

def vec_mul(size, vecA, scalar):
    return [(vecA[i] * scalar) for i in sm.range(size)]

def vec_div(size, vecA, scalar):
    return [(vecA[i] / scalar) for i in sm.range(size)]

def vec_neg(size, vecA):
    return [(-vecA[i]) for i in sm.range(size)]

def dot(size, vecA, vecB):
    dp = 0
    for i in sm.range(size):
        dp +=  vecA[i] * vecB[i]
    return dp

def magnitude(size, vecA):
    # Retain legacy arithmetic for nonfinite inputs; no new NaN/Inf policy.
    if any(math.isnan(vecA[i]) or math.isinf(vecA[i]) for i in sm.range(size)):
        return math.sqrt(sum(vecA[i] * vecA[i] for i in sm.range(size)))
    scale = max([abs(vecA[i]) for i in sm.range(size)] or [0.0])
    if scale == 0.0:
        return 0.0
    length = 0.0
    for i in sm.range(size):
        length = math.hypot(length, vecA[i] / scale)
    return length * scale

def normalize(size, vecA):
    """Normalize finite components stably; an exact zero returns zeros."""
    if any(math.isnan(vecA[i]) or math.isinf(vecA[i]) for i in sm.range(size)):
        length = magnitude(size, vecA)
        if length != 0.0:
            return [vecA[i] / length for i in sm.range(size)]
        return zero_vector(size)
    scale = max([abs(vecA[i]) for i in sm.range(size)] or [0.0])
    if scale != 0.0:
        scaled = [vecA[i] / scale for i in sm.range(size)]
        # Finite scaled components have maximum absolute value exactly 1.
        # Keep magnitude's chained hypot order without rechecking/rescaling.
        length = 0.0
        for value in scaled:
            length = math.hypot(length, value)
        return [value / length for value in scaled]
    return zero_vector(size)

def _require_nonzero(size, values):
    """Retain geometric callers' exact-zero normalization error."""
    if size != 0 and all(values[i] == 0.0 for i in sm.range(size)):
        raise ZeroDivisionError("Cannot normalize a zero geometric direction")

def maxV(size, vecA, vecB):
    return [vecA[i] if vecA[i] > vecB[i] else vecB[i] for i in sm.range(size)]

def minV(size, vecA, vecB):
    return [vecA[i] if vecA[i] < vecB[i] else vecB[i] for i in sm.range(size)]

def maxS(size, vecA):
    mScalar = vecA[0]
    for i in sm.range(size):
        if vecA[i] > mScalar:
            mScalar = vecA[i]
    return mScalar

def minS(size, vecA):
    mScalar = vecA[0]
    for i in sm.range(size):
        if vecA[i] < mScalar:
            mScalar = vecA[i]
    return mScalar

def clamp(size, value, minS, maxS):
    """Return independent clamped storage, preserving value and bound lists."""
    output = value[:]
    for i in sm.range(size):
        # Check to see if greater than max
        output[i] = maxS[i] if output[i] > maxS[i] else output[i]
        #output[i] = (output[i] > max[i]) ? max[i] : output[i]
        # Check to see is less than min
        output[i] = minS[i] if output[i] < minS[i] else output[i]
        #output[i] = (output[i] < min[i]) ? min[i] : output[i]
    return Vector(size, data=output)

def transform(size, position, matrix):
    ''' Transform a position list using a nested row-vector matrix.

    A size-by-size matrix applies the ordinary row-vector product, including
    any explicit homogeneous component. A (size+1)-by-(size+1) matrix treats
    position as a point with implicit w=1 and returns size components.
    No perspective division is performed. Neither input is mutated.
    Implicit promotion is intended for affine positions; projective results
    require projection/unprojection handling to obtain Cartesian coordinates.
    '''
    output = zero_vector(size)
    for i in sm.range(size):
        for j in sm.range(size):
            output[i] += position[j] * matrix[j][i]
        if len(matrix) == size + 1:
            output[i] += matrix[size][i]

    return output

class Vector(object):
    def __init__(self, size, data=None):
        self.size = size

        if data is None:
            self.vector = zero_vector(self.size)
        else:
            self.vector = data

    def __repr__(self):
        return 'Vector: size:{} , data:{}'.format(self.size, self.vector)

    def __add__(self, other):
        if isinstance(other, Vector):
            vecList = vec_add(self.size, self.vector, other.vector)
            return Vector(self.size, data=vecList)
        elif isinstance(other, int) or isinstance(other, float):
            vecList = s_vec_add(self.size, self.vector, other)
            return Vector(self.size, data=vecList)
        else:
            return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, Vector):
            self.vector = vec_add(self.size, self.vector, other.vector)
            return self
        elif isinstance(other, int) or isinstance(other, float):
            self.vector = s_vec_add(self.size, self.vector, other)
            return self
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Vector):
            vecList = vec_sub(self.size, self.vector, other.vector)
            return Vector(self.size, data=vecList)
        elif isinstance(other, int) or isinstance(other, float):
            vecList = s_vec_sub(self.size, self.vector, other)
            return Vector(self.size, data=vecList)
        else:
            return NotImplemented

    def __isub__(self, other):
        if isinstance(other, Vector):
            self.vector = vec_sub(self.size, self.vector, other.vector)
            return self
        elif isinstance(other, int) or isinstance(other, float):
            self.vector = s_vec_sub(self.size, self.vector, other)
            return self
        else:
            return NotImplemented

    def __mul__(self, scalar):
        if isinstance(scalar, int) or isinstance(scalar, float):
            vecList = vec_mul(self.size, self.vector, scalar)
            return Vector(self.size, data=vecList)
        else:
            return NotImplemented

    def __imul__(self, scalar):
        if isinstance(scalar, int) or isinstance(scalar, float):
            self.vector = vec_mul(self.size, self.vector, scalar)
            return self
        else:
            return NotImplemented

    def __div__(self, scalar):
        if isinstance(scalar, int) or isinstance(scalar, float):
            vecList = vec_div(self.size, self.vector, scalar)
            return Vector(self.size, data=vecList)
        else:
            return NotImplemented

    def __truediv__(self, scalar):
        if isinstance(scalar, int) or isinstance(scalar, float):
            vecList = vec_div(self.size, self.vector, scalar)
            return Vector(self.size, data=vecList)
        else:
            return NotImplemented

    def __idiv__(self, scalar):
        if isinstance(scalar, int) or isinstance(scalar, float):
            self.vector = vec_div(self.size, self.vector, scalar)
            return self
        else:
            return NotImplemented

    def __itruediv__(self, scalar):
        if isinstance(scalar, int) or isinstance(scalar, float):
            self.vector = vec_div(self.size, self.vector, scalar)
            return self
        else:
            return NotImplemented

    def __eq__(self, vecB):
        if not isinstance(vecB, Vector):
            return NotImplemented
        if self.size != vecB.size:
            return False
        for i in range(self.size):
            if self.vector[i] != vecB.vector[i]:
                return False
        return True

    def __ne__(self, vecB):
        if not isinstance(vecB, Vector):
            return NotImplemented
        if self.size != vecB.size:
            return True
        for i in range(self.size):
            if self.vector[i] != vecB.vector[i]:
                return True
        return False

    def __neg__(self):
        vecList = vec_neg(self.size, self.vector)
        return Vector(self.size, data=vecList)

    def clone(self):
        return Vector(self.size, data=self.vector[:])

    def one(self):
        self.vector = one_vector(self.size)
        return self

    def zero(self):
        self.vector = zero_vector(self.size)
        return self

    def negate(self):
        vecList = vec_neg(self.size, self.vector)
        return Vector(self.size, data=vecList)

    def maxV(self, vecB):
        return Vector(self.size, data=maxV(self.size, self.vector, vecB.vector))

    def maxS(self):
        return maxS(self.size, self.vector)

    def minV(self, vecB):
        return Vector(self.size, data=minV(self.size, self.vector, vecB.vector))

    def minS(self):
        return minS(self.size, self.vector)

    def magnitude(self):
        return magnitude(self.size, self.vector)

    def clamp(self, size, value, minS, maxS):
        ''' Returns a new clamped vector. '''
        return clamp(size, value, minS, maxS)

    def i_clamp(self, size, value, minS, maxS):
        ''' Replace receiver storage with clamped values; preserve input lists. '''
        self.vector = clamp(size, value, minS, maxS).vector
        return self

    def i_normalize(self):
        ''' Normalize the vector in place. '''
        self.vector = normalize(self.size, self.vector)
        return self

    def normalize(self):
        ''' Return a new normalized vector. '''
        vecList = normalize(self.size, self.vector)
        return Vector(self.size, data=vecList)

    def dot(self, vecB):
        ''' Return the dot product between two vectors. '''
        if isinstance(vecB, Vector):
            return dot(self.size, self.vector, vecB.vector)
        else:
            return NotImplemented

    def isInSameDirection(self, otherVec):
        ''' Return a boolean if the input vector if is in the same direction as the one it's compared against. '''
        if isinstance(otherVec, Vector):
            return self.dot(otherVec) > 0
        else:
            return NotImplemented

    def isInOppositeDirection(self, otherVec):
        ''' Return a boolean if the input vector if is in the opposite direction as the one it's compared against. '''
        if isinstance(otherVec, Vector):
            return self.dot(otherVec) < 0
        else:
            return NotImplemented

    def barycentric(self, a, b, c):
        ''' Compute barycentric coordinates (u, v, w) for point p with respect to triangle (a, b, c). '''
        # This method is shown in Christer Ericson's Real-Time Collision Detection book
        # Basically is cramer's rule to solve a linear system.
        # Returns [u, v, w] list

        v0 = b - a
        v1 = c - a
        v2 = self - a

        d00 = v0.dot(v0)
        d01 = v0.dot(v1)
        d11 = v1.dot(v1)
        d20 = v2.dot(v0)
        d21 = v2.dot(v1)

        denom = d00 * d11 - d01 * d01

        v = (d11 * d20 - d01 * d21) / denom
        w = (d00 * d21 - d01 * d20) / denom
        u = 1.0 - v - w

        return [u, v, w]

    def transform(self, position, matrix):
        ''' Return a new Vector from a position list and nested matrix.

        Uses this receiver's size; see transform for homogeneous rules.
        The receiver and inputs are preserved.
        '''
        vecList = transform(self.size, position, matrix)
        return Vector(self.size, data=vecList)

    def i_transform(self, position, matrix):
        ''' Replace this Vector with the transformed position list.

        Uses this receiver's size; see transform for homogeneous rules.
        Returns self without mutating the supplied list or nested matrix.
        '''
        self.vector = transform(self.size, position, matrix)
        return self

    # Return common components of the vector as a group
    # Vector Swizzling, similar to GLSL
    def xy(self):
        return Vector(2, [self.vector[0], self.vector[1]])

    def yz(self):
        return Vector(2, [self.vector[1], self.vector[2]])

    def xz(self):
        return Vector(2, [self.vector[0], self.vector[2]])

    def xw(self):
        return Vector(2, [self.vector[0], self.vector[3]])

    def yw(self):
        return Vector(2, [self.vector[1], self.vector[3]])

    def zw(self):
        return Vector(2, [self.vector[2], self.vector[3]])

    def xyw(self):
        return Vector(3, [self.vector[0], self.vector[1], self.vector[3]])

    def yzw(self):
        return Vector(3, [self.vector[1], self.vector[2], self.vector[3]])

    def xzw(self):
        return Vector(3, [self.vector[0], self.vector[2], self.vector[3]])

    def xyz(self):
        return Vector(3, [self.vector[0], self.vector[1], self.vector[2]])

    # 3D vector identities
    def right(self):
        return Vector(3, data=[ 1.0, 0.0, 0.0])

    def left(self):
        return Vector(3, data=[-1.0, 0.0, 0.0])

    def front(self):
        return Vector(3, data=[0.0, 0.0, -1.0])

    def back(self):
        return Vector(3, data=[ 0.0, 0.0, 1.0])

    def up(self):
        return Vector(3, data=[ 0.0, 1.0, 0.0])

    def down(self):
        return Vector(3, data=[0.0, -1.0, 0.0])
