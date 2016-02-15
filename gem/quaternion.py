import math
import six.moves as sm
from gem import vector
from gem import matrix

def quat_identity():
    ''' Returns the quaternion identity. '''
    return [1.0, 0.0, 0.0, 0.0]

def quat_add(quat, quat1):
    ''' Add two quaternions. '''
    return [quat[0] + quat1[0], quat[1] + quat1[1], quat[2] + quat1[2], quat[3] + quat1[3]]

def quat_sub(quat, quat1):
    ''' Subtract two quaternions. '''
    return [quat[0] - quat1[0], quat[1] - quat1[1], quat[2] - quat1[2], quat[3] - quat1[3]]

def quat_mul_quat(quat, quat1):
    ''' Multiply a quaternion with a quaternion. '''
    w = quat[0] * quat1[0] - quat[1] * quat1[1] - quat[2] * quat1[2] - quat[3] * quat1[3]
    x = quat[0] * quat1[1] + quat[1] * quat1[0] + quat[2] * quat1[3] - quat[3] * quat1[2]
    y = quat[0] * quat1[2] + quat[2] * quat1[0] + quat[3] * quat1[1] - quat[1] * quat1[3]
    z = quat[0] * quat1[3] + quat[3] * quat1[0] + quat[1] * quat1[2] - quat[2] * quat1[1]
    return [w, x, y, z]

def quat_mul_vect(quat, vect):
    ''' Multiply a quaternion with a vector. '''
    w = -quat[1] * vect[0] - quat[2] * vect[1] - quat[3] * vect[2]
    x =  quat[0] * vect[0] + quat[2] * vect[2] - quat[3] * vect[1]
    y =  quat[0] * vect[1] + quat[3] * vect[0] - quat[1] * vect[2]
    z =  quat[0] * vect[2] + quat[1] * vect[1] - quat[2] * vect[0]
    return [w, x, y, z]

def quat_mul_float(quat, scalar):
    ''' Multiply a quaternion with a scalar (float). '''
    return [quat[0] * scalar, quat[1] * scalar, quat[2] * scalar, quat[3] * scalar]

def quat_div_float(quat, scalar):
    ''' Divide a quaternion with a scalar (float). '''
    return [quat[0] / scalar, quat[1] / scalar, quat[2] / scalar, quat[3] / scalar]

def quat_neg(quat):
    ''' Negate the elements of a quaternion. '''
    return [-quat[0], -quat[1], -quat[2], -quat[3]]

def quat_dot(quat1, quat2):
    ''' Dot product between two quaternions. Returns a scalar. '''
    rdp= 0
    for i in sm.range(4):
        rdp += quat1[i] * quat2[i]
    return rdp

def quat_magnitude(quat):
    ''' Compute magnitude of a quaternion. Returns a scalar. '''
    rmg = 0
    for i in sm.range(4):
        rmg += quat[i] * quat[i]
    return math.sqrt(rmg)

def quat_normalize(quat):
    ''' Returns a normalized quaternion. '''
    length = quat_magnitude(quat)
    oquat = quat_identity()
    if length is not 0:
        for i in sm.range(4):
            oquat[i] = quat[i] / length
    return oquat

def quat_conjugate(quat):
    ''' Returns the conjugate of a quaternion. '''
    idquat = quat_identity()
    for i in sm.range(4):
        idquat[i] = -quat[i]
    idquat[0] = -idquat[0]
    return idquat

def quat_inverse(quat):
    ''' Returns the inverse of a quaternion. '''
    lengthSquared = quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]

    return [quat[0] / lengthSquared,
            quat[1] / lengthSquared,
            quat[2] / lengthSquared,
            quat[3] / lengthSquared]

def quat_from_axis_angle(axis, theta):
    ''' Returns a quaternion from a given axis and a angle. '''
    thetaOver2 = theta * 0.5
    sto2 = math.sin(math.radians(thetaOver2))
    cto2 = math.cos(math.radians(thetaOver2))

    quat1List = []
    if isinstance(axis, vector.Vector):
        axis.i_normalize()
        quat1List = [cto2, axis.vector[0] * sto2, axis.vector[1] * sto2, axis.vector[2] * sto2]
    elif isinstance(axis, list):
        naxis = axis.normalize()
        quat1List = (cto2, naxis[0] * sto2, naxis[1] * sto2, naxis[2] * sto2)
    else:
        return NotImplemented

    return Quaternion(data=quat1List)

def quat_rotate(origin, axis, theta):
    ''' Returns a vector that is rotated around an axis. '''
    thetaOver2 = theta * 0.5
    sinThetaOver2 = math.sin(math.radians(thetaOver2))
    cosThetaOver2 = math.cos(math.radians(thetaOver2))
    quat = Quaternion(data = [cosThetaOver2, axis[0] * sinThetaOver2, axis[1] * sinThetaOver2, axis[2] * sinThetaOver2])
    rotation = (quat * origin) * quat.conjugate()
    return vector.Vector(3, data=[rotation.data[1], rotation.data[2], rotation.data[3]])

def quat_rotate_x_from_angle(theta):
    ''' Creates a quaternion that rotates around X axis given an angle. '''
    thetaOver2 = theta * 0.5
    cto2 = math.cos(thetaOver2)
    sto2 = math.sin(thetaOver2)
    return [cto2, sto2, 0.0, 0.0]

def quat_rotate_y_from_angle(theta):
    ''' Creates a quaternion that rotates around Y axis given an angle. '''
    thetaOver2 = theta * 0.5
    cto2 = math.cos(thetaOver2)
    sto2 = math.sin(thetaOver2)
    return [cto2, 0.0, sto2, 0.0]

def quat_rotate_z_from_angle(theta):
    ''' Creates a quaternion that rotates around Z axis given an angle. '''
    thetaOver2 = theta * 0.5
    cto2 = math.cos(thetaOver2)
    sto2 = math.sin(thetaOver2)
    return [cto2, 0.0, 0.0, sto2]

def quat_rotate_from_axis_angle(axis, theta):
    ''' Creates a quaternion that rotates around an arbitary axis given an angle. '''
    thetaOver2 = theta * 0.5
    sto2 = math.sin(math.radians(thetaOver2))
    cto2 = math.cos(math.radians(thetaOver2))

    quat1List = []
    if isinstance(axis, vector.Vector):
        axis.i_normalize()
        quat1List = [cto2, axis.vector[0] * sto2, axis.vector[1] * sto2, axis.vector[2] * sto2]
    elif isinstance(axis, list):
        naxis = axis.normalize()
        quat1List = (cto2, naxis[0] * sto2, naxis[1] * sto2, naxis[2] * sto2)
    else:
        return NotImplemented

    quat1 = Quaternion(data=quat1List)
    rotation = (quat1 * axis) * quat1.conjugate()
    return rotation

def quat_rotate_vector(quat, vec):
    ''' Rotates a vector by a quaternion, returns a vector. '''
    outQuat = (quat * vec) * quat.conjugate()
    return vector.Vector(3, data=[outQuat.data[1], outQuat.data[2], outQuat.data[3]])

def quat_pow(quat, exp):
    ''' Returns a quaternion to the power of N. '''
    quatExp = Quaternion()

    if quat.data[0] is not 0.0:
        angle = math.acos(quat.data[0])
        newAngle = angle * exp
        quatExp.data[0] = math.cos(newAngle)
        divAngle = math.sin(newAngle) / math.sin(angle)
        quatExp.data[1] *= divAngle
        quatExp.data[2] *= divAngle
        quatExp.data[3] *= divAngle
    return quatExp

def quat_log(quat):
    ''' Returns the logatithm of a quaternion. '''
    alpha = math.acos(quat.data[0])
    sinAlpha = math.sin(alpha)

    outList = [1.0, 0.0, 0.0, 0.0]

    if sinAlpha > 0.0:
        outList[1] = quat.data[1] * alpha / sinAlpha
        outList[2] = quat.data[2] * alpha / sinAlpha
        outList[3] = quat.data[3] * alpha / sinAlpha
    else:
        outList = quat.data

    return outList

def quat_lerp(quat0, quat1, t):
    ''' Linear interpolation between two quaternions. '''
    k0 = 1.0 - t
    k1 = t

    output = Quaternion()
    output = (quat0 * k0) + (quat1 * k1)

    return output

def quat_slerp(quat0, quat1, t):
    ''' Spherical interpolation between two quaternions. '''
    k0 = 0.0
    k1 = 0.0

    output = Quaternion()
    quat1Neg = Quaternion()
    cosTheta = quat0.dot(quat1)

    if cosTheta < 0.0:
        quat1Neg = quat1.negate()
        cosTheta = -cosTheta
    else:
        quat1Neg = quat1

    if cosTheta > 0.999:
        k0 = 1.0 - t
        k1 = t
    else:
        theta = math.acos(cosTheta)
        oneOverSinTheta = 1.0 / math.sin(theta)
        k0 = math.sin((1.0 - t) * theta) * oneOverSinTheta
        k1 = math.sin(t * theta) * oneOverSinTheta

    output = (quat0 * k0) + (quat1Neg * k1)

    return output

def quat_slerp_no_invert(quat0, quat1, t):
    ''' Spherical interpolation between two quaternions, it does not check for theta > 90. Used by SQUAD. '''
    dotP = quat0.dot(quat1)

    output = Quaternion()

    if (dotP > -0.95) and (dotP < 0.95):
        angle = math.acos(dotP)
        k0 = math.sin(angle * (1.0 - t)) / math.sin(angle)
        k1 = math.sin(t * angle) / math.sin(angle)

        output = (quat0 * k0) + (quat1 * k1)
    else:
        output = quat_lerp(quat0, quat1, t)

    return output

def quat_squad(quat0, quat1, quat2, t):
    ''' Quaternion splines. '''
    return quat_slerp_no_invert(quat_slerp_no_invert(quat0, quat2, t), quat_slerp_no_invert(quat0, quat1, t), 2 * t(1 - t))

def quat_to_matrix(quat):
    ''' Converts a quaternion to a rotational 4x4 matrix. '''
    x2 = quat.data[1] * quat.data[1]
    y2 = quat.data[2] * quat.data[2]
    z2 = quat.data[3] * quat.data[3]
    xy = quat.data[1] * quat.data[2]
    xz = quat.data[1] * quat.data[3]
    yz = quat.data[2] * quat.data[3]
    wx = quat.data[0] * quat.data[1]
    wy = quat.data[0] * quat.data[2]
    wz = quat.data[0] * quat.data[3]

    outputMatrix = matrix.Matrix(4)

    outputMatrix.matrix[0][0] = 1.0 - 2.0 * y2 - 2.0 * z2
    outputMatrix.matrix[0][1] = 2.0 * xy + 2.0 * wz
    outputMatrix.matrix[0][2] = 2.0 * xz - 2.0 * wy
    outputMatrix.matrix[0][3] = 0.0

    outputMatrix.matrix[1][0] = 2.0 * xy - 2.0 * wz
    outputMatrix.matrix[1][1] = 1.0 - 2.0 * x2 - 2.0 * z2
    outputMatrix.matrix[1][2] = 2.0 * yz + 2.0 * wx
    outputMatrix.matrix[1][3] = 0.0

    outputMatrix.matrix[2][0] = 2.0 * xz + 2.0 * wy
    outputMatrix.matrix[2][1] = 2.0 * yz - 2.0 * wx
    outputMatrix.matrix[2][2] = 1.0 - 2.0 * x2 - 2.0 * y2
    outputMatrix.matrix[2][3] = 0.0

    return outputMatrix

class Quaternion(object):

    def __init__(self, data=None):

        if data is None:
            self.data = quat_identity()
        else:
            self.data = data

    def __add__(self, other):
        if isinstance(other, Quaternion):
            return Quaternion(quat_add(self.data, other.data))
        else:
            return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, Quaternion):
            self.data = quat_add(self.data, other.data)
            return self
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Quaternion):
            return Quaternion(quat_sub(self.data, other.data))
        else:
            return NotImplemented

    def __isub__(self, other):
        if isinstance(other, Quaternion):
            self.data = quat_sub(self.data, other.data)
            return self
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            return Quaternion(quat_mul_quat(self.data, other.data))
        elif isinstance(other, vector.Vector):
            return Quaternion(quat_mul_vect(self.data, other.vector))
        elif isinstance(other, float):
            return Quaternion(quat_mul_float(self.data, other))
        else:
            return NotImplemented

    def __imul__(self, other):
        if isinstance(other, Quaternion):
            self.data = quat_mul_quat(self.data, other.data)
            return self
        elif isinstance(other, vector.Vector):
            self.data = quat_mul_vect(self.data, other.data)
            return self
        elif isinstance(other, float):
            self.data = quat_mul_float(self.data, other)
            return self
        else:
            return NotImplemented

    def __div__(self, other):
        if isinstance(other, float):
            return Quaternion(quat_div_float(self.data, other))
        else:
            return NotImplemented

    def __idiv__(self, other):
        if isinstance(other, float):
            self.data = quat_div_float(self.data, other)
            return self
        else:
            return NotImplemented

    def i_negate(self):
        self.data = quat_neg(self.data)
        return self

    def negate(self):
        quatList = quat_neg(self.data)
        return Quaternion(quatList)

    def i_identity(self):
        self.data = quat_identity()
        return self

    def identity(self):
        quatList = quat_identity()
        return Quaternion(quatList)

    def magnitude(self):
        return quat_magnitude(self.data)

    def dot(self, quat2):
        if isinstance(quat2, Quaternion):
            return quat_dot(self.data, quat2.data)
        else:
            return NotImplemented

    def i_normalize(self):
        self.data = quat_normalize(self.data)
        return self

    def normalize(self):
        quatList = quat_normalize(self.data)
        return Quaternion(quatList)

    def i_conjugate(self):
        self.data = quat_conjugate(self.data)
        return self

    def conjugate(self):
        quatList = quat_conjugate(self.data)
        return Quaternion(quatList)

    def inverse(self):
        quatList = quat_inverse(self.data)
        return Quaternion(quatList)

    def pow(self, e):
        exponent = e
        return quat_pow(self, exponent)

    def log(self):
        return quat_log(self)

    def lerp(self, quat1, time):
        return quat_lerp(self, quat1, time)

    def slerp(self, quat1, time):
        return quat_slerp(self, quat1, time)

    def slerp_no_invert(self, quat1, time):
        return quat_slerp_no_invert(self, quat1, time)

    def squad(self, quat1, quat2, time):
        return quat_squad(self, quat1, quat2, time)

    def toMatrix(self):
        return quat_to_matrix(self)

    # The following are used for orientation and motion
    def getForward(self):
        ''' Returns the forward vector. '''
        return quat_rotate_vector(self, vector.Vector(3, data=[0.0, 0.0, 1.0]))

    def getBack(self):
        ''' Returns the backwards vector. '''
        return quat_rotate_vector(self, vector.Vector(3, data=[0.0, 0.0, -1.0]))

    def getLeft(self):
        ''' Returns the left vector. '''
        return quat_rotate_vector(self, vector.Vector(3, data=[-1.0, 0.0, 0.0]))

    def getRight(self):
        ''' Returns the right vector. '''
        return quat_rotate_vector(self, vector.Vector(3, data=[1.0, 0.0, 0.0]))

    def getUp(self):
        ''' Returns the up vector. '''
        return quat_rotate_vector(self, vector.Vector(3, data=[0.0, 1.0, 0.0]))

    def getDown(self):
        ''' Returns the down vector. '''
        return quat_rotate_vector(self, vector.Vector(3, data=[0.0, -1.0, 0.0]))

def quat_from_matrix(matrix):
    ''' Converts a 4x4 rotational matrix to quaternion. '''
    fourXSquaredMinus1 = matrix.matrix[0][0] - matrix.matrix[1][1] - matrix.matrix[2][2]
    fourYSquaredMinus1 = matrix.matrix[1][1] - matrix.matrix[0][0] - matrix.matrix[2][2]
    fourZSquaredMinus1 = matrix.matrix[2][2] - matrix.matrix[0][0] - matrix.matrix[1][1]
    fourWSquaredMinus1 = matrix.matrix[0][0] + matrix.matrix[1][1] + matrix.matrix[2][2]

    biggestIndex = 0

    fourBiggestSquaredMinus1 = fourWSquaredMinus1

    if (fourXSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 1
    elif(fourYSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 2
    elif(fourZSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 3

    biggestVal = math.sqrt(fourBiggestSquaredMinus1 + 1) * 0.5
    mult = 0.25 / biggestVal

    rquat = Quaternion()

    if biggestIndex is 0:
        rquat.data[0] = biggestVal
        rquat.data[1] = (matrix.matrix[1][2] - matrix.matrix[2][1]) * mult
        rquat.data[2] = (matrix.matrix[2][0] - matrix.matrix[0][2]) * mult
        rquat.data[3] = (matrix.matrix[0][1] - matrix.matrix[1][0]) * mult
        return rquat

    if biggestIndex is 1:
        rquat.data[0] = (matrix.matrix[1][2] - matrix.matrix[2][1]) * mult
        rquat.data[1] = biggestVal
        rquat.data[2] = (matrix.matrix[0][1] + matrix.matrix[1][0]) * mult
        rquat.data[3] = (matrix.matrix[2][0] + matrix.matrix[0][2]) * mult
        return rquat

    if biggestIndex is 2:
        rquat.data[0] = (matrix.matrix[2][0] - matrix.matrix[0][2]) * mult
        rquat.data[1] = (matrix.matrix[0][1] + matrix.matrix[1][0]) * mult
        rquat.data[2] = biggestVal
        rquat.data[3] = (matrix.matrix[1][2] + matrix.matrix[2][1]) * mult
        return rquat

    if biggestIndex is 3:
        rquat.data[0] = (matrix.matrix[0][1] - matrix.matrix[1][0]) * mult
        rquat.data[1] = (matrix.matrix[2][0] + matrix.matrix[0][2]) * mult
        rquat.data[2] = (matrix.matrix[1][2] + matrix.matrix[2][1]) * mult
        rquat.data[3] = biggestVal
        return rquat
