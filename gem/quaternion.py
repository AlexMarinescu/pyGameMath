import math
import six.moves as sm
from gem import vector
from gem import matrix
from gem import common

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
            -quat[1] / lengthSquared,
            -quat[2] / lengthSquared,
            -quat[3] / lengthSquared]

def quat_from_axis_angle(axis, theta):
    ''' Return a rotation Quaternion from a Vector/list axis and degrees.

    Normalizes a temporary axis without mutating caller data.
    '''
    thetaOver2 = theta * 0.5
    sto2 = math.sin(math.radians(thetaOver2))
    cto2 = math.cos(math.radians(thetaOver2))

    if isinstance(axis, vector.Vector):
        naxis = axis.normalize()
    elif isinstance(axis, list):
        naxis = vector.Vector(3, data=axis).normalize()
    else:
        return NotImplemented

    quat1List = [cto2, naxis.vector[0] * sto2, naxis.vector[1] * sto2, naxis.vector[2] * sto2]
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
    ''' Return the legacy Quaternion product rotating the normalized axis.

    Accepts a Vector/list axis and degrees without mutating caller data.
    The result is the pure Quaternion approximately [0, normalized_axis],
    not an axis-angle rotation constructor. Use quat_from_axis_angle to
    construct a rotation Quaternion. The numerical sandwich is retained.
    '''
    thetaOver2 = theta * 0.5
    sto2 = math.sin(math.radians(thetaOver2))
    cto2 = math.cos(math.radians(thetaOver2))

    if isinstance(axis, vector.Vector):
        naxis = axis.normalize()
    elif isinstance(axis, list):
        naxis = vector.Vector(3, data=axis).normalize()
    else:
        return NotImplemented

    quat1List = [cto2, naxis.vector[0] * sto2, naxis.vector[1] * sto2, naxis.vector[2] * sto2]
    quat1 = Quaternion(data=quat1List)
    rotation = (quat1 * naxis) * quat1.conjugate()
    return rotation

def quat_rotate_vector(quat, vec):
    ''' Rotates a vector by a quaternion, returns a vector. '''
    outQuat = (quat * vec) * quat.conjugate()
    return vector.Vector(3, data=[outQuat.data[1], outQuat.data[2], outQuat.data[3]])

def quat_pow(quat, exp):
    """Return a fresh unit-quaternion power using the principal angle.

    Nonunit inputs are unsupported. Negative identity has no unique axis,
    so only integer powers are defined for it.
    """
    w, x, y, z = quat.data
    imaginary = math.hypot(math.hypot(x, y), z)
    if imaginary == 0.0:
        if w == 0.0:
            raise ValueError("Zero quaternion has no unit-quaternion power")
        if w < 0.0:
            if exp % 1 != 0:
                raise ValueError("Negative identity requires an integer power")
            return Quaternion(data=[-1.0 if exp % 2 else 1.0, 0.0, 0.0, 0.0])
        return Quaternion()
    if exp == 0:
        return Quaternion()
    if exp == 1:
        return Quaternion(data=list(quat.data))
    angle = math.atan2(imaginary, w)
    powered_angle = angle * exp
    # Reduce the exponent first only when multiplication would overflow.
    if math.isinf(powered_angle) and not math.isinf(exp):
        powered_angle = math.fmod(exp, (2.0 * math.pi) / angle) * angle
    sine = math.sin(powered_angle)
    return Quaternion(data=[math.cos(powered_angle),
                            (x / imaginary) * sine,
                            (y / imaginary) * sine,
                            (z / imaginary) * sine])

def quat_log(quat):
    """Return a fresh [0, axis * principal angle] list for a unit quaternion.

    Nonunit inputs are unsupported; zero and negative identity raise
    ValueError. No sign canonicalization or imaginary-axis cutoff is used.
    """
    w, x, y, z = quat.data
    imaginary = math.hypot(math.hypot(x, y), z)
    if imaginary == 0.0:
        if w <= 0.0:
            raise ValueError("Quaternion logarithm has no unique imaginary axis")
        return [0.0, 0.0, 0.0, 0.0]
    angle = math.atan2(imaginary, w)
    return [0.0, (x / imaginary) * angle,
            (y / imaginary) * angle, (z / imaginary) * angle]

def quat_lerp(quat0, quat1, t):
    ''' Linear interpolation between two quaternions. '''
    k0 = 1.0 - t
    k1 = t

    output = Quaternion()
    output = (quat0 * k0) + (quat1 * k1)

    return output

def quat_slerp(quat0, quat1, t):
    """Return fresh shortest-path spherical interpolation of unit inputs.

    No input normalization or parameter clamping is performed. Equal
    orientations use the continuous limit; signs follow the first input.
    """
    end = quat1.negate() if quat0.dot(quat1) < 0.0 else quat1
    # Difference/sum norms retain angles too small for acos(dot) to resolve.
    difference = [a - b for a, b in zip(quat0.data, end.data)]
    total = [a + b for a, b in zip(quat0.data, end.data)]
    difference_norm = math.hypot(math.hypot(*difference[:2]), math.hypot(*difference[2:]))
    total_norm = math.hypot(math.hypot(*total[:2]), math.hypot(*total[2:]))
    theta = 2.0 * math.atan2(difference_norm, total_norm)
    if theta == 0.0:
        return Quaternion(data=list(quat0.data))
    denominator = math.sin(theta)
    k0 = math.sin((1.0 - t) * theta) / denominator
    k1 = math.sin(t * theta) / denominator
    return (quat0 * k0) + (end * k1)

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
    """Legacy three-control blend: quat0=start, quat2=end, quat1=control.

    Uses sign-sensitive slerp_no_invert, including its linear approximations;
    unit length is not guaranteed. See squad4 for conventional SQUAD.
    """
    a = quat_slerp_no_invert(quat0, quat2, t)
    b = quat_slerp_no_invert(quat0, quat1, t)
    return quat_slerp_no_invert(a, b, 2 * t * (1 - t))

def squad4(q0, q1, s0, s1, t):
    """Return conventional four-control SQUAD for unit quaternions.

    q0/q1 are endpoints; s0/s1 are SQUAD controls, not neighbouring
    keyframes. Uses accurate shortest-path SLERP; t is ordinarily in [0,1].
    Inputs and their storage are preserved; the result is a fresh Quaternion.
    """
    a = quat_slerp(q0, q1, t)
    b = quat_slerp(s0, s1, t)
    return quat_slerp(a, b, 2 * t * (1 - t))

def quat_to_matrix(quat):
    ''' Return a row-vector Matrix4 for a unit [w,x,y,z] quaternion.

    No implicit normalization is performed. The quaternion is preserved,
    and the returned float32 ctypes snapshot matches its matrix rows.
    '''
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

    outputMatrix.c_matrix = common.conv_list_2d(outputMatrix.matrix, common.GLfloat)
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
            self.data = quat_mul_vect(self.data, other.vector)
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

    __truediv__ = __div__
    __itruediv__ = __idiv__

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
    ''' Return a [w,x,y,z] Quaternion from a proper row-vector rotation.

    Uses the input Matrix's upper-left 3x3 block without mutating it.
    Quaternion sign is not unique; no normalization or validation is added.
    '''
    fourXSquaredMinus1 = matrix.matrix[0][0] - matrix.matrix[1][1] - matrix.matrix[2][2]
    fourYSquaredMinus1 = matrix.matrix[1][1] - matrix.matrix[0][0] - matrix.matrix[2][2]
    fourZSquaredMinus1 = matrix.matrix[2][2] - matrix.matrix[0][0] - matrix.matrix[1][1]
    fourWSquaredMinus1 = matrix.matrix[0][0] + matrix.matrix[1][1] + matrix.matrix[2][2]

    biggestIndex = 0

    fourBiggestSquaredMinus1 = fourWSquaredMinus1

    if (fourXSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 1
        fourBiggestSquaredMinus1 = fourXSquaredMinus1
    if (fourYSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 2
        fourBiggestSquaredMinus1 = fourYSquaredMinus1
    if (fourZSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 3
        fourBiggestSquaredMinus1 = fourZSquaredMinus1

    biggestVal = math.sqrt(fourBiggestSquaredMinus1 + 1) * 0.5
    mult = 0.25 / biggestVal

    rquat = Quaternion()

    if biggestIndex == 0:
        rquat.data[0] = biggestVal
        rquat.data[1] = (matrix.matrix[1][2] - matrix.matrix[2][1]) * mult
        rquat.data[2] = (matrix.matrix[2][0] - matrix.matrix[0][2]) * mult
        rquat.data[3] = (matrix.matrix[0][1] - matrix.matrix[1][0]) * mult
        return rquat

    if biggestIndex == 1:
        rquat.data[0] = (matrix.matrix[1][2] - matrix.matrix[2][1]) * mult
        rquat.data[1] = biggestVal
        rquat.data[2] = (matrix.matrix[0][1] + matrix.matrix[1][0]) * mult
        rquat.data[3] = (matrix.matrix[2][0] + matrix.matrix[0][2]) * mult
        return rquat

    if biggestIndex == 2:
        rquat.data[0] = (matrix.matrix[2][0] - matrix.matrix[0][2]) * mult
        rquat.data[1] = (matrix.matrix[0][1] + matrix.matrix[1][0]) * mult
        rquat.data[2] = biggestVal
        rquat.data[3] = (matrix.matrix[1][2] + matrix.matrix[2][1]) * mult
        return rquat

    if biggestIndex == 3:
        rquat.data[0] = (matrix.matrix[0][1] - matrix.matrix[1][0]) * mult
        rquat.data[1] = (matrix.matrix[2][0] + matrix.matrix[0][2]) * mult
        rquat.data[2] = (matrix.matrix[1][2] + matrix.matrix[2][1]) * mult
        rquat.data[3] = biggestVal
        return rquat
