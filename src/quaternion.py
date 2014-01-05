import math
import src.library.math.vector as vec
import src.library.math.matrix as mat
import src.library.math.common as com

try:
    range = xrange
except:
    pass

def Quaternion():
    ''' Easier to create a quaternion. '''
    return [1.0, 0.0, 0.0, 0.0]

def add(quat, quat1):
    ''' Add two quaternions. '''
    return [quat[0] + quat1[0], quat[1] + quat1[1], quat[2] + quat1[2], quat[3] + quat1[3]]

def sub(quat, quat1):
    ''' Substract two quaternions. '''
    return [quat[0] - quat1[0], quat[1] - quat1[1], quat[2] - quat1[2], quat[3] - quat1[3]]

def div(quat, nr):
    ''' Divide a quaternion by a number. '''
    return [quat[0] / nr, quat[1] / nr, quat[2] / nr, quat[3] / nr]

def mulQ(quat, quat1):
    ''' Multiply two quaternions and return a new one. '''
    w = quat[0] * quat1[0] - quat[1] * quat1[1] - quat[2] * quat1[2] - quat[3] * quat1[3]
    x = quat[0] * quat1[1] + quat[1] * quat1[0] + quat[2] * quat1[3] - quat[3] * quat1[2]
    y = quat[0] * quat1[2] + quat[2] * quat1[0] + quat[3] * quat1[1] - quat[1] * quat1[3]
    z = quat[0] * quat1[3] + quat[3] * quat1[0] + quat[1] * quat1[2] - quat[2] * quat1[1]
    return [w, x, y, z]

def mulV(quat, vect):
    ''' Multiply a quaternion by a vector. Returns a 3D vector. '''
    w = -quat[1] * vect[0] - quat[2] * vect[1] - quat[3] * vect[2]
    x =  quat[0] * vect[0] + quat[2] * vect[2] - quat[3] * vect[1]
    y =  quat[0] * vect[1] + quat[3] * vect[0] - quat[1] * vect[2]
    z =  quat[0] * vect[2] + quat[1] * vect[1] - quat[2] * vect[0]
    return [w, x, y, z]

def mulS(quat, nr):
    ''' Multiply a quaternion by a number. '''
    return [quat[0] / nr, quat[1] / nr, quat[2] / nr, quat[3] / nr]

def magnitude(quat):
    ''' Return the magnitude of a quaternion. '''
    result = 0.0
    for x in range(4):
        result += quat[x] * quat[x]
    return math.sqrt(result)

def normalize(quat):
    ''' Returns a normalized quaternion. '''
    length = magnitude(quat)
    output = Quaternion()
    if(length != 0):
        for x in range(4):
            output[x] = quat[x] / length
    return output

def conjugate(quat):
    ''' Returns the conjugate of a quaternion. '''
    newQ = Quaternion()
    for x in range(4):
        newQ[x] = -quat[x]
    newQ[0] = -newQ[0]
    return newQ

def inverse(quat):
    ''' Returns the inverse of a quaternion. '''
    lengthSquared = quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]

    return [quat[0] / lengthSquared,
            quat[1] / lengthSquared,
            quat[2] / lengthSquared,
            quat[3] / lengthSquared]

def negate(quat):
    ''' Negate a quaternion. '''
    return [-quat[0], -quat[1], -quat[2], -quat[3]]

def dot(quatA, quatB):
    ''' Dot product between two quaternions. Returns a scalar. '''
    result = 0
    for i in range(4):
        result += quatA[x] * quatB[x]
    return result

def rotateX(theta):
    ''' Returns a quaternion that rotates around X axis. '''
    thetaOver2 = theta * 0.5
    return [math.cos(thetaOver2), math.sin(thetaOver2), 0.0, 0.0]

def rotateY(theta):
    ''' Returns a quaternion that rotates around Y axis. '''
    thetaOver2 = theta * 0.5
    return [math.cos(thetaOver2), 0.0, math.sin(thetaOver2), 0.0]

def rotateZ(theta):
    ''' Returns a quaternion that rotates around Z axis. '''
    thetaOver2 = theta * 0.5
    return [math.cos(thetaOver2), 0.0, 0.0, math.sin(thetaOver2)]

def angleAxis(theta, axis):
    ''' Returns a quaternion from angle and axis. '''
    halfAngle = theta * 0.5
    sHAngle = math.sin(halfAngle)
    cHAngle = math.cos(halfAngle)
    nAxis = vec.normalize(axis)
    return [cHAngle, nAxis[0] * sHAngle, nAxis[1] * sHAngle, nAxis[2] * sHAngle]

def rotate(origin, axis, theta):
    ''' Returns a quaternion that rotates around a axis. '''
    thetaOver2 = theta * 0.5
    sinThetaOver2 = math.sin(math.radians(thetaOver2))
    cosThetaOver2 = math.cos(math.radians(thetaOver2))
    Q1 = [cosThetaOver2, axis[0] * sinThetaOver2, axis[1] * sinThetaOver2, axis[2] * sinThetaOver2]
    rotation = mulQ(mulV(Q1, origin), conjugate(Q1))
    return [rotation[1], rotation[2], rotation[3]]

def rotateQ(quat, axis):
    ''' Rotates a vector by a quaternion. '''
    l = len(axis)
    if l == 4:
        return mulQ(mulQ(quat, axis), conjugate(quat))
    else:
        final = mulQ(mulV(quat, axis), conjugate(quat))
        return [final[1], final[2], final[3]]

def toMatrix(quat):
    ''' Converts a quaternion to a rotational 4x4 matrix. '''
    x2 = quat[1] * quat[1]
    y2 = quat[2] * quat[2]
    z2 = quat[3] * quat[3]
    xy = quat[1] * quat[2]
    xz = quat[1] * quat[3]
    yz = quat[2] * quat[3]
    wx = quat[0] * quat[1]
    wy = quat[0] * quat[2]
    wz = quat[0] * quat[3]

    result = mat.Matrix(4)

    result[0][0] = 1.0 - 2.0 * y2 - 2.0 * z2
    result[0][1] = 2.0 * xy + 2.0 * wz
    result[0][2] = 2.0 * xz - 2.0 * wy
    result[0][3] = 0.0

    result[1][0] = 2.0 * xy - 2.0 * wz
    result[1][1] = 1.0 - 2.0 * x2 - 2.0 * z2
    result[1][2] = 2.0 * yz + 2.0 * wx
    result[1][3] = 0.0

    result[2][0] = 2.0 * xz + 2.0 * wy
    result[2][1] = 2.0 * yz - 2.0 * wx
    result[2][2] = 1.0 - 2.0 * x2  - 2.0 * y2
    result[2][3] = 0.0

    result[3][0] = 0.0
    result[3][1] = 0.0
    result[3][2] = 0.0
    result[3][3] = 1.0

    return result

def fromMatrix(m):
    ''' Converts a matrix to a quaternion. '''
    fourXSquaredMinus1 = m[0][0] - m[1][1] - m[2][2]
    fourYSquaredMinus1 = m[1][1] - m[0][0] - m[2][2]
    fourZSquaredMinus1 = m[2][2] - m[0][0] - m[1][1]
    fourWSquaredMinus1 = m[0][0] + m[1][1] + m[2][2]

    biggestIndex = 0

    fourBiggestSquaredMinus1 = fourWSquaredMinus1

    if(fourXSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 1
    elif(fourYSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 2
    elif(fourZSquaredMinus1 > fourBiggestSquaredMinus1):
        biggestIndex = 3

    biggestVal = math.sqrt(fourBiggestSquaredMinus1 + 1) * 0.5
    mult = 0.25 / biggestVal

    result = Quaternion()

    if biggestIndex == 0:
        result[0] = biggestVal
        result[1] = (m[1][2] - m[2][1]) * mult
        result[2] = (m[2][0] - m[0][2]) * mult
        result[3] = (m[0][1] - m[1][0]) * mult
        return result

    if biggestIndex == 1:
        result[0] = (m[1][2] - m[2][1]) * mult
        result[1] = biggestVal
        result[2] = (m[0][1] + m[1][0]) * mult
        result[3] = (m[2][0] + m[0][2]) * mult
        return result

    if biggestIndex == 2:
        result[0] = (m[2][0] - m[0][2]) * mult
        result[1] = (m[0][1] + m[1][0]) * mult
        result[2] = biggestVal
        result[3] = (m[1][2] + m[2][1]) * mult
        return result

    if biggestIndex == 3:
        result[0] = (m[0][1] - m[1][0]) * mult
        result[1] = (m[2][0] + m[0][2]) * mult
        result[2] = (m[1][2] + m[2][1]) * mult
        result[3] = biggestVal
        return result

def log(quat):
    ''' Return the logarithm of a quaternion. '''
    alpha = math.acos(quat[0])
    sinAlpha = math.sin(alpha)

    output = Quaternion()

    if sinAlpha > 0.0:
        output[1] = quat[1] * alpha / sinAlpha
        output[2] = quat[2] * alpha / sinAlpha
        output[3] = quat[3] * alpha / sinAlpha
    else:
        output = quat

    return output

def pow(quat, exponent):
    ''' Returns a quaternion to the power of x. '''
    quatN = Quaternion()
    if quat[0] != 0.0:
        angle = math.acos(quat[0])
        newAngle = angle * exponent
        quatN[0] = math.cos(newAngle)
        angleDiv = math.sin(newAngle) / math.sin(angle)
        quatN[1] *= angleDiv
        quatN[2] *= angleDiv
        quatN[3] *= angleDiv
    return quatN

def lerp(quat0, quat1, t):
    ''' Linear iterpolation between two quaternions. '''
    k0 = 1.0 - t
    k1 = t

    output = Quaternion()

    output[0] = quat0[0] * k0 + quat1[0] * k1
    output[1] = quat0[1] * k0 + quat1[1] * k1
    output[2] = quat0[2] * k0 + quat1[2] * k1
    output[3] = quat0[3] * k0 + quat1[3] * k1

    return output

def slerp(quat0, quat1, t):
    ''' Spherical interpolation between two quaternions. '''
    k0 = 0.0
    k1 = 0.0
    output = Quaternion()
    quat1Neg = Quaternion()
    cosTheta = dot(quat0, quat1)

    if cosTheta < 0.0:
        quat1Neg = negate(quat1)
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

    output[0] = quat0[0] * k0 + quat1Neg[0] * k1
    output[1] = quat0[1] * k0 + quat1Neg[1] * k1
    output[2] = quat0[2] * k0 + quat1Neg[2] * k1
    output[3] = quat0[3] * k0 + quat1Neg[3] * k1

    return output

def slerpNoInvert(quat0, quat1, t):
    ''' Version of SLERP, used by SQUAD, it doesn't check for theta > 90. '''
    dotP = dot(quat0, quat1)
    output = Quaternion()

    if dot > -0.95 and dot < 0.95:
        angle = math.acos(dot)
        k0 = math.sin(angle * (1.0 - t)) / math.sin(angle)
        k1 = math.sin(t * angle) / math.sin(angle)

        output[0] = quat0[0] * k0 + quat1Neg[0] * k1
        output[1] = quat0[1] * k0 + quat1Neg[1] * k1
        output[2] = quat0[2] * k0 + quat1Neg[2] * k1
        output[3] = quat0[3] * k0 + quat1Neg[3] * k1
    else:
        output = lerp(quat0, quat1, t)

    return output

def squad(quat1, q1, q2, t):
    ''' Quaternion splines. '''
    return slerpNoInvert(slerpNoInvert(quat1, q2, t), slerpNoInvert(q0, q1, t), 2 * t * (1 - t))