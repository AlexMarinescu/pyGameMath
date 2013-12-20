import math
import src.library.math.vector as vec

try:
    range = xrange
except:
    pass

def Quaternion():
    ''' Easier to create a quaternion. '''
    return [0.0 for i in range(4)]

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
    ''' Multiply two quaternion and return a new one. '''
    return [quat[0] * quat1[1] + quat[1] * quat1[0] + quat[2] * quat1[3] - quat[3] * quat1[2],
            quat[0] * quat1[2] + quat[2] * quat1[0] + quat[3] * quat1[1] - quat[1] * quat1[3],
            quat[0] * quat1[3] + quat[3] * quat1[0] + quat[1] * quat1[2] - quat[2] * quat1[1],
            quat[0] * quat1[0] - quat[1] * quat1[1] - quat[2] * quat1[2] - quat[3] * quat1[3]]

def mulV(quat, vec):
    ''' Multiply a quaternion by a vector. Returns a 3D vector. '''
    vecN = vec.normalize(vec)

    vecQuat = Quaternion()
    resQuat = Quaternion()

    vecQuat[1] = vecN[0]
    vecQuat[2] = vecN[1]
    vecQuat[3] = vecN[2]

    resQuat = mulQ(vecQuat, conjugate(quat))
    resQuat = mulQ(resQuat, quat)

    return [resQuat[1], resQuat[2], resQuat[3]]

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

def rotate(axis, theta):
    ''' Returns a quaternion that rotates around a axis. '''
    thetaOver2 = theta * 0.5
    sinThetaOver2 = math.sin(thetaOver2)
    axisN = vec.normalize(axis)
    return [math.cos(thetaOver2), axis[0] * sinThetaOver2, axis[1] * sinThetaOver2, axis[2] * sinThetaOver2]

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

    return [[1.0 - 2.0 * (y2 + z2), 2.0 * (xy + wz), 2.0 * (xz - wy), 0.0],
            [2.0 * (xy - wz), 1.0 - 2.0 * (x2 + z2), 2.0 * (yz + wz), 0.0],
            [2.0 * (xz + wy), 2.0 * (yz - wx), 1.0 - 2.0 * (x2 + y2), 0.0],
            [0.0, 0.0, 0.0, 1.0]]

def fromMatrix(m):

    s = 0.0

    q = Quaternion()

    trace = m[0][0] + m[1][1] + m[2][2]

    if trace > 0.0:

        s = math.sqrt(trace + 1.0)

        q[3] = s * 0.5

        s = 0.5 / s

        q[0] = (m[1][2] - m[2][1]) * s
        q[1] = (m[2][0] - m[0][2]) * s 
        q[3] = (m[0][1] - m[1][0]) * s 
    else:
        nxt = [1, 2, 0]

        i = 0
        j = 0
        k = 0

        if (m[1][1] > m[0][0]):
            i = 1

        if (m[2][2] > m[i][i]):
            i = 2

        j = nxt[i]
        k = nxt[j]
        s = math.sqrt((m[i][i] - (m[j][j] + m[k][k])) + 1.0)

        q[i] = s * 0.5
        s = 0.5 / s
        q[3] = (m[j][k] - m[k][j]) * s
        q[j] = (m[i][j] + m[j][i]) * s 
        q[k] = (m[i][k] + m[k][i]) * s 

    return q

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