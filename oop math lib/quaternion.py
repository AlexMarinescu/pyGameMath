# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import math
import sys

from matrix import Matrix
from vector import Vector


# Quaternion class
# w = [0]
# x = [1]
# y = [2]
# z = [3]
class Quaternion(object):
    '''Class for Quaternions'''
    def __init__(self, array):
        self.quat = array
        self.len = len(array)

        # Make sure that it will not accept more than 4 values
        if self.len > 4:
            print "Quaternion: Wrong number of values."
            sys.exit(0)

    def duplicate(self):
        '''Returns a copy of itself'''
        return Quaternion(self.quat)

    def identity(self):
        '''Returns a identity quaternion'''
        return Quaternion([1.0, 0.0, 0.0, 0.0])

    def magnitude(self):
        '''Returns the magnitude of a quaternion'''
        result = 0
        for x in xrange(self.len):
            result += self.quat[x] * self.quat[x]
        return math.sqrt(result)

    def normalize(self):
        '''Normalizes a quaternion'''
        length = self.magnitude()
        if(length != 0):
            for x in xrange(self.len):
                self.quat[x] /= length
        else:
            print "Quaternion: Division by zero."
            sys.exit(0)

    def conjugate(self):
        '''Returns the conjugate'''
        new_quat = self.quat
        for x in xrange(self.len):
            new_quat[x] = -new_quat[x]
        new_quat[0] = -new_quat[0]
        return Quaternion(new_quat)

    def inverse(self):
        '''Returns the inverse of a quaternion'''
        lengthSquared = self.quat[0] * self.quat[0] + self.quat[1] * self.quat[1] + self.quat[2] * self.quat[2] + self.quat[3] * self.quat[3]

        return Quaternion([self.quat[0] / lengthSquared,
                           self.quat[1] / lengthSquared,
                           self.quat[2] / lengthSquared,
                           self.quat[3] / lengthSquared])

    def negate(self):
        '''Negate a quaternion'''
        self.quat[0] = -self.quat[0]
        self.quat[1] = -self.quat[1]
        self.quat[2] = -self.quat[2]
        self.quat[3] = -self.quat[3]

    def dot(self, quatA):
        '''Dot quaternions'''
        result = 0
        for x in xrange(self.len):
            result += self.quat[x] * quatA.quat[x]
        return result

    def rotateX(self, theta):
        '''Rotate around X with quaternion'''
        thetaOver2 = theta * 0.5
        return Quaternion([math.cos(thetaOver2),math.sin(thetaOver2),0.0,0.0])

    def rotateY(self, theta):
        '''Rotate around Y with quaternion'''
        thetaOver2 = theta * 0.5
        return Quaternion([math.cos(thetaOver2),0.0,math.sin(thetaOver2),0.0])

    def rotateZ(self, theta):
        '''Rotate around Z with quaternion'''
        thetaOver2 = theta * 0.5
        return Quaternion([math.cos(thetaOver2),0.0,0.0,math.sin(thetaOver2)])

    def rotate(self, axis, theta):
        '''Rotate around an axis with quaternion'''
        thetaOver2 = theta * 0.5
        sinThetaOver2 = math.sin(thetaOver2)
        axis.normalize()
        return Quaternion([math.cos(thetaOver2), axis.vec[1] * sinThetaOver2,
                           axis.vec[2] * sinThetaOver2, axis.vec[3] * sinThetaOver2])

    # Convert to a 4x4 matrix
    def toMatrix(self):
        '''Convert quaternion to a matrix.'''
        x2 = self.quat[1] * self.quat[1]
        y2 = self.quat[2] * self.quat[2]
        z2 = self.quat[3] * self.quat[3]
        xy = self.quat[1] * self.quat[2]
        xz = self.quat[1] * self.quat[3]
        yz = self.quat[2] * self.quat[3]
        wx = self.quat[0] * self.quat[1]
        wy = self.quat[0] * self.quat[2]
        wz = self.quat[0] * self.quat[3]

        return Matrix(16, [[1.0 - 2.0 * (y2 - z2), 2.0 * (xy - wz), 2.0 * (xz + wy), 0.0],
                           [2.0 * (xy + wz), 1.0 - 2.0 * (x2 - z2), 2.0 * (yz - wz), 0.0],
                           [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (x2 - y2), 0.0],
                           [0.0, 0.0, 0.0, 1.0]])

    # The log of a quaternion
    def log(self):
        '''Returns the log of a quaternion'''
        alpha = math.acos(self.quat[0])
        sinAlpha = math.sin(alpha)

        output = Quaternion([0.0,0.0,0.0,0.0])

        if sinAlpha > 0.0:
            output.quat[1] = self.quat[1] * alpha / sinAlpha
            output.quat[2] = self.quat[2] * alpha / sinAlpha
            output.quat[3] = self.quat[3] * alpha / sinAlpha
        else:
            output.quat[1] = self.quat[1]
            output.quat[2] = self.quat[2]
            output.quat[3] = self.quat[3]

        return output

    # The exp of a quaternion
    def exp(self):
        '''Returns exp of a quaternion'''
        alpha = 0
        for x in xrange(self.len):
            alpha += self.quat[x] * self.quat[x]

        alpha = math.sqrt(alpha)

        sinAlpha = math.sin(alpha)
        cosAlpha = math.cos(alpha)

        output = Quaternion([cosAlpha,0.0,0.0,0.0])

        if sinAlpha > 0.0:
            output.quat[1] = self.quat[1] * sinAlpha / alpha
            output.quat[2] = self.quat[2] * sinAlpha / alpha
            output.quat[3] = self.quat[3] * sinAlpha / alpha
        else:
            output.quat[1] = self.quat[1]
            output.quat[2] = self.quat[2]
            output.quat[3] = self.quat[3]

        return output

    # Rise a quaternion to a power
    # Usefull to extract a fraction of an angular displacement.
    def pow(self, exponent):
        '''Rise a quaternion to a power'''
        if self.quat[0] != 0.0:
            # Extract the half angle
            angle = math.acos(self.quat[0])
            newAngle = angle * exponent
            # Compute a new W
            self.quat[0] = math.cos(newAngle)
            # Compute new XYZ values
            angleDiv = math.sin(newAngle) / math.sin(angle)
            self.quat[1] *= angleDiv
            self.quat[2] *= angleDiv
            self.quat[3] *= angleDiv

    # Linear Interpolation
    def lerp(self, quat0, quat1, t):
        '''Quaternion linear interpolation'''
        k0 = 1.0 - t
        k1 = t

        output = Quaternion([1.0,0.0,0.0,0.0])

        output.quat[0] = quat0.quat[0] * k0 + quat1.quat[0] * k1
        output.quat[1] = quat0.quat[1] * k0 + quat1.quat[1] * k1
        output.quat[2] = quat0.quat[2] * k0 + quat1.quat[2] * k1
        output.quat[3] = quat0.quat[3] * k0 + quat1.quat[3] * k1

        return output()


    # Quaternion Interpolation (SLERP)
    # The first two inputs are two quaternions to interpolate between.
    # The third input is the interpolation parameter.
    # slerp(q0,q1,t) = q0(q0^-1 * q1)^t
    def slerp(self, quat0, quat1, t):
        '''Quaternion  spherical interpolation'''
        output = Quaternion([1.0,0.0,0.0,0.0])

        # Compute the cosine of the angle
        cosTheta = quat0.dot(quat1)

        # If negative dot, negate one of the input quaternions for a shorter arc
        if cosTheta < 0.0:
            quat1.negate()
            cosTheta = -cosTheta

        # Interpolation params
        k0 = 0.0
        k1 = 0.0
        # Check for divions by zero
        if cosTheta > 0.999:
            # Too close, liner interpolation
            k0 = 1.0 - t
            k1 = t
        else:
            # Compute the sin of angle
            sinTheta = math.sqrt(1.0 - Theta * Theta)
            # Compute angle from the sin and cosine
            theta = math.atan2(sinTheta, cosTheta)
            # Compute inverse denominator
            oneOverSinTheta = 1.0 / sinTheta
            # Compute interpolation params
            k0 = math.sin((1.0 - t) * theta) * oneOverSinTheta
            k1 = math.sin(t * theta) * oneOverSinOmega

        # Interpolate
        output.quat[0] = quat0.quat[0] * k0 + quat1.quat[0] * k1
        output.quat[1] = quat0.quat[1] * k0 + quat1.quat[1] * k1
        output.quat[2] = quat0.quat[2] * k0 + quat1.quat[2] * k1
        output.quat[3] = quat0.quat[3] * k0 + quat1.quat[3] * k1

        return output

    # SQUAD (Quaternion Splines)
    def squad(self, quat1, q0, q1, q2, t):
        '''Quaternion splines'''
        return self.slerp(self.slerp(quat1, q2, t), self.slerp(q0, q1, t), 2*t*(1-t))


    # Overload !=
    def __ne__(self, quat):
        for i in xrange(len(self.quat)):
            if self.quat[i] != quat.quat[i]:
                return False
        return True

    # Overload ==
    def __eq__(self, quat):
        for i in xrange(len(self.quat)):
            if self.quat[i] != quat.quat[i]:
                return False
        return True

    # Overload /
    def __div__(self, input):
        return Quaternion([self.quat[0] / input, self.quat[1] / input, self.quat[2] /input, self.quat[3] / input])
    
    # Overload /=
    def __idiv__(self, input):
        return Quaternion([self.quat[0] / input, self.quat[1] / input, self.quat[2] /input, self.quat[3] / input])

    # Overload *
    def __mul__(self, input):
        # Cross Product
        if isinstance(input, Quaternion):
            return Quaternion([self.quat[0] * input.quat[1] + self.quat[1] * input.quat[0] + self.quat[2] * input.quat[3] - self.quat[3] * input.quat[2],
                               self.quat[0] * input.quat[2] + self.quat[2] * input.quat[0] + self.quat[3] * input.quat[1] - self.quat[1] * input.quat[3],
                               self.quat[0] * input.quat[3] + self.quat[3] * input.quat[0] + self.quat[1] * input.quat[2] - self.quat[2] * input.quat[1],
                               self.quat[0] * input.quat[0] - self.quat[1] * input.quat[1] - self.quat[2] * input.quat[2] - self.quat[3] * input.quat[3]])

        # Vector 3D * Quaternion
        elif isinstance(input, Vector):
            input.normalize()
            vecQuat = Quaternion([1.0, 0.0, 0.0, 0.0])
            resQuat = Quaternion([1.0, 0.0, 0.0, 0.0])

            vecQuat.quat[1] = input.vec[0]
            vecQuat.quat[2] = input.vec[1]
            vecQuat.quat[3] = input.vec[2]
            vecQuat.quat[0] = 0.0

            resQuat = vecQuat * self.conjugate()
            resQuat = self * resQuat

            return Vector([resQuat.quat[1], resQuat.quat[2], resQuat.quat[3]])

        # Multiply by a scalar
        elif isinstance(input, float):
            return Quaternion([self.quat[0] * input, self.quat[1] * input, self.quat[2] * input, self.quat[3] * input])

    # Overload *=
    def __imul__(self, input):
        # Cross Product
        if isinstance(input, Quaternion):
            return Quaternion([self.quat[0] * input.quat[1] + self.quat[1] * input.quat[0] + self.quat[2] * input.quat[3] - self.quat[3] * input.quat[2],
                               self.quat[0] * input.quat[2] + self.quat[2] * input.quat[0] + self.quat[3] * input.quat[1] - self.quat[1] * input.quat[3],
                               self.quat[0] * input.quat[3] + self.quat[3] * input.quat[0] + self.quat[1] * input.quat[2] - self.quat[2] * input.quat[1],
                               self.quat[0] * input.quat[0] - self.quat[1] * input.quat[1] - self.quat[2] * input.quat[2] - self.quat[3] * input.quat[3]])

        # Vector 3D * Quaternion
        elif isinstance(input, Vector):
            input.normalize()
            vecQuat = Quaternion([1.0, 0.0, 0.0, 0.0])
            resQuat = Quaternion([1.0, 0.0, 0.0, 0.0])

            vecQuat.quat[1] = input.vec[0]
            vecQuat.quat[2] = input.vec[1]
            vecQuat.quat[3] = input.vec[2]
            vecQuat.quat[0] = 0.0

            resQuat = vecQuat * self.conjugate()
            resQuat = self * resQuat

            return Vector([resQuat.quat[1], resQuat.quat[2], resQuat.quat[3]])

        # Multiply by a scalar
        elif isinstance(input, float):
            return Quaternion([self.quat[0] * input, self.quat[1] * input, self.quat[2] * input, self.quat[3] * input])

    # Overload +
    def __add__(self, quat):
        return Quaternion([self.quat[0] + quat.quat[0], self.quat[1] + quat.quat[1], self.quat[2] + quat.quat[2], self.quat[3] + quat.quat[3]])

    # Overload +=
    def __iadd__(self, quat):
        return Quaternion([self.quat[0] + quat.quat[0], self.quat[1] + quat.quat[1], self.quat[2] + quat.quat[2], self.quat[3] + quat.quat[3]])

    # Overload -
    def __sub__(self, quat):
        return Quaternion([self.quat[0] - quat.quat[0], self.quat[1] - quat.quat[1], self.quat[2] - quat.quat[2], self.quat[3] - quat.quat[3]])

    # Overload -=
    def __isub__(self, quat):
        return Quaternion([self.quat[0] - quat.quat[0], self.quat[1] - quat.quat[1], self.quat[2] - quat.quat[2], self.quat[3] - quat.quat[3]])

    def output(self):
        '''Prints a quaternion to screen'''
        print "Quaternion: ", self.quat

