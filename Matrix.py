# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import math
import sys

from src.engine.math.constants import MAT3x3, MAT4x4
from src.engine.math.vector import Vector, cross


# NxN Matrix Class
class Matrix(object):
    '''Clss for NxN Matrices'''
    def __init__(self, array):
        # This will be needed when returning a same type array.
        self.size = len(array)
        # Calculate once to use in loops.
        self.sizeSquared = self.size * self.size

        # Setup the most often used variables.
        self.row = self.size    # Row of the matrix
        self.col = self.size    # Column of the matrix
        self.inr = self.size    # Inner of the matrix
        self.mat = array

        # Check to see if the array is the same as the size specified.
        if (self.size * self.size) == (self.size * len(self.mat[0])):
            pass
        else:
            print "Error: The matrix is not a square matrix."
            sys.exit(0)

        # Create an output array on init to save time(dynamic).
        self.out = [[0.0 for i in xrange(self.row)] for j in xrange(self.col)]

    # Slow; use the class.mat4x4[i][j] way.
    def __getitem__(self, key):
        return self.mat[key[0]][key[1]]

    # Slow; use the class.mat4x4[i][j] way.
    def __setitem__(self, key, val):
        self.mat[key[0]][key[1]] = val

    def __add__(self, input):
        for i in xrange(self.row):
            for j in xrange(self.col):
                self.out[i][j] = self.mat[i][j] + input.mat[i][j]

    def __iadd__(self, input):
        for i in xrange(self.row):
            for j in xrange(self.col):
                self.out[i][j] = self.mat[i][j] + input.mat[i][j]

    def __sub__(self, input):
        for i in xrange(self.row):
            for j in xrange(self.col):
                self.out[i][j] = self.mat[i][j] - input.mat[i][j]

    def __isub__(self, input):
        for i in xrange(self.row):
            for j in xrange(self.col):
                self.out[i][j] = self.mat[i][j] - input.mat[i][j]

    def __mul__(self, input):
        # Multiply by Vector
        if isinstance(input, Vector):
            vectorArray = []
            for x in xrange(input.len):
                vectorArray.append(0.0)
            outputVector = Vector(vectorArray)
            for i in xrange(self.row):
                for j in xrange(self.col):
                    outputVector.vec[i] += input.vec[j] * self.mat[i][j]
            return outputVector
        # Multiply by Matrix
        elif isinstance(input, Matrix):
            for i in xrange(self.row):
                for j in xrange(self.col):
                    for k in xrange(self.inr):
                        self.out[j][i] += self.mat[i][k] * input.mat[k][j]
            return Matrix(self.out)

    def __imul__(self, input):
        # Multiply by Vector
        if isinstance(input, Vector):
            vectorArray = []
            for x in xrange(input.len):
                vectorArray.append(0.0)
            outputVector = Vector(vectorArray)
            for i in xrange(self.row):
                for j in xrange(self.col):
                    outputVector.vec[i] += input.vec[j] * self.mat[i][j]
            return outputVector
        # Multiply by Matrix
        elif isinstance(input, Matrix):
            for i in xrange(self.row):
                for j in xrange(self.col):
                    for k in xrange(self.inr):
                        self.out[j][i] += self.mat[i][k] * input.mat[k][j]
            return Matrix(self.out)

    # Common, identical functions
    def duplicate(self):
        '''Returns a duplicate of itself.'''
        return Matrix(self.mat)

    def substitue(self, array):
        '''Replace the matrix array with a different one.'''
        self.mat = array

    def identity(self):
        '''Set the matrix to be an identity.'''
        for i in xrange(self.row):
            self.mat[i][i] = 1.0

    def zero(self):
        '''Zero the matrix.'''
        for i in xrange(self.row):
            for j in xrange(self.row):
                self.mat[i][j] = 0.0

    def transpose(self):
        '''Transpose the matrix.'''
        for i in xrange(self.row):
            for j in xrange(i + 1, self.col):
                self.mat[i][j] = self.mat[j][i]

    def scale(self, vector):
        '''Scale the matrix.'''
        for i, j in enumerate(vector.vec):
            self.mat[i][i] = j

    def inverse(self):
        '''Gauss Jordan Elimination for NxN Matrix inverse.'''
        matrix = self.mat  # Make a copy of itself
        # Create a dynamic container that can hold the inverse.
        # Basically create an augumented matrix. However, they are not in the same array.
        # The augumented matrix needs to be split at the end anyway so there is no point in joining them.
        container = [[0.0 for i in xrange(self.row)] for j in xrange(self.col)]
        # Identity matrix.
        for i in xrange(self.size):
            container[i][i] = 1.0
        # Find the inverse.
        for i in xrange(self.row):
            diag_number = float(matrix[i][i])
            for x in xrange(self.row):
                # Divide row 1 by first number, row 2 by second number and so on...
                # The number happens to be on the diagonal when its incremented.
                matrix[i][x] /= diag_number
                container[i][x] /= diag_number
            # Each time a row division happens. The above stuff; the following operations needs to occur.
            for j in xrange(self.row):
                # Make sure we are not doing calculations on the row we just divided.
                if j != i:
                    # Do the same thing with the number as the first step.
                    # However, this time make sure it is the updated one and it will no longer be on the diagonal, but on a column
                    # because its location is also given by the row number that is being calculated.
                    # So its no longer [n][n] where n represents the number location on the diagonal, its [row number][n].
                    number = float(matrix[j][i])
                    for y in xrange(self.row):
                        # Here the number multiplied by the numbers of the divided row is substracted from the original matrix.
                        # Example row1[i] - forth number on the row1 * row4[i]
                        # Where row1 is the one being calculated, row4 was the divided one.
                        matrix[j][y] = matrix[j][y] - number * matrix[i][y]
                        container[j][y] = container[j][y] - number * container[i][y]

        # Return a new matrix
        return Matrix(container)

    def pivot(M):
        '''Pivot the matrix'''
        size = M.row

        # Create identity matrix.
        identity = [[0.0 for i in xrange(size)] for j in xrange(size)]

        # Identity matrix.
        for i in xrange(size):
            identity[i][i] = 1.0

        for j in xrange(size):
            row = max(xrange(j, size), key=lambda i: M.mat[i][j])

            if j != row:
                identity[j], identity[row] = identity[row], identity[j]

        return Matrix(identity)

    def LUdecomposition(A):
        '''LU matrix decomposition for NxN Matrix. Returns matU, matL.'''
        size = A.row

        # Create L and U matrix
        U = [[0.0 for i in xrange(size)] for j in xrange(size)]
        L = [[0.0 for i in xrange(size)] for j in xrange(size)]

        # Identity for L matrix
        for i in xrange(size):
            L[i][i] = 1.0

        # Pivot the matrix
        P = Matrix.pivot(A)

        A2 = P * A

        for j in xrange(size):
            for i in xrange(j + 1):
                s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
                U[i][j] = A2.mat[i][j] - s1

            for i in xrange(j, size):
                s2 = sum(U[k][j] * L[i][k] for k in xrange(j))
                L[i][j] = (A2.mat[i][j] - s2) / U[j][j]

        matU = Matrix(U)
        matL = Matrix(L)

        return matU, matL

    def determinant(self):
        '''NxN Matrix determinant.'''
        matrix = self.duplicate()

        matU, matL = Matrix.LUdecomposition(matrix)

        det = 1.0

        # Multiply the diagonal
        for i in xrange(self.row):
                det *= matU.mat[i][i]

        return det * -1.0

    def normalize(self):
        '''NxN Matrix normalization.'''
        matrix = Matrix(self.mat)
        matdet = matrix.determinant()

        if matdet == 0:
            matdet = 1.0

        row, col = matrix.row

        normalized = [[0.0 for i in xrange(row)] for j in xrange(col)]

        for x in xrange(row):
            for y in xrange(col):
                normalized[x][y] = matrix.mat[x][y] / matdet

        return Matrix(normalized)

    def output(self):
        '''General matrix console output.'''
        print "Matrix", self.row, "x", self.col, ":"

        for i in xrange(self.size):
            print self.mat[i]


def mat2x2ShearX(x):
    '''Matrix 2x2 Shear on X.'''
    container = [[1.0, x],
                 [0.0, 1.0]]
    return Matrix(container)


def mat2x2ShearY(y):
    '''Matrix 2x2 Shear on Y.'''
    container = [[1.0, y],
                 [y, 1.0]]
    return Matrix(container)


def mat2x2Rotate(angle, direction):
    '''Matrix 2x2 Roate. Direction "R" or "L".'''
    container = [[1.0, 0.0],
                 [0.0, 1.0]]

    cosAngle = math.cos(angle)
    sinAngle = math.sin(angle)

    if direction == "R":
        container[0] = cosAngle
        container[1] = sinAngle
        container[2] = -sinAngle
        container[3] = cosAngle
    elif direction == "L":
        container[0] = cosAngle
        container[1] = -sinAngle
        container[2] = sinAngle
        container[3] = cosAngle
    else:
        print "Direction undefined."
    return Matrix(container)


def perspective(fov, aspect, znear, zfar):
    '''Perspective projection matrix 4x4.'''
    y = math.tan(fov * (3.14159265358979323846 / 360))
    x = y * aspect

    if (x != 0) and (y != 0) and (zfar != znear):
        a = 1.0 / x
        b = 1.0 / y
        c = (zfar) / (znear - zfar)
        d = (zfar * znear) / (znear - zfar)

        container = [[a, 0.0, 0.0, 0.0],
                     [0.0, b, 0.0, 0.0],
                     [0.0, 0.0, c, d],
                     [0.0, 0.0, 0.0, 1.0]]
    return Matrix(container)


def lookAt(cam_loc, direction, head_pos):
    '''Matrix 4x4 lookAt function.'''
    atNew = direction - cam_loc
    atNew.normalize()

    headNew = head_pos

    xaxis = cross(atNew, headNew)
    xaxis.normalize()

    upNew = cross(xaxis, atNew)

    container = [[xaxis.vec[0], upNew.vec[0], atNew.vec[0], cam_loc.vec[0]],
                 [xaxis.vec[1], upNew.vec[1], atNew.vec[1], cam_loc.vec[1]],
                 [xaxis.vec[2], upNew.vec[2], atNew.vec[2], cam_loc.vec[2]],
                 [0.0, 0.0, 0.0, 1.0]]
    return Matrix(container)


def translate(vector, matType):
    '''Translate by a vector. Second input: MAT4x4 or MAT3x3.'''
    if matType == MAT4x4:
        container = [[1.0, 0.0, 0.0, vector.vec[0]],
                     [0.0, 1.0, 0.0, vector.vec[1]],
                     [0.0, 0.0, 1.0, vector.vec[2]],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT3x3:
        container = [[1.0, 0.0, vector.vec[0]]
                     [0.0, 1.0, vector.vec[1]]
                     [0.0, 0.0, vector.vec[2]]]
        return Matrix(container)
    else:
        print "Not defined yet."


def rotateX(theta, matType):
    '''Rotate on X-axis. Second input: MAT4x4 or MAT3x3.'''
    c = math.cos(theta)
    s = math.sin(theta)
    if matType == MAT4x4:
        container = [[1.0, 0.0, 0.0, 0.0],
                     [0.0, c, -s, 0.0],
                     [0.0, s, c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT3x3:
        container = [[1.0, 0.0, 0.0],
                     [0.0, c, -s],
                     [0.0, s, c]]
        return Matrix(container)
    else:
        print "Not defined yet."


def rotateY(theta, matType):
    '''Rotate on Y-axis. Second input: MAT4x4 or MAT3x3.'''
    c = math.cos(theta)
    s = math.sin(theta)
    if matType == MAT4x4:
        container = [[c, 0.0, s, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [-s, 0.0, c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT3x3:
        container = [[c, 0.0, s],
                     [0.0, 1.0, 0.0],
                     [-s, 0.0, c]]
        return Matrix(container)
    else:
        print "Not defined yet."


def rotateZ(theta, matType):
    '''Rotate on Z-axis. Second input: MAT4x4 or MAT3x3.'''
    c = math.cos(theta)
    s = math.sin(theta)
    if matType == MAT4x4:
        container = [[c, -s, 0.0, 0.0],
                     [s, c, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT3x3:
        container = [[c, -s, 0.0],
                     [s, c, 0.0],
                     [0.0, 0.0, 1.0]]
        return Matrix(container)
    else:
        print "Not defined yet."


def rotate(axis, theta, matType):
    '''Rotate on arbitrary axis. Second input: MAT4x4 or MAT3x3.'''
    c = math.cos(theta)
    s = math.sin(theta)

    OneMinusCos = (1.0 - c)
    axis.normalize()

    if matType == MAT4x4:
        container = [[((axis.vec[0] * axis.vec[0]) * OneMinusCos) + c,                 ((axis.vec[0] * axis.vec[1]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[0] * axis.vec[2]) * OneMinusCos) + (axis.vec[1] * s), 0.0],
                     [((axis.vec[0] * axis.vec[1]) * OneMinusCos) + (axis.vec[2] * s), ((axis.vec[1] * axis.vec[1]) * OneMinusCos) + c,                 ((axis.vec[1] * axis.vec[2]) * OneMinusCos) - (axis.vec[0] * s), 0.0],
                     [((axis.vec[0] * axis.vec[2]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[1] * axis.vec[2]) * OneMinusCos) + (axis.vec[0] * s), ((axis.vec[2] * axis.vec[2]) * OneMinusCos) + c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT3x3:
        container = [[((axis.vec[0] * axis.vec[0]) * OneMinusCos) + c,                 ((axis.vec[0] * axis.vec[1]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[0] * axis.vec[2]) * OneMinusCos) + (axis.vec[1] * s)],
                     [((axis.vec[0] * axis.vec[1]) * OneMinusCos) + (axis.vec[2] * s), ((axis.vec[1] * axis.vec[1]) * OneMinusCos) + c,                 ((axis.vec[1] * axis.vec[2]) * OneMinusCos) - (axis.vec[0] * s)],
                     [((axis.vec[0] * axis.vec[2]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[1] * axis.vec[2]) * OneMinusCos) + (axis.vec[0] * s), ((axis.vec[2] * axis.vec[2]) * OneMinusCos) + c]]
        return Matrix(container)
    else:
        print "Not defined yet."


def shearXY(x,y,matType):
    '''Shear on XY. Second input: MAT4x4 or MAT3x3.'''
    if matType == MAT3x3:
        container = [[1.0, 0.0, x],
                     [0.0, 1.0, y],
                     [0.0, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT4x4:
        container = [[1.0, 0.0, 0.0, x],
                     [0.0, 1.0, 0.0, y],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    else:
        print "Not defined yet."


def shearYZ(x, y, matType):
    '''Shear on YZ. Second input: MAT4x4 or MAT3x3.'''
    if matType == MAT3x3:
        container = [[1.0, 0.0, 0.0],
                     [x, 1.0, 0.0],
                     [y, 0.0, 1.0]]
        return Matrix(container)
    elif matType == MAT4x4:
        container = [[1.0, 0.0, 0.0, 0.0],
                     [x, 1.0, 0.0, 0.0],
                     [y, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    else:
        print "Not defined yet."


def shearXZ(x, y, matType):
    '''Shear on XZ. Second input: MAT4x4 or MAT3x3.'''
    if matType == MAT3x3:
        container = [[1.0, x, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, y, 1.0]]
        return Matrix(container)
    elif matType == MAT4x4:
        container = [[1.0, x, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, y, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(container)
    else:
        print "Not defined yet."


# Matrix Stack Class
class MatrixStack(object):
    '''A matrix strack class.'''
    def __init__(self, matrix_list):
        self.stack = matrix_list

    def __getitem__(self, key):
        return self.stack[key]

    def __setitem__(self, key, value):
        self.stack[key] = value

    def push(self, matrix):
        self.stack.append(matrix)

    def assign(self, location, matrix):
        self.stack[location] = matrix

    def pop(self, *args):
        self.stack.pop(args[0])

    def output(self):
        print self.stack
