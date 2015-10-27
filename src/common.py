import math
import six.moves as sm
import ctypes as ct

#OpenGL Types, this is needed because Matrix class is supported by OpenGL
GLfloat = ct.c_float

def convertArr(l, n):
    ''' Convert a 1D array to 2D array. '''
    return [l[i:i+n] for i in sm.range(0, len(l), n)]

# Multiply each element in the 4 length arrays eg. [a, b, c, d] x [e, f, g, h] = [axe, bxf, cxg, dxh]
def mulV4(v1, v2):
    ''' Multiply two 1D 4 Length Arrays.'''
    return [v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2], v1[3] * v2[3]]

def conv_list(listIn, cType):
    ''' Convert a python list into a ctypes array '''
    return (cType * len(listIn))(*listIn)

def conv_list_2d(listIn, cType):
    ''' Convert a python 2d list into a ctypes 2d array '''
    xlength = len(listIn)
    ylength = len(listIn[0])

    arrayOut = (cType * ylength * xlength)()

    for x in sm.range(xlength):
        for y in sm.range(ylength):
            arrayOut[x][y] = listIn[x][y]

    return arrayOut

def list_2d_to_1d(inlist):
    ''' Convert a python 2D list to python 1D list. '''
    sizeX = len(inlist[0])
    sizeY = len(inlist)

    rtnList = [None for x in sm.range(sizeY*sizeX)]

    for x in sm.range(sizeY):
        for y in sm.range(sizeX):
            rtnList[y + x * sizeX] = inlist[x][y]

    return rtnList

def convertM4to3(matrix):
    ''' Convert a 4x4 Matrix to 3x3. '''
    temp = [[0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]]
            
    for i in sm.range(3):
        for j in sm.range(3):
            temp[i][j] = matrix[i][j]
    return temp

# Sinc function
def sinc(x):
    if math.fabs(x) < 1.0e-4:
        return 1.0
    else:
        return math.sin(x) / x

# Scalar lerp
def scalarLerp(a, b, time):
    ''' Iterpolate a number over time. '''
    return a + time * (b - a)

# Returns the view port coordinates
def getViewPort(coords, width, height):
    ''' A version of glViewPort except it returns the coords. '''
    coordsN = coords.normalize()
    x = (coordsN[0] + 1) * (width / 2) + coords[0]
    y = (coordsN[1] + 1) * (height / 2) + coords[1]
    return [x,y,width,height]

# Radians to degree
def radiansToDegrees(degrees):
    return (degrees * 3.14) / 180.0

# Degrees to radians
def degreesToRadians(radians):
    return (radians * 180.0) / 3.14

# Returns the sign
def sign(x):
    if x >= 0.0:
        return 1.0
    else:
        return -1.0