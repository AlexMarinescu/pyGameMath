import src.library.math.vector as vec
import math

# Convert 1D to multidimensional array
def convertArr(l, n):
    ''' Convert a 1D array to 2D array. '''
    return [l[i:i+n] for i in range(0, len(l), n)]

# Multiply each element in the 4 length arrays eg. [a, b, c, d] x [e, f, g, h] = [axe, bxf, cxg, dxh]
def mulV4(v1, v2):
  ''' Multiply two 1D 4 Length Arrays.'''
  return [v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2], v1[3] * v2[3]]

# Converts a multidimensional array, in this case 2D, to a regular array, 1D
def convert2Dto1D(mat):
    ''' Convert a 2D array to 1D array. '''
    size = len(mat)
    output = [0.0 for i in range(size * size)]

    for i in range(size):
        for j in range(size):
            output[i * size + j] = mat[i][j]
    return output

# Conver matrix 4x4 to 3x3
def convertM4to3(matrix):
    ''' Convert a 4x4 Matrix to 3x3. '''
    temp = [[0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]]
            
    for i in range(3):
        for j in range(3):
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
    coordsN = vec.normalize(coords)
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