import math
import src.library.math.matrix as m
from src.library.math.vector import Vector, normalize, cross, dot, sub, add, neg, div, mul
from src.library.math.constants import PI
from src.library.math.common import convert2Dto1D, mulV4

def orthographic(left, right, bottom, top, near, far):
    ''' Orthographic projection matrix.'''
    a = 2.0 / (right - left)
    b = 2.0 / (top - bottom)
    c = -2.0 / (far - near)

    tx = -(right + left) / (right - left)
    ty = -(top + bottom) / (top - bottom)
    tz = -(far + near) / (far - near)

    out = [[a, 0.0, 0.0, 0.0],
           [0.0, b, 0.0, 0.0],
           [0.0, 0.0, c, 0.0],
           [tx, ty, tz, 1.0]]
    return out

def perspective(fov, aspect, znear, zfar):
    ''' Perspective projection matrix 4x4. FOVY'''
    rad = math.radians(fov)

    tanHalfFovy = math.tan(rad / 2.0)

    a = 1.0 / (aspect * tanHalfFovy)
    b = 1.0 / (tanHalfFovy)
    c = - (zfar + znear) / (zfar - znear)
    d = - (2.0 * zfar * znear) / (zfar - znear)
  
    out = [[  a, 0.0, 0.0, 0.0],
           [0.0,   b, 0.0, 0.0],
           [0.0, 0.0,   c,-1.0],
           [0.0, 0.0,   d, 0.0]]
    return out


def perspectiveX(fov, aspect, znear, zfar):
    ''' Perspective projection matrix 4x4. FOVX'''
    f = znear * math.tan((fov * PI / 360.0))

    xmax = f
    xmin = -f
    
    ymin = xmin / aspect
    ymax = xmax / aspect
    
    a = (2.0 * znear) / (xmax - xmin)
    b = (2.0 * znear) / (ymax - ymin)   
    c = -(zfar + znear) / (zfar - znear)
    d = -(2.0 * zfar * znear) / (zfar - znear)
    e = (xmax + xmin) / (xmax - xmin)
    f = (ymax + ymin) / (ymax - ymin)
    
    out = [[  a, 0.0, 0.0, 0.0],
           [0.0,   b, 0.0, 0.0],
           [  e,   f,   c,-1.0],
           [0.0, 0.0,   d, 0.0]]
    return out

def lookAt(eye, center, up):
    ''' Matrix 4x4 lookAt function.'''
    f = normalize(sub(center, eye))
    u = normalize(up)
    s = normalize(cross(f, u))
    u = cross(s, f)

    output = [[s[0], u[0], -f[0], 0.0],
              [s[1], u[1], -f[1], 0.0],
              [s[2], u[2], -f[2], 0.0],
              [-dot(s, eye), -dot(u, eye), dot(f, eye), 1.0]]
    return output

def scale(val):
    ''' Scale the matrix by a value.'''
    scale = [[val[0], 0.0, 0.0, 0.0],
             [0.0, val[1], 0.0, 0.0],
             [0.0, 0.0, val[2], 0.0],
             [0.0, 0.0, 0.0, 1.0]]
    return scale

def translate(vector):
    ''' Translate by a vector.'''
    translate = [[1.0, 0.0, 0.0, 0.0],
                 [0.0, 1.0, 0.0, 0.0],
                 [0.0, 0.0, 1.0, 0.0],
                 [vector[0], vector[1], vector[2], 1.0]]
    return translate

def rotate(axis, theta):
    ''' Rotate around an axis.'''
    c = math.cos(theta*(PI/180))
    s = math.sin(theta*(PI/180))

    OneMinusCos = (1.0 - c)

    nAxis = normalize(axis)

    x2 = nAxis[0] * nAxis[0]
    y2 = nAxis[1] * nAxis[1]
    z2 = nAxis[2] * nAxis[2]

    container = [[c + x2 * OneMinusCos, ((nAxis[1] * nAxis[0]) * OneMinusCos) + (nAxis[2] * s), ((nAxis[2] * nAxis[0]) * OneMinusCos) - (nAxis[1] * s), 0.0],
                 [((nAxis[0] * nAxis[1]) * OneMinusCos) - (nAxis[2] * s), c + y2 * OneMinusCos, ((nAxis[2] * nAxis[1]) * OneMinusCos) + (nAxis[0] * s), 0.0],
                 [((nAxis[0] * nAxis[2]) * OneMinusCos) + (nAxis[1] * s), ((nAxis[1] * nAxis[2]) * OneMinusCos) - (nAxis[0] * s), c + z2 * OneMinusCos, 0.0],
                 [0.0, 0.0, 0.0, 1.0]]
    return container

def project(obj, model, proj, viewport):
    ''' The most hacked together project code in the world. :| It works tho. :3 '''
    projM = convert2Dto1D(proj)
    modelM = convert2Dto1D(model)

    T = m.Matrix(4)
    for r in range(4):
        for c in range(4):
            T[r][c] = 0.0
            for i in range(4):
                T[r][c] += projM[r + i * 4] * modelM[i + c *4]

    result = Vector(4)

    for r in range(4):
        result[r] = dot(T[r], obj)

    rhw = 1.0 / result[3]

    return [(1 + result[0] * rhw) * viewport[2] / 2.0 + viewport[0],
            (1 + result[1] * rhw) * viewport[3] / 2.0 + viewport[1],
            (result[2] * rhw) * (1 - 0) + 0, rhw]

def unproject(winx, winy, winz, modelview, projection, viewport):
    ''' Unproject a point from the screen and return the object coordinates. '''
    m = Matrix(4)
    IN = [0.0, 0.0, 0.0, 0.0]
    objCoord = [0.0, 0.0, 0.0]

    A = m.mulM(projection, modelview)
    m = m.inverse(A)

    IN[0] = (winx - viewport[0]) / viewport[2] * 2.0 - 1.0
    IN[1] = (winy - viewport[1]) / viewport[3] * 2.0 - 1.0
    IN[2] = 2.0 * winz - 1.0
    IN[3] = 1.0

    OUT = m.mulV(m, IN)
    if(OUT[3] == 0.0):
        return [0.0, 0.0, 0.0]

    OUT[3] = 1.0 / OUT[3]
    objCoord[0] = out[0] * out[3]
    objCoord[1] = out[1] * out[3]
    objCoord[2] = out[2] * out[3]
    return objCoord

def shearXY(x, y):
    ''' Shear on XY. '''
    out = Matrix(4)
    out[0][3] = x
    out[1][3] = y
    return out

def shearYZ(y, z):
    ''' Shear on YZ. '''
    out = Matrix(4)
    out[1][0] = y
    out[2][0] = z
    return out

def shearXZ(x, z):
    ''' Shear on XZ. '''
    out = Matrix(4)
    out[0][1] = x
    out[2][1] = z
    return out

def inverse(mat):
    ''' Inverse of a 4x4 matrix.'''
    Coef00 = mat[2][2] * mat[3][3] - mat[3][2] * mat[2][3]
    Coef02 = mat[1][2] * mat[3][3] - mat[3][2] * mat[1][3]
    Coef03 = mat[1][2] * mat[2][3] - mat[2][2] * mat[1][3]

    Coef04 = mat[2][1] * mat[3][3] - mat[3][1] * mat[2][3]
    Coef06 = mat[1][1] * mat[3][3] - mat[3][1] * mat[1][3]
    Coef07 = mat[1][1] * mat[2][3] - mat[2][1] * mat[1][3]

    Coef08 = mat[2][1] * mat[3][2] - mat[3][1] * mat[2][2]
    Coef10 = mat[1][1] * mat[3][2] - mat[3][1] * mat[1][2]
    Coef11 = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]

    Coef12 = mat[2][0] * mat[3][3] - mat[3][0] * mat[2][3]
    Coef14 = mat[1][0] * mat[3][3] - mat[3][0] * mat[1][3]
    Coef15 = mat[1][0] * mat[2][3] - mat[2][0] * mat[1][3]

    Coef16 = mat[2][0] * mat[3][2] - mat[3][0] * mat[2][2]
    Coef18 = mat[1][0] * mat[3][2] - mat[3][0] * mat[1][2]
    Coef19 = mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]

    Coef20 = mat[2][0] * mat[3][1] - mat[3][0] * mat[2][1]
    Coef22 = mat[1][0] * mat[3][1] - mat[3][0] * mat[1][1]
    Coef23 = mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]

    SignA = [+1, -1, +1, -1]
    SignB = [-1, +1, -1, +1]

    Fac0 = [Coef00, Coef00, Coef02, Coef03]
    Fac1 = [Coef04, Coef04, Coef06, Coef07]
    Fac2 = [Coef08, Coef08, Coef10, Coef11]
    Fac3 = [Coef12, Coef12, Coef14, Coef15]
    Fac4 = [Coef16, Coef16, Coef18, Coef19]
    Fac5 = [Coef20, Coef20, Coef22, Coef23]

    Vec0 = [mat[1][0], mat[0][0], mat[0][0], mat[0][0]]
    Vec1 = [mat[1][1], mat[0][1], mat[0][1], mat[0][1]]
    Vec2 = [mat[1][2], mat[0][2], mat[0][2], mat[0][2]]
    Vec3 = [mat[1][3], mat[0][3], mat[0][3], mat[0][3]]

    Inv0 = mulV4(SignA, add( sub( mulV4( Vec1, Fac0), mulV4( Vec2, Fac1)), mulV4( Vec3, Fac2)))
    Inv1 = mulV4(SignB, add( sub( mulV4( Vec0, Fac0), mulV4( Vec2, Fac3)), mulV4( Vec3, Fac4)))
    Inv2 = mulV4(SignA, add( sub( mulV4( Vec0, Fac1), mulV4( Vec1, Fac3)), mulV4( Vec3, Fac5)))
    Inv3 = mulV4(SignB, add( sub( mulV4( Vec0, Fac2), mulV4( Vec1, Fac4)), mulV4( Vec2, Fac5)))

    Inverse = [Inv0, Inv1, Inv2, Inv3]

    Row0 = [Inverse[0][0], Inverse[1][0], Inverse[2][0], Inverse[3][0]]

    Determinant = dot(mat[0], Row0)

    return m.div(Inverse, Determinant)
