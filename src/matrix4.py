import math
import src.library.math.matrix as m
from src.library.math.vector import normalize, cross, dot, sub, neg, div
from src.library.math.constants import PI

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
    f = znear * math.tan((fov * PI / 360.0))
     
    # Flip signs to change dir :P
    left = -f * aspect
    right = f * aspect
    
    bottom = -f
    top = f
    
    a = (2.0 * znear) / (right - left)
    b = (2.0 * znear) / (top - bottom)
    c = -(-zfar - znear) / (znear - zfar)
    d = -(-(2.0 * zfar) * znear) / (znear - zfar)
    f = (right + left) / (right - left)
    g = (top + bottom) / (top - bottom)

    out = [[  a, 0.0, 0.0, 0.0],
           [0.0,   b, 0.0, 0.0],
           [  f,   g,   c,-1.0],
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
	translate = [[1.0, 0.0, 0.0, vector[0]],
			         [0.0, 1.0, 0.0, vector[1]],
			         [0.0, 0.0, 1.0, vector[2]],
			         [0.0, 0.0, 0.0, 1.0]]
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

def project(obj, model, proj, view):
  ''' The most hacked together project code in the world. :| '''
  matrix = m.mulM(model, proj)
  vector = m.mulV(matrix, [obj[0], obj[1], obj[2], 1.0])

  a = (obj[0] * matrix[0][3]) + (obj[1] * matrix[1][3]) + (obj[2] * matrix[2][3]) + matrix[3][3]

  if a == 0.0:
    vector = vector
  else:
    vector = div(a, vector)

  vector[0] = (vector[0] + 1) * 0.5 * view[2] + view[0]
  vector[1] = (vector[1] + 1) * 0.5 * view[3] + view[1]
  return vector

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
