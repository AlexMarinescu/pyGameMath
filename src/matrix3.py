import math
import src.library.math.matrix as m


def scale(val):
	''' Scale the matrix by a value.'''
	scale = [[val, 0.0, 0.0],
			     [0.0, val, 0.0],
			     [0.0, 0.0, val]]
	return scale

def translate(vector):
	''' Translate by a vector.'''
	translate = [[1.0, 0.0, vector[0]],
			         [0.0, 1.0, vector[1]],
			         [0.0, 0.0, vector[2]]]
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

  container = [[c + x2 * OneMinusCos, ((axis.vec[1] * nAxis[0]) * OneMinusCos) + (nAxis[2] * s), ((nAxis[2] * nAxis[0]) * OneMinusCos) - (nAxis[1] * s)],
               [((nAxis[0] * nAxis[1]) * OneMinusCos) - (nAxis[2] * s), c + y2 * OneMinusCos, ((nAxis[2] * nAxis[1]) * OneMinusCos) + (nAxis[0] * s)],
               [((nAxis[0] * nAxis[2]) * OneMinusCos) + (nAxis[1] * s), ((nAxis[1] * nAxis[2]) * OneMinusCos) - (nAxis[0] * s), c + z2 * OneMinusCos]]
  return container

def shearXY(x, y):
  ''' Shear on XY. '''
  out = Matrix(3)
  out[0][2] = x
  out[1][2] = y
  return out

def shearYZ(y, z):
  ''' Shear on YZ. '''
  out = Matrix(3)
  out[1][0] = y
  out[2][0] = z
  return out

def shearXZ(x, z):
  ''' Shear on XZ. '''
  out = Matrix(3)
  out[0][1] = x
  out[2][1] = z
  return out