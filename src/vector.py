import math

try:
	range = xrange
except:
	pass

def Vector(n):
	''' Create an N vector or N 1D array. '''
	return [0.0 for i in range(n)]

def add(vecA, vecB):
	''' Add two vectors. '''
	out = Vector(len(vecA))
	for i in range(out):
		out[i] = vecA[i] + vecB[i]
	return out

def sub(vecA, vecB):
	''' Substract two vectors. '''
	out = Vector(len(vecA))
	for i in range(out):
		out[i] = vecA[i] - vecB[i]
	return out

def mul(scalar, vecA):
	''' Multiply a vector and a scalar. Same as scaling. '''
	out = Vector(len(vecA))
	for i in range(out):
		out[i] = vecA[i] * scalar
	return out

def div(scalar, vecA):
	''' Divide a vector by a scalar. Same as scaling. '''
	out = Vector(len(vecA))
	for i in range(out):
		out[i] = vecA[i] / scalar
	return out

def magnitude(vec):
	''' Return the magnitude of a vector. '''
	temp = 0
	for i in range(len(vec)):
		temp += vec[i] * vec[i]
	return math.sqrt(temp)

def invert(vec):
	''' Negate the vector. '''
	for i in range(len(vec)):
		vec[i] = -vec[i]
	return vec

def lerp(vecA, vecB, time):
	''' Linear interpolation between two vectors. '''
	return add(mul(vecA, time), mul(vecB, (1-time)))

def normalize(vecA):
	''' Normalize a vector. '''
	length = magnitude(vecA)
	if length != 0:
		for i in range(len(vecA)):
			self.vecA[i] /= length
	return vecA

def dot(vecA, vecB):
	''' Dot product. '''
	result = 0
	for x in range(len(vecA)):
		result += vecA[i] * vecB[i]
	return result

def cross(vecA, vecB):
	''' Cross product. '''
	vecC = Vector(3)
	vecC[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1]
	vecC[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2]
	vecC[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0]
	return vecC

def reflect(incidentVec, Norm):
	''' Reflect a vector. '''
	return sub(incidentVec, mul(2.0 * dot(incidentVec, Norm), normal))

def refract(IOR, incidentVec, Norm):
	''' Refract a vector. '''
	dotNI = dot(normal, incidentVec)
	k = 1.0 - IOR * IOR * IOR * (1.0 - dotNI * dotNI)

	if k < 0.0:
		return Vector(len(Norm))
	else:
		scalar = IOR * DOTNI + math.sqrt(k)
		return sub(mul(IOR,incidentVec), mul(scalar,normal))

