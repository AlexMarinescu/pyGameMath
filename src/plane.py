from src.pycompat import *
import src.vector as vector

def Plane(*args):
	''' Define a plane using 3 vectors, return a 4 item 1D array. '''
	# Define using 3 vectors
	v1 = args[1] - args[0]
	v2 = args[2] - args[0]
	n = vector.cross(v1, v2)
	n = vector.normalize(n)
	return [n[0], n[1], n[2], vector.dot(vector.mul(n, -1.0), args[0])]

def dot(plane, vector):
	''' Return the dot product between a plane and 4D vector. '''
	return plane[0] * vector[0] + plane[1] * vector[1] + plane[2] * vector[2] + plane[3] * vector[3]

def normalize(plane):
	''' Return the normalized plane.'''
	vec = [plane[0], plane[1], plane[2]]
	vecN = vector.normalize(vec)

	length = vector.magnitude(vecN)

	if length != 0:
		return [vecN[0], vecN[1], vecN[2], plane[3] / length]
	else:
		print("Plane fail to normalize due to zero division.")
		return [0.0, 0.0, 0.0, 0.0]

def bestFitNormal(vecList):
	output = [0.0, 0.0, 0.0]
	for i in range(len(vecList)):
		output[0] += (vecList[i][2] + vecList[i + 1][2]) * (vecList[i][1] - vecList[i + 1][1])
		output[1] += (vecList[i][0] + vecList[i + 1][0]) * (vecList[i][2] - vecList[i + 1][2])
		output[2] += (vecList[i][1] + vecList[i + 1][1]) * (vecList[i][0] - vecList[i + 1][0])
	return vector.normalize(output)

def bestFitD(vecList, bestFitNormal):
	val = 0.0
	for vec in vecList:
		val += vector.dot(vec, bestFitNormal)
	return val / len(vecList)

def point_location(plane, point):
	''' Returns the location of the point. '''
	# If s > 0 then the point is on the same side as the normal. (front)
	# If s < 0 then the point is on the opposide side of the normal. (back)
	# If s = 0 then the point lies on the plane.
	s = plane[0] * point[0] + plane[1] * point[1] + plane[2] * point[2] + plane[3]

	if s > 0:
		return 1
	elif s < 0:
		return -1
	elif s == 0:
		return 0
	else:
		print("Not a clue where the point is. O_o")