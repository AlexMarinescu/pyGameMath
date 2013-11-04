import math

# Convert 1D to multidimensional array
def convertArr(l, n):
	return [l[i:i+n] for i in range(0, len(l), n)]

# Conver matrix 4x4 to 3x3
def convertM4to3(matrix):

	temp = [[0.0, 0.0, 0.0],
			[0.0, 0.0, 0.0],
			[0.0, 0.0, 0.0]]
			
	for i in range(3):
		for j in range(3):
			matrix.mat[i][j] = temp[i][j]

	return temp

# Sinc function
def sinc(x):
    if math.fabs(x) < 1.0e-4:
        return 1.0
    else:
        return math.sin(x) / x

# Scalar lerp
def scalarLerp(a, b, time):
    return a + time * (b - a)
