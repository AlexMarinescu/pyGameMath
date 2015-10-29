import six.moves as sm
from gem import matrix
from gem import vector

# Some vector operations examples
vectorA = vector.Vector(3, data=[1, 2, 3])
vectorB = vector.Vector(3, data=[4, 5, 6])
vectorC = vectorA + vectorB
print("Vector C output:" , vectorC.vector)

vectorD = vectorA * 2
print("Vector A * 2:", vectorD.vector)

print("Vector A magnitude:", vectorA.magnitude())

# Some matrix operations examples
matrixA = matrix.Matrix(4, data=[[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
matrixB = matrix.Matrix(4, data=[[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
matrixC = matrixA * matrixB

print("Matrix C output:")
for i in sm.range(matrixC.size):
	print(matrixC.matrix[i])
print("End of Matrix C output")