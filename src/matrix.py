from copy import deepcopy
from src.library.math.constants import PI

try:
	range = xrange
except:
	pass

def Matrix(size):
	'''Returns a 2D array/identity matrix of a specific size. Zero matrix.'''
	return [[1.0 if j == i else 0.0 for j in range(size)] for i in range(size)]

def add(matA, matB):
	'''Adds two matrices and returns a new one.'''
	size = len(matA)
	matC = [[0.0 for i in range(size)] for j in range(size)]
	for i in range(size):
		for j in range(size):
			matC[i][j] = matA[i][j] + matB[i][j]
	return matC

def sub(matA, matB):
	'''Substracts two matrices and retuns a new one.'''
	size = len(matA)
	matC = [[0.0 for i in range(size)] for j in range(size)]
	for i in range(size):
		for j in range(size):
			matC[i][j] = matA[i][j] - matB[i][j]
	return matC

def mulM(matA, matB):
	'''Multiply two matrices and return a new matrix.'''
	size = len(matA)
	matC = [[0.0 for i in range(size)] for j in range(size)]
	for i in range(size):
		for j in range(size):
			for k in range(size):
				matC[i][j] += matA[i][k] * matB[k][j]
	return matC

def mulV(mat, vec):
	'''Multiply a matrix and vectors and return a new vector.'''
	size = len(mat)
	vecS = len(vec)
	vecO = [0.0 for i in range(vecS)]
	for i in range(size):
		for j in range(size):
			vecO[i] += vec[j] * mat[i][j]
	return vecO

def transpose(mat):
	'''Transposes a matrix.'''
	size = len(mat)
	out = [[0.0 for i in range(size)] for j in range(size)]
	for i in range(size):
		for j in range(size):
			out[i][j]= mat[j][i]
	return out

def inverse(matrix):
  ''' Find the inverse of a NxN matrix.'''
  n = len(matrix)
  # Augumented matrix
  m = [[0.0 for j in range(n)] for i in range(n*2)]
  # Identity matrix
  i = Matrix(n)

  for cur_row in range(n):
    for cur_column in range(n):
      m[cur_column][cur_row] = matrix[cur_column][cur_row]

  for cur_row in range(n):
    for cur_column in range(n, n * 2):
      m[cur_column][cur_row] = i[(cur_column - n)][cur_row]

  return echelon(m, n)

def echelon(m, n):
  mm = m
  inverse = Matrix(n)
  for i in range(n-1):
    for cur_row in range(i, n):
      row_index_max = find_max_row(m, i, cur_row, n)
      for cur_column in range(n * 2):
        temp = m[cur_column][cur_row]
        mm[cur_column][cur_row] = mm[cur_column][row_index_max]
        mm[cur_column][row_index_max] = temp

    mm = reduce_echelon(mm, i, n)

    if mm[i][i] == 0:
      print("Matrix has no inverse")

  for i in range(n-1, 0, -1):
    mm = reverse_echelon(mm, i, n)

  for i in range(n):
    mm = canonical_form(mm, i, n)

  for cur_row in range(n):
    for cur_column in range(n):
      inverse[cur_column][cur_row] = mm[cur_column + n][cur_row]

  return inverse

def find_max_row(m, cur_column, start_row, end):
  index_max_row = start_row
  for cur_row in range(start_row, end):
    if m[cur_column][cur_row] > m[cur_column][index_max_row]:
      index_max_row = cur_row
  return index_max_row

def reduce_echelon(m, index, n):
  mt = m
  for cur_row in range(index, n-1):
    factor = mt[index][cur_row+1]/mt[index][index]
    for cur_column in range(index, n*2):
      temp = mt[cur_column][cur_row + 1] - (factor * mt[cur_column][index])
      mt[cur_column][cur_row + 1] = temp
  return mt

def reverse_echelon(m, index, n):
  mt = m
  for cur_row in range(index, 0, -1):
    factor = mt[index][cur_row-1]/mt[index][index]
    for cur_column in range(index, 2 * n):
      temp = mt[cur_column][cur_row - 1] - (factor * mt[cur_column][index])
      mt[cur_column][cur_row - 1] = temp
  return mt

def canonical_form(m, index, n):
  mt = m
  factor = mt[index][index]
  for cur_column in range(index, n*2):
    temp = mt[cur_column][index]/factor
    mt[cur_column][index] = temp
  return mt

def pivot(mat):
	'''Pivot the matrix.'''
	size = len(mat)
	# Create idenity matrix
	out = [[1.0 if j==i else 0.0 for j in range(size)] for i in range(size)]
	for j in range(size):
		row = max(range(j, size), key=lambda i: mat[i][j])
		if j != row:
			out[j], j[row] = out[row], out[j]
	return out

def LUdecompostion(mat):
	'''LU matrix decomposition for NxN Matrix. Returns U, L.'''
	size = len(mat)

	# Create L and U matrix
	lower = [[1.0 if j==i else 0.0 for j in range(size)] for i in range(size)]
	upper = [[1.0 if j==i else 0.0 for j in range(size)] for i in range(size)]

	# Pivot the matrix
	pivotM = pivot(mat)
	a2 = mulM(pivot, mat)
	for j in range(size):
		for i in range(j + 1):
			itemSum = sum(upper[k][j] * lower[i][k] for k in range(i))
			upper[i][j] = mat[i][j] - itemSum
		for i in range(j, size):
			itemSum = sum(upper[k][j] * lower[i][k] for k in range(j))
			if upper[j][i] == 0.0:
				upper[j][j] = 1.0
			lower[i][j] = (mat[i][j] - itemSum) / upper[j][j]
	return upper, lower

def determinant(mat):
	'''NxN Matrix determinant.'''
	U, L = LUdecompostion(mat)
	size = len(mat)
	det = 1.0
	# Multiply the diagonal
	for i in range(size):
		det *= U[i][i]
	return det * -1.0

def normalize(mat):
	'''NxN Matrix normalization.'''
	size = len(mat)
	det = determinant(mat)
	if det == 0.0:
		det = 1.0
	out = [[0.0 for i in range(size)] for j in range(size)]
	for x in range(size):
		for y in range(size):
			out[x][y] = mat[x][y] / det
	return out