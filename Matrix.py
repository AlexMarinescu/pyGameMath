import math
import sys
import Vector as vec

# NxN Matrix Class
class Matrix(object):
    
    def __init__(self, size, array):
        # This will be needed when returning a same type array.
        self.fullsize = size

        # Calculate once to use in loops.
        square_of_size = int(math.sqrt(size))

        # Setup the most often used variables.
        self.row = square_of_size
        self.col = square_of_size
        self.mat = array

        # Calculate the array length once, it might be needed in certain places.
        self.array_len = len(array)

        # Check to see if the array is the same as the size specified.
        if self.array_len != square_of_size:
            print "The array doesn't match the size specified."
            sys.exit(0)

        # Create an output array on init to save time(dynamic).
        self.out = []
        for x in xrange(self.row):
            self.out.append([])
            for y in xrange(self.col):
                self.out[x].append(0.0)

    # Slow; use the class.mat4x4[i][j] way.
    def __getitem__(self, key):
        return self.mat[key[0]][key[1]]

    # Slow; use the class.mat4x4[i][j] way.
    def __setitem__(self, key, val):
        self.mat[key[0]][key[1]] = val

    def __add__(self, input):
        for i in xrange(self.row):
            for j in xrange(self.col):
                self.out[i][j] = self.mat[i][j] + input.mat[i][j]

    def __sub__(self, input):
        for i in xrange(self.row):
            for j in xrange(self.col):
                self.out[i][j] = self.mat[i][j] - input.mat[i][j]
        
    def __mul__(self, input):     
        for i in xrange(self.row):
            for j in xrange(self.col):
                for k in xrange(self.row):
                    self.out[j][i] = self.out[j][i] + self.mat[i][k] * input.mat[k][j]
                                        
        return Matrix(self.fullsize, self.out)

    # Common, identical functions
    def sub(self, array):
        self.mat = array

    def identity(self):
        for i in xrange(self.row):
            self.mat[i][i] = 1.0

    def zero(self):
        for i in xrange(self.row):
            for j in xrange(self.row):
                self.mat[i][j] = 0.0
                
    def transpose(self):
        for i in xrange(self.row):
            for j in xrange(i+1,self.col):
                self.mat[i][j] = self.mat[j][i]

    def scale(self, vector):
        for i,j in enumerate(vector):
            self.mat[i][i] = j
            
    # This is an implementation of the Gauss Jordan Elimination to find a nxn matrix inverse.           
    def inverse(self):
        # Create a backup of the original matrix because by the end the original will be converted to identity.
        matrix_backup = self.mat
        # Create a dynamic container that can hold the inverse.
        # Basically create an augumented matrix. However, they are not in the same array.
        # The augumented matrix needs to be split at the end anyway so there is no point in joining them.
        container = []
        for x in xrange(self.row):
            container.append([])
            for y in xrange(self.row):
                container[x].append(0.0)
        # Identity matrix.
        for i in xrange(self.row):
            container[i][i] = 1.0
        # Find the inverse.
        for i in xrange(self.row):
            diag_number = float(self.mat[i][i])
            for x in xrange(self.row):
                # Divide row 1 by first number, row 2 by second number and so on... 
                # The number happens to be on the diagonal when its incremented.
                self.mat[i][x] /= diag_number
                container[i][x] /= diag_number
            # Each time a row division happens. The above stuff; the following operations needs to occur.
            for j in xrange(self.row):
                # Make sure we are not doing calculations on the row we just divided.
                if j != i:
                    # Do the same thing with the number as the first step. 
                    # However, this time make sure it is the updated one and it will no longer be on the diagonal, but on a column
                    # because its location is also given by the row number that is being calculated. 
                    # So its no longer [n][n] where n represents the number location on the diagonal, its [row number][n].
                    number = float(self.mat[j][i])
                    for y in xrange(self.row):
                        # Here the number multiplied by the numbers of the divided row is substracted from the original matrix.
                        # Example row1[i] - forth number on the row1 * row4[i]
                        # Where row1 is the one being calculated, row4 was the divided one.
                        self.mat[j][y] = self.mat[j][y] - number * self.mat[i][y]
                        container[j][y] = container[j][y] - number * container[i][y]
        # Restore the original matrix.
        self.mat = matrix_backup
        # Return a new matrix
        return Matrix(self.fullsize, container)

    def output(self):
        print "Matrix",self.row,"x",self.col,":"
        print self.mat

# Matrix 2x2 specific functions
def mat2x2ShearX(x):
    container = [[1.0,  x],
                 [0.0,1.0]]
    return Matrix(4,container)

def mat2x2ShearY(y):
    container = [[1.0,  y],
                 [  y,1.0]]
    return Matrix(4,container)

def mat2x2Rotate(angle, direction):
    container = [[1.0, 0.0],
                 [0.0, 1.0]]
    if direction == "R":
        container[0] = cos(angle)
        container[1] = sin(angle)
        container[2] =-sin(angle)
        container[3] = cos(angle)
    elif direction == "L":
        container[0] = cos(angle)
        container[1] =-sin(angle)
        container[2] = sin(angle)
        container[3] = cos(angle)
    else:
        print "Direction undefined."
    return Matrix(4, container)

    
# Matrix 4x4 specific functions
def perspective(fov, aspect, znear, zfar):
    y = tan(fov * (3.14159265358979323846 / 360))
    x = y * aspect

    if (x != 0) and (y != 0) and (zfar != znear):
        a = 1.0 / x
        b = 1.0 / y
        c = (zfar)/(znear-zfar)
        d = (zfar*znear)/(znear-zfar)

        container = [[  a, 0.0, 0.0, 0.0],
                     [0.0,   b, 0.0, 0.0],
                     [0.0, 0.0,   c,   d],
                     [0.0, 0.0, 0.0, 1.0]]
    return Matrix(16, container)

def lookAt(cam_loc, direction, head_pos):
    atNew = direction - cam_loc
    atNew.normalize()

    headNew = head_pos

    xaxis = vec.cross(atNew, headNew)
    xaxis.normalize()

    upNew = vec.cross(xaxis, atNew)

    container = [[xaxis.vec[0], upNew.vec[0], atNew.vec[0], cam_loc.vec[0]],
                 [xaxis.vec[1], upNew.vec[1], atNew.vec[1], cam_loc.vec[1]],
                 [xaxis.vec[2], upNew.vec[2], atNew.vec[2], cam_loc.vec[2]],
                 [         0.0,          0.0,          0.0,             1.0]]
    return Matrix(16, container)

# 3x3 and 4x4 Matrix functions
def translate(vector, matType):
    if matType == "4x4":
        container = [[1.0,0.0,0.0,vector.vec[0]],
                     [0.0,1.0,0.0,vector.vec[1]],
                     [0.0,0.0,1.0,vector.vec[2]],
                     [0.0,0.0,0.0,      1.0    ]]
        return Matrix(16, container)
    elif matType == "3x3":
        container = [[1.0,0.0,vector.vec[0]]
                     [0.0,1.0,vector.vec[1]]
                     [0.0,0.0,vector.vec[2]]]
        return Matrix(9, container)
    else:
        print "Not defined yet."
        
def roateX(theta, matType):
    c = cos(theta)
    s = sin(theta)   
    if matType == "4x4":
        container = [[1.0, 0.0, 0.0, 0.0],
                     [0.0,   c,  -s, 0.0],
                     [0.0,   s,   c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    elif matType == "3x3":
        container = [[1.0, 0.0, 0.0],
                     [0.0,   c,  -s],
                     [0.0,   s,   c]]
        return Matrix(9, container)
    else:
        print "Not defined yet."
        
def roateY(theta, matType):
    c = cos(theta)
    s = sin(theta)   
    if matType == "4x4":
        container = [[  c, 0.0,   s, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [ -s, 0.0,   c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    elif matType == "3x3":
        container = [[  c, 0.0,   s],
                     [0.0, 1.0, 0.0],
                     [ -s, 0.0,   c]]
        return Matrix(9, container)
    else:
        print "Not defined yet."
        
def roateZ(theta, matType):
    c = cos(theta)
    s = sin(theta)   
    if matType == "4x4":
        container = [[  c,  -s, 0.0, 0.0],
                     [  s,   c, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    elif matType == "3x3":
        container = [[  c,  -s, 0.0],
                     [  s,   c, 0.0],
                     [0.0, 0.0, 1.0]]
        return Matrix(9, container)
    else:
        print "Not defined yet."

# Roate around an arbitary axis.
def rotate(axis, theta, matType):
    c = cos(theta)
    s = sin(theta)
 
    OneMinusCos = (1.0 - c)
    axis.normalize()
    
    if matType == "4x4":
        container = [[((axis.vec[0] * axis.vec[0]) * OneMinusCos) + c,                 ((axis.vec[0] * axis.vec[1]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[0] * axis.vec[2]) * OneMinusCos) + (axis.vec[1] * s), 0.0],
                     [((axis.vec[0] * axis.vec[1]) * OneMinusCos) + (axis.vec[2] * s), ((axis.vec[1] * axis.vec[1]) * OneMinusCos) + c,                 ((axis.vec[1] * axis.vec[2]) * OneMinusCos) - (axis.vec[0] * s), 0.0],
                     [((axis.vec[0] * axis.vec[2]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[1] * axis.vec[2]) * OneMinusCos) + (axis.vec[0] * s), ((axis.vec[2] * axis.vec[2]) * OneMinusCos) + c, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    elif matType == "3x3":
        container = [[((axis.vec[0] * axis.vec[0]) * OneMinusCos) + c,                 ((axis.vec[0] * axis.vec[1]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[0] * axis.vec[2]) * OneMinusCos) + (axis.vec[1] * s)],
                     [((axis.vec[0] * axis.vec[1]) * OneMinusCos) + (axis.vec[2] * s), ((axis.vec[1] * axis.vec[1]) * OneMinusCos) + c,                 ((axis.vec[1] * axis.vec[2]) * OneMinusCos) - (axis.vec[0] * s)],
                     [((axis.vec[0] * axis.vec[2]) * OneMinusCos) - (axis.vec[2] * s), ((axis.vec[1] * axis.vec[2]) * OneMinusCos) + (axis.vec[0] * s), ((axis.vec[2] * axis.vec[2]) * OneMinusCos) + c]]
        return Matrix(9, container)
    else:
        print "Not defined yet."

def shearXY(x,y,matType):
    if matType == "3x3":
        container = [[1.0, 0.0,  x],
                     [0.0, 1.0,  y],
                     [0.0, 0.0, 1.0]]
        return Matrix(9, container)
    elif matType == "4x4":
        container = [[1.0, 0.0, 0.0,   x],
                     [0.0, 1.0, 0.0,   y],
                     [0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    else:
        print "Not defined yet."

def shearYZ(x,y,matType):
    if matType == "3x3":
        container = [[1.0, 0.0, 0.0],
                     [  x, 1.0, 0.0],
                     [  y, 0.0, 1.0]]
        return Matrix(9, container)
    elif matType == "4x4":
        container = [[1.0, 0.0, 0.0, 0.0],
                     [  x, 1.0, 0.0, 0.0],
                     [  y, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    else:
        print "Not defined yet."
        
def shearXZ(x,y,matType):
    if matType == "3x3":
        container = [[1.0,   x, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0,   y, 1.0]]
        return Matrix(9, container)
    elif matType == "4x4":
        container = [[1.0,   x, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0,   y, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0]]
        return Matrix(16, container)
    else:
        print "Not defined yet."
        
# Matrix Stack Class
class MatrixStack(object):

    def __init__(self,matrix_list):
        self.stack = matrix_list

    def __getitem__(self, key):
        return self.stack[key]

    def __setitem__(self, key, value):
        self.stack[key] = value

    def push(self, matrix):
        self.stack.append(matrix)

    def assign(self, location, matrix):
        self.stack[location] = matrix

    def pop(self,*args):
        self.stack.pop(args[0])

    def output(self):
        print self.stack