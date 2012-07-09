import math
import sys

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

    def output(self):
        print "Matrix",self.row,"x",self.col,":"
        print self.mat

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
