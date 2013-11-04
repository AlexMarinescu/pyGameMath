import src.library.math.vector as vec
import src.library.math.matrix as mat
import src.library.math.quaternion as quat


# Ray Class
class Ray(object):
    def __init__(self, startVector, dirVector):
        ''' Initiated a ray with the start vector and direction vector. '''
        self.start = startVector
        self.dir = dirVector
        self.distance = vec.magnitude(self.dir)
        self.dir = vec.normalize(self.dir)
        # The end is used when intersections are added so
        # we can know where the ray stops.
        self.end = vec.Vector(3)

    def duplicate(self):
        ''' Make a duplicate of the ray. '''
        return Ray(self.start, self.dir)

    def roateUsingMatrix(self, matrix):
        ''' Rotate the ray using a matrix. '''
        self.start = mat.mulV(matrix, self.start)
        self.dir = mat.mulV(matrix, self.dir)
        self.dir = vec.normalize(self.dir)

    def rotateUsingQuaternion(self, quat1):
        ''' Rotate the ray using a quaternion. '''
        self.start = quat.mulV(quat1, self.start)
        self.dir = quat.mulV(quat1, self.dir)
        self.dir = vec.normalize(self.dir)

    def translate(self, matrix):
        ''' Translate the ray using a matrix. '''
        self.dir = mat.mulV(matrix, self.dir)
        self.distance = vec.magnitude(self.dir)
        self.dir = vec.normalize(self.dir)

    def output(self):
        ''' Show information regarding the ray's behaviour. '''
        print ("Ray:")
        print ("Start:")
        print (self.start)
        print ("Dir:")
        print (self.dir)
