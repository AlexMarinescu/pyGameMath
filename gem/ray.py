import gem.vector as vec

# Ray Class
class Ray(object):
    def __init__(self, startVector, dirVector):
        ''' Initiated a ray with the start vector and direction vector. '''
        self.start = startVector
        self.dir = dirVector
        self.distance = self.dir.magnitude()
        self.dir.i_normalize()
        # The end is used when intersections are added so
        # we can know where the ray stops.
        self.end = vec.Vector(3)

    def duplicate(self):
        ''' Make a duplicate of the ray. '''
        return Ray(self.start, self.dir)

    def roateUsingMatrix(self, matrix):
        ''' Rotate the ray using a matrix. '''
        self.start = matrix * self.start
        self.dir = matrix * self.dir
        self.dir.i_normalize()

    def rotateUsingQuaternion(self, quat1):
        ''' Rotate the ray using a quaternion. '''
        self.start = quat1 * self.start
        self.dir = quat1 * self.dir
        self.dir.i_normalize()

    def translate(self, matrix):
        ''' Translate the ray using a matrix. '''
        self.dir = matrix * self.dir
        self.distance = self.dir.magnitude()
        self.dir.i_normalize()

    def output(self):
        ''' Show information regarding the ray's behaviour. '''
        print ("Ray:")
        print ("Start:")
        print (self.start.vector)
        print ("Dir:")
        print (self.dir.vector)
