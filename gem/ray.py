import gem.vector as vec
import gem.quaternion as quat

# Ray Class
class Ray(object):
    def __init__(self, startVector, dirVector):
        ''' Initiated a ray with the start vector and direction vector. '''
        self.start = startVector
        self.dir = dirVector
        self.distance = self.dir.magnitude()
        vec._require_nonzero(self.dir.size, self.dir.vector)
        self.dir.i_normalize()
        # The end is used when intersections are added so
        # we can know where the ray stops.
        self.end = vec.Vector(3)

    def duplicate(self):
        """Copy stored geometry and intersection state without normalization."""
        result = object.__new__(Ray)
        result.start = self.start.clone()
        result.dir = self.dir.clone()
        result.end = self.end.clone()
        result.distance = self.distance
        return result

    def roateUsingMatrix(self, matrix):
        ''' Rotate the ray using a matrix. '''
        self.start = matrix * self.start
        self.dir = matrix * self.dir
        vec._require_nonzero(self.dir.size, self.dir.vector)
        self.dir.i_normalize()

    def rotateUsingQuaternion(self, quat1):
        """Rotate about the coordinate origin using a unit Quaternion.

        Distance and stored intersection state are preserved.
        """
        self.start = quat.quat_rotate_vector(quat1, self.start)
        self.dir = quat.quat_rotate_vector(quat1, self.dir)
        vec._require_nonzero(self.dir.size, self.dir.vector)
        self.dir.i_normalize()

    def translate(self, matrix):
        """Apply a pure translation, preserving distance and intersection state.

        Vector3/Matrix4 positions and directions receive local w=1/w=0.
        No perspective division or general operator promotion is performed.
        """
        if matrix.size == 4 and self.start.size == 3 and self.dir.size == 3:
            position = matrix * vec.Vector(4, self.start.vector + [1.0])
            direction = matrix * vec.Vector(4, self.dir.vector + [0.0])
            self.start = vec.Vector(3, position.vector[:3])
            self.dir = vec.Vector(3, direction.vector[:3])
        else:
            self.start = matrix * self.start
            self.dir = matrix * self.dir

    def output(self):
        ''' Show information regarding the ray's behaviour. '''
        print ("Ray:")
        print ("Start:")
        print (self.start.vector)
        print ("Dir:")
        print (self.dir.vector)
