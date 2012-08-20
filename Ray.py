from Vector import*
from Matrix import*

# Ray Class
class Ray(object):
    
    def __init__(self, startVector, dirVector):
        self.start = startVector
        self.dir = dirVector
        self.distance = self.dir.magnitude()
        self.dir.normalize()
        # The end is used when intersections are added so we can know where the ray stops.
        self.end = Vector([0.0,0.0,0.0])
        
    def duplicate(self):
        return Ray(self.start, self.dir)
        
    def roateUsingMatrix(self, matrix):
        self.start = matrix * self.start
        self.dir   = matrix * self.dir
        self.dir.normalize()
        
    def rotateUsingQuaternion(self, quat):
        self.start = quat * self.start
        self.dir   = quat * self.dir
        self.dir.normalize()
        
    def translate(self, matrix):
        self.dir = matrix * self.dir
        self.distance = self.dir.magnitude()
        self.dir.normalize()
        
    def output(self):
        print "Ray:"
        print "Start:"
        self.start.output()
        print "Dir:"
        self.dir.output()