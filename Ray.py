from Vector import*
from Matrix import*

# Ray Class
class Ray(object):
    
    def __init__(self, startVector, dirVector):
        self.start = startVector
        self.dir = dirVector
        self.dir.normalize()
        
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
        self.dir.normalize()
        
    def output(self):
        print "Ray:"
        print "Start:"
        self.start.output()
        print "Dir:"
        self.dir.output()