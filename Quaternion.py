import math
from Matrix import*
from Vector import*

# Quaternion class
class Quaternion(object):
    def __init__(self,array):
        self.quat = array
        
    def identity(self):
        self.quat = [1.0, 0.0, 0.0, 0.0]
        self.len = len(self.quat)
        
    def magnitude(self):
        result = 0
        for x in xrange(self.len):
            result += self.quat[x] * self.quat[x]
        return math.sqrt(result)
        
    def normalize(self):
        length = self.magnitude()
        if(length != 0):
            for x in xrange(self.len):
                self.quat[x] /= length
                
    def conjugate(self):
        new_quat = self.quat
        for x in xrange(self.len):
            new_quat[x] = -new_quat[x]
        new_quat[0] = -new_quat[0]
        return Quaternion(new_quat)
    
    def dot(self, quatA):
        result = 0
        for x in xrange(self.len):
            result += self.quat[x] * quatA.quat[x]
        return result
        
    def rotateX(self, theta):
        thetaOver2 = theta * 0.5
        return Quaternion([cos(thetaOver2),sin(thetaOver2),0.0,0.0])
    
    def rotateY(self, theta):
        thetaOver2 = theta * 0.5
        return Quaternion([cos(thetaOver2),0.0,sin(thetaOver2),0.0])
        
    def rotateZ(self, theta):
        thetaOver2 = theta * 0.5
        return Quaternion([cos(thetaOver2),0.0,0.0,sin(thetaOver2)])
        
    def rotate(self, axis, theta):
        thetaOver2 = theta * 0.5
        sinThetaOver2 = sin(thetaOver2)
        axis.normalize()
        return Quaternion([cos(thetaOver2), axis.vec[1] * sinThetaOver2,
                           axis.vec[2] * sinThetaOver2, axis.vec[3] * sinThetaOver2])
    
    # Conver to a 4x4 matrix
    def toMatrix(self):
        x2 = self.quat[1] * self.quat[1]
        y2 = self.quat[2] * self.quat[2]
        z2 = self.quat[3] * self.quat[3]
        xy = self.quat[1] * self.quat[2]
        xz = self.quat[1] * self.quat[3]
        yz = self.quat[2] * self.quat[3]
        wx = self.quat[0] * self.quat[1]
        wy = self.quat[0] * self.quat[2]
        wz = self.quat[0] * self.quat[3]
        
        return Matrix(16, [[1.0 - 2.0 * (y2 * z2), 2.0 * (xy - wz), 2.0 * (xz + wy), 0.0],
                           [2.0 * (xy + wz), 1.0 - 2.0 * (x2 + z2), 2.0 * (yz - wz), 0.0],
                           [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (x2 + y2), 0.0],
                           [0.0, 0.0, 0.0, 1.0]]
    
    # Cross Product
    def __mul__(self, quat):
        return Quaternion([self.quat[0] * quat.quat[1] + self.quat[1] * quat.quat[0] + self.quat[2] * quat.quat[3] - self.quat[3] * quat.quat[2],
                           self.quat[0] * quat.quat[2] + self.quat[2] * quat.quat[0] + self.quat[3] * quat.quat[1] - self.quat[1] * quat.quat[3],
                           self.quat[0] * quat.quat[3] + self.quat[3] * quat.quat[0] + self.quat[1] * quat.quat[2] - self.quat[2] * quat.quat[1],
                           self.quat[0] * quat.quat[0] - self.quat[1] * quat.quat[1] - self.quat[2] * quat.quat[2] - self.quat[3] * quat.quat[3]])
                           
    # Vector 3D * Quaternion
    def __mul__(self, vec):
        vec.normalize()
        vecQuat = Quaternion([1.0, 0.0, 0.0, 0.0])
        resQuat = Quaternion([1.0, 0.0, 0.0, 0.0])
        
        vecQuat.quat[1] = vec.vec[0]
        vecQuat.quat[2] = vec.vec[1]
        vecQuat.quat[3] = vec.vec[2]
        vecQuat.quat[0] = 0.0
        
        resQuat = vecQuat * self.conjugate()
        resQuat = self * resQuat
        
        return Vector([resQuat.quat[1], resQuat.quat[2], resQuat.quat[3]])
    
    def __imul__(self, quat):
        return Quaternion([self.quat[0] * quat.quat[1] + self.quat[1] * quat.quat[0] + self.quat[2] * quat.quat[3] - self.quat[3] * quat.quat[2],
                           self.quat[0] * quat.quat[2] + self.quat[2] * quat.quat[0] + self.quat[3] * quat.quat[1] - self.quat[1] * quat.quat[3],
                           self.quat[0] * quat.quat[3] + self.quat[3] * quat.quat[0] + self.quat[1] * quat.quat[2] - self.quat[2] * quat.quat[1],
                           self.quat[0] * quat.quat[0] - self.quat[1] * quat.quat[1] - self.quat[2] * quat.quat[2] - self.quat[3] * quat.quat[3]])
                                   
    def __add__(self, quat):
        return Quaternion([self.quat[0] + quat.quat[0], self.quat[1] + quat.quat[1], self.quat[2] + quat.quat[2], self.quat[3] + quat.quat[3]])

    def __iadd__(self, quat):
        return Quaternion([self.quat[0] + quat.quat[0], self.quat[1] + quat.quat[1], self.quat[2] + quat.quat[2], self.quat[3] + quat.quat[3]])
 
    def __sub__(self, quat):
        return Quaternion([self.quat[0] - quat.quat[0], self.quat[1] - quat.quat[1], self.quat[2] - quat.quat[2], self.quat[3] - quat.quat[3]])
        
    def __isub__(self, quat):
        return Quaternion([self.quat[0] - quat.quat[0], self.quat[1] - quat.quat[1], self.quat[2] - quat.quat[2], self.quat[3] - quat.quat[3]])
        