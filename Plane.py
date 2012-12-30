from Math import*
from Vector import*
from Point import*

# I know that I'm missing some shit in this "plane" class but bleh, for now the current stuff will do. :3
class Plane(object):
    # A, B, C make up the normal of the plane
    def __init__(self, *args):
        # Define using 3 vectors.
        if isinstance(args[0], Vector) and isinstance(args[1], Vector) and isinstance(args[2], Vector):
            v1 = args[1] - args[0]
            v2 = args[2] - args[0]
            n = cross(v1, v2)
            n.normalize()
            self.A = n.vec[0]
            self.B = n.vec[1]
            self.C = n.vec[2]
            self.D = dot(n.scale(-1.0), args[0])
        # Define using a point and a normal.
        elif isintance(args[0], Point) and isinstance(args[1], Vector):
            self.A = args[1].vec[0]
            self.B = args[1].vec[1]
            self.C = args[1].vec[2]
            self.D = -dot(args[1], args[0].point)
        # Define by input A,B,C,D manually.
        elif isinstance(args[0], float) and isinstance(args[1], float) and isinstance(args[2], float):
            self.A = A
            self.B = B
            self.C = C
            self.D = D
        else:
            print "FFFUUU.... Cannot define plane."
            
    def duplicate(self):
        return Plane(self.A, self.B, self.C, self.D)
            
    # Dot product for the plane, same as vector 4D dot product.
    def dot(self, vector):
        return self.A * vector.vec[0] + self.B * vector.vec[1] + self.C * vector.vec[2] + self.D * vector.vec[3]
        
    def normalize(self):
        # Normalize the normal. -L-
        temp_vec = Vector([self.A, self.B, self.C])
        temp_vec.normalize()
        
        # Calculate D so everything is balanced. Ax + By + Cz = -D
        # Divide it by the magnitude of the normal.
        len = temp_vec.magnitude()
        if len != 0:
            new_D = self.D / len
        else:
            "Dividing by zero... LOL"
        
        # Return the normalized plane
        return Plane(temp_vec.vec[0], temp_vec.vec[1], temp_vec.vec[2], new_D)
		
    def bestFitNormal(self, vecList):
        output = Vector([0.0,0.0,0.0])
        
        for i,j in enumerate(vecList):
            output.vec[0] += (vecList[i].vec[2] + vecList[i+1].vec[2]) * (vecList[i].vec[1] - vecList[i+1].vec[1])
            output.vec[1] += (vecList[i].vec[0] + vecList[i+1].vec[0]) * (vecList[i].vec[2] - vecList[i+1].vec[2])
            output.vec[3] += (vecList[i].vec[1] + vecList[i+1].vec[1]) * (vecList[i].vec[0] - vecList[i+1].vec[0])
            
        output.normalize()
        return output
        
    def bestFitD(self, vecList, bestFitNormal):
        val = 0.0
        
        for i,j in enumerate(vecList):
            val = dot(j, bestFitNormal)
        
        return val / len(vecList)
        
    def distanceToPoint(thePoint):
        return dot(thePoint, Vector([self.A, self.B, self.C])) - self.D
        
    def point_location(self, point):
        # If s > 0 then the point is on the same side as the normal. (front)
        # If s < 0 then the point is on the opposide side of the normal. (back)
        # If s = 0 then the point lies on the plane.
        s = self.A * point.vec[0] + self.B * point.vec[1] + self.C * point.vec[2] + self.D
        
        if s > 0:
            point.location = "FRONT"
        elif s < 0:
            potin.location = "BEHIND"
        elif s == 0:
            point.location = "ON_PLANE"
        else
            print "No clue where the point is. :O"
            
    def output(self):
        print "Plane: ", "A: ", self.A, " B: ", self.B, " C: ", self.C, " D: ", self.D