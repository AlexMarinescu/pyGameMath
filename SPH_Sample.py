import math
import sys
from Vector import*

class SPHSample (object):
    
    def __init__(self, theta, phi, dir, val):
        
        # Spherical coordinates
        self.theta = theta
        self.phi = phi
        
        # Value of SH function at this point
        self.val = val
        
        if isinstance(dir, Vector):
            # Direction (Vector3D)
            self.dir = dir
        else:
            self.dir = Vector([0.0, 0.0, 0.0])
            
    # Generate the samples
    def GenerateSamples(self, sqrtNumSamples, numBands):
    
        numSamples = sqrtNumSamples * sqrtNumSamples
        numFunctions = numBands * numBands
        
        # Create the array to hold the samples
        samples = []
        for x in xrange(numSamples):
            samples[x].append(SPHSample(0.0, 0.0, Vector([0.0, 0.0, 0.0]), 0.0))
            
        for i in xrange(sqrtNumSamples):
            for j in xrange(sqrtNumSamples):
                # The code 
        
        # Return the samples array
        return samples