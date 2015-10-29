import random
import six.moves as sm
import math
from gem import vector
from gem.experimental import sph


# Spherical Harmonics Sample Class
class SPHSample (object):
    def __init__(self, theta, phi, dirc, sampleNumber):
        # Spherical coordinates
        self.theta = theta
        self.phi = phi

        # Values of SH function at this point
        self.values = []
        for _ in sm.range(sampleNumber):
            self.values.append(0.0)

        # Direction (Vector3D)
        if isinstance(dirc, vector.Vector):
            self.dir = dirc
        else:
            self.dir = vector.Vector(3, data=[0.0, 0.0, 0.0])


# Generate the samples
def GenerateSamples(sqrtNumSamples, numBands):
    x = 0
    y = 0
    theta = 0
    phi = 0
    ii = 0
    index = 0

    numSamples = sqrtNumSamples * sqrtNumSamples
    numFunctions = numBands * numBands
    invertedNumSamples = 1.0 / sqrtNumSamples

    # Create the array to hold the samples
    samples = []
    for x in sm.range(numSamples):
        samples.append(SPHSample(0.0, 0.0, vector.Vector(3, data=[0.0, 0.0, 0.0]), numFunctions))

    print ("Generating Samples...")
    # Loop through a grid of numSamples X numSamples
    for i in sm.range(sqrtNumSamples):
        for j in sm.range(sqrtNumSamples):
            x = (i + random.random()) * invertedNumSamples
            y = (j + random.random()) * invertedNumSamples

            # Spherical Angles
            theta = 2.0 * math.acos(math.sqrt(1.0 - x))
            phi = 2.0 * math.PI * y

            samples[ii].theta = theta
            samples[ii].phi = phi

            # Vec is an array and dir is a vector
            samples[ii].dir.vec = [math.sin(theta) * math.cos(phi),
                                   math.sin(theta) * math.sin(phi),
                                   math.cos(theta)]

            # Calculate SH coefficients of current sample
            for l in sm.range(numBands):
                for m in sm.range(-l, l + 1):
                    index = l * (l + 1) + m
                    samples[ii].values[index] = sph.SPH(l, m, theta, phi)
            ii += 1

    # Return the samples array
    return samples
