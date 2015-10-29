import pprint
import math
import struct
import six.moves as sm

from gem.common import sinc


# Calculate Spherical Harmonics coefficients using irradiance maps
class SPH_IrradianceMapCoeff(object):

    def __init__(self, fileU, width, height):
        # File details
        self.file = fileU
        self.width = width
        self.height = height
        # Alocate space for the image
        self.hdr = [[[0 for _ in sm.range(3)] for _ in sm.range(self.width)] for _ in sm.range(self.height)]
        # Alocate space for the coefficients
        self.coeffs = [[0 for _ in sm.range(3)] for _ in sm.range(9)]
        # Load the image
        self.load()

    # Load the file
    def load(self):
        with open(self.file, 'rb') as f:
            for i in sm.range(self.width):
                for j in sm.range(self.height):
                    # The depth is 3 (R,G,B)
                    for k in sm.range(3):
                        # Decode the 32 bit binary into float point number
                        # 4 bytes
                        self.hdr[i][j][k] = struct.unpack('f', f.read(4))[0]

        # Calculate coefficients
        self.calculateCoefficients()

    def calculateCoefficients(self):
        for i in sm.range(self.width):
            for j in sm.range(self.height):
                v = (self.width / 2.0 - i) / (self.width / 2.0)
                u = (j - self.width / 2.0) / (self.width / 2.0)
                r = math.sqrt(u * u + v * v)
                # Only withn a unit circle
                if r <= 1.0:
                    theta = math.pi * r
                    phi = math.atan2(v, u)

                    x = math.sin(theta) * math.cos(phi)
                    y = math.sin(theta) * math.sin(phi)
                    z = math.cos(theta)

                    domega = (2 * math.pi / self.width) * (2.0 * math.pi / self.width) * sinc(theta)

                    self.updateCoefficients(self.hdr[i][j], domega, x, y, z)

    def updateCoefficients(self, hdr, domega, x, y, z):
        for col in sm.range(3):

            c = 0.282095
            self.coeffs[0][col] += hdr[col] * c * domega

            c = 0.488603
            self.coeffs[1][col] += hdr[col] * (c * y) * domega
            self.coeffs[2][col] += hdr[col] * (c * z) * domega
            self.coeffs[3][col] += hdr[col] * (c * x) * domega

            c = 1.092548
            self.coeffs[4][col] += hdr[col] * (c * x * y) * domega
            self.coeffs[5][col] += hdr[col] * (c * y * z) * domega
            self.coeffs[7][col] += hdr[col] * (c * x * z) * domega

            c = 0.315392
            self.coeffs[6][col] += hdr[col] * (c * (3 * z * z - 1)) * domega

            c = 0.546274
            self.coeffs[8][col] += hdr[col] * (c * (x * x - y * y)) * domega

    def output(self):
        pprint.pprint(self.coeffs)

# Example
#test_irradiance_map = SPH_IrradianceMapCoeff("uffizi_probe.float", 1500, 1500)
#test_irradiance_map.output()
