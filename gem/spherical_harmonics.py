"""Real spherical harmonics, angular-probe projection and diffuse convolution.

Canonical coefficients use Condon–Shortley phase and l*(l+1)+m indexing.
RGB radiance coefficients are distinct from cosine-convolved irradiance.
The legacy probe class retains its historical nine-polynomial basis.
"""
import math
import pprint
import random
import struct
from gem.vector import Vector

from gem.legendre import Legendre


# n! where n >= 0
def Factorial(n):
    """Historical factorial helper; supported domain is nonnegative integers."""
    if n <= 1:
        return 1

    result = n

    while n > 1:
        n -= 1
        result *= n

    return result

# Normalization constant for a Spherical Harmonic function
def K(l, m):
    """SH normalization for integer degrees with 0 <= m <= l."""
    K = ((2.0 * l + 1.0) * Factorial(l - m)) / ((4.0 * math.pi) * Factorial(l + m))
    return math.sqrt(K)


# Sample a Spherical Harmonic function Y(l, m) at a point on the unit sphere
def SPH(l, m, theta, phi):
    """Real SH with Condon–Shortley phase, polar theta and azimuth phi."""
    root2 = math.sqrt(2.0)
    if m == 0:
        return K(l, 0) * Legendre(l, m, math.cos(theta)).run()
    elif m > 0:
        return root2 * K(l, m) * math.cos(m * phi) * Legendre(l, m, math.cos(theta)).run()
    elif m < 0:
        return root2 * K(l, -m) * math.sin(-m * phi) * Legendre(l, -m, math.cos(theta)).run()
    else:
        print ("WTF... The m is ...")
        return 0

class SPHSample (object):
    def __init__(self, theta, phi, dirc, sampleNumber):
        # Spherical coordinates
        self.theta = theta
        self.phi = phi

        # Values of SH function at this point
        self.values = []
        for _ in range(sampleNumber):
            self.values.append(0.0)

        # Direction (Vector3D)
        if isinstance(dirc, Vector):
            self.dir = dirc
        else:
            self.dir = Vector(3, data=[0.0, 0.0, 0.0])



def GenerateSamples(sqrtNumSamples, numBands):
    """Generate jittered equal-solid-angle strata; seed random for repeatability.

    Each sample has integration weight 4*pi/sqrtNumSamples**2.
    Constructor direction ownership is unchanged. Legacy invalid-argument
    arithmetic is not replaced with a generalized validation policy.
    """
    inverse = 1.0 / sqrtNumSamples
    samples = []
    for i in range(sqrtNumSamples):
        for j in range(sqrtNumSamples):
            u = (i+random.random())*inverse
            v = (j+random.random())*inverse
            theta = 2.0*math.acos(math.sqrt(1.0-u))
            phi = 2.0*math.pi*v
            direction = Vector(3, [math.sin(theta)*math.cos(phi),
                                   math.sin(theta)*math.sin(phi), math.cos(theta)])
            sample = SPHSample(theta, phi, direction, numBands*numBands)
            for l in range(numBands):
                for m in range(-l,l+1):
                    sample.values[l*(l+1)+m] = SPH(l,m,theta,phi)
            samples.append(sample)
    return samples


def _positive_integer(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError(name+" must be a positive integer")
    return value


def _rgb(values):
    try:
        valid = len(values) == 3 and all(math.isfinite(v) for v in values)
    except (TypeError, ValueError, OverflowError):
        valid = False
    if not valid:
        raise ValueError("RGB values must contain three finite components")
    return list(values)


def _coefficient_bands(coefficients):
    if not coefficients:
        raise ValueError("coefficients must contain complete nonempty bands")
    bands = math.isqrt(len(coefficients))
    if bands*bands != len(coefficients):
        raise ValueError("coefficients must contain complete bands")
    for value in coefficients:
        _rgb(value)
    return bands


def _basis(bands, theta, phi):
    return [SPH(l,m,theta,phi) for l in range(bands) for m in range(-l,l+1)]


def project_radiance(samples, radiances, weights=None):
    """Project RGB radiances at SPHSamples into fresh canonical coefficients.

    Omitted weights assume uniform sphere samples with weight 4*pi/N.
    Explicit weights are finite nonnegative solid angles in steradians.
    Samples must contain matching precomputed basis values.
    """
    if not samples or len(samples) != len(radiances):
        raise ValueError("require matching nonempty samples and radiances")
    count = len(samples[0].values)
    bands = math.isqrt(count)
    if count == 0 or bands*bands != count:
        raise ValueError("sample values must contain complete nonempty bands")
    rgb = [_rgb(value) for value in radiances]
    if weights is None:
        weights = [4*math.pi/len(samples)]*len(samples)
    try:
        valid_weights = len(weights) == len(samples) and all(math.isfinite(w) and w >= 0 for w in weights)
    except (TypeError, ValueError, OverflowError):
        valid_weights = False
    if not valid_weights:
        raise ValueError("weights must be matching finite nonnegative solid angles")
    if any(len(sample.values) != count or not all(math.isfinite(v) for v in sample.values)
           for sample in samples):
        raise ValueError("sample basis values must be matching and finite")
    return [[math.fsum(color[channel]*sample.values[index]*weight
                      for sample,color,weight in zip(samples,rgb,weights))
             for channel in range(3)] for index in range(count)]


def _angular_pixels(hdr):
    if not hdr or not hdr[0]:
        raise ValueError("angular probe must be a nonempty rectangular RGB image")
    height, width = len(hdr), len(hdr[0])
    try:
        rectangular = all(len(row) == width for row in hdr)
    except TypeError:
        rectangular = False
    if not rectangular:
        raise ValueError("angular probe rows must have equal widths")
    for row in hdr:
        for pixel in row:
            _rgb(pixel)
    area = 4*math.pi*math.pi/(width*height)
    for row in range(height):
        v = 1-2*(row+0.5)/height
        for column in range(width):
            u = 2*(column+0.5)/width-1
            radius = math.hypot(u,v)
            if radius <= 1:
                theta = math.pi*radius
                phi = math.atan2(v,u)
                sinc = math.sin(theta)/theta if theta else 1.0
                yield hdr[row][column], area*sinc, theta, phi


def project_angular_probe(hdr, numBands=3):
    """Project row-major RGB angular-disk pixels using center quadrature.

    A rectangular image stretches the disk independently along each axis.
    No weight renormalization is performed. Finite pixel integration is an
    approximation; latitude-longitude and mirrored-ball maps are unsupported.
    """
    bands = _positive_integer(numBands, "numBands")
    # Compensated accumulation uses memory proportional to coefficient count.
    totals = [[0.0]*3 for _ in range(bands*bands)]
    corrections = [[0.0]*3 for _ in range(bands*bands)]
    for color,weight,theta,phi in _angular_pixels(hdr):
        for index,value in enumerate(_basis(bands,theta,phi)):
            for channel in range(3):
                term = color[channel]*weight*value-corrections[index][channel]
                updated = totals[index][channel]+term
                corrections[index][channel] = (updated-totals[index][channel])-term
                totals[index][channel] = updated
    return totals


def reconstruct(coefficients, direction):
    """Reconstruct RGB from canonical coefficients at a caller-supplied unit direction.

    Accept Vector3 or XYZ triples. Unit length is a prerequisite; no implicit
    normalization or cosine convolution is performed.
    """
    bands = _coefficient_bands(coefficients)
    values = direction.vector if isinstance(direction, Vector) else direction
    try:
        valid = len(values) == 3 and all(math.isfinite(v) for v in values) and any(values)
    except (TypeError, ValueError, OverflowError):
        valid = False
    if not valid:
        raise ValueError("direction must be a finite nonzero unit XYZ triple")
    x,y,z = values
    if not -1 <= z <= 1:
        raise ValueError("unit direction Z must be within [-1,1]")
    basis = _basis(bands, math.acos(z), math.atan2(y,x))
    return [math.fsum(coefficient[channel]*value for coefficient,value in zip(coefficients,basis))
            for channel in range(3)]


def convolve_diffuse(radiance_coefficients):
    """Return separate irradiance coefficients for one to three SH bands.

    Cosine convolution factors are pi, 2*pi/3 and pi/4. Reconstruction does
    not apply convolution again. Lambertian outgoing radiance additionally
    multiplies irradiance by albedo/pi; that is a rendering operation.
    """
    bands = _coefficient_bands(radiance_coefficients)
    if bands > 3:
        raise ValueError("diffuse convolution supports at most three bands")
    factors = [math.pi, 2*math.pi/3, math.pi/4]
    return [[v*factors[math.isqrt(index)] for v in coefficient]
            for index,coefficient in enumerate(radiance_coefficients)]


def legacy_to_canonical(coefficients):
    """Convert nine legacy probe radiance coefficients, including rounded scales."""
    if _coefficient_bands(coefficients) != 3:
        raise ValueError("legacy probe conversion requires nine RGB coefficients")
    first = math.sqrt(3/(4*math.pi))/0.488603
    cross = math.sqrt(15/(4*math.pi))/1.092548
    ratios = [math.sqrt(1/(4*math.pi))/0.282095, -first, first, -first,
              cross, -cross, math.sqrt(5/(16*math.pi))/0.315392,
              -cross, math.sqrt(15/(16*math.pi))/0.546274]
    return [[value*ratios[index] for value in coefficient]
            for index,coefficient in enumerate(coefficients)]

class SPH_IrradianceMapCoeff(object):
    """Raw native-endian RGB float32 angular probe, retaining the legacy basis.

    coeffs are nine radiance coefficients, not cosine-convolved irradiance.
    Use legacy_to_canonical before canonical reconstruction or convolution.
    """
    def __init__(self, fileU, width, height):
        self.file = fileU
        self.width = _positive_integer(width, "width")
        self.height = _positive_integer(height, "height")
        self.hdr = []
        self.coeffs = [[0.0]*3 for _ in range(9)]
        self.load()

    def load(self):
        required = self.width*self.height*3*4
        with open(self.file, 'rb') as handle:
            data = handle.read(required)
        if len(data) != required:
            raise ValueError("probe file does not contain enough RGB float32 data")
        values = struct.unpack('%df' % (self.width*self.height*3), data)
        self.hdr = [[list(values[(row*self.width+column)*3:(row*self.width+column+1)*3])
                     for column in range(self.width)] for row in range(self.height)]
        self.calculateCoefficients()

    def calculateCoefficients(self):
        self.coeffs = [[0.0]*3 for _ in range(9)]
        for color,weight,theta,phi in _angular_pixels(self.hdr):
            self.updateCoefficients(color,weight,math.sin(theta)*math.cos(phi),
                                    math.sin(theta)*math.sin(phi),math.cos(theta))

    def updateCoefficients(self, hdr, domega, x, y, z):
        for col in range(3):

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
