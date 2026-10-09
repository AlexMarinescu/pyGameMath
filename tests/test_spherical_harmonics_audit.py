"""Independent SH audit: Cartesian polynomials, Rodrigues and exact quadrature.

No gem basis, matrix conversion or quaternion rotation is used as an oracle.
Confirmed near-pole losses are strict expected failures, not repaired here.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import copy
import math
import random
import struct

import pytest

from gem import spherical_harmonics as sh
from gem.quaternion import Quaternion
from gem.vector import Vector


def cartesian_basis(direction):
    x, y, z = direction
    a = math.sqrt(3/(4*math.pi))
    b = math.sqrt(15/(4*math.pi))
    return [1/math.sqrt(4*math.pi), -a*y, a*z, -a*x, b*x*y, -b*y*z,
            math.sqrt(5/(16*math.pi))*(3*z*z-1), -b*x*z,
            math.sqrt(15/(16*math.pi))*(x*x-y*y)]


def direction(theta, phi):
    return [math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)]


@lru_cache(maxsize=None)
def rodrigues_coefficients(degree, order=0):
    # Differentiate the explicit Rodrigues polynomial. No recurrence.
    return tuple((degree-2*k-order, Fraction(
        (-1)**k*math.factorial(2*degree-2*k),
        2**degree*math.factorial(k)*math.factorial(degree-k)*math.factorial(degree-2*k-order)))
        for k in range((degree-order)//2+1))


def normalization_reference(degree, order):
    with localcontext() as ctx:
        ctx.prec = 180
        return ((Decimal(2*degree+1)*Decimal(math.factorial(degree-order)))
                / (4*Decimal.from_float(math.pi)*Decimal(math.factorial(degree+order)))).sqrt()


def spherical_reference(degree, order, theta, phi):
    x = Fraction(math.cos(theta))
    m = abs(order)
    derivative = sum(c*x**p for p, c in rodrigues_coefficients(degree, m))
    base = 1-x*x
    with localcontext() as ctx:
        ctx.prec = 180
        value = Decimal(derivative.numerator)/Decimal(derivative.denominator)
        if m:
            value *= (Decimal(base.numerator)/Decimal(base.denominator)).sqrt()**m
            value *= Decimal((-1)**m)*Decimal(2).sqrt()
            value *= Decimal.from_float(math.cos(m*phi) if order > 0 else math.sin(m*phi))
        return float(value*normalization_reference(degree, m))


def decimal_polynomial(coefficients, x):
    return sum(Decimal(c.numerator)/Decimal(c.denominator)*x**p if p else
               Decimal(c.numerator)/Decimal(c.denominator) for p, c in coefficients)


@lru_cache(maxsize=None)
def gauss_nodes(count=12):
    # Newton solves explicit Rodrigues polynomials in Decimal, independently
    # of gem.Legendre; weights integrate polynomials through degree 2N-1.
    polynomial = rodrigues_coefficients(count)
    derivative = tuple((p-1, c*p) for p, c in polynomial if p)
    nodes = []
    with localcontext() as ctx:
        ctx.prec = 80
        for i in range(count):
            x = Decimal.from_float(math.cos(math.pi*(i+.75)/(count+.5)))
            for _ in range(30):
                step = decimal_polynomial(polynomial, x)/decimal_polynomial(derivative, x)
                x -= step
                if abs(step) < Decimal('1e-65'):
                    break
            else:
                raise AssertionError('independent Gauss node did not converge')
            d = decimal_polynomial(derivative, x)
            nodes.append((float(x), float(2/((1-x*x)*d*d))))
    return tuple(sorted(nodes))


@lru_cache(maxsize=None)
def sphere_quadrature():
    rows = []
    for z, weight in gauss_nodes():
        r = math.sqrt((1-z)*(1+z))
        for j in range(32):
            phi = 2*math.pi*(j+.5)/32
            rows.append(([r*math.cos(phi), r*math.sin(phi), z], weight*2*math.pi/32))
    return rows


def orientation(axis, angle):
    n = math.hypot(*axis)
    return Quaternion([math.cos(angle/2)]+[v/n*math.sin(angle/2) for v in axis])


def rodrigues_rotate(point, axis, angle):
    a = [v/math.hypot(*axis) for v in axis]
    x, y, z = point
    cross = [a[1]*z-a[2]*y, a[2]*x-a[0]*z, a[0]*y-a[1]*x]
    dot = math.fsum(u*v for u, v in zip(a, point))
    c, s = math.cos(angle), math.sin(angle)
    return [c*v+s*w+(1-c)*dot*u for v, w, u in zip(point, cross, a)]


@pytest.mark.parametrize('degree', [0, 1, 2, 3, 8, 12, 20, 32, 64, 85])
@pytest.mark.parametrize('x', [-.875, -.25, 0., .25, .875])
def test_rodrigues_normalization_and_negative_orders(degree, x):
    theta, phi = math.acos(x), .73
    for m in sorted({0, min(1, degree), degree//2, degree}):
        assert sh.K(degree, m) == pytest.approx(float(normalization_reference(degree, m)), rel=4e-15, abs=0)
        for order in ({m, -m} if m else {0}):
            expected = spherical_reference(degree, order, theta, phi)
            actual = sh.SPH(degree, order, theta, phi)
            assert actual == pytest.approx(expected, rel=8e-12, abs=3e-13)


@pytest.mark.parametrize('theta,phi', [(0, .7), (math.pi, -1.2), (math.pi/2, 0),
    (math.pi/2, math.pi/2), (.73, 1.27), (2.3, -.9)])
def test_low_band_cartesian_signs_and_index(theta, phi):
    expected = cartesian_basis(direction(theta, phi))
    assert [l*(l+1)+m for l in range(3) for m in range(-l, l+1)] == list(range(9))
    actual = [sh.SPH(l, m, theta, phi) for l in range(3) for m in range(-l, l+1)]
    assert actual == pytest.approx(expected, abs=3e-15)
    assert sh._basis(3, theta, phi) == pytest.approx(expected, abs=3e-15)


@pytest.mark.parametrize('degree', [0, 1, 2, 3, 6, 12, 32, 64, 85])
def test_addition_theorem_parity_and_poles(degree):
    for theta, phi in [(.71, .93), (math.pi/2, -.53), (2.4, 1.13)]:
        values = [sh.SPH(degree, m, theta, phi) for m in range(-degree, degree+1)]
        assert math.fsum(v*v for v in values) == pytest.approx((2*degree+1)/(4*math.pi), rel=2e-13)
        antipodal = [sh.SPH(degree, m, math.pi-theta, phi+math.pi) for m in range(-degree, degree+1)]
        assert antipodal == pytest.approx([(-1)**degree*v for v in values], abs=1e-13)
    for theta, sign in [(0., 1), (math.pi, (-1)**degree)]:
        assert sh.SPH(degree, 0, theta, .3) == pytest.approx(sign*math.sqrt((2*degree+1)/(4*math.pi)))
        assert all(sh.SPH(degree, m, theta, .3) == 0 for m in range(-degree, degree+1) if m)


def test_independent_gauss_fourier_orthonormality_through_degree_six():
    terms = [(l, m) for l in range(7) for m in range(-l, l+1)]
    values = [[sh.SPH(l, m, math.acos(d[2]), math.atan2(d[1], d[0])) for l, m in terms]
              for d, _ in sphere_quadrature()]
    weights = [w for _, w in sphere_quadrature()]
    assert math.fsum(weights) == pytest.approx(4*math.pi, abs=4e-15)
    for i in range(len(terms)):
        for j in range(i+1):
            integral = math.fsum(row[i]*row[j]*w for row, w in zip(values, weights))
            assert integral == pytest.approx(float(i == j), abs=8e-15)


def analytical_signal(d):
    x, y, z = d
    return [2+.3*x+.2*x*y, 3-.7*y+.8*z*z, 4+.4*z+.5*(x*x-y*y)]


def signal_coefficients():
    a = math.sqrt(3/(4*math.pi)); b = math.sqrt(15/(4*math.pi))
    d = math.sqrt(5/(16*math.pi)); e = math.sqrt(15/(16*math.pi))
    c = [[0.]*3 for _ in range(9)]
    c[0] = [v*math.sqrt(4*math.pi) for v in [2, 3+.8/3, 4]]
    c[3][0] = -.3/a; c[4][0] = .2/b; c[1][1] = .7/a
    c[6][1] = .8/(3*d); c[2][2] = .4/a; c[8][2] = .5/e
    return c


def analytical_irradiance(n):
    x, y, z = n
    return [2*math.pi+2*math.pi*.3*x/3+math.pi*.2*x*y/4,
            math.pi*(3+.8/3)-2*math.pi*.7*y/3+math.pi*.8*(z*z-1/3)/4,
            4*math.pi+2*math.pi*.4*z/3+math.pi*.5*(x*x-y*y)/4]


@pytest.mark.parametrize('bands', [1, 3, 5, 7])
def test_projection_of_low_degree_rgb_and_storage(bands):
    samples, colors, weights = [], [], []
    for d, weight in sphere_quadrature():
        sample = sh.SPHSample(math.acos(d[2]), math.atan2(d[1], d[0]), Vector(3, d[:]), bands*bands)
        # Independent polynomial values, not gem.SPH, drive the integration.
        sample.values = [spherical_reference(l, m, sample.theta, sample.phi)
                         for l in range(bands) for m in range(-l, l+1)]
        samples.append(sample); colors.append(analytical_signal(d)); weights.append(weight)
    saved = copy.deepcopy([(s.dir.vector, s.values) for s in samples])
    stores = [(s.dir.vector, s.values) for s in samples]
    original = copy.deepcopy(colors); w_original = weights[:]
    expected = signal_coefficients()[:bands*bands]+[[0.]*3 for _ in range(max(0, bands*bands-9))]
    for _ in range(2):
        result = sh.project_radiance(samples, colors, weights)
        for row, reference in zip(result, expected):
            assert row == pytest.approx(reference, abs=2e-14)
        assert len({id(row) for row in result}) == bands*bands
        result[0][0] = -123
        assert colors == original and weights == w_original
    assert [(s.dir.vector, s.values) for s in samples] == saved
    assert all(s.dir.vector is a and s.values is b for s, (a, b) in zip(samples, stores))


@pytest.mark.parametrize('n', [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (.3,.4,math.sqrt(.75))])
def test_reconstruction_and_diffuse_direct_hemisphere_integral(n):
    c = signal_coefficients(); saved = copy.deepcopy(c)
    assert sh.reconstruct(c, n) == pytest.approx(analytical_signal(n), abs=3e-15)
    irradiance = sh.convolve_diffuse(c)
    expected = analytical_irradiance(n)
    assert sh.reconstruct(irradiance, Vector(3, list(n))) == pytest.approx(expected, abs=6e-15)
    # Integrate L(d)*(n.d) over a parameterized hemisphere, independently.
    axis = [0, 1, 0] if abs(n[1]) < .9 else [1, 0, 0]
    u = [axis[1]*n[2]-axis[2]*n[1], axis[2]*n[0]-axis[0]*n[2], axis[0]*n[1]-axis[1]*n[0]]
    u = [v/math.hypot(*u) for v in u]
    v = [n[1]*u[2]-n[2]*u[1], n[2]*u[0]-n[0]*u[2], n[0]*u[1]-n[1]*u[0]]
    terms = [[] for _ in range(3)]
    for z, w in gauss_nodes():
        z = (z+1)/2; r = math.sqrt(1-z*z)
        for j in range(32):
            phi = 2*math.pi*(j+.5)/32
            d = [z*a+r*(math.cos(phi)*b+math.sin(phi)*c) for a,b,c in zip(n,u,v)]
            for channel, value in enumerate(analytical_signal(d)):
                terms[channel].append(value*z*w*math.pi/32)
    assert [math.fsum(t) for t in terms] == pytest.approx(expected, abs=8e-15)
    assert c == saved and all(a is not b for a,b in zip(c, irradiance))


@pytest.mark.parametrize('axis', [(1,0,0), (0,1,0), (0,0,1), (1,-2,3)])
@pytest.mark.parametrize('angle', [0., math.pi/2, -math.pi/2, math.pi, .73])
@pytest.mark.parametrize('count', [1, 4, 9])
def test_rotation_covariance_energy_ownership_and_sign(axis, angle, count):
    rng = random.Random(473001+count)
    c = [[rng.uniform(-2,2) for _ in range(3)] for _ in range(count)]
    saved = copy.deepcopy(c); q = orientation(axis, angle); qstore = q.data; qsaved = q.data[:]
    out = sh.rotate_coefficients(c, q)
    for d in [(1,0,0), (0,1,0), (0,0,1), (.3,.4,math.sqrt(.75))]:
        b = cartesian_basis(rodrigues_rotate(d, axis, -angle))[:count]
        expected = [math.fsum(row[k]*value for row,value in zip(c,b)) for k in range(3)]
        assert sh.reconstruct(out, d) == pytest.approx(expected, abs=9e-15)
    for start, stop in [(0,1), (1,4), (4,9)]:
        for channel in range(3):
            assert math.hypot(*(row[channel] for row in out[start:stop])) == pytest.approx(
                math.hypot(*(row[channel] for row in c[start:stop])), rel=3e-15, abs=1e-15)
    inverse = sh.rotate_coefficients(out, orientation(axis, -angle))
    negative = sh.rotate_coefficients(c, Quaternion([-v for v in q.data]))
    for a,b,ref in zip(inverse,negative,c):
        assert a == pytest.approx(ref, abs=7e-15)
    for a,b in zip(out, negative):
        assert a == b
    assert c == saved and q.data is qstore and q.data == qsaved
    assert all(a is not b for a,b in zip(out, c))
    out[0][0] = 999
    assert c == saved


@pytest.mark.parametrize('scale', [1e-280, 1., 1e280])
@pytest.mark.parametrize('rgb', [False, True])
def test_rotation_mixed_channels_scales_and_noncommuting_composition(scale, rgb):
    base = [(-1)**i*(i+.25)/8 for i in range(9)]
    c = [[scale*v, scale*v/8, 0.] for v in base] if rgb else [scale*v for v in base]
    saved = copy.deepcopy(c)
    a = orientation((1,0,0), .7); b = orientation((0,1,0), -.9)
    out = sh.rotate_coefficients(sh.rotate_coefficients(c, a), b)
    reverse = sh.rotate_coefficients(sh.rotate_coefficients(c, b), a)
    for d in [(.3,.4,math.sqrt(.75)), (1,0,0)]:
        original = rodrigues_rotate(rodrigues_rotate(d, (0,1,0), .9), (1,0,0), -.7)
        expected = math.fsum(v*y for v,y in zip(base,cartesian_basis(original)))
        actual = math.fsum((row[0] if rgb else row)/scale*y for row,y in zip(out,cartesian_basis(d)))
        assert actual == pytest.approx(expected, abs=6e-15)
    flat = lambda values: [v for row in values for v in row] if rgb else values
    assert max(abs(x/scale-y/scale) for x,y in zip(flat(out),flat(reverse))) > .1
    if rgb:
        assert all(row[2] == 0 for row in out)
        assert [row[0]/scale for row in out] == pytest.approx([8*row[1]/scale for row in out], abs=1e-15)
    assert c == saved


@pytest.mark.parametrize('shape', [(3,3), (5,3), (8,4)])
def test_angular_pixel_area_signs_and_legacy_conversion(tmp_path, shape):
    width, height = shape
    image = [[[.125*(1+r), .25*(1+c), .5] for c in range(width)] for r in range(height)]
    saved = copy.deepcopy(image)
    expected = [[[] for _ in range(3)] for _ in range(9)]
    for row in range(height):
        for col in range(width):
            u = (2*col+1-width)/width; v = (height-2*row-1)/height
            radius = math.hypot(u,v)
            if radius > 1: continue
            theta = math.pi*radius
            d = [math.sin(theta)*u/radius, math.sin(theta)*v/radius, math.cos(theta)] if radius else [0,0,1]
            weight = 4*math.pi**2/(width*height)*(math.sin(theta)/theta if radius else 1.)
            for index,y in enumerate(cartesian_basis(d)):
                for channel in range(3): expected[index][channel].append(image[row][col][channel]*weight*y)
    expected = [[math.fsum(terms) for terms in row] for row in expected]
    canonical = sh.project_angular_probe(image)
    for actual, ref in zip(canonical,expected): assert actual == pytest.approx(ref, abs=5e-15)
    raw = tmp_path/'probe.float'
    raw.write_bytes(struct.pack(f'{width*height*3}f', *(v for row in image for pixel in row for v in pixel)))
    probe = sh.SPH_IrradianceMapCoeff(str(raw), width, height)
    original = copy.deepcopy(probe.coeffs)
    for _ in range(2):
        probe.calculateCoefficients()
        converted = sh.legacy_to_canonical(probe.coeffs)
        for actual,ref in zip(converted, expected): assert actual == pytest.approx(ref, abs=7e-15)
    assert probe.coeffs == original and image == saved


@pytest.mark.parametrize('bad', [[], [[1,2,3]]*2, [[1,2]], [[1,2,math.nan]], [[1,2,math.inf]]])
def test_checked_coefficient_layout_and_nonfinite_domains(bad):
    for function in [lambda: sh.reconstruct(bad,[0,0,1]), lambda: sh.convolve_diffuse(bad),
                     lambda: sh.rotate_coefficients(bad,Quaternion())]:
        with pytest.raises(ValueError): function()


def test_zero_energy_channel_aliases_and_independent_results():
    shared = [0.,0.,0.]; c = [shared]*9
    q = orientation((1,2,3), .8)
    for result in [sh.rotate_coefficients(c,q), sh.convolve_diffuse(c), sh.legacy_to_canonical(c)]:
        assert result == [[0.]*3 for _ in range(9)]
        assert len({id(row) for row in result}) == 9
        result[0][0] = 1
        assert c == [[0.]*3 for _ in range(9)]
    assert sh.reconstruct(c,[0,0,1]) == [0.,0.,0.]


@pytest.mark.parametrize('theta,phi,l,m', [
    (1e-9, 0., 1, 1), (1e-12, .7, 1, -1), (1e-100, 0., 1, 1),
    (math.pi-1e-9, 0., 1, 1), (math.pi-1e-12, .7, 2, -1), (1e-9, .3, 2, 2)])
@pytest.mark.defect('4G3-A01: near-pole transverse information lost through cos(theta)')
def test_near_pole_basis_retains_transverse_components(theta, phi, l, m):
    expected = cartesian_basis(direction(theta, phi))[l*(l+1)+m]
    assert expected != 0 and math.isfinite(expected)
    assert sh.SPH(l,m,theta,phi) == pytest.approx(expected, rel=3e-14, abs=0)


@pytest.mark.parametrize('d,index', [([1e-9,0.,1.],3), ([0.,1e-9,1.],1),
                                     ([1e-9,0.,-1.],7), ([0.,1e-9,-1.],5)])
@pytest.mark.defect('4G3-A01: reconstruction discards supplied transverse components')
def test_near_pole_reconstruction_retains_unit_direction(d, index):
    # Unit sphere coordinates rounded independently to binary64; hypot is 1.
    assert math.hypot(*d) == 1.
    c = [[0.]*3 for _ in range(9)]; c[index] = [1.,2.,4.]
    saved = copy.deepcopy(c); position = Vector(3,d[:]); storage = position.vector
    expected = cartesian_basis(d)[index]
    out = sh.reconstruct(c, position)
    assert c == saved and position.vector is storage and position.vector == d
    assert out == pytest.approx([expected, 2*expected, 4*expected], rel=3e-14, abs=0)


@pytest.mark.parametrize('light_axis', [(1,0,0), (0,1,0), (.8,0,.6)])
def test_directional_projection_and_truncated_cosine_kernel(light_axis):
    sample = sh.SPHSample(0., 0., Vector(3,list(light_axis)), 9)
    sample.values = cartesian_basis(light_axis)
    projected = sh.project_radiance([sample], [[1.,2.,4.]], [1.])
    convolved = sh.convolve_diffuse(projected)
    for n in [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (.3,.4,math.sqrt(.75))]:
        u = math.fsum(a*b for a,b in zip(n,light_axis))
        # Addition theorem and exact cosine-kernel factors, no SH oracle.
        kernel = .25+.5*u+5*(3*u*u-1)/32
        assert sh.reconstruct(convolved,n) == pytest.approx([kernel,2*kernel,4*kernel], abs=3e-15)
    # The truncated kernel at an aligned axis is 17/16, not exact max(u,0).
    assert sh.reconstruct(convolved,light_axis)[0] == pytest.approx(17/16, abs=2e-15)


def test_sharp_hdr_golden_against_continuous_zonal_integrals():
    import json
    from pathlib import Path
    from examples.hdr_sh import reference
    manifest = json.loads((Path(__file__).resolve().parents[1]/'examples/output/visualization.json').read_text())
    # For g(u)=max(u,0)^32, Funk-Hecke gives c_lm=2pi*I_l*Y_lm(axis),
    # I_l=int_0^1 u^32 P_l(u)du. These first moments are exact rationals.
    integrals = [Fraction(1,33), Fraction(1,34), (Fraction(3,35)-Fraction(1,33))/2]
    peaks = [45.,30.,18.]; ambient = [.08,.06,.04]
    original = []
    for name,axis in [('sh_original',[.8,0.,.6]), ('sh_rotated',[0.,.8,.6])]:
        radiance = [[2*math.pi*float(integrals[math.isqrt(i)])*y*p for p in peaks]
                    for i,y in enumerate(cartesian_basis(axis))]
        for channel,value in enumerate(ambient): radiance[0][channel] += value*math.sqrt(4*math.pi)
        stored = manifest['radiance_coefficients'] if name == 'sh_original' else manifest['rotated_radiance_coefficients']
        for actual,expected in zip(stored,radiance):
            # Integration, not floating roundoff, dominates at 128x64.
            assert actual == pytest.approx(expected, abs=5e-4)
        coefficients = [[v*[math.pi,2*math.pi/3,math.pi/4][math.isqrt(i)] for v in row]
                        for i,row in enumerate(radiance)]
        for pixel in manifest['images'][name]['reference_pixels']:
            b = cartesian_basis(pixel['normal'])
            expected = [math.fsum(row[k]*y for row,y in zip(coefficients,b)) for k in range(3)]
            assert pixel['irradiance_linear_rgb'] == pytest.approx(expected, abs=9e-4)
            lambertian = [v*.65/math.pi for v in expected]
            assert pixel['lambertian_linear_rgb'] == pytest.approx(lambertian, abs=2e-4)
            assert reference.reflected_radiance(expected,[.65]*3) == pytest.approx(lambertian, abs=1e-15)
        original.append(radiance)
    assert original[0][3][0] < 0 and original[1][1][0] < 0
    assert original[0][1][0] == original[1][3][0] == 0


@pytest.mark.parametrize('resolution', [16, 32, 64])
def test_angular_and_latlong_constant_quadrature_with_declared_error(resolution):
    from examples.hdr_sh import environment
    width, height = resolution, resolution//2
    color = [1.,2.,4.]
    image = [[color[:] for _ in range(width)] for _ in range(height)]
    for mapping, project in [('angular', sh.project_angular_probe), ('latlong', environment.project_latlong)]:
        terms = [[] for _ in range(9)]
        for row in range(height):
            for col in range(width):
                if mapping == 'angular':
                    u = (2*col+1-width)/width; v = (height-2*row-1)/height
                    r = math.hypot(u,v)
                    if r > 1: continue
                    theta = math.pi*r
                    d = [math.sin(theta)*u/r,math.sin(theta)*v/r,math.cos(theta)] if r else [0,0,1]
                    weight = 4*math.pi**2/(width*height)*(math.sin(theta)/theta if r else 1)
                else:
                    theta = math.pi*(row+.5)/height; phi = 2*math.pi*(col+.5)/width
                    d = direction(theta,phi)
                    weight = 2*math.pi/width*(math.cos(math.pi*row/height)-math.cos(math.pi*(row+1)/height))
                for index,y in enumerate(cartesian_basis(d)): terms[index].append(weight*y)
        expected = [[math.fsum(t)*v for v in color] for t in terms]
        c = project(image)
        for actual, ref in zip(c,expected): assert actual == pytest.approx(ref, abs=4e-14)
        for n in [(1,0,0), (0,1,0), (0,0,1)]:
            b = cartesian_basis(n)
            prediction = [math.fsum(row[ch]*y*[math.pi,2*math.pi/3,math.pi/4][math.isqrt(i)]
                                   for i,(row,y) in enumerate(zip(expected,b))) for ch in range(3)]
            actual = sh.reconstruct(sh.convolve_diffuse(c),n)
            assert actual == pytest.approx(prediction, abs=4e-14)
            # Independently predicted quadrature error, separated from
            # roundoff. A coarse image is not an exact constant integral.
            observed_error = max(abs(x-math.pi*v) for x,v in zip(actual,color))
            reference_error = max(abs(x-math.pi*v) for x,v in zip(prediction,color))
            assert observed_error == pytest.approx(reference_error, rel=2e-11, abs=1e-14)


@pytest.mark.parametrize('bad', [[], [-1.], [math.nan], [math.inf], [0.,1.]])
def test_projection_weight_validation(bad):
    sample = sh.SPHSample(0.,0.,Vector(3,[0,0,1]),1); sample.values = [1/math.sqrt(4*math.pi)]
    with pytest.raises(ValueError): sh.project_radiance([sample],[[1,2,3]],bad)


@pytest.mark.parametrize('values', [[], [1,2], [math.nan], [math.inf]])
def test_projection_basis_validation(values):
    sample = sh.SPHSample(0.,0.,Vector(3,[0,0,1]),len(values)); sample.values = values
    with pytest.raises(ValueError): sh.project_radiance([sample],[[1,2,3]])


@pytest.mark.parametrize('count', [1,4,9])
def test_orientation_actual_norm_boundary_and_temporary_normalization(count):
    c = [i+.25 for i in range(count)]
    axis = [0,0,1]; angle = 2*math.atan2(.8,.6)
    for sign in [-1,1]:
        scale = 1+sign*1e-12
        def norm(s): return math.hypot(.6*s,.8*s)
        for _ in range(20):
            if abs(norm(scale)-1) <= 1e-12: break
            scale = math.nextafter(scale,1.)
        else: raise AssertionError('boundary not found')
        rejected = math.nextafter(scale, math.inf if sign>0 else 0.)
        while abs(norm(rejected)-1) <= 1e-12:
            rejected = math.nextafter(rejected, math.inf if sign>0 else 0.)
        q = Quaternion([.6*scale,0,0,.8*scale]); saved = q.data[:]; storage = q.data
        out = sh.rotate_coefficients(c,q)
        for d in [(1,0,0),(.3,.4,math.sqrt(.75))]:
            expected = math.fsum(x*y for x,y in zip(c,cartesian_basis(rodrigues_rotate(d,axis,-angle))))
            assert math.fsum(x*y for x,y in zip(out,cartesian_basis(d))) == pytest.approx(expected, abs=1e-14)
        assert q.data is storage and q.data == saved
        with pytest.raises(ValueError): sh.rotate_coefficients(c,Quaternion([.6*rejected,0,0,.8*rejected]))


def test_constant_integrations_converge_without_forced_renormalization():
    from examples.hdr_sh import environment
    for project in [sh.project_angular_probe, environment.project_latlong]:
        errors = []
        for width in [16,32,64]:
            image = [[[1.,2.,4.] for _ in range(width)] for _ in range(width//2)]
            irradiance = sh.convolve_diffuse(project(image))
            errors.append(max(abs(x-math.pi*v) for d in [(1,0,0),(0,1,0),(0,0,1)]
                              for x,v in zip(sh.reconstruct(irradiance,d),[1,2,4])))
        assert errors == sorted(errors,reverse=True) and errors[-1] < errors[0]/8


@pytest.mark.parametrize('scale', [1e-200, 1., 1e200])
def test_projection_compensated_cancellation_and_channel_independence(scale):
    samples = [sh.SPHSample(0.,0.,Vector(3,[0,0,1]),1) for _ in range(3)]
    value = 1/math.sqrt(4*math.pi)
    for sample in samples: sample.values = [value]
    colors = [[1e16*scale, scale, 0.], [scale, 2*scale, 0.], [-1e16*scale, 4*scale, 0.]]
    expected = [value*scale, 7*value*scale, 0.]
    actual = sh.project_radiance(samples,colors,[1.,1.,1.])[0]
    assert actual == pytest.approx(expected, rel=3e-15, abs=0)
    zero = sh.project_radiance(samples,colors,[0.,0.,0.])
    assert zero == [[0.,0.,0.]]


@pytest.mark.parametrize('count', [1,4,9,16,25,49])
def test_reconstruction_channels_complete_bands_and_repeatability(count):
    c = [[(i+.25)*1e-200, (-1)**i*(i+.5)*1e200, 0.] for i in range(count)]
    saved = copy.deepcopy(c); d = [.3,.4,math.sqrt(.75)]
    basis = [spherical_reference(l,m,math.acos(d[2]),math.atan2(d[1],d[0]))
             for l in range(math.isqrt(count)) for m in range(-l,l+1)]
    expected = [math.fsum(row[k]*y for row,y in zip(c,basis)) for k in range(3)]
    for _ in range(2):
        output = sh.reconstruct(c,d)
        for k in range(3):
            budget = 2e-14*math.fsum(abs(row[k]*y) for row,y in zip(c,basis))
            assert abs(output[k]-expected[k]) <= budget
        output[0] = 999
        assert c == saved and d == [.3,.4,math.sqrt(.75)]


def test_generate_samples_independent_values_global_rng_and_ownership():
    state = random.getstate()
    try:
        random.seed(473002)
        samples = sh.GenerateSamples(6,5)
        random.seed(473002)
        again = sh.GenerateSamples(6,5)
        assert [(s.theta,s.phi,s.dir.vector,s.values) for s in samples] == [
            (s.theta,s.phi,s.dir.vector,s.values) for s in again]
        for index,sample in enumerate(samples):
            i,j = divmod(index,6)
            assert i/6 <= (1-sample.dir.vector[2])/2 < (i+1)/6
            assert j/6 <= sample.phi/(2*math.pi) < (j+1)/6
            assert sample.dir.vector == pytest.approx(direction(sample.theta,sample.phi), abs=1e-15)
            values = [spherical_reference(l,m,sample.theta,sample.phi) for l in range(5) for m in range(-l,l+1)]
            assert sample.values == pytest.approx(values, abs=4e-15)
        assert len({id(s.values) for s in samples}) == len({id(s.dir.vector) for s in samples}) == 36
        v = Vector(3,[1,0,0]); sample = sh.SPHSample(0,0,v,4)
        assert sample.dir is v
        sample.dir.vector[0] = 2
        assert v.vector == [2,0,0]  # Documented container reference ownership.
    finally:
        random.setstate(state)


def test_higher_band_rotation_and_convolution_rejected_without_basis_expansion():
    c = [[0.,0.,0.] for _ in range(16)]
    with pytest.raises(ValueError): sh.rotate_coefficients(c,Quaternion())
    with pytest.raises(ValueError): sh.convolve_diffuse(c)
    assert sh.reconstruct(c,[0,0,1]) == [0.,0.,0.]
