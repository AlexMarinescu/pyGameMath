"""Independent near-pole references: Cartesian polynomials and Rodrigues.

Error budgets count binary64 ulps, including unavoidable subnormal rounding.
No core SH or quaternion/matrix rotation is used to derive expected values.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
import math

import pytest

from gem import spherical_harmonics as sh
from gem.legendre import Legendre
from gem.quaternion import Quaternion
from gem.vector import Vector
from tests.test_spherical_harmonics_audit import rodrigues_coefficients, normalization_reference


def decimal_basis(d):
    with localcontext() as context:
        context.prec = 180
        x,y,z = map(Decimal.from_float, map(float,d))
        pi = Decimal.from_float(math.pi)
        a = (3/(4*pi)).sqrt(); b = (15/(4*pi)).sqrt()
        values = [1/(4*pi).sqrt(),-a*y,a*z,-a*x,b*x*y,-b*y*z,
                  (5/(16*pi)).sqrt()*(3*z*z-1),-b*x*z,(15/(16*pi)).sqrt()*(x*x-y*y)]
        return [float(v) for v in values]


def cosine_rounding_boundary():
    lo,hi = 0.,1e-7
    for _ in range(100):
        mid = (lo+hi)/2
        if mid == lo or mid == hi: break
        if math.cos(mid) == 1: lo = mid
        else: hi = mid
    assert math.nextafter(lo,math.inf) == hi
    assert math.cos(lo) == 1 > math.cos(hi)
    return lo,hi


BOUNDARY = cosine_rounding_boundary()
ANGLES = [0., 2**-10, 2**-20, *BOUNDARY, math.nextafter(BOUNDARY[1],math.inf),
          1e-9, 1e-12, 1e-100, 1e-150, 1e-160, 1e-308, 1e-320,
          20*math.ulp(0.), math.nextafter(math.pi,0.), math.pi-1e-9, math.pi]


@pytest.mark.parametrize('theta', ANGLES)
@pytest.mark.parametrize('phi', [0., .3, -.7])
def test_low_band_angle_boundaries_and_polar_endpoint(theta,phi):
    sine = 0. if theta in (0.,math.pi) else math.sin(theta)
    d = [sine*math.cos(phi),sine*math.sin(phi),math.cos(theta)]
    expected = decimal_basis(d)
    actual = [sh.SPH(l,m,theta,phi) for l in range(3) for m in range(-l,l+1)]
    for a,b in zip(actual,expected):
        assert abs(a-b) <= 16*math.ulp(b)
    assert sh._basis(3,theta,phi) == actual
    # Cached normalization never owns or reuses a sampled basis array.
    result = sh._basis(3,theta,phi); result[0] = 999
    assert sh._basis(3,theta,phi) == actual


@pytest.mark.parametrize('t', [1e-7, 1e-9, 1e-50, 1e-150, 1e-300, 1e-320])
@pytest.mark.parametrize('pole', [-1.,1.])
@pytest.mark.parametrize('xy', [(1.,0.),(-1.,0.),(0.,1.),(0.,-1.),(.8,-.6)])
def test_reconstruction_signed_transverse_components_and_ownership(t,pole,xy):
    x,y = t*xy[0],t*xy[1]
    d = [x,y,pole*math.sqrt(1-x*x-y*y)]
    assert math.hypot(*d) == pytest.approx(1.,abs=2e-16)
    expected = decimal_basis(d)
    position = Vector(3,d[:]); storage = position.vector
    for index in range(9):
        coefficients = [[0.]*3 for _ in range(9)]
        coefficients[index] = [1.,2.,-4.]
        original = [row[:] for row in coefficients]
        out = sh.reconstruct(coefficients,position)
        # Form channel expectations before rounding to retain subnormal terms.
        with localcontext() as context:
            context.prec = 180
            reference = [float(Decimal.from_float(expected[index])*Decimal(v)) for v in [1,2,-4]]
        for actual,target in zip(out,reference): assert abs(actual-target) <= 20*math.ulp(target)
        assert sh.reconstruct(coefficients,d) == out
        out[0] = 999
        assert coefficients == original
    assert position.vector is storage and position.vector == d


@pytest.mark.parametrize('xy', [(1e-9,1e-100), (-1e-100,1e-9), (1e-300,-1e-300),
    (math.ulp(0.),0.), (0.,-20*math.ulp(0.)), (40*math.ulp(0.),-30*math.ulp(0.))])
@pytest.mark.parametrize('pole', [-1.,1.])
def test_mixed_scale_and_subnormal_transverse_directions(xy,pole):
    d = [*xy,pole]
    assert math.hypot(*d) == 1.
    expected = decimal_basis(d)
    for index in [1,3,4,5,7,8]:
        c = [[0.]*3 for _ in range(9)]; c[index] = [1.,0.,0.]
        actual = sh.reconstruct(c,d)
        assert abs(actual[0]-expected[index]) <= 16*math.ulp(expected[index])
        assert actual[1:] == [0.,0.]


@pytest.mark.parametrize('theta', [1e-9, math.pi-1e-9])
@pytest.mark.parametrize('degree', [3,8,12,32])
def test_shared_associated_seed_with_unchanged_higher_degree_recurrence(theta,degree):
    x = Fraction(math.cos(theta)); sine = math.sin(theta); phi = .73
    for m in [1,2,3]:
        derivative = sum(c*x**p for p,c in rodrigues_coefficients(degree,m))
        for order in [m,-m]:
            with localcontext() as context:
                context.prec = 180
                value = Decimal(derivative.numerator)/Decimal(derivative.denominator)
                value *= Decimal.from_float(sine)**m*Decimal((-1)**m)*Decimal(2).sqrt()
                value *= normalization_reference(degree,m)
                value *= Decimal.from_float(math.cos(m*phi) if order>0 else math.sin(m*phi))
                expected = float(value)
            assert sh.SPH(degree,order,theta,phi) == pytest.approx(expected,rel=2e-14,abs=0)


@pytest.mark.parametrize('theta', [0.,math.pi])
def test_exact_poles_and_public_legendre_scratch_contract(theta):
    for degree in range(13):
        for order in range(-degree,degree+1):
            actual = sh.SPH(degree,order,theta,.73)
            if order: assert actual == 0.
            else:
                sign = 1 if theta == 0. else (-1)**degree
                expected = sign*math.sqrt((2*degree+1)/(4*math.pi))
                assert abs(actual-expected) <= 2*math.ulp(expected)
    p = Legendre(8,3,.25)
    p.calculatePML(8)
    saved = p.__dict__.copy()
    for _ in range(2): sh.SPH(8,3,theta,.73)
    assert p.__dict__ == saved
    assert p.run() == p.PML


@pytest.mark.parametrize('theta', [Fraction(1,2),Decimal('.5'),0,1])
def test_existing_numeric_angle_representations(theta):
    assert sh._basis(3,theta,.7) == sh._basis(3,float(theta),.7)
    assert sh.SPH(2,-1,theta,.7) == sh.SPH(2,-1,float(theta),.7)


def test_near_pole_rotation_convolution_projection_and_channel_scaling():
    d = [1e-9,-2e-9,1.]
    c = [[0.]*3 for _ in range(9)]; c[3] = [1e6,2e6,0.]
    original = [row[:] for row in c]
    q = Quaternion([math.sqrt(.5),0,0,math.sqrt(.5)])
    rotated = sh.rotate_coefficients(c,q)
    a = math.sqrt(3/(4*math.pi))
    expected = [-a*d[1]*1e6,-a*d[1]*2e6,0.]
    assert sh.reconstruct(rotated,d) == pytest.approx(expected,rel=2e-15,abs=1e-18)
    irradiance = sh.convolve_diffuse(rotated)
    assert sh.reconstruct(irradiance,d) == pytest.approx([x*2*math.pi/3 for x in expected],rel=2e-15,abs=1e-18)
    sample = sh.SPHSample(1e-9,.3,Vector(3,[1e-9,0,1]),9)
    sample.values = sh._basis(3,1e-9,.3)
    projected = sh.project_radiance([sample],[[1.,2.,0.]],[.5])
    sine = math.sin(1e-9)
    basis = decimal_basis([sine*math.cos(.3),sine*math.sin(.3),math.cos(1e-9)])
    for actual,b in zip(projected,basis): assert actual == pytest.approx([b*.5,b,0.],rel=3e-15,abs=0)
    assert c == original
