import math
import random
import struct
import pytest
from gem import vector
from gem.experimental import bezier, legendre, sph, sph_sample, sph_object, sph_irradiance_map


def V(*xs):
    return vector.Vector(len(xs),list(xs))


@pytest.mark.parametrize('t',[0,0.25,0.5,0.75,1])
def test_quadratic_bezier_scalar(t):
    assert bezier.quadraticBezierPoint(t,0,1,2) == pytest.approx(2*t)


@pytest.mark.parametrize('t',[0,0.25,0.5,0.75,1])
@pytest.mark.defect('E01')
def test_cubic_bezier_scalar(t):
    assert bezier.cubicBezierPoint(t,0,1,2,3) == pytest.approx(3*t)


@pytest.mark.defect('E02')
def test_bezier_vectors():
    assert bezier.quadraticBezierPoint(0.5,V(0,0),V(1,1),V(2,0)).vector == [1,0.5]


@pytest.mark.defect('E03')
def test_bezier_path_count_python3():
    path = bezier.BezierPath()
    path.setControlPoints([0,1,2,3])
    assert isinstance(path.curveCount,int)
    assert path.getDrawingPoints()


@pytest.mark.defect('E04')
def test_bezier_recursive_sampling():
    path = bezier.BezierPath()
    path.setControlPoints([V(0,0),V(1,0),V(2,0),V(3,0)])
    assert path.findDrawingPoints(0)


@pytest.mark.defect('E04')
def test_bezier_sample_points():
    path = bezier.BezierPath()
    path.samplePoints([V(0,0),V(1,0),V(2,0)],0.01,1,0.5)
    assert path.controlPoints


@pytest.mark.parametrize('l,m,x,expected',[(0,0,0.2,1),(1,0,0.2,0.2),(1,1,0,-1),(2,0,0.5,-0.125),(2,1,0.5,-1.5*math.sqrt(0.75)),(2,2,0.5,2.25)])
def test_legendre_low_orders(l,m,x,expected):
    assert legendre.Legendre(l,m,x).run() == pytest.approx(expected)


@pytest.mark.parametrize('l,x',[(3,0.2),(4,0.5),(5,-0.4)])
@pytest.mark.defect('E05')
def test_legendre_recurrence(l,x):
    a,b = 1.0,x
    for n in range(2,l+1):
        a,b = b,((2*n-1)*x*b-(n-1)*a)/n
    assert legendre.Legendre(l,0,x).run() == pytest.approx(b)


@pytest.mark.defect('E05')
def test_legendre_repeatability():
    p = legendre.Legendre(2,2,0.5)
    assert p.run() == p.run()


@pytest.mark.parametrize('n',range(8))
def test_factorial(n):
    assert sph.Factorial(n) == math.factorial(n)


def test_spherical_harmonics_low_orders():
    assert sph.SPH(0,0,0.2,0.3) == pytest.approx(1/math.sqrt(4*math.pi))
    assert sph.SPH(1,0,0,0) == pytest.approx(math.sqrt(3/(4*math.pi)))
    assert sph.SPH(1,1,math.pi/2,0) == pytest.approx(-math.sqrt(3/(4*math.pi)))
    assert sph.SPH(1,-1,math.pi/2,math.pi/2) == pytest.approx(-math.sqrt(3/(4*math.pi)))


@pytest.mark.defect('E06')
def test_generate_samples(monkeypatch):
    monkeypatch.setattr(random,'random',lambda:0.5)
    samples = sph_sample.GenerateSamples(2,3)
    assert len(samples) == 4
    for sample in samples:
        assert sample.dir.magnitude() == pytest.approx(1)
        assert len(sample.values) == 9
        assert sample.values[0] == pytest.approx(1/math.sqrt(4*math.pi))


@pytest.mark.defect('E07')
def test_generate_object_coefficients():
    vertex = sph_object.SPHVertex(V(0,0,0),V(0,0,1))
    sample = sph_sample.SPHSample(0,0,V(0,0,1),1)
    sample.values[0] = 1/math.sqrt(4*math.pi)
    obj = sph_object.SPHObject([0],[vertex])
    sph_object.GenereateCoeffs(1,1,[sample],[obj])
    assert vertex.unshadowedCoeffs == pytest.approx([math.sqrt(4*math.pi)])


@pytest.mark.defect('E08')
def test_rectangular_irradiance_file(tmp_path):
    path = tmp_path/'probe.float'
    path.write_bytes(struct.pack('12f',*([1.0]*12)))
    result = sph_irradiance_map.SPH_IrradianceMapCoeff(str(path),2,1)
    assert len(result.hdr) == 1
    assert len(result.hdr[0]) == 2


def test_irradiance_coefficient_known_answer(tmp_path):
    path = tmp_path/'probe.float'
    path.write_bytes(struct.pack('12f',*([1.0]*12)))
    result = sph_irradiance_map.SPH_IrradianceMapCoeff(str(path),2,2)
    result.coeffs = [[0.0]*3 for _ in range(9)]
    result.updateCoefficients([1,2,3],1,0,0,1)
    assert result.coeffs[0] == pytest.approx([0.282095,0.564190,0.846285])
    assert result.coeffs[2] == pytest.approx([0.488603,0.977206,1.465809])
    assert result.coeffs[6] == pytest.approx([0.630784,1.261568,1.892352])
    assert result.coeffs[8] == [0,0,0]


@pytest.mark.parametrize('l', [0,1,2,pytest.param(3,marks=pytest.mark.defect('E05')),pytest.param(4,marks=pytest.mark.defect('E05'))])
def test_spherical_harmonics_addition_theorem(l):
    # Real orthonormal harmonics: sum_m Y_lm(theta,phi)^2 = (2l+1)/(4*pi).
    theta,phi = 0.73,1.27
    assert sum(sph.SPH(l,m,theta,phi)**2 for m in range(-l,l+1)) == pytest.approx((2*l+1)/(4*math.pi))


def test_spherical_harmonics_orthonormality():
    # Midpoint quadrature in z=cos(theta),phi; deterministic and pure Python.
    pairs = [(0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),(2,0),(2,1),(2,2)]
    gram = [[0.0]*len(pairs) for _ in pairs]
    nz,np = 120,24
    weight = 2/nz*2*math.pi/np
    for i in range(nz):
        theta = math.acos(-1+(i+0.5)*2/nz)
        for j in range(np):
            phi = (j+0.5)*2*math.pi/np
            values = [sph.SPH(l,m,theta,phi) for l,m in pairs]
            for a in range(len(pairs)):
                for b in range(a+1):
                    gram[a][b] += values[a]*values[b]*weight
    for a in range(len(pairs)):
        for b in range(a+1):
            assert gram[a][b] == pytest.approx(int(a==b),abs=4e-4)
