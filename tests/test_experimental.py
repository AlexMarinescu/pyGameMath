import math
import struct
import pytest
from gem import vector
from gem import spherical_harmonics as sph
from gem import spherical_harmonics as sph_irradiance_map


def V(*xs):
    return vector.Vector(len(xs),list(xs))


@pytest.mark.parametrize('n',range(8))
def test_factorial(n):
    assert sph.Factorial(n) == math.factorial(n)


def test_spherical_harmonics_low_orders():
    assert sph.SPH(0,0,0.2,0.3) == pytest.approx(1/math.sqrt(4*math.pi))
    assert sph.SPH(1,0,0,0) == pytest.approx(math.sqrt(3/(4*math.pi)))
    assert sph.SPH(1,1,math.pi/2,0) == pytest.approx(-math.sqrt(3/(4*math.pi)))
    assert sph.SPH(1,-1,math.pi/2,math.pi/2) == pytest.approx(-math.sqrt(3/(4*math.pi)))



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


@pytest.mark.parametrize('l', [0,1,2,3,4])
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
