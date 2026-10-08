import math
from fractions import Fraction
import pytest
from gem import legendre

@pytest.mark.parametrize('l,m,x,expected',[(0,0,0.2,1),(1,0,0.2,0.2),(1,1,0,-1),(2,0,0.5,-0.125),(2,1,0.5,-1.5*math.sqrt(0.75)),(2,2,0.5,2.25)])
def test_legendre_low_orders(l,m,x,expected):
    assert legendre.Legendre(l,m,x).run() == pytest.approx(expected)


@pytest.mark.parametrize('l,x',[(3,0.2),(4,0.5),(5,-0.4)])
def test_legendre_recurrence(l,x):
    a,b = 1.0,x
    for n in range(2,l+1):
        a,b = b,((2*n-1)*x*b-(n-1)*a)/n
    assert legendre.Legendre(l,0,x).run() == pytest.approx(b)


def test_legendre_repeatability():
    p = legendre.Legendre(2,2,0.5)
    assert p.run() == p.run()




def rodrigues(l,m,x):
    # Differentiate the explicit Rodrigues polynomial, using exact rational
    # coefficients. This reference does not use the three-term recurrence.
    value=Fraction(0)
    for k in range((l-m)//2+1):
        power=l-2*k
        coefficient=Fraction((-1)**k*math.factorial(2*l-2*k),
                             2**l*math.factorial(k)*math.factorial(l-k)*math.factorial(power-m))
        value += coefficient*x**(power-m)
    return (-1)**m * float(value) * (1-float(x)**2)**(m/2)


@pytest.mark.parametrize('l',range(13))
@pytest.mark.parametrize('x',[Fraction(-1),Fraction(-3,5),Fraction(0),Fraction(1,5),Fraction(3,5),Fraction(1)])
def test_all_orders_independent_rodrigues(l,x):
    for m in range(l+1):
        expected=rodrigues(l,m,x)
        assert legendre.Legendre(l,m,float(x)).run()==pytest.approx(expected,rel=3e-12,abs=1e-9)


@pytest.mark.parametrize('l',range(1,13))
def test_recurrence_parity_and_boundary(l):
    for m in range(l):
        x=.37
        previous=legendre.Legendre(l-1,m,x).run()
        current=legendre.Legendre(l,m,x).run()
        following=legendre.Legendre(l+1,m,x).run()
        assert (l-m+1)*following==pytest.approx((2*l+1)*x*current-(l+m)*previous,rel=3e-14,abs=1e-8)
        assert legendre.Legendre(l,m,-x).run()==pytest.approx((-1)**(l+m)*current)
    assert legendre.Legendre(l,0,1).run()==1
    assert legendre.Legendre(l,0,-1).run()==(-1)**l


@pytest.mark.parametrize('l,m',[(0,0),(1,0),(4,2),(12,7),(12,12)])
def test_run_preserves_all_state_and_helper_initialization(l,m):
    p=legendre.Legendre(l,m,.2)
    p.P,p.PM1,p.PML=123.,456.,789.
    saved=p.__dict__.copy()
    expected=rodrigues(l,m,Fraction(1,5))
    for _ in range(3):
        assert p.run()==pytest.approx(expected)
        assert p.__dict__==saved
    p.mGreaterThan0(); first=p.P
    p.mGreaterThan0(); assert p.P==first
    p.calculatePM1(); first=p.PM1
    p.calculatePM1(); assert p.PM1==first
    p.calculatePML(max(l,m+2)); first=p.PML
    p.calculatePML(max(l,m+2)); assert p.PML==first
    assert p.PML==pytest.approx(rodrigues(max(l,m+2),m,Fraction(1,5)))


@pytest.mark.parametrize('l',range(13))
@pytest.mark.parametrize('theta,phi',[(.2,-.3),(.73,1.27),(2.2,3.1),(0,0),(math.pi,.9)])
def test_spherical_addition_theorem_high_orders(l,theta,phi):
    from gem.experimental import sph
    assert sum(sph.SPH(l,m,theta,phi)**2 for m in range(-l,l+1))==pytest.approx((2*l+1)/(4*math.pi),rel=2e-13)


def test_phase_normalization_and_compatibility():
    from gem.experimental.legendre import Legendre
    from gem.experimental import sph
    assert Legendre is legendre.Legendre is sph.Legendre
    assert Legendre(1,1,0).run()==-1
    assert Legendre(2,2,0).run()==3
    assert sph.SPH(1,1,math.pi/2,0)==pytest.approx(-math.sqrt(3/(4*math.pi)))
    assert sph.SPH(1,-1,math.pi/2,math.pi/2)==pytest.approx(-math.sqrt(3/(4*math.pi)))
    assert Legendre(3,0,2).run()==17 # ordinary polynomial outside [-1,1]
