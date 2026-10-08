import math
import pytest
from gem.quaternion import Quaternion
from gem import spherical_harmonics as sh


def q(axis,angle):
    norm=math.hypot(*axis);s=math.sin(angle/2)/norm
    return Quaternion([math.cos(angle/2)]+[v*s for v in axis])


def rodrigues(v,axis,angle):
    norm=math.hypot(*axis);a=[x/norm for x in axis];c=math.cos(angle);s=math.sin(angle)
    cross=[a[1]*v[2]-a[2]*v[1],a[2]*v[0]-a[0]*v[2],a[0]*v[1]-a[1]*v[0]]
    dot=sum(x*y for x,y in zip(a,v))
    return [c*x+s*y+(1-c)*dot*z for x,y,z in zip(v,cross,a)]


def basis(d):
    x,y,z=d;a=math.sqrt(3/(4*math.pi));b=math.sqrt(15/(4*math.pi))
    return [1/math.sqrt(4*math.pi),-a*y,a*z,-a*x,b*x*y,-b*y*z,
            math.sqrt(5/(16*math.pi))*(3*z*z-1),-b*x*z,math.sqrt(15/(16*math.pi))*(x*x-y*y)]


@pytest.mark.parametrize('axis',[(1,0,0),(0,1,0),(0,0,1),(1,2,-3)])
@pytest.mark.parametrize('angle',[math.pi/2,-math.pi/2,math.pi,.73])
@pytest.mark.parametrize('count',[1,4,9])
def test_directional_identity_energy_and_inverse(axis,angle,count):
    c=[(-1)**i*(i+.3) for i in range(count)];orientation=q(axis,angle)
    saved=orientation.data[:];storage=orientation.data
    out=sh.rotate_coefficients(c,orientation)
    for direction in [(1,0,0),(0,1,0),(0,0,1),(.3,.4,math.sqrt(.75))]:
        original_direction=rodrigues(direction,axis,-angle)
        expected=sum(x*y for x,y in zip(c,basis(original_direction)))
        assert sum(x*y for x,y in zip(out,basis(direction)))==pytest.approx(expected,abs=3e-14)
    for start,end in [(0,1),(1,4),(4,9)]:
        assert sum(x*x for x in out[start:end])==pytest.approx(sum(x*x for x in c[start:end]),rel=2e-14)
    assert out[0]==c[0]
    assert sh.rotate_coefficients(out,q(axis,-angle))==pytest.approx(c,abs=2e-14)
    assert orientation.data is storage and orientation.data==saved
    assert out is not c


def test_active_z_axis_feature_and_l2_known_answer():
    c=[0.]*9;c[3]=-1;c[8]=1
    out=sh.rotate_coefficients(c,q((0,0,1),math.pi/2))
    expected=[0.]*9;expected[1]=-1;expected[8]=-1
    assert out==pytest.approx(expected,abs=1e-15)


def test_rgb_identity_composition_and_diffuse():
    c=[[i+.2,(-1)**i*(i+1),.5*i] for i in range(9)];saved=[r[:] for r in c]
    a=q((1,0,0),.7);b=q((0,1,0),-.9)
    sequential=sh.rotate_coefficients(sh.rotate_coefficients(c,a),b)
    combined=sh.rotate_coefficients(c,b*a)
    reverse=sh.rotate_coefficients(c,a*b)
    assert any(abs(x-y)>1e-3 for r,s in zip(combined,reverse) for x,y in zip(r,s))
    for x,y in zip(sequential,combined):assert x==pytest.approx(y,abs=2e-14)
    identity=sh.rotate_coefficients(c,Quaternion())
    for x,y in zip(c,identity):assert x==pytest.approx(y,abs=1e-15);assert x is not y
    for channel in range(3):
        scalar=sh.rotate_coefficients([r[channel] for r in c],a)
        assert [r[channel] for r in sh.rotate_coefficients(c,a)]==pytest.approx(scalar)
    left=sh.rotate_coefficients(sh.convolve_diffuse(c),a)
    right=sh.convolve_diffuse(sh.rotate_coefficients(c,a))
    for x,y in zip(left,right):assert x==pytest.approx(y,abs=1e-14)
    assert c==saved
    negative=Quaternion([-v for v in a.data])
    for x,y in zip(sh.rotate_coefficients(c,a),sh.rotate_coefficients(c,negative)):assert x==pytest.approx(y)


def test_norm_boundary_adjacent_binary64_values():
    # 1+1e-12 rounds just outside the tolerance. Locate adjacent floats on
    # each side of the actual representable boundary, rather than assume it.
    upper=1+1e-12
    while upper-1>1e-12:upper=math.nextafter(upper,1)
    lower=1-1e-12
    while 1-lower>1e-12:lower=math.nextafter(lower,1)
    for accepted,rejected in [(upper,math.nextafter(upper,math.inf)),(lower,math.nextafter(lower,0))]:
        assert abs(accepted-1)<=1e-12<abs(rejected-1)
        assert sh.rotate_coefficients([1,2,3,4],Quaternion([accepted,0,0,0]))==[1,2,3,4]
        with pytest.raises(ValueError):sh.rotate_coefficients([1],Quaternion([rejected,0,0,0]))


@pytest.mark.parametrize('orientation',[None,[1,0,0,0],Quaternion([0,0,0,0]),Quaternion([2,0,0,0]),Quaternion([math.nan,0,0,0]),Quaternion([math.inf,0,0,0]),Quaternion([1,0,0])])
def test_invalid_orientations(orientation):
    with pytest.raises(ValueError):sh.rotate_coefficients([1],orientation)


@pytest.mark.parametrize('coefficients',[[],[1,2],[[1,2]],[[1,2,3],1,1,1],[math.nan],[math.inf]])
def test_invalid_coefficients(coefficients):
    with pytest.raises(ValueError):sh.rotate_coefficients(coefficients,Quaternion())


def test_legacy_conversion_is_explicit():
    legacy=[[i+1.,2*i+1.,3*i+1.] for i in range(9)]
    canonical=sh.legacy_to_canonical(legacy);orientation=q((1,2,3),.8)
    rotated=sh.rotate_coefficients(canonical,orientation)
    d=[.3,.4,math.sqrt(.75)];inverse=rodrigues(d,(1,2,3),-.8)
    assert sh.reconstruct(rotated,d)==pytest.approx(sh.reconstruct(canonical,inverse),abs=3e-14)
    assert legacy==[[i+1.,2*i+1.,3*i+1.] for i in range(9)]
