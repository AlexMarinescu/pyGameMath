"""Independent angle, tangent-plane and ownership interpolation checks."""
import math
import itertools
import pytest
from gem import quaternion as q


def axis_quat(degrees):
    angle = math.radians(degrees)/2
    return q.Quaternion([math.cos(angle), 0, 0, math.sin(angle)])


def norm(values):
    return math.sqrt(sum(c*c for c in values))


@pytest.mark.parametrize('degrees', [0, 0.001, 1, 5.125, 5.126, 90, 180])
@pytest.mark.parametrize('t', [0, .13, .5, .87, 1])
def test_slerp_known_axis(degrees, t):
    result = axis_quat(0).slerp(axis_quat(degrees),t)
    assert result.data == pytest.approx(axis_quat(degrees*t).data, abs=1e-12, rel=0)
    assert norm(result.data) == pytest.approx(1, abs=1e-12, rel=0)


@pytest.mark.parametrize('t', [0, .2, .5, .8, 1])
def test_slerp_tangent_plane_reference(t):
    # Orthogonal unit quaternions: exact interpolation on their great circle.
    a=[.5,.5,.5,.5]; b=[.5,-.5,.5,-.5]
    expected=[x*math.cos(t*math.pi/2)+y*math.sin(t*math.pi/2) for x,y in zip(a,b)]
    assert q.quat_slerp(q.Quaternion(a),q.Quaternion(b),t).data == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize('t', [0.0, .5, .8, 1.0])
def test_legacy_squad_known_axis(t):
    # All three nested blends use spherical branches at these interior t.
    # q0=0, q1=180, q2=90 spatial degrees.
    expected_degrees=90*t + 180*t*t*(1-t)
    assert axis_quat(0).squad(axis_quat(180),axis_quat(90),t).data == pytest.approx(axis_quat(expected_degrees).data,abs=1e-12)


@pytest.mark.parametrize('degrees', [0, 30, 180])
def test_legacy_identical_and_endpoints(degrees):
    a=axis_quat(degrees); b=axis_quat(70); c=axis_quat(100)
    assert a.squad(b,c,0.0).data == pytest.approx(a.data)
    assert a.squad(b,c,1.0).data == pytest.approx(c.data)
    assert a.squad(a,a,.37).data == pytest.approx(a.data)


def test_legacy_linear_branch_and_sign_policy():
    a=axis_quat(0); b=axis_quat(1); c=axis_quat(2)
    # At t=.5, the retained linear branches give .5*a+.25*b+.25*c.
    expected=[.5*x+.25*y+.25*z for x,y,z in zip(a.data,b.data,c.data)]
    assert a.squad(b,c,.5).data == pytest.approx(expected)
    assert norm(expected) < 1
    out=a.squad(b,c,.5)
    assert a.negate().squad(b.negate(),c.negate(),.5).data == pytest.approx([-x for x in out.data])
    assert a.squad(b.negate(),c,.5).data != pytest.approx(out.data)


@pytest.mark.parametrize('t', [0,.1,.25,.5,.75,.9,1])
def test_squad4_independent_axis_reference(t):
    # Endpoint blend angle=90t; control blend angle=30+90t.
    # Final spatial angle=90t+60t(1-t), with no shortest-path wrap.
    result=q.squad4(axis_quat(0),axis_quat(90),axis_quat(30),axis_quat(120),t)
    assert result.data == pytest.approx(axis_quat(90*t+60*t*(1-t)).data,abs=1e-12,rel=0)
    assert norm(result.data) == pytest.approx(1,abs=1e-12,rel=0)


def test_squad4_midpoint_and_control_influence():
    a=axis_quat(0); b=axis_quat(90)
    out=q.squad4(a,b,axis_quat(30),axis_quat(120),.5)
    assert out.data == pytest.approx([math.sqrt(3)/2,0,0,.5],abs=1e-12)
    assert out.data != pytest.approx(q.squad4(a,b,a,b,.5).data)


@pytest.mark.parametrize('signs', list(itertools.product([-1,1],repeat=4)))
@pytest.mark.parametrize('t', [0,.23,.5,1])
def test_squad4_sign_equivalence(signs,t):
    controls=[axis_quat(d) for d in [0,90,30,120]]
    signed=[q.Quaternion([sign*c for c in quat.data]) for sign,quat in zip(signs,controls)]
    out=q.squad4(*signed,t)
    expected=axis_quat(90*t+60*t*(1-t))
    assert abs(sum(x*y for x,y in zip(out.data,expected.data))) == pytest.approx(1,abs=1e-12,rel=0)


@pytest.mark.parametrize('data', [[1,0,0,0],[.5,.5,.5,.5],[-1,0,0,0]])
def test_squad4_identical(data):
    a=q.Quaternion(data)
    assert q.squad4(a,a,a,a,.43).data == pytest.approx(data)


@pytest.mark.parametrize('function', [q.quat_slerp,q.quat_squad,'squad4'])
def test_storage_preservation(function):
    controls=[axis_quat(d) for d in [0,90,30,120]]
    storage=[c.data for c in controls]; snapshots=[v[:] for v in storage]
    if function == 'squad4':
        result=q.squad4(*controls,.37)
    else:
        result=function(*controls[:2 if function is q.quat_slerp else 3],.37)
    assert isinstance(result,q.Quaternion)
    assert all(result is not c and result.data is not c.data for c in controls)
    result.data[0]=123
    assert all(c.data is s and s==v for c,s,v in zip(controls,storage,snapshots))


@pytest.mark.parametrize('seed', range(10))
def test_unit_length_and_continuity_deterministic(seed):
    import random
    rng=random.Random(seed)
    controls=[]
    for _ in range(4):
        values=[rng.uniform(-1,1) for _ in range(4)]
        length=norm(values)
        controls.append(q.Quaternion([c/length for c in values]))
    for i in range(101):
        t=i/100
        out=q.squad4(*controls,t)
        assert norm(out.data) == pytest.approx(1,abs=1e-12,rel=0)
    # Local continuity away from possible shortest-path branch ties.
    a=q.squad4(*controls,.371);b=q.squad4(*controls,.37100001)
    assert abs(sum(x*y for x,y in zip(a.data,b.data))) > 1-1e-12


@pytest.mark.parametrize('sign', [-1,1])
def test_slerp_equal_and_sign_equivalent_limits(sign):
    a=q.Quaternion([.5,.5,.5,.5]);b=q.Quaternion([sign*x for x in a.data])
    assert a.slerp(b,.3).data == pytest.approx(a.data,abs=1e-12)


def test_slerp_tiny_angle_and_extrapolation():
    a=q.Quaternion();b=q.Quaternion([1,1e-200,0,0])
    assert a.slerp(b,.5).data[1] == pytest.approx(5e-201,rel=1e-12,abs=0)
    assert axis_quat(0).slerp(axis_quat(90),1.5).data == pytest.approx(axis_quat(135).data)


def test_squad4_orthogonal_midpoint_reference():
    # Two orthogonal great-circle midpoints are themselves orthogonal.
    controls=[q.Quaternion(v) for v in [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]]
    assert q.squad4(*controls,.5).data == pytest.approx([.5,.5,.5,.5],abs=1e-12,rel=0)


def test_shortest_path_half_turn_tie_retains_supplied_sign():
    a=q.Quaternion();b=q.Quaternion([0,1,0,0])
    assert a.slerp(b,.5).data == pytest.approx([math.sqrt(.5),math.sqrt(.5),0,0])
    assert a.slerp(b.negate(),.5).data == pytest.approx([math.sqrt(.5),-math.sqrt(.5),0,0])
