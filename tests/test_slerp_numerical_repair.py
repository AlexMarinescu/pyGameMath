"""Independent endpoint, binary64 boundary and ordinary rotation references."""
from fractions import Fraction as F
import math
import random
import struct

import pytest

from gem import quaternion as q
from gem.vector import Vector

U = math.ulp(0.)
MIN_NORMAL = float.fromhex('0x1p-1022')


def bits(data):
    return [struct.pack('!d',x) for x in data]


def exact_rotation(data, point):
    # Expanded Hamilton sandwich using exact represented inputs, not gem.
    w,x,y,z = map(F,data)
    rows = [[w*w+x*x-y*y-z*z,2*(x*y+w*z),2*(x*z-w*y)],
            [2*(x*y-w*z),w*w-x*x+y*y-z*z,2*(y*z+w*x)],
            [2*(x*z+w*y),2*(y*z-w*x),w*w-x*x-y*y+z*z]]
    return [float(sum((F(point[i])*rows[i][j] for i in range(3)),F())) for j in range(3)]


@pytest.mark.parametrize('data', [
    [1.,U,0.,-0.], [1.,0.,-U,U], [0.,1.,U,-0.],
    [.5,-.5,.5,.5], [-.5,.5,.5,.5], [1.,1e-200,-1e-200,0.]])
@pytest.mark.parametrize('t', [0.,1.])
@pytest.mark.parametrize('sign', [-1.,1.])
def test_exact_endpoint_components_storage_and_independent_rotations(data,t,sign):
    selected = [sign*x for x in data]
    other = [.5,.5,-.5,.5]
    a,b = (selected,other) if t == 0. else (other,selected)
    controls = [q.Quaternion(values[:]) for values in (a,b,[0.,1.,0.,0.],[0.,0.,1.,0.])]
    storages = [c.data for c in controls]
    saved = [bits(c.data) for c in controls]
    flip = sum(F(x)*F(y) for x,y in zip(a,b)) < 0
    expected = a if t == 0. else [-x for x in b] if flip else b
    for out in (q.quat_slerp(*controls[:2],t),controls[0].slerp(controls[1],t),q.squad4(*controls,t)):
        assert type(out) is q.Quaternion
        assert bits(out.data) == bits(expected)
        assert all(out is not c and out.data is not c.data for c in controls)
        # No absolute tolerance may swallow the representable transverse value.
        actual = q.quat_rotate_vector(out,Vector(3,[0.,1.,0.])).vector
        reference = exact_rotation(expected,[0.,1.,0.])
        for x,y in zip(actual,reference):
            assert abs(x-y) <= 2*math.ulp(y)
        out.data[0] = 17.
    assert all(c.data is s and bits(s) == old for c,s,old in zip(controls,storages,saved))


@pytest.mark.parametrize('size', [U,2*U,3*U,5*U,7*U,16*U,
    math.nextafter(MIN_NORMAL,0.),MIN_NORMAL,math.nextafter(MIN_NORMAL,math.inf),1e-300,1e-200])
@pytest.mark.parametrize('sign', [-1.,1.])
def test_near_zero_angle_weights_against_exact_linear_limit(size,sign):
    # atan(size)=size+O(size**3), and sin(t*atan(size))=t*size+O(size**3).
    # Every omitted correction here is below the minimum subnormal grid.
    endpoint = q.Quaternion([1.,sign*size,0.,0.])
    for t in (0.,.125,.375,.5,.75,math.nextafter(1.,0.),1.):
        out = q.quat_slerp(q.Quaternion(),endpoint,t)
        expected = float(F(t)*F(sign*size))
        if size < MIN_NORMAL/2:
            assert out.data[1] == expected
        else:
            # Normal-angle sine divisions and the final product round separately.
            assert abs(out.data[1]-expected) <= 2*math.ulp(expected)
        assert out.data[2:] == [0.,0.]
        assert abs(out.data[0]-1.) <= 2*math.ulp(1.)
        assert abs(math.hypot(*out.data)-1.) <= 2*math.ulp(1.)


@pytest.mark.parametrize('first,last', [(U,3*U),(-U,3*U),(5*U,-3*U)])
@pytest.mark.parametrize('t', [.125,.375,.5,.75])
def test_nonidentity_controls_with_subnormal_transverse_axes(first,last,t):
    a,b = q.Quaternion([0.,1.,first,0.]),q.Quaternion([0.,1.,last,0.])
    out = q.quat_slerp(a,b,t)
    expected = float((1-F(t))*F(first)+F(t)*F(last))
    # Two weighted products can round individually before their sum, at most
    # one minimum-subnormal ulp in these exact dyadic configurations.
    assert abs(out.data[2]-expected) <= U
    assert out.data[:2] == [0.,1.] and out.data[3] == 0.
    assert a.data == [0.,1.,first,0.] and b.data == [0.,1.,last,0.]


def test_unrepresentable_midpoint_is_not_promoted_to_a_false_nonzero_value():
    assert q.quat_slerp(q.Quaternion(),q.Quaternion([1.,U,0.,0.]),.5).data == [1.,0.,0.,0.]
    assert q.quat_slerp(q.Quaternion(),q.Quaternion([1.,U,0.,0.]),.75).data == [1.,U,0.,0.]
    assert q.quat_slerp(q.Quaternion(),q.Quaternion([1.,3*U,0.,0.]),.375).data == [1.,U,0.,0.]


@pytest.mark.parametrize('t',[F(1,8),F(3,8),F(3,4)])
def test_fraction_parameters_keep_the_slerp_float_weight_dispatch(t):
    # The spherical path produces float weights even for a Fraction parameter.
    # Quaternion component LERP has narrower operand dispatch; it is not an
    # interchangeable fallback for these already accepted SLERP parameters.
    for size in (U,3*U,64*U,1e-200):
        out = q.quat_slerp(q.Quaternion(),q.Quaternion([1.,size,0.,0.]),t)
        expected = float(t*F(size))
        assert abs(out.data[1]-expected) <= math.ulp(expected)
        assert abs(out.data[0]-1.) <= math.ulp(1.)


@pytest.mark.parametrize('signs', [(1,1),(-1,1),(1,-1),(-1,-1)])
def test_subnormal_squad_controls_use_independent_dyadic_angle_polynomial(signs):
    # Endpoint angle (quaternion half-angle) is 64*u*t; control angle is
    # (32+64*t)*u. At these dyadic t every weighted term is representable.
    controls = [q.Quaternion([sign,sign*n*U,0.,0.])
                for n,sign in zip((0,64,32,96),(signs[0],signs[1],signs[0],signs[1]))]
    saved = [c.data[:] for c in controls]
    for t in (0.,.25,.5,.75,1.):
        out = q.squad4(*controls,t)
        reference = float((64*F(t)+64*F(t)*(1-F(t)))*F(U))
        assert out.data == [float(signs[0]),signs[0]*reference,0.,0.]
        assert all(out.data is not c.data for c in controls)
    assert [c.data for c in controls] == saved


def test_squad_subnormal_products_have_a_derived_rounding_budget():
    controls = [q.Quaternion([1.,n*U,0.,0.]) for n in (0,8,4,12)]
    out = q.squad4(*controls,.25)
    # Exact nested angle is 3.5*u. In the final blend, weighted terms are
    # 1.25*u and 2.25*u; separately rounding them yields 3*u. The sum's
    # error is bounded by the two half-ulp product errors, not forced to
    # equal the single-rounded 4*u reference.
    exact = F(7,2)*F(U)
    assert abs(F(out.data[1])-exact) <= F(U)
    assert out.data[0] == 1. and out.data[2:] == [0.,0.]


@pytest.mark.parametrize('seed', range(12))
def test_general_spherical_paths_and_matrix_rotations_use_independent_tangent_reference(seed):
    rng = random.Random(475100+seed)
    first = [rng.uniform(-1,1) for _ in range(4)]
    length = math.hypot(*first); first = [x/length for x in first]
    tangent = [rng.uniform(-1,1) for _ in range(4)]
    projection = math.fsum(a*b for a,b in zip(first,tangent))
    tangent = [b-projection*a for a,b in zip(first,tangent)]
    length = math.hypot(*tangent); tangent = [x/length for x in tangent]
    angle = (.03,.4,.9,1.3)[seed%4]
    last = [math.cos(angle)*a+math.sin(angle)*b for a,b in zip(first,tangent)]
    if seed%2: last = [-x for x in last]
    controls = q.Quaternion(first[:]),q.Quaternion(last[:])
    for t in (0.,math.nextafter(0.,1.),.17,.5,.83,math.nextafter(1.,0.),1.):
        expected = [math.cos(t*angle)*a+math.sin(t*angle)*b for a,b in zip(first,tangent)]
        out = q.quat_slerp(*controls,t)
        assert out.data == pytest.approx(expected,rel=0.,abs=2e-15)
        point = [1.,-2.,3.]
        reference = exact_rotation(expected,point)
        assert q.quat_rotate_vector(out,Vector(3,point)).vector == pytest.approx(reference,rel=0.,abs=2e-14)
        assert (out.toMatrix()*Vector(4,point+[0.])).vector == pytest.approx(reference+[0.],rel=0.,abs=2e-14)
        assert abs(math.hypot(*out.data)-1.) <= 4*math.ulp(1.)


@pytest.mark.parametrize('sign', [-1.,1.])
def test_exact_half_turn_tie_still_retains_supplied_endpoint_sign(sign):
    first,last = q.Quaternion(),q.Quaternion([0.,sign,0.,0.])
    assert first.slerp(last,1.).data == last.data
    out = first.slerp(last,.5)
    assert out.data == pytest.approx([math.sqrt(.5),sign*math.sqrt(.5),0.,0.],rel=0.,abs=2e-16)


def test_legacy_nonunit_interpolation_and_parameter_policies_are_not_redefined():
    # Supported component LERP and legacy no-invert/SQUAD branches remain
    # unnormalized. General nonunit accurate SLERP is still unsupported.
    a,b,c = [q.Quaternion(x) for x in ([2.,0.,0.,0.],[2.,1.,0.,0.],[2.,0.,1.,0.])]
    assert q.quat_lerp(a,b,.5).data == [2.,.5,0.,0.]
    assert q.quat_slerp_no_invert(a,b,.5).data == [2.,.5,0.,0.]
    assert q.quat_squad(a,b,c,.5).data == [2.,.25,.25,0.]
    end = q.Quaternion([math.cos(math.pi/4),0.,0.,math.sin(math.pi/4)])
    assert q.quat_slerp(q.Quaternion(),end,1.5).data == pytest.approx(
        [math.cos(3*math.pi/8),0.,0.,math.sin(3*math.pi/8)],rel=0.,abs=3e-16)
    assert q.quat_slerp(q.Quaternion(),q.Quaternion(),float('nan')).data == [1.,0.,0.,0.]
