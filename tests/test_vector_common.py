import ctypes
import math
import random
import pytest
from gem import common, vector


def V(*values):
    return vector.Vector(len(values), list(values))


@pytest.mark.parametrize('size', [0, 1, 2, 3, 4, 8])
def test_fresh_storage(size):
    a, b = vector.Vector(size), vector.Vector(size)
    assert a.vector == [0] * size
    a.one()
    assert b.vector == [0] * size
    if size:
        c = a.clone()
        c.vector[0] = 99
        assert a.vector[0] == 1


@pytest.mark.parametrize('seed', range(20))
def test_vector_properties(seed):
    rng = random.Random(seed)
    a, b, c = [V(*(rng.uniform(-10, 10) for _ in range(3))) for _ in range(3)]
    assert ((a+b)-b).vector == pytest.approx(a.vector)
    assert a.dot(b+c) == pytest.approx(a.dot(b)+a.dot(c))
    cross = vector.cross(a, b)
    assert cross.dot(a) == pytest.approx(0, abs=1e-11)
    assert cross.dot(b) == pytest.approx(0, abs=1e-11)
    assert a.normalize().magnitude() == pytest.approx(1)
    assert (a / 2).vector == pytest.approx((a*0.5).vector)
    assert vector.lerp(a, b, 0).vector == pytest.approx(a.vector)
    assert vector.lerp(a, b, 1).vector == pytest.approx(b.vector)


@pytest.mark.parametrize('op', ['eq', 'ne'])
def test_equality_all_components(op):
    a, b = V(1, 2, 3), V(1, 9, 3)
    assert (a == b) is False if op == 'eq' else (a != b) is True


@pytest.mark.parametrize('size', [1, 2, 3, 4, 8])
def test_equality_each_component_exact(size):
    values = [1.0] * size
    a = V(*values)
    assert (a == V(*values)) is True
    assert (a != V(*values)) is False
    for index in range(size):
        changed = list(values)
        changed[index] += 1e-12
        b = V(*changed)
        assert (a == b) is False
        assert (b == a) is False
        assert (a != b) is True
        assert (b != a) is True


def test_equality_unsupported_operand():
    a = V(1, 2, 3)
    assert a.__eq__([1, 2, 3]) is NotImplemented
    assert a.__ne__([1, 2, 3]) is NotImplemented


@pytest.mark.contract_question('V01-dimension-policy')
def test_equality_dimensions():
    assert (V(1, 2) == V(1, 2, 3)) is False


@pytest.mark.contract_question('V01-empty-policy')
def test_empty_equality():
    assert (V() == V()) is True


def test_zero_normalization():
    assert V(0, 0, 0).normalize().vector == [0, 0, 0]


@pytest.mark.parametrize('scale', [1e200, 1e-200])
def test_extreme_magnitude(scale):
    assert V(scale, scale).magnitude() == pytest.approx(math.hypot(scale, scale), rel=1e-14, abs=0)


def test_known_geometry():
    assert vector.cross(V(1,0,0), V(0,1,0)).vector == [0,0,1]
    assert vector.reflect(V(1,-1,0), V(0,1,0)).vector == [1,1,0]
    assert V(0.25,0.25,0).barycentric(V(0,0,0), V(1,0,0), V(0,1,0)) == [0.5,0.25,0.25]
    assert V(1,2,3,4).yzw().vector == [2,3,4]
    assert V(1,2,3).isInSameDirection(V(1,0,0))
    assert common.list_2d_to_1d([[1,2],[3,4]]) == [1,2,3,4]
    assert common.convertArr([1,2,3,4], 2) == [[1,2],[3,4]]
    assert list(common.conv_list([1,2], ctypes.c_double)) == [1,2]
    assert common.mulV4([1,2,3,4],[2,2,2,2]) == [2,4,6,8]
    assert common.convertM4to3([[1,2,3,4]]*4) == [[1,2,3]]*3
    assert common.sinc(0) == 1
    assert common.sinc(0.5) == pytest.approx(math.sin(0.5)/0.5)
    assert common.scalarLerp(2,6,0.25) == 3
    assert vector.toAngle([0,1]) == pytest.approx(math.pi/2)
    assert vector.lperp([1,0]).vector == [0,1]
    assert vector.rperp([1,0]).vector == [0,-1]


@pytest.mark.parametrize('function,value,expected', [
    (common.radiansToDegrees, math.pi, 180),
    (common.degreesToRadians, 180, math.pi)])
def test_angle_conversion(function, value, expected):
    assert function(value) == pytest.approx(expected)


def test_refraction_normal_incidence():
    assert vector.refract(0.5, V(0,-1,0), V(0,1,0)).vector == pytest.approx([0,-1,0])


def test_refraction_critical_angle():
    # eta=1.5, sin(theta)=0.6: k=0.19 > 0, not total internal reflection.
    out = vector.refract(1.5, V(0.6,-0.8,0), V(0,1,0))
    assert out.vector == pytest.approx([0.9,-math.sqrt(0.19),0])


def test_total_internal_reflection():
    assert vector.refract(2.0, V(0.8,-0.6,0), V(0,1,0)).vector == [0,0,0]


def test_transform_identity():
    assert vector.transform(3, [2,3,4], [[1,0,0],[0,1,0],[0,0,1]]) == [2,3,4]


@pytest.mark.contract_question('V04-value-list-ownership')
def test_clamp_preserves_input():
    values = [-2, 2, 10]
    assert vector.clamp(3, values, [0]*3, [5]*3).vector == [0,2,5]
    assert values == [-2,2,10]


@pytest.mark.defect('C02')
def test_viewport_vector():
    assert len(common.getViewPort(V(1,1), 100, 100)) == 4


def test_quaternion_extreme_magnitude():
    from gem.quaternion import Quaternion
    assert Quaternion([1e200]*4).magnitude() == pytest.approx(2e200)
