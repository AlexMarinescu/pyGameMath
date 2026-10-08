import math
import random
import pytest
from gem import matrix, quaternion as q, vector
from .helpers import assert_matrix


def V(*xs):
    return vector.Vector(len(xs), list(xs))


@pytest.mark.parametrize('seed', range(20))
def test_quaternion_algebra(seed):
    rng = random.Random(seed)
    a,b,c = [q.Quaternion([rng.uniform(-2,2) for _ in range(4)]) for _ in range(3)]
    assert ((a*b)*c).data == pytest.approx((a*(b*c)).data)
    assert (a*q.Quaternion()).data == a.data
    assert (a*a.conjugate()).data == pytest.approx([a.magnitude()**2,0,0,0], abs=1e-14)
    assert a.normalize().magnitude() == pytest.approx(1)
    assert a.conjugate().conjugate().data == a.data


@pytest.mark.parametrize('values,expected', [
    ([1,2,3,4], [1/30,-2/30,-3/30,-4/30]),
    ([0,1,0,0], [0,-1,0,0]),
    ([2,0,0,0], [0.5,0,0,0]),
    ([0,0,2,0], [0,0,-0.5,0]),
    ([0,0,0,-4], [0,0,0,0.25]),
    ([1,0,0,0], [1,0,0,0]),
    ([-1,2,-3,4], [-1/30,-2/30,3/30,-4/30]),
])
def test_inverse_identity(values, expected):
    original = list(values)
    quat = q.Quaternion(values)
    result = quat.inverse()
    assert result is not quat
    assert q.quat_inverse(values) == pytest.approx(expected)
    assert result.data == pytest.approx(expected)
    assert (quat*result).data == pytest.approx([1,0,0,0], abs=1e-14)
    assert (result*quat).data == pytest.approx([1,0,0,0], abs=1e-14)
    assert values == original


@pytest.mark.parametrize('seed', range(20))
def test_inverse_general_quaternion(seed):
    rng = random.Random(seed)
    quat = q.Quaternion([rng.uniform(-2,2) for _ in range(4)])
    result = quat.inverse()
    assert (quat*result).data == pytest.approx([1,0,0,0], abs=1e-14)
    assert (result*quat).data == pytest.approx([1,0,0,0], abs=1e-14)


@pytest.mark.parametrize('axis', [[1,0,0],[0,1,0],[0,0,1],[1,2,3]])
@pytest.mark.parametrize('angle', [0,45,90,135])
def test_matrix_rotation_agreement(axis,angle):
    rot = q.quat_from_axis_angle(V(*axis),angle)
    m = matrix.Matrix(4).rotate(V(*axis),angle)
    assert_matrix(rot.toMatrix().matrix, m.matrix, abs=1e-14)
    out = q.quat_rotate_vector(rot,V(2,3,4))
    assert out.vector == pytest.approx((m*V(2,3,4,0)).vector[:3])
    assert out.magnitude() == pytest.approx(V(2,3,4).magnitude())


@pytest.mark.parametrize('axis', [[1,0,0],[0,1,0],[0,0,1],[1,2,3]])
@pytest.mark.defect('Q02')
def test_half_turn_matrix_roundtrip(axis):
    rot = q.quat_from_axis_angle(V(*axis),180)
    recovered = q.quat_from_matrix(rot.toMatrix())
    assert abs(recovered.dot(rot)) == pytest.approx(1)


def test_small_rotation_matrix_roundtrip():
    rot = q.quat_from_axis_angle(V(1,2,3),45)
    assert abs(q.quat_from_matrix(rot.toMatrix()).dot(rot)) == pytest.approx(1)


@pytest.mark.defect('Q03')
def test_power_one():
    rot = q.quat_from_axis_angle(V(0,0,1),90)
    assert rot.pow(1).data == pytest.approx(rot.data)


@pytest.mark.defect('Q03')
def test_identity_power():
    assert q.Quaternion().pow(0.5).data == [1,0,0,0]


@pytest.mark.defect('Q04')
def test_unit_log():
    rot = q.quat_from_axis_angle(V(0,0,1),90)
    assert rot.log() == pytest.approx([0,0,0,math.pi/4])


@pytest.mark.defect('Q05')
def test_squad_identical():
    rot = q.Quaternion()
    assert rot.squad(rot,rot,0.5).data == [1,0,0,0]


@pytest.mark.parametrize('function',[q.quat_from_axis_angle,q.quat_rotate_from_axis_angle])
@pytest.mark.defect('Q06')
def test_list_axis(function):
    assert isinstance(function([0,0,1],90),q.Quaternion)


@pytest.mark.defect('Q07')
def test_rotation_matrix_ctypes_snapshot():
    m = q.quat_from_axis_angle(V(0,0,1),90).toMatrix()
    assert_matrix([list(row) for row in m.c_matrix],m.matrix,abs=1e-6)


@pytest.mark.defect('Q08')
def test_inplace_vector_multiply():
    rot = q.Quaternion()
    rot *= V(1,2,3)
    assert rot.data == [0,1,2,3]


@pytest.mark.defect('Q09')
def test_python3_quaternion_division():
    assert (q.Quaternion([2,4,6,8])/2.0).data == [1,2,3,4]


@pytest.mark.parametrize('t', [0.0,0.25,0.5,0.75,1.0])
def test_slerp_known_answer(t):
    a,b = q.Quaternion(),q.quat_from_axis_angle(V(0,0,1),90)
    assert a.slerp(b,t).data == pytest.approx([math.cos(t*math.pi/4),0,0,math.sin(t*math.pi/4)])
    assert a.slerp(b.negate(),t).data == pytest.approx(a.slerp(b,t).data)


@pytest.mark.contract_question('Q10-unit-accuracy')
def test_slerp_nearby_unit_length():
    a,b = q.Quaternion(),q.quat_from_axis_angle(V(0,0,1),1)
    assert a.slerp(b,0.5).magnitude() == pytest.approx(1,abs=1e-12)


@pytest.mark.parametrize('function', [q.quat_rotate_x_from_angle,q.quat_rotate_y_from_angle,q.quat_rotate_z_from_angle])
def test_axis_helpers_use_radians(function):
    assert function(math.pi)[0] == pytest.approx(0,abs=1e-15)


@pytest.mark.defect('N01')
def test_zero_quaternion_normalization():
    assert q.Quaternion([0,0,0,0]).normalize().data == [1,0,0,0]


def test_quaternion_vector_product_is_not_rotation():
    assert (q.Quaternion()*V(1,2,3)).data == [0,1,2,3]
