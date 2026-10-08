"""Additional public API characterization; policy choices are not called bugs."""
import math
import pytest
from gem import common, matrix, plane, quaternion as q, vector, ray
from .helpers import assert_matrix


def V(*xs):
    return vector.Vector(len(xs),list(xs))


@pytest.mark.parametrize('size',[2,3,4])
def test_inplace_matrix_methods(size):
    axis = V(0,0,1) if size > 2 else V(0,0)
    a = matrix.Matrix(size)
    assert a.i_rotate(axis,30) is a
    b = matrix.Matrix(size).rotate(axis,30)
    assert_matrix(a.matrix,b.matrix)
    scale = V(2,3,4)
    assert a.i_scale(scale) is a
    assert_matrix(a.matrix,b.scale(scale).matrix)
    if size > 2:
        a.i_inverse()
        assert_matrix(a.matrix,b.scale(scale).inverse().matrix)
    assert_matrix([list(row) for row in a.c_matrix],a.matrix,rel=1e-6,abs=1e-6)


@pytest.mark.parametrize('name,args',[('shearYZ',(0.2,0.3)),('shearXZ',(0.2,0.3))])
@pytest.mark.parametrize('size',[3,4])
def test_inplace_shear(name,args,size):
    a,b = matrix.Matrix(size),matrix.Matrix(size)
    getattr(a,'i_'+name)(*args)
    assert_matrix(a.matrix,getattr(b,name)(*args).matrix)


def test_inplace_vector_methods():
    v = V(1,2,3)
    v += V(2,3,4)
    v -= 1
    v *= 2
    v /= 2
    assert v.vector == [2,4,6]
    assert v.i_normalize() is v
    assert v.magnitude() == pytest.approx(1)
    assert V(1,2,3).maxV(V(2,1,4)).vector == [2,2,4]
    assert V(1,2,3).minV(V(2,1,4)).vector == [1,1,3]
    assert V(1,2,3).minS() == 1
    assert V(1,2,3).maxS() == 3
    assert V(1,2,3,4).xy().vector == [1,2]
    assert V(1,2,3,4).xz().vector == [1,3]
    assert V(1,2,3,4).yz().vector == [2,3]
    assert V(1,2,3,4).xw().vector == [1,4]
    assert V(1,2,3,4).yw().vector == [2,4]
    assert V(1,2,3,4).zw().vector == [3,4]
    assert V(1,2,3,4).xyw().vector == [1,2,4]
    assert V(1,2,3,4).xzw().vector == [1,3,4]
    assert V(1,2,3,4).xyz().vector == [1,2,3]
    assert (-V(1,2,3)).vector == V(1,2,3).negate().vector


def test_inplace_quaternion_methods():
    a = q.Quaternion([1,2,3,4])
    a += q.Quaternion()
    a -= q.Quaternion()
    a *= 2.0
    assert a.data == [2,4,6,8]
    assert a.i_normalize() is a
    assert a.magnitude() == pytest.approx(1)
    assert a.i_conjugate() is a
    assert a.i_negate() is a
    assert a.i_identity() is a
    assert a.data == a.identity().data == [1,0,0,0]
    a *= q.Quaternion([0,1,0,0])
    assert a.data == [0,1,0,0]


def test_current_boundary_errors():
    with pytest.raises(ZeroDivisionError):
        V(1,2)/0
    with pytest.raises(ZeroDivisionError):
        V(0,0,0).barycentric(V(0,0,0),V(1,1,1),V(2,2,2))
    with pytest.raises(ZeroDivisionError):
        q.Quaternion([0,0,0,0]).inverse()
    with pytest.raises(NotImplementedError):
        matrix.Matrix(5).inverse()
    with pytest.raises(NotImplementedError):
        matrix.Matrix(5).det()
    with pytest.raises(ValueError):
        matrix.Matrix(3)*matrix.Matrix(2)
    with pytest.raises(TypeError):
        matrix.Matrix(3).rotate([0,0,1],90)
    assert (V(1,2) == object()) is False


@pytest.mark.defect('R01')
def test_duplicate_aliases():
    r = ray.Ray(V(1,2,3),V(0,0,5))
    copy = r.duplicate()
    copy.start.vector[0] = 99
    assert r.start.vector[0] == 1


def test_ray_constructor_mutates_direction_characterization():
    direction = V(0,0,5)
    ray.Ray(V(0,0,0),direction)
    assert direction.vector == [0,0,1]  # Existing ownership convention, not changed.


def test_forward_convention_mismatch_characterization():
    assert vector.Vector(3).front().vector == [0,0,-1]
    assert q.Quaternion().getForward().vector == [0,0,1]
    assert q.Quaternion().getBack().vector == [0,0,-1]
    assert q.Quaternion().getLeft().vector == [-1,0,0]
    assert q.Quaternion().getRight().vector == [1,0,0]
    assert q.Quaternion().getUp().vector == [0,1,0]
    assert q.Quaternion().getDown().vector == [0,-1,0]


def test_angle_units_characterization():
    # Preserve this mixed-unit API until a compatibility policy is approved.
    assert_matrix(matrix.rotate_origin2(math.pi/2), matrix.rotate3([0,0,1],90),abs=1e-14)
    assert q.quat_rotate(V(1,0,0),[0,0,1],90).vector == pytest.approx([0,1,0],abs=1e-14)


@pytest.mark.contract_question('Q11-return-semantics')
def test_arbitrary_axis_helper_returns_rotation_quaternion():
    rotation = q.quat_rotate_from_axis_angle(V(0,0,1),90)
    expected = q.quat_from_axis_angle(V(0,0,1),90)
    assert rotation.data == pytest.approx(expected.data)


@pytest.mark.defect('M07')
def test_legacy_inplace_division_ctypes_sync():
    a = matrix.Matrix(2,[[2,0],[0,4]])
    a.__idiv__(2.0)
    assert_matrix([list(row) for row in a.c_matrix],a.matrix)


def test_sinc_small_argument_error_bound():
    x = 1e-5
    assert abs(common.sinc(x)-math.sin(x)/x) < 2e-11


@pytest.mark.parametrize('size',[3,4])
def test_ill_conditioned_inverse_characterization(size):
    # No promise of backward stability: quantify an ordinary small-scale case.
    a = matrix.Matrix(size, [[(1e-12 if i == 0 else i+1) if i == j else 0 for j in range(size)] for i in range(size)])
    assert_matrix((a*a.inverse()).matrix,matrix.identity(size))
