import math
import random
import pytest
from gem import matrix, vector
from .helpers import determinant, inverse, assert_matrix


def V(*xs):
    return vector.Vector(len(xs), list(xs))


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('seed', range(20))
def test_determinant_and_inverse(size, seed):
    rng = random.Random(seed)
    a = [[rng.randint(-3,3) + (10 if i == j else 0) for j in range(size)] for i in range(size)]
    m = matrix.Matrix(size, a)
    assert m.det() == determinant(a)
    assert_matrix(m.inverse().matrix, inverse(a))
    assert_matrix((m*m.inverse()).matrix, matrix.identity(size), abs=1e-12)
    assert_matrix((m.inverse()*m).matrix, matrix.identity(size), abs=1e-12)


@pytest.mark.parametrize('values,expected', [
    ([[1,2],[3,4]], [[-2,1],[1.5,-0.5]]),
    ([[4,7],[2,6]], [[0.6,-0.7],[-0.2,0.4]]),
    ([[2,0],[0,-4]], [[0.5,0],[0,-0.25]]),
    ([[0,2],[4,0]], [[0,0.25],[0.5,0]]),
])
def test_inverse2_known_answer(values, expected):
    original = [row[:] for row in values]
    assert_matrix(matrix.inverse2(values), expected)
    m = matrix.Matrix(2, values)
    result = m.inverse()
    assert result is not m
    assert_matrix(result.matrix, expected)
    assert_matrix((m*result).matrix, matrix.identity(2), abs=1e-12)
    assert_matrix((result*m).matrix, matrix.identity(2), abs=1e-12)
    assert values == original
    assert m.i_inverse() is m
    assert_matrix(m.matrix, expected)
    assert_matrix([list(row) for row in m.c_matrix], expected, abs=1e-7)


@pytest.mark.parametrize('size', [2,3,4])
def test_singular_inverse_current_contract(size):
    with pytest.raises(ZeroDivisionError):
        matrix.Matrix(size, matrix.zero_matrix(size)).inverse()


@pytest.mark.parametrize('size', [2,3,4,5])
def test_matrix_algebra(size):
    a = matrix.Matrix(size, [[i*size+j for j in range(size)] for i in range(size)])
    assert_matrix((a*matrix.Matrix(size)).matrix, a.matrix)
    assert_matrix(a.transpose().transpose().matrix, a.matrix)
    b = matrix.Matrix(size)
    b.i_transpose()
    assert_matrix(b.matrix, matrix.identity(size))


def test_storage_and_row_vector_convention():
    m = matrix.Matrix(2, [[1,2],[3,4]])
    assert (m*V(5,6)).vector == [23,34]
    assert [list(row) for row in m.c_matrix] == [[1,2],[3,4]]
    t = matrix.Matrix(4).translate(V(2,3,4))
    assert (t*V(1,2,3,1)).vector == [3,5,7,1]


@pytest.mark.parametrize('angle', [-180, -90, 0, 45, 90, 180])
def test_rotations(angle):
    r = matrix.Matrix(3).rotate(V(0,0,1), angle)
    assert (r*V(1,0,0)).vector == pytest.approx([math.cos(math.radians(angle)),math.sin(math.radians(angle)),0], abs=1e-15)
    assert_matrix((r*r.transpose()).matrix, matrix.identity(3), abs=1e-14)
    assert r.det() == pytest.approx(1)


def test_scalar_division_helper():
    assert_matrix(matrix.matrix_div([[2,4],[6,8]], 2), [[1,2],[3,4]])


def test_python3_division():
    assert_matrix((matrix.Matrix(2, [[2,4],[6,8]])/2.0).matrix, [[1,2],[3,4]])


def test_inplace_translate_matches_returning():
    a, b = matrix.Matrix(4), matrix.Matrix(4)
    a.i_translate(V(2,3,4))
    assert_matrix(a.matrix, b.translate(V(2,3,4)).matrix)


@pytest.mark.defect('M05')
def test_shear_xy3():
    assert len(matrix.shearXY3(1,2)) == 3


@pytest.mark.defect('M06')
def test_rotate2_fixed_pivot():
    r = matrix.Matrix(3, matrix.rotate2([2,3], 90))
    assert (r*V(2,3,1)).vector == pytest.approx([2,3,1])


def test_matrix2_rotation():
    r = matrix.Matrix(2).rotate(V(0,0), 90)
    assert (r*V(1,0)).vector == pytest.approx([0,1], abs=1e-14)


@pytest.mark.parametrize('f', [matrix.perspective, matrix.perspectiveX])
def test_projection_depth(f):
    p = f(90, 1, 1, 10)
    for z, expected in [(-1,-1),(-10,1)]:
        clip = p*V(0,0,z,1)
        assert clip.vector[2]/clip.vector[3] == pytest.approx(expected)


def test_orthographic_and_lookat():
    p = matrix.orthographic(-2,2,-3,3,1,10)
    assert (p*V(-2,-3,-1,1)).vector == pytest.approx([-1,-1,-1,1])
    eye = V(1,2,3)
    view = matrix.lookAt(eye,V(1,2,2),V(0,1,0))
    assert (view*V(1,2,3,1)).vector == [0,0,0,1]
    assert_matrix(matrix.lookAt(V(0,0,0),V(0,0,-1),V(0,1,0)).matrix, matrix.identity(4))


@pytest.mark.parametrize('use_objects', [False,True])
def test_project_center(use_objects):
    m = matrix.Matrix(4)
    arg = m if use_objects else m.matrix
    assert matrix.project(V(0,0,0,1),arg,arg,[0,0,100,100]).vector[:3] == [50,50,0.5]


def test_unproject_identity():
    assert matrix.unproject(50,50,0.5,matrix.Matrix(4),matrix.Matrix(4),[0,0,100,100]).vector == [0,0,0]


def test_unproject_noncommuting_transforms():
    model = matrix.Matrix(4).translate(V(2,0,0))
    projection = matrix.orthographic(-4,4,-4,4,1,10)
    point = V(1,1,-3,1)
    clip = projection*(model*point)
    ndc = [x/clip.vector[3] for x in clip.vector[:3]]
    window = [(ndc[0]+1)*50,(ndc[1]+1)*50,(ndc[2]+1)/2]
    assert matrix.unproject(*window,model,projection,[0,0,100,100]).vector == pytest.approx(point.vector[:3])


def test_scale_and_composition_order():
    t = matrix.Matrix(4).translate(V(2,3,4))
    s = matrix.Matrix(4).scale(V(2,2,2))
    v = V(1,1,1,1)
    assert ((t*s)*v).vector == (s*(t*v)).vector == [6,8,10,1]


@pytest.mark.parametrize('scale',[1e-100,1e100])
@pytest.mark.defect('N03')
def test_inverse4_extreme_uniform_scale(scale):
    a = matrix.Matrix(4,[[scale if i==j else 0 for j in range(4)] for i in range(4)])
    assert_matrix(a.inverse().matrix,[[1/scale if i==j else 0 for j in range(4)] for i in range(4)],abs=0)


def test_matrix5_scale_characterization():
    # NxN multiplication exists, but scaling treats all indices >=3 as homogeneous.
    scaled = matrix.Matrix(5).scale(V(2,3,4,5,6))
    assert [scaled.matrix[i][i] for i in range(5)] == [2,3,4,1,1]
