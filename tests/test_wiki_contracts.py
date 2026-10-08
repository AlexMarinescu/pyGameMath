"""Unambiguous historical wiki requirements, with documented example corrections.

Frozen evidence: audit/wiki-snapshot at wiki commit 715e5039c75e080814a12e957f5148c35cdf8bda.
No malformed example is silently treated as the desired library contract.
"""
import pytest
from gem import matrix, vector
from .helpers import assert_matrix


def V(*values):
    return vector.Vector(len(values),list(values))


@pytest.mark.parametrize('size',[1,2,3,4,8])
def test_wiki_vector_nd_construction(size):
    """Vector page lines 65–85: default zero and supplied values."""
    assert vector.Vector(size).vector == [0]*size
    values = list(range(1,size+1))
    assert vector.Vector(size,data=values).vector == values
    assert vector.Vector(size,values).vector == values
    assert vector.Vector(size).one().vector == [1]*size
    assert vector.Vector(size).zero().vector == [0]*size


@pytest.mark.parametrize('op,expected',[
    ('add',[6,7,13]),('sub',[4,-1,3]),('mul',[10,6,16]),
    ('div',[1/3,4/3,5/3]),('neg',[-5,-3,-8]),('cross',[-17,-17,17]),
])
def test_wiki_vector_operator_examples(op,expected):
    """Vector page lines 92–119: unrounded mathematical answers."""
    a,b = V(5,3,8),V(1,4,5)
    operations = {'add':lambda:a+b,'sub':lambda:a-b,'mul':lambda:a*2.0,
                  'div':lambda:b/3.0,'neg':lambda:-a,'cross':lambda:vector.cross(a,b)}
    result = operations[op]()
    assert isinstance(result,vector.Vector)
    assert result.size == 3
    assert result.vector == pytest.approx(expected)
    assert a.vector == [5,3,8] and b.vector == [1,4,5]


def test_wiki_vector_dot_example():
    assert V(5,3,8).dot(V(1,4,5)) == 57


@pytest.mark.parametrize('name,expected',[
    ('right',[1,0,0]),('left',[-1,0,0]),('front',[0,0,-1]),
    ('back',[0,0,1]),('up',[0,1,0]),('down',[0,-1,0])])
def test_wiki_vector_axis_identities(name,expected):
    assert getattr(vector.Vector(3),name)().vector == expected


@pytest.mark.parametrize('name,indices',[
    ('xy',[0,1]),('yz',[1,2]),('xz',[0,2]),('xw',[0,3]),
    ('yw',[1,3]),('zw',[2,3]),('xyw',[0,1,3]),('yzw',[1,2,3]),('xzw',[0,2,3])])
def test_wiki_swizzle_contracts(name,indices):
    v = V(11,22,33,44)
    result = getattr(v,name)()
    assert result.size == len(indices)
    assert result.vector == [v.vector[i] for i in indices]
    result.vector[0] = 99
    assert v.vector == [11,22,33,44]


def test_wiki_clone_and_normalization_example():
    # Correct only the documented `data[...]` typo to `data=[...]`.
    a = vector.Vector(3,data=[2.0,4.0,8.0])
    b = a.clone()
    assert b is not a and b.vector is not a.vector
    n = a.normalize()
    assert n is not a and a.vector == [2,4,8]
    b.i_normalize()
    assert n.vector == pytest.approx([0.218218,0.436436,0.872872],abs=1e-6)
    assert b.vector == pytest.approx(n.vector)
    assert a.vector == [2,4,8]


@pytest.mark.parametrize('name,expected',[
    ('maxV',[5,4,8]),('minV',[1,3,5])])
def test_wiki_vector_component_extrema(name,expected):
    assert getattr(V(5,3,8),name)(V(1,4,5)).vector == expected


def test_wiki_scalar_extrema_are_unary():
    # maxS/minS descriptions mistakenly say "between two vectors";
    # documented signatures take no second vector.
    assert V(5,3,8).maxS() == 8
    assert V(5,3,8).minS() == 3
    assert V(3,4,0).magnitude() == 5


@pytest.mark.parametrize('size',[1,2,3,4,8])
@pytest.mark.parametrize('t',[0.0,0.5,1.0])
def test_wiki_lerp_returns_input_dimension(size,t):
    a,b = vector.Vector(size,[1]*size),vector.Vector(size,[3]*size)
    out = vector.lerp(a,b,t)
    assert out.size == size and out.vector == [1+2*t]*size
    assert a.vector == [1]*size and b.vector == [3]*size


def test_wiki_clamp_returning_vs_inplace_receiver():
    """General non-i/i rule: do not infer ownership of separate `value` list."""
    receiver = V(9,9,9)
    out = receiver.clamp(3,[-2,2,8],[0]*3,[5]*3)
    assert isinstance(out,vector.Vector) and out is not receiver
    assert out.vector == [0,2,5] and receiver.vector == [9,9,9]
    receiver.i_clamp(3,[-2,2,8],[0]*3,[5]*3)
    assert receiver.vector == [0,2,5]


@pytest.mark.parametrize('size',[2,3,4,8])
def test_wiki_matrix_explicit_constructors_and_identities(size):
    data=matrix.identity(size)
    assert_matrix(matrix.Matrix(size,data=data).matrix,data)
    assert_matrix(matrix.Matrix(size,data=matrix.zero_matrix(size)).matrix,[[0]*size for _ in range(size)])


def test_matrix_default_identity_historical_characterization():
    # Wiki prints zero, but both contemporaneous 2015 and current source use identity.
    # Preserve the established behavior rather than making a contradictory example normative.
    assert_matrix(matrix.Matrix(4).matrix,matrix.identity(4))


@pytest.mark.parametrize('size',[1,2,3,4,8])
def test_wiki_dimension_matched_multiplication(size):
    a,b = matrix.Matrix(size,data=matrix.identity(size)),matrix.Matrix(size,data=matrix.identity(size))
    v=vector.Vector(size,list(range(size)))
    product=a*b
    assert isinstance(product,matrix.Matrix) and product.size==size
    out=a*v
    assert isinstance(out,vector.Vector) and out.size==size
    assert out.vector==v.vector


def test_wiki_matrix_scaling_example_corrected():
    # Matrix page line 103: fix `indetity` and choose size 4 to match its 4x4 data/output.
    a = matrix.Matrix(4,data=matrix.identity(4))
    s = V(2,3,4)
    scaled = a.scale(s)
    expected = [[2,0,0,0],[0,3,0,0],[0,0,4,0],[0,0,0,1]]
    assert scaled is not a
    assert_matrix(scaled.matrix,expected)
    assert_matrix(a.matrix,matrix.identity(4))
    a.i_scale(s)
    # Wiki prints identity afterward, contradicting its in-place prose and determinant 24.
    assert_matrix(a.matrix,expected)
    assert a.det()==24


def test_wiki_python3_matrix_division_example():
    a=matrix.Matrix(4,data=matrix.identity(4))  # fix `indetity` typo only
    result=a/2.0
    assert isinstance(result,matrix.Matrix)
    assert_matrix(result.matrix,[[0.5 if i==j else 0 for j in range(4)] for i in range(4)])


@pytest.mark.parametrize('size',[2,3,4])
def test_wiki_inverse_returning_contract(size):
    a=matrix.Matrix(size,[[i+2 if i==j else 0 for j in range(size)] for i in range(size)])
    before=[row[:] for row in a.matrix]
    if size==2:
        # Value defect M01 has separate regressions; here check return type/receiver ownership.
        out=a.inverse()
    else:
        out=a.inverse()
        assert_matrix((a*out).matrix,matrix.identity(size))
    assert isinstance(out,matrix.Matrix) and out is not a
    assert_matrix(a.matrix,before)


@pytest.mark.parametrize('name,args',[
    ('scale',(V(2,3,4),)),('rotate',(V(0,0,1),30)),('translate',(V(2,3,4,0),)),
    ('transpose',()),('inverse',()),('shearXY',(0.2,0.3)),('shearYZ',(0.2,0.3)),('shearXZ',(0.2,0.3))])
def test_wiki_matrix_returning_and_inplace_ownership(name,args):
    # Use supported inputs. In-place Vector3 translation's M04 defect is tested elsewhere.
    source=matrix.Matrix(4)
    before=[row[:] for row in source.matrix]
    out=getattr(source,name)(*args)
    assert isinstance(out,matrix.Matrix) and out is not source
    assert_matrix(source.matrix,before)
    target=matrix.Matrix(4)
    getattr(target,'i_'+name)(*args)
    assert_matrix(target.matrix,out.matrix)


@pytest.mark.parametrize('name,args',[
    ('orthographic',(-2,2,-3,3,1,10)),('perspective',(90,2,1,10)),
    ('perspectiveX',(90,2,1,10)),('lookAt',(V(0,0,0),V(0,0,-1),V(0,1,0)))])
def test_wiki_projection_function_return_types(name,args):
    out=getattr(matrix,name)(*args)
    assert isinstance(out,matrix.Matrix) and out.size==4
    assert len(out.matrix)==4 and all(len(row)==4 for row in out.matrix)


def test_wiki_horizontal_vs_vertical_fov():
    # Page names vertical FOVY and horizontal FOVX, not the angle unit.
    py,px = matrix.perspective(90,2,1,10),matrix.perspectiveX(90,2,1,10)
    assert py.matrix[0][0] == pytest.approx(0.5)
    assert py.matrix[1][1] == pytest.approx(1)
    assert px.matrix[0][0] == pytest.approx(1)
    assert px.matrix[1][1] == pytest.approx(2)


def test_wiki_project_returns_three_components():
    # Lists match the source's initial flattening code. No window-depth policy is asserted.
    out=matrix.project(V(0,0,0,1),matrix.identity(4),matrix.identity(4),[0,0,100,100])
    assert isinstance(out,vector.Vector)
    assert out.size==len(out.vector)==3


def test_wiki_unproject_returns_object_vector():
    out=matrix.unproject(50,50,0.5,matrix.Matrix(4),matrix.Matrix(4),[0,0,100,100])
    assert isinstance(out,vector.Vector) and out.size==3
    assert out.vector==[0,0,0]
