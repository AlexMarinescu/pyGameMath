"""Independent ray geometry and compatibility regressions."""
import math
import pytest
from gem import ray, vector, matrix, quaternion


def V(*values):
    return vector.Vector(len(values),list(values))


def geometric_tip(r):
    # Derived test geometry, deliberately separate from intersection state.
    return [p+d*r.distance for p,d in zip(r.start.vector,r.dir.vector)]


@pytest.mark.parametrize('distance', [0.25,5,17])
@pytest.mark.parametrize('end', [[0,0,0],[7,-8,9]])
def test_duplicate_exact_state_and_independent_storage(distance,end):
    r=ray.Ray(V(1,2,3),V(0,0,5))
    r.distance=distance
    r.dir.vector=[.2,.3,.4]  # Copy stored state, without re-normalization.
    r.end=V(*end)
    old=[v.vector for v in [r.start,r.dir,r.end]]
    snapshots=[v[:] for v in old]
    duplicate=r.duplicate()
    assert isinstance(duplicate,ray.Ray) and duplicate is not r
    assert duplicate.distance == distance
    for original,copy,storage,expected in zip([r.start,r.dir,r.end],
                                              [duplicate.start,duplicate.dir,duplicate.end],old,snapshots):
        assert isinstance(copy,vector.Vector) and copy.size==original.size
        assert copy is not original and copy.vector is not storage
        assert copy.vector == expected and original.vector is storage
        copy.vector[0]=99
        assert original.vector == expected
    r.end.vector[1]=123
    assert duplicate.end.vector[1] == end[1]


def test_constructor_ownership_and_placeholder_are_retained():
    start=V(1,2,3);direction=V(0,0,5)
    old_start=start.vector
    r=ray.Ray(start,direction)
    assert r.start is start and r.dir is direction
    assert start.vector is old_start and start.vector==[1,2,3]
    assert direction.vector==[0,0,1] and r.distance==5
    assert isinstance(r.end,vector.Vector) and r.end.vector==[0,0,0]
    assert r.end is not start and r.end is not direction


ROTATIONS = [
    ([1,0,0,0], [1,2,3], [2,-3,6], [3,-1,9]),
    ([math.sqrt(.5),0,0,math.sqrt(.5)], [-2,1,3], [3,2,6], [1,3,9]),
    ([0,0,1,0], [-1,2,-3], [-2,-3,-6], [-3,-1,-9]),
    ([.5,.5,.5,.5], [3,1,2], [6,2,-3], [9,3,-1]),
    ([math.sqrt(.5),0,0,-math.sqrt(.5)], [2,-1,3], [-3,-2,6], [-1,-3,9]),
]


@pytest.mark.parametrize('values,start,displacement,tip',ROTATIONS)
@pytest.mark.parametrize('end', [[0,0,0],[8,-7,6]])
def test_quaternion_rotations_known_geometry(values,start,displacement,tip,end):
    source_start=V(1,2,3);source_dir=V(2,-3,6)
    r=ray.Ray(source_start,source_dir);r.end=V(*end)
    q=quaternion.Quaternion(values[:]);q_storage=q.data
    state=r.end;state_storage=state.vector
    assert r.rotateUsingQuaternion(q) is None
    assert isinstance(r.start,vector.Vector) and isinstance(r.dir,vector.Vector)
    assert r.start.vector == pytest.approx(start,abs=1e-14)
    assert r.dir.vector == pytest.approx([v/7 for v in displacement],abs=1e-14)
    assert r.dir.magnitude()==pytest.approx(1,abs=1e-14)
    assert r.distance==7
    assert geometric_tip(r)==pytest.approx(tip,abs=1e-14)
    assert q.data is q_storage and q.data==values
    assert source_start.vector==[1,2,3]
    assert source_dir.vector==pytest.approx([2/7,-3/7,6/7])
    assert r.end is state and state.vector is state_storage and state.vector==end


@pytest.mark.parametrize('rows,start,displacement,tip',[
    ([[1,0,0],[0,1,0],[0,0,1]], [1,2,3],[2,-3,6],[3,-1,9]),
    ([[0,1,0],[-1,0,0],[0,0,1]], [-2,1,3],[3,2,6],[1,3,9]),
    ([[-1,0,0],[0,1,0],[0,0,-1]], [-1,2,-3],[-2,-3,-6],[-3,-1,-9]),
    ([[0,1,0],[0,0,1],[1,0,0]], [3,1,2],[6,2,-3],[9,3,-1]),
])
@pytest.mark.parametrize('end_values',[[0,0,0],[8,-7,6]])
def test_existing_matrix3_rotation(rows,start,displacement,tip,end_values):
    r=ray.Ray(V(1,2,3),V(2,-3,6));r.end=V(*end_values);end=r.end;end_storage=end.vector
    m=matrix.Matrix(3,[row[:] for row in rows]);storage=m.matrix
    assert r.roateUsingMatrix(m) is None
    assert r.start.vector==pytest.approx(start)
    assert r.dir.vector==pytest.approx([c/7 for c in displacement])
    assert r.distance==7 and r.end is end and end.vector is end_storage and end.vector==end_values
    assert geometric_tip(r)==pytest.approx(tip)
    assert m.matrix is storage and m.matrix==rows


@pytest.mark.parametrize('offset',[[0,0,0],[2,3,4],[-5,2,-7],[.25,-.5,.75]])
@pytest.mark.parametrize('end',[[0,0,0],[8,-7,6]])
def test_translation_known_geometry_and_storage(offset,end):
    source_start=V(1,2,3);source_dir=V(2,-3,6)
    r=ray.Ray(source_start,source_dir);r.end=V(*end)
    old_end=r.end;old_end_storage=old_end.vector
    before=r.dir.vector[:]
    rows=[[1,0,0,0],[0,1,0,0],[0,0,1,0],list(offset)+[1]]
    m=matrix.Matrix(4,[row[:] for row in rows]);storage=m.matrix
    export=[list(row) for row in m.c_matrix]
    assert r.translate(m) is None
    assert r.start.vector==pytest.approx([1+offset[0],2+offset[1],3+offset[2]])
    assert isinstance(r.start,vector.Vector) and r.start.size==3
    assert isinstance(r.dir,vector.Vector) and r.dir.size==3
    assert r.dir.vector==before and r.distance==7
    assert geometric_tip(r)==pytest.approx([3+offset[0],-1+offset[1],9+offset[2]])
    assert r.end is old_end and r.end.vector is old_end_storage and old_end.vector==end
    assert source_start.vector==[1,2,3]
    assert source_dir.vector==before
    assert m.matrix is storage and m.matrix==rows
    assert [list(row) for row in m.c_matrix]==export


def test_noncommuting_translation_and_rotation():
    rows=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[2,-3,4,1]]
    m=matrix.Matrix(4,rows)
    q=quaternion.Quaternion([math.sqrt(.5),0,0,math.sqrt(.5)])
    a=ray.Ray(V(1,2,3),V(2,-3,6));b=a.duplicate()
    a.translate(m);a.rotateUsingQuaternion(q)
    b.rotateUsingQuaternion(q);b.translate(m)
    assert a.start.vector==pytest.approx([1,3,7])
    assert b.start.vector==pytest.approx([0,-2,7])
    assert geometric_tip(a)==pytest.approx([4,5,13])
    assert geometric_tip(b)==pytest.approx([3,0,13])
    assert a.dir.vector==pytest.approx(b.dir.vector) and a.distance==b.distance==7


def test_no_general_matrix_vector_promotion_added():
    with pytest.raises(IndexError):
        matrix.Matrix(4)*V(1,2,3)
