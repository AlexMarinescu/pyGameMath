"""Decimal/Fraction known answers for finite binary64 numerical operations."""
import math
import sys
from decimal import Decimal, localcontext
from fractions import Fraction
import pytest
from gem import vector, quaternion, matrix, plane, ray, common
from .helpers import inverse as fraction_inverse, assert_matrix


def V(values):return vector.Vector(len(values),values[:])


VALUES = [[0,0,0],[3,4,0],[1e200,-2e200,2e200],[1e-200,-2e-200,2e-200],
          [1e300,1e-300,-1e300],[1e-300,1e-308,-2e-300],
          [math.ldexp(1.,-1074)]*3,[1e308,1e308,1e308],
          [sys.float_info.max]*3,[1,math.ldexp(1.,-52),-1]]


@pytest.mark.parametrize('values',VALUES)
@pytest.mark.parametrize('kind',['vector','quaternion'])
def test_decimal_length_and_normalization(values,kind):
    data=values[:] if kind=='vector' else values[:]+[0]
    with localcontext() as ctx:
        ctx.prec=800
        exact=[Decimal.from_float(float(c)) for c in data]
        length=sum(c*c for c in exact).sqrt()
        expected_length=float(length)
        expected=[float(c/length) for c in exact] if length else ([0]*len(data) if kind=='vector' else [1,0,0,0])
    obj=V(data) if kind=='vector' else quaternion.Quaternion(data[:])
    storage=obj.vector if kind=='vector' else obj.data
    result=obj.normalize()
    output=result.vector if kind=='vector' else result.data
    assert isinstance(result,type(obj)) and result is not obj
    assert output is not storage
    if math.isinf(expected_length):assert math.isinf(obj.magnitude())
    else:assert obj.magnitude()==pytest.approx(expected_length,rel=2e-15,abs=0)
    assert output==pytest.approx(expected,rel=2e-15,abs=math.ldexp(1.,-1074))
    assert result.magnitude()==pytest.approx(0 if not length and kind=='vector' else 1,rel=2e-15,abs=0)
    assert storage==data
    receiver=obj
    assert obj.i_normalize() is receiver
    assert (obj.vector if kind=='vector' else obj.data)==output
    assert storage==data


@pytest.mark.parametrize('size',[0,1,2,3,4,7])
def test_zero_vector_dimensions(size):
    data=[-0.0]*size
    assert vector.normalize(size,data)==[0]*size
    assert vector.magnitude(size,data)==0


BASES={3:[[2,1,0],[0,2,1],[0,0,2]],4:[[2,1,0,1],[0,3,1,0],[0,0,2,1],[1,0,0,3]]}


@pytest.mark.parametrize('size',[3,4])
@pytest.mark.parametrize('scale',[1e-300,1e-200,1e-100,1,1e100,1e200,1e300])
def test_fraction_scaled_inverse_known_answers(size,scale):
    base=BASES[size]
    rows=[[scale*c for c in row] for row in base]
    expected=[[float(c) for c in row] for row in fraction_inverse(rows)]
    original=[row[:] for row in rows]
    wrapper=matrix.Matrix(size,rows)
    exported=[list(row) for row in wrapper.c_matrix]
    result=wrapper.inverse()
    assert_matrix(result.matrix,expected,rel=3e-14,abs=0)
    assert_matrix((wrapper*result).matrix,matrix.identity(size),rel=0,abs=3e-14)
    assert_matrix((result*wrapper).matrix,matrix.identity(size),rel=0,abs=3e-14)
    assert wrapper.matrix is rows and rows==original
    assert [list(row) for row in wrapper.c_matrix]==exported
    assert result is not wrapper and result.matrix is not rows
    for output,row in zip(result.c_matrix,result.matrix):
        for actual,value in zip(output,row):
            # Public exports remain float32, including finite binary64 overflow.
            import ctypes
            assert actual==ctypes.c_float(value).value
    receiver=wrapper
    assert wrapper.i_inverse() is receiver
    assert_matrix(wrapper.matrix,expected,rel=3e-14,abs=0)
    for output,row in zip(wrapper.c_matrix,wrapper.matrix):
        for actual,value in zip(output,row):
            import ctypes
            assert actual==ctypes.c_float(value).value


@pytest.mark.parametrize('size',[3,4])
def test_diagonal_mixed_exponents(size):
    diagonal=[1e-300,-1e300,1e-200,1e200][:size]
    rows=[[diagonal[i] if i==j else 0 for j in range(size)] for i in range(size)]
    result=matrix.Matrix(size,rows).inverse()
    assert_matrix(result.matrix,[[1/diagonal[i] if i==j else 0 for j in range(size)] for i in range(size)],rel=3e-15,abs=0)
    assert_matrix((matrix.Matrix(size,rows)*result).matrix,matrix.identity(size),abs=2e-15,rel=0)
    assert_matrix((result*matrix.Matrix(size,rows)).matrix,matrix.identity(size),abs=2e-15,rel=0)


@pytest.mark.parametrize('size',[3,4])
@pytest.mark.parametrize('scale',[1e-200,1,1e200])
def test_singular_matrices_keep_exception(size,scale):
    rows=[[float(i+j+1)*scale for j in range(size)] for i in range(size)]
    # Duplicate row gives an exact represented singularity.
    rows[1]=rows[0][:]
    with pytest.raises(ZeroDivisionError):matrix.Matrix(size,rows).inverse()
    with pytest.raises(ZeroDivisionError):matrix.Matrix(size).scale(V([0]*size)).inverse()


@pytest.mark.parametrize('size',[3,4])
def test_unrepresentable_inverse_coefficients(size):
    tiny=math.ldexp(1.,-1074)
    diagonal=[tiny,-tiny]+[1.]*(size-2)
    rows=[[diagonal[i] if i==j else 0 for j in range(size)] for i in range(size)]
    out=matrix.Matrix(size,rows).inverse().matrix
    assert out[0][0]==float('inf') and out[1][1]==-float('inf')
    assert all(out[i][j]==0 for i in range(size) for j in range(size) if i!=j)


def test_ordinary_singular_integer_cancellation():
    with pytest.raises(ZeroDivisionError):matrix.inverse3([[1,2,3],[4,5,6],[7,8,9]])


def test_degenerate_core_callers_retain_errors():
    z=V([0,0,0]);x=V([1,0,0])
    calls=[lambda:quaternion.quat_from_axis_angle(z,90),
           lambda:quaternion.quat_rotate_from_axis_angle([0,0,0],90),
           lambda:matrix.rotate3([0,0,0],90),lambda:matrix.rotate4([0,0,0],90),
           lambda:matrix.lookAt(z,z,V([0,1,0])),
           lambda:matrix.lookAt(z,x,z),lambda:matrix.lookAt(z,x,x),
           lambda:plane.Plane().fromPoints(z,x,V([2,0,0])),
           lambda:plane.Plane().bestFitNormal([]),lambda:ray.Ray(x,z),
           lambda:common.getViewPort(z,100,100)]
    for call in calls:
        with pytest.raises(ZeroDivisionError):call()
    for method,arg in [('roateUsingMatrix',matrix.Matrix(3,matrix.zero_matrix(3))),
                       ('rotateUsingQuaternion',quaternion.Quaternion([0,0,0,0]))]:
        r=ray.Ray(V([1,0,0]),V([1,0,0]))
        with pytest.raises(ZeroDivisionError):getattr(r,method)(arg)


@pytest.mark.parametrize('values',[[float('inf'),0,0],[float('inf'),float('nan'),0],[float('nan'),1,0]])
def test_legacy_nonfinite_norm_paths(values):
    expected=math.sqrt(sum(c*c for c in values))
    for data,cls in [(values,lambda a:V(a)),(values+[0],quaternion.Quaternion)]:
        obj=cls(data[:]);actual=obj.magnitude()
        assert math.isnan(actual) if math.isnan(expected) else actual==expected
        result=obj.normalize();output=result.vector if hasattr(result,'vector') else result.data
        for c,value in zip(data,output):
            reference=c/expected
            assert math.isnan(value) if math.isnan(reference) else value==reference
