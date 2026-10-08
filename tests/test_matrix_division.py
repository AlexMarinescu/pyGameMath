"""Phase 2B: elementwise division, explicit adjugate, and ctypes consistency."""
import pytest
from gem import matrix
from .helpers import inverse, assert_matrix


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('scalar', [2.0, -4.0, 0.5])
def test_nonsymmetric_scalar_division(size, scalar):
    values = [[i*size+j+1 for j in range(size)] for i in range(size)]
    original = [row[:] for row in values]
    expected = [[value/scalar for value in row] for row in values]
    source = matrix.Matrix(size, values)
    assert_matrix(matrix.matrix_div(values, scalar), expected)
    # The helper has always accepted integer divisors, unlike the wrapper.
    assert_matrix(matrix.matrix_div(values, 2), [[value/2 for value in row] for row in values])
    result = source / scalar
    legacy = source.__div__(scalar)
    for divided in [result, legacy]:
        assert isinstance(divided, matrix.Matrix)
        assert divided.size == size
        assert divided is not source
        assert_matrix(divided.matrix, expected)
        assert_matrix([list(row) for row in divided.c_matrix], expected, abs=1e-6)
        assert all(new is not old for new, old in zip(divided.matrix, values))
    assert values == original
    assert source.matrix == original
    alias = source
    source /= scalar
    assert source is alias
    assert_matrix(source.matrix, expected)
    assert_matrix([list(row) for row in source.c_matrix], expected, abs=1e-6)
    assert values == original
    legacy_inplace = matrix.Matrix(size, [row[:] for row in values])
    assert legacy_inplace.__idiv__(scalar) is legacy_inplace
    assert_matrix(legacy_inplace.matrix, expected)
    assert_matrix([list(row) for row in legacy_inplace.c_matrix], expected, abs=1e-6)


@pytest.mark.parametrize('values,expected', [
    ([[1,2],[3,4]], [[-2,1],[1.5,-0.5]]),
    ([[1,2,0],[0,1,3],[0,0,1]], [[1,-2,6],[0,1,-3],[0,0,1]]),
    ([[1,2,0,0],[0,1,3,0],[0,0,1,4],[0,0,0,1]],
     [[1,-2,6,-24],[0,1,-3,12],[0,0,1,-4],[0,0,0,1]]),
])
def test_nonsymmetric_inverse_and_division(values, expected):
    size = len(values)
    original = [row[:] for row in values]
    source = matrix.Matrix(size, values)
    result = source.inverse()
    assert_matrix(result.matrix, expected)
    assert_matrix(result.matrix, inverse(values))
    assert_matrix((source*result).matrix, matrix.identity(size), abs=1e-12)
    assert_matrix((result*source).matrix, matrix.identity(size), abs=1e-12)
    scaled = source / 2.0
    scaled_inverse = scaled.inverse()
    assert_matrix(scaled_inverse.matrix, [[2*x for x in row] for row in expected])
    assert_matrix((scaled*scaled_inverse).matrix, matrix.identity(size), abs=1e-12)
    assert_matrix((scaled_inverse*scaled).matrix, matrix.identity(size), abs=1e-12)
    assert source.i_inverse() is source
    assert_matrix(source.matrix, expected)
    assert_matrix([list(row) for row in source.c_matrix], expected, abs=1e-6)
    assert values == original


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('method', ['__div__', '__truediv__', '__idiv__', '__itruediv__'])
def test_zero_division_preserves_receiver(size, method):
    values = [[i*size+j+1 for j in range(size)] for i in range(size)]
    source = matrix.Matrix(size, values)
    snapshot = [list(row) for row in source.c_matrix]
    with pytest.raises(ZeroDivisionError):
        getattr(source, method)(0.0)
    assert source.matrix == values
    assert [list(row) for row in source.c_matrix] == snapshot


@pytest.mark.parametrize('method', ['__div__', '__truediv__', '__idiv__', '__itruediv__'])
@pytest.mark.parametrize('operand', [2, '2', None, object()])
def test_division_unsupported_operand(method, operand):
    source = matrix.Matrix(2)
    assert getattr(source, method)(operand) is NotImplemented
    assert source.matrix == matrix.identity(2)
    assert [list(row) for row in source.c_matrix] == matrix.identity(2)
