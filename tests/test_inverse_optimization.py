"""Independent Fraction/integer references for the optimized scaled inverses."""
import ctypes
from fractions import Fraction
import math
import random

import pytest

from gem import matrix
from .helpers import assert_matrix, determinant, inverse


@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('seed', range(12))
@pytest.mark.parametrize('scale', [1e-300, 1e-150, 1, 1e150, 1e300])
def test_well_conditioned_fraction_reference(size, seed, scale):
    rng = random.Random(seed)
    rows = [[((4 if i == j else 0) + rng.uniform(-.5, .5))*scale
             for j in range(size)] for i in range(size)]
    original = [row[:] for row in rows]
    expected = inverse(rows)  # Exact rational elimination of represented inputs.
    receiver = matrix.Matrix(size, rows)
    raw = getattr(matrix, 'inverse'+str(size))(rows)
    assert_matrix(raw, expected, rel=3e-14, abs=0)
    for result in [raw, receiver.inverse().matrix]:
        assert_matrix(matrix.matrix_multiply(rows, result), matrix.identity(size), rel=0, abs=3e-14)
        assert_matrix(matrix.matrix_multiply(result, rows), matrix.identity(size), rel=0, abs=3e-14)
    returned = receiver.inverse()
    assert returned is not receiver and returned.matrix is not rows
    assert rows == original and receiver.matrix is rows
    assert receiver.i_inverse() is receiver
    assert rows == original and receiver.matrix is not rows
    assert_matrix(receiver.matrix, expected, rel=3e-14, abs=0)
    for obj in [returned, receiver]:
        for exported, row in zip(obj.c_matrix, obj.matrix):
            assert list(exported) == [ctypes.c_float(v).value for v in row]


@pytest.mark.parametrize('seed', range(30))
def test_reduced_integer_determinant_exact(seed):
    rng = random.Random(seed)
    rows = [[rng.randint(-20, 20) << rng.randrange(0, 1100) for _ in range(4)] for _ in range(4)]
    assert matrix._det4_exact(rows) == determinant(rows)
    rows[-1] = rows[0][:]
    assert matrix._det4_exact(rows) == determinant(rows) == 0


@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('exponent', [-900, 0, 900])
def test_exact_singular_linear_combination(size, exponent):
    rows = [[float(i == j) for j in range(size)] for i in range(size)]
    rows[-1] = [rows[0][j] + 2*rows[1][j] for j in range(size)]
    rows = [[math.ldexp(v, exponent) for v in row] for row in rows]
    assert determinant([[Fraction(v) for v in row] for row in rows]) == 0
    with pytest.raises(ZeroDivisionError):getattr(matrix, 'inverse'+str(size))(rows)


@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('exponent', [-900, 0, 900])
def test_near_singular_binary_known_answer(size, exponent):
    rows = matrix.identity(size)
    rows[0][:2] = [1., 1.]
    rows[1][:2] = [1., 1. + 2**-30]
    rows = [[math.ldexp(v, exponent) for v in row] for row in rows]
    expected = inverse(rows)
    actual = getattr(matrix, 'inverse'+str(size))(rows)
    assert_matrix(actual, expected, rel=3e-14, abs=0)
    for product in [matrix.matrix_multiply(rows, actual), matrix.matrix_multiply(actual, rows)]:
        assert_matrix(product, matrix.identity(size), rel=0, abs=3e-14)


@pytest.mark.parametrize('size', [3, 4])
def test_mixed_row_exponents_nonsymmetric(size):
    exponents = [-300, 300, -200, 200][:size]
    rows = [[math.ldexp(float(2 if i==j else 1 if j==i+1 else 0), exponents[i])
             for j in range(size)] for i in range(size)]
    expected = [[math.ldexp((-1)**(j-i)/2**(j-i+1), -exponents[j]) if j>=i else 0
                 for j in range(size)] for i in range(size)]
    actual = getattr(matrix, 'inverse'+str(size))(rows)
    assert_matrix(actual, expected, rel=3e-14, abs=0)
    for product in [matrix.matrix_multiply(rows, actual), matrix.matrix_multiply(actual, rows)]:
        assert_matrix(product, matrix.identity(size), rel=0, abs=3e-14)


@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('value,error', [(1j, TypeError), ('x', TypeError), (10**500, OverflowError)])
def test_invalid_values_keep_exception_class(size, value, error):
    rows = matrix.identity(size);rows[0][0] = value
    with pytest.raises(error):getattr(matrix, 'inverse'+str(size))(rows)


@pytest.mark.parametrize('size', [3, 4])
def test_short_row_keeps_index_error(size):
    rows = matrix.identity(size);rows[0] = [1.]
    with pytest.raises(IndexError):getattr(matrix, 'inverse'+str(size))(rows)


@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_legacy_kernel_output(size, value):
    rows = matrix.identity(size);rows[0][0] = value
    result = getattr(matrix, 'inverse'+str(size))(rows)
    for i, row in enumerate(result):
        for j, actual in enumerate(row):
            if math.isnan(value) or (i > 0 and j > 0):assert math.isnan(actual)
            else:
                assert actual == 0
                # Legacy division gives reciprocal sign at [0,0].
                if i == j == 0:assert math.copysign(1, actual) == math.copysign(1, value)
