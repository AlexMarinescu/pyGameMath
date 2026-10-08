import copy
import ctypes
import math
import random

import pytest

from gem import matrix, quaternion, vector
from .helpers import assert_matrix, determinant


def V(*values):
    return vector.Vector(len(values), list(values))


def assert_equivalent(actual, expected):
    assert isinstance(actual, quaternion.Quaternion)
    assert len(actual.data) == 4
    dot = sum(a * b for a, b in zip(actual.data, expected))
    assert abs(dot) == pytest.approx(1, abs=1e-13)
    sign = 1 if dot >= 0 else -1
    assert actual.data == pytest.approx([sign * value for value in expected], abs=1e-13)


def assert_export(m):
    assert isinstance(m, matrix.Matrix) and m.size == 4
    assert [list(row) for row in m.c_matrix] == [
        [ctypes.c_float(value).value for value in row] for row in m.matrix]


def rodrigues(axis, angle, point):
    # Independent geometric vector rotation: c*v + s*(n cross v) + (1-c)*n*(n dot v).
    length = math.sqrt(sum(value * value for value in axis))
    nx, ny, nz = [value / length for value in axis]
    x, y, z = point
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    cross = [ny * z - nz * y, nz * x - nx * z, nx * y - ny * x]
    dot = nx * x + ny * y + nz * z
    return [c * value + s * perpendicular + (1 - c) * n * dot
            for value, perpendicular, n in zip(point, cross, [nx, ny, nz])]


def rotation_rows(axis, angle):
    return [rodrigues(axis, angle, basis) + [0] for basis in
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]] + [[0, 0, 0, 1]]


SQRT_HALF = math.sqrt(0.5)
@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('rows,expected', [
    ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 0, 0, 0]),
    ([[1, 0, 0], [0, 0, 1], [0, -1, 0]], [SQRT_HALF, SQRT_HALF, 0, 0]),
    ([[0, 0, -1], [0, 1, 0], [1, 0, 0]], [SQRT_HALF, 0, SQRT_HALF, 0]),
    ([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], [SQRT_HALF, 0, 0, SQRT_HALF]),
    ([[0, -1, 0], [1, 0, 0], [0, 0, 1]], [SQRT_HALF, 0, 0, -SQRT_HALF]),
    ([[1, 0, 0], [0, -1, 0], [0, 0, -1]], [0, 1, 0, 0]),
    ([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], [0, 0, 1, 0]),
    ([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], [0, 0, 0, 1]),
    ([[-6 / 7, 2 / 7, 3 / 7], [2 / 7, -3 / 7, 6 / 7], [3 / 7, 6 / 7, 2 / 7]],
     [0, 1 / math.sqrt(14), 2 / math.sqrt(14), 3 / math.sqrt(14)]),
    # 120 degrees about (1,1,1): X -> Y -> Z -> X.
    ([[0, 1, 0], [0, 0, 1], [1, 0, 0]], [0.5, 0.5, 0.5, 0.5]),
])
def test_literal_rotation_known_answers(rows, expected, size):
    source = copy.deepcopy(rows) if size == 3 else [row[:] + [0] for row in rows] + [[0, 0, 0, 1]]
    m = matrix.Matrix(size, source)
    original = copy.deepcopy(source)
    old_export = [list(row) for row in m.c_matrix]
    recovered = quaternion.quat_from_matrix(m)
    assert_equivalent(recovered, expected)
    assert recovered.magnitude() == pytest.approx(1, abs=1e-13)
    assert source == original
    assert [list(row) for row in m.c_matrix] == old_export
    for values in [expected, [-value for value in expected]]:
        rot = quaternion.Quaternion(list(values))
        result = rot.toMatrix()
        assert_matrix([row[:3] for row in result.matrix[:3]], rows, abs=1e-13)
        assert_export(result)
        assert rot.data == values


@pytest.mark.parametrize('dominant,components', [
    ('w', [4, 1, 2, 3]), ('x', [1, 4, 2, 3]),
    ('y', [1, 2, 4, 3]), ('z', [1, 2, 3, 4]),
])
def test_each_largest_component_branch(dominant, components):
    # All four components are nonzero; X can exceed W while Y or Z is still larger.
    norm = math.sqrt(sum(value * value for value in components))
    expected = [value / norm for value in components]
    angle = math.degrees(2 * math.acos(expected[0]))
    rows = rotation_rows(expected[1:], angle)
    assert 'wxyz'[max(range(4), key=lambda i: abs(expected[i]))] == dominant
    recovered = quaternion.quat_from_matrix(matrix.Matrix(4, rows))
    assert_equivalent(recovered, expected)
    assert_matrix(recovered.toMatrix().matrix, rows, abs=1e-13)


@pytest.mark.parametrize('axis', [[3, 1, 2], [1, 3, 2], [1, 2, 3], [-2, 1, -3]])
@pytest.mark.parametrize('angle', [35, 90, 151, 180, -180, 179.999999, -151])
def test_arbitrary_axis_round_trips_and_rotation(axis, angle):
    length = math.sqrt(sum(value * value for value in axis))
    half = math.radians(angle) / 2
    expected = [math.cos(half)] + [value / length * math.sin(half) for value in axis]
    rows = rotation_rows(axis, angle)
    original = copy.deepcopy(rows)
    m = matrix.Matrix(4, rows)
    recovered = quaternion.quat_from_matrix(m)
    assert_equivalent(recovered, expected)
    rot = quaternion.Quaternion(expected[:])
    converted = rot.toMatrix()
    assert_matrix(converted.matrix, rows, abs=1e-13)
    assert_equivalent(quaternion.quat_from_matrix(converted), expected)
    assert_matrix(recovered.toMatrix().matrix, rows, abs=1e-13)
    assert_matrix(rot.negate().toMatrix().matrix, rows, abs=1e-13)
    for i in range(3):
        for j in range(3):
            assert sum(converted.matrix[i][k] * converted.matrix[j][k]
                       for k in range(3)) == pytest.approx(1 if i == j else 0, abs=1e-13)
    assert determinant([row[:3] for row in converted.matrix[:3]]) == pytest.approx(1, abs=1e-13)
    assert converted.det() == pytest.approx(1, abs=1e-13)
    point = V(2, -3, 4)
    answer = rodrigues(axis, angle, point.vector)
    assert quaternion.quat_rotate_vector(rot, point).vector == pytest.approx(answer, abs=1e-13)
    assert (converted * V(2, -3, 4, 0)).vector == pytest.approx(answer + [0], abs=1e-13)
    assert rows == original and rot.data == expected and point.vector == [2, -3, 4]
    assert_export(converted)
    assert_export(recovered.toMatrix())


@pytest.mark.parametrize('seed', range(12))
def test_seeded_sign_equivalence_and_export_ownership(seed):
    rng = random.Random(2900 + seed)
    components = [rng.uniform(-2, 2) for _ in range(4)]
    length = math.sqrt(sum(value * value for value in components))
    values = [value / length for value in components]
    source = values[:]
    rot = quaternion.Quaternion(source)
    a = quaternion.quat_to_matrix(rot)
    b = rot.toMatrix()
    assert a is not b and a.matrix is not b.matrix and a.c_matrix is not b.c_matrix
    assert_matrix(a.matrix, b.matrix)
    assert_export(a)
    assert_export(b)
    assert_equivalent(quaternion.quat_from_matrix(a), values)
    assert_equivalent(quaternion.quat_from_matrix(rot.negate().toMatrix()), values)
    assert source == values
    original = copy.deepcopy(b.matrix)
    a.matrix[0][0] = 99
    a.c_matrix[1][1] = 88
    assert b.matrix == original
    assert_export(b)
    assert source == values


@pytest.mark.parametrize('values,rows', [
    ([0, 0, 0, 0], [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    ([2, 0, 0, 1], [[-1, 4, 0, 0], [-4, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
])
def test_to_matrix_retains_nonunit_and_zero_numerics(values, rows):
    # Characterize unchanged formulas; these are not promises of proper rotations.
    rot = quaternion.Quaternion(values[:])
    result = rot.toMatrix()
    assert_matrix(result.matrix, rows)
    assert_export(result)
    assert rot.data == values


def test_from_matrix_preserves_existing_rotation_block_extraction():
    rows = [[0, 1, 0, 4], [-1, 0, 0, -5], [0, 0, 1, 6], [7, -8, 9, 2]]
    original = copy.deepcopy(rows)
    m = matrix.Matrix(4, rows)
    old_export = [list(row) for row in m.c_matrix]
    recovered = quaternion.quat_from_matrix(m)
    assert_equivalent(recovered, [SQRT_HALF, 0, 0, SQRT_HALF])
    assert rows == original and [list(row) for row in m.c_matrix] == old_export
