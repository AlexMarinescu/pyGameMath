import copy
import ctypes
import math
import random

import pytest

from gem import matrix, vector
from .helpers import assert_matrix


def V(*values):
    return vector.Vector(len(values), list(values))


def assert_export(m):
    assert [list(row) for row in m.c_matrix] == [
        [ctypes.c_float(value).value for value in row] for row in m.matrix]


@pytest.mark.parametrize('pivot', [(2, 3), (-4, 1), (0, 0)])
@pytest.mark.parametrize('angle,offset', [
    (0, (1, 2)), (90, (-2, 1)), (-90, (2, -1)), (180, (-1, -2)),
])
def test_pivot_rotation_known_answers(pivot, angle, offset):
    point = list(pivot)
    rows = matrix.rotate2(point, angle)
    assert len(rows) == 3 and all(len(row) == 3 for row in rows)
    r = matrix.Matrix(3, rows)
    assert (r * V(*pivot, 1)).vector == pytest.approx([*pivot, 1], abs=1e-14)
    position = V(pivot[0] + 1, pivot[1] + 2, 1)
    assert (r * position).vector == pytest.approx(
        [pivot[0] + offset[0], pivot[1] + offset[1], 1], abs=1e-14)
    assert position.vector == [pivot[0] + 1, pivot[1] + 2, 1]
    # Homogeneous directions rotate without the pivot's affine offset.
    assert (r * V(1, 2, 0)).vector == pytest.approx([*offset, 0], abs=1e-14)
    assert point == list(pivot)
    assert_export(r)


def test_pivot_rotation_eighth_turn_known_answer():
    r = matrix.Matrix(3, matrix.rotate2([2, -3], 45))
    assert (r * V(4, -3, 1)).vector == pytest.approx(
        [2 + math.sqrt(2), -3 + math.sqrt(2), 1], abs=1e-14)


@pytest.mark.parametrize('seed', range(10))
def test_pivot_rotation_analytic_properties(seed):
    rng = random.Random(2700 + seed)
    px, py, x, y = [rng.uniform(-6, 6) for _ in range(4)]
    angle = rng.uniform(-180, 180)
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    # Rotate the relative point before adding back the pivot.
    expected = [px + (x - px) * c - (y - py) * s,
                py + (x - px) * s + (y - py) * c, 1]
    rows = matrix.rotate2([px, py], angle)
    r = matrix.Matrix(3, rows)
    result = (r * V(x, y, 1)).vector
    assert result == pytest.approx(expected, abs=1e-13)
    assert math.hypot(result[0] - px, result[1] - py) == pytest.approx(math.hypot(x - px, y - py))
    assert (r * V(px, py, 1)).vector == pytest.approx([px, py, 1], abs=1e-13)
    reverse = matrix.Matrix(3, matrix.rotate2([px, py], -angle))
    assert_matrix((r * reverse).matrix, matrix.identity(3), abs=1e-13)


def test_pivot_rotation_noncommuting_translation_and_inplace_composition():
    r = matrix.Matrix(3, matrix.rotate2([2, 3], 90))
    t = matrix.Matrix(3).translate(V(4, -1))
    point = V(3, 3, 1)
    # Translate first: (7,2) rotates about (2,3) to (3,8).
    assert ((t * r) * point).vector == pytest.approx([3, 8, 1], abs=1e-14)
    # Rotate first: (2,4), then translate to (6,3).
    assert ((r * t) * point).vector == pytest.approx([6, 3, 1], abs=1e-14)
    for a, b in [(t, r), (r, t)]:
        original = copy.deepcopy(a.matrix)
        result = a * b
        in_place = matrix.Matrix(3, copy.deepcopy(a.matrix))
        assert in_place.__imul__(b) is in_place
        assert_matrix(in_place.matrix, result.matrix)
        assert a.matrix == original
        assert_export(result)
        assert_export(in_place)
    assert point.vector == [3, 3, 1]


@pytest.mark.parametrize('angle,expected', [(90, (-3, 2)), (-90, (3, -2))])
def test_matrix_rotation_dispatch_preserved(angle, expected):
    # Matrix2 cannot encode a translated pivot; its historical argument is ignored.
    axis2 = V(5, -7)
    m = matrix.Matrix(2)
    result = m.rotate(axis2, angle)
    assert (result * V(2, 3)).vector == pytest.approx(expected, abs=1e-14)
    assert m.i_rotate(axis2, angle) is m
    assert_matrix(m.matrix, result.matrix)
    assert axis2.vector == [5, -7]
    assert_export(m)
    for size in [3, 4]:
        axis = V(0, 0, 1)
        m = matrix.Matrix(size)
        result = m.rotate(axis, angle)
        values = [2, 3, 4] + ([1] if size == 4 else [])
        assert (result * V(*values)).vector == pytest.approx([*expected, *values[2:]], abs=1e-14)
        assert m.i_rotate(axis, angle) is m
        assert_matrix(m.matrix, result.matrix)
        assert axis.vector == [0, 0, 1]
        assert_export(m)


def sheared(values, plane, a, b):
    # Independent coordinate mapping, not matrix-index arithmetic.
    x, y, z = values[:3]
    if plane == 'XY':
        return [x, y, z + a * x + b * y] + values[3:]
    if plane == 'YZ':
        return [x + a * y + b * z, y, z] + values[3:]
    return [x, y + a * x + b * z, z] + values[3:]


@pytest.mark.parametrize('plane', ['XY', 'YZ', 'XZ'])
@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('factors', [(0, 0), (1, 2), (-2, 3), (0.5, -0.25), (-1, -2)])
@pytest.mark.parametrize('nonidentity', [False, True])
def test_shear_coordinate_oracles_and_methods(plane, size, factors, nonidentity):
    a, b = factors
    helper = getattr(matrix, 'shear' + plane + str(size))
    rows = helper(a, b)
    # The images of row basis vectors independently specify the matrix.
    expected_rows = [sheared(row, plane, a, b) for row in matrix.identity(size)]
    assert rows == expected_rows
    assert len(rows) == size and all(len(row) == size for row in rows)
    transform = matrix.Matrix(size, rows)
    for w in ([0, 1, -2, 0.5] if size == 4 else [None]):
        values = [2, -3, 4] + ([] if w is None else [w])
        v = V(*values)
        assert (transform * v).vector == sheared(values, plane, a, b)
        assert v.vector == values
    source = ([[2, -1, 3], [4, 5, -2], [-3, 2, 6]] if size == 3 else
              [[2, -1, 3, 4], [4, 5, -2, -1], [-3, 2, 6, 2], [7, -8, 9, 3]])
    if not nonidentity:
        source = matrix.identity(size)
    original = copy.deepcopy(source)
    m = matrix.Matrix(size, source)
    expected = [sheared(row, plane, a, b) for row in source]
    old_export = [list(row) for row in m.c_matrix]
    result = getattr(m, 'shear' + plane)(a, b)
    assert result is not m
    assert_matrix(result.matrix, expected)
    assert m.matrix == source == original
    assert [list(row) for row in m.c_matrix] == old_export
    assert_export(result)
    assert getattr(m, 'i_shear' + plane)(a, b) is m
    assert_matrix(m.matrix, expected)
    assert_export(m)
    assert source == original


@pytest.mark.parametrize('size', [3, 4])
def test_shear_planes_do_not_commute(size):
    xy = matrix.Matrix(size).shearXY(1, 2)
    yz = matrix.Matrix(size).shearYZ(3, -1)
    values = [2, 3, 4] + ([1] if size == 4 else [])
    suffix = values[3:]
    assert ((xy * yz) * V(*values)).vector == [-1, 3, 12] + suffix
    assert ((yz * xy) * V(*values)).vector == [7, 3, 17] + suffix
    result = matrix.Matrix(size).shearXY(1, 2)
    assert result.i_shearYZ(3, -1) is result
    assert_matrix(result.matrix, (xy * yz).matrix)
    assert_export(result)


def test_shear_translation_composition_order():
    t = matrix.Matrix(4).translate(V(2, -3, 4))
    s = matrix.Matrix(4).shearXY(1, 2)
    point = V(1, 2, 3, 1)
    assert ((t * s) * point).vector == [3, -1, 8, 1]
    assert ((s * t) * point).vector == [3, -1, 12, 1]
    assert point.vector == [1, 2, 3, 1]


@pytest.mark.parametrize('seed', range(6))
def test_shear_seeded_properties(seed):
    rng = random.Random(2800 + seed)
    a, b = rng.uniform(-3, 3), rng.uniform(-3, 3)
    values = [rng.uniform(-6, 6) for _ in range(4)]
    for plane in ['XY', 'YZ', 'XZ']:
        for size in [3, 4]:
            m = getattr(matrix.Matrix(size), 'shear' + plane)(a, b)
            expected = sheared(values[:size], plane, a, b)
            assert (m * V(*values[:size])).vector == pytest.approx(expected, abs=1e-13)
            assert m.det() == pytest.approx(1, abs=1e-13)
            undo = getattr(matrix.Matrix(size), 'shear' + plane)(-a, -b)
            assert_matrix((m * undo).matrix, matrix.identity(size), abs=1e-13)
            assert_export(m)
