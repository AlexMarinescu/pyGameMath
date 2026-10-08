import copy
import ctypes
import random

import pytest

from gem import matrix, vector
from .helpers import assert_matrix


def V(*values):
    return vector.Vector(len(values), list(values))


def assert_export(m):
    expected = [[ctypes.c_float(value).value for value in row] for row in m.matrix]
    assert [list(row) for row in m.c_matrix] == expected


@pytest.mark.parametrize('offset', [(2, -3, 4), (-5, 6, -7), (0, 0, 0)])
@pytest.mark.parametrize('homogeneous_offset', [False, True])
@pytest.mark.parametrize('rows', [
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    [[2, 0, 0, 0], [0, -3, 0, 0], [0, 0, 4, 0], [5, -6, 7, 1]],
    [[1, 2, 3, 4], [5, 6, 7, 8], [-2, 3, -4, 5], [6, -7, 8, 9]],
])
def test_translation4_known_answers(rows, homogeneous_offset, offset):
    source = copy.deepcopy(rows)
    original = copy.deepcopy(rows)
    m = matrix.Matrix(4, source)
    old_export = [list(row) for row in m.c_matrix]
    # The fourth offset component is historically ignored.
    displacement = V(*(offset + (37,) if homogeneous_offset else offset))
    expected = [[row[j] + row[3] * offset[j] for j in range(3)] + [row[3]]
                for row in original]
    result = m.translate(displacement)
    assert result is not m
    assert_matrix(result.matrix, expected)
    assert m.matrix == source == original
    assert [list(row) for row in m.c_matrix] == old_export
    assert_export(result)
    assert m.i_translate(displacement) is m
    assert_matrix(m.matrix, expected)
    assert_export(m)
    assert source == original
    assert displacement.vector == list(offset) + ([37] if homogeneous_offset else [])


@pytest.mark.parametrize('seed', range(12))
def test_translation4_seeded_properties(seed):
    rng = random.Random(2200 + seed)
    offset = [rng.uniform(-8, 8) for _ in range(3)]
    rows = [[rng.uniform(-4, 4) for _ in range(4)] for _ in range(4)]
    expected = [[row[j] + row[3] * offset[j] for j in range(3)] + [row[3]]
                for row in rows]
    m = matrix.Matrix(4, copy.deepcopy(rows))
    assert_matrix(m.translate(V(*offset)).matrix, expected)
    m.i_translate(V(*offset))
    assert_matrix(m.matrix, expected)
    assert_export(m)
    m.i_translate(V(*[-value for value in offset]))
    assert_matrix(m.matrix, rows, abs=1e-13)
    assert_export(m)


@pytest.mark.parametrize('values,offset,expected', [
    ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (2, -3),
     [[1, 0, 0], [0, 1, 0], [2, -3, 1]]),
    ([[2, 1, 0], [-3, 4, 0], [5, 6, 1]], (-2, 3),
     [[2, 1, 0], [-3, 4, 0], [3, 9, 1]]),
    ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (2, -3, 4),
     [[1, 0, 0], [0, 1, 0], [2, -3, 4]]),
])
def test_translation3_existing_variants(values, offset, expected):
    # Vector3 offsets retain the legacy translate3 semantics, not affine 3D translation.
    m = matrix.Matrix(3, copy.deepcopy(values))
    assert_matrix(m.translate(V(*offset)).matrix, expected)
    assert m.matrix == values
    assert m.i_translate(V(*offset)) is m
    assert_matrix(m.matrix, expected)
    assert_export(m)


def test_translation2_remains_unsupported():
    m = matrix.Matrix(2)
    for method in [m.translate, m.i_translate]:
        with pytest.raises(NotImplementedError):
            method(V(2, -3))
    assert m.matrix == [[1, 0], [0, 1]]
    assert_export(m)


@pytest.mark.parametrize('size', [1, 2, 3, 4, 8])
def test_transform_identity_and_ownership(size):
    position = [i - 3.5 for i in range(size)]
    original = position[:]
    rows = matrix.identity(size)
    receiver = vector.Vector(size, [99] * size)
    previous = receiver.vector
    raw = vector.transform(size, position, rows)
    assert raw == original
    assert raw is not position
    result = receiver.transform(position, rows)
    assert result is not receiver
    assert result.size == size
    assert result.vector == original
    assert receiver.vector is previous
    assert previous == [99] * size
    assert receiver.i_transform(position, rows) is receiver
    assert receiver.vector == original
    assert receiver.vector is not position
    assert position == original
    assert rows == matrix.identity(size)


@pytest.mark.parametrize('size', [1, 2, 3, 4, 8])
def test_transform_nonsymmetric_basis_oracle(size):
    rows = [[10 * i + j + 1 for j in range(size)] for i in range(size)]
    original = copy.deepcopy(rows)
    for i in range(size):
        position = [0] * size
        position[i] = 1
        # A row basis vector selects the corresponding matrix row.
        assert vector.transform(size, position, rows) == rows[i]
    assert vector.transform(size, [0] * size, rows) == [0] * size
    assert rows == original


@pytest.mark.parametrize('size', [1, 2, 3, 4, 8])
@pytest.mark.parametrize('seed', range(8))
def test_transform_linear_properties(size, seed):
    rng = random.Random(2300 + 10 * size + seed)
    a = [rng.randint(-5, 5) for _ in range(size)]
    b = [rng.randint(-5, 5) for _ in range(size)]
    rows = [[rng.randint(-4, 4) for _ in range(size)] for _ in range(size)]
    original = copy.deepcopy(rows)
    # Independent sum of scaled row basis images, using exact integer arithmetic.
    expected = [sum(value * row[j] for value, row in zip(a, rows))
                for j in range(size)]
    assert vector.transform(size, a, rows) == expected
    ta, tb = vector.transform(size, a, rows), vector.transform(size, b, rows)
    assert vector.transform(size, [x + y for x, y in zip(a, b)], rows) == [
        x + y for x, y in zip(ta, tb)]
    assert vector.transform(size, [3 * x for x in a], rows) == [3 * x for x in ta]
    assert rows == original


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('angle,expected_xy', [
    (0, (2, 3)), (90, (-3, 2)), (-90, (3, -2)), (180, (-2, -3)),
])
def test_transform_known_rotations(size, angle, expected_xy):
    axis = V(0, 0) if size == 2 else V(0, 0, 1)
    m = matrix.Matrix(size).rotate(axis, angle)
    position = [2, 3] + ([4] if size >= 3 else []) + ([1] if size == 4 else [])
    expected = list(expected_xy) + position[2:]
    assert vector.transform(size, position, m.matrix) == pytest.approx(expected, abs=1e-14)
    assert position == [2, 3] + ([4] if size >= 3 else []) + ([1] if size == 4 else [])


@pytest.mark.parametrize('size', [2, 3, 4])
def test_inplace_transform_receiver_alias(size):
    # Reading self.vector must finish before its replacement.
    m = matrix.Matrix(size, [[i * size + j + 1 for j in range(size)] for i in range(size)])
    receiver = vector.Vector(size, [0] * size)
    receiver.vector[0] = 1
    old_values = receiver.vector
    assert receiver.i_transform(receiver.vector, m.matrix) is receiver
    assert receiver.vector == m.matrix[0]
    assert receiver.vector is not old_values
    assert old_values == [1] + [0] * (size - 1)


@pytest.mark.parametrize('size', [1, 2, 3, 4, 8])
def test_transform_implicit_homogeneous_identity(size):
    position = [i - 2 for i in range(size)]
    assert vector.transform(size, position, matrix.identity(size + 1)) == position


@pytest.mark.parametrize('size', [2, 3])
@pytest.mark.parametrize('offset_sign', [-1, 1])
def test_transform_translation_and_homogeneous_weights(size, offset_sign):
    offset = [offset_sign * value for value in [2, -3, 4][:size]]
    position = [1, 2, -5][:size]
    m = matrix.Matrix(size + 1).translate(V(*offset))
    original = copy.deepcopy(m.matrix)
    expected = [value + delta for value, delta in zip(position, offset)]
    assert vector.transform(size, position, m.matrix) == expected
    receiver = vector.Vector(size)
    assert receiver.transform(position, m.matrix).vector == expected
    assert receiver.vector == [0] * size
    assert receiver.i_transform(position, m.matrix) is receiver
    assert receiver.vector == expected
    for w in [0, 1, -2, 0.5]:
        homogeneous = position + [w]
        explicit = [value + w * delta for value, delta in zip(position, offset)] + [w]
        assert vector.transform(size + 1, homogeneous, m.matrix) == explicit
        assert vector.Vector(size + 1).transform(homogeneous, m.matrix).vector == explicit
        in_place = vector.Vector(size + 1)
        assert in_place.i_transform(homogeneous, m.matrix) is in_place
        assert in_place.vector == explicit
        assert homogeneous == position + [w]
    assert position == [1, 2, -5][:size]
    assert m.matrix == original
    assert_export(m)


@pytest.mark.parametrize('size', [2, 3])
def test_noncommuting_translation_rotation(size):
    offset = [2, -3, 4][:size]
    # Use a literal quarter-turn to avoid the unrelated pivot-rotation API.
    rotation_rows = matrix.identity(size + 1)
    rotation_rows[0][:2] = [0, 1]
    rotation_rows[1][:2] = [-1, 0]
    rotation = matrix.Matrix(size + 1, rotation_rows)
    translation = matrix.Matrix(size + 1).translate(V(*offset))
    position = [1, 2, 3][:size]
    translated_then_rotated = [1, 3] + ([7] if size == 3 else [])
    rotated_then_translated = [0, -2] + ([7] if size == 3 else [])
    tr, rt = translation * rotation, rotation * translation
    assert vector.transform(size, position, tr.matrix) == translated_then_rotated
    assert vector.transform(size, position, rt.matrix) == rotated_then_translated
    assert vector.transform(size, vector.transform(size, position, translation.matrix),
                            rotation.matrix) == translated_then_rotated
    assert vector.transform(size, vector.transform(size, position, rotation.matrix),
                            translation.matrix) == rotated_then_translated
    for m, expected in [(tr, translated_then_rotated), (rt, rotated_then_translated)]:
        assert vector.transform(size + 1, position + [1], m.matrix) == expected + [1]
        assert (m * V(*(position + [1]))).vector == expected + [1]
        assert_export(m)
    in_place = matrix.Matrix(size + 1, copy.deepcopy(rotation.matrix))
    assert in_place.i_translate(V(*offset)) is in_place
    assert_matrix(in_place.matrix, rt.matrix)
    assert_export(in_place)
    assert position == [1, 2, 3][:size]


@pytest.mark.parametrize('seed', range(12))
def test_transform_affine_analytic_properties(seed):
    rng = random.Random(2400 + seed)
    x, y, z = [rng.randint(-5, 5) for _ in range(3)]
    sx, sy, sz = [rng.randint(1, 5) for _ in range(3)]
    tx, ty, tz = [rng.randint(-5, 5) for _ in range(3)]
    # Scale, quarter-turn about +Z, then translate, derived analytically.
    rows = [[0, sx, 0, 0], [-sy, 0, 0, 0], [0, 0, sz, 0], [tx, ty, tz, 1]]
    position = [x, y, z]
    expected = [-sy * y + tx, sx * x + ty, sz * z + tz]
    result = vector.transform(3, position, rows)
    assert result == expected
    assert vector.transform(4, position + [1], rows) == expected + [1]
    assert vector.transform(4, position + [0], rows) == [-sy * y, sx * x, sz * z, 0]
    midpoint = [value / 2 for value in position]
    assert vector.transform(3, midpoint, rows) == [
        (value + delta) / 2 for value, delta in zip(result, [tx, ty, tz])]
    assert position == [x, y, z]


@pytest.mark.parametrize('offset_w', [0, 2, -4, -6])
def test_transform_computes_projective_w_without_division(offset_w):
    rows = [[1, 0, 0, 2], [0, 1, 0, 0], [0, 0, 1, 0], [5, -6, 7, offset_w]]
    original = copy.deepcopy(rows)
    for x in [0, 2, -1]:
        for w in [0, 1, 2, -0.5]:
            position = [x, 3, 4, w]
            expected = [x + 5 * w, 3 - 6 * w, 4 + 7 * w, 2 * x + offset_w * w]
            assert vector.transform(4, position, rows) == expected
            receiver = vector.Vector(4)
            assert receiver.transform(position, rows).vector == expected
            assert receiver.i_transform(position, rows) is receiver
            assert receiver.vector == expected
            assert position == [x, 3, 4, w]
        # Implicit coordinates omit computed output w, even if zero/nonunit.
        assert vector.transform(3, [x, 3, 4], rows) == [x + 5, -3, 11]
    assert rows == original


def test_general_matrix_operator_does_not_promote_vector3():
    m = matrix.Matrix(4).translate(V(2, -3, 4))
    position = V(1, 2, 3)
    # Preserve the existing operator's dimension precondition and failure.
    with pytest.raises(IndexError):
        m * position
    assert position.vector == [1, 2, 3]
    assert_export(m)
