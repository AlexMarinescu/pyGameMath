import copy
import random

import pytest

from gem import matrix, vector
from .helpers import inverse


def V(*values):
    return vector.Vector(len(values), list(values))


def argument(rows, wrapped):
    return matrix.Matrix(4, copy.deepcopy(rows)) if wrapped else copy.deepcopy(rows)


def snapshot(value):
    if isinstance(value, matrix.Matrix):
        return (copy.deepcopy(value.matrix), [list(row) for row in value.c_matrix])
    return copy.deepcopy(value)


def assert_vector3(result, expected):
    assert isinstance(result, vector.Vector)
    assert result.size == len(result.vector) == 3
    assert result.vector == pytest.approx(expected, abs=1e-12)


# Literal windows are derived from the viewport and analytic camera equations.
# They are used as unprojection inputs without calling project first.
@pytest.mark.parametrize('model_wrapped,proj_wrapped', [
    (False, False), (False, True), (True, False), (True, True),
])
@pytest.mark.parametrize('kind,position,viewport,window', [
    ('identity', (0, 0, 0, 1), (0, 0, 100, 100), (50, 50, 0.5)),
    ('identity', (-1, 1, -1, 1), (-20, 30, 320, 160), (-20, 190, 0)),
    ('identity', (1, -1, 1, 1), (-20, 30, 320, 160), (300, 30, 1)),
    ('identity', (2, -2, 3, 1), (-20, 30, 320, 160), (460, -50, 2)),
    ('orthographic', (-2, -3, -2, 1), (-20, 30, 320, 160), (-20, 30, 0)),
    ('orthographic', (6, 5, -10, 1), (-20, 30, 320, 160), (300, 190, 1)),
    ('orthographic', (2, 1, -6, 1), (-20, 30, 320, 160), (140, 110, 0.5)),
    ('perspective', (0, 0, -1, 1), (10, -20, 400, 200), (210, 80, 0)),
    ('perspective', (0, 0, -9, 1), (10, -20, 400, 200), (210, 80, 1)),
    ('perspective', (1, 1, -2, 1), (10, -20, 400, 200), (260, 130, 9 / 16)),
    ('model-perspective', (1, 2, -2, 1), (-30, 20, 600, 300), (270, 120, 15 / 16)),
    ('model-orthographic', (1, 2, -2, 1), (-30, 20, 600, 300), (120, 57.5, 0.5)),
])
def test_independent_known_answers(kind, position, viewport, window,
                                  model_wrapped, proj_wrapped):
    model_rows = matrix.identity(4)
    if kind.startswith('model-'):
        # +Z quarter-turn followed by translation (2, -3, -4).
        model_rows = [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [2, -3, -4, 1]]
    if kind.endswith('orthographic'):
        projection_rows = matrix.orthographic(-2, 6, -3, 5, 2, 10).matrix
    elif kind.endswith('perspective'):
        projection_rows = matrix.perspective(90, 2, 1, 9).matrix
    else:
        projection_rows = matrix.identity(4)
    model = argument(model_rows, model_wrapped)
    projection = argument(projection_rows, proj_wrapped)
    obj = V(*position)
    view = list(viewport)
    before = (snapshot(model), snapshot(projection), obj.vector[:], view[:])
    projected = matrix.project(obj, model, projection, view)
    assert_vector3(projected, window)
    unprojected = matrix.unproject(*window, model, projection, view)
    assert_vector3(unprojected, position[:3])
    assert projected is not obj and unprojected is not obj
    assert before == (snapshot(model), snapshot(projection), obj.vector, view)
    projected.vector[0] = 999
    unprojected.vector[0] = 888
    assert before == (snapshot(model), snapshot(projection), obj.vector, view)


@pytest.mark.parametrize('function,expected', [
    (matrix.perspective, (260, 130, 9 / 16)),
    (matrix.perspectiveX, (310, 180, 9 / 16)),
])
def test_vertical_and_horizontal_perspective_known_answers(function, expected):
    projection = function(90, 2, 1, 9)
    identity = matrix.Matrix(4)
    viewport = [10, -20, 400, 200]
    assert_vector3(matrix.project(V(1, 1, -2, 1), identity, projection, viewport), expected)
    assert_vector3(matrix.unproject(*expected, identity, projection, viewport), [1, 1, -2])


@pytest.mark.parametrize('kind', ['orthographic', 'perspective'])
@pytest.mark.parametrize('seed', range(8))
def test_analytic_camera_properties_and_round_trips(kind, seed):
    rng = random.Random(2500 + seed)
    x, y, z = rng.uniform(-2, 2), rng.uniform(-2, 2), rng.uniform(-8, -3)
    tx, ty, tz = rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-2, 0)
    model_rows = [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [tx, ty, tz, 1]]
    ex, ey, ez = -y + tx, x + ty, z + tz
    viewport = [rng.randint(-50, 50), rng.randint(-50, 50), 640, 240]
    near, far = 2, 14
    if kind == 'orthographic':
        proj = matrix.orthographic(-4, 6, -3, 5, near, far)
        fractions = [(ex + 4) / 10, (ey + 3) / 8, (-ez - near) / (far - near)]
    else:
        proj = matrix.perspective(90, 2, near, far)
        # tan(90 degrees / 2) = 1; depth follows the analytic pinhole mapping.
        fractions = [(ex / (-2 * ez) + 1) / 2, (ey / (-ez) + 1) / 2,
                     far / (far - near) + far * near / ((far - near) * ez)]
    expected = [viewport[0] + viewport[2] * fractions[0],
                viewport[1] + viewport[3] * fractions[1], fractions[2]]
    model = argument(model_rows, seed % 2 == 0)
    projection = proj if seed % 3 == 0 else copy.deepcopy(proj.matrix)
    obj = V(x, y, z, 1)
    result = matrix.project(obj, model, projection, viewport)
    assert_vector3(result, expected)
    # Independently calculated windows prevent cancelling errors in the inverse.
    assert_vector3(matrix.unproject(*expected, model, projection, viewport), [x, y, z])
    assert_vector3(matrix.unproject(*result.vector, model, projection, viewport), [x, y, z])


@pytest.mark.parametrize('scale', [-3, 0.5, 2])
def test_explicit_homogeneous_scale(scale):
    obj = V(scale, scale, -2 * scale, scale)
    projection = matrix.perspective(90, 2, 1, 9)
    viewport = [10, -20, 400, 200]
    assert_vector3(matrix.project(obj, matrix.Matrix(4), projection, viewport), [260, 130, 9 / 16])
    assert obj.vector == [scale, scale, -2 * scale, scale]


def test_direction_input_uses_supplied_zero_w():
    # With perspective projection, a direction can have nonzero clip W.
    model = matrix.Matrix(4).translate(V(5, -6, -7))
    projection = matrix.perspective(90, 2, 1, 9)
    assert_vector3(matrix.project(V(1, 1, -2, 0), model, projection, [0, 0, 400, 200]),
                   [250, 150, 9 / 8])


def test_project_does_not_promote_vector3():
    model, projection = matrix.Matrix(4), matrix.Matrix(4)
    position = V(1, 2, 3)
    before = (snapshot(model), snapshot(projection), position.vector[:])
    with pytest.raises(IndexError):
        matrix.project(position, model, projection, [0, 0, 100, 100])
    assert before == (snapshot(model), snapshot(projection), position.vector)


@pytest.mark.parametrize('model_wrapped,proj_wrapped', [
    (False, False), (False, True), (True, False), (True, True),
])
def test_zero_w_and_singular_contracts(model_wrapped, proj_wrapped):
    model = argument(matrix.identity(4), model_wrapped)
    swap = argument([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]],
                    proj_wrapped)
    viewport = [0, 0, 100, 100]
    original = (snapshot(model), snapshot(swap), viewport[:])
    obj = V(0, 2, 3, 1)
    with pytest.raises(ZeroDivisionError):
        matrix.project(obj, model, swap, viewport)
    sentinel = matrix.unproject(50, 50, 0.5, model, swap, viewport)
    assert_vector3(sentinel, [0, 0, 0])
    again = matrix.unproject(50, 50, 0.5, model, swap, viewport)
    assert again is not sentinel and again.vector is not sentinel.vector
    assert_vector3(matrix.unproject(75, 50, 0.5, model, swap, viewport), [2, 0, 0])
    assert obj.vector == [0, 2, 3, 1]
    assert original == (snapshot(model), snapshot(swap), viewport)

    singular_rows = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]
    singular = argument(singular_rows, proj_wrapped)
    # Projection need not invert its input; this finite collapse is valid.
    assert_vector3(matrix.project(V(1, -1, 7, 1), model, singular, viewport), [100, 0, 0.5])
    with pytest.raises(ZeroDivisionError):
        matrix.unproject(50, 50, 0.5, model, singular, viewport)
    assert snapshot(singular) == snapshot(argument(singular_rows, proj_wrapped))


@pytest.mark.parametrize('seed', range(6))
def test_general_projective_matrix_against_independent_inverse(seed):
    rng = random.Random(2600 + seed)
    # Diagonally dominant integer matrices are invertible at ordinary scales.
    rows = [[rng.randint(-2, 2) + (12 if i == j else 0) for j in range(4)]
            for i in range(4)]
    inverse_rows = inverse(rows)  # Independent Fraction Gauss-Jordan oracle.
    position = [rng.randint(-3, 3) for _ in range(3)] + [1]
    clip = [sum(position[i] * rows[i][j] for i in range(4)) for j in range(4)]
    assert clip[3] != 0
    ndc = [value / clip[3] for value in clip[:3]]
    viewport = [-13, 27, 320, 180]
    window = [-13 + (ndc[0] + 1) * 160, 27 + (ndc[1] + 1) * 90, (ndc[2] + 1) / 2]
    projected = matrix.project(V(*position), matrix.identity(4), rows, viewport)
    assert_vector3(projected, window)
    output = [sum(([ndc[0], ndc[1], ndc[2], 1])[i] * inverse_rows[i][j]
                  for i in range(4)) for j in range(4)]
    expected = [value / output[3] for value in output[:3]]
    assert_vector3(matrix.unproject(*window, matrix.identity(4), rows, viewport), expected)
    assert expected == pytest.approx(position[:3], abs=1e-12)
