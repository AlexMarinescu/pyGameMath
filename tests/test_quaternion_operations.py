import math
import random

import pytest

from gem import quaternion as q, vector


def V(*values):
    return vector.Vector(len(values), list(values))


@pytest.mark.parametrize('function', [q.quat_from_axis_angle, q.quat_rotate_from_axis_angle])
@pytest.mark.parametrize('wrapped', [False, True])
@pytest.mark.parametrize('axis', [[2, 0, 0], [0, -3, 0], [0, 0, 4], [1, 2, 2], [-2, 1, -2]])
@pytest.mark.parametrize('angle', [0, 90, -90, 180])
def test_axis_known_answers_and_ownership(function, wrapped, axis, angle):
    source = axis[:]
    argument = vector.Vector(3, source) if wrapped else source
    old_storage = argument.vector if wrapped else argument
    length = math.sqrt(sum(value * value for value in axis))
    unit = [value / length for value in axis]
    if function is q.quat_from_axis_angle:
        half = math.radians(angle) / 2
        expected = [math.cos(half)] + [value * math.sin(half) for value in unit]
    else:
        # Preserve the legacy result: a unit axis rotated about itself is unchanged.
        expected = [0] + unit
    result = function(argument, angle)
    assert isinstance(result, q.Quaternion)
    assert isinstance(result.data, list) and len(result.data) == 4
    assert result.data == pytest.approx(expected, abs=1e-14)
    assert source == axis
    assert (argument.vector if wrapped else argument) is old_storage
    result.data[1] = 99
    assert source == axis


@pytest.mark.parametrize('seed', range(6))
def test_axis_representation_and_scale_properties(seed):
    rng = random.Random(3000 + seed)
    axis = [rng.uniform(-5, 5) for _ in range(3)]
    angle = rng.uniform(-180, 180)
    norm = math.sqrt(sum(value * value for value in axis))
    expected = [math.cos(math.radians(angle) / 2)] + [
        value / norm * math.sin(math.radians(angle) / 2) for value in axis]
    for scale in [0.5, 3]:
        scaled = [scale * value for value in axis]
        original = scaled[:]
        result = q.quat_from_axis_angle(scaled, angle)
        assert result.data == pytest.approx(expected, abs=1e-14)
        assert scaled == original
        wrapped = V(*scaled)
        assert q.quat_from_axis_angle(wrapped, angle).data == pytest.approx(expected, abs=1e-14)
        assert wrapped.vector == original
        assert q.quat_rotate_from_axis_angle(scaled, angle).data == pytest.approx(
            [0] + [value / norm for value in axis], abs=1e-14)
        assert scaled == original


@pytest.mark.parametrize('function', [q.quat_from_axis_angle, q.quat_rotate_from_axis_angle])
def test_axis_unsupported_and_zero_behavior(function):
    for axis in [(0, 0, 1), None, 1, 'z']:
        assert function(axis, 90) is NotImplemented
    for argument in [[0, 0, 0], V(0, 0, 0)]:
        old = argument.vector if isinstance(argument, vector.Vector) else argument
        with pytest.raises(ZeroDivisionError):
            function(argument, 90)
        assert (argument.vector if isinstance(argument, vector.Vector) else argument) is old
        assert old == [0, 0, 0]


@pytest.mark.parametrize('helper,axis,expected', [
    (q.quat_rotate_x_from_angle, [1, 0, 0], [math.sqrt(0.5), math.sqrt(0.5), 0, 0]),
    (q.quat_rotate_y_from_angle, [0, 1, 0], [math.sqrt(0.5), 0, math.sqrt(0.5), 0]),
    (q.quat_rotate_z_from_angle, [0, 0, 1], [math.sqrt(0.5), 0, 0, math.sqrt(0.5)]),
])
def test_axis_angle_units_remain_distinct(helper, axis, expected):
    assert helper(math.pi / 2) == pytest.approx(expected)
    assert q.quat_from_axis_angle(axis, 90).data == pytest.approx(expected)


@pytest.mark.parametrize('components,values,expected', [
    ([1, 2, 3, 4], [5, 6, 7], [-56, 2, 12, 4]),
    ([-2, 3, -1, 4], [-5, 2, 6], [-7, -4, -42, -11]),
    ([1, 0, 0, 0], [1, -2, 3], [0, 1, -2, 3]),
    ([0, 0, 0, 1], [1, 0, 0], [0, 0, 1, 0]),
    ([2, -3, 4, 5], [0, 0, 0], [0, 0, 0, 0]),
])
def test_vector_product_known_answers(components, values, expected):
    source = components[:]
    rot = q.Quaternion(source)
    v = V(*values)
    result = rot * v
    assert isinstance(result, q.Quaternion) and result is not rot
    assert result.data == expected
    assert rot.data is source and source == components
    receiver = rot
    rot *= v
    assert rot is receiver
    assert isinstance(rot, q.Quaternion) and rot.data == expected
    assert source == components and rot.data is not source
    assert v.vector == values


@pytest.mark.parametrize('seed', range(8))
def test_vector_product_hamilton_oracle(seed):
    rng = random.Random(3100 + seed)
    w, x, y, z = [rng.randint(-5, 5) for _ in range(4)]
    a, b, c = [rng.randint(-5, 5) for _ in range(3)]
    # q*(0,v) = (-u dot v, w*v + u cross v), without the rotation sandwich.
    expected = [-sum(u * v for u, v in zip([x, y, z], [a, b, c])),
                w * a + y * c - z * b, w * b + z * a - x * c, w * c + x * b - y * a]
    rot = q.Quaternion([w, x, y, z])
    result = rot * V(a, b, c)
    assert result.data == expected
    assert rot.__imul__(V(a, b, c)) is rot
    assert rot.data == expected


@pytest.mark.parametrize('components,divisor,expected', [
    ([2, 4, 6, 8], 2.0, [1, 2, 3, 4]),
    ([1, -2, 3, -4], -0.5, [-2, 4, -6, 8]),
    ([0, 0, 0, 0], 3.0, [0, 0, 0, 0]),
])
def test_division_known_answers_and_ownership(components, divisor, expected):
    source = components[:]
    rot = q.Quaternion(source)
    for method in [rot.__div__, rot.__truediv__]:
        result = method(divisor)
        assert isinstance(result, q.Quaternion) and result is not rot
        assert result.data == expected
        assert result.data is not source
    assert (rot / divisor).data == expected
    assert rot.data is source and source == components
    legacy = q.Quaternion(components[:])
    assert legacy.__idiv__(divisor) is legacy
    assert legacy.data == expected
    receiver = rot
    rot /= divisor
    assert rot is receiver and rot.data == expected
    assert source == components and rot.data is not source


@pytest.mark.parametrize('operand', [2, True, 2 + 0j, '2', [2], V(2, 0, 0), q.Quaternion(), None])
def test_division_retains_float_only_acceptance(operand):
    source = [1, 2, 3, 4]
    rot = q.Quaternion(source)
    for method in [rot.__div__, rot.__idiv__, rot.__truediv__, rot.__itruediv__]:
        assert method(operand) is NotImplemented
    with pytest.raises(TypeError):
        rot / operand
    receiver = rot
    with pytest.raises(TypeError):
        rot /= operand
    assert rot is receiver and rot.data is source and source == [1, 2, 3, 4]


@pytest.mark.parametrize('divisor', [0.0, -0.0])
def test_division_zero_preserves_receiver(divisor):
    source = [1, -2, 3, -4]
    rot = q.Quaternion(source)
    for method in [rot.__div__, rot.__idiv__, rot.__truediv__, rot.__itruediv__]:
        with pytest.raises(ZeroDivisionError):
            method(divisor)
        assert rot.data is source and source == [1, -2, 3, -4]
    with pytest.raises(ZeroDivisionError):
        rot / divisor
    receiver = rot
    with pytest.raises(ZeroDivisionError):
        rot /= divisor
    assert rot is receiver and rot.data is source


def test_division_accepts_existing_float_subclasses():
    class Float(float):
        pass
    rot = q.Quaternion([2, 4, 6, 8])
    assert (rot / Float(2)).data == [1, 2, 3, 4]
    assert rot.__itruediv__(Float(2)) is rot
    assert rot.data == [1, 2, 3, 4]
