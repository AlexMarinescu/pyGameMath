"""Exact equality and clamp ownership contracts independent of implementation."""
import pytest
from gem.vector import Vector, clamp


def V(values):
    return Vector(len(values), list(values))


@pytest.mark.parametrize('size', [0, 2, 3, 4])
def test_exact_equality_and_inequality(size):
    values = [float(i + 1) for i in range(size)]
    a, b = V(values), V(values)
    assert (a == b) is True and (b == a) is True
    assert (a != b) is False and (b != a) is False
    for i in range(size):
        changed = values[:]; changed[i] += 1e-12
        c = V(changed)
        assert (a == c) is False and (c == a) is False
        assert (a != c) is True and (c != a) is True
    assert a.vector == values and b.vector == values


@pytest.mark.parametrize('first', [0, 2, 3, 4])
@pytest.mark.parametrize('second', [0, 2, 3, 4])
def test_dimensions_both_operand_orders(first, second):
    a, b = Vector(first), Vector(second)
    assert (a == b) is (first == second)
    assert (b == a) is (first == second)
    assert (a != b) is (first != second)
    assert (b != a) is (first != second)


@pytest.mark.parametrize('other', [None, [], [1, 2], 1, 'vector'])
def test_unsupported_operand_protocol(other):
    a = V([1, 2])
    assert a.__eq__(other) is NotImplemented
    assert a.__ne__(other) is NotImplemented
    assert (a == other) is False and (a != other) is True


@pytest.mark.parametrize('size', [0, 2, 3, 4])
def test_clamp_fresh_storage_and_receiver_mutation(size):
    values = [-2, 2, 10, -8][:size]
    lower, upper = [0]*size, [5]*size
    receiver = Vector(size, values)
    original = receiver.vector[:]
    old_storage = receiver.vector
    for result in [clamp(size, values, lower, upper), receiver.clamp(size, values, lower, upper)]:
        assert result.vector == [0, 2, 5, 0][:size]
        assert result.vector is not values and result.vector is not old_storage
        assert values == original and receiver.vector == original
        if size:result.vector[0] = 99
        assert values == original
    assert receiver.i_clamp(size, values, lower, upper) is receiver
    assert receiver.vector == [0, 2, 5, 0][:size]
    assert receiver.vector is not old_storage
    assert old_storage == original and values == original


def test_clamp_shared_value_and_bound_lists():
    values = [-2, 2, 10]
    a, b = Vector(3, values), Vector(3, values)
    a.i_clamp(3, values, [0]*3, [5]*3)
    assert a.vector == [0, 2, 5]
    assert b.vector == values == [-2, 2, 10]
    # Aliased bounds remain inputs, not output scratch space.
    result = clamp(3, values, values, [5]*3)
    assert result.vector == [-2, 2, 10] and values == [-2, 2, 10]
    assert clamp(3, values, [0]*3, values).vector == [0, 2, 10]
    assert values == [-2, 2, 10]


@pytest.mark.parametrize('components,width,height,expected', [
    ([3, 4], 100, 200, [83, 184, 100, 200]),
    ([-3, -4], 100, 200, [17, 16, 100, 200]),
    ([0, 2], 80, 40, [40, 42, 80, 40]),
    ([3, 4, 12], 130, 260, [83, 174, 130, 260]),
    ([1, 2, 2, 4], 100, 50, [61, 37, 100, 50]),
    ([3, 4], 0, -10, [3, -5, 0, -10]),
])
def test_viewport_literal_known_answers(components, width, height, expected):
    from gem.common import getViewPort
    coords = V(components)
    storage = coords.vector
    first = getViewPort(coords, width, height)
    second = getViewPort(coords, width, height)
    assert first == pytest.approx(expected)
    assert second == pytest.approx(expected)
    assert first is not second
    first[0] = 999
    assert coords.vector is storage and coords.vector == components


@pytest.mark.parametrize('size', [2, 3, 4])
def test_zero_viewport_preserves_historical_error(size):
    from gem.common import getViewPort
    coords = Vector(size)
    with pytest.raises(ZeroDivisionError):
        getViewPort(coords, 100, 50)
    assert coords.vector == [0]*size
