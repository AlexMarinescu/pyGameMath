"""Independent unit-quaternion power and principal-logarithm regressions."""
import math
import pytest
from gem.quaternion import Quaternion, quat_pow, quat_log


def hamilton(a, b):
    w, x, y, z = a
    s, u, v, t = b
    return [w*s-x*u-y*v-z*t, w*u+x*s+y*t-z*v,
            w*v-x*t+y*s+z*u, w*t+x*v-y*u+z*s]


AXES = [(1, 0, 0), (0, 1, 0), (0, 0, -1), (1/3, -2/3, 2/3)]


@pytest.mark.parametrize('axis', AXES)
@pytest.mark.parametrize('angle', [math.pi/4, math.pi/2, 2*math.pi/3])
@pytest.mark.parametrize('exponent', [-3, -1, 0, 1, 2, 4])
def test_integer_powers_against_hamilton_products(axis, angle, exponent):
    data = [math.cos(angle)] + [c*math.sin(angle) for c in axis]
    operand = data if exponent >= 0 else [data[0]] + [-c for c in data[1:]]
    expected = [1, 0, 0, 0]
    for _ in range(abs(exponent)):
        expected = hamilton(expected, operand)
    assert Quaternion(data=data).pow(exponent).data == pytest.approx(expected, abs=2e-14)


@pytest.mark.parametrize('axis', AXES)
@pytest.mark.parametrize('exponent', [-0.5, 0.5, 1.5])
def test_fractional_known_angles(axis, exponent):
    # Quaternion angle pi/2 is a spatial half-turn.
    data = [0] + list(axis)
    expected = [math.cos(exponent*math.pi/2)] + [c*math.sin(exponent*math.pi/2) for c in axis]
    assert Quaternion(data=data).pow(exponent).data == pytest.approx(expected)


@pytest.mark.parametrize('exponent', [-4, -1, -0.5, 0, 0.5, 1, 8])
def test_identity_powers(exponent):
    assert Quaternion().pow(exponent).data == [1, 0, 0, 0]


@pytest.mark.parametrize('exponent', [-4, -3, 0, 1, 2, 3.0])
def test_negative_identity_integer_parity(exponent):
    assert Quaternion(data=[-1, 0, 0, 0]).pow(exponent).data == [(-1)**exponent, 0, 0, 0]


@pytest.mark.parametrize('exponent', [-0.5, 0.5, 1.5])
def test_negative_identity_fractional_rejected(exponent):
    with pytest.raises(ValueError):
        Quaternion(data=[-1, 0, 0, 0]).pow(exponent)


@pytest.mark.parametrize('axis', AXES)
@pytest.mark.parametrize('angle', [math.pi/4, math.pi/2, 3*math.pi/4])
def test_log_known_angles(axis, angle):
    data = [math.cos(angle)] + [c*math.sin(angle) for c in axis]
    assert Quaternion(data=data).log() == pytest.approx([0] + [c*angle for c in axis])


def test_identity_and_negative_identity_log():
    assert Quaternion().log() == [0, 0, 0, 0]
    with pytest.raises(ValueError):
        Quaternion(data=[-1, 0, 0, 0]).log()


@pytest.mark.parametrize('exponent', [-1, 0, 0.5, 1])
def test_zero_power_rejected(exponent):
    with pytest.raises(ValueError):
        Quaternion(data=[0, 0, 0, 0]).pow(exponent)


def test_zero_log_rejected():
    with pytest.raises(ValueError):
        Quaternion(data=[0, 0, 0, 0]).log()


@pytest.mark.parametrize('scale', [1e-20, 1e-200, 1e-320])
@pytest.mark.parametrize('scalar', [1, -1])
def test_tiny_imaginary_direction(scale, scalar):
    data = [scalar, scale, -2*scale, 2*scale]
    q = Quaternion(data=data)
    angle = 3*scale if scalar > 0 else math.pi
    assert q.log() == pytest.approx([0, angle/3, -2*angle/3, 2*angle/3], rel=2e-3, abs=0)
    expected = ([1, scale/2, -scale, scale] if scalar > 0 else [0, 1/3, -2/3, 2/3])
    assert q.pow(0.5).data == pytest.approx(expected, rel=2e-3, abs=1e-16 if scalar < 0 else 0)


def test_sign_is_not_canonicalized():
    c = math.sqrt(0.5)
    q = Quaternion(data=[-c, 0, 0, -c])
    assert q.log() == pytest.approx([0, 0, 0, -3*math.pi/4])
    assert q.pow(0.5).data == pytest.approx([math.cos(3*math.pi/8), 0, 0, -math.sin(3*math.pi/8)])


@pytest.mark.parametrize('data', [[1, 0, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0]])
def test_fresh_power_storage(data):
    original = data[:]
    q = Quaternion(data=data)
    for exponent in [0, 1, 2]:
        result = quat_pow(q, exponent)
        assert isinstance(result, Quaternion)
        assert result is not q and result.data is not data
        if exponent == 1:
            assert result.data == original
        result.data[0] = 12
    assert q.data is data and data == original


@pytest.mark.parametrize('data', [[1, 0, 0, 0], [0, 1, 0, 0]])
def test_fresh_log_storage(data):
    original = data[:]
    q = Quaternion(data=data)
    result = quat_log(q)
    assert type(result) is list and len(result) == 4 and result is not data
    assert result == q.log()
    result[0] = 12
    assert q.data is data and data == original


def test_large_finite_exponent_has_finite_unit_result():
    # The principal angle times this finite exponent exceeds float range.
    data = [-math.sqrt(0.5), math.sqrt(0.5), 0, 0]
    result = Quaternion(data=data).pow(1e308).data
    assert all(math.isfinite(c) for c in result)
    assert sum(c*c for c in result) == pytest.approx(1)
