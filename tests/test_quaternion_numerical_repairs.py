"""Independent range references and ordinary-path compatibility for 4G-1R."""
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import math
import random
import struct

import pytest

from gem.quaternion import Quaternion, quat_mul_quat
from .test_core_algebra_audit import hamilton


MIN_NORMAL = math.ldexp(1.0, -1022)
TINY = math.ldexp(1.0, -1074)


def decimal_axis(values):
    with localcontext() as context:
        context.prec = 800
        entries = [Decimal.from_float(x) for x in values]
        length = sum(x*x for x in entries).sqrt()
        return [float(x/length) for x in entries]


@pytest.mark.parametrize('scale', [TINY, math.ldexp(1.0, -1060), 1e-320,
                                   MIN_NORMAL/4, MIN_NORMAL, 1e-200])
@pytest.mark.parametrize('components', [(1, 1, 0), (1, 2, -3), (0, 1, -1)])
@pytest.mark.parametrize('exponent', [-0.5, 0.5, 1.5])
def test_scaled_axis_fractional_powers_and_log(scale, components, exponent):
    data = [-1.0] + [component*scale for component in components]
    axis = decimal_axis(data[1:])
    # The principal angle differs from pi by at most 4e-200 in these fixtures.
    expected = [math.cos(math.pi*exponent)] + [x*math.sin(math.pi*exponent) for x in axis]
    q = Quaternion(data)
    original = data[:]
    power, logarithm = q.pow(exponent), q.log()
    assert power.data == pytest.approx(expected, rel=3e-15, abs=3e-15)
    assert math.hypot(*power.data) == pytest.approx(1.0, abs=3e-15)
    assert logarithm == pytest.approx([0.0] + [x*math.pi for x in axis], rel=3e-15, abs=0)
    assert power is not q and power.data is not data and logarithm is not data
    assert q.data is data and data == original


@pytest.mark.parametrize('scale', [TINY, 4*TINY, 1e-320, MIN_NORMAL/4,
                                   math.nextafter(MIN_NORMAL, 0), MIN_NORMAL,
                                   math.nextafter(MIN_NORMAL, math.inf)])
def test_positive_identity_adjacent_representable_results(scale):
    data = [1.0, scale, -scale, 0.0]
    q = Quaternion(data)
    # log([1,v]) = [0,v] and its square root = [1,v/2] to the
    # representable accuracy here. A subnormal result may round by one ulp.
    assert q.log() == pytest.approx([0.0] + data[1:], rel=3e-15, abs=TINY)
    assert q.pow(0.5).data == pytest.approx([1.0, scale/2, -scale/2, 0.0], rel=3e-15, abs=TINY)
    assert q.pow(0).data == [1, 0, 0, 0]
    assert q.pow(1).data == data and q.pow(1).data is not data


CYCLIC = [[0.0] + [float(sign) if i == axis else 0.0 for i in range(3)]
          for axis in range(3) for sign in [-1, 1]]
CYCLIC += [list(values) for values in itertools.product([-0.5, 0.5], repeat=4)]


def test_exact_control_set_is_closed_under_hamilton_products():
    controls = CYCLIC + [[1.0, 0.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0]]
    members = set(tuple(data) for data in controls)
    assert len(members) == 24
    for a, b in itertools.product(controls, repeat=2):
        exact = tuple(hamilton(a, b))
        assert exact in members
        assert quat_mul_quat(a, b) == [float(x) for x in exact]


@pytest.mark.parametrize('data', CYCLIC)
@pytest.mark.parametrize('exponent', [-10**16-1, 10**16+1, 2**200+1,
                                     -2**200+1, float(2**53-1), 1e308])
def test_exact_cyclic_controls_against_fraction_period(data, exponent):
    # The six unit basis controls have period 4; the sixteen half-component
    # controls have period 6 or 3. Exact arithmetic proves the period first.
    identity = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    product = identity
    period = None
    for count in range(1, 7):
        product = hamilton(product, data)
        if product == identity:
            period = count
            break
    assert period is not None
    assert 12 % period == 0
    expected = identity
    for _ in range(int(exponent) % period):
        expected = hamilton(expected, data)
    q = Quaternion(data)
    original = data[:]
    result = q.pow(exponent)
    assert result.data == [float(x) for x in expected]
    assert sum(x*x for x in result.data) == 1.0
    assert result is not q and result.data is not data
    assert q.data is data and data == original


def legacy_principal_power(data, exponent):
    """Pre-repair formula, only for compatibility checks of unchanged paths."""
    w, x, y, z = data
    length = math.hypot(math.hypot(x, y), z)
    angle = math.atan2(length, w)*exponent
    sine = math.sin(angle)
    return [math.cos(angle), x/length*sine, y/length*sine, z/length*sine]


def bits(values):
    return [struct.pack('!d', x) for x in values]


@pytest.mark.parametrize('seed', range(16))
def test_ordinary_general_powers_and_log_preserve_previous_bits(seed):
    rng = random.Random(471000+seed)
    axis = [rng.uniform(-2, 2) for _ in range(3)]
    length = math.sqrt(sum(x*x for x in axis))
    angle = rng.uniform(0.1, 2.9)
    data = [math.cos(angle)] + [x/length*math.sin(angle) for x in axis]
    q = Quaternion(data)
    imaginary = math.hypot(math.hypot(*data[1:3]), data[3])
    principal = math.atan2(imaginary, data[0])
    for exponent in [-16, -1, -0.5, 0.5, 2, 16]:
        assert bits(q.pow(exponent).data) == bits(legacy_principal_power(data, exponent))
    assert bits(q.log()) == bits([0.0] + [x/imaginary*principal for x in data[1:]])


@pytest.mark.parametrize('data', [[2.0, 1.0, 0.0, 0.0], [-2.0, 3.0, -4.0, 5.0],
                                 [0.5, 0.5, 0.5, 0.25], [0.0, 2.0, 0.0, 0.0]])
def test_nonunit_legacy_paths_are_not_replaced_with_general_powers(data):
    q = Quaternion(data)
    for exponent in [-3, -0.5, 0.5, 2, 10**16]:
        assert bits(q.pow(exponent).data) == bits(legacy_principal_power(data, exponent))
    assert q.pow(1).data == data and q.pow(1).data is not data
    assert q.pow(0).data == [1, 0, 0, 0]
    # These remain unsupported inputs; no new normalization or rejection policy.
    assert q.data is data


@pytest.mark.parametrize('index', range(4))
@pytest.mark.parametrize('toward', [0.0, math.inf])
def test_nextafter_neighbors_do_not_enter_exact_cyclic_path(index, toward):
    data = [0.5]*4
    data[index] = math.nextafter(0.5, toward)
    q = Quaternion(data)
    for exponent in [-2, 2, 10**16]:
        assert bits(q.pow(exponent).data) == bits(legacy_principal_power(data, exponent))


@pytest.mark.parametrize('exponent', [-3.0, -0.5, 0.0, 0.5, 1.0, 2.0])
def test_cyclic_fractional_principal_branch_and_ownership(exponent):
    data = [0.5, -0.5, 0.5, -0.5]
    q = Quaternion(data)
    if exponent % 1:
        assert bits(q.pow(exponent).data) == bits(legacy_principal_power(data, exponent))
    else:
        expected = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
        for _ in range(int(exponent) % 6):
            expected = hamilton(expected, data)
        assert q.pow(exponent).data == [float(x) for x in expected]
    assert q.pow(exponent).data is not data


def test_exact_zero_negative_identity_and_nonfinite_exponents_keep_behavior():
    zero = Quaternion([0.0]*4)
    for exponent in [0, 1, -1, 0.5]:
        with pytest.raises(ValueError):
            zero.pow(exponent)
    with pytest.raises(ValueError):
        Quaternion([-1, 0, 0, 0]).log()
    with pytest.raises(ValueError):
        Quaternion([-1, 0, 0, 0]).pow(0.5)
    assert Quaternion([-1, 0, 0, 0]).pow(-10**16-1).data == [-1, 0, 0, 0]
    for data in [[0, 1, 0, 0], [0.6, 0.8, 0, 0]]:
        with pytest.raises(ValueError):
            Quaternion(data).pow(math.inf)
        assert all(math.isnan(x) for x in Quaternion(data).pow(math.nan).data)


@pytest.mark.parametrize('data', [[0, 1, 0, 0], [0.6, 0.8, 0, 0]])
def test_existing_exponent_protocol_and_conversion_errors(data):
    q = Quaternion(data)
    for exponent in [Decimal(2), '2', 2j, None]:
        with pytest.raises(TypeError):
            q.pow(exponent)
    with pytest.raises(OverflowError):
        q.pow(10**400)
    assert q.pow(Fraction(1, 2)).data == q.pow(0.5).data
