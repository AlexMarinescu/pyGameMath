"""Independent range and ownership regressions for Phase 4G-2R."""
from fractions import Fraction
import math
import sys

import pytest

from gem import bezier
from gem.legendre import Legendre
from gem.vector import Vector
from .test_curves_legendre_audit import bernstein, polynomial_coefficients


def evaluate(degree, controls, t):
    function = bezier.quadraticBezierPoint if degree == 2 else bezier.cubicBezierPoint
    return function(t, *controls)


def rounded(value):
    try:
        return float(value)
    except OverflowError:
        return math.copysign(float('inf'), 1 if value > 0 else -1)


def ordinary_reference(degree, x):
    return rounded(sum(c*Fraction(x)**power for power, c in polynomial_coefficients(degree).items()))


def assert_rounded(actual, expected, ulps=4):
    # Absolute ulp budgets also cover subnormal answers and exact zero.
    assert actual == pytest.approx(expected, rel=0, abs=ulps*math.ulp(expected))


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('dimension', [None, 0, 1, 2, 3, 4, 8])
def test_tiny_parameters_representable_and_unrepresentable_terms(degree, dimension):
    for t in [math.ldexp(1., -1074), 1e-300, 1e-200, 1e-170, 1e-160,
              1e-155, 1e-154, 1e-110, 1e-105, 1e-103, 5e-103, 1e-80]:
        for signed_t in [t, -t]:
            values = [0.]*degree + [1e300]
            expected = float(bernstein(values, signed_t))
            if dimension is None:
                assert_rounded(evaluate(degree, values, signed_t), expected)
                continue
            rows = [[p*(-1)**j for j in range(dimension)] for p in values]
            controls = [Vector(dimension, row) for row in rows]
            snapshot = [row[:] for row in rows]
            result = evaluate(degree, controls, signed_t)
            assert type(result) is Vector and result.size == dimension
            for j, actual in enumerate(result.vector):
                assert_rounded(actual, expected*(-1)**j)
            assert rows == snapshot and all(p.vector is row for p, row in zip(controls, rows))
            assert all(result.vector is not row for row in rows)
            if dimension:
                result.vector[0] = 999
                assert rows == snapshot


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('dimension', [None, 2, 3])
def test_every_bernstein_term_with_large_and_small_controls(degree, dimension):
    for index in range(degree+1):
        for t, control in [(1e-200, 1e308), (-1e-110, -1e308),
                           (1e-114, 1e-210), (1e-110, 1e-100),
                           (math.ldexp(1., -1074), 1e308)]:
            values = [0.]*(degree+1); values[index] = control
            expected = float(bernstein(values, t))
            if dimension is None:
                assert_rounded(evaluate(degree, values, t), expected)
            else:
                controls = [Vector(dimension, [p*(-1)**j for j in range(dimension)]) for p in values]
                result = evaluate(degree, controls, t)
                for j, actual in enumerate(result.vector):
                    assert_rounded(actual, expected*(-1)**j)


@pytest.mark.parametrize('degree,bound', [(2, math.ldexp(1., -511)), (3, math.ldexp(1., -340))])
def test_parameter_range_boundary_and_cancellation(degree, bound):
    for t in [math.nextafter(bound, 0.), bound, math.nextafter(bound, math.inf)]:
        values = [0.]*degree + [1e300]
        assert_rounded(evaluate(degree, values, t), float(bernstein(values, t)))
    t = 1e-200 if degree == 2 else 1e-150
    tail = float(bernstein([0.]*degree+[1e300], t))
    # Cancellation may expose rounding; use the absolute weighted input scale
    # instead of a relative comparison with an almost-zero polynomial.
    for p0 in [-tail, math.nextafter(-tail, 0.), math.nextafter(-tail, -math.inf)]:
        values = [p0]+[0.]*(degree-1)+[1e300]
        expected = float(bernstein(values, t))
        scale = float(bernstein([abs(p) for p in values], t))
        assert evaluate(degree, values, t) == pytest.approx(expected, rel=0, abs=8*math.ulp(scale))


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('t', [-1.5, -.25, 0., .125, .5, 1., 1.75])
def test_ordinary_dyadic_extrapolation_and_endpoint_ownership(degree, t):
    controls = [Vector(3, [(-1)**i*(i+1), 2*i-3, i*i]) for i in range(degree+1)]
    original = [p.vector[:] for p in controls]
    result = evaluate(degree, controls, t)
    expected = [float(bernstein([row[j] for row in original], t)) for j in range(3)]
    assert result.vector == expected
    assert all(result is not p and result.vector is not p.vector for p in controls)
    assert [p.vector for p in controls] == original


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('t', [0., 1., 1e-200, -1e-200])
def test_custom_vector_multiplication_dispatch_is_preserved(degree, t):
    events = []
    class Custom(Vector):
        def __mul__(self, weight):
            events.append(weight)
            return super().__mul__(weight)
    controls = [Custom(2, [float(i), -float(i)]) for i in range(degree+1)]
    evaluate(degree, controls, t)
    u = 1-t
    expected = [u*u, 2*u*t, t*t] if degree == 2 else [u*u*u, 3*u*u*t, 3*u*t*t, t*t*t]
    assert events == expected


@pytest.mark.parametrize('degree', [2, 3])
def test_fraction_dispatch_and_native_mismatched_dimensions(degree):
    t = Fraction(1, 10**200)
    values = [Fraction(0)]*degree+[Fraction(10**300)]
    assert evaluate(degree, values, t) == bernstein(values, t)
    controls = [Vector(2, [0.,0.])] + [Vector(3, [0.,0.,7.]) for _ in range(degree)]
    result = evaluate(degree, controls, 1e-200)
    assert result.size == 2 and result.vector == [0.,0.]


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('nonfinite', [math.inf, -math.inf, math.nan])
def test_nonfinite_controls_keep_historical_weight_arithmetic(degree, nonfinite):
    # A zero materialized tail weight times infinity historically produces NaN.
    values = [0.]*degree+[nonfinite]
    assert math.isnan(evaluate(degree, values, 1e-200))
    vectors = [Vector(2, [p,p]) for p in values]
    assert all(math.isnan(x) for x in evaluate(degree, vectors, 1e-200).vector)


@pytest.mark.parametrize('degree', [2, 3])
def test_large_integer_and_malformed_storage_exceptions_remain(degree):
    with pytest.raises(OverflowError):
        evaluate(degree, [0.]*degree+[10**1000], 1e-200)
    controls = [Vector(2, [0.,0.]) for _ in range(degree+1)]
    controls[-1].vector.pop()
    with pytest.raises(IndexError):
        evaluate(degree, controls, 1e-200)


@pytest.mark.parametrize('degree,x', [(2, 1e154), (3, 4e102), (4, 8e76),
                                    (8, 1e38), (12, 1.3*math.ldexp(1.,84))])
@pytest.mark.parametrize('sign', [-1, 1])
def test_ordinary_scaled_recurrence_exact_polynomial(degree, x, sign):
    x *= sign
    expected = ordinary_reference(degree, x)
    assert math.isfinite(expected)
    polynomial = Legendre(degree, 0, x)
    polynomial.P, polynomial.PM1, polynomial.PML = 123., -456., 789.
    snapshot = polynomial.__dict__.copy()
    actual = polynomial.run()
    assert actual == pytest.approx(expected, rel=4e-15, abs=0)
    assert polynomial.__dict__ == snapshot and polynomial.run() == actual
    assert polynomial.calculatePML(degree) is None
    assert polynomial.P == 1. and polynomial.PM1 == x and polynomial.PML == actual


@pytest.mark.parametrize('degree', [2, 3, 4])
@pytest.mark.parametrize('sign', [-1, 1])
def test_legendre_at_finite_output_boundary(degree, sign):
    leading = polynomial_coefficients(degree)[degree]
    # Establish the exact polynomial result on either side of binary64's range;
    # the boundary estimate is only for choosing inputs, not the oracle.
    x = math.exp((math.log(sys.float_info.max)-math.log(float(leading)))/degree)
    for factor in [.9999999999999, 1., 1.0000000000001]:
        argument = sign*x*factor
        expected = ordinary_reference(degree, argument)
        actual = Legendre(degree, 0, argument).run()
        if math.isfinite(expected):
            assert actual == pytest.approx(expected, rel=4e-15)
        else:
            assert actual == expected


@pytest.mark.parametrize('degree', [63, 64, 65])
@pytest.mark.parametrize('x', [-2., -1.25, 1.25, 2., math.nextafter(2., math.inf)])
def test_conservative_recurrence_range_bound(degree, x):
    assert Legendre(degree,0,x).run() == pytest.approx(ordinary_reference(degree,x), rel=3e-14)


@pytest.mark.parametrize('x', [-10**154, 10**154])
def test_ordinary_integer_argument_with_finite_answer(x):
    polynomial = Legendre(2,0,x)
    assert polynomial.run() == pytest.approx(ordinary_reference(2,x), rel=4e-15)
    assert polynomial.x == x and type(polynomial.x) is int


@pytest.mark.parametrize('degree,order,x', [(2,1,.5), (12,7,-.375), (32,12,.875)])
def test_associated_phase_and_scratch_contract(degree, order, x):
    from .test_curves_legendre_audit import legendre_reference
    polynomial = Legendre(degree, order, x)
    original = polynomial.__dict__.copy()
    assert polynomial.run() == pytest.approx(legendre_reference(degree, order, x), rel=3e-13)
    assert polynomial.__dict__ == original
    assert Legendre(degree, order, -x).run() == (-1)**(degree+order)*polynomial.run()


def test_legendre_historical_invalid_and_nonfinite_behavior():
    with pytest.raises(ValueError):
        Legendre(2,1,1.01).run()
    with pytest.raises(TypeError):
        Legendre(2,0.,2.).run()
    with pytest.raises(TypeError):
        Legendre(2,1.5,.2).run()
    assert Legendre(0,0,math.inf).run() == 1.
    assert Legendre(2,0,math.inf).run() == math.inf
    assert math.isnan(Legendre(3,0,math.inf).run())
    assert math.isnan(Legendre(2,0,math.nan).run())
