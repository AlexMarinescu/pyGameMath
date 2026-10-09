"""Independent Phase 4G-2 references; numerical findings have been repaired.

Fraction Bernstein sums and differentiated Rodrigues coefficients are oracles.
Private subdivision helpers are tested as implementation details, not new APIs.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
import math
import random

import pytest

from gem import bezier
from gem.legendre import Legendre
from gem.vector import Vector
from gem.quaternion import Quaternion, squad4


def bernstein(controls, t):
    t = Fraction(t)
    degree = len(controls) - 1
    return sum(Fraction(p) * math.comb(degree, i) * t**i
               * (1-t)**(degree-i) for i, p in enumerate(controls))


def polynomial_coefficients(degree, order=0):
    """Differentiate the explicit Rodrigues polynomial, without recurrence."""
    result = {}
    for k in range((degree-order)//2+1):
        power = degree-2*k-order
        result[power] = Fraction(
            (-1)**k * math.factorial(2*degree-2*k),
            2**degree * math.factorial(k) * math.factorial(degree-k)
            * math.factorial(power))
    return result


def legendre_reference(degree, order, x):
    x = Fraction(x)  # Use the exact represented argument, not its decimal label.
    derivative = sum(c*x**power for power, c in polynomial_coefficients(degree, order).items())
    with localcontext() as context:
        context.prec = 180
        value = Decimal(derivative.numerator)/Decimal(derivative.denominator)
        base = 1-x*x
        factor = (Decimal(base.numerator)/Decimal(base.denominator)).sqrt()**order if order else Decimal(1)
        return float((-1)**order * value * factor)


def segment_distance(point, start, end):
    """Exact rational projection with a Decimal square root, not gem geometry."""
    delta = [Fraction(b)-Fraction(a) for a, b in zip(start, end)]
    offset = [Fraction(p)-Fraction(a) for p, a in zip(point, start)]
    denominator = sum(d*d for d in delta)
    t = max(Fraction(0), min(Fraction(1), sum(p*d for p, d in zip(offset, delta))/denominator)) if denominator else 0
    squared = sum((p-t*d)**2 for p, d in zip(offset, delta))
    with localcontext() as context:
        context.prec = 100
        return float((Decimal(squared.numerator)/Decimal(squared.denominator)).sqrt())


def evaluate(degree, controls, t):
    function = bezier.quadraticBezierPoint if degree == 2 else bezier.cubicBezierPoint
    return function(t, *controls)


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('dimension', [1, 2, 3, 4])
@pytest.mark.parametrize('seed', range(4))
def test_bernstein_affine_reversal_and_storage(degree, dimension, seed):
    rng = random.Random(472000+seed)
    coords = [[rng.randint(-32, 32)/8 for _ in range(dimension)] for _ in range(degree+1)]
    controls = [Vector(dimension, row) for row in coords]
    saved = [row[:] for row in coords]
    # Non-diagonal affine mapping, calculated without gem matrix operations.
    rows = [[Fraction((i+2*j)%5-2, 4) for j in range(dimension)] for i in range(dimension)]
    offsets = [Fraction(i-2, 2) for i in range(dimension)]
    mapped = [Vector(dimension, [float(sum(rows[j][k]*Fraction(row[k]) for k in range(dimension)) + offsets[j])
                                for j in range(dimension)]) for row in coords]
    for t in [-.5, 0, .125, .5, .875, 1, 1.5]:
        exact = [bernstein([p[j] for p in coords], t) for j in range(dimension)]
        result = evaluate(degree, controls, t)
        assert result.vector == [float(x) for x in exact]
        assert evaluate(degree, controls[::-1], 1-t).vector == result.vector
        expected_mapped = [float(sum(rows[j][k]*exact[k] for k in range(dimension))+offsets[j])
                           for j in range(dimension)]
        assert evaluate(degree, mapped, t).vector == expected_mapped
        assert all(result is not p and result.vector is not p.vector for p in controls)
        result.vector[0] = 12345
        assert coords == saved and all(p.vector is row for p, row in zip(controls, coords))


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('scale', [1e-250, 1., 1e250])
@pytest.mark.parametrize('coincident', [False, True])
def test_evaluation_representable_uniform_scales(degree, scale, coincident):
    values = [scale*(1 if coincident else (-1)**i*(i+1)) for i in range(degree+1)]
    for t in [0., .125, .5, .875, 1.]:
        expected = float(bernstein(values, t))
        # Cancellation can have an exact zero answer despite nonzero rounded
        # products. Bound absolute roundoff by the weighted input scale.
        weighted_scale = float(bernstein([abs(p) for p in values], t))
        assert evaluate(degree, values, t) == pytest.approx(
            expected, rel=3e-15, abs=16*2**-53*weighted_scale)


@pytest.mark.parametrize('degree', [2, 3])
def test_exact_scalar_fraction_arithmetic_fallback(degree):
    controls = [Fraction((-1)**i*(i+1), 7) for i in range(degree+1)]
    for t in [Fraction(-1, 3), Fraction(0), Fraction(2, 7), Fraction(1), Fraction(4, 3)]:
        assert evaluate(degree, controls, t) == bernstein(controls, t)


@pytest.mark.parametrize('degree,t', [(2, 1e-200), (3, 1e-150)])
@pytest.mark.parametrize('dimension', [None, 2, 3])
def test_tiny_parameter_large_control_representable_result(degree, t, dimension):
    values = [0.]*degree + [1e300]
    expected = float(bernstein(values, t))
    assert math.isfinite(expected) and expected > 1e-200
    if dimension is None:
        actual = evaluate(degree, values, t)
    else:
        data = [[p*(-1)**j for j in range(dimension)] for p in values]
        controls = [Vector(dimension, row) for row in data]
        saved = [row[:] for row in data]
        result = evaluate(degree, controls, t)
        assert all(result.vector is not row for row in data)
        assert data == saved and all(p.vector is row for p, row in zip(controls, data))
        actual = result.vector[0]
        # Each channel has an independent sign and the same reference scale.
        assert result.vector == pytest.approx([expected*(-1)**j for j in range(dimension)], rel=3e-14, abs=0)
    assert actual == pytest.approx(expected, rel=3e-14, abs=0)


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('t', [0., .25, .5, .75, 1.])
def test_derivative_identity_of_existing_evaluators(degree, t):
    # Five-point differentiation is exact for degree <= 3 in real arithmetic.
    # No public derivative helper exists; this checks the evaluated polynomial.
    controls = [-2., 3., -1., 4.][:degree+1]
    h = 1/32
    derivative = (evaluate(degree, controls, t-2*h) - 8*evaluate(degree, controls, t-h)
                  + 8*evaluate(degree, controls, t+h) - evaluate(degree, controls, t+2*h))/(12*h)
    expected = float(bernstein([degree*(b-a) for a, b in zip(controls, controls[1:])], t))
    assert derivative == pytest.approx(expected, rel=0, abs=2e-13)


@pytest.mark.parametrize('dimension', [1, 2, 3])
@pytest.mark.parametrize('t', [0., .125, .5, .875, 1.])
def test_private_split_independent_restricted_controls(dimension, t):
    coords = [tuple((i*i-2*i+j)*(-1)**j for j in range(dimension)) for i in range(4)]
    saved = coords[:]
    left, right = bezier._split(coords, t)
    expected_left = [tuple(float(bernstein([p[j] for p in coords[:i+1]], t))
                           for j in range(dimension)) for i in range(4)]
    expected_right = [tuple(float(bernstein([p[j] for p in coords[i:]], t))
                            for j in range(dimension)) for i in range(4)]
    assert left == expected_left and right == expected_right and coords == saved
    for s in [0., .25, .5, 1.]:
        for j in range(dimension):
            assert bernstein([p[j] for p in left], s) == bernstein([p[j] for p in coords], t*s)
            assert bernstein([p[j] for p in right], s) == bernstein([p[j] for p in coords], t+(1-t)*s)


@pytest.mark.parametrize('scale', [1e-250, 1., 1e250])
@pytest.mark.parametrize('polygon', [
    [(0.,), (-2.,), (6.,), (4.,)],
    [(0., 0.), (3., 4.), (-3., -4.), (0., 0.)],
    [(0., 0.), (-2., 1.), (6., -1.), (4., 0.)],
    [(0., 0., 0.), (1., 3., 4.), (3., -3., -4.), (4., 0., 0.)],
])
def test_flatness_and_chord_distance_rational_projection(scale, polygon):
    controls = [tuple(x*scale for x in p) for p in polygon]
    expected = [segment_distance(p, controls[0], controls[-1]) for p in controls[1:-1]]
    assert bezier._flatness(controls) == pytest.approx(max(expected), rel=4e-15, abs=0)
    for p, distance in zip(controls[1:-1], expected):
        assert bezier._chord_distance(p, controls[0], controls[-1]) == pytest.approx(distance, rel=4e-15, abs=0)


@pytest.mark.parametrize('dimension', [2, 3])
def test_adaptive_samples_dense_bernstein_error_and_parameter_order(dimension):
    coords = [[i/4, [0., 2., -2., 0.][i]] + ([-i/2] if dimension == 3 else []) for i in range(4)]
    source = [Vector(dimension, row) for row in coords]
    path = bezier.BezierPath(); path.setControlPoints(source)
    counts, errors = [], []
    for tolerance in [.04, .01]:
        path.minimum_sqr_distance = tolerance*tolerance
        samples = path.findDrawingPoints(0)
        counts.append(len(samples))
        assert samples[0].vector == coords[0] and samples[-1].vector == coords[-1]
        assert all(p.vector[0] < q.vector[0] for p, q in zip(samples, samples[1:]))
        for sample in samples:
            t = Fraction(sample.vector[0])/Fraction(3, 4)
            assert sample.vector == pytest.approx([float(bernstein([p[j] for p in coords], t))
                                                  for j in range(dimension)], rel=0, abs=2e-14)
            assert all(sample.vector is not row for row in coords)
        # Since x is monotonic, use the unique bracketing chord, independently
        # of the sampler's stack and flatness calculation.
        error, index = 0., 0
        for i in range(257):
            t = Fraction(i, 256)
            point = [float(bernstein([p[j] for p in coords], t)) for j in range(dimension)]
            while index+2 < len(samples) and samples[index+1].vector[0] < point[0]:
                index += 1
            error = max(error, segment_distance(point, samples[index].vector, samples[index+1].vector))
        assert error <= tolerance
        errors.append(error)
    assert counts[1] >= counts[0] and errors[1] <= errors[0]
    assert all(p.vector is row for p, row in zip(source, coords))


@pytest.mark.parametrize('dimension', [2, 3])
@pytest.mark.parametrize('interval', [(0., 0.), (.25, .25), (1., 1.), (0., 1.), (.125, .875)])
def test_subinterval_samples_independent_polynomial_and_ownership(dimension, interval):
    # x=3t/4, y=6t(1-t)(1-2t), z=-3t/2.
    coords = [[i/4, [0., 2., -2., 0.][i]] + ([-i/2] if dimension == 3 else []) for i in range(4)]
    controls = [Vector(dimension, row) for row in coords]
    path = bezier.BezierPath(); path.setControlPoints(controls)
    path.minimum_sqr_distance = 1e-4
    a, b = interval
    def reference(t):
        return [.75*t, 6*t*(1-t)*(1-2*t)] + ([-1.5*t] if dimension == 3 else [])
    start, end = Vector(dimension, reference(a)), Vector(dimension, reference(b))
    points = [start, end]
    count = path.findDrawingPointsAdded(0, a, b, points, 1)
    assert count == len(points)-2
    assert points[0] is start and points[-1] is end
    if a == b:
        assert count == 0
    else:
        assert all(p.vector[0] < q.vector[0] for p, q in zip(points, points[1:]))
    for point in points[1:-1]:
        assert point.vector == pytest.approx(reference(point.vector[0]/.75), rel=0, abs=3e-14)
        assert all(point.vector is not row for row in coords)
    assert all(p.vector is row for p, row in zip(controls, coords))


@pytest.mark.parametrize('dimension', [2, 3])
def test_interpolated_noncollinear_controls_and_geometric_join(dimension):
    # Interior tangent follows (3,4[,0])/5. Adjacent distances are 3 and 4.
    coords = [[0., 0.], [3., 0.], [3., 4.]]
    if dimension == 3:
        coords = [row+[2.] for row in coords]
    source = [Vector(dimension, row) for row in coords]
    snapshot = [row[:] for row in coords]
    path = bezier.BezierPath()
    assert path.interpolate(source, .5) is None
    expected = [[0, 0], [1.5, 0], [2.1, -1.2], [3, 0], [4.2, 1.6], [3, 2], [3, 4]]
    if dimension == 3:
        expected = [row+[2] for row in expected]
    for actual, reference in zip(path.controlPoints, expected):
        assert actual.vector == pytest.approx(reference, rel=0, abs=2e-15)
        assert all(actual.vector is not row for row in coords)
    assert len(path.controlPoints) == 7 and path.curveCount == 2
    # Same tangent direction, different parameter speeds: G1, not promised C1.
    incoming = [3*(b-a) for a, b in zip(expected[2], expected[3])]
    outgoing = [3*(b-a) for a, b in zip(expected[3], expected[4])]
    assert outgoing == pytest.approx([4*x/3 for x in incoming], rel=0, abs=2e-15)
    assert coords == snapshot


@pytest.mark.parametrize('dimension', [None, 2, 3])
def test_builder_repeated_coincident_controls_and_returned_storage(dimension):
    point = 2. if dimension is None else Vector(dimension, [2.]*dimension)
    source = [point]*3
    path = bezier.BezierPath()
    path.samplePoints(source, 0., 1., .5)
    first = path.controlPoints
    expected = [2.]*7 if dimension is None else [[2.]*dimension]*7
    def coordinates():
        return path.controlPoints if dimension is None else [p.vector for p in path.controlPoints]
    assert coordinates() == expected
    path.samplePoints(source, 0., 1., .5)
    assert path.controlPoints is not first and coordinates() == expected
    for i in range(2):
        points = path.findDrawingPoints(i)
        assert len(points) == 2
        if dimension is not None:
            assert points[0].vector is not points[1].vector
            points[0].vector[0] = 99
            assert point.vector == [2.]*dimension and coordinates() == expected


@pytest.mark.parametrize('axis', [1, 2, 3])
@pytest.mark.parametrize('t', [0., .25, .5, .75, 1.])
def test_existing_squad4_spline_same_axis_polynomial(axis, t):
    # Positive same-axis controls stay within one shortest-path branch.
    # Half-angles: endpoints 0, pi/3; controls pi/6, pi/4.
    angles = [0., math.pi/3, math.pi/6, math.pi/4]
    inputs = []
    for angle in angles:
        data = [math.cos(angle), 0., 0., 0.]
        data[axis] = math.sin(angle)
        inputs.append(Quaternion(data))
    saved = [q.data[:] for q in inputs]; storage = [q.data for q in inputs]
    blend = 2*t*(1-t)
    half_angle = (1-blend)*t*math.pi/3 + blend*((1-t)*math.pi/6+t*math.pi/4)
    expected = [math.cos(half_angle), 0., 0., 0.]
    expected[axis] = math.sin(half_angle)
    result = squad4(*inputs, t)
    assert result.data == pytest.approx(expected, rel=0, abs=3e-15)
    assert math.hypot(*result.data) == pytest.approx(1., rel=0, abs=3e-15)
    assert all(q.data is data for q, data in zip(inputs, storage))
    assert [q.data for q in inputs] == saved
    assert all(result.data is not data for data in storage)


@pytest.mark.parametrize('axis', [1, 2, 3])
@pytest.mark.parametrize('t', [0., .5, .75, 1.])
def test_existing_legacy_squad_spline_known_spherical_branches(axis, t):
    # q0=identity, control half-angle pi/4, endpoint half-angle pi/2.
    # At .5/.75 all nontrivial blends use the documented spherical branch;
    # the final parameter is zero at endpoints. No normalization is added.
    inputs = []
    for angle in [0., math.pi/4, math.pi/2]:
        data = [math.cos(angle), 0., 0., 0.]
        data[axis] = math.sin(angle)
        inputs.append(Quaternion(data))
    storage = [q.data for q in inputs]; saved = [data[:] for data in storage]
    half_angle = t*math.pi/2 - 2*t*(1-t)*t*math.pi/4
    expected = [math.cos(half_angle), 0., 0., 0.]
    expected[axis] = math.sin(half_angle)
    result = inputs[0].squad(inputs[1], inputs[2], t)
    assert result.data == pytest.approx(expected, rel=0, abs=3e-15)
    assert all(q.data is data for q, data in zip(inputs, storage))
    assert [q.data for q in inputs] == saved
    assert all(result.data is not data for data in storage)


@pytest.mark.parametrize('degree', [0, 1, 2, 3, 7, 12, 20, 32, 64, 128])
@pytest.mark.parametrize('x', [-1., -.75, -.25, 0., .375, .875, 1.])
def test_legendre_all_selected_orders_rodrigues_and_state(degree, x):
    for order in sorted({0, min(1, degree), degree//2, degree}):
        polynomial = Legendre(degree, order, x)
        polynomial.P, polynomial.PM1, polynomial.PML = -123., 456., -789.
        snapshot = polynomial.__dict__.copy()
        expected = legendre_reference(degree, order, x)
        # Near zeros use a scale from the same degree/order at x=0 or 1/4,
        # not an absolute tolerance that would hide tiny associated results.
        scale = max(abs(legendre_reference(degree, order, 0.)),
                    abs(legendre_reference(degree, order, .25)), abs(expected))
        assert polynomial.run() == pytest.approx(expected, rel=2e-12, abs=scale*3e-14)
        assert polynomial.run() == polynomial.run()
        assert polynomial.__dict__ == snapshot
        assert Legendre(degree, order, -x).run() == (-1)**(degree+order)*polynomial.run()


@pytest.mark.parametrize('degree', [1, 2, 12, 64, 256, 1000])
def test_high_degree_endpoints_and_near_endpoint_reference(degree):
    assert Legendre(degree, 0, 1.).run() == 1.
    assert Legendre(degree, 0, -1.).run() == (-1)**degree
    # Exact finite Taylor expansion about 1, avoiding cancellation of the
    # expanded Rodrigues polynomial at a high-degree endpoint.
    x = math.nextafter(1., 0.)
    delta = Fraction(x)-1
    expected = sum(Fraction(math.factorial(degree+k),
                            2**k*math.factorial(k)**2*math.factorial(degree-k))*delta**k
                   for k in range(degree+1))
    # Endpoint recurrence roundoff accumulates with degree. This exploratory
    # O(l^2 ulp) budget is not a public high-degree accuracy guarantee.
    assert Legendre(degree, 0, x).run() == pytest.approx(
        float(expected), rel=0, abs=2*degree*(degree+1)*math.ulp(1.))


@pytest.mark.parametrize('degree,order', [(0,0), (1,0), (1,1), (2,1), (5,2), (8,0), (8,8)])
def test_unnormalized_associated_integral(degree, order):
    # Exact orthogonality normalization; Simpson is independent of recurrence.
    count = 4096
    values = [Legendre(degree, order, -1+2*i/count).run()**2 for i in range(count+1)]
    integral = (values[0]+values[-1] + 4*math.fsum(values[1:-1:2])
                + 2*math.fsum(values[2:-1:2]))*(2/count)/3
    expected = float(Fraction(2*math.factorial(degree+order),
                              (2*degree+1)*math.factorial(degree-order)))
    assert integral == pytest.approx(expected, rel=2e-8, abs=1e-14)


@pytest.mark.parametrize('degree,order,x', [(0,0,.2), (1,1,-.75), (7,3,.375), (32,12,-.25)])
def test_scratch_helper_target_degree_repeatability(degree, order, x):
    polynomial = Legendre(degree, order, x)
    for target in [order, order+1, max(degree, order+2)]:
        for _ in range(2):
            assert polynomial.calculatePML(target) is None
            assert polynomial.PML == pytest.approx(legendre_reference(target, order, x), rel=3e-13)
            assert (polynomial.l, polynomial.m, polynomial.x) == (degree, order, x)


@pytest.mark.parametrize('degree,order', [(2,0), (3,1), (12,5), (32,12), (64,0)])
def test_associated_recurrence_with_independent_right_hand_side(degree, order):
    x = -.375
    expected = ((2*degree+1)*x*legendre_reference(degree, order, x)
                - (degree+order)*legendre_reference(degree-1, order, x))
    actual = (degree-order+1)*Legendre(degree+1, order, x).run()
    assert actual == pytest.approx(expected, rel=3e-12, abs=abs(expected)*3e-14)


@pytest.mark.parametrize('degree', range(9))
@pytest.mark.parametrize('x', [-2., -1.25, 1.25, 2.])
def test_ordinary_extrapolation_independent_polynomial(degree, x):
    assert Legendre(degree, 0, x).run() == pytest.approx(legendre_reference(degree, 0, x), rel=3e-14)


@pytest.mark.parametrize('x', [-1e300, -1e-300, -math.ldexp(1., -1074),
                              math.ldexp(1., -1074), 1e-300, 1e300])
def test_low_degree_extreme_arguments_with_representable_answers(x):
    assert Legendre(0, 0, x).run() == 1.
    assert Legendre(1, 0, x).run() == x


@pytest.mark.parametrize('degree,order', [(1,1), (2,1), (12,7), (32,12)])
@pytest.mark.parametrize('x', [math.nextafter(-1., 0.), math.nextafter(1., 0.)])
def test_associated_near_pole_factored_seed(degree, order, x):
    assert Legendre(degree, order, x).run() == pytest.approx(
        legendre_reference(degree, order, x), rel=3e-14, abs=0)


@pytest.mark.parametrize('x', [-1e154, 1e154])
def test_low_degree_extrapolation_avoids_intermediate_overflow(x):
    expected = legendre_reference(2, 0, x)
    assert math.isfinite(expected) and expected > 1e308
    polynomial = Legendre(2, 0, x)
    snapshot = polynomial.__dict__.copy()
    actual = polynomial.run()
    assert polynomial.__dict__ == snapshot
    assert actual == pytest.approx(expected, rel=3e-15)


def test_existing_domain_and_compatibility_characterization():
    from gem.experimental.bezier import BezierPath, quadraticBezierPoint
    from gem.experimental.legendre import Legendre as legacy
    assert BezierPath is bezier.BezierPath and quadraticBezierPoint is bezier.quadraticBezierPoint
    assert legacy is Legendre
    with pytest.raises(ValueError):
        Legendre(2, 1, 1.01).run()
    with pytest.raises(TypeError):
        Legendre(4, 1.5, .2).run()
    invalid = Legendre(0, 1, .2)
    invalid.PML = 123.
    assert invalid.run() == 123.  # Historical behavior, not a new invalid policy.
